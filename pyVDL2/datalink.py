import numpy as np
import logging
from scipy import signal
from textwrap import wrap

from . import common
from .librs import RSCodec


max_len = 0x3FFF
max_len_corrected = 0x1FFF


class DataLink(object):
    """docstring for DataLink."""

    def __init__(self):
        super(DataLink, self).__init__()
        # linear-feedback shift register
        self.rscodec = RSCodec(nsym=6, nsize=255, fcr=120, prim=0x187)

    def reset_lfsr(self):
        self.lfsr = [1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1]

    def getlen(self, header_bitstr):
        # 3 bits of reserved symbol (000)
        # 17 bits transmission length
        # 5 bits header FEC

        header0 = header_bitstr[:25]
        logging.debug("FEC: {}".format(header0[20:25]))

        length0 = int(header0[3:20][::-1], 2)
        logging.debug(
            "Transmission length (before correction): {}".format(length0)
        )

        if length0 > max_len:
            logging.debug("Message length seems too long, should not decode.")
            length = None

        header, syndrome = self.header_correction(header0)

        length = int(header[3:20][::-1], 2)
        logging.debug("FEC syndrome: {}".format(syndrome))
        logging.debug(
            "Transmission length (after correction): {}".format(length)
        )

        # drop frame with very large length
        if syndrome > 0 and length > max_len_corrected:
            logging.debug("Message length seems too long, should not decode.")
            length = None

        return length

    def decode_data(self, length, databin, strout=True):
        n_block = int(np.ceil(length / (249 * 8)))

        logging.debug("Length of data: {} bits".format(length))
        logging.debug("Number of blocks: {}".format(n_block))

        # --- De-interleavering ---
        databin = databin + "0" * (8 - len(databin) % 8)

        data_octs = [d8[::-1] for d8 in wrap(databin, 8)]

        de_inter_matrix = np.ones((n_block, 255), dtype=int)

        last_block_content_bits = length % (249 * 8)
        last_block_content_octs = int(np.ceil(last_block_content_bits / 8))

        if last_block_content_octs >= 68:
            de_inter_matrix[-1][last_block_content_octs:249] = 0
            last_block_rs_erase = []
        elif last_block_content_octs >= 31:
            de_inter_matrix[-1][last_block_content_octs:249] = 0
            de_inter_matrix[-1][249 + 4 :] = 0
            last_block_rs_erase = [224, 225]
        elif last_block_content_octs >= 3:
            de_inter_matrix[-1][last_block_content_octs:249] = 0
            de_inter_matrix[-1][249 + 2 :] = 0
            last_block_rs_erase = [222, 223, 224, 225]
        else:
            de_inter_matrix[-1][last_block_content_octs:] = 0
            last_block_rs_erase = None

        i = 0
        for c in range(255):
            for r in range(n_block):
                if de_inter_matrix[r, c] == 0:
                    continue
                de_inter_matrix[r, c] = int(data_octs[i], 2)
                i += 1

        blocks = []
        for r in range(n_block):
            block = bytes(de_inter_matrix[r].tolist())
            blocks.append(block)

        erasures = [[]] * n_block
        erasures[-1] = last_block_rs_erase

        block_content_sizes = [249] * n_block
        block_content_sizes[-1] = last_block_content_octs

        # Contruct original AVLC frame

        error = False
        msgints = []

        for i, (block, size, erasure) in enumerate(
            zip(blocks, block_content_sizes, erasures)
        ):
            logging.debug("Before RS check:\n" + common.dumphex(block))
            if erasure is None:
                corrected = block
                msgints.extend(corrected[:size])
            else:
                try:
                    corrected = self.rscodec.decode(block, erase_pos=erasure)
                    logging.debug("RS check, block {}, success".format(i + 1))
                    logging.debug(
                        "After RS check:\n" + common.dumphex(corrected)
                    )
                    msgints.extend(corrected[:size])
                except Exception as err:
                    logging.debug(
                        "RS check, block {}, error: {}".format(i + 1, err)
                    )
                    error = True
                    pass

        if error:
            return None

        else:
            msg_stuffed = "".join(
                [format(i, "08b")[::-1] for i in msgints[1:]]
            )
            msg_unstuff = msg_stuffed[: length - 8].replace("111110", "11111")
            message = "".join([b8[::-1] for b8 in wrap(msg_unstuff, 8)])

            logging.info(
                "AVLC message (after un-stuffing):\n" + common.dumphex(message)
            )
            message = common.bin2hex(message)

            return message

    def descramble(self, bitstream):
        # descramble raw binary bit string
        scrambled = [int(i) for i in bitstream]
        descrambled = []
        i = 0
        while i < len(scrambled):
            bit = self.lfsr[0] ^ self.lfsr[-1]
            self.lfsr = [bit] + self.lfsr[0:-1]
            descrambled.append(scrambled[i] ^ bit)
            i += 1

        descrambled = "".join(map(str, descrambled))
        return descrambled

    def header_correction(self, headerbits):
        def parity(x):
            p = 0
            while x:
                p = 0 if p is 1 else 1
                x = x & (x - 1)
            return p

        headerint = int(headerbits, 2)
        syndrome = 0
        row = 0

        for i, h in enumerate(common.fec_H):
            row = headerint & h
            syndrome |= parity(row) << (5 - 1 - i)

        headerint ^= common.fec_syndrome_table[syndrome]

        headerbits_corr = format(headerint, "025b")
        return headerbits_corr, syndrome
