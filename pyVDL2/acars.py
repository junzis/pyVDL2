from . import common
import string
from textwrap import wrap


def decode_acars(hexstr):
    databin = common.hex2bin(hexstr)

    # replace first bit in each octect with "0"
    strchunks = ["0" + d8[1:] for d8 in wrap(databin, 8)]
    acars_frame = [int(i, 2).to_bytes(1, "big").decode() for i in strchunks]
    # print(acars_frame)

    regid = "".join(acars_frame[1:8])

    if acars_frame[8] == "\x06":
        ack = "ACK"
    elif acars_frame[8] == "\x15":
        ack = "NAK"
    else:
        ack = acars_frame[8]

    label = "".join(acars_frame[9:11])
    label = label.replace("\x7f", "DEL")

    block_id = acars_frame[11]

    # acars_frame[12] is '\x02', STX (start of text) in ascii code

    msg_num = "".join(acars_frame[13:17])

    flight = "".join(acars_frame[17:23])

    msg = "".join(acars_frame[23:])
    msg_printable = "".join(filter(lambda x: x in string.printable, msg))

    res = {
        "mode": acars_frame[0],
        "regid": regid,
        "ack": ack,
        "label": label,
        "block_id": block_id,
        "msg_num": msg_num,
        "flight": flight,
        "msg": msg_printable,
    }

    return res
