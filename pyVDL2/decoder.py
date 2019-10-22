import numpy as np
import time
import logging
from scipy import signal
from .datalink import DataLink
from . import common

# pre-compute constant terms for linear regression in sync search errors
lr_x = np.arange(16) - np.mean(np.arange(16))
lr_denom = np.sum(lr_x ** 2)
err_threshold = 4


class VDL2Decoder(object):
    """docstring for VDL2Decoder."""

    def __init__(self, sample_rate, channel_freq):
        super(VDL2Decoder, self).__init__()
        self.sample_rate = sample_rate
        self.channel_freq = channel_freq

        # samples per symbol
        self.sps = sample_rate // common.symbol_rate

        # Chebyshev low pass filter
        self.sos = signal.cheby2(2, 45, 0.22, "low", output="sos")

        # steps used for searching, to reduce computation load
        self.signal_search_step = self.sps * 10
        self.sync_search_step = max(self.sps // 20, 1)
        self.minimum_noise = 0.02

        self.datalink = DataLink()

    def reset_sync_params(self):
        self.sync_start = None
        self.dph_init = None
        self.ppm_err = None
        self.prev_dphi = None
        self.ph_errs = [1e6, 1e6, 1e6]

    def process(self, iqdata, t0=None):
        """Process raw iqdata."""

        if t0 is None:
            t0 = time.time()

        self.re = iqdata[:-1:2]
        self.im = iqdata[1::2]

        self.lfre = signal.sosfilt(self.sos, self.re)
        self.lfim = signal.sosfilt(self.sos, self.im)

        self.phase = np.arctan2(self.lfim, self.lfre)
        self.magnitude = np.sqrt(self.lfim ** 2 + self.lfre ** 2)

        self.idx_curr = 0

        n = len(self.phase)
        self.reset_sync_params()

        self.noise_floor = 0.01
        self.search_sync_count = 0

        while self.idx_curr < n - self.signal_search_step:
            next_mag = self.magnitude[
                self.idx_curr : self.idx_curr + self.signal_search_step
            ]

            # skipping if signal is not higher than twice of noise floor
            if np.mean(next_mag) < self.noise_floor * 2:
                self.idx_curr += self.signal_search_step
                self.search_sync_count = 0
                continue

            # skipping if the header can not be synced within 20 symbols
            elif (
                self.search_sync_count > self.sps * 20 / self.sync_search_step
            ):
                self.idx_curr += self.sps
                continue

            else:
                self.search_sync_count += 1

            # find the sync point
            errvec = np.empty(common.training_len)
            unwrap = 0
            prev_err = (
                self.phase[self.idx_curr + self.sps]
                - common.training_phases[0]
            )
            errvec[0] = prev_err

            for i in range(1, common.training_len):
                cur_err = (
                    self.phase[self.idx_curr + (i + 1) * self.sps]
                    - common.training_phases[i]
                )
                errdiff = cur_err - prev_err

                if errdiff > np.pi:
                    unwrap -= 2 * np.pi
                elif errdiff < -np.pi:
                    unwrap += 2 * np.pi

                errvec[i] = cur_err + unwrap

                prev_err = cur_err

            # compute the squared sum of phase error
            lr_y = errvec - np.mean(errvec)
            freq_err = np.sum(lr_x * lr_y) / lr_denom
            self.ph_errs[0] = np.sum((lr_y - freq_err * lr_x) ** 2)

            # found the sync pattern
            if (
                self.ph_errs[1] < err_threshold
                and self.ph_errs[0] > self.ph_errs[1]
            ):
                # find the sync point where the error is the smallest
                ls = self._find_min_shift(self.sync_search_step, self.ph_errs)
                self.sync_start = self.idx_curr - ls
                self.dph_init = self.prev_dphi
                self.ppm_err = (
                    common.symbol_rate
                    * self.dph_init
                    / (2 * np.pi * self.channel_freq)
                    * 1e6
                )

                # start demodulating and decoding
                msg = self.demod_and_decode()

                if msg is not None:
                    t = t0 + self.sync_start / self.sample_rate
                    yield t, self.ppm_err, msg

                # update noise_floor
                self.noise_floor = np.mean(
                    self.magnitude[
                        self.sync_start
                        - 2 * self.signal_search_step : self.sync_start
                        - self.signal_search_step
                    ]
                )

                # reset sync parameters for next messsage
                self.reset_sync_params()

                self.search_sync_count = 0

            self.ph_errs[2] = self.ph_errs[1]
            self.ph_errs[1] = self.ph_errs[0]
            self.prev_dphi = freq_err

            self.idx_curr += self.sync_search_step

    def demod_and_decode(self):
        header_start = self.sync_start + 16 * self.sps
        self.idx_curr = header_start

        logging.debug("[-- New chunk --]")
        logging.debug("Data starts at: {}".format(header_start))
        logging.debug("Initial phase shift: {}".format(self.dph_init))
        logging.debug("PPM error: {}".format(self.ppm_err))

        header_phases = self.phase[
            header_start : header_start + 10 * self.sps : self.sps
        ]

        header = []
        for i in range(1, len(header_phases)):
            dph = header_phases[i] - header_phases[i - 1] - self.dph_init
            if dph < 0:
                dph += 2 * np.pi
            elif dph > 2 * np.pi:
                dph -= 2 * np.pi
            dph = dph / (np.pi / 4)
            i = int(round(dph)) % 8
            header.append(common.graycode[i])

        logging.debug("Headr:{}".format("".join(map(str, header))))

        headerbin = "".join([format(c, "03b") for c in header])

        # reset LFSR for new message
        self.datalink.reset_lfsr()
        headerbin = self.datalink.descramble(headerbin)
        length = self.datalink.getlen(headerbin)  # number of bit

        if length is None:
            return None
        else:
            data_size = self.calc_data_symbol_size(length)
            datastart = header_start + (10 - 1) * self.sps
            dataend = header_start + (10 + data_size + 1) * self.sps

            data_phases = self.phase[datastart : dataend : self.sps]

            data = []
            for i in range(1, len(data_phases)):
                dph = data_phases[i] - data_phases[i - 1] - self.dph_init
                if dph < 0:
                    dph += 2 * np.pi
                elif dph > 2 * np.pi:
                    dph -= 2 * np.pi
                dph = dph / (np.pi / 4)
                i = int(round(dph)) % 8
                data.append(common.graycode[i])

            logging.debug("Data:\n{}".format(common.dumpdemod(data)))

            databin = "".join([format(c, "03b") for c in data])

            # join with left over from header octects
            databin = headerbin[-2:] + self.datalink.descramble(databin)

            msg_decoded = self.datalink.decode_data(length, databin)

            if msg_decoded is not None:
                self.idx_curr = dataend

            return msg_decoded

    def calc_data_symbol_size(self, length):
        full_block_octets = int(length / 8 / 249) * 255
        last_block_octets = (length / 8) % 249

        if last_block_octets <= 2:
            total_octets = full_block_octets + last_block_octets
        elif 3 < last_block_octets <= 30:
            total_octets = full_block_octets + last_block_octets + 2
        elif 31 < last_block_octets <= 67:
            total_octets = full_block_octets + last_block_octets + 4
        else:
            total_octets = full_block_octets + last_block_octets + 6

        n_symbols = int(np.ceil(total_octets * 8 / 3))

        return n_symbols

    def _find_min_shift(self, d, ys):
        # return the x = min(y) from a parabolic function
        # func: y = a(x-b)^2+c
        # points: (0, y1), (-d, y2) (-2d, y3)
        # b always < 0
        y3, y2, y1 = ys
        b = -d / 2 * (y3 - y2) / (y1 + y3)
        left_shift = int(round(-b))
        return left_shift
