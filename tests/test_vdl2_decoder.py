import os
import numpy as np
import logging
from pyVDL2 import VDL2Decoder
from pyVDL2 import common

# logging.basicConfig(level=logging.DEBUG)


dir = os.path.dirname(os.path.realpath(__file__))

sample = dir + "/data/sample_iqdata_rtlsdr.dat"

iqdata = ((np.fromfile(sample, dtype="uint8") - 127.5)) / 255


decoder = VDL2Decoder(sample_rate=1050000, channel_freq=136975000)

count = 0
for t, ppm, msg in decoder.process(iqdata, t0=0):
    count += 1

    print()
    print("[message:{} | Timestamp:{}]----------".format(count, round(t, 6)))
    print(common.dumphex(msg))
    print()
