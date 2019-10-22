import os
import json
from pyVDL2.avlc import decode_avlc

frames = []

dir = os.path.dirname(os.path.realpath(__file__))

with open(dir + "/data/sample_avlc_frame.txt") as f:
    frame = "".join(f.readlines())
    frame = frame.replace(" ", "").replace("\n", "")

    decoded = decode_avlc(frame)
    print(json.dumps(decoded, indent=2))
