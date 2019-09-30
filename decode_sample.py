import pyModeS as pms
from avlc import decode_avlc
import json

frames = []
ts = []
with open("data/sample.txt") as f:
    buffer = ""
    t = ""
    for line in f:
        if "[2019-" in line:
            t = line[1:20]

        elif "|" in line:
            buffer += line.split("|")[0]

        elif line == "\n":
            if not buffer:
                continue
            else:
                frame = buffer.replace(" ", "")
                frames.append(frame.upper())
                ts.append(t)
                buffer = ""
                t = ""
        else:
            continue


label_stats = {}
for t, frame in zip(ts, frames):
    decoded = decode_avlc(frame)
    if "general_format" in decoded:
        print("-" * 30, t, "-" * 30)
        print(json.dumps(decoded, indent=2))
