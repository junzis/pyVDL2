import numpy as np
from textwrap import wrap
from collections.abc import Iterable


training_phases = [
    0 * np.pi / 4,
    3 * np.pi / 4,
    -3 * np.pi / 4,
    1 * np.pi / 4,
    1 * np.pi / 4,
    2 * np.pi / 4,
    0 * np.pi / 4,
    4 * np.pi / 4,
    -3 * np.pi / 4,
    4 * np.pi / 4,
    -2 * np.pi / 4,
    3 * np.pi / 4,
    1 * np.pi / 4,
    -2 * np.pi / 4,
    -3 * np.pi / 4,
    0 * np.pi / 4,
]

training_len = 16

graycode = [0, 1, 3, 2, 6, 7, 5, 4]

symbol_rate = 10500

fec_H = [
    0b0000000011111111111110000,
    0b0011111100001111111101000,
    0b1100011100110000111100100,
    0b1101101101010011001100010,
    0b0110100111100101010100001,
]

fec_syndrome_table = [
    0b0000000000000000000000000,
    0b0000000000000000000000001,
    0b0000000000000000000000010,
    0b0100000000000000000000100,
    0b0000000000000000000000100,
    0b0100000000000000000000010,
    0b1000000000000000000000000,
    0b0100000000000000000000000,
    0b0000000000000000000001000,
    0b0010000000000000000000000,
    0b0001000000000000000000000,
    0b0000100000000000000000000,
    0b0000010000000000000000000,
    0b1000100000000000000000000,
    0b0000001000000000000000000,
    0b0000000100000000000000000,
    0b0000000000000000000010000,
    0b0000000010000000000000000,
    0b0100000000100000000000000,
    0b0000000001000000000000000,
    0b0100000001000000000000000,
    0b0000000000100000000000000,
    0b0000000000010000000000000,
    0b1000000010000000000000000,
    0b0000000000001000000000000,
    0b0000000000000100000000000,
    0b0000000000000010000000000,
    0b0000000000000001000000000,
    0b0000000000000000100000000,
    0b0000000000000000010000000,
    0b0000000000000000001000000,
    0b0000000000000000000100000,
]


def hex2bin(hexstr):
    """Convert a hexdecimal string to binary string, with zero fillings."""
    n_bits = len(hexstr) * 4
    binstr = "{i:0{n}b}".format(i=int(hexstr, 16), n=n_bits)
    return binstr


def bin2hex(data):
    if not isinstance(data, Iterable):
        raise Exception("Data is not iterable: {}".format(data))

    if isinstance(data, str) and set(data) == {"0", "1"}:
        binstr = data
    elif set(data) == {0, 1} or set(data) == {"0", "1"}:
        binstr = "".join(map(str, data))
    else:
        raise Exception("Can not convert this to binary: {}".format(data))

    if len(binstr) % 8 > 0:
        binstr = binstr + "0" * (8 - len(binstr) % 8)

    hexstr = "".join([format(int(b8, 2), "02x") for b8 in wrap(binstr, 8)])

    return hexstr


def dumphex(data):

    if isinstance(data, Iterable) and set(data).issubset({0, 1}):
        data = "".join(map(str, data))
        if len(data) % 8 > 0:
            data = data + "0" * (8 - len(data) % 8)
        data = [int(b8, 2) for b8 in wrap(data, 8)]

    elif isinstance(data, Iterable) and set(data).issubset({"0", "1"}):
        if len(data) % 8 > 0:
            data = data + "0" * (8 - len(data) % 8)
        data = [int(b8, 2) for b8 in wrap(data, 8)]

    elif isinstance(data, Iterable) and set(data).issubset(
        set("0123456789abcdefABCDEF")
    ):
        if len(data) % 2 > 0:
            data = data + "0"
        data = [int(h2, 16) for h2 in wrap(data, 2)]

    res = ""
    for i, b in enumerate(data):
        if i > 0 and i % 10 == 0:
            res += " "

        if i > 0 and i % 20 == 0:
            res += "\n"

        if isinstance(b, bytes):
            b = ord(b)

        # print(b)
        hexstr = format(b, "02x")
        res = res + hexstr + " "

    return res


def dumpdemod(data):
    if isinstance(data[0], int):
        intstr = "".join(map(str, data))
    elif isinstance(data[0], str):
        intstr = data

    res = ""
    for i, s in enumerate(intstr):
        if i > 0 and i % 5 == 0:
            res += " "

        if i > 0 and i % 50 == 0:
            res += "\n"

        res = res + s

    return res
