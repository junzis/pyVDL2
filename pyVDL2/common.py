import numpy as np


def hex2bin(hexstr):
    """Convert a hexdecimal string to binary string, with zero fillings."""
    num_of_bits = len(hexstr) * 4
    binstr = bin(int(hexstr, 16))[2:].zfill(int(num_of_bits))
    return binstr


def hex2int(hexstr):
    """Convert a hexdecimal string to integer."""
    return int(hexstr, 16)


def int2hex(n):
    """Convert a integer to hexadecimal string."""
    # strip 'L' for python 2
    return hex(n)[2:].rjust(6, "0").upper().rstrip("L")


def bin2int(binstr):
    """Convert a binary string to integer."""
    return int(binstr, 2)


def bin2hex(hexstr):
    """Convert a hexdecimal string to integer."""
    return int2hex(bin2int(hexstr))


def bin2np(binstr):
    """Convert a binary string to numpy array."""
    return np.array([int(i) for i in binstr])


def np2bin(npbin):
    """Convert a binary numpy array to string."""
    return np.array2string(npbin, separator="")[1:-1]


def df(msg):
    """Decode Downlink Format vaule, bits 1 to 5."""
    msgbin = hex2bin(msg)
    return min(bin2int(msgbin[0:5]), 24)
