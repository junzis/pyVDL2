from . import common
from .acars import decode_acars


def decode_address(hexaddr):
    """Deocde address from the hexadecimal string.

    Args:
        hexaddr (String): 16-char (8-byte or 64-bit) hexadecimal string

    Returns:
        dict: docoded source and destination addresses

    """

    def _decode(hex8):
        abin = common.hex2bin(hex8)  # 32 bits
        newbin = abin[24:31] + abin[16:23] + abin[8:15] + abin[0:6]
        newbin = newbin[::-1]

        addr_type_bits = newbin[:3]
        if addr_type_bits == "111":
            addr = "FFFFFF"
            addr_type = "all stations"
            addr_comment = "broadcast"
        elif addr_type_bits == "001":
            addr = common.bin2hex(newbin[3:])
            addr_type = "aircraft"
            addr_comment = "ICAO address"
        elif addr_type_bits == "100":
            addr = common.bin2hex(newbin[3:])
            addr_type = "ground station"
            addr_comment = "ICAO administered"
        elif addr_type_bits == "101":
            addr = common.bin2hex(newbin[3:])
            addr_type = "ground station"
            addr_comment = "ICAO delegated"
        else:
            addr = common.bin2hex(newbin[3:])
            addr_type = "reserved"
            addr_comment = "future use"

        return {"address": addr, "type": addr_type, "comment": addr_comment}

    return {"source": _decode(hexaddr[8:]), "destination": _decode(hexaddr[:8])}


def decode_linkcontrol(hexlc):
    """Deocde link control information from the hexadecimal string.

    Args:
        hexlc (String): 2-char (1-byte or 8-bit) hexadecimal string

    Returns:
        dict: link control information

    """
    binstr = common.hex2bin(hexlc)[::-1]
    if binstr[0] == "0":
        linkcontrol = {
            "raw": binstr,
            "type": "information",
            "send_sequence": int(binstr[1:4], 2),
            "require_response": bool(int(binstr[4])),
            "receive_sequence": int(binstr[5:8], 2),
        }
    else:
        if binstr[1] == "0":
            linkcontrol = {
                "raw": binstr,
                "type": "supervisory",
                "function_bits": binstr[2:4],
                "require_response": bool(int(binstr[4])),
                "receive_sequence": int(binstr[5:8], 2),
            }
        else:
            linkcontrol = {
                "raw": binstr,
                "type": "unnumbered",
                "function_bits_1": binstr[2:4],
                "require_response": bool(int(binstr[4])),
                "function_bits_2": binstr[5:8],
            }

    return linkcontrol


def decode_avlc(hexframe):
    """Deocde AVLC information from the ALVC frame in hexadecimal formats.

    Args:
        hexframe (String): hexadecimal string

    Returns:
        dict: decoded VDL2 information

    """
    if len(hexframe) < 18:
        return None

    addr = decode_address(hexframe[:16])
    link_control = decode_linkcontrol(hexframe[16:18])

    info = hexframe[18:-12]
    fcs = hexframe[-12:-4]
    flag = hexframe[-4:]

    general_format = "N/A"
    information = info

    if info is None or len(info) < 4:
        pass

    elif common.hex2bin(info[0]) == "1111":
        # ACARS over AVLC
        data = info[6:]
        general_format = "ACARS"
        information = decode_acars(data)

    elif common.hex2bin(info[0]) == "0001" and len(info) > 6:
        # X.25 (ISO 8208)
        general_format_identifier = common.hex2bin(info[0])
        logical_channel_group_number = common.hex2bin(info[1])
        logical_channel_number = common.hex2bin(info[2:4])
        packet_data_fields = common.hex2bin(info[4:6])

        if packet_data_fields == "00001011":
            packet_type = "call request"
        elif packet_data_fields == "00001111":
            packet_type = "call accepted"
        elif packet_data_fields == "00010011":
            packet_type = "clear request"
        elif packet_data_fields == "00010111":
            packet_type = "clear confirmation"

        elif packet_data_fields[-1] == "0":
            packet_type = "data"
        elif packet_data_fields == "00100011":
            packet_type = "interupt"
        elif packet_data_fields == "00100111":
            packet_type = "interupt confirmation"

        elif packet_data_fields[-5:] == "00001":
            packet_type = "receive ready"
        elif packet_data_fields[-5:] == "00101":
            packet_type = "receive not ready"
        elif packet_data_fields[-5:] == "01001":
            packet_type = "reject"
        elif packet_data_fields == "00011011":
            packet_type = "reset request"
        elif packet_data_fields == "00011111":
            packet_type = "reset confirmation"

        elif packet_data_fields == "11111011":
            packet_type = "restart request"
        elif packet_data_fields == "11111111":
            packet_type = "restart confirmation"

        elif packet_data_fields == "11110001":
            packet_type = "diagnostic"

        else:
            packet_type = "unknown"

        general_format = "X.25"
        information = (
            {
                "packet_type": packet_type,
                "info": info,
                "raw_header_info": {
                    "general_format_identifier": general_format_identifier,
                    "logical_channel_group_number": logical_channel_group_number,
                    "logical_channel_number": logical_channel_number,
                    "packet_data_fields": packet_data_fields,
                },
            },
        )

    elif link_control["type"] == "unnumbered":
        if (
            link_control["function_bits_1"] == "11"
            and link_control["function_bits_2"] == "101"
        ):
            general_format = "XID"
        if (
            link_control["function_bits_1"] == "00"
            and link_control["function_bits_2"] == "001"
        ):
            general_format = "TEST"

        information = info

    decoded = {
        "address": addr,
        "link_control": link_control,
        "general_format": general_format,
        "information": information,
        "fcs": fcs,
        "flag": flag,
    }

    return decoded
