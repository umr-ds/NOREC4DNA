import struct
from typing import Any, Callable, List, Optional, Union

from bitarray import bitarray
from norec4dna.helper import calc_crc
from reedsolo import RSCodec

try:
    missing_crclib = False
except ImportError:
    missing_crclib = True

rscodec: Optional[RSCodec] = None
number_r_symbols: int = 2
ErrorCorrectionCallable = Callable[..., bytes]


def nocode(txt: bytes, *args: Any) -> bytes:
    return txt


def crc32(txt: Union[bytes, bytearray], crc_len_format: str = "B") -> bytes:
    crc = calc_crc(bytes(txt), crc_len_format)
    packed = struct.pack("<" + str(len(txt)) + "s" + crc_len_format, txt, crc)
    return packed


def crc32_decode(txt: bytes, crc_len_format="B") -> bytes:
    crc_len = -struct.calcsize("<" + crc_len_format)
    crc = struct.unpack("<" + crc_len_format, txt[crc_len:])[0]
    payload = txt[:crc_len]
    calced_crc = calc_crc(payload, crc_len_format)
    assert crc == calced_crc, "CRC-Error - " + str(hex(crc)) + " != " + str(hex(calced_crc))
    return payload


def reed_solomon_encode(txt: Union[bytes, str, bytearray], number_repair_symbols: int = 2) -> bytes:
    global rscodec, number_r_symbols
    if rscodec is None or number_r_symbols != number_repair_symbols:
        rscodec = RSCodec(number_repair_symbols)
        number_r_symbols = number_repair_symbols
    return bytes(rscodec.encode(txt))


def reed_solomon_decode(txt: Union[bytes, str, bytearray], number_repair_symbols: int = 2) -> bytes:
    global rscodec, number_r_symbols
    if rscodec is None or number_r_symbols != number_repair_symbols:
        rscodec = RSCodec(number_repair_symbols)
        number_r_symbols = number_repair_symbols
    decoded = rscodec.decode(txt)
    # reedsolo returns tuple of arrays, we need bytes
    if isinstance(decoded, tuple):
        return bytes(decoded[0])
    return bytes(decoded)


def dna_reed_solomon_encode(
    txt: bytes,
    number_repair_symbols: int = 2,
    c_exp: int = 2,
    prim_poly: int = 7,
) -> bytes:
    """
    Warning: dna_reed_solomon_* requires a custom RSCodec version with support for custom c_exp and nsize.
    """
    tmp_rscodec = RSCodec(number_repair_symbols, c_exp=c_exp, prim=prim_poly, nsize=2**c_exp - 1)
    dec_list = bits_to_dec(txt)
    return bytes(tmp_rscodec.encode(dec_list))


def dna_reed_solomon_decode(
    txt: bytes,
    number_repair_symbols: int = 2,
    c_exp: int = 2,
    prim_poly: int = 7,
) -> bytes:
    """
    Warning: dna_reed_solomon_* requires a custom RSCodec version with support for custom c_exp and nsize.
    """
    tmp_rscodec = RSCodec(number_repair_symbols, c_exp=c_exp, prim=prim_poly, nsize=2**c_exp - 1)
    decoded = tmp_rscodec.decode(txt)
    # reedsolo returns tuple of arrays, we need bytes
    if isinstance(decoded, tuple):
        return bytes(decoded[0])
    return bytes(decoded)


# Helper functions.
def bits_to_bytes(decoded_bytes: Union[bytearray, List[int], bytes]) -> bytes:
    decoded_bits = dec_to_bits(
        list(decoded_bytes) if isinstance(decoded_bytes, (bytearray, bytes)) else decoded_bytes
    )
    return bitarray(decoded_bits).tobytes()


def str_to_bits(input_string: str) -> str:
    return "".join(format(x, "08b") for x in input_string)


def bytes_to_bits(input_bytes: bytes) -> str:
    bits = bitarray(endian="big")
    bits.frombytes(input_bytes)
    return bits.to01()


def bits_to_dec(input_string: bytes) -> List[int]:
    translation = {"00": 0, "01": 1, "10": 2, "11": 3}
    input_bits = bytes_to_bits(input_string)
    if len(input_bits) % 2 == 1:
        print(input_string)
        print(input_bits)
    return [translation[input_bits[i : i + 2]] for i in range(0, len(input_bits), 2)]


def dec_to_bits(decoded_bytes: Union[List[int], bytearray, bytes]) -> str:
    translation = {0: "00", 1: "01", 2: "10", 3: "11"}
    byte_list = (
        list(decoded_bytes) if isinstance(decoded_bytes, (bytearray, bytes)) else decoded_bytes
    )
    return "".join([translation[bits] for bits in byte_list])


def get_error_correction_decode(e_correction: str, repair_symbols: int) -> Callable[[bytes], bytes]:
    def reed_solomon_decode_with_symbols(data: bytes) -> bytes:
        return reed_solomon_decode(data, repair_symbols)

    def dna_reed_solomon_decode_with_symbols(data: bytes) -> bytes:
        return dna_reed_solomon_decode(data, repair_symbols)

    if e_correction == "nocode":
        error_correction: Callable[[bytes], bytes] = nocode
    elif e_correction == "crc":
        error_correction = crc32_decode
    elif e_correction == "reedsolomon":
        if repair_symbols != 2:
            error_correction = reed_solomon_decode_with_symbols
        else:
            error_correction = reed_solomon_decode
    elif e_correction == "dna_reedsolomon":
        if repair_symbols != 2:
            error_correction = dna_reed_solomon_decode_with_symbols
        else:
            error_correction = dna_reed_solomon_decode
    else:
        raise NotImplementedError(
            "Selected Error Correction not supported, choose: 'nocode', 'crc', 'reedsolomon' or 'dna_reedsolomon'"
        )
    return error_correction


def get_error_correction_encode(e_correction: str, repair_symbols: int) -> Callable[[bytes], bytes]:
    def reed_solomon_encode_with_symbols(data: bytes) -> bytes:
        return reed_solomon_encode(data, repair_symbols)

    def dna_reed_solomon_encode_with_symbols(data: bytes) -> bytes:
        return dna_reed_solomon_encode(data, repair_symbols)

    if e_correction == "nocode":
        error_correction: Callable[[bytes], bytes] = nocode
    elif e_correction == "crc":
        error_correction = crc32
    elif e_correction == "reedsolomon":
        if repair_symbols != 2:
            error_correction = reed_solomon_encode_with_symbols
        else:
            error_correction = reed_solomon_encode
    elif e_correction == "dna_reedsolomon":
        if repair_symbols != 2:
            error_correction = dna_reed_solomon_encode_with_symbols
        else:
            error_correction = dna_reed_solomon_encode
    else:
        raise NotImplementedError(
            "Selected Error Correction not supported, choose: 'nocode', 'crc', 'reedsolomon' or 'dna_reedsolomon'"
        )
    return error_correction


_ERROR_CORRECTION_NAMES = {
    "nocode": "nocode",
    "crc32": "crc",
    "reed_solomon_encode": "reedsolomon",
    "dna_reed_solomon_encode": "dna_reedsolomon",
}


def _get_callable_name(error_correction_func: Callable[..., Any]) -> Optional[str]:
    code = getattr(error_correction_func, "__code__", None)
    if code is None:
        return None
    return code.co_name


def _get_lambda_closure_name(error_correction_func: Callable[..., Any]) -> Optional[str]:
    if _get_callable_name(error_correction_func) != "<lambda>":
        return None
    for cell in getattr(error_correction_func, "__closure__", ()) or ():
        try:
            cell_func = cell.cell_contents
            cell_name = _get_callable_name(cell_func)
            if cell_name in _ERROR_CORRECTION_NAMES:
                return cell_name
        except (ValueError, AttributeError):
            continue
    return None


def get_error_correction_name(error_correction_func: Callable[..., Any]) -> str:
    """
    Get the error correction name string from an error correction function.

    Args:
        error_correction_func: The error correction function (e.g., nocode, crc32, reed_solomon_encode)

    Returns:
        The error correction name string ('nocode', 'crc', 'reedsolomon', or 'dna_reedsolomon')

    Example:
        >>> get_error_correction_name(nocode)
        'nocode'
        >>> get_error_correction_name(crc32)
        'crc'
    """
    func_name = _get_callable_name(error_correction_func)
    if func_name in _ERROR_CORRECTION_NAMES:
        return _ERROR_CORRECTION_NAMES[func_name]

    closure_name = _get_lambda_closure_name(error_correction_func)
    if closure_name in _ERROR_CORRECTION_NAMES:
        return _ERROR_CORRECTION_NAMES[closure_name]

    return "nocode"
