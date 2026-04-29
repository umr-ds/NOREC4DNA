import io
import math
import os
import typing

import numpy as np

from .distributions.OnlineDistribution import OnlineDistribution
from .Encoder import Encoder
from .ErrorCorrection import get_error_correction_decode, get_error_correction_encode
from .HeaderChunk import HeaderChunk
from .helper.quaternary2Bin import tranlate_quat_to_byte
from .OnlineDecoder import OnlineDecoder
from .OnlineEncoder import OnlineEncoder
from .Packet import Packet
from .rules.FastDNARules import FastDNARules

# INPUT_FILE = "Dorn"
OVERHEAD = 0.2
INSERT_HEADER = True
DROP_UPPER_BOUND = (
    1.0  # decreasing this value will drop more packets but ensure all rules are followed
)
NUMBER_OF_CHUNKS_IN_PACKET = False

# either set  NUMBER_OF_CHUNKS or CHUNK_SIZE !
CHUNK_SIZE = 75
# the decoder needs to know the number of chunks: if CHUNK_SIZE was used, enter the number of chunks here
# (you could store the NUMBER OF CHUNKS in each packet header but this would increase the overhead)
# alternatively, one could infer the number of chunks from the number of encoded packets and the expected overhead
# OR one could brtueforce the number of chunks (there are only a few possible values...)
NUMBER_OF_CHUNKS = None

# if the header is too small for the data, consider renaming the file or decreasing
# the size of the checksum (the last chunk len field can not be changed in this tool...):
CHECKSUM_LEN_STR = "H"

# For Online encoding, one should set the epsilon and quality parameters and derive the number of chunks from the file size..
# EPS = 0.19
QUALITY = 8
SEED = 2

ERROR_CORRECTION = "reedsolomon"  # packet level error correction
REPAIR_SYMBOLS = 2  # number of symbols / bytes for each packet
error_correction_func = get_error_correction_encode(ERROR_CORRECTION, REPAIR_SYMBOLS)
error_correction_func_dec = get_error_correction_decode(ERROR_CORRECTION, REPAIR_SYMBOLS)

DIST = OnlineDistribution  # in general this should be left unchanged

DNA_RULES = FastDNARules()

READ_ALL = True
NULL_IS_TERMINATOR = False  # should only be set
DECODER_NUM_CHUNK_LEN_FORMAT = (
    ""  # should only be set if the number of chunks is stored in each packet
)
SEED_LEN_FORMAT = "I"


def _s_from_eps(eps: float) -> typing.Optional[int]:
    if not (0.0 < eps < 1.0):
        return None
    num = math.log((eps * eps) / 4.0)
    den = math.log(1.0 - eps / 2.0)
    if den >= 0.0:  # should be negative for 0<eps<1
        return None
    s = math.ceil(num / den)
    return s if s > 0 else None


def _p1_from_eps_s(eps: float, s: int) -> float:
    # p1 = 1 - (1 + 1/s) / (1 + eps)
    return 1.0 - (1.0 + 1.0 / s) / (1.0 + eps)


def _evaluate_online_epsilon(
    eps: float, k_value: int, s_max: int, p1_min: float
) -> typing.Optional[typing.Tuple[float, int, int, int, float]]:
    s_value = _s_from_eps(eps)
    if s_value is None or s_value > s_max:
        return None

    extra = math.ceil(eps * k_value)
    if extra < s_value:
        return None

    p1 = _p1_from_eps_s(eps, s_value)
    if not (0.0 <= p1 <= 1.0) or p1 < p1_min:
        return None
    return eps, extra, k_value, s_value, p1


def best_epsilon_online(
    file_size_bytes: int,
    droplet_len_bytes: int,
    # tuning knobs
    s_max_frac: float = 0.05,  # require s <= s_max_frac * K
    p1_min: float = 0.2,  # require degree-1 mass >= p1_min
    eps_min: float = 0.02,  # search lower bound
    eps_max: float = 0.50,  # search upper bound
    steps: int = 2000,  # grid resolution for the search
) -> typing.Tuple[float, int, int, int, float]:
    """
    Find the smallest epsilon in [eps_min, eps_max] such that:
      - s = ceil(log(eps^2/4)/log(1 - eps/2)) is defined and positive
      - s <= s_max_frac * K
      - epsilon*K >= s  (enough budget for the late checks)
      - p1 >= p1_min, where p1 = 1 - (1 + 1/s)/(1 + eps)

    Returns:
      (epsilon, extra_symbols, K, s, p1)
    """
    if droplet_len_bytes <= 0 or file_size_bytes <= 0:
        raise ValueError("file_size_bytes and droplet_len_bytes must be positive")

    K = math.ceil(file_size_bytes / droplet_len_bytes)
    if K <= 0:
        raise ValueError("Derived K must be positive")

    s_max = max(1, int(math.floor(s_max_frac * K)))

    # Monotone grid from small to larger eps (s decreases with increasing eps).
    for i in range(steps + 1):
        eps = eps_min + (eps_max - eps_min) * (i / steps)
        best = _evaluate_online_epsilon(eps, K, s_max, p1_min)
        if best is not None:
            return best

    s_value = _s_from_eps(eps_max)
    if s_value is None:
        raise RuntimeError("Unable to derive s from the provided eps_max; widen search bounds.")
    extra = math.ceil(eps_max * K)
    p1 = _p1_from_eps_s(eps_max, s_value)
    return eps_max, extra, K, s_value, p1


def encode(string_file_name):
    global NUMBER_OF_CHUNKS
    # we operate on files, not raw bits, thus we only use the file name and load the data from disk,
    # if required we could change this...
    if NUMBER_OF_CHUNKS is None:
        NUMBER_OF_CHUNKS = Encoder.get_number_of_chunks_for_file_with_chunk_size(
            string_file_name, chunk_size=CHUNK_SIZE, insert_header=INSERT_HEADER
        )
    encoder = OnlineEncoder(
        string_file_name,
        NUMBER_OF_CHUNKS,
        DIST(eps=EPS),
        insert_header=INSERT_HEADER,
        epsilon=EPS,
        quality=QUALITY,
        rules=DNA_RULES,
        error_correction=error_correction_func,
        number_of_chunks_len_format=DECODER_NUM_CHUNK_LEN_FORMAT,
        check_block_number_len_format=SEED_LEN_FORMAT,
        save_number_of_chunks_in_packet=NUMBER_OF_CHUNKS_IN_PACKET,
        drop_upper_bound=DROP_UPPER_BOUND,
        checksum_len_str=CHECKSUM_LEN_STR,
    )
    encoder.set_overhead_limit(OVERHEAD)
    encoder.encode_to_packets()
    return [x.get_dna_struct(True) for x in encoder.encodedPackets], encoder


def _iter_online_packets(decoder: OnlineDecoder, list_of_dna_strings):
    for dna_str in list_of_dna_strings:
        new_pack = decoder.parse_raw_packet(
            io.BytesIO(tranlate_quat_to_byte(dna_str)).read(),
            crc_len_format=CHECKSUM_LEN_STR,
            number_of_chunks_len_format=DECODER_NUM_CHUNK_LEN_FORMAT,
            check_block_number_len_format=SEED_LEN_FORMAT,
        )
        if new_pack is not None and new_pack != "CORRUPT":
            yield new_pack


def _ensure_online_header_chunk(decoder: OnlineDecoder, gepp, row_index: int) -> None:
    if not INSERT_HEADER or decoder.headerChunk is not None:
        return
    header_row = gepp.result_mapping[0]
    decoder.headerChunk = HeaderChunk(
        Packet(gepp.b[header_row], {0}, decoder.number_of_chunks, read_only=True),
        checksum_len_format=CHECKSUM_LEN_STR,
    )


def _write_online_row(
    decoder: OnlineDecoder, gepp, row_index: int, output_buffer: io.BytesIO
) -> bool:
    if row_index < 0:
        output_buffer.write(b"\x00" * len(gepp.b[row_index][0]))
        return False

    _ensure_online_header_chunk(decoder, gepp, row_index)
    if row_index == 0 and INSERT_HEADER:
        return False

    if decoder.number_of_chunks - 1 == row_index and INSERT_HEADER:
        header_chunk = decoder.headerChunk
        if header_chunk is None:
            raise RuntimeError("Header chunk missing during Online decode")
        output_buffer.write(gepp.b[row_index][0][0 : header_chunk.get_last_chunk_length()])
        return False

    if NULL_IS_TERMINATOR:
        splitter = gepp.b[row_index].tobytes().decode().split("\x00")
        output_buffer.write(splitter[0].encode())
        return len(splitter) > 1

    output_buffer.write(gepp.b[row_index].tobytes())
    return False


def decode(string_file_name, list_of_dna_strings):
    # make sure that the dist is freshly initialized...
    assert NUMBER_OF_CHUNKS is not None
    decoder = OnlineDecoder(
        string_file_name,
        error_correction=error_correction_func_dec,
        use_headerchunk=INSERT_HEADER,
        static_number_of_chunks=NUMBER_OF_CHUNKS,
    )
    decoder.read_all_before_decode = READ_ALL

    for new_pack in _iter_online_packets(decoder, list_of_dna_strings):
        decoder.input_new_packet(new_pack)

    decoder.solve()
    gepp = decoder.GEPP
    if gepp is None:
        raise RuntimeError("Decoder GEPP was not initialized")
    __byte_io = io.BytesIO()
    with __byte_io as f:
        for x in gepp.result_mapping:
            if _write_online_row(decoder, gepp, x, f):
                break
        # convert the byte array __ByteIO to a numpy bool array
        numpy_boolean_array = np.unpackbits(np.frombuffer(f.getvalue(), dtype=np.uint8))
        return numpy_boolean_array


if __name__ == "__main__":
    for INPUT_FILE in [
        ".INFILES/Dorn",
        ".INFILES/sleeping_beauty",
        "README.md",
        ".INFILES/logo.jpg",
        ".INFILES/data_1mb.test",
        ".INFILES/data_2mb.test",
    ]:
        NUMBER_OF_CHUNKS = None
        file_size = os.path.getsize(INPUT_FILE)
        EPS, extra, K, s, p1 = best_epsilon_online(file_size, CHUNK_SIZE)
        res, encoder = encode(INPUT_FILE)
        try:
            with open(INPUT_FILE, "rb") as f:
                org = np.unpackbits(np.frombuffer(f.read(), dtype=np.uint8))
                decoded = decode(None, res)
                print(f"{INPUT_FILE}, {np.all(np.equal(org, decoded))}")
        except Exception:
            print(f"{INPUT_FILE}, False")
