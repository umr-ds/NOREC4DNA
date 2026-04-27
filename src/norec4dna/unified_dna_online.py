import io
import math
import os
import typing

import numpy as np
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.Packet import Packet

from . import (
    Encoder,
    OnlineDecoder,
    OnlineDistribution,
    OnlineEncoder,
    get_error_correction_decode,
    get_error_correction_encode,
)
from .helper.quaternary2Bin import tranlate_quat_to_byte
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
    best = None
    for i in range(steps + 1):
        eps = eps_min + (eps_max - eps_min) * (i / steps)
        s = _s_from_eps(eps)
        if s is None:
            continue

        # Basic feasibility: s must be well below K
        if s > s_max:
            continue

        # Budget feasibility: need at least s extras
        extra = math.ceil(eps * K)
        if extra < s:
            continue

        # Degree-1 mass constraint for peeling stability
        p1 = _p1_from_eps_s(eps, s)
        if not (0.0 <= p1 <= 1.0):
            continue
        if p1 < p1_min:
            continue

        # First feasible eps is the minimal one due to monotone search
        best = (eps, extra, K, s, p1)
        break

    if best is None:
        # If nothing feasible was found, fall back to eps_max and report what it yields.
        eps = eps_max
        s = _s_from_eps(eps)
        if s is None:
            raise RuntimeError("Unable to derive s from the provided eps_max; widen search bounds.")
        extra = math.ceil(eps * K)
        p1 = _p1_from_eps_s(eps, s)
        return eps, extra, K, s, p1

    return best


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


def decode(string_file_name, list_of_dna_strings):
    # make sure that the dist is freshly initialized...
    decoder = OnlineDecoder(
        string_file_name,
        error_correction=error_correction_func_dec,
        use_headerchunk=INSERT_HEADER,
        static_number_of_chunks=NUMBER_OF_CHUNKS,
    )
    decoder.read_all_before_decode = READ_ALL

    for dna_str in list_of_dna_strings:
        new_pack = decoder.parse_raw_packet(
            io.BytesIO(tranlate_quat_to_byte(dna_str)).read(),
            crc_len_format=CHECKSUM_LEN_STR,
            number_of_chunks_len_format=DECODER_NUM_CHUNK_LEN_FORMAT,
            check_block_number_len_format=SEED_LEN_FORMAT,
        )
        if new_pack is not None and new_pack != "CORRUPT":
            decoder.input_new_packet(new_pack)

    decoder.solve()
    __byte_io = io.BytesIO()
    with __byte_io as f:
        for x in decoder.GEPP.result_mapping:
            if x < 0:
                f.write(b"\x00" * len(decoder.GEPP.b[x][0]))
                dirty = True
                continue
            if INSERT_HEADER and decoder.headerChunk is None:
                header_row = decoder.GEPP.result_mapping[0]
                decoder.headerChunk = HeaderChunk(
                    Packet(
                        decoder.GEPP.b[header_row], {0}, decoder.number_of_chunks, read_only=True
                    ),
                    checksum_len_format=CHECKSUM_LEN_STR,
                )
            if 0 != x or not INSERT_HEADER:
                if decoder.number_of_chunks - 1 == x and INSERT_HEADER:
                    output = decoder.GEPP.b[x][0][0 : decoder.headerChunk.get_last_chunk_length()]
                    f.write(output)
                else:
                    if NULL_IS_TERMINATOR:
                        splitter: str = decoder.GEPP.b[x].tostring().decode().split("\x00")
                        output = splitter[0].encode()
                        f.write(output)
                        if len(splitter) > 1:
                            break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                    else:
                        output = decoder.GEPP.b[x]
                        f.write(output)
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
        except:
            print(f"{INPUT_FILE}, False")
