import functools
import io

import numpy as np
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import (
    ErlichZielinskiRobustSolitonDistribution,
)
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.Packet import Packet

from . import (
    Encoder,
    LTDecoder,
    LTEncoder,
    get_error_correction_decode,
    get_error_correction_encode,
)
from .helper.quaternary2Bin import tranlate_quat_to_byte
from .rules.DNARules_ErlichZielinski import DNARules_ErlichZielinski

# INPUT_FILE = "data_2mb.test"
OVERHEAD = 0.9
INSERT_HEADER = True
IMPLICIT_MODE = True  # should be left True
NUMBER_OF_CHUNKS_IN_PACKET = False

# either set  NUMBER_OF_CHUNKS or CHUNK_SIZE !
CHUNK_SIZE = 75
# the decoder needs to know the number of chunks: if CHUNK_SIZE was used, enter the number of chunks here
# (you could store the NUMBER OF CHUNKS in each packet header but this would increase the overhead)
# alternatively, one could infer the number of chunks from the number of encoded packets and the expected overhead
# OR one could bruteforce the number of chunks (there are only a few possible values...)
NUMBER_OF_CHUNKS = None

SEED = 2
DIST = ErlichZielinskiRobustSolitonDistribution
dist_func = functools.partial(DIST, seed=SEED)

ERROR_CORRECTION = "reedsolomon"  # packet level error correction
REPAIR_SYMBOLS = 2  # number of symbols / bytes for each packet
error_correction_func = get_error_correction_encode(ERROR_CORRECTION, REPAIR_SYMBOLS)
error_correction_func_dec = get_error_correction_decode(ERROR_CORRECTION, REPAIR_SYMBOLS)

# THESE RULES MUST BE ADAPTED TO THE ERROR CHANNEL! "FastDNARules()" should be used and modified accordingly!
DNA_RULES = DNARules_ErlichZielinski()

DROP_UPPER_BOUND = 0.9

NUMBER_OF_CHUNKS_LEN_STR = "H"
ID_LEN_STR = "I"
USED_PACKETS_LEN_STR = "H"
CHECKSUM_LEN_STR = "H"

DEGREE_LEN_STR = "I"

# instead of the Gaussian elimination based decoder, the belief-propagation decoder could be used to speed things up.
# however: this might require a higher overhead and is thus discouraged for dna data storage!
# additionally we would have to adapt the decode-logic used in _this_ script
DECODER_CLASS = LTDecoder

READ_ALL = True
NULL_IS_TERMINATOR = False  # should only be set
NUM_CHUNK_LEN_FORMAT = ""  # should only be set if the number of chunks is stored in each packet

RAISE_ON_UNSOLVED = True  # if set to false, the decoder will return partial results if a full decode is not possible.


def encode(string_file_name, numpy_boolean_array):
    global NUMBER_OF_CHUNKS
    # we operate on files, not raw bits, thus we only use the file name and load the data from disk,
    # if required we could change this...
    if NUMBER_OF_CHUNKS is None:
        NUMBER_OF_CHUNKS = Encoder.get_number_of_chunks_for_file_with_chunk_size(
            string_file_name, chunk_size=CHUNK_SIZE, insert_header=INSERT_HEADER
        )
    encoder = LTEncoder(
        string_file_name,
        NUMBER_OF_CHUNKS,
        dist_func(NUMBER_OF_CHUNKS),
        insert_header=INSERT_HEADER,
        rules=DNA_RULES,
        error_correction=error_correction_func,
        number_of_chunks_len_format=NUM_CHUNK_LEN_FORMAT,
        id_len_format=ID_LEN_STR,
        used_packets_len_format="H",
        save_number_of_chunks_in_packet=NUMBER_OF_CHUNKS_IN_PACKET,
        implicit_mode=IMPLICIT_MODE,
        drop_upper_bound=DROP_UPPER_BOUND,
        pseudo_decoder=None,
    )  # ,
    # checksum_len_str=CHECKSUM_LEN_STR)
    encoder.set_overhead_limit(OVERHEAD)
    encoder.encode_to_packets()
    return [x.get_dna_struct(True) for x in encoder.encodedPackets], encoder


def decode(string_file_name, list_of_dna_strings):
    # make sure that the dist is freshly initialized...
    decoder = DECODER_CLASS(
        string_file_name,
        error_correction=error_correction_func_dec,
        use_headerchunk=INSERT_HEADER,
        static_number_of_chunks=NUMBER_OF_CHUNKS,
        implicit_mode=IMPLICIT_MODE,
        dist=dist_func(NUMBER_OF_CHUNKS),
    )
    decoder.read_all_before_decode = READ_ALL

    for dna_str in list_of_dna_strings:
        new_pack = decoder.parse_raw_packet(
            io.BytesIO(tranlate_quat_to_byte(dna_str)).read(),
            crc_len_format=CHECKSUM_LEN_STR,
            number_of_chunks_len_format=NUM_CHUNK_LEN_FORMAT,
            degree_len_format=DEGREE_LEN_STR,
            seed_len_format=ID_LEN_STR,
        )
        if new_pack is not None and new_pack != "CORRUPT":
            decoder.input_new_packet(new_pack)

    solved = decoder.solve()
    if RAISE_ON_UNSOLVED and not solved:
        raise RuntimeError(
            "Could not solve the system of equations. A partial recovery might be possible!"
        )
    if not solved:
        print("Could not solve the system of equations. A partial recovery will be performed:")
    __byte_io = io.BytesIO()
    with __byte_io as f:
        for x in decoder.GEPP.result_mapping:
            if x < 0:
                f.write(b"\x00" * len(decoder.GEPP.b[x][0]))
                continue
            if INSERT_HEADER and decoder.headerChunk is None:
                decoder.headerChunk = HeaderChunk(
                    Packet(decoder.GEPP.b[0], {0}, decoder.number_of_chunks, read_only=True)
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
        res, encoder = encode(INPUT_FILE, None)
        try:
            with open(INPUT_FILE, "rb") as f:
                org = np.unpackbits(np.frombuffer(f.read(), dtype=np.uint8))
                decoded = decode(None, res)
                print(f"{INPUT_FILE}, {np.all(np.equal(org, decoded))}")
        except Exception as ex:
            print(ex)
            print(f"{INPUT_FILE}, False")
