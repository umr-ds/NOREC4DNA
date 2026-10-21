import pytest

from norec4dna.Encoder import Encoder
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.Packet import Packet


@pytest.mark.parametrize("last_chunk_format", ["B", "H", "I"])
@pytest.mark.parametrize(
    "checksum_format, checksum",
    [(None, None), ("", None), ("B", 0), ("H", 1), ("I", 1),
     ("I", 0x12003456), ("I", 0x12345678), ("Q", 1)],
)
@pytest.mark.parametrize("as_array", [False, True])
def test_header_filename_round_trip(last_chunk_format, checksum_format, checksum, as_array):
    encoder = Encoder("archive.bin", 2, None, chunk_size=64)
    encoder.chunks = [b"payload"]
    data = encoder.encode_header_info(
        checksum=checksum, checksum_len_str=checksum_format,
        last_chunk_len_format=last_chunk_format,
    )
    if not as_array:
        data = data.tobytes()
    packet = Packet(data, {0}, 2, read_only=True)

    header = HeaderChunk(packet, last_chunk_format, checksum_format)

    assert header.get_file_name() == b"archive.bin"
    assert header.get_last_chunk_length() == len(b"payload")
    assert header.checksum == checksum
