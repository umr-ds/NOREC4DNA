import filecmp
import os
import shutil
from math import ceil

import pytest
from norec4dna.helper.bin2Quaternary import string2QUATS
from norec4dna.helper.quaternary2Bin import quats_to_bytes
from norec4dna.ReedSolomonSuite import ReedSolomonDecoder, ReedSolomonEncoder, get_file_size

file = "logo.jpg"
file_dec = "RS_logo.jpg.DECODED"
out_dir = "RS_logo.jpg"


def _flip_base(base: str) -> str:
    transitions = {"A": "T", "T": "G", "G": "C", "C": "A"}
    return transitions.get(base, base)


def _bytes_from_quats(dna_data: str) -> bytes:
    dna_data_temp = b""
    for index in range(0, len(dna_data), 4):
        try:
            dna_data_temp += quats_to_bytes(dna_data[index : index + 4])
        except ValueError:
            continue
    return dna_data_temp


def _mutate_packet_file(path: str, flip_bases: int) -> None:
    with open(path, "rb+") as f:
        dna_data = list("".join(string2QUATS(f.read())))
        if not dna_data:
            return
        for index in range(16, 16 + flip_bases):
            dna_data[index] = _flip_base(dna_data[index])
        f.seek(0)
        f.write(_bytes_from_quats("".join(dna_data)))
        f.truncate()


# testing the ReedSolomonSuite with different parameters
@pytest.mark.parametrize("overhead", [0.10])
@pytest.mark.parametrize("chunksize", [100])
@pytest.mark.parametrize("headerchunk", [True])
@pytest.mark.parametrize("flip_bases", [1])
def test_suite(overhead, chunksize, headerchunk, flip_bases):
    file_size = get_file_size(file)
    number_of_chunks = ceil(1.0 * file_size / chunksize) + 1 if headerchunk else 0
    encoder = ReedSolomonEncoder(file, number_of_chunks, 0, overhead)
    encoder.save_packets(True)
    # manipulating some packets
    for fi in os.listdir(out_dir):
        _mutate_packet_file(out_dir + "/" + fi, flip_bases)
    decoder = ReedSolomonDecoder(out_dir)
    decoder.decodeFolder()
    assert os.path.exists(file_dec)
    assert filecmp.cmp(file_dec, file)
    os.remove(file_dec)
    shutil.rmtree(out_dir)
