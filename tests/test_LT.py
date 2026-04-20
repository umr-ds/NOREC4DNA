#!/usr/bin/python
# -*- coding: latin-1 -*-
import filecmp
import os
import shutil
from pathlib import Path

import pytest
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import (
    ErlichZielinskiRobustSolitonDistribution,
)
from norec4dna.distributions.IdealSolitonDistribution import IdealSolitonDistribution
from norec4dna.distributions.RobustSolitonDistribution import RobustSolitonDistribution
from norec4dna.Encoder import Encoder
from norec4dna.ErrorCorrection import reed_solomon_decode, reed_solomon_encode
from norec4dna.LTBPDecoder import LTBPDecoder
from norec4dna.LTDecoder import LTDecoder
from norec4dna.LTEncoder import LTEncoder
from norec4dna.rules.DNARules_ErlichZielinski import DNARules_ErlichZielinski
from norec4dna.rules.FastDNARules import FastDNARules

# Get the directory containing this test file
TEST_DIR = Path(__file__).parent.absolute()
# Get the current working directory - decoder saves files here
CWD = Path.cwd()

file = str(TEST_DIR / "logo.jpg")
out_dir = str(TEST_DIR / "LT_logo.jpg")
cmp_file = str(TEST_DIR / "cmp_logo.jpg")

file2 = str(TEST_DIR / "Dorn")
out_dir2 = str(TEST_DIR / "LT_Dorn")
cmp_file2 = str(TEST_DIR / "cmp_dorn")


@pytest.fixture(autouse=True)
def run_between_tests():
    for _file in [file, file2]:
        if os.path.exists(_file):
            os.remove(_file)
        shutil.copy(cmp_file, _file)


@pytest.mark.parametrize("as_dna", [False, True])
@pytest.mark.parametrize("decoder_instance", [LTDecoder, LTBPDecoder])
@pytest.mark.parametrize("distribution", ["robust", "ideal", "ErlichZielinski"])
@pytest.mark.parametrize("use_header", [True, False])
@pytest.mark.parametrize("implicit_mode", [True, False])
def test_suite(as_dna, decoder_instance, distribution, use_header, implicit_mode):
    try:
        os.remove(file)
    except:
        print("Not deleting, File did not exists")
    shutil.copyfile(cmp_file, file)
    chunksize = 200
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks)
    if distribution == "robust":
        dist = RobustSolitonDistribution(S=number_of_chunks, delta=0.2, seed=2)
    elif distribution == "ideal":
        dist = IdealSolitonDistribution(S=number_of_chunks, seed=2)
    else:
        dist = ErlichZielinskiRobustSolitonDistribution(k=number_of_chunks, delta=0.2, seed=2)
    rules = FastDNARules() if as_dna else None
    encoder = LTEncoder(
        file,
        number_of_chunks,
        dist,
        chunk_size=chunksize,
        pseudo_decoder=pseudo_decoder,
        rules=rules,
        insert_header=use_header,
        number_of_chunks_len_format="H",
        id_len_format="H",
        used_packets_len_format="H",
        implicit_mode=implicit_mode,
    )
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
        pseudo_decoder.is_decoded() and pseudo_decoder.getSolvedCount() == encoder.number_of_chunks
    )
    assert os.path.exists(out_dir)
    decoder = decoder_instance(
        out_dir, use_headerchunk=use_header, dist=dist, implicit_mode=implicit_mode
    )
    decoder.decodeFolder(
        number_of_chunks_len_format="H", seed_len_format="H", degree_len_format="H"
    )
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file)
    decoder.saveDecodedFile(print_to_output=False)
    # Decoder saves files:
    # - When use_header=False: Creates DEC_ + basename of output folder in CWD
    # - When use_header=True: Reads full path from header chunk and saves there
    if not use_header:
        out_file = str(CWD / "DEC_LT_logo.jpg")
    else:
        # Header contains full path, decoder saves to that exact path
        out_file = str(TEST_DIR / "logo.jpg")
    assert os.path.exists(out_file) and filecmp.cmp(out_file, cmp_file)
    if decoder_instance == LTBPDecoder:
        # since ApproxDecoder defines an upper bound Gauss-Decoder MUST be able to decode!
        decoder = LTDecoder(
            out_dir, use_headerchunk=use_header, dist=dist, implicit_mode=implicit_mode
        )
        decoder.decodeFolder(
            number_of_chunks_len_format="H", seed_len_format="H", degree_len_format="H"
        )
        assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
        os.remove(out_file)
        decoder.saveDecodedFile(print_to_output=False)
        assert os.path.exists(out_file) and filecmp.cmp(out_file, cmp_file)
    shutil.rmtree(out_dir)


def test_erlich_zielinski_dnarules():
    try:
        os.remove(str(f"{TEST_DIR}/DEC_{os.path.basename(out_dir2)}"))
    except:
        print("Not deleting, File did not exists")
    shutil.copyfile(cmp_file2, file2)
    chunksize = 75
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file2, chunksize)
    pseudo_decoder = LTDecoder.pseudo_decoder(number_of_chunks)
    dist = ErlichZielinskiRobustSolitonDistribution(k=number_of_chunks, delta=0.2, seed=2)
    rules = DNARules_ErlichZielinski()
    encoder = LTEncoder(
        file2,
        number_of_chunks,
        dist,
        error_correction=reed_solomon_encode,
        chunk_size=chunksize,
        rules=rules,
        insert_header=False,
        number_of_chunks_len_format="H",
        id_len_format="H",
        used_packets_len_format="H",
        pseudo_decoder=pseudo_decoder,
        implicit_mode=True,
        save_number_of_chunks_in_packet=False,
    )
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=True)
    assert (
        pseudo_decoder.is_decoded() and pseudo_decoder.getSolvedCount() == encoder.number_of_chunks
    )
    assert os.path.exists(out_dir2)
    decoder = LTDecoder(
        out_dir2,
        use_headerchunk=False,
        dist=dist,
        implicit_mode=True,
        static_number_of_chunks=encoder.number_of_chunks,
        error_correction=reed_solomon_decode,
    )
    decoder.decodeFolder(
        number_of_chunks_len_format="H", seed_len_format="H", degree_len_format="H"
    )
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    decoder.saveDecodedFile(print_to_output=False, null_is_terminator=True)
    # Decoder saves files to current working directory (where pytest runs from)
    out_file2 = str(CWD / "DEC_LT_Dorn")
    assert os.path.exists(out_file2) and filecmp.cmp(out_file2, cmp_file2)
    os.remove(out_file2)


def test_size_shrink():
    try:
        os.remove(file)
    except:
        print("Not deleting, File did not exists")
    shutil.copyfile(cmp_file, file)
    chunksize = 200
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    dist = RobustSolitonDistribution(S=number_of_chunks, delta=0.2, seed=2)
    encoder1 = LTEncoder(
        file,
        number_of_chunks,
        dist,
        chunk_size=chunksize,
        rules=None,
        insert_header=False,
        number_of_chunks_len_format="H",
        id_len_format="H",
        used_packets_len_format="H",
        implicit_mode=False,
    )
    number_of_chunks = encoder1.number_of_chunks
    encoder1.prepareEncoder()
    packet1 = encoder1.create_new_packet()

    encoder2 = LTEncoder(
        file,
        number_of_chunks,
        dist,
        chunk_size=chunksize,
        rules=None,
        insert_header=False,
        number_of_chunks_len_format="H",
        id_len_format="H",
        used_packets_len_format="H",
        implicit_mode=True,
    )
    encoder2.prepareEncoder()
    packet2 = encoder2.create_new_packet()
    assert len(packet1.get_dna_struct(True)) > len(packet2.get_dna_struct(True))

    encoder3 = LTEncoder(
        file,
        number_of_chunks,
        dist,
        chunk_size=chunksize,
        rules=None,
        insert_header=False,
        number_of_chunks_len_format="I",
        id_len_format="I",
        used_packets_len_format="I",
        implicit_mode=False,
    )
    encoder3.prepareEncoder()
    packet3 = encoder3.create_new_packet()

    encoder4 = LTEncoder(
        file,
        number_of_chunks,
        dist,
        chunk_size=chunksize,
        rules=None,
        insert_header=False,
        number_of_chunks_len_format="I",
        id_len_format="I",
        used_packets_len_format="I",
        implicit_mode=True,
    )
    encoder4.prepareEncoder()
    packet4 = encoder4.create_new_packet()
    assert len(packet3.get_dna_struct(True)) > len(packet4.get_dna_struct(True))
    assert len(packet3.get_dna_struct(True)) > len(packet1.get_dna_struct(True))
    assert len(packet4.get_dna_struct(True)) > len(packet2.get_dna_struct(True))
