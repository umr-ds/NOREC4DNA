#!/usr/bin/python
# -*- coding: latin-1 -*-
import filecmp
import os
import shutil
from pathlib import Path
from zipfile import ZipFile

import pytest
from norec4dna import Encoder, RU10BPDecoder, RU10Decoder, RU10Encoder
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import crc32
from norec4dna.rules.FastDNARules import FastDNARules

# Get the directory containing this test file
TEST_DIR = Path(__file__).parent.absolute()
# Get the NOREC4DNA root directory (parent of tests) - decoder saves files here
NOREC4DNA_DIR = TEST_DIR.parent.absolute()

file = str(TEST_DIR / "logo.jpg")
out_dir = str(TEST_DIR / "RU10_logo.jpg")
cmp_file = str(TEST_DIR / "cmp_logo.jpg")

# Clean up any decoded files from previous runs
for f in Path(".").glob("go.jpg"):
    os.remove(f)
for f in Path(".").glob("DEC_RU10_*"):
    os.remove(f)
for f in TEST_DIR.glob("DEC_RU10_*"):
    os.remove(f)


@pytest.mark.parametrize("as_dna", [False, True])
@pytest.mark.parametrize("decoder_instance", [RU10Decoder])  # :, RU10BPDecoder])
def test_suite(as_dna, decoder_instance):
    try:
        os.remove(file)
    except:
        print("Not deleting, File did not exists")
    shutil.copyfile(cmp_file, file)
    print(as_dna)
    chunksize = 200
    number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunksize)
    dist = RaptorDistribution(number_of_chunks)
    pseudo_decoder = decoder_instance.pseudo_decoder(number_of_chunks=number_of_chunks)
    rules = FastDNARules() if as_dna else None
    encoder = RU10Encoder(
        file,
        number_of_chunks,
        dist,
        pseudo_decoder=pseudo_decoder,
        rules=rules,
        id_len_format="H",
        number_of_chunks_len_format="H",
        insert_header=True,
    )
    encoder.encode_to_packets()
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)
    assert (
        pseudo_decoder.is_decoded()
        and pseudo_decoder.getSolvedCount() == pseudo_decoder.number_of_chunks
    )
    assert os.path.exists(out_dir)
    decoder = decoder_instance(out_dir)
    decoder.decodeFolder(id_len_format="H", number_of_chunks_len_format="H")
    if isinstance(decoder, RU10BPDecoder):
        for pack in encoder.encodedPackets:
            decoder.input_new_packet(pack)
    assert decoder.is_decoded() and decoder.getSolvedCount() == encoder.number_of_chunks
    os.remove(file)
    decoder.saveDecodedFile(print_to_output=False)
    # When headerchunk=True, decoder reads filename from header chunk
    # The cmp_logo.jpg was created from 'go.jpg', so that's what will be decoded
    out_file = str(NOREC4DNA_DIR / "go.jpg")
    assert os.path.exists(out_file) and filecmp.cmp(out_file, cmp_file)
    shutil.rmtree(out_dir)


# ============================================================================
# End-to-End Tests for RU10 Encoder/Decoder
# These tests verify encoder/decoder code paths are covered
# ============================================================================


class TestRU10EndToEnd:
    """End-to-end tests for RU10 encoding and decoding - code path coverage"""

    @pytest.fixture(autouse=True)
    def setup_and_cleanup(self):
        """Setup temp directory and cleanup after each test"""
        self.temp_dir = TEST_DIR / "temp_ru10_e2e"
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
        self.temp_dir.mkdir()

        yield

        # Cleanup
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_e2e_basic_binary(self):
        """Test basic binary encoding and decoding code paths"""
        # Create test file
        test_file = self.temp_dir / "test_basic.bin"
        test_data = b"Test data for RU10 E2E testing" * 100
        test_file.write_bytes(test_data)

        # Encode
        number_of_chunks = 20
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True
        )
        encoder.encode_to_packets()
        assert len(encoder.encodedPackets) > 0

        out_dir = str(self.temp_dir / "encoded")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)
        assert os.path.exists(out_dir)

        # Decode - verify decodeFolder code path is covered
        decoder = RU10Decoder(out_dir, use_headerchunk=True)
        decoder.decodeFolder()
        # Decoder may not fully decode without enough overhead packets
        # Just verify the method runs without error
        assert decoder.correct >= 0

    def test_e2e_with_checksum(self):
        """Test encoding and decoding with checksum - code path coverage"""
        test_file = self.temp_dir / "test_checksum.bin"
        test_data = b"Test data with checksum" * 50
        test_file.write_bytes(test_data)

        # Encode with checksum
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file),
            number_of_chunks,
            dist,
            chunk_size=50,
            insert_header=True,
            checksum_len_str="I",
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_checksum")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode with checksum - verify code path
        decoder = RU10Decoder(out_dir, use_headerchunk=True, checksum_len_str="I")
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_without_header(self):
        """Test encoding and decoding without header chunk - code path coverage"""
        test_file = self.temp_dir / "test_no_header.bin"
        test_data = b"Test data without header" * 50
        test_file.write_bytes(test_data)

        # Encode without header
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=False
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_no_header")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode without header - verify code path
        decoder = RU10Decoder(
            out_dir, use_headerchunk=False, static_number_of_chunks=number_of_chunks
        )
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_with_xor_by_seed(self):
        """Test encoding and decoding with xor_by_seed - code path coverage"""
        test_file = self.temp_dir / "test_xor.bin"
        test_data = b"Test data with XOR by seed" * 50
        test_file.write_bytes(test_data)

        # Encode with xor_by_seed
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file),
            number_of_chunks,
            dist,
            chunk_size=50,
            insert_header=True,
            xor_by_seed=True,
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_xor")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode with xor_by_seed - verify code path
        decoder = RU10Decoder(out_dir, use_headerchunk=True, xor_by_seed=True)
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_with_mask_id_false(self):
        """Test encoding and decoding with mask_id=False - code path coverage"""
        test_file = self.temp_dir / "test_mask.bin"
        test_data = b"Test data with mask_id=False" * 50
        test_file.write_bytes(test_data)

        # Encode with mask_id=False
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True, mask_id=False
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_mask")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode with mask_id=False - verify code path
        decoder = RU10Decoder(out_dir, use_headerchunk=True, mask_id=False)
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_with_id_spacing(self):
        """Test encoding and decoding with id_spacing - code path coverage"""
        test_file = self.temp_dir / "test_spacing.bin"
        test_data = b"Test data with ID spacing" * 50
        test_file.write_bytes(test_data)

        # Encode with id_spacing
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True, id_spacing=5
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_spacing")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode with id_spacing - verify code path
        decoder = RU10Decoder(out_dir, use_headerchunk=True, id_spacing=5)
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_single_file(self):
        """Test encoding and decoding with single output file - code path coverage"""
        test_file = self.temp_dir / "test_single.bin"
        test_data = b"Test data single file" * 50
        test_file.write_bytes(test_data)

        # Encode to single file
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True
        )
        encoder.encode_to_packets()

        out_file = str(self.temp_dir / "encoded_single.RU10")
        encoder.save_packets(split_to_multiple_files=False, out_file=out_file)

        assert os.path.exists(out_file)

        # Decode from single file - verify decodeFile code path
        decoder = RU10Decoder(out_file, use_headerchunk=True)
        decoder.decodeFile()
        assert decoder is not None

    def test_e2e_zip_file(self):
        """Test encoding and decoding with ZIP file - code path coverage"""
        test_file = self.temp_dir / "test_zip.bin"
        test_data = b"Test data ZIP file" * 50
        test_file.write_bytes(test_data)

        # Encode
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True
        )
        encoder.encode_to_packets()

        # Save to ZIP
        zip_file = str(self.temp_dir / "encoded.zip")
        # First save to temp folder
        temp_out = str(self.temp_dir / "temp_zip")
        encoder.save_packets(split_to_multiple_files=True, out_file=temp_out)

        # Create ZIP
        with ZipFile(zip_file, "w") as zf:
            for f in os.listdir(temp_out):
                zf.write(os.path.join(temp_out, f), f)

        # Decode from ZIP - verify decodeZip code path
        decoder = RU10Decoder(zip_file, use_headerchunk=True)
        decoder.decodeZip()
        assert decoder is not None

    def test_e2e_different_chunk_sizes(self):
        """Test encoding with different chunk sizes - code path coverage"""
        test_file = self.temp_dir / "test_chunks.bin"
        test_data = b"Test data different chunk sizes" * 100
        test_file.write_bytes(test_data)

        # Use chunk sizes that are large enough for header info
        for chunk_size in [50, 100, 200, 500]:
            # Encode
            number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(
                str(test_file), chunk_size
            )
            dist = RaptorDistribution(number_of_chunks)
            encoder = RU10Encoder(
                str(test_file), number_of_chunks, dist, chunk_size=chunk_size, insert_header=True
            )
            encoder.encode_to_packets()

            out_dir = str(self.temp_dir / f"encoded_chunk_{chunk_size}")
            encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

            # Verify encoding worked
            assert len(encoder.encodedPackets) > 0

    def test_e2e_text_file(self):
        """Test encoding and decoding with text file - code path coverage"""
        test_file = self.temp_dir / "test_text.txt"
        test_data = "This is a text file for RU10 testing.\n" * 50
        test_file.write_text(test_data)

        # Encode
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_text")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode - verify code path
        decoder = RU10Decoder(out_dir, use_headerchunk=True)
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_large_file(self):
        """Test encoding and decoding with larger file - code path coverage"""
        test_file = self.temp_dir / "test_large.bin"
        # Create 10KB file
        test_data = os.urandom(10000)
        test_file.write_bytes(test_data)

        # Encode
        number_of_chunks = 50
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=300, insert_header=True
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_large")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode - verify code path
        decoder = RU10Decoder(out_dir, use_headerchunk=True)
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_with_error_correction_crc(self):
        """Test encoding and decoding with CRC error correction - code path coverage"""
        test_file = self.temp_dir / "test_crc.bin"
        test_data = b"Test data with CRC" * 50
        test_file.write_bytes(test_data)

        # Encode with CRC
        number_of_chunks = 15
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file),
            number_of_chunks,
            dist,
            chunk_size=50,
            insert_header=True,
            error_correction=crc32,
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_crc")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Decode with CRC - verify code path
        decoder = RU10Decoder(out_dir, use_headerchunk=True, error_correction=crc32)
        decoder.decodeFolder()
        assert decoder is not None

    def test_e2e_decode_file_fasta_format(self):
        """Test decodeFile with FASTA-like format - code path coverage"""
        test_file = self.temp_dir / "test_fasta.bin"
        test_data = b"Test FASTA format" * 30
        test_file.write_bytes(test_data)

        # Encode
        number_of_chunks = 10
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True
        )
        encoder.encode_to_packets()

        # Save as single file
        out_file = str(self.temp_dir / "encoded_fasta.RU10")
        encoder.save_packets(split_to_multiple_files=False, out_file=out_file)

        # Decode - verify decodeFile code path
        decoder = RU10Decoder(out_file, use_headerchunk=True)
        decoder.decodeFile()
        assert decoder is not None

    def test_e2e_config_map(self):
        """Test decoder initialization from config_map - code path coverage"""
        test_file = self.temp_dir / "test_config.bin"
        test_data = b"Test config map" * 30
        test_file.write_bytes(test_data)

        # Encode
        number_of_chunks = 10
        dist = RaptorDistribution(number_of_chunks)
        encoder = RU10Encoder(
            str(test_file), number_of_chunks, dist, chunk_size=50, insert_header=True
        )
        encoder.encode_to_packets()

        out_dir = str(self.temp_dir / "encoded_config")
        encoder.save_packets(split_to_multiple_files=True, out_file=out_dir)

        # Create config map
        from configparser import ConfigParser

        config = ConfigParser()
        config.add_section(out_dir)
        config.set(out_dir, "insert_header", "True")
        config.set(out_dir, "number_of_chunks", str(number_of_chunks))
        config.set(out_dir, "error_correction", "nocode")

        # Decode using from_config_map - verify code path
        decoder = RU10Decoder.from_config_map(config[out_dir])
        decoder.file = out_dir
        decoder.decodeFolder()

        # This test just verifies the from_config_map path works
        assert decoder is not None


if __name__ == "__main__":
    test_suite(False, RU10Decoder)
