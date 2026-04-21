#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
Coverage improvement tests for helper.py module.
"""

import os
import shutil
from pathlib import Path

import pytest
from norec4dna.helper.helper import (
    base_str_to_int,
    calc_crc,
    calc_file_crc,
    cluster_and_remove_index,
    crc_algo_from_str,
    fasta_cluster_and_remove_index,
    find_ceil_power_of_four,
    merge_folder_content,
    merge_parts,
    number_to_base_str,
    split_file,
    split_first,
    xor_with_seed,
)

# Get the directory containing this test file
TEST_DIR = Path(__file__).parent.absolute()
TEMP_DIR = TEST_DIR / "temp_test_helper"


@pytest.fixture(autouse=True)
def setup_and_cleanup():
    """Setup and cleanup for each test"""
    # Setup
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR)
    TEMP_DIR.mkdir()

    yield

    # Cleanup
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR)
    if (TEST_DIR / "cluster_out").exists():
        shutil.rmtree(TEST_DIR / "cluster_out")


class TestSplitFile:
    """Test split_file function"""

    def test_split_file_basic(self):
        """Test basic file splitting"""
        # Create a test file in current directory (split_file creates subdir there)
        test_file = TEMP_DIR / "test_split.txt"
        test_content = b"0123456789" * 100  # 1000 bytes
        test_file.write_bytes(test_content)

        # Change to TEMP_DIR for split_file to work correctly
        import os

        old_cwd = os.getcwd()
        try:
            os.chdir(str(TEMP_DIR))
            # Split into 5 parts
            result = split_file("test_split.txt", 5)

            assert len(result) == 5
            assert all(os.path.exists(f) for f in result)

            # Verify content
            merged = b""
            for f in result:
                merged += open(f, "rb").read()
            assert merged == test_content
        finally:
            os.chdir(old_cwd)


class TestFindCeilPowerOfFour:
    """Test find_ceil_power_of_four function"""

    def test_find_ceil_power_of_four_various(self):
        """Test with various inputs"""
        assert find_ceil_power_of_four(1) == 0  # 4^0 = 1
        assert find_ceil_power_of_four(4) == 1  # 4^1 = 4
        assert find_ceil_power_of_four(5) == 2  # 4^2 = 16
        assert find_ceil_power_of_four(16) == 2
        assert find_ceil_power_of_four(17) == 3  # 4^3 = 64
        assert find_ceil_power_of_four(64) == 3
        assert find_ceil_power_of_four(100) == 4  # 4^4 = 256


class TestNumberToBaseStr:
    """Test number_to_base_str function"""

    def test_number_to_base_str_basic(self):
        """Test basic number to base string conversion"""
        assert number_to_base_str(0, 4) == "AAAA"
        assert number_to_base_str(1, 4) == "AAAC"
        assert number_to_base_str(4, 4) == "AACA"

    def test_number_to_base_str_roundtrip(self):
        """Test roundtrip conversion"""
        for num in [0, 1, 10, 100, 1000]:
            base_str = number_to_base_str(num, 8)
            result = base_str_to_int(base_str)
            assert result == num


class TestBaseStrToInt:
    """Test base_str_to_int function"""

    def test_base_str_to_int_basic(self):
        """Test basic base string to int conversion"""
        assert base_str_to_int("AAAA") == 0
        assert base_str_to_int("AAAC") == 1
        assert base_str_to_int("AACA") == 4
        assert base_str_to_int("TTTT") == 255  # 3 + 3*4 + 3*16 + 3*64 = 255


class TestClusterAndRemoveIndex:
    """Test cluster_and_remove_index function"""

    def test_cluster_and_remove_index_start(self):
        """Test clustering with index at start"""
        import os

        # Create test DNA files
        dna_folder = TEMP_DIR / "dna_files"
        dna_folder.mkdir()

        # The function reads the index from the DNA content (first 4 chars)
        # AAAA=0, AAAC=1, AACA=4 in base-4
        files_content = [
            ("test1.DNA", "AAAAACGT"),  # Index 0 (AAAA at start)
            ("test2.DNA", "AAACTGCA"),  # Index 1 (AAAC at start)
            ("test3.DNA", "AACAGGGG"),  # Index 4 (AACA at start)
        ]

        for filename, content in files_content:
            (dna_folder / filename).write_text(content)

        # Change to TEMP_DIR for cluster_out to be created there
        old_cwd = os.getcwd()
        try:
            os.chdir(str(TEMP_DIR))
            # Cluster
            folders, last_folder = cluster_and_remove_index("start", 4, str(dna_folder))

            # Should have 3 different indices (0, 1, 4)
            assert len(folders) == 3
            assert os.path.exists("cluster_out/0")
            assert os.path.exists("cluster_out/1")
            assert os.path.exists("cluster_out/4")
        finally:
            os.chdir(old_cwd)


class TestFastaClusterAndRemoveIndex:
    """Test fasta_cluster_and_remove_index function"""

    def test_fasta_cluster_basic(self):
        """Test FASTA clustering"""
        # Create test FASTA file
        fasta_file = TEMP_DIR / "test.fasta"
        fasta_content = """>seq1
AAAATTTT
>seq2
AAACGGGG
>seq3
AAAACCCC
"""
        fasta_file.write_text(fasta_content)

        # Cluster
        folders, last_folder = fasta_cluster_and_remove_index("start", 4, str(fasta_file))

        assert len(folders) == 2
        assert os.path.exists("cluster_out/0.fasta")
        assert os.path.exists("cluster_out/1.fasta")


class TestMergeFolderContent:
    """Test merge_folder_content function"""

    def test_merge_folder_basic(self):
        """Test basic folder merging"""
        # Create source folder structure
        src_base = TEMP_DIR / "src"
        src_base.mkdir()

        folder1 = src_base / "folder1"
        folder1.mkdir()
        (folder1 / "file1.txt").write_text("content1")

        folder2 = src_base / "folder2"
        folder2.mkdir()
        (folder2 / "file2.txt").write_text("content2")

        # Create dest folder
        dest = TEMP_DIR / "dest"

        # Merge
        merge_folder_content(
            str(src_base), str(dest), append_folder_name=True, clear_dest_folder=True
        )

        assert (dest / "folder1_file1.txt").exists()
        assert (dest / "folder2_file2.txt").exists()


class TestSplitFirst:
    """Test split_first function"""

    def test_split_first_fasta(self):
        """Test split_first with FASTA file"""
        assert split_first("test_5.fasta") == "5"
        assert split_first("data_100.fasta") == "100"

    def test_split_first_other(self):
        """Test split_first with other file types"""
        assert split_first("test.txt") == "txt"
        assert split_first("data.bin") == "bin"


class TestMergeParts:
    """Test merge_parts function"""

    def test_merge_parts_basic(self):
        """Test merging file parts"""
        base_name = str(TEMP_DIR / "merged")

        # Create parts
        parts_content = [b"part0", b"part1", b"part2"]
        for i, content in enumerate(parts_content):
            with open(f"{base_name}.{i}", "wb") as f:
                f.write(content)

        # Get filenames
        filenames = [f"{base_name}.{i}" for i in range(3)]

        # Merge
        merge_parts(filenames, remove_tmp_on_success=True)

        # Verify merged file
        with open(base_name, "rb") as f:
            merged = f.read()

        assert merged == b"".join(parts_content)
        # Verify parts were removed
        assert not any(os.path.exists(f) for f in filenames)


class TestXorWithSeed:
    """Test xor_with_seed function"""

    def test_xor_with_seed_basic(self):
        """Test basic XOR with seed"""
        data = b"test data"
        seed = 42

        # XOR the data
        xored = xor_with_seed(data, seed)

        # XOR again should give original (with same seed)
        # Note: This test verifies the function runs without error
        # The actual XOR behavior depends on numpy random
        assert len(xored) == len(data)
        assert xored != data  # Should be different


class TestCrcAlgoFromStr:
    """Test crc_algo_from_str function"""

    def test_crc_algo_various(self):
        """Test CRC algorithm selection"""
        algo_b = crc_algo_from_str("B")
        algo_h = crc_algo_from_str("H")
        algo_i = crc_algo_from_str("I")

        # Test they work
        data = b"test"
        assert isinstance(algo_b(data, 0), int)
        assert isinstance(algo_h(data, 0), int)
        assert isinstance(algo_i(data, 0), int)

    def test_crc_algo_invalid(self):
        """Test invalid CRC length string"""
        with pytest.raises(ValueError):
            crc_algo_from_str("X")


class TestCalcFileCrc:
    """Test calc_file_crc function"""

    def test_calc_file_crc_various(self):
        """Test CRC calculation for files"""
        test_file = TEMP_DIR / "crc_test.txt"
        test_file.write_bytes(b"test content")

        crc_i = calc_file_crc(str(test_file), "I")
        crc_h = calc_file_crc(str(test_file), "H")
        crc_b = calc_file_crc(str(test_file), "B")

        assert isinstance(crc_i, int)
        assert isinstance(crc_h, int)
        assert isinstance(crc_b, int)
        assert crc_i != crc_h != crc_b  # Different CRC lengths give different results


class TestCalcCrc:
    """Test calc_crc function"""

    def test_calc_crc_bytes(self):
        """Test CRC calculation for bytes"""
        data = b"test content"
        crc = calc_crc(data, "I")
        assert isinstance(crc, int)

    def test_calc_crc_bytearray(self):
        """Test CRC calculation for bytearray"""
        data = bytearray(b"test content")
        crc = calc_crc(data, "I")
        assert isinstance(crc, int)

    def test_calc_crc_string(self):
        """Test CRC calculation for string"""
        data = "test content"
        crc = calc_crc(data, "I")
        assert isinstance(crc, int)

    def test_calc_crc_file_io(self):
        """Test CRC calculation for file IO"""
        import io

        data = io.BytesIO(b"test content")
        crc = calc_crc(data, "I")
        assert isinstance(crc, int)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
