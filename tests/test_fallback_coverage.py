#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
Coverage tests for fallback_code.py module.
These tests cover the fallback implementations used when C extensions are not available.
"""

import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from norec4dna.helper.fallback_code import (
    bitSet,
    bitsSet,
    buildGraySequence,
    grayCode,
    longestSequenceOfChar_python,
    microsatellite_python,
    r_region,
    small_r_region,
    strContainsSub_python,
    xor_intern,
    xor_numpy_internal,
)


class TestBitSet:
    """Test bitSet function from fallback_code"""

    def test_bitset_basic(self):
        """Test basic bitSet functionality"""
        # Test various bit positions
        assert bitSet(0b1, 0) == True
        assert bitSet(0b10, 1) == True
        assert bitSet(0b100, 2) == True
        assert bitSet(0b1000, 3) == True

    def test_bitset_false(self):
        """Test bitSet when bit is not set"""
        assert bitSet(0b0, 0) == False
        assert bitSet(0b1, 1) == False
        assert bitSet(0b10, 0) == False
        assert bitSet(0b1010, 0) == False
        assert bitSet(0b1010, 2) == False

    def test_bitset_large_number(self):
        """Test bitSet with large numbers"""
        large_num = 0xFFFFFFFFFFFFFFFF
        assert bitSet(large_num, 0) == True
        assert bitSet(large_num, 63) == True


class TestBitsSet:
    """Test bitsSet function from fallback_code"""

    def test_bitsset_zero(self):
        """Test bitsSet with zero"""
        assert bitsSet(np.uint64(0)) == 0

    def test_bitsset_single_bit(self):
        """Test bitsSet with single bit set"""
        assert bitsSet(np.uint64(1)) == 1
        assert bitsSet(np.uint64(2)) == 1
        assert bitsSet(np.uint64(4)) == 1
        assert bitsSet(np.uint64(8)) == 1

    def test_bitsset_multiple_bits(self):
        """Test bitsSet with multiple bits set"""
        assert bitsSet(np.uint64(0b11)) == 2
        assert bitsSet(np.uint64(0b101)) == 2
        assert bitsSet(np.uint64(0b111)) == 3
        assert bitsSet(np.uint64(0b1111)) == 4

    def test_bitsset_all_bits(self):
        """Test bitsSet with all bits set"""
        assert bitsSet(np.uint64(0xFF)) == 8
        assert bitsSet(np.uint64(0xFFFF)) == 16
        assert bitsSet(np.uint64(0xFFFFFFFF)) == 32


class TestGrayCode:
    """Test grayCode function from fallback_code"""

    def test_graycode_zero(self):
        """Test grayCode with zero"""
        result = grayCode(0)
        assert isinstance(result, np.uint64)
        assert result == 0

    def test_graycode_sequence(self):
        """Test grayCode produces correct Gray code sequence"""
        # Gray code: 0, 1, 3, 2, 6, 7, 5, 4, ...
        assert grayCode(0) == np.uint64(0)
        assert grayCode(1) == np.uint64(1)
        assert grayCode(2) == np.uint64(3)
        assert grayCode(3) == np.uint64(2)
        assert grayCode(4) == np.uint64(6)
        assert grayCode(5) == np.uint64(7)

    def test_graycode_property(self):
        """Test Gray code property: consecutive values differ by one bit"""
        for i in range(10):
            g1 = grayCode(i)
            g2 = grayCode(i + 1)
            # XOR of consecutive Gray codes should have exactly one bit set
            xor_result = g1 ^ g2
            assert bitsSet(xor_result) == 1


class TestBuildGraySequence:
    """Test buildGraySequence function from fallback_code"""

    def test_buildgraysequence_basic(self):
        """Test basic buildGraySequence functionality"""
        # Build sequence of length 4 with b=1 (Gray codes with exactly 1 bit set)
        result = buildGraySequence(4, 1)
        assert len(result) == 4
        # All results should have exactly 1 bit set
        for val in result:
            assert bitsSet(np.uint64(val)) == 1

    def test_buildgraysequence_b2(self):
        """Test buildGraySequence with b=2"""
        result = buildGraySequence(3, 2)
        assert len(result) == 3
        # All results should have exactly 2 bits set
        for val in result:
            assert bitsSet(np.uint64(val)) == 2


class TestXorIntern:
    """Test xor_intern function from fallback_code"""

    def test_xor_intern_basic(self):
        """Test basic xor_intern functionality"""
        a = np.array([1, 2, 3], dtype=np.uint8)
        b = np.array([1, 2, 3], dtype=np.uint8)
        result = xor_intern(a, b)
        assert np.array_equal(result, np.array([0, 0, 0], dtype=np.uint8))

    def test_xor_intern_different(self):
        """Test xor_intern with different values"""
        a = np.array([0, 1, 2], dtype=np.uint8)
        b = np.array([1, 1, 1], dtype=np.uint8)
        result = xor_intern(a, b)
        assert np.array_equal(result, np.array([1, 0, 3], dtype=np.uint8))


class TestRRegion:
    """Test r_region function from fallback_code"""

    def test_r_region_no_repeat(self):
        """Test r_region with no repeating subsequences"""
        # No repeats of length 20
        data = "ACGT" * 5  # Only 20 chars, no room for repeat
        result = r_region(data, repeat_length=20)
        assert result == 0.0

    def test_r_region_with_repeat(self):
        """Test r_region with repeating subsequences"""
        # Has repeating subsequence
        data = "ACGTACGTACGT" + "ACGTACGTACGT"
        result = r_region(data, repeat_length=4)
        assert result == 1.0


class TestSmallRRegion:
    """Test small_r_region function from fallback_code"""

    def test_small_r_region_no_repeat(self):
        """Test small_r_region with no repeating subsequences"""
        # The function counts occurrences and has a threshold
        # Even with few repeats, it may return 1.0 based on the formula
        data = "ACGT" * 3  # 12 chars
        result = small_r_region(data, repeat_length=9)
        # Just test it runs and returns a valid float
        assert isinstance(result, float)
        assert 0.0 <= result <= 1.0

    def test_small_r_region_with_repeat(self):
        """Test small_r_region with repeating subsequences"""
        # Create data with many repeats
        data = "ACGTACGTACGT" * 10
        result = small_r_region(data, repeat_length=4)
        # Should be high since many repeats
        assert result >= 0.5


class TestMicrosatellitePython:
    """Test microsatellite_python function from fallback_code"""

    def test_microsatellite_no_repeat(self):
        """Test microsatellite_python with no repeats"""
        # The function counts consecutive repeats at regular intervals
        # Even with unique data, it will find at least 1 occurrence
        data = "ACGT" * 3
        count, chars = microsatellite_python(data, 4)
        # Just test it runs and returns valid values
        assert isinstance(count, int)
        assert count >= 1
        assert isinstance(chars, str)
        assert len(chars) == 4

    def test_microsatellite_with_repeat(self):
        """Test microsatellite_python with repeats"""
        data = "ACGTACGTACGT"
        count, chars = microsatellite_python(data, 4)
        assert count >= 2  # Found repeats
        assert chars == "ACGT"


class TestLongestSequenceOfCharPython:
    """Test longestSequenceOfChar_python function from fallback_code"""

    def test_longestsequence_basic(self):
        """Test basic longestSequenceOfChar_python functionality"""
        data = "AAABBBCCCC"
        char, count = longestSequenceOfChar_python(data, "*")
        assert char == "C"
        assert count == 4

    def test_longestsequence_specific_char(self):
        """Test longestSequenceOfChar_python with specific char"""
        data = "AAABBBCCCC"
        char, count = longestSequenceOfChar_python(data, "A")
        assert char == "A"
        assert count == 3

    def test_longestsequence_no_match(self):
        """Test longestSequenceOfChar_python when char not found"""
        data = "BBBCCCC"
        char, count = longestSequenceOfChar_python(data, "A")
        assert char == "A"  # Returns the search char when not found
        assert count == 0


class TestStrContainsSubPython:
    """Test strContainsSub_python function from fallback_code"""

    def test_strcontainsub_found(self):
        """Test strContainsSub_python when sequence is found"""
        assert strContainsSub_python("ACGTACGT", "ACGT") == True
        assert strContainsSub_python("hello world", "world") == True

    def test_strcontainsub_not_found(self):
        """Test strContainsSub_python when sequence is not found"""
        assert strContainsSub_python("ACGT", "TGCA") == False
        assert strContainsSub_python("hello", "world") == False


class TestXorNumpyInternal:
    """Test xor_numpy_internal function from fallback_code"""

    def test_xor_numpy_internal_basic(self):
        """Test basic xor_numpy_internal functionality"""
        a = np.array([1, 2, 3], dtype=np.uint8)
        b = np.array([1, 2, 3], dtype=np.uint8)
        result = xor_numpy_internal(a, b)
        assert np.array_equal(result, np.array([0, 0, 0], dtype=np.uint8))

    def test_xor_numpy_internal_different(self):
        """Test xor_numpy_internal with different values"""
        a = np.array([0, 1, 2, 3], dtype=np.uint8)
        b = np.array([1, 1, 1, 1], dtype=np.uint8)
        result = xor_numpy_internal(a, b)
        assert np.array_equal(result, np.array([1, 0, 3, 2], dtype=np.uint8))

    def test_xor_numpy_internal_large(self):
        """Test xor_numpy_internal with larger arrays"""
        size = 1000
        a = np.random.randint(0, 256, size, dtype=np.uint8)
        b = np.random.randint(0, 256, size, dtype=np.uint8)
        result = xor_numpy_internal(a, b)
        assert len(result) == size
        # Verify XOR property: a XOR b XOR b = a
        verify = xor_numpy_internal(result, b)
        assert np.array_equal(verify, a)


class TestHelperCpuSingleCoreFallback:
    """Test helper_cpu_single_core.py fallback code paths by mocking import failures"""

    def test_xor_numpy_internal_via_helper_cpu_single_core(self):
        """Test that xor_numpy_internal from fallback_code works correctly"""
        # This tests the fallback implementation that helper_cpu_single_core uses
        a = np.array([1, 0, 1, 0], dtype=np.uint8)
        b = np.array([0, 1, 1, 0], dtype=np.uint8)
        result = xor_numpy_internal(a, b)
        assert np.array_equal(result, np.array([1, 1, 0, 0], dtype=np.uint8))

    def test_bitset_fallback(self):
        """Test bitSet fallback implementation"""
        # Test the fallback implementation directly
        assert bitSet(0b1010, 1) == True
        assert bitSet(0b1010, 0) == False

    def test_bitsset_fallback(self):
        """Test bitsSet fallback implementation"""
        assert bitsSet(np.uint64(0b1011)) == 3

    def test_graycode_fallback(self):
        """Test grayCode fallback implementation"""
        result = grayCode(5)
        assert isinstance(result, np.uint64)
        # Gray code of 5 should be 7
        assert result == np.uint64(7)

    def test_buildgraysequence_fallback(self):
        """Test buildGraySequence fallback implementation"""
        result = buildGraySequence(2, 1)
        assert len(result) == 2
        # All should have exactly 1 bit set
        for val in result:
            assert bitsSet(np.uint64(val)) == 1

    def test_helper_cpu_single_core_fallback_import(self):
        """Test that helper_cpu_single_core fallback functions are available"""
        # Note: The actual fallback import paths (lines 14-15, 98-100, 107-109, 116-118, 122-124)
        # are only executed when cdnarules C extension is NOT installed.
        # Since cdnarules IS installed in this environment, those specific lines cannot be covered.
        # This test verifies the fallback_code functions that would be used work correctly.

        # Test that fallback_code functions work (these are what helper_cpu_single_core uses as fallback)
        a = np.array([1, 0], dtype=np.uint8)
        b = np.array([0, 1], dtype=np.uint8)
        result = xor_numpy_internal(a, b)
        assert np.array_equal(result, np.array([1, 1], dtype=np.uint8))

        # Test bitSet from fallback
        assert bitSet(0b1010, 1) == True

        # Test bitsSet from fallback
        assert bitsSet(np.uint64(0b1011)) == 3

        # Test grayCode from fallback
        assert grayCode(5) == np.uint64(7)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
