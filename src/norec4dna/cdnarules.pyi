"""
Type stub file for cdnarules C extension module.

This module provides optimized DNA sequence processing functions implemented in C
with full type safety and AVX2/SSE2 vectorization support.

Author: DR4DNA Team
License: MIT
"""

from typing import Any, Tuple

import numpy as np

UInt64Array = np.ndarray[Any, np.dtype[np.uint64]]
UInt8Array = np.ndarray[Any, np.dtype[np.uint8]]
BoolArray = np.ndarray[Any, np.dtype[np.bool_]]
IntpArray = np.ndarray[Any, np.dtype[np.intp]]

# ============================================================================
# DNA Sequence Analysis Functions
# ============================================================================

def microsatellite(text: str, lengthToLookFor: int) -> Tuple[int, str]:
    """
    Find the maximum microsatellite repeat of given length in text.

    A microsatellite is a repeating sequence of DNA bases.

    Args:
        text: DNA sequence string containing only A, C, G, T
        lengthToLookFor: Length of subsequence to search for (must be positive)

    Returns:
        Tuple of (count, subsequence) where:
            - count: Number of times the most frequent subsequence repeats
            - subsequence: The most frequently repeating subsequence

    Raises:
        ValueError: If lengthToLookFor is not positive

    Example:
        >>> microsatellite("ACGACGAAGAAGAAGAAGAG", 2)
        (4, 'AG')
    """
    ...

def longestSequenceOfChar(text: str, character_to_look_for: str = "*") -> Tuple[str, int]:
    """
    Find the longest consecutive sequence of a character in text.

    Args:
        text: DNA sequence string
        character_to_look_for: Character to search for, or "*" for any character

    Returns:
        Tuple of (character, length) where:
            - character: The character with the longest sequence
            - length: Length of the longest sequence

    Raises:
        ValueError: If character_to_look_for is not a single character

    Example:
        >>> longestSequenceOfChar("AAACCCGGGTTT", "*")
        ('A', 3)
        >>> longestSequenceOfChar("AAACCCGGGTTT", "G")
        ('G', 3)
    """
    ...

def repeatRegion(text: str, lengthToLookFor: int) -> int:
    """
    Check if a region of given length is repeated within text.

    Detects duplicate subsequences (including overlapping repeats).

    Args:
        text: DNA sequence string
        lengthToLookFor: Length of region to search for (must be positive)

    Returns:
        1 if any repeat is found, 0 otherwise

    Raises:
        ValueError: If lengthToLookFor is not positive

    Example:
        >>> repeatRegion("ACGTACGT", 4)
        1
        >>> repeatRegion("ACGT", 4)
        0
    """
    ...

def smallRepeatRegion(text: str, lengthToLookFor: int) -> float:
    """
    Calculate error value based on repeat density in text.

    Returns a score indicating the likelihood of sequencing errors
    due to repeat regions.

    Args:
        text: DNA sequence string
        lengthToLookFor: Length of region to search for (must be positive)

    Returns:
        Error value between 0.0 and 1.0:
            - 0.0: No repeats detected
            - 1.0: High repeat density (likely sequencing errors)

    Raises:
        ValueError: If lengthToLookFor is not positive

    Example:
        >>> smallRepeatRegion("AAAAAAAAAA", 2)
        1.0
    """
    ...

def strContainsSub(text: str, substr: str) -> bool:
    """
    Check if substr is present in text.

    Args:
        text: DNA sequence string to search in
        substr: Subsequence to search for

    Returns:
        True if substr is found in text, False otherwise

    Example:
        >>> strContainsSub("ACGTACGT", "CGT")
        True
        >>> strContainsSub("ACGT", "CGTA")
        False
    """
    ...

def gc_content(text: str) -> float:
    """
    Calculate GC content percentage of DNA sequence.

    GC content is the percentage of guanine (G) and cytosine (C) bases
    in a DNA sequence. Important for assessing sequence stability.

    Args:
        text: DNA sequence string (case-insensitive)

    Returns:
        Percentage of G and C bases (0.0 to 100.0)

    Example:
        >>> gc_content("GGCC")
        100.0
        >>> gc_content("AATT")
        0.0
        >>> gc_content("ACGT")
        50.0
    """
    ...

# ============================================================================
# Bit Manipulation Functions
# ============================================================================

def bitSet(x: int, b: int) -> bool:
    """
    Check if bit b is set in integer x.

    Args:
        x: Integer value (treated as 64-bit unsigned)
        b: Bit position to check (0-63, where 0 is least significant)

    Returns:
        True if bit b is set (1), False if clear (0)

    Raises:
        ValueError: If b is not in range 0-63

    Example:
        >>> bitSet(0b1010, 1)
        True
        >>> bitSet(0b1010, 0)
        False
    """
    ...

def bitsSet(x: int) -> int:
    """
    Count the number of bits set to 1 in integer x (population count).

    Uses hardware POPCNT instruction when available for optimal performance.

    Args:
        x: Integer value (treated as 64-bit unsigned)

    Returns:
        Number of bits set to 1 (0-64)

    Example:
        >>> bitsSet(0b1011)
        3
        >>> bitsSet(0xFFFF)
        16
    """
    ...

def grayCode(x: int) -> int:
    """
    Convert integer to Gray code.

    Gray code is a binary numeral system where two successive values
    differ in only one bit. Useful for error correction and digital communications.

    Args:
        x: Integer value to convert

    Returns:
        Gray code representation of x

    Example:
        >>> grayCode(0)
        0
        >>> grayCode(1)
        1
        >>> grayCode(2)
        3
        >>> grayCode(3)
        2
    """
    ...

def buildGraySequence(length: int, b: int) -> UInt64Array:
    """
    Generate Gray code sequence of given length with exactly b bits set.

    Args:
        length: Number of Gray codes to generate (must be positive)
        b: Number of bits that should be set in each code (0-64)

    Returns:
        Numpy array of uint64 Gray codes with exactly b bits set

    Raises:
        ValueError: If length is not positive
        MemoryError: If insufficient memory for result array

    Note:
        May return fewer than 'length' values if not enough Gray codes
        exist with exactly b bits set.

    Example:
        >>> buildGraySequence(5, 2)
        array([ 3,  5,  6,  9, 10], dtype=uint64)
    """
    ...

# ============================================================================
# DNA Encoding/Decoding Functions
# ============================================================================

def getQUAT(bit1: bool, bit2: bool) -> str:
    """
    Convert two bits to DNA base using standard encoding.

    Encoding mapping:
        - (0, 0) → 'A' (Adenine)
        - (0, 1) → 'C' (Cytosine)
        - (1, 0) → 'G' (Guanine)
        - (1, 1) → 'T' (Thymine)

    Args:
        bit1: First bit (most significant)
        bit2: Second bit (least significant)

    Returns:
        Single character string: 'A', 'C', 'G', or 'T'

    Example:
        >>> getQUAT(False, False)
        'A'
        >>> getQUAT(True, True)
        'T'
    """
    ...

def byte2QUATS(byte: int) -> str:
    """
    Convert a byte to 4-character DNA quaternary representation.

    Encodes 8 bits as 4 DNA bases (2 bits per base).

    Args:
        byte: Integer byte value (0-255)

    Returns:
        String of 4 DNA bases representing the byte

    Example:
        >>> byte2QUATS(0b11100100)
        'TGA'
    """
    ...

# ============================================================================
# Array Operations
# ============================================================================

def xorArray(X: UInt8Array, Y: UInt8Array) -> UInt8Array:
    """
    Perform element-wise XOR on two uint8 numpy arrays.

    Uses AVX2/SSE2 vectorization when available for optimal performance.

    Args:
        X: First input array (1D, dtype=uint8)
        Y: Second input array (1D, dtype=uint8, same length as X)

    Returns:
        Numpy array of uint8 with element-wise XOR results

    Raises:
        ValueError: If arrays have different lengths or wrong dtype

    Example:
        >>> import numpy as np
        >>> X = np.array([0b1010, 0b1100], dtype=np.uint8)
        >>> Y = np.array([0b1001, 0b1010], dtype=np.uint8)
        >>> xorArray(X, Y)
        array([ 3,  6], dtype=uint8)
    """
    ...

# ============================================================================
# Gaussian Elimination Functions
# ============================================================================

def elimination(
    A: BoolArray,
    b: UInt8Array,
    packet_mapping: IntpArray,
    chunk_to_used_packets: BoolArray,
) -> bool:
    """
    Perform Gaussian elimination with partial pivoting on matrix system.

    Solves the system Ax = b using XOR-based Gaussian elimination over GF(2).
    Modifies A, b, packet_mapping, and chunk_to_used_packets in-place.

    Args:
        A: Boolean matrix of shape (rows, cols) representing the system
        b: uint8 data matrix of shape (rows, data_width) containing payload data
        packet_mapping: Integer array of shape (rows,) mapping rows to packet indices
        chunk_to_used_packets: Square boolean matrix of shape (n, n) where
                              n >= max(rows, cols), tracks chunk-to-packet mappings

    Returns:
        True if system was successfully solved (all columns have pivots),
        False if system is unsolvable (rank-deficient)

    Raises:
        ValueError: If array dimensions are incompatible

    Note:
        This function uses XOR operations instead of standard arithmetic,
        making it suitable for binary linear systems in DNA storage decoding.

    Example:
        >>> import numpy as np
        >>> A = np.array([[True, False], [True, True]], dtype=bool)
        >>> b = np.array([[1], [0]], dtype=np.uint8)
        >>> pm = np.array([0, 1], dtype=np.intp)
        >>> c2p = np.eye(2, dtype=bool)
        >>> elimination(A, b, pm, c2p)
        True
    """
    ...

def elimination_with_first_row(
    A: BoolArray,
    b: UInt8Array,
    packet_mapping: IntpArray,
    chunk_to_used_packets: BoolArray,
    first_row_idx: int = -1,
) -> bool:
    """
    Perform Gaussian elimination with optional first row pivot selection.

    Extended version of elimination() that allows specifying a particular
    row to use as the first pivot. Useful when certain packets are known
    to be more reliable.

    Args:
        A: Boolean matrix of shape (rows, cols)
        b: uint8 data matrix of shape (rows, data_width)
        packet_mapping: Integer array of shape (rows,)
        chunk_to_used_packets: Square boolean matrix of shape (n, n)
        first_row_idx: Index of row to use as first pivot.
                      If -1 or if A[first_row_idx, 0] is False, uses default pivoting.

    Returns:
        True if system was successfully solved, False otherwise

    Raises:
        ValueError: If first_row_idx is out of bounds or arrays incompatible

    Note:
        If first_row_idx is specified and valid, that row is used as the
        initial pivot and all other rows with a 1 in column 0 are XORed
        with it. Standard elimination then proceeds from column 1.

    Example:
        >>> # Use row 2 as the first pivot
        >>> elimination_with_first_row(A, b, pm, c2p, first_row_idx=2)
        True
    """
    ...

# ============================================================================
# Module Information
# ============================================================================

__doc__ = """
cdnarules - Optimized DNA Sequence Processing Module

This C extension module provides high-performance functions for DNA sequence
analysis and fountain code decoding. Key features include:

- **Vectorized Operations**: AVX2/SSE2 acceleration for XOR operations
- **Hardware Optimization**: POPCNT instruction for bit counting
- **Memory Efficiency**: Aligned memory allocation for SIMD operations
- **Type Safety**: Full type stub support for static type checkers

Functions are organized into categories:
1. DNA sequence analysis (microsatellite detection, GC content)
2. Bit manipulation (Gray codes, population count)
3. Array operations (XOR with SIMD acceleration)
4. Linear algebra (Gaussian elimination over GF(2))

All functions handle error conditions gracefully and raise appropriate
Python exceptions for invalid inputs.
"""

__author__ = "DR4DNA Team"
__license__ = "MIT"
__version__ = "1.0.0"
