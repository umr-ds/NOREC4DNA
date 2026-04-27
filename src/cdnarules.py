"""Compatibility wrapper for the packaged ``norec4dna.cdnarules`` extension."""

from norec4dna.cdnarules import (
    bitSet,
    bitsSet,
    buildGraySequence,
    byte2QUATS,
    elimination,
    elimination_with_first_row,
    gc_content,
    getQUAT,
    grayCode,
    longestSequenceOfChar,
    microsatellite,
    repeatRegion,
    smallRepeatRegion,
    strContainsSub,
    xorArray,
)

__all__ = [
    "bitSet",
    "bitsSet",
    "buildGraySequence",
    "byte2QUATS",
    "elimination",
    "elimination_with_first_row",
    "gc_content",
    "getQUAT",
    "grayCode",
    "longestSequenceOfChar",
    "microsatellite",
    "repeatRegion",
    "smallRepeatRegion",
    "strContainsSub",
    "xorArray",
]
