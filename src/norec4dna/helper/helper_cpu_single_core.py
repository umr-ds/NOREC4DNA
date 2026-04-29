#!/usr/bin/python
# -*- coding: latin-1 -*-
from __future__ import annotations

import typing
from functools import reduce
from random import random
from typing import TYPE_CHECKING, Any, Union, overload

import numpy
from crccheck.crc import Crc8Lte as crc8
from crccheck.crc import Crc16, Crc32, Crc64

from .fallback_code import buildGraySequence as fallback_build_gray_sequence

UInt8Array = Any
Int64Array = Any
BoolArray = Any

try:
    from cdnarules import xorArray as xor_numpy_uint8_internal
except ImportError:
    xor_numpy_uint8_internal = None

if TYPE_CHECKING:
    from norec4dna.Packet import Packet


@overload
def xor_numpy(p1: UInt8Array, p2: UInt8Array) -> UInt8Array: ...


@overload
def xor_numpy(p1: Int64Array, p2: Int64Array) -> Int64Array: ...


@overload
def xor_numpy(p1: BoolArray, p2: BoolArray) -> BoolArray: ...


@overload
def xor_numpy(
    p1: Union[bytes, bytearray, UInt8Array],
    p2: Union[bytes, bytearray, UInt8Array],
) -> UInt8Array: ...


def xor_numpy(
    p1: Union[bytes, bytearray, UInt8Array, Int64Array, BoolArray],
    p2: Union[bytes, bytearray, UInt8Array, Int64Array, BoolArray],
) -> Union[UInt8Array, Int64Array, BoolArray]:
    if (isinstance(p2, numpy.ndarray) and isinstance(p1, numpy.ndarray)) and (
        (p1.dtype == numpy.uint8 and p2.dtype == numpy.uint8)
        or (p1.dtype == numpy.int64 and p2.dtype == numpy.int64)
        or (p1.dtype == bool and p2.dtype == bool)
    ):
        return numpy.bitwise_xor(p1, p2)
    n_p1 = numpy.frombuffer(p1, dtype="uint8")
    n_p2 = numpy.frombuffer(p2, dtype="uint8")
    if xor_numpy_uint8_internal is not None:
        return xor_numpy_uint8_internal(n_p1, n_p2)
    return numpy.bitwise_xor(n_p1, n_p2)


def listXOR(plist: list) -> Any:
    return reduce(xor_numpy, plist)


def logical_xor(plist: typing.List[typing.List[bool]]) -> BoolArray:
    return numpy.logical_xor.reduce(plist)


def xor_pakets(packet1: str, packet2: str) -> list:
    assert len(packet1) == len(packet2)
    a = [a ^ b for (a, b) in zip(bytes(packet1, "utf-8"), bytes(packet2, "utf-8"))]
    return a


def should_drop_packet(
    rules: Any, packet: "Packet", upper_bound: float = 1.0, limit_only: bool = True
) -> bool:
    rand = upper_bound * random()  # create number from [0, upper_bound)
    drop_chance = rules.apply_all_rules(packet)
    if isinstance(drop_chance, list):
        drop_chance = drop_chance[0]
    packet.set_error_prob(drop_chance)
    # print(str(rand) + " , " + str(drop_chance))
    # drop packet if rand bigger than the drop_chance for this Packet.
    return (drop_chance > upper_bound) if limit_only else (drop_chance > rand)


def calc_crc(data: bytes, crc_len_format: str = "B") -> int:
    if crc_len_format == "L":
        return Crc32.calc(data)
    elif crc_len_format == "H":
        return Crc16.calc(data)
    elif crc_len_format == "B":
        return crc8.calc(data)
    elif crc_len_format == "Q":
        return Crc64.calc(data)
    else:
        raise ValueError("Unknown crc_len_format: " + str(crc_len_format))


T = typing.TypeVar("T", int, bytes, typing.Literal[12])


def xor_mask(
    data: T,
    len_format: str = "I",
    mask: int = 0b11111001110000110110111110011100,
    enabled: bool = True,
) -> T:
    if not enabled:
        return data
    if len_format == "B":
        return data
    if len_format == "H":
        mask = 0b1111100111000011
    if len_format == "I":
        mask = 0b11111001110000110110111110011100
    if len_format == "Q":
        mask = 0b1111100111000011011011111001110011111001110000110110111110011100
    with numpy.errstate(over="ignore"):
        return numpy.bitwise_xor(data, mask)


try:
    from cdnarules import bitSet as bitSet_c

    def bitSet(x: int, b: int) -> bool:
        return bitSet_c(int(x), int(b))

except ImportError:
    print("BitSet - C Module failed to load, falling back to slow mode")

try:
    from cdnarules import bitsSet as bitsSet_c

    def bitsSet(x: int) -> int:
        return bitsSet_c(int(x))

except ImportError:
    print("BitsSet - C Module failed to load, falling back to slow mode")

try:
    from cdnarules import grayCode as grayCode_c

    def grayCode(x: int) -> typing.Any:
        return numpy.uint64(grayCode_c(int(x)))

except ImportError:
    print("Gray-Code - C Module failed to load, falling back to slow mode")

try:
    from cdnarules import buildGraySequence as cdnarules_build_gray_sequence

    buildGraySequence = cdnarules_build_gray_sequence
except ImportError:
    print("Graysequence - C Module failed to load, falling back to slow mode")
    buildGraySequence = fallback_build_gray_sequence
