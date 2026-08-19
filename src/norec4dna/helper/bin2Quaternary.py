#!/usr/bin/python
# -*- coding: latin-1 -*-
"""Helpers for converting binary payloads to quaternary DNA bases."""

from __future__ import annotations

import os
import typing
from importlib import import_module
from typing import List, Union


def _fallback_get_quat(bit1: bool, bit2: bool) -> str:
    if not (bit1 or bit2):
        return "A"
    if (not bit1) and bit2:
        return "C"
    if bit1 and (not bit2):
        return "G"
    if bit1 and bit2:
        return "T"
    print("ERROR, this might never happen.")
    return "E"


def _fallback_byte_to_quats(byte: int) -> str:
    res = ""
    if not isinstance(byte, int):
        byte = byte[0]
    if not isinstance(byte, str):
        byt = iter(bin(byte)[2:].rjust(8, "0"))
    else:
        byt = iter(bin(ord(byte))[2:].rjust(8, "0"))
    for x, y in zip(byt, byt):
        res += _fallback_get_quat(str2bool(x), str2bool(y))
    return res


def _load_quat_helpers() -> (
    typing.Tuple[typing.Callable[[int], str], typing.Callable[[bool, bool], str]]
):
    try:
        cdnarules = import_module("cdnarules")
        return cdnarules.byte2QUATS, cdnarules.getQUAT
    except ImportError:
        print("C Module failed to load, falling back to slow mode")
        return _fallback_byte_to_quats, _fallback_get_quat


byte2QUATS, getQUAT = _load_quat_helpers()
_BYTE2QUATS_LOOKUP: typing.Tuple[str, ...] = tuple(byte2QUATS(i) for i in range(256))


def bin2Quaternary(filename: str) -> None:
    with open(filename, "rb") as f:
        with open(filename + ".quat", "w") as ff:
            byte = int(f.read(1))
            while byte != b"":
                ff.write(byte2QUATS(byte))
                byte = int(f.read(1))


def string2QUATS(text: typing.Union[str, bytes]) -> List[str]:
    if isinstance(text, str):
        text = text.encode("utf-8")
    lookup = _BYTE2QUATS_LOOKUP
    return [lookup[x] for x in text]


def str2bool(s: str) -> bool:
    return s == "1"


def quads2dna(quads: Union[List[int], bytes]) -> str:
    translation = {0: "A", 1: "C", 2: "G", 3: "T"}
    return "".join(translation[x] for x in quads)


def main() -> None:
    folder = "RU10_b_lq.webm"
    for filename in os.listdir(folder):
        if filename.endswith(".RU10"):
            bin2Quaternary(folder + "/" + filename)


if __name__ == "__main__":
    # main()
    filename = "vergleich bzgl ACGT-verteilung/mit dna rules/RU10_logo.jpg/0.RU10"  # ""RU10_b_lq.webm/1.RU10"  # "raptor.pdf"
    x = [
        0,
        0,
        0,
        1,
        2,
        3,
        1,
        1,
        1,
    ]
    print(quads2dna(x))

    def bitstring_to_bytes(s: Union[str, bytes, bytearray]) -> bytes:
        v = int(s, 2)
        b = bytearray()
        while v:
            b.append(v & 0xFF)
            v >>= 8
        return bytes(b[::-1])

    s = (
        "0010110001101001010111100010101100001101101011100111110001010001"
        "00000011111110110001011010010011101110001111101010110100"
    )
    print("".join(string2QUATS(bitstring_to_bytes(s))))
