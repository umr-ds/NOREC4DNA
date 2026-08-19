#!/usr/bin/python
# -*- coding: latin-1 -*-
"""Rule parsing helpers with optional acceleration via the cdnarules extension."""

from collections import Counter
from re import compile, search
from typing import Any, Callable, Dict, List, Pattern, Tuple

from ..helper.fallback_code import (
    longestSequenceOfChar_python,
    microsatellite_python,
    strContainsSub_python,
)

cdnarules: Any = None
try:
    import cdnarules
except ModuleNotFoundError:
    print("C Module failed to load, falling back to slow mode")


def microsatellite(text: str, length_to_look_for: int) -> Tuple[int, str]:
    if cdnarules is None:
        return microsatellite_python(text, length_to_look_for)
    try:
        return cdnarules.microsatellite(text, length_to_look_for)
    except AttributeError:
        return microsatellite_python(text, length_to_look_for)


def longestSequenceOfChar(text: str, char_x: str = "*") -> Tuple[str, int]:
    if cdnarules is None:
        return longestSequenceOfChar_python(text, char_x)
    try:
        return cdnarules.longestSequenceOfChar(text, char_x)
    except AttributeError:
        return longestSequenceOfChar_python(text, char_x)


def strContainsSub(text: str, sequence: str) -> bool:
    return sequence in text


debug = False


# @jit
def switch(name: str) -> Callable[[str, Any, Any], int]:
    switcher: Dict[str, Callable[[str, Any, Any], int]] = {
        "longestSequenceOfChar": (
            lambda x, y, z: 1 if int(z) <= longestSequenceOfChar(x, y)[1] else 0
        ),
        "strContainsSub": (lambda x, y, z: 1 if y in x else 0),
        "strContainsSubRegex": (lambda x, y, z: 1 if strContainsSubRegex(x, y) else 0),
        "strContainsIllegalChars": (lambda x, y, z: 1 if strContainsIllegalChars(x, y) else 0),
        "charCountBiggerEqualThanX": (
            lambda x, y, z: 1 if int(z) <= charCountBiggerEqualThanX(x, y) else 0
        ),
        "microsatelliteLongerThanX": (
            lambda x, y, z: 1 if int(z) <= microsatellite(x, int(y))[0] else 0
        ),
        "len": (lambda x, y, z: 1 if int(y) <= length(x) else 0),
        "gcContent": (lambda x, y, z: 1 if float(y) >= gc_content(x) else 0),
        "gcContentLQ": (lambda x, y, z: 1 if float(y) <= gc_content(x) else 0),
        "*": (lambda x, y, z: 1),
    }
    handler = switcher.get(name)
    if handler is None:
        raise KeyError(f"Unknown rule: {name}")
    return handler


_cdnarules_gc = getattr(cdnarules, "gc_content", None) if cdnarules is not None else None


def gc_content(text: str) -> float:
    if _cdnarules_gc is not None:
        return float(_cdnarules_gc(text))
    return gc_content_python(text)


def gc_content_python(text: str) -> float:
    n = len(text)
    if n == 0:
        return 0.0
    return ((text.count("G") + text.count("C")) / n) * 100


def iupac_replace(sequence: str) -> Pattern[str]:
    iupac_regex = {
        "M": "[AC]",
        "R": "[AG]",
        "W": "[AT]",
        "S": "[CG]",
        "Y": "[CT]",
        "K": "[GT]",
        "V": "[ACG]",
        "H": "[ACT]",
        "D": "[AGT]",
        "B": "[CGT]",
        "X": "[ACGT]",
        "N": "[ACGT]",
    }
    for i, j in iupac_regex.items():
        sequence = sequence.replace(i, j)
    if debug:
        print(sequence)
    return compile(sequence)


def strContainsSubRegex(text: str, sequence: str) -> bool:
    iupac_seq = iupac_replace(sequence)
    res = search(iupac_seq, text)
    if debug:
        print(res)
    return bool(res)


def strContainsIllegalChars(text: str, allowed_chars: str) -> int:
    for cha in text:
        if cha not in allowed_chars:
            return 1
    return 0


def charCountBiggerEqualThanX(text: str, cha: str) -> int:
    res = text.count(cha)
    if debug:
        print(res)
    return res


def length(text: str) -> int:
    res = len(text)
    if debug:
        print(res)
    return res


def executeRule(rule_kind: str, data: str) -> int:
    if "(" not in rule_kind:
        rule_kind += "(*,0)"
    name, params = rule_kind.split("(")
    params = params.split(")")[0].replace(" ", "")
    if params == "":
        params = "*,0"
    if "," not in params:
        params += ",0"
    p1, p2 = params.split(",")
    if p2 == "":
        p2 = "0"
    return switch(name)(data, p1, p2)


def shouldDrop(data: str, rules: List[Tuple[str, float]]) -> float:
    drop_chance = 0.0
    for rule in rules:
        rule_kind, drop_prob = rule
        drop_chance += executeRule(rule_kind, data) * drop_prob
    return min(1.0, drop_chance)


def shouldDropMax(data: str, rules: List[Tuple[str, float]]) -> float:
    drop_chance = 0.0
    for rule in rules:
        rule_kind, drop_prob = rule
        tmp = executeRule(rule_kind, data) * drop_prob
        if drop_chance < tmp:
            drop_chance = tmp
    return min(1.0, drop_chance)


def shouldDropMin(data: str, rules: List[Tuple[str, float]]) -> float:
    drop_chance = 1.0
    for rule in rules:
        rule_kind, drop_prob = rule
        tmp = executeRule(rule_kind, data) * drop_prob
        if drop_chance > tmp:
            drop_chance = tmp
    return min(1.0, drop_chance)


if __name__ == "__main__":
    print(microsatellite("ACGACGAAGAAGAAGAAGAGTAGTAGAAGA", 2))

    in_text = "AAAGCCGAGAGAATTTTTTCACAAAAAAAAAAAAAGTTATAAATCCAATCA" * 10000
    rules1 = [
        ("*", 0.1),
        ("len(15)", 0.2),
        ("strContainsIllegalChars(ACGT)", 1.0),
        ("charCountBiggerEqualThanX(A, 15)", 0.5),
        ("microsatelliteLongerThanX(2,10)", 0.01),
    ]
    print(
        shouldDrop(in_text, rules1)
        + shouldDropMax(
            in_text,
            [
                ("longestSequenceOfChar(*,2)", 0.001),
                ("longestSequenceOfChar(*,3)", 0.005),
                ("longestSequenceOfChar(*,4)", 0.01),
                ("longestSequenceOfChar(*,5)", 0.05),
                ("longestSequenceOfChar(*,6)", 0.1),
                ("longestSequenceOfChar(*,7)", 0.9),
                ("longestSequenceOfChar(*,8)", 1.0),
            ],
        )
    )
    print(
        shouldDrop(
            in_text,
            [
                ("gcContent(50)", 0.001),
                ("gcContent(40)", 0.01),
                ("gcContent(30)", 0.02),
                ("gcContent(20)", 0.03),
                ("gcContent(10)", 0.04),
                ("gcContent(0)", 0.05),
            ],
        )
    )
    print(
        shouldDrop(
            in_text,
            [
                ("strContainsSub(TATAAA)", 0.01),
                ("strContainsSub(TTGACA)", 0.05),
            ],
        )
    )
    print(
        shouldDrop(
            in_text,
            [
                ("strContainsSubRegex(CANYYY)", 0.01),
                ("strContainsSubRegex(ANCCAATCA)", 0.01),
            ],
        )
    )
