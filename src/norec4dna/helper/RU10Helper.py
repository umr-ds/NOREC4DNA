#!/usr/bin/python
# -*- coding: latin-1 -*-
from __future__ import annotations

import typing
from functools import lru_cache
from math import ceil, floor, log, pow, sqrt

import numpy as np

from ..distributions.RaptorDistribution import RaptorDistribution

BoolArray = typing.Any


class RandomWithRandInt(typing.Protocol):
    def seed(self, seed: int) -> None: ...

    def randint(self, low: int, high: typing.Optional[int] = None) -> int: ...


int63 = int(pow(2, 63) - 1)
int31 = int(pow(2, 31) - 1)


_PACKET_NUMBERS_CACHE: typing.Dict[typing.Tuple[int, int, bool, typing.Optional[int]], typing.List[int]] = {}
_TRIPLE_CACHE: typing.Dict[typing.Tuple[int, int], typing.Tuple[int, int, int]] = {}


def choose_packet_numbers(
    number_of_chunks: int,
    code_block_index: int,
    dist: RaptorDistribution,
    systematic: bool = False,
    max_l: typing.Optional[int] = None,
) -> typing.List[int]:
    cache_key = (number_of_chunks, code_block_index, systematic, max_l)
    cached = _PACKET_NUMBERS_CACHE.get(cache_key)
    if cached is not None:
        return list(cached)

    if systematic:
        d, a, b = systematic_ru10_triple_generator(number_of_chunks, code_block_index, dist)
    else:
        d, a, b = ru10_triple_generator(number_of_chunks, code_block_index, dist, max_l)
    if max_l is None:
        total_intermediate_symbols, _, _ = intermediate_symbols(number_of_chunks, dist)
    else:
        total_intermediate_symbols = max_l
    lprime: int = int(np.uint32(dist.smallestPrimeGreaterOrEqual(total_intermediate_symbols)))

    if d > total_intermediate_symbols:
        d = total_intermediate_symbols
    total_intermediate_symbols = np.uint32(total_intermediate_symbols)
    indices: typing.List[int] = [0] * d
    while b >= total_intermediate_symbols:
        b = (b + a) % lprime

    indices[0] = b

    for idx in range(1, d):
        b = (b + a) % lprime
        while b >= total_intermediate_symbols:
            b = (b + a) % lprime
        indices[idx] = b
    res = sorted(indices)
    _PACKET_NUMBERS_CACHE[cache_key] = res
    return list(res)


@lru_cache(maxsize=None)
def intermediate_symbols(k, dist) -> typing.Tuple[int, int, int]:
    # X is the smallest positive integer such that X*(X-1) >= 2*K
    x = int(floor(sqrt(2 * np.float64(k))))
    if x < 1:
        x = 1

    while (x * (x - 1)) < (2 * k):
        x += 1
    s = int(ceil(0.01 * np.float64(k))) + x
    s = dist.smallestPrimeGreaterOrEqual(s)
    h = int(floor(log(np.float64(s) + np.float64(k)) / log(4)))
    while dist.centerBinomial(h) < k + s:
        h += 1
    return k + s + h, s, h


def ru10_triple_generator(
    k: int, x: int, dist: RaptorDistribution, max_l: typing.Optional[int] = None
) -> typing.Tuple[int, int, int]:
    if max_l is None:
        total_intermediate_symbols, _, _ = intermediate_symbols(k, dist)
    else:
        total_intermediate_symbols = max_l
    lprime = dist.smallestPrimeGreaterOrEqual(total_intermediate_symbols)
    cache_key = (lprime, x)
    cached = _TRIPLE_CACHE.get(cache_key)
    if cached is not None:
        return cached
    rng: RandomWithRandInt = np.random
    rng.seed(x)
    v = np.uint32(r_int63(rng) % 1048576)
    a = np.uint32(1 + (r_int63(rng) % (lprime - 1)))
    b = np.uint32(r_int63(rng) % lprime)
    d = dist.deg(int(v))
    res = (d, int(a), int(b))
    _TRIPLE_CACHE[cache_key] = res
    return res


def systematic_ru10_triple_generator(
    k: int, x: int, dist: RaptorDistribution
) -> typing.Tuple[int, int, int]:
    l, _, _ = intermediate_symbols(k, dist)
    lprime = dist.smallestPrimeGreaterOrEqual(l)
    q = 65521  # largest prime < 2 ^ 16
    jk = dist.systematicIndextable[k]

    a = int(53591 + jk * 997) % int(q)
    b = 10267 * (jk + 1) % q
    y = int(b + (x * a)) % q
    v: int = dist.raptor_rand(y, 0, 1048576)
    d: int = dist.deg(v)
    a: int = 1 + dist.raptor_rand(y, 1, int(lprime - 1))
    b: int = dist.raptor_rand(y, 2, int(lprime))
    return d, a, b


def r_int63(rng: RandomWithRandInt) -> int:
    return rng.randint(0, int31)


def from_true_false_list(
    tf_list: typing.Union[typing.List[bool], BoolArray],
) -> typing.List[int]:
    if isinstance(tf_list, np.ndarray):
        return np.flatnonzero(tf_list).tolist()
    return [i for i, x in enumerate(tf_list) if x]


if __name__ == "__main__":
    print(choose_packet_numbers(500, 123, RaptorDistribution(500)), False, None)
