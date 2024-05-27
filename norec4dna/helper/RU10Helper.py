#!/usr/bin/python
# -*- coding: latin-1 -*-
import typing
import numpy as np
from functools import lru_cache
from math import ceil, floor, pow, sqrt, log

from norec4dna.distributions import RaptorDistribution

int63 = int(pow(2, 63) - 1)
int31 = int(pow(2, 31) - 1)


def choose_packet_numbers(number_of_chunks: int, code_block_index: int, dist: RaptorDistribution,
                          systematic: bool = False, max_l: int = None) -> typing.List[int]:
    if systematic:
        d, a, b = systematic_ru10_triple_generator(number_of_chunks, code_block_index, dist)
    else:
        d, a, b = ru10_triple_generator(number_of_chunks, code_block_index, dist, max_l)
    if max_l is None:
        l, _, _ = intermediate_symbols(number_of_chunks, dist)
    else:
        l = max_l
    lprime: np.uint32 = np.uint32(dist.smallestPrimeGreaterOrEqual(l))

    if d > l:
        d = l
    l = np.uint32(l)
    indices: typing.List[int] = [0] * d
    while b >= l:
        b = (b + a) % lprime

    indices[0] = b

    for idx in range(1, d):
        b = (b + a) % lprime
        while b >= l:
            b = (b + a) % lprime
        indices[idx] = b
    return sorted(indices)


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


def ru10_triple_generator(k: int, x: int, dist: RaptorDistribution, max_l: typing.Optional[int] = None) \
        -> typing.Tuple[int, np.uint32, np.uint32]:
    if max_l is None:
        l, _, _ = intermediate_symbols(k, dist)
    else:
        l = max_l
    lprime = dist.smallestPrimeGreaterOrEqual(l)
    rng: np.random = np.random
    rng.seed(x)
    v = np.uint32(r_int63(rng) % 1048576)
    a = np.uint32(1 + (r_int63(rng) % (lprime - 1)))
    b = np.uint32(r_int63(rng) % lprime)
    d = dist.deg(v)
    return d, a, b


def systematic_ru10_triple_generator(k: int, x: int, dist: RaptorDistribution) -> typing.Tuple[int, int, int]:
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


def r_int63(rng: np.random) -> int:
    return rng.randint(0, int31)


def from_true_false_list_old(tf_list: typing.List[bool]) -> typing.List[int]:
    return [i for i, x in enumerate(tf_list) if x]


def from_true_false_list(tf_list: typing.List[bool]) -> typing.List[int]:
    return np.nonzero(tf_list)[0].tolist()


if __name__ == '__main__':
    # generate a list of 1000 True/False values:
    true_false_list = [np.random.choice([True, False]) for _ in range(1000)]
    print(res1 := from_true_false_list(true_false_list))
    print(res2 := from_true_false_list_old(true_false_list))
    print(res1 == res2)
    print(choose_packet_numbers(500, 123, RaptorDistribution.RaptorDistribution(500)), False, None)
