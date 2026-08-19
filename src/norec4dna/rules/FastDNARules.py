import copy
import math
from collections.abc import Callable, Iterable
from functools import partial
from typing import List, Optional, Set, Tuple

from ..helper.fallback_code import r_region, small_r_region
from .RuleParser import (
    charCountBiggerEqualThanX,
    gc_content,
    length,
    longestSequenceOfChar,
    microsatellite,
    strContainsIllegalChars,
    strContainsSub,
    strContainsSubRegex,
)

try:
    from cdnarules import repeatRegion as _repeat_region
    from cdnarules import smallRepeatRegion as _small_repeat_region
except Exception:
    print("C Module failed to load, falling back to slow mode")
    _repeat_region = None
    _small_repeat_region = None


def rRegion(data: str, repeat_length: int = 20) -> float:
    if _repeat_region is None:
        return r_region(data, repeat_length)
    return float(_repeat_region(data, repeat_length))


def smallrRegion(data: str, repeat_length: int = 9) -> float:
    if _small_repeat_region is None:
        return small_r_region(data, repeat_length)
    return float(_small_repeat_region(data, repeat_length))


_undes_motifs = [
    ("CTCGTAGACTGCGTACCA", 1.01),
    ("GACGATGAGTCCTGAGTA", 1.01),
    ("CTGTCTCTTATACACATCT", 1.01),
    ("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", 1.01),
    ("GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", 1.01),
]

undes_motifs = [
    # Lox sites.
    ("ATAACTTCGTATAGCATACATTATACGAAGTTAT", 1.01),
    ("ATAACTTCGTATAGCATACATTATACGAACGGTA", 1.01),
    ("TACCGTTCGTATAGCATACATTATACGAAGTTAT", 1.01),
    ("TACCGTTCGTATAGCATACATTATACGAACGGTA", 1.01),
    ("TACCGTTCGTATATGGTATTATATACGAAGTTAT", 1.01),
    ("TACCGTTCGTATATTCTATCTTATACGAAGTTAT", 1.01),
    ("TACCGTTCGTATAGGATACTTTATACGAAGTTAT", 1.01),
    ("TACCGTTCGTATATACTATACTATACGAAGTTAT", 1.01),
    ("TACCGTTCGTATACTATAGCCTATACGAAGTTAT", 1.01),
    ("ATAACTTCGTATATGGTATTATATACGAACGGTA", 1.01),
    ("ATAACTTCGTATAGTATACCTTATACGAAGTTAT", 1.01),
    # Twister Adapters:
    ("GAAGTGCCATTCCGCCTGACCT", 1.0),  # Twister 5' Adapter
    ("AGGCTAGGTGGAGGCTCAGTG", 1.0),  # Twister 3' Adapter
]


def gc_error_calculation(gc_percentage):
    return (
        100
        + (175 * gc_percentage) / 6
        - (121 * gc_percentage**2) / 72
        + (gc_percentage**3) / 36
        - (gc_percentage**4) / 7200
    ) / 100


def fs_gc_error_calculation(gc_percentage):
    return 1.0 if gc_percentage > 60 or gc_percentage < 40 else 0.0


def ts_gc_error_calculation(gc_percentage):
    return 1.0 if gc_percentage > 70 or gc_percentage < 30 else 0.0


def gc_strict_calculation(gc_percentage):
    return (
        100
        + 49970.8 * gc_percentage
        - 2582 * gc_percentage**2
        + 41.6458 * gc_percentage**3
        - 0.208229 * gc_percentage**4
    ) / 100


def strict_homopolymers():
    return [0.0, 0.0, 0.2, 0.5, 0.8, 1.0]


# x = |Homopolymer| ; x > 3 -> 100% ; x <= 3 -> 0%
def three_strict_homopolymers():
    return [0.0, 0.0, 0.0, 0.0, 1.0]


def four_strict_homopolymers():
    return [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]


def lax_increasing_homopolymers():
    return [0.0, 0.0, 0.01, 0.05, 0.4, 0.7, 0.9, 1.0]


def lax_homopolymers():
    res = [0.0, 0.0]
    for x in strict_homopolymers():
        res.append(x)
    return res


class FastDNARules:
    def __init__(self, active_rules=None, extra_motifs=None):
        self.nineteen_mers: Set[str] = set()
        self.tmp_nineteen_mers: Set[str] = set()
        # Extra forbidden motifs supplied at construction time or via add_forbidden_sequences().
        # Each entry is a (dna_sequence, error_probability) tuple; error_probability >= 1.01
        # guarantees the packet is rejected (same convention as the built-in undesired motifs).
        self.extra_motifs: List[Tuple[str, float]] = list(extra_motifs) if extra_motifs else []
        if active_rules is None:
            self.active_rules: List[Callable[[str], float]] = [
                # FastDNARules.a_permutation,
                # FastDNARules.t_permutation,
                # FastDNARules.c_permutation,
                # FastDNARules.g_permutation,
                # FastDNARules.dinucleotid_runs,
                # FastDNARules.homopolymers,
                partial(FastDNARules.homopolymers, probs=four_strict_homopolymers()),
                # FastDNARules.overall_gc_content,
                # To change the GC error function:
                partial(FastDNARules.overall_gc_content, calc_func=fs_gc_error_calculation),
                # FastDNARules.windowed_gc_content,
                partial(FastDNARules.windowed_gc_content, calc_func=ts_gc_error_calculation),
                #  FastDNARules.long_strands,
                #  FastDNARules.illegal_symbols,
                # FastDNARules.trinucleotid_runs,
                # FastDNARules.random_permutations,
                self._motif_search,  # bound method so extra_motifs are included
                # FastDNARules.motif_regex_search,
                # FastDNARules.repeatRegion,
                # FastDNARules.smallRepeatRegion,
                # self.check_and_add_mers
            ]
        else:
            self.active_rules = active_rules

    def add_forbidden_sequences(self, sequences: Iterable[str], error_prob: float = 1.01) -> None:
        """Register additional DNA sequences that must not appear in any encoded packet.

        Sequences are added to the instance-level extra_motifs list and are checked
        by ``_motif_search`` alongside the built-in undesired-motif catalogue.
        An ``error_prob`` of 1.01 (the default) guarantees the packet is rejected
        because ``generate_new_packets`` only accepts packets with ``error_prob < 1.0``.

        Args:
            sequences: Iterable of DNA strings (uppercase A/C/G/T) to forbid.
            error_prob: Error probability assigned to each sequence (default 1.01).
        """
        for seq in sequences:
            self.extra_motifs.append((seq, error_prob))

    def _motif_search(self, data: str) -> float:
        """Instance wrapper around ``motif_search`` that also checks ``extra_motifs``."""
        drop = FastDNARules.motif_search(data)
        if drop >= 1.0 or not self.extra_motifs:
            return min(1.0, drop)
        for motif, prob in self.extra_motifs:
            if motif in data:
                drop += prob
                if drop >= 1.0:
                    return 1.0
        return drop

    def check_and_add_mers(self, data: str, length: int = 19) -> float:
        chunks = [data[i : i + length] for i in range(0, len(data), length)]
        self.tmp_nineteen_mers.clear()
        res = 0.0
        for chunk in chunks:
            if chunk in self.nineteen_mers or chunk in self.tmp_nineteen_mers:
                res = 0.7
                break
            else:
                self.tmp_nineteen_mers.add(chunk)
        for chunk in chunks:
            self.nineteen_mers.add(chunk)
        return res

    @staticmethod
    def repeatRegion(data, repeat_length=20):
        return rRegion(data, repeat_length)

    @staticmethod
    def smallRepeatRegion(data, repeat_length=9):
        """
        :param data:
        :param repeat_length:
        :return:
        """
        return smallrRegion(data, repeat_length)

    # @staticmethod
    def apply_all_rules(self, packet):
        """
        Apply all rules to the DNA sequence of the given packet and returns the summed error probability.
        :param packet:
        :return:
        """
        try:
            dna_data = packet.get_dna_struct(True)
        except Exception:
            dna_data = packet
        res_arr = [x(dna_data) for x in self.active_rules]
        return sum(res_arr)

    # @staticmethod
    def apply_all_rules_with_data(self, packet):
        """

        :param packet:
        :return:
        """
        try:
            dna_data = packet.get_dna_struct(True)
        except AttributeError:
            dna_data = packet
        res_arr = [x(dna_data) for x in self.active_rules]
        return sum(res_arr), res_arr, packet

    _DEFAULT_HOMOPOLYMER_PROBS: Tuple[float, ...] = (0.0, 0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0)

    @staticmethod
    def homopolymers(data: str, probs: Optional[List[float]] = None) -> float:
        """Calculates the dropchance based on homopolymer run length."""
        p = probs if probs is not None else FastDNARules._DEFAULT_HOMOPOLYMER_PROBS
        homopolymer_length = longestSequenceOfChar(data, "*")[1]
        return min(1.0, p[homopolymer_length] if homopolymer_length < len(p) else 1.0)

    @staticmethod
    def dinucleotid_runs(data):
        """
        Calculates the dropchance based on the chance of dinucleotid microsatellites to mutate. Microsatellites are
        repeats of two to six units that have a higher chance to mutate than regular polymers.
        The dropchances are additive in this case (@shouldDrop()).
        Function y = 0.00002x^2 - 0.0001x
        :param data: The DNA sequence to check for microsatellites with length 2.
        :return: The dropchance based on the occurence of microsatellites with length 2.
        """
        _len = microsatellite(data, 2)[0]
        return max(min(1.0, 0.00002 * _len**2 - 0.0001 * _len), 0.0)

    @staticmethod
    def trinucleotid_runs(data):
        """
        Calculates the dropchance based on the chance of trinucleotid microsatellites to mutate. Microsatellites are
        repeats of two to six units that have a higher chance to mutate than regular polymers.
        The dropchances are additive in this case (@shouldDrop()).
        Function y = 0.00002x^2 - 0.0001x
        :param data: The DNA sequence to check for microsatellites with length 3.
        :return: The dropchance based on the occurence of microsatellites with length 3.
        """
        _len = microsatellite(data, 3)[0]
        return max(min(1.0, 0.00002 * _len**2 - 0.0001 * _len), 0.0)

    @staticmethod
    def long_strands(data):
        """
        Calculates the dropchance based on the length of the sequence since longer strands have a higher chance to mutate.
        Function y = 0.00144363 * exp(0.0140971 * x), x > 117
        :param data: The sequence to check for the strand length.
        :return: The dropchance based on the strand length.
        """
        _len = length(data)
        if 117 <= _len:
            dropchance = 0.00144363 * math.exp(0.0140971 * _len)
        else:
            dropchance = 0.0
        return min(dropchance, 0.2)

    @staticmethod
    def random_permutations(data):
        """
        Calculates the dropchance for simulated random mutations in the DNA-Data.
        :param data: The sequence to simulate random mutations for.
        :return: The dropchance based on random mutations.
        """
        return 0.02

    @staticmethod
    def illegal_symbols(data):
        """
        Checks the DNA data for illegal symbols and returns a dropchance of 1.0 if the sequence contains them.
        :param data: The sequence to check for illegal symbols.
        :return: The dropchance based on the occurence of illegal symbols (0.0 or 1.0).
        """
        if strContainsIllegalChars(data, "ACGT"):
            return 1.0
        return 0.0

    @staticmethod
    def ch_cg_permutation(data, ch):
        """
        Calculates the dropchance based on the number of occurrences of a nucleotide, G and C in this case since their
        chance to mutate is generally higher.
        These errors are additiv.
        :param data: The sequence to check for the number of occurences of a nucleotide.
        :param ch: The nucleotide to check for.
        :return: The dropchance based on the number of occurences of the given nucleotide.
        """
        occ = charCountBiggerEqualThanX(data, ch)
        dropchance = math.floor(occ / 20) * 0.002
        return min(dropchance, 0.02)

    @staticmethod
    def ch_at_permutation(data, ch):
        """
        Calculates the dropchance based on the number of occurrences of a nucleotide, A and T in this case since their
        chance to mutate is generally lower.
        These errors are additiv.
        :param data: The sequence to check for the number of occurences of a nucleotide.
        :param ch: The nucleotide to check for.
        :return: The dropchance based on the number of occurences of the given nucleotide.
        """
        occ = charCountBiggerEqualThanX(data, ch)
        dropchance = math.floor(occ / 20) * 0.001
        return min(dropchance, 0.01)

    @staticmethod
    def a_permutation(data):
        """
        Calls @ch_at_permutation for 'A'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_at_permutation(data, "A")

    @staticmethod
    def t_permutation(data):
        """
        Calls @ch_at_permutation for 'T'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_at_permutation(data, "T")

    @staticmethod
    def g_permutation(data):
        """
        Calls @ch_cg_permutation for 'G'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_cg_permutation(data, "G")

    @staticmethod
    def c_permutation(data):
        """
        Calls @ch_cg_permutation for 'C'
        :param data: Sequence
        :return: Dropchance
        """
        return FastDNARules.ch_cg_permutation(data, "C")

    @staticmethod
    def overall_gc_content(data, calc_func=gc_error_calculation):
        """
        Calculates the dropchance of the sequence based on the GC-content since GC-rich DNA is more stable.
        :param calc_func: function to calculate the error-prob based on the gc-content %
        :param data: The sequence to check for the GC-content.
        :return: The Dropchance based on the GC-content.
        """
        return max(0.0, min(1.0, calc_func(gc_content(data))))

    @staticmethod
    def windowed_gc_content(data, window_size=50, calc_func=gc_error_calculation):
        chunks = [
            FastDNARules.overall_gc_content(data[i : i + window_size], calc_func=calc_func)
            for i in range(0, len(data), window_size)
        ]
        return min(1.0, max(chunks))

    _UNDES_MOTIFS_CATALOGUE: Tuple[Tuple[str, float], ...] = (
        # Promoter recognition motif (Euk).
        ("TATAAA", 1.01),
        # Promoter recognition motifs (Prok).
        ("TTGACA", 1.05),
        ("TGTATAATG", 1.05),
        # Polyadenylation signals (Euk).
        ("AATAAA", 1.01),
        ("TTGTGTGTTG", 1.01),
        # Lox sites.
        ("ATAACTTCGTATAGCATACATTATACGAAGTTAT", 1.01),
        ("ATAACTTCGTATAGCATACATTATACGAACGGTA", 1.01),
        ("TACCGTTCGTATAGCATACATTATACGAAGTTAT", 1.01),
        ("TACCGTTCGTATAGCATACATTATACGAACGGTA", 1.01),
        ("TACCGTTCGTATATGGTATTATATACGAAGTTAT", 1.01),
        ("TACCGTTCGTATATTCTATCTTATACGAAGTTAT", 1.01),
        ("TACCGTTCGTATAGGATACTTTATACGAAGTTAT", 1.01),
        ("TACCGTTCGTATATACTATACTATACGAAGTTAT", 1.01),
        ("TACCGTTCGTATACTATAGCCTATACGAAGTTAT", 1.01),
        ("ATAACTTCGTATATGGTATTATATACGAACGGTA", 1.01),
        ("ATAACTTCGTATAGTATACCTTATACGAAGTTAT", 1.01),
        # Lox site spacers not covered by the Lox sites.
        ("AGGTATGC", 1.01),
        ("TTGTATGG", 1.01),
        ("GGATAGTA", 1.01),
        ("GTGTATTT", 1.01),
        ("GGTTACGG", 1.01),
        ("TTTTAGGT", 1.01),
        ("GTACACAT", 1.01),
        # Restriction enzyme recognition motifs.
        # BpiI
        ("GAAGAC", 1.01),
        # inverse BpiI
        ("CTTCTG", 1.01),
        # BsaI
        ("GGTCTC", 1.01),
        # inverse BsaI
        ("CCAGAG", 1.01),
        ("CGTCTC", 1.01),
        ("GCGATG", 1.01),
        ("GCTCTTC", 1.01),
        # Oligo Adapters.
        ("CTCGTAGACTGCGTACCA", 1.01),
        ("GACGATGAGTCCTGAGTA", 1.01),
        # 5' extensions.
        ("GGTTCCACGTAAGCTTCC", 1.01),
        ("GCGATTACCCTGTACACC", 1.01),
        ("GCCAGTACATCAATTGCC", 1.01),
        # Twister Adapters:
        ("GAAGTGCCATTCCGCCTGACCT", 1.0),  # Twister 5' Adapter
        ("AGGCTAGGTGGAGGCTCAGTG", 1.0),  # Twister 3' Adapter
    )

    @staticmethod
    def motif_search(data: str) -> float:
        """Searches for undesired motifs with given error probabilities."""
        dropchance = 0.0
        for motif_seq, prob in FastDNARules._UNDES_MOTIFS_CATALOGUE:
            if motif_seq in data:
                dropchance += prob
                if dropchance >= 1.0:
                    return 1.0
        return min(1.0, dropchance)

    @staticmethod
    def motif_regex_search(data):
        """
        Searches for undesired motifs with given error probabilities.
        :param data:
        :return:
        """
        undes_motifs = [
            # Promoter recognition motif (Euk).
            ("CANYYY", 0.01),
            ("ANCCAATCA", 0.01),
            ("KGGGCGGRRY", 0.01),
            ("KRGGCGKRRY", 0.01),
            # Promoter recognition motifs (Prok).
            ("AAAWWTWTTTTNNNAAA", 0.05),
            # Ribosomal binding site (Euk).
            ("RCCACCATGG", 0.05),
            # Ribosomal binding site (Prok).
            ("AGGAGGACAGCTAUG", 0.05),
            # Lox sites.
            ("ATAACTTCGTATAGTAYACATTATACGAAGTTAT", 0.01),
        ]
        dropchance = 0.0
        for motif in undes_motifs:
            if strContainsSubRegex(data, motif[0]):
                dropchance += motif[1]
        return min(1.0, dropchance)

    @staticmethod
    def add_complementary(input_lst):
        def revert(in_chr):
            return {"G": "C", "C": "G", "T": "A", "A": "T"}[in_chr]

        tmp = copy.deepcopy(input_lst)

        for x in input_lst:
            tmp.append("".join([revert(chr) for chr in x]))
        return tmp

    @staticmethod
    def add_reverse_complementary(input_lst):
        def revert(in_chr):
            return {"G": "C", "C": "G", "T": "A", "A": "T"}[in_chr]

        tmp = copy.deepcopy(input_lst)

        for x in input_lst:
            tmp.append(("".join([revert(chr[0]) for chr in x[0][::-1]]), x[1]))
        return tmp

    @staticmethod
    def add_reverse(input_lst):
        tmp = copy.deepcopy(input_lst)
        for x in input_lst:
            tmp.append(x[::-1])
        return tmp

    @staticmethod
    def simple_motif_search(data):
        return (
            1.0
            if any(
                x in data
                for x in FastDNARules.add_complementary(
                    [
                        "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
                        "ATAACTTCGTATAGCATACATTATACGAACGGTA",
                        "TACCGTTCGTATAGCATACATTATACGAAGTTAT",
                        "TACCGTTCGTATAGCATACATTATACGAACGGTA",
                        "TACCGTTCGTATATGGTATTATATACGAAGTTAT",
                        "TACCGTTCGTATATTCTATCTTATACGAAGTTAT",
                        "TACCGTTCGTATAGGATACTTTATACGAAGTTAT",
                        "TACCGTTCGTATATACTATACTATACGAAGTTAT",
                        "TACCGTTCGTATACTATAGCCTATACGAAGTTAT",
                        "ATAACTTCGTATATGGTATTATATACGAACGGTA",
                        "ATAACTTCGTATAGTATACCTTATACGAAGTTAT",
                        "ATAACTTCGTATAGTATACATTATACGAAGTTAT",
                        "ATAACTTCGTATAGTACACATTATACGAAGTTAT",
                        "GCATACAT",
                        "TGGTATTA",
                        "TTCTATCT",
                        "GGATACTT",
                        "TACTATAC",
                        "CTATAGCC",
                        "AGGTATGC",
                        "TTGTATGG",
                        "GGATAGTA",
                        "GTGTATTT",
                        "GGTTACGG",
                        "TTTTAGGT",
                        "GTATACCT",
                        "GTACACAT",
                        "GAAGAC",
                        "CTTCTG",
                        "GGTCTC",
                        "CCAGAG",
                    ]
                )
            )
            else 0.0
        )


if __name__ == "__main__":
    x = FastDNARules()
    x.motif_search("TTGACA")
    # print(x.add_reverse_complementary([("CTCGTAGACTGCGTACCA", 1.01)]))
    # print(x.check_and_add_mers("AAAAGAGAGAGAGAGAGAGCCCCCCCCCCCCCCCCCCCAACAGAGAGAGAGAGAGAG", 19))
    # print(x.check_and_add_mers("CCCCCCCCCCCCCCCCCCC", 19))
    print(
        x.homopolymers(
            "GGTCTCGCAAGTTACGTGTCTATTTAGCGCGGCATATCACAGCGGCGGTACGCATAACAGTTTACA"
            "GGGAAAGTAGATCATCAGGCGTGGCTAGGGAGCGCGTGTCCTCATTTGTTGAGGAGACGCTAAAGC"
            "ACCCGGGTAGTAAATATCTGAACATGGGGGGG",
            probs=three_strict_homopolymers(),
        )
    )
    print(
        x.overall_gc_content(
            "GGTCTCGCAAGTTACGTGTCTATTTAGCGCGGCATATCACAGCGGCGGTACGCATAACAGTTTACA"
            "GGGAAAGTAGATCATCAGGCGTGGCTAGGGAGCGCGTGTCCTCATTTGTTGAGGAGACGCTAAAGC"
            "ACCCGGGTAGTAAATATCTGAACATGGGGGGG"
        )
    )
    # print(x.motif_search("ATGGTACGCAAGTCTACGAG"))
