from __future__ import annotations

import struct
import typing
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
from reedsolo import ReedSolomonError

from .distributions.Distribution import Distribution
from .distributions.RaptorDistribution import RaptorDistribution
from .helper import logical_xor
from .helper.helper import xor_with_seed
from .helper.RU10Helper import choose_packet_numbers, from_true_false_list
from .RU10Packet import RU10Packet

BoolArray = np.ndarray[Any, np.dtype[np.bool_]]


class RU10Shared:
    """Methods shared verbatim by :class:`RU10Decoder` and :class:`RU10BPDecoder`.

    Both RU10 decoders implement the same packet-format parsing and Raptor
    aux-block bookkeeping, but they differ in the solver (Gaussian elimination
    vs. belief propagation).  Keeping the shared logic here (instead of
    delegating ``RU10Decoder.method(self, ...)`` from the BP decoder) gives both
    classes a single, correctly-typed implementation.
    """

    number_of_chunks: int
    s: int
    h: int
    use_method: bool
    ldpcANDhalf: Dict[int, Any]
    _distribution: Optional[Distribution]

    @property
    def distribution(self) -> Optional[Distribution]:
        return self._distribution

    @distribution.setter
    def distribution(self, val: Optional[Distribution]) -> None:
        self._distribution = val
    correct: int
    corrupt: int
    id_spacing: int
    config_map: Optional[Any]
    error_correction: Any
    xor_by_seed: bool
    static_number_of_chunks: Optional[int]

    # Subclasses provide these (they differ between GEPP and BP):
    def _update_packet_header(
        self,
        len_data: Tuple[Any, ...],
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> int:  # pragma: no cover - abstract
        raise NotImplementedError

    def _ensure_distribution_ready(self) -> None:  # pragma: no cover - abstract
        raise NotImplementedError

    # ── packet parsing helpers ─────────────────────────────────────────────

    def _packet_method_chunks(self, data: bytes) -> Tuple[bytes, List[int]]:
        if not self.use_method:
            return data, []
        method_data = bin(data[-1])[2:].rjust(8, "0")
        return data[:-1], self._chunk_list_from_method_data(method_data)

    def _chunk_list_from_method_data(self, method_data: str) -> List[int]:
        if method_data.startswith("00"):
            return [chunk for chunk in range(0, self.number_of_chunks + 1) if chunk % 2 == 0]
        if method_data.startswith("01"):
            return [chunk for chunk in range(0, self.number_of_chunks + 1) if chunk % 2 != 0]
        if method_data.startswith("10"):
            return self._window_chunk_list(method_data, 30)
        if method_data.startswith("11"):
            return self._window_chunk_list(method_data, 40)
        raise RuntimeError("Not a valid start:", method_data)

    def _window_chunk_list(self, method_data: str, window_size: int) -> List[int]:
        window = int(method_data[2:], 2)
        start = window * (window_size - 10)
        return [
            chunk for chunk in range(start, start + window_size) if chunk <= self.number_of_chunks
        ]

    def _choose_used_packets(
        self, unxored_id: int, chunk_lst: Optional[List[int]]
    ) -> typing.Sequence[int]:
        distribution = typing.cast(RaptorDistribution, self.distribution)
        if self.use_method and chunk_lst is not None:
            numbers = choose_packet_numbers(
                len(chunk_lst), unxored_id, distribution, systematic=False, max_l=len(chunk_lst)
            )
            return [chunk_lst[i] for i in numbers]
        return choose_packet_numbers(
            self.number_of_chunks,
            unxored_id,
            distribution,
            systematic=False,
        )

    # ── seed spacing / appended version fields ─────────────────────────────

    def parse_raw_packet(
        self,
        packet_input: bytes,
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "L",
        packet_len_format: str = "I",
        id_len_format: str = "L",
        dna_str: Optional[str] = None,
    ) -> Optional[Union[RU10Packet, typing.Literal["CORRUPT"]]]:
        """
        Creates a RU10 packet from a raw given packet. Also checks if the packet is corrupted. If any method was used to
        create packets from specific chunks, set self.use_method = True. This will treat the last byte of the raw packet
        data as the byte that contains the information about the used method ("even", "odd", "window_30 + window" or
        "window_40 + window". See RU10Encoder.create_new_packet_from_chunks for further information.
        :param dna_str: dna string for reference in the packet object
        :param packet_input: A raw packet
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: RU10Packet or an error message
        """
        struct_str = "<" + number_of_chunks_len_format + id_len_format
        struct_len = struct.calcsize(struct_str)
        try:
            packet = self.error_correction(packet_input)
        except (AssertionError, ReedSolomonError, ValueError):
            self.corrupt += 1
            return "CORRUPT"
        header = typing.cast(bytes, packet)[:struct_len]
        data = typing.cast(bytes, packet)[struct_len:]
        data, chunk_lst = self._packet_method_chunks(data)
        len_data = struct.unpack(struct_str, header)
        unxored_id = self._update_packet_header(
            len_data, number_of_chunks_len_format, id_len_format
        )
        if self.xor_by_seed:
            data = xor_with_seed(data, unxored_id)
        self._ensure_distribution_ready()
        used_packets = self._choose_used_packets(unxored_id, chunk_lst)
        res = RU10Packet(
            data,
            used_packets,
            self.number_of_chunks,
            unxored_id,
            read_only=True,
            packet_len_format=packet_len_format,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            id_len_format=id_len_format,
            save_number_of_chunks_in_packet=self.static_number_of_chunks is None,
        )
        res.dna_data = dna_str
        res.packed_used_packets = packet  # without error correction
        res.packed_struct = packet_input  # with error correction
        return res

    def revert_seed_spacing(self, dna_str: str, id_len_format: str) -> str:
        struct_len = struct.calcsize(id_len_format) * 4
        if self.id_spacing > 0 and struct_len > 0:
            res = ""
            input_str = list(dna_str)
            i = 0
            while len(res) < struct_len:
                res += input_str[i]
                input_str[i] = " "
                i += self.id_spacing + 1
            input_str = "".join(input_str)
            input_str = input_str.replace(" ", "")
            res += input_str
            return res
        return dna_str

    def strip_appended_version_fields(self, dna_str: str) -> str:
        """Strip appended version fields from DNA string if enabled in config."""
        if self.config_map is None:
            return dna_str

        append_version_fields = self.config_map.get("append_version_fields", False)
        if isinstance(append_version_fields, str):
            append_version_fields = append_version_fields.lower() in ("true", "1", "yes", "on")

        if not append_version_fields:
            return dna_str

        try:
            version_bits = int(self.config_map.get("version_bits", 8))
            chunk_idx_bits = int(self.config_map.get("chunk_idx_bits", 8))
            algo_id_bits = int(self.config_map.get("algo_id_bits", 3))
            iter_pos_bits = int(self.config_map.get("iter_pos_bits", 5))
            append_version_include_magic = self.config_map.get("append_version_include_magic", True)
            if isinstance(append_version_include_magic, str):
                append_version_include_magic = append_version_include_magic.lower() in (
                    "true",
                    "1",
                    "yes",
                    "on",
                )

            total_bits = version_bits + chunk_idx_bits + algo_id_bits + iter_pos_bits
            total_bytes = (total_bits + 7) // 8
            total_bases = total_bytes * 4

            if append_version_include_magic:
                magic_string = self.config_map.get("magic_string", "GAGCCAGTGAGTCGTA")
                total_bases += len(magic_string)

            if total_bases > 0 and len(dna_str) >= total_bases:
                return dna_str[:-total_bases]
        except Exception:
            pass
        return dna_str

    # ── aux-block index helpers (used by removeAndXorAuxPackets_from_indices) ──

    def _normalize_packet_index_set(
        self, packet_indices: Union[Set[int], List[int], np.ndarray]
    ) -> Set[int]:
        if isinstance(packet_indices, np.ndarray) and packet_indices.dtype == bool:
            return set(from_true_false_list(packet_indices.tolist()))
        return set(packet_indices)

    def _split_packet_indices(
        self, packet_set: Set[int]
    ) -> Tuple[Set[int], Set[int], Set[int]]:
        systematic_indices: Set[int] = set()
        ldpc_indices: Set[int] = set()
        half_indices: Set[int] = set()
        for idx in packet_set:
            if 0 <= idx < self.number_of_chunks:
                systematic_indices.add(idx)
            elif self.number_of_chunks <= idx < self.number_of_chunks + self.s:
                ldpc_indices.add(idx - self.number_of_chunks)
            elif self.number_of_chunks + self.s <= idx < self.number_of_chunks + self.s + self.h:
                half_indices.add(idx - self.number_of_chunks - self.s)
        return systematic_indices, ldpc_indices, half_indices

    def _systematic_result(self, systematic_indices: Set[int]) -> np.ndarray:
        result = np.zeros(self.number_of_chunks, dtype=bool)
        for chunk_idx in systematic_indices:
            result[chunk_idx] = True
        return result

    def _pad_bool_arrays(self, arrays: List[List[bool]], target_size: int) -> List[List[bool]]:
        padded_arrays: List[List[bool]] = []
        for arr in arrays:
            if len(arr) < target_size:
                padded = [False] * target_size
                padded[: len(arr)] = arr
                padded_arrays.append(padded)
            else:
                padded_arrays.append(arr)
        return padded_arrays

    def _half_packet_arrays(self, half_indices: Set[int]) -> Tuple[List[List[bool]], int]:
        half_arrays: List[List[bool]] = []
        half_arr_size = 0
        for half_idx in half_indices:
            ldpc_half_idx = self.s + half_idx
            if ldpc_half_idx in self.ldpcANDhalf:
                half_arr = self.ldpcANDhalf[ldpc_half_idx].get_bool_array_used_and_ldpc_packets()
                half_arr_size = max(half_arr_size, len(half_arr))
                half_arrays.append(half_arr)
        return half_arrays, half_arr_size

    def _systematic_and_ldpc_array(
        self, systematic_indices: Set[int], ldpc_indices: Set[int], half_arr_size: int
    ) -> np.ndarray:
        systematic_and_ldpc = np.zeros(half_arr_size, dtype=bool)
        for chunk_idx in systematic_indices:
            if chunk_idx < half_arr_size:
                systematic_and_ldpc[chunk_idx] = True
        for ldpc_idx in ldpc_indices:
            if self.number_of_chunks + ldpc_idx < half_arr_size:
                systematic_and_ldpc[self.number_of_chunks + ldpc_idx] = True
        return systematic_and_ldpc

    def _extract_half_result(self, combined: np.ndarray) -> Tuple[np.ndarray, Set[int]]:
        result = np.zeros(self.number_of_chunks, dtype=bool)
        new_ldpc_from_half: Set[int] = set()
        for idx in from_true_false_list(combined):
            if 0 <= idx < self.number_of_chunks:
                result[idx] = True
            elif self.number_of_chunks <= idx < self.number_of_chunks + self.s:
                new_ldpc_from_half.add(idx - self.number_of_chunks)
        return result, new_ldpc_from_half

    def _process_half_indices(
        self, half_indices: Set[int], systematic_indices: Set[int], ldpc_indices: Set[int]
    ) -> Tuple[np.ndarray, Set[int]]:
        result = self._systematic_result(systematic_indices)
        new_ldpc_from_half: Set[int] = set()
        if not half_indices:
            return result, new_ldpc_from_half
        half_arrays, half_arr_size = self._half_packet_arrays(half_indices)
        if not half_arrays:
            return result, new_ldpc_from_half
        half_result = logical_xor(self._pad_bool_arrays(half_arrays, half_arr_size))
        systematic_and_ldpc = self._systematic_and_ldpc_array(
            systematic_indices, ldpc_indices, half_arr_size
        )
        combined = logical_xor([half_result.tolist(), systematic_and_ldpc.tolist()])
        return self._extract_half_result(combined)

    def _process_ldpc_indices(self, ldpc_indices: Set[int], result: np.ndarray) -> np.ndarray:
        aux_list: List[List[bool]] = []
        for ldpc_idx in ldpc_indices:
            if ldpc_idx in self.ldpcANDhalf:
                aux_arr = self.ldpcANDhalf[ldpc_idx].get_bool_array_used_packets()
                if aux_arr is not None and len(aux_arr) < self.number_of_chunks:
                    padded = [False] * self.number_of_chunks
                    padded[: len(aux_arr)] = aux_arr
                    aux_list.append(padded)
                elif aux_arr is not None:
                    aux_list.append(aux_arr)
        if not aux_list:
            return result
        return logical_xor([logical_xor(aux_list).tolist(), result.tolist()])

    def removeAndXorAuxPackets_from_indices(
        self, packet_indices: Union[Set[int], List[int]]
    ) -> np.ndarray:
        """
        Removes auxpackets (LDPC and Half) from a list/set of packet indices to get the chunk composition.

        This is the equivalent of removeAndXorAuxPackets but operates on packet indices directly
        without requiring a RU10Packet object.

        Packet index interpretation:
        - Indices 0 to number_of_chunks-1: systematic packets (each contains exactly one chunk)
        - Indices number_of_chunks to number_of_chunks+s-1: LDPC packets (indexed in ldpcANDhalf as 0 to s-1)
        - Indices number_of_chunks+s to number_of_chunks+s+h-1: Half packets (indexed in ldpcANDhalf as s to s+h-1)

        Returns:
            Boolean numpy array where index i=True means chunk i is in the result after aux removal
        """
        packet_set = self._normalize_packet_index_set(packet_indices)
        if not packet_set:
            return np.zeros(self.number_of_chunks, dtype=bool)
        systematic_indices, ldpc_indices, half_indices = self._split_packet_indices(packet_set)
        result, new_ldpc_from_half = self._process_half_indices(
            half_indices, systematic_indices, ldpc_indices
        )
        ldpc_to_process = new_ldpc_from_half if half_indices else ldpc_indices
        return self._process_ldpc_indices(ldpc_to_process, result)
