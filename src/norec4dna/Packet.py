from __future__ import annotations

import struct
import typing
from typing import Any, List, Optional, Set, Union

import numpy as np

from .ErrorCorrection import ErrorCorrectionCallable, nocode
from .helper import xor_mask, xor_numpy
from .helper.bin2Quaternary import quads2dna, string2QUATS

UInt8Array = np.ndarray[Any, np.dtype[np.uint8]]
BoolArray = np.ndarray[Any, np.dtype[np.bool_]]


def interleave_spacing(input_str: str, spacing: int, spacing_length: int) -> str:
    if spacing <= 0 or spacing_length <= 0:
        return input_str

    left, right = input_str[:spacing_length], input_str[spacing_length:]

    interleaved_str = ""
    i = 0
    j = 0

    while i < len(left) and j < len(right):
        if j % spacing == 0:
            interleaved_str += left[i]
            i += 1
        interleaved_str += right[j]
        j += 1

    # If there are remaining characters in left or right
    interleaved_str += left[i:] + right[j:]

    return interleaved_str


def deinterleave_spacing(interleaved_str: str, spacing: int, spacing_length: int) -> str:
    if spacing <= 0 or spacing_length <= 0:
        return interleaved_str

    left_chars: List[str] = []
    right_chars: List[str] = []

    i = 0
    left_count = 0
    right_count = 0

    while i < len(interleaved_str):
        # Alle `spacing` Zeichen von rechts wird ein Zeichen von links eingefügt
        if right_count % spacing == 0 and left_count < spacing_length:
            left_chars.append(interleaved_str[i])
            left_count += 1
            i += 1

        if i < len(interleaved_str):
            right_chars.append(interleaved_str[i])
            right_count += 1
            i += 1

    return "".join(left_chars) + "".join(right_chars)


class Packet:
    def __init__(
        self,
        data: Union[bytes, UInt8Array],
        used_packets: Set[int],
        total_number_of_chunks: int,
        read_only: bool = False,
        seed: int = 0,
        implicit_mode: bool = True,
        error_correction: ErrorCorrectionCallable = nocode,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        used_packets_len_format: str = "I",
        id_len_format: str = "I",
        save_number_of_chunks_in_packet: bool = True,
        prepend: str = "",
        append: str = "",
    ):
        self.data: Union[bytes, UInt8Array] = data
        self.error_correction: ErrorCorrectionCallable = error_correction
        self.total_number_of_chunks: int = total_number_of_chunks
        self.used_packets: Set[int] = set()
        self.set_used_packets(used_packets)
        self.internal_hash: Optional[int] = None
        self.degree: int = 0  # just a stub
        self.update_degree()
        self.dna_data: Optional[str] = None
        self.id: int = seed
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.used_packets_len_format: str = used_packets_len_format
        self.id_len_format: str = id_len_format
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.error_prob: Optional[int] = None
        self.implicit_mode: bool = implicit_mode
        self.did_change: bool = False
        self.packed_data: Optional[bytes] = None
        self.packed_used_packets: bytes
        self.packed: bytes
        if not read_only:
            self.packed_used_packets = self.prepare_and_pack()
            self.packed = self.calculate_packed_data()
        self.prepend = prepend
        self.append = append

    @classmethod
    def from_packet(cls, packet: "Packet", pseudo: bool = False) -> "Packet":
        if not pseudo:
            data = packet.get_data()
        else:
            data = b""
        used_packets = packet.get_used_packets()
        number_of_packets = packet.get_total_number_of_chunks()
        res = cls(data, used_packets, number_of_packets)
        res.error_correction = packet.get_error_correction()
        res.dna_data = None
        return res

    def get_org_class(self) -> str:
        return self.__module__.split(".")[1]

    def set_used_packets(self, used_packets: Set[int]) -> None:
        self.used_packets = used_packets
        self.internal_hash = None

    def prepare_and_pack(self) -> bytes:
        # Format = Highest possible Packetnumber for this file, number of used Packets for this
        # File and the Indizies of the used Packets
        struct_str = (
            "<"
            + (self.number_of_chunks_len_format if self.save_number_of_chunks_in_packet else "")
            + (self.used_packets_len_format if not self.implicit_mode else "")
            + self.id_len_format
        )
        if self.save_number_of_chunks_in_packet:
            if self.implicit_mode:
                return struct.pack(
                    struct_str,
                    xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                    xor_mask(self.id, self.id_len_format),
                )
            else:
                return struct.pack(
                    struct_str,
                    xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                    xor_mask(len(self.used_packets), self.used_packets_len_format),
                    xor_mask(self.id, self.id_len_format),
                )
        else:
            if self.implicit_mode:
                return struct.pack(struct_str, xor_mask(self.id, self.id_len_format))
            else:
                return struct.pack(
                    struct_str,
                    xor_mask(len(self.used_packets), self.used_packets_len_format),
                    xor_mask(self.id, self.id_len_format),
                )

    def calculate_packed_data(self) -> bytes:
        # Laenge des Packets + UsedPackets + Data + crc
        self.packed_data = struct.pack("<" + str(len(self.data)) + "s", bytes(self.data))
        assert self.packed_data is not None
        payload = struct.pack(
            "<" + str(len(self.packed_used_packets)) + "s" + str(len(self.packed_data)) + "s",
            self.packed_used_packets,
            self.packed_data,
        )
        return self.error_correction(payload)

    def get_struct(self, split_to_multiple_files: bool) -> bytes:
        packed = self.packed
        if not split_to_multiple_files:
            return struct.pack(
                "<" + self.packet_len_format + str(len(packed)) + "s", len(packed), packed
            )
        else:
            return packed

    def get_dna_struct(
        self,
        split_to_multiple_files: bool,
        spacing: int = 0,
        spacing_length: int = 0,
        recalculate: bool = False,
    ) -> str:
        if recalculate:
            self.packed_used_packets = self.prepare_and_pack()
            self.packed = self.calculate_packed_data()
            self.dna_data = None
        if self.dna_data is None and self.error_correction.__name__ == "dna_reed_solomon_encode":
            self.dna_data = (
                self.prepend + quads2dna(self.get_struct(split_to_multiple_files)) + self.append
            )
        elif self.dna_data is None:
            self.dna_data = (
                self.prepend
                + interleave_spacing(
                    "".join(string2QUATS(self.get_struct(split_to_multiple_files))),
                    spacing,
                    spacing_length,
                )
                + self.append
            )
        assert self.dna_data is not None, "Should never happen"
        self.internal_hash = None  # enforce recalculation of hash
        return self.dna_data

    def get_data(self) -> Union[bytes, UInt8Array]:
        return self.data

    def set_data(self, data: Union[bytes, UInt8Array]) -> None:
        self.data = data

    def get_used_packets(self) -> Set[int]:
        return self.used_packets

    def get_bool_array_used_packets(self) -> BoolArray:
        return np.array(
            [x in self.used_packets for x in range(self.total_number_of_chunks)], dtype=bool
        )

    def set_bool_array_used_packet(self, b_array: List[bool]) -> None:
        assert len(b_array) == self.total_number_of_chunks, "Problem"
        self.set_used_packets({k for k in range(self.total_number_of_chunks) if b_array[k]})

    def get_total_number_of_chunks(self) -> int:
        return self.total_number_of_chunks

    def get_error_correction(self) -> ErrorCorrectionCallable:
        return self.error_correction

    def update_degree(self) -> None:
        self.degree = len(self.used_packets)

    def get_degree(self) -> int:
        return self.degree

    def remove_packets(self, packet_set: Set[int]) -> None:
        self.set_used_packets(self.used_packets.difference(packet_set))
        self.update_degree()
        self.did_change = True  # CRC is no longer valid

    def xor_and_remove_packet(self, packet: "Packet") -> None:
        self.remove_packets(packet.get_used_packets())
        self.data = xor_numpy(self.data, packet.get_data())
        self.did_change = True  # CRC is no longer valid

    def set_error_prob(self, error_prob: Optional[int] = None) -> None:
        self.error_prob = error_prob

    def __str__(self) -> str:
        return (
            "< "
            + "Id"
            + str(self.id)
            + "used_packets: "
            + str(self.used_packets)
            + " , Data: "
            + str(self.data)
            + " , Error Correction: "
            + str(self.error_correction)
            + " >"
        )

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other: Any) -> bool:
        # if self.error_prob is not None and other.error_prob is not None:
        #    return self.error_prob == other.error_prob
        # else:
        if not isinstance(other, self.__class__):
            return False
        self_data = self.data.tobytes() if isinstance(self.data, np.ndarray) else bytes(self.data)
        other_data = (
            other.data.tobytes() if isinstance(other.data, np.ndarray) else bytes(other.data)
        )
        return hash(self) == hash(other) and self_data == other_data

    def __lt__(self, other: "Packet") -> bool:
        if self.error_prob is not None and other.error_prob is not None:
            return self.error_prob < other.error_prob
        else:
            return min(self.get_used_packets()) < min(other.get_used_packets())

    def __hash__(self) -> int:
        if self.internal_hash is None:
            normalized_data = (
                self.data.tobytes() if isinstance(self.data, np.ndarray) else bytes(self.data)
            )
            self.internal_hash = hash(
                str(self.total_number_of_chunks)
                + str(self.id)
                + str("" if self.error_prob is None else self.error_prob)
                + self.__module__
                + (normalized_data.hex() if self.dna_data is None else self.dna_data)
            )
        assert self.internal_hash is not None
        return self.internal_hash


class ParallelPacket:
    def __init__(
        self,
        used_packets: Set[int],
        total_number_of_chunks: int,
        p_id: int,
        data: Union[bytes, UInt8Array],
        dna_data: typing.Optional[str],
        packed: bytes,
        error_prob: typing.Optional[int],
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
        save_number_of_chunks_in_packet: bool,
        org_class: str = "",
        prepend: str = "",
        append: str = "",
        org_hash: Optional[int] = None,
    ):
        self.used_packets: Set[int] = used_packets
        self.total_number_of_chunks: int = total_number_of_chunks
        self.id: int = p_id
        self.datadata: Union[bytes, UInt8Array] = data
        self.dna_data: typing.Optional[str] = dna_data
        self.packed: bytes = packed
        self.error_prob: typing.Optional[int] = error_prob
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.id_len_format: str = id_len_format
        self.safe_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.org_class: str = org_class.split(".")[1]
        self.bool_arrayused_packets: List[bool] = [
            x in self.used_packets for x in range(0, self.total_number_of_chunks)
        ]
        self.prepend: str = prepend
        self.append: str = append
        self.calculated_hash: Optional[int] = org_hash

    @classmethod
    def from_packet(cls, packet: Packet) -> "ParallelPacket":
        return ParallelPacket(
            packet.used_packets,
            packet.total_number_of_chunks,
            packet.id,
            packet.data,
            packet.dna_data,
            packet.packed,
            packet.error_prob,
            packet.packet_len_format,
            packet.crc_len_format,
            packet.number_of_chunks_len_format,
            packet.id_len_format,
            packet.save_number_of_chunks_in_packet,
            org_class=packet.__module__,
            prepend=packet.prepend,
            append=packet.append,
            org_hash=hash(packet),
        )

    def get_org_class(self) -> str:
        return self.org_class

    def __hash__(self) -> typing.Optional[int]:

        if self.calculated_hash is None:
            self.calculated_hash = hash(
                f"{self.total_number_of_chunks}{self.id}{self.error_prob}{self.org_class}{self.dna_data}"
            )
        return self.calculated_hash

    def __eq__(self, other: Any) -> bool:
        # if self.error_prob is not None and other.error_prob is not None:
        #    return self.error_prob == other.error_prob
        # else:
        return hash(self) == hash(other)

    def __lt__(self, other: "ParallelPacket") -> bool:
        if self.error_prob is not None and other.error_prob is not None:
            return self.error_prob < other.error_prob
        else:
            return min(self.get_used_packets()) < min(other.get_used_packets())

    def __gt__(self, other: "ParallelPacket") -> bool:
        if self.error_prob is not None and other.error_prob is not None:
            return self.error_prob > other.error_prob
        else:
            return min(self.get_used_packets()) > min(other.get_used_packets())

    def get_dna_struct(self, split_to_multiple_files: bool) -> str:
        if self.dna_data is None:
            raise RuntimeError("DNA-Data should not be None!")
        return self.dna_data

    def get_used_packets(self) -> Set[int]:
        return self.used_packets

    def get_struct(self, split_to_multiple_files: bool) -> bytes:
        packed = self.packed
        if not split_to_multiple_files:
            return struct.pack(
                "<" + self.packet_len_format + str(len(packed)) + "s", len(packed), packed
            )
        else:
            return packed
