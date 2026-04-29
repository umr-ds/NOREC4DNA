#!/usr/bin/python
# -*- coding: latin-1 -*-
from __future__ import annotations

import copy
import logging
import struct
import typing
from importlib import import_module, util

import numpy as np

from .distributions.RaptorDistribution import RaptorDistribution
from .ErrorCorrection import nocode
from .helper import xor_mask
from .helper.helper import xor_with_seed
from .helper.RU10Helper import intermediate_symbols
from .Packet import Packet


class _FallbackBitArray:
    def __init__(self, *, bin: str):
        self._bits = bin

    def tobytes(self) -> bytes:
        padded_bits = self._bits.ljust(((len(self._bits) + 7) // 8) * 8, "0")
        return bytes(
            int(padded_bits[index : index + 8], 2) for index in range(0, len(padded_bits), 8)
        )


_bitstring = import_module("bitstring") if util.find_spec("bitstring") is not None else None
BitArray = _bitstring.BitArray if _bitstring is not None else _FallbackBitArray


class RU10Packet(Packet):
    def __init__(
        self,
        data,
        used_packets: typing.Collection[int],
        total_number_of_chunks,
        id,
        dist=None,
        read_only=False,
        error_correction=nocode,
        packet_len_format="I",
        crc_len_format="L",
        number_of_chunks_len_format="L",
        id_len_format="L",
        save_number_of_chunks_in_packet=True,
        method=None,
        window=None,
        prepend="",
        append="",
        xor_by_seed=False,
        mask_id=True,
        id_spacing=0,
        contains_meta_or_version=False,
    ):
        self.id: int = id
        self.bool_arrayused_packets: typing.Optional[np.ndarray] = None
        self.total_number_of_chunks: int = total_number_of_chunks
        self.data: bytes = data
        self.used_packets: typing.Optional[typing.Iterable[int]] = None
        self.internal_hash: typing.Optional[int] = None
        self.set_used_packets(used_packets)
        self.degree: int = len(used_packets)
        self.dna_data: typing.Optional[str] = None
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.id_len_format: str = id_len_format
        self.error_correction: typing.Callable[[typing.Any], typing.Any] = error_correction
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.error_prob: typing.Optional[float] = None
        self.xor_by_seed = xor_by_seed
        self.mask_id = mask_id
        if id_spacing < 0:
            id_spacing = 0
        self.id_spacing = id_spacing
        self.id_spacing_length = struct.calcsize(id_len_format) * 4
        if method:
            self.method: typing.Optional[str] = method
            self.window: typing.Optional[int] = window
            self.packedMethod: typing.Optional[bytes] = self.packMethod()
        else:
            self.packedMethod = None
        if dist is None:
            self.dist = RaptorDistribution(total_number_of_chunks)
        else:
            self.dist = dist
        _, self.s, self.h = intermediate_symbols(total_number_of_chunks, self.dist)
        self.prepend = prepend
        self.append = append
        self.packed_used_packets: typing.Optional[bytes]
        self.packed: typing.Optional[bytes]
        if not read_only and (len(self.data) > 0 or self.data != ""):
            self.packed_used_packets = self.prepare_and_pack()
            self.packed = self.calculate_packed_data()
            self.get_dna_struct(True, self.id_spacing, self.id_spacing_length)
        else:
            self.packed_used_packets = None
            self.packed = None
        self.packed_struct: typing.Optional[bytes] = None
        self.contains_meta_or_version = contains_meta_or_version
        # super().__init__(data, used_packets, total_number_of_chunks, read_only, error_correction=error_correction)

    def get_packet_header_size(self) -> int:
        size = 0
        if self.save_number_of_chunks_in_packet:
            size += struct.calcsize(self.number_of_chunks_len_format)
        size += struct.calcsize(self.id_len_format)
        return size

    def set_used_packets(self, u_packets: typing.Collection[int]):
        self.used_packets = u_packets
        tmp_lst = np.zeros(self.total_number_of_chunks, dtype=bool)
        valid_indices = np.array(u_packets)[np.array(u_packets) < self.total_number_of_chunks]
        if len(u_packets) > 0:
            tmp_lst[valid_indices] = True
        else:
            logging.warning(
                "Degenerated Packet! - No valid indices found for used packets: " + str(u_packets)
            )
        self.internal_hash = hash(np.packbits(tmp_lst).tobytes())
        self.bool_arrayused_packets = tmp_lst
        self.update_degree()

    def prepare_and_pack(self) -> bytes:
        # Format = Highest possible Packetnumber for this file,
        # number of used Packets for this File and the seed for the Indices of the used Packets
        if self.save_number_of_chunks_in_packet:
            return struct.pack(
                "<" + self.number_of_chunks_len_format + self.id_len_format,
                xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                xor_mask(self.id, self.id_len_format, enabled=self.mask_id),
            )
        else:
            return struct.pack(
                "<" + self.id_len_format,
                xor_mask(self.id, self.id_len_format, enabled=self.mask_id),
            )

    def packMethod(self) -> bytes:
        assert self.method is not None, "method must not be None!"
        if "window" not in self.method:
            if self.method == "even":
                data = BitArray(bin="00101010").tobytes()
            elif self.method == "odd":
                data = BitArray(bin="01101010").tobytes()
            else:
                raise RuntimeError("Unknown method: ", self.method)
        else:
            if self.method == "window_30":
                data_str = "10"
            elif self.method == "window_40":
                data_str = "01"
            else:
                raise RuntimeError("Unknown method: ", self.method)
            assert self.window, "window must not be None!"
            bin_win = bin(self.window)[2:]
            while len(data_str) + len(bin_win) < 8:
                data_str += "0"
            data_str += bin_win
            data = BitArray(bin=data_str).tobytes()
        return data

    def calculate_packed_data(self) -> bytes:
        # size of the packets + UsedPackets + Data + crc
        self.packed_data: bytes = struct.pack("<" + str(len(self.data)) + "s", bytes(self.data))
        if self.xor_by_seed:
            self.packed_data: bytes = xor_with_seed(self.packed_data, self.id)
        assert self.packed_used_packets is not None
        if self.packedMethod:
            payload = struct.pack(
                "<"
                + str(len(self.packed_used_packets))
                + "s"
                + str(len(self.packed_data))
                + "s"
                + str(len(self.packedMethod))
                + "s",  # method data
                self.packed_used_packets,
                self.packed_data,
                self.packedMethod,
            )
        else:
            # i = 0
            # payload = b""
            # for fragment in self.packed_used_packets:
            #    payload += fragment.to_bytes(1, "little") + self.packed_data[i:i + self.id_spacing]
            #    i += self.id_spacing
            # payload += self.packed_data[i:]
            # self.packed_used_packets = ""
            payload = struct.pack(
                "<" + str(len(self.packed_used_packets)) + "s" + str(len(self.packed_data)) + "s",
                self.packed_used_packets,
                self.packed_data,
            )
        return self.error_correction(
            payload
        )  # proxy payload through dynamic error correction / detection

    def getId(self) -> int:
        return self.id

    def setId(self, id: int):
        self.id = id

    @classmethod
    def from_packet(cls, packet: "RU10Packet", pseudo: bool = False) -> "RU10Packet":
        if not pseudo:
            data = packet.get_data()
        else:
            data = ""
        used_packets = packet.get_used_packets()
        number_of_packets = packet.get_total_number_of_chunks()
        res = cls(data, used_packets, number_of_packets, packet.getId())
        res.error_correction = packet.get_error_correction()
        return res

    def get_number_of_half_blocks(self) -> int:
        return self.h

    def get_number_of_ldpc_blocks(self) -> int:
        return self.s

    def get_bool_array_used_packets(self) -> typing.Optional[typing.List[bool]]:
        return (
            self.bool_arrayused_packets.tolist()
            if self.bool_arrayused_packets is not None
            else None
        )

    def get_bool_array_all_used_packets(self) -> typing.List[bool]:
        used = self.used_packets or []
        return [
            x in used
            for x in range(
                self.total_number_of_chunks
                + self.get_number_of_ldpc_blocks()
                + self.get_number_of_half_blocks()
            )
        ]

    def get_bool_array_used_and_ldpc_packets(self) -> typing.List[bool]:
        # speedup candidate
        u_bound = self.total_number_of_chunks + self.get_number_of_ldpc_blocks()
        tmp_lst = np.full((1, u_bound), False)
        for x in self.used_packets or []:
            if x < u_bound:
                tmp_lst[0, x] = True
        res = tmp_lst[0]
        del tmp_lst
        return res

    def get_bool_array_ldpc_packets(self) -> typing.List[bool]:
        # speedup candidate
        return [
            x in (self.used_packets or [])
            for x in range(
                self.total_number_of_chunks,
                self.total_number_of_chunks + self.get_number_of_ldpc_blocks(),
            )
        ]

    def get_bool_array_half_packets(self) -> typing.List[bool]:
        return [
            x in (self.used_packets or [])
            for x in range(
                self.total_number_of_chunks + self.get_number_of_ldpc_blocks(),
                self.total_number_of_chunks
                + self.get_number_of_ldpc_blocks()
                + self.get_number_of_half_blocks(),
            )
        ]

    def get_bool_array_repair_packets(self) -> typing.List[bool]:
        return [
            x in (self.used_packets or [])
            for x in range(
                self.total_number_of_chunks,
                self.total_number_of_chunks
                + self.get_number_of_ldpc_blocks()
                + self.get_number_of_half_blocks(),
            )
        ]

    def __str__(self) -> str:
        return "< used_packets: " + str(self.used_packets) + " , Data: " + str(self.data) + " >"

    # copy method: create a deep copy of the packet:
    def copy(self) -> "RU10Packet":
        new_packet = RU10Packet(
            data=bytes(self.data),
            used_packets=(
                list(copy.deepcopy(self.used_packets)) if self.used_packets is not None else []
            ),
            total_number_of_chunks=self.total_number_of_chunks,
            id=self.id,
            dist=self.dist,
            read_only=False,
            error_correction=self.error_correction,
            packet_len_format=self.packet_len_format,
            crc_len_format=self.crc_len_format,
            number_of_chunks_len_format=self.number_of_chunks_len_format,
            id_len_format=self.id_len_format,
            save_number_of_chunks_in_packet=self.save_number_of_chunks_in_packet,
            method=self.method if hasattr(self, "method") else None,
            window=self.window if hasattr(self, "window") else None,
            prepend=self.prepend,
            append=self.append,
            xor_by_seed=self.xor_by_seed,
            mask_id=self.mask_id,
            id_spacing=self.id_spacing,
        )
        return new_packet


if __name__ == "__main__":
    print("This class must not be called by itself.")
