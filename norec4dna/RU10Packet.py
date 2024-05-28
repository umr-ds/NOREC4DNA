#!/usr/bin/python
# -*- coding: latin-1 -*-
import logging
import struct
import typing
import numpy as np

from norec4dna.helper import xor_mask
from norec4dna.Packet import Packet
from bitstring import BitArray
from norec4dna.helper.helper import xor_with_seed
from norec4dna.helper.RU10Helper import intermediate_symbols
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import nocode


class RU10Packet(Packet):
    def __init__(self, data, used_packets: typing.Collection[int], total_number_of_chunks, id, dist=None,
                 read_only=False,
                 error_correction=nocode, packet_len_format="I", crc_len_format="L", number_of_chunks_len_format="L",
                 id_len_format="L", save_number_of_chunks_in_packet=True, method=None, window=None, prepend="",
                 append="", xor_by_seed=False, mask_id=True, id_spacing=0):
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
        if not read_only and (len(self.data) > 0 or self.data != ""):
            self.packed_used_packets = self.prepare_and_pack()
            self.packed = self.calculate_packed_data()
            self.get_dna_struct(True, self.id_spacing, self.id_spacing_length)
        else:
            self.packed_used_packets = None
            self.packed = None
        # super().__init__(data, used_packets, total_number_of_chunks, read_only, error_correction=error_correction)

    def set_used_packets_old(self, u_packets):
        self.used_packets = u_packets
        tmp_lst = np.zeros(self.total_number_of_chunks, dtype=bool)
        valid_indices = np.array(u_packets)[np.array(u_packets) < self.total_number_of_chunks]
        if len(u_packets) > 0:
            tmp_lst[valid_indices] = True
        else:
            logging.warning("Degenerated Packet! - No valid indices found for used packets: " + str(u_packets))
        self.internal_hash = hash(np.packbits(tmp_lst).tobytes())
        self.bool_arrayused_packets = tmp_lst
        self.update_degree()

    def set_used_packets(self, u_packets):
        self.used_packets = u_packets

        # Create boolean array directly
        #self.bool_arrayused_packets = np.zeros(self.total_number_of_chunks, dtype=bool)

        # Set True only for valid indices in u_packets
        #valid_indices = np.array(u_packets)[np.array(u_packets) < self.total_number_of_chunks]
        #self.bool_arrayused_packets[valid_indices] = True

        if len(u_packets) > 0:
            # Compute hash directly from boolean array
            try:
                self.internal_hash = hash(np.packbits(self.used_packets).tobytes())
            except:
                self.internal_hash = 0
            self.update_degree()
        else:
            logging.warning("Degenerated Packet! - No valid indices found for used packets: " + str(u_packets))

    def prepare_and_pack(self) -> bytes:
        # Format = Highest possible Packetnumber for this file,
        # number of used Packets for this File and the seed for the Indices of the used Packets
        if self.save_number_of_chunks_in_packet:
            return struct.pack("<" + self.number_of_chunks_len_format + self.id_len_format,
                               xor_mask(self.total_number_of_chunks, self.number_of_chunks_len_format),
                               xor_mask(self.id, self.id_len_format, enabled=self.mask_id))
        else:
            return struct.pack("<" + self.id_len_format, xor_mask(self.id, self.id_len_format, enabled=self.mask_id))

    def packMethod(self) -> bytes:
        if "window" not in self.method:
            if self.method == "even":
                data = BitArray(bin='00101010').tobytes()
            elif self.method == "odd":
                data = BitArray(bin='01101010').tobytes()
            else:
                raise RuntimeError("Unknown method: ", self.method)
        else:
            if self.method == "window_30":
                data_str = '10'
            elif self.method == "window_40":
                data_str = '01'
            else:
                raise RuntimeError("Unknown method: ", self.method)
            bin_win = bin(self.window)[2:]
            while len(data_str) + len(bin_win) < 8:
                data_str += '0'
            data_str += bin_win
            data = BitArray(bin=data_str).tobytes()
        return data

    def calculate_packed_data(self) -> bytes:
        # size of the packets + UsedPackets + Data + crc
        self.packed_data = struct.pack("<" + str(len(self.data)) + "s", bytes(self.data))
        if self.xor_by_seed:
            self.packed_data = xor_with_seed(self.packed_data, self.id)
        if self.packedMethod:
            payload = struct.pack(
                "<" + str(len(self.packed_used_packets)) + "s" + str(len(self.packed_data)) + "s" + str(
                    len(self.packedMethod)) + "s",  # method data
                self.packed_used_packets, self.packed_data, self.packedMethod)
        else:
            # i = 0
            # payload = b""
            # for fragment in self.packed_used_packets:
            #    payload += fragment.to_bytes(1, "little") + self.packed_data[i:i + self.id_spacing]
            #    i += self.id_spacing
            # payload += self.packed_data[i:]
            # self.packed_used_packets = ""
            payload = struct.pack("<" + str(len(self.packed_used_packets)) + "s" + str(len(self.packed_data)) + "s",
                                  self.packed_used_packets, self.packed_data)
        return self.error_correction(payload)  # proxy payload through dynamic error correction / detection

    def getId(self) -> int:
        return self.id

    def setId(self, id: int):
        self.id = id

    @classmethod
    def from_packet(cls, packet: 'RU10Packet', pseudo: bool = False) -> 'RU10Packet':
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
        if self.bool_arrayused_packets is None:
            self.bool_arrayused_packets = np.zeros(self.total_number_of_chunks, dtype=bool)

            # Set True only for valid indices in u_packets
            valid_indices = np.array(self.used_packets)[np.array(self.used_packets) < self.total_number_of_chunks]
            self.bool_arrayused_packets[valid_indices] = True
        return self.bool_arrayused_packets

    def get_valid_used_packets(self) -> typing.Set[int]:
        return set(filter(lambda x: x < self.total_number_of_chunks, self.used_packets))

    def get_bool_array_all_used_packets(self) -> typing.List[bool]:
        return [x in self.used_packets for x in
                range(
                    self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks())]

    def get_bool_array_used_and_ldpc_packets(self) -> typing.List[bool]:
        # speedup candidate
        u_bound = self.total_number_of_chunks + self.get_number_of_ldpc_blocks()
        tmp_lst = np.full((1, u_bound), False)
        for x in self.used_packets:
            if x < u_bound:
                tmp_lst[0, x] = True
        res = tmp_lst[0]
        del tmp_lst
        return res

    def get_used_and_ldpc_packets(self) -> typing.Set[int]:
        u_bound = self.total_number_of_chunks + self.get_number_of_ldpc_blocks()
        # get all number from self.use_packets smaller than u_bound:
        return set(filter(lambda x: x < u_bound, self.used_packets))

    def get_bool_array_ldpc_packets(self) -> typing.List[bool]:
        # speedup candidate
        return [x in self.used_packets for x in
                range(self.total_number_of_chunks, self.total_number_of_chunks + self.get_number_of_ldpc_blocks(), )]

    def get_ldpc_packets(self) -> typing.Set[int]:
        return set(filter(
            lambda x: self.total_number_of_chunks <= x < self.total_number_of_chunks + self.get_number_of_ldpc_blocks(),
            self.used_packets))

    def get_bool_array_half_packets_old(self) -> typing.List[bool]:
        return [x in self.used_packets for x in
                range(self.total_number_of_chunks + self.get_number_of_ldpc_blocks(),
                      self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks(), )]

    def get_bool_array_half_packets(self) -> typing.List[bool]:
        used_packets_set = set(self.used_packets)
        ldpc_blocks = self.get_number_of_ldpc_blocks()
        half_blocks = self.get_number_of_half_blocks()
        start = self.total_number_of_chunks + ldpc_blocks
        end = start + half_blocks
        return [x in used_packets_set for x in range(start, end)]

    def get_half_packets(self) -> typing.Set[int]:
        return set(filter(lambda
                              x: self.total_number_of_chunks + self.get_number_of_ldpc_blocks() <= x < self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks(),
                          self.used_packets))

    def get_bool_array_repair_packets_old(self) -> typing.List[bool]:
        return [x in self.used_packets for x in
                range(self.total_number_of_chunks,
                      self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks(), )]

    def get_bool_array_repair_packets(self) -> typing.List[bool]:
        used_packets_set = set(self.used_packets)
        start = self.total_number_of_chunks
        end = start + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks()
        return [x in used_packets_set for x in range(start, end)]
    def get_repair_packets(self):
        return set(filter(lambda
                              x: self.total_number_of_chunks <= x < self.total_number_of_chunks + self.get_number_of_ldpc_blocks() + self.get_number_of_half_blocks(),
                          self.used_packets))

    def __str__(self) -> str:
        return "< used_packets: " + str(self.used_packets) + " , Data: " + str(self.data) + " >"


if __name__ == "__main__":
    print("This class must not be called by itself.")
