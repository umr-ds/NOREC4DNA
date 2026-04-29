#!/usr/bin/python
# -*- coding: latin-1 -*-
"""Belief-propagation decoder for LT-coded files and DNA strands."""

from __future__ import annotations

import os
import struct
import time
import typing
from collections import deque

import numpy as np
from reedsolo import ReedSolomonError

from .BPDecoder import BPDecoder
from .DecodePacket import DecodePacket
from .distributions.Distribution import Distribution
from .ErrorCorrection import crc32, nocode
from .HeaderChunk import HeaderChunk
from .helper import calc_crc, xor_mask
from .helper.quaternary2Bin import quat_file_to_bin
from .Packet import Packet

if typing.TYPE_CHECKING:
    from .OnlinePacket import OnlinePacket
    from .RU10Packet import RU10Packet


class LTBPDecoder(BPDecoder):
    def __init__(
        self,
        file: typing.Optional[str] = None,
        error_correction: typing.Callable[[bytes], bytes] = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: typing.Optional[int] = None,
        implicit_mode: bool = True,
        dist: typing.Optional[Distribution] = None,
    ):
        super().__init__(file, error_correction, use_headerchunk, static_number_of_chunks)
        self.implicit_mode: bool = implicit_mode
        self.use_headerchunk: bool = use_headerchunk
        self.file: typing.Optional[str] = file
        self.decodedPackets: typing.Dict[int, Packet] = {}
        self.degreeToPacket: typing.Dict[int, typing.Set[Packet]] = {}
        if file is not None:
            self.isFolder: bool = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        self.number_of_chunks: int = 1000000
        self.headerChunk: typing.Optional[HeaderChunk] = None
        self.queue: deque = deque()
        self.error_correction: typing.Callable[[bytes], bytes] = error_correction
        self.static_number_of_chunks: typing.Optional[int] = static_number_of_chunks
        self.dist: typing.Optional[Distribution] = (
            dist  # if implicit_mode is True, dist MUST be != None
        )

    def _store_decoded_packet(self, packet: Packet) -> None:
        if packet.get_degree() != 1:
            return
        [chunk_index] = packet.get_used_packets()
        if chunk_index == 0 and self.use_headerchunk:
            if self.headerChunk is None:
                self.headerChunk = HeaderChunk(packet)
            return
        self.decodedPackets[chunk_index] = packet

    def _update_decode_formats(
        self, number_of_chunks_len_format: str, degree_len_format: str
    ) -> typing.Tuple[str, str]:
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = ""
        if self.implicit_mode:
            degree_len_format = ""
        return number_of_chunks_len_format, degree_len_format

    def _open_folder_packet(self, file_name: str) -> None:
        assert self.file is not None, "file must be set before opening folder packets"
        self.EOF = False
        file_path = self.file + "/" + file_name
        self.f = quat_file_to_bin(file_path) if file_name.endswith("DNA") else open(file_path, "rb")

    def _decode_folder_packet(
        self,
        file_name: str,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        degree_len_format: str,
        seed_len_format: str,
        last_chunk_len_format: str,
    ) -> bool:
        if not (file_name.endswith(".LT") or file_name.endswith("DNA")):
            return False
        self._open_folder_packet(file_name)
        new_pack = self.getNextValidPacket(
            True,
            packet_len_format=packet_len_format,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            degree_len_format=degree_len_format,
            seed_len_format=seed_len_format,
            last_chunk_len_format=last_chunk_len_format,
        )
        return new_pack is not None and self.input_new_packet(new_pack)

    def _get_sorted_packets(self) -> typing.List[Packet]:
        return [self.decodedPackets[chunk_index] for chunk_index in sorted(self.decodedPackets)]

    def _ensure_header_chunk(self, last_chunk_len_format: str) -> None:
        if not self.use_headerchunk or self.headerChunk is not None:
            return
        for packet in self.degreeToPacket.get(1, set()):
            if packet.get_used_packets() == {0}:
                self.headerChunk = HeaderChunk(packet, last_chunk_len_format=last_chunk_len_format)
                break

    def _default_output_file_name(self) -> str:
        assert self.file is not None, "Cannot save file: file is None"
        return "DEC_" + os.path.basename(self.file.split("\x00")[0])

    def _resolve_output_file_name(self) -> str:
        file_name = self._default_output_file_name()
        if self.headerChunk is None:
            return file_name
        try:
            header_file_name = self.headerChunk.get_file_name()
            return (
                header_file_name.decode("utf-8")
                if isinstance(header_file_name, bytes)
                else header_file_name
            )
        except (UnicodeDecodeError, AttributeError) as exc:
            print(
                f"Warning: Could not decode filename from header chunk ({exc}), using default filename"
            )
            return file_name

    def _packet_output_bytes(
        self, decoded: Packet, null_is_terminator: bool
    ) -> typing.Tuple[bytes, bool]:
        if (
            self.number_of_chunks - 1 in decoded.get_used_packets()
            and self.use_headerchunk
            and self.headerChunk is not None
        ):
            output = decoded.get_data()[0 : self.headerChunk.get_last_chunk_length()]
            return output if isinstance(output, bytes) else output.tobytes(), False
        data = decoded.get_data()
        data_bytes = data if isinstance(data, bytes) else data.tobytes()
        if not null_is_terminator:
            return data_bytes, False
        splitter = data_bytes.decode().split("\x00")
        return splitter[0].encode(), len(splitter) > 1

    def _read_packet_data(
        self, file_handle: typing.IO[bytes], from_multiple_files: bool, packet_len_format: str
    ) -> typing.Tuple[bytes, int]:
        if from_multiple_files:
            packet = file_handle.read()
            return packet, len(packet)
        packet_len_raw = file_handle.read(struct.calcsize("<" + packet_len_format))
        packet_len = struct.unpack("<" + packet_len_format, packet_len_raw)[0]
        packet = file_handle.read(int(packet_len))
        return packet, int(packet_len)

    def _retry_read_packet(
        self,
        from_multiple_files: bool,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        degree_len_format: str,
        seed_len_format: str,
        last_chunk_len_format: str,
    ) -> typing.Optional[Packet]:
        return self.getNextValidPacket(
            from_multiple_files,
            packet_len_format=packet_len_format,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            degree_len_format=degree_len_format,
            seed_len_format=seed_len_format,
            last_chunk_len_format=last_chunk_len_format,
        )

    def _decode_crc_packet(self, packet: bytes, crc_len_format: str) -> typing.Optional[int]:
        crc_len = -struct.calcsize("<" + crc_len_format)
        payload = packet[:crc_len]
        crc: typing.Union[int, typing.Any] = struct.unpack("<" + crc_len_format, packet[crc_len:])[
            0
        ]
        calced_crc = calc_crc(payload)
        if crc == calced_crc:
            return crc_len
        print("[-] CRC-Error - " + str(hex(crc)) + " != " + str(hex(calced_crc)))
        self.corrupt += 1
        return None

    def _convert_xor_value(self, value: typing.Any) -> int:
        if isinstance(value, np.ndarray):
            return int(value.item())
        if isinstance(value, (int, float)):
            return int(value)
        return int(value)

    def _decode_lt_header(
        self,
        len_data: typing.Tuple[typing.Any, ...],
        number_of_chunks_len_format: str,
        degree_len_format: str,
        seed_len_format: str,
    ) -> typing.Tuple[int, typing.Set[int]]:
        degree: typing.Optional[int] = None
        if self.static_number_of_chunks is None:
            if self.implicit_mode:
                number_of_chunks_raw, seed_raw = len_data
            else:
                number_of_chunks_raw, degree_raw, seed_raw = len_data
                degree = int(degree_raw)
            self.number_of_chunks = self._convert_xor_value(
                xor_mask(number_of_chunks_raw, number_of_chunks_len_format)
            )
        else:
            if self.implicit_mode:
                (seed_raw,) = len_data
            else:
                degree_raw, seed_raw = len_data
                degree = int(degree_raw)
        seed_value = self._convert_xor_value(xor_mask(seed_raw, seed_len_format))
        if degree is None:
            if self.dist is not None:
                self.dist.set_seed(seed_value)
                degree = int(self.dist.getNumber())
            else:
                degree = 1
        else:
            degree = self._convert_xor_value(xor_mask(degree_raw, degree_len_format))
        return seed_value, self.choose_packet_numbers(degree, seed=seed_value)

    def decodeFolder(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        degree_len_format: str = "I",
        seed_len_format: str = "I",
        last_chunk_len_format: str = "I",
    ) -> typing.Optional[int]:
        decoded: bool = False
        self.EOF: bool = False
        if self.file is None:
            return None
        number_of_chunks_len_format, degree_len_format = self._update_decode_formats(
            number_of_chunks_len_format, degree_len_format
        )
        for file_name in os.listdir(self.file):
            decoded = self._decode_folder_packet(
                file_name,
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                degree_len_format,
                seed_len_format,
                last_chunk_len_format,
            )
            if decoded:
                break
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        if self.f is not None:
            self.f.close()
        if not decoded and self.EOF:
            print("Unable to retrieve file from chunks. Too many errors?")
            return -1

    def decodeFile(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        degree_len_format: str = "I",
        seed_len_format: str = "I",
        last_chunk_len_format: str = "I",
    ) -> typing.Optional[int]:
        decoded: bool = False
        self.EOF: bool = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        if self.implicit_mode:
            degree_len_format = ""
        while not (decoded or self.EOF):
            new_pack: typing.Optional[Packet] = self.getNextValidPacket(
                False,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                degree_len_format=degree_len_format,
                seed_len_format=seed_len_format,
                last_chunk_len_format=last_chunk_len_format,
            )
            if new_pack is None:
                break
            decoded = self.input_new_packet(new_pack)
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        if self.f is not None:
            self.f.close()
        if not decoded and self.EOF:
            print("Unable to retrieve file from chunks. Too many errors?")
            return -1

    def decodeZip(self, *args: typing.Any, **kwargs: typing.Any) -> typing.Optional[typing.Any]:
        raise RuntimeError("Not implemented for BP-based decoders in the current version!")

    def input_new_packet(self, packet: Packet) -> bool:
        """Used for easy PseudoDecode"""
        if not isinstance(packet, DecodePacket):
            decode_packet: Packet = DecodePacket.from_packet(packet)
            self.addPacket(decode_packet)
            return self.updatePackets(decode_packet)
        self.addPacket(packet)
        return self.updatePackets(packet)

    def addPacket(self, packet: Packet) -> None:
        if (packet.get_degree() not in self.degreeToPacket) or (
            not isinstance(self.degreeToPacket[packet.get_degree()], set)
        ):
            self.degreeToPacket[packet.get_degree()] = set()
        self.number_of_chunks = packet.get_total_number_of_chunks()
        self.degreeToPacket[packet.get_degree()].add(packet)

    def updatePackets(self, packet: Packet) -> bool:
        if packet.get_degree() == 1:
            self._store_decoded_packet(packet)
            if self.is_decoded():
                return True
        self.queue.append(packet)
        finished: bool = False
        while len(self.queue) > 0 and not finished:
            finished = self.reduceAll(self.queue.popleft())
        return finished

    def compareAndReduce(self, packet: Packet, other: Packet) -> typing.Union[bool, int]:
        if self.file is None:  # In case of PseudoDecode: DO NOT REALLY COMPUTE XOR
            packet.remove_packets(other.get_used_packets())
        else:
            packet.xor_and_remove_packet(other)
        degree = packet.get_degree()
        if (degree not in self.degreeToPacket) or (
            not isinstance(self.degreeToPacket[degree], set)
        ):
            self.degreeToPacket[degree] = set()
        if degree == 1:
            self._store_decoded_packet(packet)
        self.degreeToPacket[degree].add(packet)
        if self.is_decoded():
            return True
        self.queue.append(packet)
        return degree

    """def reduceAll(self, packet: Packet) -> bool:
        # looup all packets for this to solve with ( when this packet has a subset of used Packets)
        fin: bool = False

        lookup: typing.List[int] = [i for i in self.degreeToPacket.keys() if packet.get_degree() < i]
        for i in lookup:
            if not isinstance(self.degreeToPacket[i], set):
                self.degreeToPacket[i] = set()
            for p in self.degreeToPacket[i].copy():
                p_used = p.get_used_packets()
                pack_used = packet.get_used_packets()
                if len(pack_used) < len(p_used) and pack_used.issubset(p_used):
                    self.degreeToPacket[i].remove(p)
                    degree = self.compareAndReduce(p, packet)
                    if isinstance(degree, bool) and degree:
                        return degree
        degree = packet.get_degree()
        lookup = [i for i in self.degreeToPacket.keys() if packet.get_degree() > i]
        for i in lookup:
            if not isinstance(self.degreeToPacket[i], set):
                self.degreeToPacket[i] = set()
            for p in self.degreeToPacket[i].copy():
                p_used = p.get_used_packets()
                pack_used = packet.get_used_packets()
                if len(pack_used) > len(p_used) and p_used.issubset(pack_used):
                    try:
                        self.degreeToPacket[degree].remove(packet)
                        degree = self.compareAndReduce(packet, p)
                        if isinstance(degree, bool) and degree:
                            return degree
                    except Exception:
                        continue
        return fin or self.is_decoded()"""

    def is_decoded(self) -> bool:
        solved_chunks = len(self.decodedPackets)
        required_chunks = self.number_of_chunks - (1 if self.use_headerchunk else 0)
        if solved_chunks < required_chunks:
            return False
        return not self.use_headerchunk or self.headerChunk is not None

    def removeAndXorAuxPackets(
        self, packet: typing.Union[Packet, "OnlinePacket", "RU10Packet"]
    ) -> typing.List[bool]:
        """Override parent class method - LTBPDecoder doesn't use aux packets"""
        return []

    def getSolvedCount(self) -> int:
        return len(self.decodedPackets) + (1 if self.headerChunk is not None else 0)

    def getNextValidPacket(
        self,
        from_multiple_files: bool = False,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        degree_len_format: str = "I",
        seed_len_format: str = "I",
        last_chunk_len_format: str = "I",
    ) -> typing.Optional[Packet]:
        if self.f is None:
            raise RuntimeError("Input file not open")
        file_handle = self.f
        packet, packet_len = self._read_packet_data(
            file_handle, from_multiple_files, packet_len_format
        )
        if not packet or not packet_len:  # EOF
            self.EOF: bool = True
            file_handle.close()
            return None
        crc_len: typing.Optional[int]
        if self.error_correction.__code__.co_name == crc32.__code__.co_name:
            crc_len = self._decode_crc_packet(packet, crc_len_format)
            if crc_len is None:
                return self._retry_read_packet(
                    from_multiple_files,
                    packet_len_format,
                    crc_len_format,
                    number_of_chunks_len_format,
                    degree_len_format,
                    seed_len_format,
                    last_chunk_len_format,
                )
        else:
            crc_len = None
            try:
                packet = self.error_correction(packet)
            except (AssertionError, ReedSolomonError, ValueError):
                self.corrupt += 1
                return self._retry_read_packet(
                    from_multiple_files,
                    packet_len_format,
                    crc_len_format,
                    number_of_chunks_len_format,
                    degree_len_format,
                    seed_len_format,
                    last_chunk_len_format,
                )

        struct_str = "<" + number_of_chunks_len_format + degree_len_format + seed_len_format
        struct_len = struct.calcsize(struct_str)
        len_data = struct.unpack(struct_str, packet[0:struct_len])
        _, used_packets = self._decode_lt_header(
            len_data, number_of_chunks_len_format, degree_len_format, seed_len_format
        )
        data = packet[struct_len:crc_len]

        self.correct += 1
        res = DecodePacket(
            np.frombuffer(data, dtype=np.uint8),
            used_packets,
            error_correction=self.error_correction,
            number_of_chunks=self.number_of_chunks,
        )
        if used_packets.issubset({0}) and self.headerChunk is None and self.use_headerchunk:
            self.headerChunk = HeaderChunk(res)
        return res

    def choose_packet_numbers(self, degree: int, seed: int = 0) -> typing.Set[int]:
        assert degree <= self.number_of_chunks
        res: typing.Set[int] = set()
        rng = np.random
        rng.seed(int(seed))
        for _ in range(0, degree):
            tmp = int(rng.choice(range(0, self.number_of_chunks)))
            while tmp in res:
                tmp = int(rng.choice(range(0, self.number_of_chunks)))
            res.add(tmp)
        return res

    def saveDecodedFile(
        self,
        last_chunk_len_format: str = "I",
        null_is_terminator: bool = False,
        print_to_output: bool = True,
    ) -> None:
        assert self.is_decoded(), "Can not save File: Unable to reconstruct."
        assert self.file is not None, "Cannot save file: file is None"
        self._ensure_header_chunk(last_chunk_len_format)
        file_name = self._resolve_output_file_name()
        output_concat: bytes = b""
        with open(file_name, "wb") as f:
            for decoded in self._get_sorted_packets():
                output_bytes, stop_writing = self._packet_output_bytes(decoded, null_is_terminator)
                output_concat += output_bytes
                f.write(output_bytes)
                if stop_writing:
                    break

        print("Saved file as '" + str(file_name) + "'")
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))


if __name__ == "__main__":
    filename = "LT_logo.jpg"
    p_start = time.time()
    x = LTBPDecoder(filename)
    x.decode()
    p_end = time.time() - p_start
    print(
        "LT_Approx_Decode,"
        + str(x.number_of_chunks)
        + ","
        + str(len(x.decodedPackets))
        + ","
        + str(p_end)
        + ","
        + filename
    )
    x.saveDecodedFile()
