from __future__ import annotations

import argparse
import logging
import os
import struct
from io import BytesIO
from math import ceil, floor
from typing import Any, Callable, List, Optional, Union

import numpy as np
from PIL import Image
from reedsolo import ReedSolomonError

from .BPDecoder import BPDecoder
from .distributions.Distribution import Distribution
from .distributions.RaptorDistribution import RaptorDistribution
from .ErrorCorrection import get_error_correction_decode, nocode
from .HeaderChunk import HeaderChunk
from .helper import bitSet, buildGraySequence, logical_xor, xor_mask
from .helper.quaternary2Bin import quad_file_to_bytes, quat_file_to_bin, tranlate_quat_to_byte
from .helper.RU10Helper import choose_packet_numbers, from_true_false_list, intermediate_symbols
from .Packet import Packet
from .RU10IntermediatePacket import RU10IntermediatePacket
from .RU10Packet import RU10Packet

logger = logging.getLogger(__name__)

IntArray = np.ndarray[Any, np.dtype[np.int_]]


class RU10BPDecoder(BPDecoder):
    def __init__(
        self,
        file: Optional[str] = None,
        error_correction: Callable[..., Any] = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: Optional[int] = None,
        use_method: bool = False,
    ):
        super().__init__()
        self.file: Optional[str] = file
        self.use_method: bool = use_method
        self.f: Optional[Any] = None
        self.auxBlocks: dict[int, RU10IntermediatePacket] = {}
        self.isFolder: bool = False
        if file is not None:
            self.isFolder = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(file, "rb")
        self.number_of_chunks: int = 1000000
        self.s: int = -1
        self.h: int = -1
        self.error_correction: Callable[..., Any] = error_correction
        self.use_headerchunk: bool = use_headerchunk
        self.static_number_of_chunks: Optional[int] = static_number_of_chunks
        self.dist: Optional[Distribution] = None

    def _update_decode_number_format(self, number_of_chunks_len_format: str) -> str:
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            return ""
        return number_of_chunks_len_format

    def _open_folder_packet(self, filepath: str, is_dna_file: bool) -> bool:
        if not is_dna_file:
            self.f = open(filepath, "rb")
            return True
        loader = (
            quad_file_to_bytes
            if self.error_correction.__name__ == "dna_reed_solomon_decode"
            else quat_file_to_bin
        )
        try:
            self.f = loader(filepath)
        except TypeError:
            logger.warning("skipping CORRUPT file - contains illegal character(s)")
            self.corrupt += 1
            return False
        return True

    def _decode_folder_packet(
        self,
        filename: str,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> bool:
        if not (filename.endswith(".RU10") or filename.endswith("DNA")):
            return False
        assert self.file is not None, "file must be set before decoding folder packets"
        self.EOF = False
        filepath = self.file + "/" + filename
        if not self._open_folder_packet(filepath, filename.endswith("DNA")):
            return False
        new_pack = self.getNextValidPacket(
            True,
            packet_len_format=packet_len_format,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            id_len_format=id_len_format,
        )
        return new_pack is not None and self.input_new_packet(new_pack)

    def _default_output_file_name(self) -> str:
        return "DEC_" + os.path.basename(self.file) if self.file is not None else "RU10.BIN"

    def _resolve_output_file_name(self) -> str:
        file_name = self._default_output_file_name()
        if self.headerChunk is None:
            return file_name
        header_file_name = self.headerChunk.get_file_name()
        resolved = (
            header_file_name.decode("utf-8")
            if isinstance(header_file_name, bytes)
            else header_file_name
        )
        return resolved.split("\x00")[0]

    def _should_write_decoded_packet(self, chunk_number: int) -> bool:
        return chunk_number != 0 or not self.use_headerchunk or self.number_of_chunks - 1 == 0

    def _packet_output_bytes(
        self, decoded: Packet, chunk_number: int, null_is_terminator: bool
    ) -> tuple[bytes, bool]:
        if (
            self.number_of_chunks - 1 == chunk_number
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

    def _open_single_file_input(self) -> bool:
        assert self.file is not None, "file must be set before decoding"
        if not self.file.lower().endswith("dna"):
            return True
        try:
            if self.f is not None:
                self.f.close()
            self.f = quat_file_to_bin(self.file)
        except TypeError:
            logger.warning("skipping CORRUPT file - contains illegal character(s)")
            self.corrupt += 1
            return False
        return True

    def _open_fasta_input(self) -> None:
        assert self.file is not None, "file must be set before decoding"
        if self.f is not None:
            self.f.close()
        self.f = open(self.file, "r")

    def _read_fasta_entry(self) -> Optional[tuple[str, str, str]]:
        assert self.f is not None, "fasta file handle must be open"
        line = self.f.readline()
        if not line:
            self.EOF = True
            return None
        try:
            error_prob, seed = line[1:].replace("\n", "").split("_")
        except ValueError:
            error_prob, seed = "0", "0"
        line = self.f.readline()
        if not line:
            self.EOF = True
            return None
        return error_prob, seed, line.replace("\n", "")

    def _decode_fasta_packet(
        self,
        dna_str: str,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> Optional[RU10Packet]:
        return self.parse_raw_packet(
            BytesIO(tranlate_quat_to_byte(dna_str)).read(),
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            packet_len_format=packet_len_format,
            id_len_format=id_len_format,
        )

    def _decode_fasta_file(
        self,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> bool:
        self._open_fasta_input()
        decoded = False
        while not (decoded or self.EOF):
            entry = self._read_fasta_entry()
            if entry is None:
                break
            _, _, dna_str = entry
            new_pack = self._decode_fasta_packet(
                dna_str,
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
            if new_pack is not None:
                decoded = self.input_new_packet(new_pack)
        return decoded

    def _decode_binary_file(
        self,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> bool:
        decoded = False
        while not (decoded or self.EOF):
            new_pack = self.getNextValidPacket(
                False,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                id_len_format=id_len_format,
            )
            if new_pack is None:
                break
            decoded = self.input_new_packet(new_pack)
        return decoded

    def _packet_method_data(self, data: bytes) -> tuple[bytes, Optional[list[int]]]:
        if not self.use_method:
            return data, None
        method_data = bin(data[-1])[2:].rjust(8, "0")
        return data[:-1], self._chunk_list_from_method_data(method_data)

    def _chunk_list_from_method_data(self, method_data: str) -> list[int]:
        if method_data.startswith("00"):
            return [chunk for chunk in range(0, self.number_of_chunks + 1) if chunk % 2 == 0]
        if method_data.startswith("01"):
            return [chunk for chunk in range(0, self.number_of_chunks + 1) if chunk % 2 != 0]
        if method_data.startswith("10"):
            return self._window_chunk_list(method_data, 30)
        if method_data.startswith("11"):
            return self._window_chunk_list(method_data, 40)
        raise RuntimeError("Invalid method_data: %s" % method_data)

    def _window_chunk_list(self, method_data: str, window_size: int) -> list[int]:
        window = int(method_data[2:], 2)
        start = window * (window_size - 10)
        return [
            chunk for chunk in range(start, start + window_size) if chunk <= self.number_of_chunks
        ]

    def _update_packet_header(
        self, len_data: tuple[Any, ...], number_of_chunks_len_format: str, id_len_format: str
    ) -> int:
        if self.static_number_of_chunks is None:
            self.number_of_chunks = xor_mask(len_data[0], number_of_chunks_len_format)
            return xor_mask(len_data[1], id_len_format)
        return xor_mask(len_data[0], id_len_format)

    def _ensure_distribution_ready(self) -> None:
        if self.dist is None:
            self.dist = RaptorDistribution(self.number_of_chunks)
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.dist)
        if self.correct == 0:
            self.createAuxBlocks()
        self.correct += 1

    def _choose_used_packets(self, unxored_id: int, chunk_lst: Optional[list[int]]) -> set[int]:
        assert isinstance(self.dist, RaptorDistribution)
        if chunk_lst is None:
            return set(
                choose_packet_numbers(
                    self.number_of_chunks, unxored_id, self.dist, systematic=False
                )
            )
        numbers = choose_packet_numbers(
            len(chunk_lst), unxored_id, self.dist, systematic=False, max_l=len(chunk_lst)
        )
        return {chunk_lst[index] for index in numbers}

    def decodeFolder(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "I",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
    ) -> Optional[int]:
        """
        Decodes the information from a folder if self.file represents a folder and the packets were saved
        in multiple files and prints the number of decoded and corrupted packets.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: -1 if the decoding wasn't successful
        """
        if self.file is None:
            logger.error("Error: No file specified")
            return -1
        decoded = False
        self.EOF = False
        number_of_chunks_len_format = self._update_decode_number_format(number_of_chunks_len_format)
        for filename in os.listdir(self.file):
            decoded = self._decode_folder_packet(
                filename,
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
            if decoded:
                break
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets: %s", self.corrupt)
        if self.f is not None:
            self.f.close()
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve File from Chunks. Too many errors?")
            return -1

    def decodeFile(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
    ) -> Optional[int]:
        """
        Decodes the information from a file if self.file represents a file and the packets were saved in a single file.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: -1 if the decoding wasn't successful
        """
        if self.file is None:
            logger.error("Error: No file specified")
            return -1
        self.EOF = False
        decoded = False
        if not self._open_single_file_input():
            decoded = False
        number_of_chunks_len_format = self._update_decode_number_format(number_of_chunks_len_format)
        if self.file.lower().endswith("fasta"):
            decoded = self._decode_fasta_file(
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
        else:
            decoded = self._decode_binary_file(
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets : %s", self.corrupt)
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve file from chunks. Too many errors?")
            return -1

    def getNumberOfLDPCBlocks(self):
        return self.s

    def getNumberOfHalfBlocks(self):
        return self.h

    def getNumberOfRepairBlocks(self):
        return self.getNumberOfHalfBlocks() + self.getNumberOfLDPCBlocks()

    def removeAndXorAuxPackets(self, packet: RU10Packet) -> List[bool]:
        """
        Removes auxpackets (LDCP and Half) from a given packet to get the packets data.
        :param packet: Packet to remove auxpackets from
        :return: The data without the auxpackets
        """
        aux_mapping = self.getHalfPacketListFromPacket(packet)  # Enthaelt Data + LDPC Nummern
        aux_mapping.append(packet.get_bool_array_used_and_ldpc_packets())
        xored_list = logical_xor(aux_mapping)
        tmp = set(from_true_false_list(xored_list))  # Nur noch Data + LDPC sind vorhanden
        if self.debug:
            logger.debug("%s", tmp)
        tmp = type(packet)("", tmp, self.number_of_chunks, packet.id, packet.dist, read_only=True)
        aux_mapping = self.getAuxPacketListFromPacket(tmp)
        aux_mapping.append(tmp.get_bool_array_used_packets())  # [-len(self.auxBlocks):])
        return logical_xor(aux_mapping).tolist()

    def input_new_packet(self, packet: RU10Packet):
        """
        Removes auxpackets (LDPC and Half) and adds the remaining data to the GEPP matrix.
        :param packet: A Packet to add to the GEPP matrix
        :return: True: If solved. False: Else.
        """
        if self.auxBlocks == {} and self.dist is None:  # self.isPseudo and
            self.dist = RaptorDistribution(self.number_of_chunks)
            self.number_of_chunks = packet.get_total_number_of_chunks()
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.dist)
            self.createAuxBlocks()
        # we need to do it twice sine half symbols may contain ldpc symbols (which by definition are repair codes.)
        if self.debug:
            logger.debug("----")
            logger.debug("Id = %s", packet.id)
            logger.debug("%s", packet.used_packets)
        removed = self.removeAndXorAuxPackets(packet)
        if self.debug:
            logger.debug("%s", from_true_false_list(removed))
            logger.debug("%s", packet.get_error_correction())
            logger.debug("----")
        packet.set_used_packets(set(from_true_false_list(removed)))
        if self.count:
            for i in range(len(removed)):
                if i in self.counter.keys():
                    if removed[i]:
                        self.counter[i] += 1
                else:
                    self.counter[i] = 1
        self.addPacket(packet)
        return self.updatePackets(packet)

    def createAuxBlocks(self):
        """
        Reconstructs the auxblocks to be able to remove them afterwards.
        :return:
        """
        assert (
            self.number_of_chunks is not None
        ), "createAuxBlocks can only be called AFTER first Packet"
        if self.debug:
            logger.debug(
                "We should have %s LDPC-Blocks, %s Half-Blocks and %s normal Chunks (including 1 HeaderChunk)",
                self.getNumberOfLDPCBlocks(),
                self.getNumberOfHalfBlocks(),
                self.number_of_chunks,
            )
        for i in range(0, self.getNumberOfRepairBlocks()):
            self.repairBlockNumbers[i] = set()
        i = 0
        for group in self.generateIntermediateBlocksFormat(self.number_of_chunks):
            for elem in group:
                self.repairBlockNumbers[i] = set(elem)
                i += 1
        # XOR all Chunks into the corresponding AUX-Block
        for aux_number in self.repairBlockNumbers.keys():
            self.auxBlocks[aux_number] = RU10IntermediatePacket(
                "",
                self.repairBlockNumbers[aux_number],
                total_number_of_chunks=self.number_of_chunks,
                id=aux_number,
                dist=self.dist,
            )  # # We will add the Data once we have it.
            if self.debug:
                logger.debug("%s : %s", aux_number, self.auxBlocks[aux_number].used_packets)
        # Correct

    def getAuxPacketListFromPacket(self, packet: RU10Packet):
        """
        Creates a list for a packet with information about whether auxpackets have been used for that packet.
        :param packet: The packet to check.
        :return: Information about used auxpackets.
        """
        res = []
        aux_used_packets = packet.get_bool_array_repair_packets()
        for i in range(len(aux_used_packets)):
            if aux_used_packets[i]:
                res.append((self.auxBlocks[i].get_bool_array_used_packets()))

        return res

    def getHalfPacketListFromPacket(self, packet: RU10Packet) -> List[List[bool]]:
        """
        Generates a list of halfpackets from a packet.
        :param packet: The packet to get the list from
        :return: List of halfpackets
        """
        res: List[List[bool]] = []
        aux_used_packets = packet.get_bool_array_half_packets()
        for i in range(len(aux_used_packets)):
            if aux_used_packets[i]:
                res.append(
                    (
                        self.auxBlocks[
                            packet.get_number_of_ldpc_blocks() + i
                        ].get_bool_array_used_and_ldpc_packets()
                    )
                )
        return res

    def solve(self):
        if self.use_headerchunk and self.headerChunk is None:
            self.decodeHeader()
        # Decoder.solve(self)
        super(self.__class__, self).solve()

    def is_decoded(self) -> bool:
        return self.getSolvedCount() >= self.number_of_chunks

    def getSolvedCount(self) -> int:
        return len(self.decodedPackets) + (1 if self.headerChunk is not None else 0)

    def getNextValidPacket(
        self,
        from_multiple_files: bool = False,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
    ) -> Optional[RU10Packet]:
        """
        Takes a raw packet from a file and calls @parse_raw_packet to get a RU10 packet. If the packet is corrupt the
        next one will be taken.
        :param from_multiple_files: True: The packets were saved in multiple files. False: Packets were saved in one file.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: RU10Packet
        """
        assert self.f is not None
        if not from_multiple_files:
            packet_len = self.f.read(struct.calcsize("<" + packet_len_format))
            try:
                packet_len = struct.unpack("<" + packet_len_format, packet_len)[0]
                packet = self.f.read(int(packet_len))
            except struct.error:
                return None
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF = True
            self.f.close()
            return None
        res = self.parse_raw_packet(
            packet,
            crc_len_format=crc_len_format,
            packet_len_format=packet_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            id_len_format=id_len_format,
        )
        if res is None:
            res = self.getNextValidPacket(
                from_multiple_files,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                id_len_format=id_len_format,
            )
        return res

    def parse_raw_packet(
        self,
        packet: bytes,
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "L",
        packet_len_format: str = "I",
        id_len_format: str = "L",
    ) -> Optional[RU10Packet]:
        """
        Creates a RU10 packet from a raw given packet. Also checks if the packet is corrupted. If any method was used to
        create packets from specific chunks, set self.use_method = True. This will treat the last byte of the raw packet
        data as the byte that contains the information about the used method ("even", "odd", "window_30 + window" or
        "window_40 + window". See RU10Encoder.create_new_packet_from_chunks for further information.
        :param packet: A raw packet
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: RU10Packet or an error message
        """
        struct_str = "<" + number_of_chunks_len_format + id_len_format
        struct_len = struct.calcsize(struct_str)
        """"
        if self.error_correction.__code__.co_name == crc32.__code__.co_name:
            crc_len = -struct.calcsize("<" + crc_len_format)
            payload = packet[:crc_len]
            crc = struct.unpack("<" + crc_len_format, packet[crc_len:])[0]
            calced_crc = calc_crc(payload)
            if crc != calced_crc:  # If the Packet is corrupt, try next one
                logger.warning("CRC-Error - %s != %s", hex(crc), hex(calced_crc))
                self.corrupt += 1
                return "CORRUPT"

        else:
        """
        try:
            packet = self.error_correction(packet)
        except (AssertionError, ReedSolomonError, ValueError):
            self.corrupt += 1
            return None

        data = packet[struct_len:]
        data, chunk_lst = self._packet_method_data(data)
        len_data = struct.unpack(struct_str, packet[0:struct_len])
        unxored_id = self._update_packet_header(
            len_data, number_of_chunks_len_format, id_len_format
        )
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
        return res

    def generateIntermediateBlocksFormat(self, number_of_chunks: int) -> List[List[List[int]]]:
        """
        Generates the format of the intermediate blocks from the number of used chunks.
        :param number_of_chunks: The number of used chunks.
        :return:
        """
        compositions: List[List[int]] = [[] for _ in range(self.s)]
        for i in range(0, number_of_chunks):
            a = 1 + (int(floor(np.float64(i) / np.float64(self.s))) % (self.s - 1))
            b = int(i % self.s)
            compositions[b].append(i)
            b = (b + a) % self.s
            compositions[b].append(i)
            b = (b + a) % self.s
            compositions[b].append(i)

        hprime: int = int(ceil(np.float64(self.h) / 2))
        m = buildGraySequence(number_of_chunks + self.s, hprime)
        hcompositions: List[List[int]] = [[] for _ in range(self.h)]
        for i in range(0, self.h):
            hcomposition = []
            for j in range(0, number_of_chunks + self.s):
                gray_row = m[j]
                gray_value = (
                    int(gray_row.item()) if isinstance(gray_row, np.ndarray) else int(gray_row)
                )
                if bitSet(gray_value, i):
                    hcomposition.append(j)
            hcompositions[i] = hcomposition
        res = [compositions, hcompositions]
        return res

    def decodeHeader(self, last_chunk_len_format: str = "I") -> None:
        if self.headerChunk is not None or 1 not in self.degreeToPacket.keys():
            return  # Header already set
        for decoded in self.degreeToPacket[1]:
            if decoded.get_used_packets().issubset({0}):
                self.headerChunk = HeaderChunk(decoded, last_chunk_len_format=last_chunk_len_format)
                return

    def decodeZip(self, *args: Any, **kwargs: Any) -> Optional[Any]:
        raise RuntimeError("Not implemented for BP-based decoders in the current version!")

    def saveDecodedFile(
        self,
        last_chunk_len_format: str = "I",
        null_is_terminator: bool = False,
        print_to_output: bool = True,
        return_file_name: bool = False,
    ) -> Union[bytes, str]:
        """
        Saves the file - if decoded. The filename is either taken from the headerchunk or generated based on the input
        filename.
        :param return_file_name: if set to true, this function will return the filename under which the file as been saved
        :param last_chunk_len_format: Format of the last chunk length
        :param null_is_terminator: True: The file is handled as null-terminated C-String.
        :param print_to_output: True: Result we be printed to the command line.
        :return:
        """
        assert self.is_decoded(), "Can not save File: Unable to reconstruct."
        if self.use_headerchunk:
            self.decodeHeader()
        file_name = self._resolve_output_file_name()
        output_concat = b""
        with open(file_name, "wb") as f:
            for decoded in sorted(self.decodedPackets):
                [num] = decoded.get_used_packets()
                if not self._should_write_decoded_packet(num):
                    continue
                output_bytes, stop_writing = self._packet_output_bytes(
                    decoded, num, null_is_terminator
                )
                output_concat += output_bytes
                f.write(output_bytes)
                if stop_writing:
                    break
        logger.info("Saved file as '%s'", file_name)
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))
        if return_file_name:
            return file_name
        return output_concat

    def mode_1_bmp_decode(self, last_chunk_len_format: str = "I"):
        dec_out = self.saveDecodedFile(
            last_chunk_len_format=last_chunk_len_format,
            null_is_terminator=False,
            print_to_output=False,
        )
        assert isinstance(dec_out, bytes), "saveDecodedFile did not return a bytes object!"
        return self.bytes_to_bitmap(dec_out)

    def bytes_to_bitmap(self, img_byt: bytes):
        width, height = (
            int(struct.unpack(">H", img_byt[:2])[0]),
            int(struct.unpack(">H", img_byt[2:4])[0]),
        )
        unpack = (
            np.unpackbits(
                np.frombuffer(img_byt, dtype=np.uint8, count=int((width * height) / 8), offset=4)
            )
            .reshape(height, width)
            .transpose()
        )
        flip_bits: IntArray = np.logical_not(unpack).astype(int)
        new_img = self.draw_img(flip_bits, width, height)
        assert self.file is not None, "filename must be known, not None!"
        tmp_file_name = os.path.basename(self.file) + ".bmp"
        file_name = "DEC_" + tmp_file_name if self.file is not None else "RU10.BIN.bmp"
        new_img.save(file_name)
        return file_name

    @staticmethod
    def draw_img(unpacked_flipped_bits: IntArray, width: int, height: int) -> Image.Image:
        new_img = Image.new("1", (width, height))
        pixels = new_img.load()
        if pixels is None:
            raise RuntimeError("Failed to load image pixel buffer")

        for i in range(new_img.size[0]):
            for j in range(new_img.size[1]):
                pixels[i, j] = int(unpacked_flipped_bits[i, j])
        return new_img


def main(
    file: str,
    num_of_chunks: int,
    err_correction: Callable[[bytes], bytes] = nocode,
    insert_header: bool = False,
    mode_1_bmp: bool = False,
):
    x = RU10BPDecoder(
        file,
        use_headerchunk=insert_header,
        error_correction=err_correction,
        static_number_of_chunks=num_of_chunks,
    )
    x.decode(id_len_format="I", number_of_chunks_len_format="I")
    x.saveDecodedFile(null_is_terminator=False, print_to_output=False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file / folder to Decode")
    parser.add_argument(
        "--error_correction",
        metavar="error_correction",
        type=str,
        required=False,
        default="nocode",
        help="Error Correction Method to use; possible values: \
                                    nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)",
    )
    parser.add_argument("--insert_header", required=False, action="store_true", default=False)
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=True, type=int)
    parser.add_argument(
        "--repair_symbols",
        metavar="repair_symbols",
        type=int,
        required=False,
        default=2,
        help="number of repair symbols for ReedSolomon (default=2)",
    )
    parser.add_argument("--as_mode_1_bmp", required=False, action="store_true")
    args = parser.parse_args()
    _file = args.filename
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _mode_1_bmp = args.as_mode_1_bmp
    _number_of_chunks = args.number_of_chunks
    _error_correction = get_error_correction_decode(args.error_correction, _repair_symbols)
    logger.info("File / Folder to decode: %s", _file)
    main(_file, _number_of_chunks, _error_correction, _insert_header, _mode_1_bmp)
    logger.info("Decoding finished.")
