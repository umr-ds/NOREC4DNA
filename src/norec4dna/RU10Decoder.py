#!/usr/bin/python
# -*- coding: latin-1 -*-
from __future__ import annotations

import argparse
import io
import logging
import os
import struct
import typing
from configparser import SectionProxy
from io import BytesIO
from math import ceil, floor
from zipfile import ZipFile

import numpy as np
from PIL import Image
from typing_extensions import Callable

from .Decoder import Decoder
from .distributions.Distribution import Distribution
from .distributions.RaptorDistribution import RaptorDistribution
from .ErrorCorrection import get_error_correction_decode, nocode
from .GEPP import GEPP, GEPP_intern
from .HeaderChunk import HeaderChunk
from .helper import bitSet, buildGraySequence, calc_file_crc, logical_xor, xor_mask
from .helper.quaternary2Bin import (
    quad_file_to_bytes,
    quat_file_to_bin,
    tranlate_quat_to_byte,
)
from .helper.RU10Helper import (
    from_true_false_list,
    intermediate_symbols,
)
from .Packet import Packet
from .RU10IntermediatePacket import RU10IntermediatePacket
from .RU10Packet import RU10Packet
from .RU10Shared import RU10Shared

DEBUG = False
logger = logging.getLogger(__name__)

BoolArray = np.ndarray[typing.Any, np.dtype[np.bool_]]


class RU10Decoder(RU10Shared, Decoder):
    def __init__(
        self,
        file: typing.Optional[str] = None,
        error_correction: Callable[[bytes], bytes] = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: typing.Optional[int] = None,
        use_method: bool = False,
        checksum_len_str: typing.Optional[str] = None,
        xor_by_seed: bool = False,
        mask_id: bool = True,
        id_spacing: int = 0,
        config_map: typing.Optional[typing.Any] = None,
    ):
        self.debug = False
        super().__init__()
        if checksum_len_str is None:
            self.checksum_len_str = ""
        if not use_headerchunk and (checksum_len_str != "" and checksum_len_str is not None):
            logger.warning(
                "[Warning] Header-checksums are only supported with headerchunks! Checksum from config file will be ignored!"
            )
        self.checksum_len_str = checksum_len_str
        self.isPseudo: bool = False
        self.file: typing.Optional[str] = file
        self.degreeToPacket: dict = {}
        self.use_method: bool = use_method
        self.xor_by_seed = xor_by_seed
        self.mask_id = mask_id
        if self.file is not None:
            self.isFolder = os.path.isdir(self.file)
            self.isZip = self.file.endswith(".zip")
            if not self.isFolder:
                self.f = open(self.file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        if static_number_of_chunks is not None:
            self.number_of_chunks = static_number_of_chunks
        self.headerChunk: typing.Optional[HeaderChunk] = None
        self.GEPP: typing.Optional[GEPP_intern] = None
        self.pseudoCount: int = 0
        self.ldpcANDhalf: typing.Dict[int, RU10IntermediatePacket] = {}
        self.repairBlockNumbers: dict = {}
        self.s: int = -1
        self.h: int = -1
        self.distribution: typing.Optional[Distribution] = None
        self.EOF: bool = False
        self.counter: typing.Dict[int, int] = {}
        self.count: bool = True
        self.error_correction: typing.Callable[[bytes], bytes] = error_correction
        self.use_headerchunk: bool = use_headerchunk
        self.static_number_of_chunks: typing.Optional[int] = static_number_of_chunks
        self.id_spacing = id_spacing
        self.packets = []
        self.config_map = config_map

    def _update_number_of_chunks_format(self, number_of_chunks_len_format: str) -> str:
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            return ""
        return number_of_chunks_len_format

    def _open_dna_stream(self, filepath: str) -> bool:
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

    def _open_folder_packet(self, filename: str) -> bool:
        if self.file is None:
            raise ValueError("self.file must be set before opening packets")
        filepath = self.file + "/" + filename
        if filename.endswith("DNA"):
            return self._open_dna_stream(filepath)
        self.f = open(filepath, "rb")
        return True

    def _decode_input_packet(
        self,
        from_multiple_files: bool,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> typing.Optional[RU10Packet]:
        new_pack = self.getNextValidPacket(
            from_multiple_files,
            packet_len_format=packet_len_format,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            id_len_format=id_len_format,
        )
        if new_pack is None or new_pack == "CORRUPT":
            return None
        return new_pack

    def _record_packet(self, packet: RU10Packet) -> bool:
        self.packets.append(packet)
        return self.input_new_packet(packet)

    def _decode_zip_entry(
        self,
        archive: ZipFile,
        name: str,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> bool:
        self.f = io.BytesIO(archive.read(name))
        new_pack = self._decode_input_packet(
            True,
            packet_len_format,
            crc_len_format,
            number_of_chunks_len_format,
            id_len_format,
        )
        if hasattr(self, "f"):
            self.f.close()
        return new_pack is not None and self.input_new_packet(new_pack)

    def _sorted_zip_namelist(self, archive: ZipFile) -> list[str]:
        namelist = archive.namelist()
        try:
            split_names = [name.split("_") for name in namelist]
            sorted_names = sorted(split_names, key=lambda parts: float(parts[1]), reverse=False)
            return [parts[0] + "_" + parts[1] for parts in sorted_names]
        except Exception:
            return namelist

    def _folder_candidate_files(self) -> list[str]:
        if self.file is None:
            raise ValueError("self.file must be set for decodeFolder")
        return [
            filename
            for filename in os.listdir(self.file)
            if filename.endswith(".RU10") or filename.endswith("DNA")
        ]

    def _decode_folder_entry(
        self,
        filename: str,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> bool:
        self.EOF = False
        if not self._open_folder_packet(filename):
            return False
        new_pack = self._decode_input_packet(
            True,
            packet_len_format,
            crc_len_format,
            number_of_chunks_len_format,
            id_len_format,
        )
        return new_pack is not None and self.input_new_packet(new_pack)

    def _open_file_input(self) -> None:
        if self.file is None:
            raise ValueError("self.file must be set for decodeFile")
        if self.file.lower().endswith("dna") and not self._open_dna_stream(self.file):
            return

    def _open_fasta_input(self) -> None:
        if self.file is None:
            raise ValueError("self.file must be set for decodeFile")
        self.f.close()
        self.f = open(self.file, "r")

    def _read_fasta_entry(self) -> typing.Optional[tuple[str, str, str]]:
        line = self.f.readline()
        if not line:
            self.EOF = True
            return None
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        try:
            error_prob, seed = line[1:].replace("\n", "").split("_")
        except ValueError:
            error_prob, seed = "0", "0"
        line = self.f.readline()
        if not line:
            self.EOF = True
            return None
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        return error_prob, seed, line.replace("\n", "")

    def _decode_fasta_packet(
        self,
        dna_str: str,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> typing.Union[RU10Packet, str, None]:
        reverted_dna_str = self.revert_seed_spacing(dna_str, id_len_format)
        return self.parse_raw_packet(
            BytesIO(tranlate_quat_to_byte(reverted_dna_str)).read(),
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            packet_len_format=packet_len_format,
            id_len_format=id_len_format,
            dna_str=dna_str,
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
            # Strip appended version fields if enabled
            dna_str_stripped = self.strip_appended_version_fields(dna_str)
            try:
                new_pack = self._decode_fasta_packet(
                    dna_str_stripped,
                    packet_len_format,
                    crc_len_format,
                    number_of_chunks_len_format,
                    id_len_format,
                )
            except Exception:
                new_pack = "CORRUPT"
                self.packets.append(dna_str_stripped)
            if new_pack != "CORRUPT":
                assert isinstance(new_pack, RU10Packet), "new_pack must be RU10Packet"
                decoded = self._record_packet(new_pack)
                if self.progress_bar is not None:
                    self.progress_bar.update(self.correct, Corrupt=self.corrupt)
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
            new_pack = self._decode_input_packet(
                False,
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
            if new_pack is None:
                break
            decoded = self._record_packet(new_pack)
        return decoded

    def _can_solve_now(self) -> bool:
        return (
            self.GEPP is not None
            and self.GEPP.isPotentionallySolvable()
            and not self.read_all_before_decode
        )

    def _finalize_decode(self, decoded: bool, eof_message: str) -> typing.Union[int, bool]:
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets: %s", self.corrupt)
        if self.GEPP is None:
            logger.warning("No Packet was correctly decoded. Check your configuration.")
            return -1
        if self._can_solve_now():
            decoded = self.GEPP.solve()
        if not decoded and self.EOF:
            logger.warning(eof_message)
            return -1
        return decoded

    def _initialize_distribution(self, packet: RU10Packet) -> None:
        if len(self.ldpcANDhalf) != 0 or self.distribution is not None:
            return
        self.distribution = RaptorDistribution(self.number_of_chunks)
        self.number_of_chunks = packet.get_total_number_of_chunks()
        _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.distribution)
        self.createAuxBlocks()
        self.progress_bar = self.create_progress_bar(
            int(self.number_of_chunks + 0.02 * self.number_of_chunks)
        )

    def _update_counter(self, removed: BoolArray) -> None:
        if not self.count:
            return
        for i in range(len(removed)):
            if i in self.counter.keys():
                if removed[i]:
                    self.counter[i] += 1
            else:
                self.counter[i] = 1

    def _store_removed_packet(self, removed: BoolArray, packet: RU10Packet) -> None:
        packet_data = np.frombuffer(packet.get_data(), dtype="uint8")
        if self.GEPP is None:
            self.GEPP = GEPP(np.array([removed], dtype=bool), packet_data)
            return
        if self.GEPP.b is not None and self.GEPP.b.size > 0 and len(packet_data) != self.GEPP.b.shape[-1]:
            return
        self.GEPP.addRow(np.array(removed, dtype=bool), packet_data)




    def _update_packet_header(
        self,
        len_data: tuple[typing.Any, ...],
        number_of_chunks_len_format: str,
        id_len_format: str,
    ) -> int:
        if self.static_number_of_chunks is None:
            self.number_of_chunks = xor_mask(len_data[0], number_of_chunks_len_format)
            return xor_mask(len_data[1], id_len_format, enabled=self.mask_id)
        return xor_mask(len_data[0], id_len_format, enabled=self.mask_id)

    def _ensure_distribution_ready(self) -> None:
        if self.distribution is None:
            self.distribution = RaptorDistribution(self.number_of_chunks)
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.distribution)
            self.progress_bar = self.create_progress_bar(
                int(self.number_of_chunks + 0.02 * self.number_of_chunks)
            )
        if self.correct == 0:
            self.createAuxBlocks()
        self.correct += 1


    def _resolve_output_file_name(self) -> str:
        file_name = "DEC_" + os.path.basename(self.file) if self.file is not None else "RU10.BIN"
        if self.headerChunk is None:
            return file_name.split("\x00")[0]
        try:
            header_file_name = self.headerChunk.get_file_name()
            resolved = (
                header_file_name.decode("utf-8")
                if isinstance(header_file_name, bytes)
                else header_file_name
            )
            return resolved.split("\x00")[0]
        except Exception as ex:
            logger.warning("%s", ex)
            return file_name.split("\x00")[0]

    def _write_output_chunk(
        self, gepp: GEPP_intern, x: int, null_is_terminator: bool
    ) -> tuple[bytes, bool, bool]:
        if x < 0:
            return b"\x00" * len(gepp.b[x][0]), False, True
        if self.number_of_chunks - 1 == x and self.use_headerchunk:
            assert self.headerChunk is not None
            output = gepp.b[x][0][0 : self.headerChunk.get_last_chunk_length()]
            output_bytes = output if isinstance(output, bytes) else output.tobytes()
            return output_bytes, False, False
        if null_is_terminator:
            splitter = gepp.b[x].tobytes().decode().split("\x00")
            return splitter[0].encode(), len(splitter) > 1, False
        output = gepp.b[x]
        output_bytes = output if isinstance(output, bytes) else output.tobytes()
        return output_bytes, False, False

    def _validate_checksum(self, file_name: str, ignore_crc: bool) -> None:
        if self.checksum_len_str is None or self.checksum_len_str == "":
            return
        decoded_crc = calc_file_crc(file_name, self.checksum_len_str)
        assert self.headerChunk is not None
        if self.headerChunk.checksum != decoded_crc:
            logger.warning("Decoded CRC: %s", decoded_crc)
            logger.warning("Header CRC: %s", self.headerChunk.checksum)
            if not ignore_crc:
                raise ValueError(
                    "Checksum of decoded file does not match checksum in header chunk!",
                    file_name,
                )










    @staticmethod
    def from_config_map(config_map: SectionProxy) -> "RU10Decoder":
        """
        missing (set elsewhere :
        id_len_format = decode_conf.get("id_len_format")
        number_of_chunks_len_format = decode_conf.get("number_of_chunks_len_format", "I")
        crc_len_format = decode_conf.get("crc_len_format", "L")
        """
        return RU10Decoder(
            file=config_map.name,
            error_correction=get_error_correction_decode(
                config_map.get("error_correction", "nocode"), config_map.getint("repair_symbols", 2)
            ),
            use_headerchunk=config_map.getboolean("insert_header", True),
            static_number_of_chunks=config_map.getint("number_of_chunks", None),
            use_method=config_map.getboolean(
                "method", False
            ),  # default to False as it is rarely used.
            checksum_len_str=config_map.get("checksum_len_str", ""),
            xor_by_seed=config_map.getboolean("xor_by_seed", False),
            mask_id=config_map.getboolean("mask_id", True),
            id_spacing=config_map.getint("id_spacing", 0),
            config_map=config_map,
        )

    def decodeZip(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "I",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
        store_parsed_packets: bool = False,
        *args,
        **kwargs,
    ):
        if hasattr(self, "f"):
            self.f.close()
        decoded = False
        self.EOF = False
        if self.file is None:
            raise ValueError("self.file must be set for decodeZip")
        number_of_chunks_len_format = self._update_number_of_chunks_format(
            number_of_chunks_len_format
        )
        archive = ZipFile(self.file, "r")
        for name in self._sorted_zip_namelist(archive):
            decoded = self._decode_zip_entry(
                archive,
                name,
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
            if decoded:
                break
            if self.progress_bar is not None:
                self.progress_bar.update(self.correct, Corrupt=self.corrupt)
        return self._finalize_decode(
            decoded, "Unable to retrieve File from Chunks. Too many errors?"
        )

    def decodeFolder(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "I",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
        store_parsed_packets: bool = False,
        *args,
        **kwargs,
    ):
        """
        Decodes the information from a folder if self.file represents a folder and the packets were saved
        in multiple files and prints the number of decoded and corrupted packets.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: -1 if the decoding wasn't successful
        """
        decoded = False
        self.EOF = False
        if self.file is None:
            raise ValueError("self.file must be set for decodeFolder")
        number_of_chunks_len_format = self._update_number_of_chunks_format(
            number_of_chunks_len_format
        )
        for file_in_folder in self._folder_candidate_files():
            decoded = self._decode_folder_entry(
                file_in_folder,
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
            if decoded:
                break
            if self.progress_bar is not None:
                self.progress_bar.update(self.correct, Corrupt=self.corrupt)
        if hasattr(self, "f"):
            self.f.close()
        return self._finalize_decode(
            decoded, "Unable to retrieve File from Chunks. Too many errors?"
        )



    def decodeFile(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
        store_parsed_packets=False,
        *args,
        **kwargs,
    ):
        """
        Decodes the information from a file if self.file represents a file and the packets were saved in a single file.
        :param packet_len_format: Format of the packet length
        :param crc_len_format:  Format of the crc length
        :param number_of_chunks_len_format: Format of the number of chunks length
        :param id_len_format: Format of the ID length
        :return: -1 if the decoding wasn't successful
        """
        decoded = False
        self.EOF = False
        if self.file is None:
            raise ValueError("self.file must be set for decodeFile")
        self._open_file_input()
        number_of_chunks_len_format = self._update_number_of_chunks_format(
            number_of_chunks_len_format
        )
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
        if self._can_solve_now():
            assert self.GEPP is not None
            decoded = self.GEPP.solve()
        if not decoded and self.EOF and not self.read_all_before_decode:
            logger.warning("Unable to retrieve file from chunks. Too many errors??")
            return -1
        return decoded
        # self.f.close()

    def getNumberOfLDPCBlocks(self):
        return self.s

    def getNumberOfHalfBlocks(self):
        return self.h

    def getNumberOfRepairBlocks(self):
        return self.getNumberOfHalfBlocks() + self.getNumberOfLDPCBlocks()

    def input_new_packet(self, packet: RU10Packet) -> bool:
        """
        Removes auxpackets (LDPC and Half) and adds the remaining data to the GEPP matrix.
        :param packet: A Packet to add to the GEPP matrix
        :return: True: If solved. False: Else.
        """
        self._initialize_distribution(packet)
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
        self._update_counter(removed)
        self._store_removed_packet(removed, packet)
        assert self.GEPP is not None
        if (
            self.isPseudo or not self.read_all_before_decode
        ) and self.GEPP.isPotentionallySolvable():
            # and self.GEPP.n % 5 == 0:  # Nur alle 5 Packete versuch starten
            if self.debug:
                logger.debug("current size: %s", self.GEPP.n)
            return self.GEPP.solve(partial=False)
        return False

    # Correct
    def removeAndXorAuxPackets(self, packet: RU10Packet) -> BoolArray:
        """
        Removes auxpackets (LDCP and Half) from a given packet to get the packets data.
        :param packet: Packet to remove auxpackets from
        :return: The data without the auxpackets
        """
        aux_mapping = self.getHalfPacketListFromPacket(packet)  # Enthaelt Data + LDPC Nummern
        aux_mapping.append(packet.get_bool_array_used_and_ldpc_packets())
        xored_list = logical_xor(aux_mapping)
        del aux_mapping
        tmp = from_true_false_list(xored_list)  # Nur noch Data + LDPC sind vorhanden
        if self.debug:
            logger.debug("%s", tmp)
        tmp = RU10Packet("", tmp, self.number_of_chunks, packet.id, packet.dist, read_only=True)
        aux_mapping = self.getAuxPacketListFromPacket(tmp)
        bool_used = tmp.get_bool_array_used_packets()
        if bool_used is not None:
            aux_mapping.append(bool_used)  # [-len(self.auxBlocks):])
        else:
            logger.warning("Aux packets for current packet returned None!")
        res = logical_xor(aux_mapping)
        del tmp, aux_mapping
        return res


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
                "We should have "
                + str(self.getNumberOfLDPCBlocks())
                + " LDPC-Blocks, "
                + str(self.getNumberOfHalfBlocks())
                + " Half-Blocks and "
                + str(self.number_of_chunks)
                + " normal Chunks (including 1 HeaderChunk)"
            )
        for i in range(0, self.getNumberOfRepairBlocks()):
            self.repairBlockNumbers[i] = set()
        i = 0
        for group in self.generateIntermediateBlocksFormat(self.number_of_chunks):
            for elem in group:
                self.repairBlockNumbers[i] = elem
                i += 1
        # XOR all Chunks into the corresponding AUX-Block
        for aux_number in self.repairBlockNumbers.keys():
            self.ldpcANDhalf[aux_number] = RU10IntermediatePacket(
                "",
                self.repairBlockNumbers[aux_number],
                total_number_of_chunks=self.number_of_chunks,
                id=aux_number,
                dist=self.distribution,
            )
            if self.debug:
                logger.debug("%s : %s", aux_number, self.ldpcANDhalf[aux_number].used_packets)

    # Correct
    def getAuxPacketListFromPacket(self, packet: RU10Packet) -> typing.List[typing.List[bool]]:
        """
        Creates a list for a packet with information about whether auxpackets have been used for that packet.
        :param packet: The packet to check.
        :return: Information about used auxpackets.
        """
        res: typing.List[typing.List[bool]] = []
        aux_used_packets = packet.get_bool_array_repair_packets()
        for i in range(len(aux_used_packets)):
            if aux_used_packets[i]:
                tmp = self.ldpcANDhalf[i].get_bool_array_used_packets()
                if tmp is not None:
                    res.append(tmp)
                else:
                    logger.warning("Aux-List was None!")
        return res

    def getHalfPacketListFromPacket(self, packet: RU10Packet) -> typing.List[typing.List[bool]]:
        """
        Generates a list of halfpackets from a packet.
        :param packet: The packet to get the list from
        :return: List of halfpackets
        """
        res: typing.List[typing.List[bool]] = []
        aux_used_packets = packet.get_bool_array_half_packets()
        for i in range(len(aux_used_packets)):
            if aux_used_packets[i]:
                res.append(
                    (
                        self.ldpcANDhalf[
                            packet.get_number_of_ldpc_blocks() + i
                        ].get_bool_array_used_and_ldpc_packets()
                    )
                )
        return res

    def solve(self, partial=False) -> bool:
        """
        Calls GEPP.solve()
        :return: True: If GEPP was able to solve the matrix. False: Else.
        """
        assert self.GEPP is not None, "GEPP must not be None!"
        return self.GEPP.solve(partial=partial)

    def getSolvedCount(self) -> int:
        if self.GEPP is None:
            return 0
        return self.GEPP.getSolvedCount()

    def is_decoded(self) -> bool:
        """
        Checks if the data is decoded.
        :return: True: If the decoding was successfull. False: Else.
        """
        return (
            self.GEPP is not None and self.GEPP.isPotentionallySolvable() and self.GEPP.isSolved()
        )

    def getNextValidPacket(
        self,
        from_multiple_files: bool = False,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
    ) -> typing.Optional[RU10Packet]:
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
        if not from_multiple_files:
            packet_len_bytes = typing.cast(
                bytes, self.f.read(struct.calcsize("<" + packet_len_format))
            )
            try:
                packet_len = struct.unpack("<" + packet_len_format, packet_len_bytes)[0]
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
        # cast to bytes if packet is of type str:
        if isinstance(packet, str):
            packet = packet.encode("utf-8")
        res = self.parse_raw_packet(
            packet,
            crc_len_format=crc_len_format,
            packet_len_format=packet_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            id_len_format=id_len_format,
        )
        if res == "CORRUPT" and not from_multiple_files:
            return self.getNextValidPacket(
                from_multiple_files,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                id_len_format=id_len_format,
            )
        if res == "CORRUPT":
            return None
        return res


    def generateIntermediateBlocksFormat(
        self, number_of_chunks: int
    ) -> typing.List[typing.List[typing.List[int]]]:
        """
        Generates the format of the intermediate blocks from the number of used chunks.
        :param number_of_chunks: The number of used chunks.
        :return:
        """
        compositions: typing.List[typing.List[int]] = [[] for _ in range(self.s)]
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
        hcompositions: typing.List[typing.List[int]] = [[] for _ in range(self.h)]
        for i in range(0, self.h):
            hcomposition = []
            for j in range(0, number_of_chunks + self.s):
                if bitSet(int(m[j]), int(i)):
                    hcomposition.append(j)
            hcompositions[i] = hcomposition
        res = [compositions, hcompositions]
        return res

    def populate_header_chunk(self, last_chunk_len_str: typing.Optional[str] = None):
        if last_chunk_len_str is None:
            if self.config_map is None:
                last_chunk_len_str = "I"
            else:
                last_chunk_len_str = self.config_map.get("last_chunk_len_str", "I")
        assert last_chunk_len_str is not None
        assert self.GEPP is not None, "GEPP must be set before populating Header!"
        if self.use_headerchunk:
            header_row = self.GEPP.result_mapping[0]
            if header_row >= 0:
                self.headerChunk = HeaderChunk(
                    Packet(
                        self.GEPP.b[header_row].tobytes(),
                        {0},
                        self.number_of_chunks,
                        read_only=True,
                    ),
                    last_chunk_len_format=last_chunk_len_str,
                    checksum_len_format=self.checksum_len_str or "",
                )

    def saveDecodedFile(
        self,
        last_chunk_len_format: str = "I",
        null_is_terminator: bool = False,
        print_to_output: bool = True,
        return_file_name=False,
        partial_decoding: bool = True,
        ignore_crc=False,
    ) -> typing.Union[bytes, str]:
        """
        Saves the file - if decoded. The filename is either taken from the headerchunk or generated based on the input
        filename.
        :param partial_decoding: perform partial decoding if full decoding failed, missing parts will be filled with "\x00"
        :param return_file_name: if set to true, this function will return the filename under which the file as been saved
        :param last_chunk_len_format: Format of the last chunk length
        :param null_is_terminator: True: The file is handled as null-terminated C-String.
        :param print_to_output: True: Result we be printed to the command line.
        :return:
        """
        assert (
            self.is_decoded() or partial_decoding
        ), "Can not save File: Unable to reconstruct. You may try saveDecodedFile(partial_decoding=True)"
        if partial_decoding:
            self.solve(partial=True)
        dirty = False
        self.populate_header_chunk(last_chunk_len_str=last_chunk_len_format)
        if self.GEPP is None:
            raise RuntimeError("GEPP not initialized")
        gepp = self.GEPP
        file_name = self._resolve_output_file_name()
        output_concat = b""
        with open(file_name, "wb") as f:
            for x in gepp.result_mapping:
                if x == 0 and self.use_headerchunk:
                    continue
                output_bytes, stop_writing, marked_dirty = self._write_output_chunk(
                    gepp, x, null_is_terminator
                )
                output_concat += output_bytes
                f.write(output_bytes)
                dirty = dirty or marked_dirty
                if stop_writing:
                    break
        logger.info("Saved file as '%s'", file_name)
        self._validate_checksum(file_name, ignore_crc)
        if dirty:
            logger.warning(
                "Some parts could not be restored, file WILL contain sections with \\x00 !"
            )
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))
        if self.progress_bar is not None:
            self.progress_bar.update(self.number_of_chunks, Corrupt=self.corrupt)
        if return_file_name:
            return file_name
        return output_concat

    def mode_1_bmp_decode(self, last_chunk_len_format: str = "I"):
        dec_out = self.saveDecodedFile(
            last_chunk_len_format=last_chunk_len_format,
            null_is_terminator=False,
            print_to_output=False,
        )
        return self.bytes_to_bitmap(typing.cast(bytes, dec_out))

    def bytes_to_bitmap(self, img_byt: bytes):
        width, height = struct.unpack(">H", img_byt[:2])[0], struct.unpack(">H", img_byt[2:4])[0]
        unpack = (
            np.unpackbits(
                np.frombuffer(img_byt, dtype=np.uint8, count=int((width * height) / 8), offset=4)
            )
            .reshape(height, width)
            .transpose()
        )
        flip_bits = np.logical_not(unpack).astype(int)
        new_img = self.draw_img(flip_bits, width, height)
        tmp_file_name = (
            os.path.basename(self.file) + ".bmp" if self.file is not None else "RU10.BIN.bmp"
        )
        file_name = "DEC_" + tmp_file_name if self.file is not None else "RU10.BIN.bmp"
        new_img.save(file_name)
        return file_name

    @staticmethod
    def draw_img(
        unpacked_flipped_bits: np.ndarray[typing.Any, np.dtype[typing.Any]], width: int, height: int
    ) -> "Image.Image":
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
    number_of_chunks: int,
    error_correction: typing.Callable = nocode,
    insert_header: bool = False,
    mode_1_bmp: bool = False,
    _header_crc_str: typing.Optional[str] = None,
    xor_by_seed=False,
    _id_spacing=0,
):
    logger.info("Pure Gauss-Mode")
    x = RU10Decoder(
        file,
        use_headerchunk=insert_header,
        error_correction=error_correction,
        static_number_of_chunks=number_of_chunks,
        checksum_len_str=_header_crc_str,
        xor_by_seed=xor_by_seed,
        id_spacing=_id_spacing,
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
    parser.add_argument(
        "--header_crc_str", metavar="header_crc_str", required=False, type=str, default=""
    )
    parser.add_argument(
        "--as_mode_1_bmp",
        required=False,
        action="store_true",
        help="convert to a header-less B/W BMP format. (use only for image/bmp input)",
    )
    parser.add_argument("--xor_by_seed", required=False, action="store_true")
    parser.add_argument("--id_spacing", metavar="id_spacing", required=False, type=int, default=0)
    args = parser.parse_args()
    _file = args.filename
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _mode_1_bmp = args.as_mode_1_bmp
    _number_of_chunks = args.number_of_chunks
    _error_correction = get_error_correction_decode(args.error_correction, _repair_symbols)
    _header_crc_str = args.header_crc_str
    _xor_by_seed = args.xor_by_seed
    _id_spacing = args.id_spacing
    logger.info("File / Folder to decode: %s", _file)
    main(
        _file,
        _number_of_chunks,
        _error_correction,
        _insert_header,
        _mode_1_bmp,
        _header_crc_str,
        _xor_by_seed,
        _id_spacing,
    )
    logger.info("Decoding finished.")