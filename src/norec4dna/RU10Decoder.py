#!/usr/bin/python
# -*- coding: latin-1 -*-
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
from numpy.typing import NDArray
from PIL import Image
from typing_extensions import Callable

from .Decoder import Decoder
from .distributions.Distribution import Distribution
from .distributions.RaptorDistribution import RaptorDistribution
from .ErrorCorrection import get_error_correction_decode, nocode
from .GEPP import GEPP, GEPP_intern
from .HeaderChunk import HeaderChunk
from .helper import bitSet, buildGraySequence, calc_file_crc, logical_xor, xor_mask
from .helper.helper import xor_with_seed
from .helper.quaternary2Bin import (
    quad_file_to_bytes,
    quat_file_to_bin,
    tranlate_quat_to_byte,
)
from .helper.RU10Helper import (
    choose_packet_numbers,
    from_true_false_list,
    intermediate_symbols,
)
from .Packet import Packet
from .RU10IntermediatePacket import RU10IntermediatePacket
from .RU10Packet import RU10Packet

DEBUG = False
logger = logging.getLogger(__name__)


class RU10Decoder(Decoder):
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
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        if self.file is None:
            raise ValueError("self.file must be set for decodeZip")
        archive = ZipFile(self.file, "r")
        namelist = archive.namelist()
        try:
            nam = [x.split("_") for x in namelist]
            sorted_by_second = sorted(nam, key=lambda tup: float(tup[1]), reverse=False)
            namelist = [x[0] + "_" + x[1] for x in sorted_by_second]
        except Exception:
            pass
        for name in namelist:
            self.f = io.BytesIO(archive.read(name))
            new_pack = self.getNextValidPacket(
                True,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                id_len_format=id_len_format,
            )
            if hasattr(self, "f"):
                self.f.close()
            if new_pack is None:
                break
                # koennte durch input_new_packet ersetzt werden:
                # self.addPacket(new_pack)
                if new_pack != "CORRUPT":
                    decoded = self.input_new_packet(new_pack)
                else:
                    if DEBUG:
                        logger.debug("Packet with name=%s corrupt.", name)
            if decoded:
                break
            if self.progress_bar is not None:
                self.progress_bar.update(self.correct, Corrupt=self.corrupt)
            ##
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets: %s", self.corrupt)

        if self.GEPP is None:
            logger.warning("No Packet was correctly decoded. Check your configuration.")
            return -1
        if (
            self.GEPP is not None
            and self.GEPP.isPotentionallySolvable()
            and not self.read_all_before_decode
        ):
            decoded = self.GEPP.solve()
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve File from Chunks. Too many errors?")
            return -1
        return decoded

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
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        for file_in_folder in os.listdir(self.file):
            if file_in_folder.endswith(".RU10") or file_in_folder.endswith("DNA"):
                self.EOF = False
                if file_in_folder.endswith("DNA"):
                    if self.error_correction.__name__ == "dna_reed_solomon_decode":
                        try:
                            self.f = quad_file_to_bytes(self.file + "/" + file_in_folder)
                        except TypeError:
                            logger.warning("skipping CORRUPT file - contains illegal character(s)")
                            self.corrupt += 1
                            continue
                    else:
                        try:
                            self.f = quat_file_to_bin(self.file + "/" + file_in_folder)
                        except TypeError:
                            logger.warning("skipping CORRUPT file - contains illegal character(s)")
                            self.corrupt += 1
                            continue
                else:
                    self.f = open(self.file + "/" + file_in_folder, "rb")
                new_pack = self.getNextValidPacket(
                    True,
                    packet_len_format=packet_len_format,
                    crc_len_format=crc_len_format,
                    number_of_chunks_len_format=number_of_chunks_len_format,
                    id_len_format=id_len_format,
                )
                if new_pack is not None and new_pack != "CORRUPT":
                    # koennte durch input_new_packet ersetzt werden:
                    # self.addPacket(new_pack)
                    decoded = self.input_new_packet(new_pack)
                if decoded:
                    break
            if self.progress_bar is not None:
                self.progress_bar.update(self.correct, Corrupt=self.corrupt)
            ##
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets: %s", self.corrupt)
        if hasattr(self, "f"):
            self.f.close()
        if self.GEPP is None:
            logger.warning("No Packet was correctly decoded. Check your configuration.")
            return -1
        if (
            self.GEPP is not None
            and self.GEPP.isPotentionallySolvable()
            and not self.read_all_before_decode
        ):
            decoded = self.GEPP.solve()
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve File from Chunks. Too many errors?")
            return -1
        return decoded

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
        if self.file.lower().endswith("dna"):
            try:
                self.f.close()
                self.f = quat_file_to_bin(self.file)
            except TypeError:
                logger.warning("skipping CORRUPT file - contains illegal character(s)")
                self.corrupt += 1
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        if self.file.lower().endswith("fasta"):
            self.f.close()
            self.f = open(self.file, "r")
            raw_packet_list = []
            while not (decoded or self.EOF):
                line = self.f.readline()
                if not line:
                    self.EOF = True
                    break
                try:
                    error_prob, seed = line[1:].replace("\n", "").split("_")
                except:
                    error_prob, seed = "0", "0"
                line = self.f.readline()
                if not line:
                    self.EOF = True
                    break
                dna_str = line.replace("\n", "")
                # un-space the dna string:
                reverted_dna_str = self.revert_seed_spacing(dna_str, id_len_format)

                raw_packet_list.append((error_prob, seed, dna_str))
                try:
                    new_pack = self.parse_raw_packet(
                        BytesIO(tranlate_quat_to_byte(reverted_dna_str)).read(),
                        crc_len_format=crc_len_format,
                        number_of_chunks_len_format=number_of_chunks_len_format,
                        packet_len_format=packet_len_format,
                        id_len_format=id_len_format,
                        dna_str=dna_str,
                    )

                except Exception:
                    new_pack = "CORRUPT"
                    # TODO: as the metadata trick might invalidate the checksum, we do not want to throw this packet away!
                    # after reverting the metadata changes, we may parse the DNA string by calling:
                    # self.input_new_packets(self.parse_raw_packet(BytesIO(translate_quat_to_byte(repaired_dna)...)
                    # (separate and check for exception)!
                    self.packets.append(dna_str)
                if new_pack != "CORRUPT":
                    assert isinstance(new_pack, RU10Packet), "new_pack must be RU10Packet"
                    decoded = self.input_new_packet(new_pack)
                    self.packets.append(new_pack)
                    if self.progress_bar is not None:
                        self.progress_bar.update(self.correct, Corrupt=self.corrupt)
        else:
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
                # koennte durch input_new_packet ersetzt werden:
                # self.addPacket(new_pack)
                if new_pack != "CORRUPT":
                    self.packets.append(new_pack)
                    decoded = self.input_new_packet(new_pack)
                #
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets : %s", self.corrupt)
        if (
            self.GEPP is not None
            and self.GEPP.isPotentionallySolvable()
            and not self.read_all_before_decode
        ):
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
        if len(self.ldpcANDhalf) == 0 and self.distribution is None:  # self.isPseudo and
            self.distribution = RaptorDistribution(self.number_of_chunks)
            self.number_of_chunks = packet.get_total_number_of_chunks()
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.distribution)
            self.createAuxBlocks()
            self.progress_bar = self.create_progress_bar(
                int(self.number_of_chunks + 0.02 * self.number_of_chunks)
            )
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
        if self.count:
            for i in range(len(removed)):
                if i in self.counter.keys():
                    if removed[i]:
                        self.counter[i] += 1
                else:
                    self.counter[i] = 1
        if self.GEPP is None:
            self.GEPP = GEPP(
                np.array([removed], dtype=bool),
                np.frombuffer(packet.get_data(), dtype="uint8"),
            )
        else:
            self.GEPP.addRow(
                np.array(removed, dtype=bool),
                np.frombuffer(packet.get_data(), dtype="uint8"),
            )
        if (
            self.isPseudo or not self.read_all_before_decode
        ) and self.GEPP.isPotentionallySolvable():
            # and self.GEPP.n % 5 == 0:  # Nur alle 5 Packete versuch starten
            if self.debug:
                logger.debug("current size: %s", self.GEPP.n)
            return self.GEPP.solve(partial=False)
        return False

    # Correct
    def removeAndXorAuxPackets(self, packet: RU10Packet) -> NDArray[np.bool_]:
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

    def removeAndXorAuxPackets_from_indices(
        self, packet_indices: typing.Union[typing.Set[int], typing.List[int]]
    ) -> np.ndarray:
        """
        Removes auxpackets (LDPC and Half) from a list/set of packet indices to get the chunk composition.

        This is the equivalent of removeAndXorAuxPackets but operates on packet indices directly
        without requiring a RU10Packet object.

        Packet index interpretation:
        - Indices 0 to number_of_chunks-1: systematic packets (each contains exactly one chunk)
        - Indices number_of_chunks to number_of_chunks+s-1: LDPC packets (indexed in ldpcANDhalf as 0 to s-1)
        - Indices number_of_chunks+s to number_of_chunks+s+h-1: Half packets (indexed in ldpcANDhalf as s to s+h-1)

        Args:
            self: Decoder instance with ldpcANDhalf, distribution, and number_of_chunks initialized
            packet_indices: Set/list/array of packet indices to process

        Returns:
            Boolean numpy array where index i=True means chunk i is in the result after aux removal
        """
        import numpy as np

        from .helper.helper import logical_xor
        from .helper.RU10Helper import from_true_false_list

        # Convert input to set of indices
        if isinstance(packet_indices, np.ndarray):
            if packet_indices.dtype == bool:
                # Boolean array - extract True indices
                packet_indices = from_true_false_list(packet_indices.tolist())
            # else: integer array - use as-is
        packet_set = set(packet_indices)

        if not packet_set:
            return np.zeros(self.number_of_chunks, dtype=bool)

        # Separate packet types
        systematic_indices = set()
        ldpc_indices = set()
        half_indices = set()

        for idx in packet_set:
            if 0 <= idx < self.number_of_chunks:
                # Systematic packet - directly represents a chunk
                systematic_indices.add(idx)
            elif self.number_of_chunks <= idx < self.number_of_chunks + self.s:
                # LDPC packet
                ldpc_indices.add(idx - self.number_of_chunks)
            elif self.number_of_chunks + self.s <= idx < self.number_of_chunks + self.s + self.h:
                # Half packet
                half_indices.add(idx - self.number_of_chunks - self.s)

        # Build result starting with systematic packets (each represents one chunk)
        result = np.zeros(self.number_of_chunks, dtype=bool)
        for chunk_idx in systematic_indices:
            result[chunk_idx] = True

        # Process Half packets first (they contain data + LDPC)
        # Track newly discovered LDPC indices from half packet contents
        new_ldpc_from_half = set()

        if half_indices:
            half_list = []
            half_arr_size = 0

            for half_idx in half_indices:
                ldpc_half_idx = self.s + half_idx  # Index in ldpcANDhalf
                if ldpc_half_idx in self.ldpcANDhalf:
                    half_arr = self.ldpcANDhalf[
                        ldpc_half_idx
                    ].get_bool_array_used_and_ldpc_packets()
                    half_arr_size = max(half_arr_size, len(half_arr))
                    half_list.append(half_arr)

            if half_list:
                # Pad all arrays to the same size
                padded_half_list = []
                for arr in half_list:
                    if len(arr) < half_arr_size:
                        arr_padded = np.zeros(half_arr_size, dtype=bool)
                        arr_padded[: len(arr)] = arr
                        padded_half_list.append(arr_padded)
                    else:
                        padded_half_list.append(arr)

                # XOR all half packets together
                half_result = logical_xor(padded_half_list)

                # Add systematic AND original LDPC packets to the XOR (to match removeAndXorAuxPackets behavior)
                systematic_and_ldpc_arr = np.zeros(half_arr_size, dtype=bool)
                for chunk_idx in systematic_indices:
                    if chunk_idx < half_arr_size:
                        systematic_and_ldpc_arr[chunk_idx] = True
                # Also mark original LDPC indices
                for ldpc_idx in ldpc_indices:
                    if self.number_of_chunks + ldpc_idx < half_arr_size:
                        systematic_and_ldpc_arr[self.number_of_chunks + ldpc_idx] = True

                half_list_with_sys = [half_result, systematic_and_ldpc_arr]
                combined = logical_xor(half_list_with_sys)

                # Extract data chunks and LDPC indices from combined result
                data_and_ldpc = from_true_false_list(combined)
                result = np.zeros(self.number_of_chunks, dtype=bool)
                for idx in data_and_ldpc:
                    if 0 <= idx < self.number_of_chunks:
                        result[idx] = True
                    elif self.number_of_chunks <= idx < self.number_of_chunks + self.s:
                        # LDPC index from half packet content - track for processing
                        new_ldpc_from_half.add(idx - self.number_of_chunks)

        # Process LDPC packets:
        # - If half packets existed: only process NEW LDPC discovered from half contents (original LDPC already XORed)
        # - If no half packets: process all original LDPC
        ldpc_to_process = new_ldpc_from_half if half_indices else ldpc_indices
        if ldpc_to_process:
            aux_list = []
            for ldpc_idx in ldpc_to_process:
                if ldpc_idx in self.ldpcANDhalf:
                    aux_arr = self.ldpcANDhalf[ldpc_idx].get_bool_array_used_packets()
                    # Ensure correct size
                    if aux_arr is not None and len(aux_arr) < self.number_of_chunks:
                        aux_arr_padded = np.zeros(self.number_of_chunks, dtype=bool)
                        aux_arr_padded[: len(aux_arr)] = aux_arr
                        aux_list.append(aux_arr_padded)
                    elif aux_arr is not None:
                        aux_list.append(aux_arr)

            if aux_list:
                # XOR all LDPC packets together with current result
                aux_result = logical_xor(aux_list)
                combined_list = [aux_result, result]
                result = logical_xor(combined_list)

        return result

    def createAuxBlocks(self):
        """
        Reconstructs the auxblocks to be able to remove them afterwards.
        :return:
        """
        assert self.number_of_chunks is not None, (
            "createAuxBlocks can only be called AFTER first Packet"
        )
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
        return self.GEPP.solve(partial=partial)

    def getSolvedCount(self) -> int:
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
            except:
                return None
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF = True
            try:
                self.f.close()
            except:
                return None
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
            res = self.getNextValidPacket(  # type: ignore[assignment]
                from_multiple_files,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                id_len_format=id_len_format,
            )
        return res  # type: ignore[return-value]

    def parse_raw_packet(
        self,
        packet_input: bytes,
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "L",
        packet_len_format: str = "I",
        id_len_format: str = "L",
        dna_str: typing.Optional[str] = None,
    ) -> typing.Union[RU10Packet, str]:
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
        except:
            self.corrupt += 1
            return "CORRUPT"
        header = typing.cast(bytes, packet)[:struct_len]
        data = typing.cast(bytes, packet)[struct_len:]
        chunk_lst = []
        if self.use_method:
            method_data = bin(data[-1])[2:]
            while len(method_data) < 8:
                method_data = "0" + method_data
            data = data[:-1]
            if method_data.startswith("00"):
                chunk_lst = [ch for ch in range(0, self.number_of_chunks + 1) if ch % 2 == 0]
            elif method_data.startswith("01"):
                chunk_lst = [ch for ch in range(0, self.number_of_chunks + 1) if ch % 2 != 0]
            elif method_data.startswith("10"):
                window = int(method_data[2:], 2)
                window_size = 30
                start = window * (window_size - 10)
                chunk_lst = [
                    ch for ch in range(start, start + window_size) if ch <= self.number_of_chunks
                ]
            elif method_data.startswith("11"):
                window = int(method_data[2:], 2)
                window_size = 40
                start = window * (window_size - 10)
                chunk_lst = [
                    ch for ch in range(start, start + window_size) if ch <= self.number_of_chunks
                ]
            else:
                raise RuntimeError("Not a valid start:", method_data)
        len_data = struct.unpack(struct_str, header)
        if self.static_number_of_chunks is None:
            self.number_of_chunks = xor_mask(len_data[0], number_of_chunks_len_format)
            unxored_id = xor_mask(len_data[1], id_len_format, enabled=self.mask_id)
        else:
            unxored_id = xor_mask(len_data[0], id_len_format, enabled=self.mask_id)
        if self.xor_by_seed:
            data = xor_with_seed(data, unxored_id)
        if self.distribution is None:
            self.distribution = RaptorDistribution(self.number_of_chunks)
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.distribution)
            self.progress_bar = self.create_progress_bar(
                int(self.number_of_chunks + 0.02 * self.number_of_chunks)
            )

        if self.correct == 0:
            self.createAuxBlocks()
        self.correct += 1
        if self.use_method:
            numbers = choose_packet_numbers(
                len(chunk_lst),
                unxored_id,
                typing.cast(RaptorDistribution, self.distribution),  # type: ignore[arg-type]
                systematic=False,
                max_l=len(chunk_lst),
            )
            used_packets = [chunk_lst[i] for i in numbers]
        else:
            used_packets = choose_packet_numbers(
                self.number_of_chunks,
                unxored_id,
                typing.cast(RaptorDistribution, self.distribution),
                systematic=False,  # type: ignore[arg-type]
            )
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
        assert self.is_decoded() or partial_decoding, (
            "Can not save File: Unable to reconstruct. You may try saveDecodedFile(partial_decoding=True)"
        )
        if partial_decoding:
            self.solve(partial=True)
        dirty = False
        self.populate_header_chunk(last_chunk_len_str=last_chunk_len_format)
        file_name = "DEC_" + os.path.basename(self.file) if self.file is not None else "RU10.BIN"
        output_concat = b""
        if self.headerChunk is not None:
            try:
                file_name = self.headerChunk.get_file_name().decode("utf-8")
            except Exception as ex:
                logger.warning("%s", ex)
        file_name = file_name.split("\x00")[0]
        assert self.GEPP is not None
        with open(file_name, "wb") as f:
            for x in self.GEPP.result_mapping:
                if x < 0:
                    f.write(b"\x00" * len(self.GEPP.b[x][0]))
                    dirty = True
                    continue
                if 0 != x or not self.use_headerchunk:
                    if self.number_of_chunks - 1 == x and self.use_headerchunk:
                        assert self.headerChunk is not None
                        output = self.GEPP.b[x][0][0 : self.headerChunk.get_last_chunk_length()]
                        output_concat += output.tobytes()
                        f.write(output)
                    else:
                        if null_is_terminator:
                            splitter = self.GEPP.b[x].tostring().decode().split("\x00")
                            output = splitter[0].encode()
                            if type(output) == bytes:
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
                            if len(splitter) > 1:
                                break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                        else:
                            output = self.GEPP.b[x]
                            if isinstance(output, bytes):
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
        logger.info("Saved file as '%s'", file_name)
        if self.checksum_len_str is not None and self.checksum_len_str != "":
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
        unpacked_flipped_bits: NDArray[typing.Any], width: int, height: int
    ) -> "Image.Image":
        new_img = Image.new("1", (width, height))
        pixels = new_img.load()

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
