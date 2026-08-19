from __future__ import annotations

import argparse
import io
import logging
import os
import struct
from configparser import SectionProxy
from io import BytesIO
from math import ceil, floor
from typing import Any, Callable, Dict, List, Optional, Set, Union
from zipfile import ZipFile

import numpy as np
from PIL import Image

from .BPDecoder import BPDecoder
from .distributions.Distribution import Distribution
from .distributions.RaptorDistribution import RaptorDistribution
from .ErrorCorrection import get_error_correction_decode, nocode
from .GEPP import GEPP, GEPP_intern
from .HeaderChunk import HeaderChunk
from .helper import bitSet, buildGraySequence, calc_file_crc, xor_mask
from .helper.quaternary2Bin import quad_file_to_bytes, quat_file_to_bin, tranlate_quat_to_byte
from .helper.RU10Helper import intermediate_symbols
from .Packet import Packet
from .RU10IntermediatePacket import RU10IntermediatePacket
from .RU10Packet import RU10Packet
from .RU10Shared import RU10Shared

logger = logging.getLogger(__name__)

IntArray = np.ndarray[Any, np.dtype[np.int_]]


class RU10BPDecoder(RU10Shared, BPDecoder):
    def __init__(
        self,
        file: Optional[str] = None,
        error_correction: Callable[..., Any] = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: Optional[int] = None,
        use_method: bool = False,
        checksum_len_str: Optional[str] = None,
        xor_by_seed: bool = False,
        mask_id: bool = True,
        id_spacing: int = 0,
        config_map: Optional[Any] = None,
    ):
        super().__init__(
            file=file,
            error_correction=error_correction,
            use_headerchunk=use_headerchunk,
            static_number_of_chunks=static_number_of_chunks,
            use_method=use_method,
        )
        if checksum_len_str is None:
            self.checksum_len_str = ""
        else:
            self.checksum_len_str = checksum_len_str
        self.xor_by_seed = xor_by_seed
        self.mask_id = mask_id
        self.id_spacing = id_spacing
        self.config_map = config_map
        self.GEPP: Optional[GEPP_intern] = None
        self.packets: List[RU10Packet] = []
        self._gepp_finalized: bool = False
        self.file: Optional[str] = file
        self.use_method: bool = use_method
        self.f: Optional[Any] = None
        self.auxBlocks: dict[int, RU10IntermediatePacket] = {}
        self.ldpcANDhalf: dict[int, RU10IntermediatePacket] = self.auxBlocks
        self.isFolder: bool = False
        if file is not None:
            self.isFolder = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(file, "rb")
        self.number_of_chunks: int = 1000000
        if static_number_of_chunks is not None:
            self.number_of_chunks = static_number_of_chunks
        self.s: int = -1
        self.h: int = -1
        self.error_correction: Callable[..., Any] = error_correction
        self.use_headerchunk: bool = use_headerchunk
        self.static_number_of_chunks: Optional[int] = static_number_of_chunks
        self.dist: Optional[Distribution] = None
        self.counter: Dict[int, int] = {}
        self.count: bool = True
        self._solved_cache: set[int] = set()
        self._solved_cache_len: int = -1

    @property
    def distribution(self) -> Optional[Distribution]:
        return self.dist

    @distribution.setter
    def distribution(self, val: Optional[Distribution]) -> None:
        self.dist = val

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
            return file_name.split("\x00")[0]
        try:
            header_file_name = self.headerChunk.get_file_name()
            raw_name = (
                header_file_name.decode("utf-8", errors="replace")
                if isinstance(header_file_name, bytes)
                else str(header_file_name)
            )
            resolved_clean = raw_name.split("\x00")[0].strip()
            return resolved_clean if resolved_clean else file_name.split("\x00")[0]
        except Exception as ex:
            logger.warning("%s", ex)
            return file_name.split("\x00")[0]

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
            # Strip appended version fields if enabled (mirrors RU10Decoder)
            dna_str_stripped = self.strip_appended_version_fields(dna_str)
            new_pack = self._decode_fasta_packet(
                dna_str_stripped,
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                id_len_format,
            )
            if new_pack is not None:
                decoded = self.input_new_packet(new_pack)
        if not decoded:
            decoded = self.solve(partial=True)
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
        if not decoded:
            decoded = self.solve(partial=True)
        return decoded


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
        if not number_of_chunks_len_format:
            return xor_mask(len_data[0], id_len_format, enabled=self.mask_id)
        if self.static_number_of_chunks is None:
            self.number_of_chunks = xor_mask(len_data[0], number_of_chunks_len_format)
            return xor_mask(len_data[1], id_len_format, enabled=self.mask_id)
        return xor_mask(len_data[0], id_len_format, enabled=self.mask_id)

    def _ensure_distribution_ready(self) -> None:
        if self.dist is None:
            self.dist = RaptorDistribution(self.number_of_chunks)
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.dist)
        if self.correct == 0:
            self.createAuxBlocks()
        self.correct += 1


    def decodeFolder(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "I",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
        store_parsed_packets: bool = False,
        **kwargs: Any,
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
        store_parsed_packets: bool = False,
        **kwargs: Any,
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

    def _remove_and_xor_aux_packets_set(self, packet: RU10Packet) -> Set[int]:
        """O(degree) aux removal returning the set of data-chunk indices.

        Half packets reference data + LDPC indices; LDPC blocks reference data
        indices only.  Used by the peeling hot path; the boolean-mask
        ``removeAndXorAuxPackets`` is kept for the multi-version layer.
        """
        used_set: Set[int] = set(packet.get_used_packets() or [])
        ldpc_offset = self.number_of_chunks
        # Step 1: remove half packets (they reference data + LDPC indices).
        for i, used_half in enumerate(packet.get_bool_array_half_packets()):
            if used_half:
                block = self.auxBlocks[self.s + i]
                used_set.symmetric_difference_update(block.get_used_packets() or [])
        # Step 2: remove the LDPC blocks still referenced *after* half removal.
        ldpc_used = [idx - ldpc_offset for idx in used_set if ldpc_offset <= idx < ldpc_offset + self.s]
        for ldpc_idx in ldpc_used:
            block = self.auxBlocks[ldpc_idx]
            used_set.symmetric_difference_update(block.get_used_packets() or [])
        return {u for u in used_set if u < self.number_of_chunks}

    def removeAndXorAuxPackets(self, packet: RU10Packet) -> List[bool]:
        """
        Removes auxpackets (LDCP and Half) from a given packet to get the packets data.
        :param packet: Packet to remove auxpackets from
        :return: The data without the auxpackets
        """
        used_set = self._remove_and_xor_aux_packets_set(packet)
        mask = np.zeros(self.number_of_chunks, dtype=bool)
        for u in used_set:
            mask[u] = True
        return mask.tolist()

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
        used_set = self._remove_and_xor_aux_packets_set(packet)
        packet.set_used_packets(used_set)

        packet_data = packet.get_data()
        if self.GEPP is not None and self.GEPP.b is not None and self.GEPP.b.size > 0:
            expected_len = self.GEPP.b.shape[-1]
            if len(packet_data) != expected_len:
                return False

        self.packets.append(packet)
        self._store_removed_packet(used_set, packet)
        return self.updatePackets(packet)

    def _store_removed_packet(self, removed: Any, packet: RU10Packet) -> None:
        # BP peeling does not need the accumulated packet matrix: solved chunks
        # are tracked in solved_symbols / decodedPackets, and the chunk-indexed
        # GEPP is materialised once at the end (see _materialize_gepp_b).  We
        # only keep a width placeholder so callers can read GEPP.b.shape[-1]
        # (chunk size) during feeding.
        packet_data = np.frombuffer(packet.get_data(), dtype="uint8")
        self._update_counter(removed)
        if self.GEPP is None or (
            self.GEPP.b is not None and self.GEPP.b.shape[-1] != len(packet_data)
        ):
            mask = np.zeros(self.number_of_chunks, dtype=bool)
            for u in removed:
                if u < self.number_of_chunks:
                    mask[u] = True
            self.GEPP = GEPP(np.array([mask], dtype=bool), packet_data)

    def _update_counter(self, removed: Any) -> None:
        if not self.count:
            return
        for i in removed:
            self.counter[i] = self.counter.get(i, 0) + 1

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

    def _sync_gepp_to_decoded_packets(self) -> None:
        """If GEPP has solved rows, copy any missing chunks into self.decodedPackets."""
        if self.GEPP is None or getattr(self.GEPP, "result_mapping", None) is None:
            return
        existing_chunks: set[int] = set()
        for pkt in self.decodedPackets:
            existing_chunks.update(pkt.get_used_packets())

        b_arr = self.GEPP.b
        res_map = self.GEPP.result_mapping
        for i in range(self.number_of_chunks):
            if i not in existing_chunks and i < len(res_map):
                row_idx = int(res_map[i][0]) if getattr(res_map[i], "ndim", 0) > 0 else int(res_map[i])
                if 0 <= row_idx < len(b_arr):
                    row_data = b_arr[row_idx]
                    pkt_bytes = row_data.tobytes() if isinstance(row_data, np.ndarray) else bytes(row_data)
                    pkt = RU10Packet(
                        pkt_bytes,
                        {i},
                        self.number_of_chunks,
                        i,
                        read_only=True,
                    )
                    self.decodedPackets.add(pkt)

    def _inactivation_solve(self, partial: bool = False) -> None:
        """Solve symbols the BP peeling could not resolve via a residual GEPP.

        After belief propagation stalls, every received packet has been XOR-
        reduced: symbols solved during peeling were removed from its
        ``used_packets`` and XORed out of its data, so each remaining packet
        references only unsolved symbols and its data is the XOR of exactly
        those symbols.  Building a GEPP system over the *unsolved columns only*
        and solving it recovers the remaining chunks (RFC 5053 "inactivation
        decoding").  Restricting the system to the unsolved symbol set keeps
        the solve O(r·u²) instead of O(r·K²), which is near-linear whenever
        peeling resolves the bulk of the pool.
        """
        if self.GEPP is None:
            return
        unsolved = [
            i for i in range(self.number_of_chunks) if i not in self.solved_symbols
        ]
        if not unsolved:
            return
        col_of = {sym: c for c, sym in enumerate(unsolved)}

        # Every received packet (after peeling) is an equation over the
        # remaining unsolved symbols.  ``symbol_to_packets`` only ever contains
        # packets that are already in ``self.packets``, so iterating
        # ``self.packets`` alone covers every equation (no dedup needed).
        residual_rows: List[np.ndarray] = []
        residual_data: List[np.ndarray] = []
        for pkt in self.packets:
            used = pkt.get_used_packets()
            if not used:
                continue
            cols = [col_of[u] for u in used if u in col_of]
            if not cols:
                continue
            mask = np.zeros(len(unsolved), dtype=bool)
            mask[cols] = True
            
            d_raw = pkt.get_data()
            if isinstance(d_raw, (bytes, bytearray)):
                d_arr = np.frombuffer(d_raw, dtype=np.uint8).copy()
            else:
                d_arr = np.asarray(d_raw, dtype=np.uint8).copy()
            # XOR out all symbols that were solved during BP peeling
            for u in used:
                if u in self.solved_symbols and u not in col_of:
                    s_obj = self.solved_symbols[u]
                    s_data = s_obj.get_data() if hasattr(s_obj, "get_data") else s_obj
                    if isinstance(s_data, (bytes, bytearray)):
                        s_arr = np.frombuffer(s_data, dtype=np.uint8)
                    else:
                        s_arr = np.asarray(s_data, dtype=np.uint8)
                    n = min(len(d_arr), len(s_arr))
                    if n:
                        d_arr[:n] ^= s_arr[:n]
            residual_rows.append(mask)
            residual_data.append(d_arr)
        if not residual_rows:
            return
        rgepp = GEPP(np.array(residual_rows, dtype=bool), np.array(residual_data, dtype=np.uint8))
        try:
            rgepp.solve(partial=True)
        except Exception:
            return
        res_map = rgepp.result_mapping
        for c, sym in enumerate(unsolved):
            if sym in self.solved_symbols or c >= len(res_map):
                continue
            row_idx = int(res_map[c][0]) if getattr(res_map[c], "ndim", 0) > 0 else int(res_map[c])
            if 0 <= row_idx < len(rgepp.b):
                data = rgepp.b[row_idx]
                pkt = RU10Packet(
                    bytes(data),
                    {sym},
                    self.number_of_chunks,
                    sym,
                    read_only=True,
                )
                self.solved_symbols[sym] = pkt
                self.decodedPackets.add(pkt)
                self.queue.append(pkt)

    def _materialize_gepp_b(self) -> None:
        """Make ``GEPP.b`` chunk-indexed so the multi-version layer can read it.

        The DR4DNA multi-version coder/decoder reads chunk data straight from
        ``GEPP.b[chunk_id]`` (and via ``result_mapping``) rather than through
        the BP solver's symbol bookkeeping.  After solving, rebuild ``GEPP`` as
        an identity system whose row ``i`` holds the content of chunk ``i``.
        Chunks that are NOT yet solved keep ``result_mapping == -1`` so callers
        can tell partial from complete reconstruction.
        """
        if self.GEPP is None:
            return
        n = self.number_of_chunks
        csize = int(self.GEPP.b.shape[-1])
        new_b = np.zeros((n, csize), dtype=np.uint8)
        res_map = np.full((n, 1), -1, dtype=np.int32)
        for i in range(n):
            obj = self.solved_symbols.get(i)
            if obj is None:
                continue
            data = obj.get_data() if hasattr(obj, "get_data") else obj
            if isinstance(data, (bytes, bytearray)):
                row = np.frombuffer(bytes(data), dtype=np.uint8)
            else:
                row = np.asarray(data, dtype=np.uint8).ravel()
            new_b[i, : min(csize, row.size)] = row[:csize]
            res_map[i, 0] = i
        self.GEPP.A = np.identity(n, dtype=bool)
        self.GEPP.b = new_b
        self.GEPP.result_mapping = res_map
        self.GEPP.packet_mapping = np.arange(n, dtype=np.intp)
        self.GEPP.n = n
        self.GEPP.m = n
        self.GEPP.chunk_to_used_packets = np.identity(n, dtype=bool)

    def _finalize_gepp(self) -> None:
        """Materialise the chunk-indexed GEPP once the solver has finished."""
        self._materialize_gepp_b()
        self._sync_gepp_to_decoded_packets()
        self._gepp_finalized = True

    def _peel(self) -> bool:
        # Near-linear peeling; materialise GEPP as soon as the pool is solved so
        # callers reading GEPP.b/result_mapping always see correct data.
        res = super()._peel()
        if res and not self._gepp_finalized:
            self._finalize_gepp()
        return res

    def solve(self, partial: bool = False) -> bool:
        # 1. Peeling (belief propagation): solve as many symbols as possible.
        self._peel()
        # 2. Inactivation decoding: GEPP on the residual (unsolved) system.
        if not self.is_decoded():
            self._inactivation_solve(partial=partial)
            self._peel()
        # 3. Materialize a chunk-indexed GEPP for the multi-version layer.
        self._finalize_gepp()
        if self.use_headerchunk and self.headerChunk is None:
            self.decodeHeader()
        return self.is_decoded()

    def is_decoded(self) -> bool:
        return len(self.getSolvedChunkIds()) >= self.number_of_chunks

    def getSolvedCount(self) -> int:
        return len(self.getSolvedChunkIds())

    def getSolvedChunkIds(self) -> set[int]:
        """Return the set of chunk ids that have been uniquely solved."""
        dec_len = len(self.decodedPackets)
        sym_len = len(self.solved_symbols) if hasattr(self, "solved_symbols") and self.solved_symbols else 0
        hdr_flag = self.headerChunk is not None
        gepp_res = self.GEPP.result_mapping if (self.GEPP is not None and getattr(self.GEPP, "result_mapping", None) is not None) else None
        gepp_id = id(gepp_res) if gepp_res is not None else None

        cache_state = (dec_len, sym_len, hdr_flag, gepp_id)
        if getattr(self, "_solved_cache_state", None) == cache_state:
            return self._solved_cache

        chunks: set[int] = set()
        for pkt in self.decodedPackets:
            chunks.update(pkt.get_used_packets())
        if sym_len > 0:
            chunks.update(self.solved_symbols.keys())
        if hdr_flag:
            chunks.add(0)
        if gepp_res is not None:
            for cid, row in enumerate(gepp_res):
                row_idx = int(row[0]) if getattr(row, "ndim", 0) > 0 else int(row)
                if row_idx >= 0:
                    chunks.add(cid)
        self._solved_cache = chunks
        self._solved_cache_state = cache_state
        return chunks

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
        packet_input: bytes,
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "L",
        packet_len_format: str = "I",
        id_len_format: str = "L",
        dna_str: Optional[str] = None,
    ) -> Optional[RU10Packet]:
        res = super().parse_raw_packet(
            packet_input,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            packet_len_format=packet_len_format,
            id_len_format=id_len_format,
            dna_str=dna_str,
        )
        if res == "CORRUPT" or not isinstance(res, RU10Packet):
            return None
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

    @staticmethod
    def from_config_map(config_map: SectionProxy) -> "RU10BPDecoder":
        return RU10BPDecoder(
            file=config_map.name,
            error_correction=get_error_correction_decode(
                config_map.get("error_correction", "nocode"), config_map.getint("repair_symbols", 2)
            ),
            use_headerchunk=config_map.getboolean("insert_header", True),
            static_number_of_chunks=config_map.getint("number_of_chunks", None),
            use_method=config_map.getboolean("method", False),
            checksum_len_str=config_map.get("checksum_len_str", ""),
            xor_by_seed=config_map.getboolean("xor_by_seed", False),
            mask_id=config_map.getboolean("mask_id", True),
            id_spacing=config_map.getint("id_spacing", 0),
            config_map=config_map,
        )

    def populate_header_chunk(self, last_chunk_len_str: Optional[str] = None):
        if last_chunk_len_str is None:
            if self.config_map is None:
                last_chunk_len_str = "I"
            else:
                last_chunk_len_str = self.config_map.get("last_chunk_len_str", "I")
        assert last_chunk_len_str is not None
        if self.headerChunk is None:
            self.decodeHeader(last_chunk_len_format=last_chunk_len_str)
        if self.headerChunk is None and 0 in self.decodedPackets:
            pkt = next(p for p in self.decodedPackets if p.get_used_packets() == {0})
            self.headerChunk = HeaderChunk(
                pkt,
                last_chunk_len_format=last_chunk_len_str,
                checksum_len_format=self.checksum_len_str or "",
            )
        if self.headerChunk is None and self.GEPP is not None:
            try:
                if not self.GEPP.isSolved():
                    self.GEPP.solve(partial=True)
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
            except Exception:
                pass











    def decodeHeader(self, last_chunk_len_format: str = "") -> None:
        if not last_chunk_len_format:
            if self.config_map is None:
                last_chunk_len_format = "I"
            else:
                last_chunk_len_format = str(self.config_map.get("last_chunk_len_str", "I"))
        if (
            self.headerChunk is not None
            and self.headerChunk.last_chunk_len_format == last_chunk_len_format
        ):
            return  # Header already parsed with the requested format
        if getattr(self, "solved_symbols", None) is not None and 0 in self.solved_symbols:
            self.headerChunk = HeaderChunk(
                self.solved_symbols[0],
                last_chunk_len_format=last_chunk_len_format,
                checksum_len_format=self.checksum_len_str or "",
            )
            return
        for pkt in self.decodedPackets:
            if pkt.get_used_packets() == {0}:
                self.headerChunk = HeaderChunk(
                    pkt,
                    last_chunk_len_format=last_chunk_len_format,
                    checksum_len_format=self.checksum_len_str or "",
                )
                return
        if self.GEPP is not None and getattr(self.GEPP, "result_mapping", None) is not None:
            res_map = self.GEPP.result_mapping
            if 0 < len(res_map):
                row_idx = int(res_map[0][0]) if getattr(res_map[0], "ndim", 0) > 0 else int(res_map[0])
                if 0 <= row_idx < len(self.GEPP.b):
                    row_data = self.GEPP.b[row_idx]
                    pkt_bytes = row_data.tobytes() if isinstance(row_data, np.ndarray) else bytes(row_data)
                    pkt = RU10Packet(pkt_bytes, {0}, self.number_of_chunks, 0, read_only=True)
                    self.headerChunk = HeaderChunk(
                        pkt,
                        last_chunk_len_format=last_chunk_len_format,
                        checksum_len_format=self.checksum_len_str or "",
                    )
                    return

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
        new_pack = self.getNextValidPacket(
            True,
            packet_len_format=packet_len_format,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            id_len_format=id_len_format,
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

    def decodeZip(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "I",
        number_of_chunks_len_format: str = "I",
        id_len_format: str = "I",
        store_parsed_packets: bool = False,
        *args: Any,
        **kwargs: Any,
    ) -> Optional[int]:
        if hasattr(self, "f") and self.f is not None:
            self.f.close()
        decoded = False
        self.EOF = False
        if self.file is None:
            raise ValueError("self.file must be set for decodeZip")
        number_of_chunks_len_format = self._update_decode_number_format(
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
        archive.close()
        if not self.is_decoded():
            self.solve(partial=True)
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets : %s", self.corrupt)
        if not self.is_decoded() and self.EOF:
            logger.warning("Unable to retrieve File from Chunks. Too many errors?")
            return -1
        return 1 if self.is_decoded() else 0

    def saveDecodedFile(
        self,
        last_chunk_len_format: str = "I",
        null_is_terminator: bool = False,
        print_to_output: bool = True,
        return_file_name: bool = False,
        partial_decoding: bool = True,
        ignore_crc: bool = False,
    ) -> Union[bytes, str]:
        """
        Saves the file - if decoded. The filename is either taken from the headerchunk or generated based on the input
        filename.
        :param partial_decoding: perform partial decoding if full decoding failed, missing parts will be filled with "\\x00"
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
        if self.use_headerchunk:
            self.decodeHeader(last_chunk_len_format=last_chunk_len_format)
        file_name = self._resolve_output_file_name()
        output_concat = b""
        with open(file_name, "wb") as f:
            chunk_size = getattr(self, "chunk_size", 0) or (
                self.GEPP.b.shape[1] if self.GEPP is not None and hasattr(self.GEPP, "b") and self.GEPP.b.ndim > 1 else 0
            )
            last_chunk_len = (
                self.headerChunk.get_last_chunk_length()
                if self.headerChunk is not None and hasattr(self.headerChunk, "get_last_chunk_length")
                else chunk_size
            )
            missing_chunks = 0
            for num in range(1, self.number_of_chunks):
                if not self._should_write_decoded_packet(num):
                    continue
                if num in self.solved_symbols:
                    decoded = self.solved_symbols[num]
                    output_bytes, stop_writing = self._packet_output_bytes(
                        decoded, num, null_is_terminator
                    )
                elif self.GEPP is not None and hasattr(self.GEPP, "result_mapping") and num < len(self.GEPP.result_mapping):
                    row_idx = int(self.GEPP.result_mapping[num][0]) if getattr(self.GEPP.result_mapping[num], "ndim", 0) > 0 else int(self.GEPP.result_mapping[num])
                    if 0 <= row_idx < len(self.GEPP.b):
                        chunk_data = bytes(self.GEPP.b[row_idx])
                        if num == self.number_of_chunks - 1 and last_chunk_len and last_chunk_len > 0:
                            chunk_data = chunk_data[:last_chunk_len]
                        output_bytes = chunk_data
                        stop_writing = False
                    else:
                        missing_chunks += 1
                        output_bytes = b"\x00" * (chunk_size or 1)
                        stop_writing = False
                else:
                    missing_chunks += 1
                    output_bytes = b"\x00" * (chunk_size or 1)
                    stop_writing = False

                output_concat += output_bytes
                f.write(output_bytes)
                if stop_writing:
                    break
        if missing_chunks:
            logger.warning(
                "Decode was NOT fully finished: %d/%d data chunks could not be recovered "
                "and were written as 0x00 padding.",
                missing_chunks,
                self.number_of_chunks - 1,
            )
        logger.info("Saved file as '%s'", file_name)
        self._validate_checksum(file_name, ignore_crc)
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))
        if return_file_name:
            return file_name
        return output_concat

    def _validate_checksum(self, file_name: str, ignore_crc: bool) -> None:
        """Validate the decoded file against the checksum stored in the header."""
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