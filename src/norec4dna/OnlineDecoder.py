#!/usr/bin/python
# -*- coding: latin-1 -*-
import argparse
import logging
import os
import struct
from io import BytesIO
from math import ceil
from typing import Any, Callable, Dict, Optional, Set, Union

import numpy as np
from norec4dna.Decoder import Decoder
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.ErrorCorrection import crc32, get_error_correction_decode, nocode
from norec4dna.GEPP import GEPP, GEPP_intern
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.helper import calc_crc, calc_file_crc, logical_xor, xor_mask
from norec4dna.helper.quaternary2Bin import quat_file_to_bin, tranlate_quat_to_byte
from norec4dna.OnlineAuxPacket import OnlineAuxPacket
from norec4dna.OnlinePacket import OnlinePacket
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


class OnlineDecoder(Decoder):
    def __init__(
        self,
        file: Optional[str] = None,
        error_correction: Callable[[bytes], bytes] = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: Optional[int] = None,
        read_all: bool = True,
        checksum_len_str: Optional[str] = None,
        config_map: Optional[Any] = None,
    ):
        super().__init__(file)
        self.checksum_len_str: str = checksum_len_str if checksum_len_str is not None else ""
        if not use_headerchunk and (
            self.checksum_len_str != "" and self.checksum_len_str is not None
        ):
            raise Exception("Header-checksums are only supported with headerchunks.")
        self.debug: bool = False
        self.isPseudo: bool = False
        self.file: Optional[str] = file
        self.decodedPackets: Set[OnlinePacket] = set()
        self.degreeToPacket: Dict[int, Set[OnlinePacket]] = {}
        self.f: Optional[Any] = None
        self.isFolder: bool = False
        if file is not None:
            self.isFolder = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        self.rng: Any = np.random
        self.number_of_chunks: int = 1000000
        if static_number_of_chunks is not None:
            self.number_of_chunks = static_number_of_chunks
        self.headerChunk: Optional[HeaderChunk] = None
        self.auxBlockNumbers: Dict[int, Set[int]] = dict()
        self.auxBlocks: Dict[int, OnlineAuxPacket] = dict()
        self.GEPP: Optional[GEPP_intern] = None
        self.dist: Optional[OnlineDistribution] = None
        self.read_all_before_decode: bool = read_all
        self.numberOfDecodedAuxBlocks: int = 0
        self.do_count: bool = True
        self.counter: Dict[int, int] = dict()
        self.error_correction: Callable[[bytes], bytes] = error_correction
        self.use_headerchunk: bool = use_headerchunk
        self.static_number_of_chunks: Optional[int] = static_number_of_chunks
        self.EOF: bool = False
        self.quality: int = 0
        self.epsilon: float = 0.0
        self.config_map: Optional[Any] = config_map

    def decodeFolder(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        quality_len_format: str = "I",
        epsilon_len_format: str = "f",
        check_block_number_len_format: str = "I",
    ):
        assert self.file is not None, "file must be set before calling decodeFolder!"
        decoded: bool = False
        self.EOF: bool = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        for dir_file in os.listdir(self.file):
            if dir_file.endswith(".ONLINE") or dir_file.endswith("DNA"):
                self.EOF = False
                if dir_file.endswith("DNA"):
                    self.f = quat_file_to_bin(self.file + "/" + dir_file)
                else:
                    self.f = open(self.file + "/" + dir_file, "rb")
                new_pack = self.getNextValidPacket(
                    True,
                    packet_len_format=packet_len_format,
                    crc_len_format=crc_len_format,
                    number_of_chunks_len_format=number_of_chunks_len_format,
                    quality_len_format=quality_len_format,
                    epsilon_len_format=epsilon_len_format,
                    check_block_number_len_format=check_block_number_len_format,
                )
                if new_pack is not None:
                    decoded = self.input_new_packet(new_pack)
                if decoded:
                    break
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets : %s", self.corrupt)
        if self.GEPP is not None and self.GEPP.isPotentionallySolvable():
            decoded = self.GEPP.solve()
        if hasattr(self, "f"):
            self.f.close()
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve file from chunks. Too many errors?")
            return -1

    def decodeFile(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        quality_len_format: str = "I",
        epsilon_len_format: str = "f",
        check_block_number_len_format: str = "I",
    ):
        assert self.file is not None, "file must be set before calling decodeFile!"
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        decoded: bool = False
        self.EOF: bool = False
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
                raw_packet_list.append((error_prob, seed, dna_str))
                new_pack = self.parse_raw_packet(
                    BytesIO(tranlate_quat_to_byte(dna_str)).read(),
                    crc_len_format=crc_len_format,
                    number_of_chunks_len_format=number_of_chunks_len_format,
                    epsilon_len_format=epsilon_len_format,
                    quality_len_format=quality_len_format,
                    check_block_number_len_format=check_block_number_len_format,
                )
                if new_pack is not None:
                    decoded = self.input_new_packet(new_pack)
                if self.progress_bar is not None:
                    self.progress_bar.update(self.correct, Corrupt=self.corrupt)
            else:
                while not (decoded or self.EOF):
                    new_pack = self.getNextValidPacket(
                        False,
                        packet_len_format=packet_len_format,
                        crc_len_format=crc_len_format,
                        number_of_chunks_len_format=number_of_chunks_len_format,
                        quality_len_format=quality_len_format,
                        epsilon_len_format=epsilon_len_format,
                        check_block_number_len_format=check_block_number_len_format,
                    )
                    if new_pack is None:
                        break
                    decoded = self.input_new_packet(new_pack)
                    ##
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets : %s", self.corrupt)
        if self.GEPP.isPotentionallySolvable():
            return self.GEPP.solve()
        self.f.close()
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve file from chunks. Too many errors?")
            return -1

    def decodeZip(self, *args: Any, **kwargs: Any) -> Optional[Any]:
        raise NotImplementedError("Not implemented for OnlineDecoder!")

    def createAuxBlocks(self) -> None:
        assert (
            self.number_of_chunks is not None
        ), "createAuxBlocks can only be called AFTER first Packet"
        self.rng.seed(self.number_of_chunks)
        if self.debug:
            logger.debug(
                "We should have %s Aux-Blocks and %s normal Chunks (+ 1 HeaderChunk)",
                self.getNumberOfAuxBlocks(),
                self.number_of_chunks,
            )
        for i in range(0, self.getNumberOfAuxBlocks()):
            self.auxBlockNumbers[i] = set()
        for chunk_no in range(
            0, self.number_of_chunks
        ):  # + (1 if self.use_headerchunk else 0)):  # + 1 for HeaderChunk
            # Insert this Chunk into quality different Aux-Packets
            for i in range(0, self.quality):
                # uniform choose a number of aux blocks
                aux_no = self.rng.randint(0, self.getNumberOfAuxBlocks())
                self.auxBlockNumbers[aux_no].add(chunk_no)

        # XOR all Chunks into the corresponding AUX-Block
        for aux_number in self.auxBlockNumbers.keys():
            self.auxBlocks[aux_number] = OnlineAuxPacket(
                b"",
                self.auxBlockNumbers[aux_number],
                aux_number=aux_number,
                total_number_of_chunks=self.number_of_chunks,
            )  # , numberOfAuxPackets=self.getNumberOfAuxBlocks()) # We will add the Data once we have it.

    def getAuxPacketListFromPacket(self, packet: OnlinePacket) -> NDArray[np.bool_]:
        """Return a 2D boolean array (rows = aux-blocks included by packet).

        Each row is the boolean mask returned by
        `self.auxBlocks[i].get_bool_array_used_packets()` for the aux blocks
        that are marked as used in the given packet. If no aux blocks are
        referenced by the packet an empty array with shape (0, number_of_chunks)
        is returned.
        """
        aux_used_packets = packet.getBoolArrayAuxPackets()
        rows = []
        # collect rows for aux blocks that are used
        for i, aux in enumerate(aux_used_packets):
            if aux:
                arr = self.auxBlocks[i].get_bool_array_used_packets()
                rows.append(np.asarray(arr, dtype=bool))

        if not rows:
            # no aux rows -> return empty 2D array with appropriate width
            return np.zeros((0, self.number_of_chunks), dtype=np.bool_)

        # stack into a 2D numpy array where each row corresponds to one aux block
        return np.vstack(rows).astype(np.bool_)

    def removeAndXorAuxPackets(self, packet: OnlinePacket) -> NDArray[np.bool_]:
        aux_mapping = self.getAuxPacketListFromPacket(packet)
        packet_row = np.asarray(packet.get_bool_array_used_packets(), dtype=bool)
        # combine aux rows with the packet's own row as the last row
        if aux_mapping.size == 0:
            combined = packet_row[np.newaxis, :]
        else:
            combined = np.vstack((aux_mapping, packet_row))

        return logical_xor(combined)

    def input_new_packet(self, packet: OnlinePacket) -> bool:
        if self.isPseudo and self.auxBlocks == {}:
            self.number_of_chunks = packet.get_total_number_of_chunks()
            self.quality = packet.getQuality()
            self.epsilon = round(packet.getEpsilon(), 6)
            self.dist = OnlineDistribution(self.epsilon)
            self.createAuxBlocks()
        removed: NDArray[np.bool_] = self.removeAndXorAuxPackets(packet)
        if self.do_count:
            for i in range(len(removed)):
                if i in self.counter.keys():
                    if bool(removed[i]):
                        self.counter[i] += 1
                else:
                    self.counter[i] = 1
        if self.GEPP is None:
            self.GEPP = GEPP(
                removed[np.newaxis, :],
                np.frombuffer(packet.get_data(), dtype="uint8"),
            )
        else:
            self.GEPP.addRow(
                self.removeAndXorAuxPackets(packet),
                np.frombuffer(packet.get_data(), dtype="uint8"),
            )
        if (
            self.isPseudo
            and not self.read_all_before_decode
            and (self.GEPP.isPotentionallySolvable() and self.GEPP.n % 25 == 0)
        ):
            if self.debug:
                logger.debug("current size: %s", self.GEPP.n)
            return self.GEPP.solve(partial=False)
        return False

    def solve(self, partial: bool = False) -> bool:
        assert self.GEPP is not None, "GEPP must not be None!"
        return self.GEPP.solve(partial)

    def getSolvedCount(self) -> int:
        return self.GEPP.getSolvedCount()

    def is_decoded(self) -> bool:
        return (
            self.GEPP is not None and self.GEPP.isPotentionallySolvable() and self.GEPP.isSolved()
        )

    def getNextValidPacket(
        self,
        from_multiple_files: bool = False,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        quality_len_format: str = "I",
        epsilon_len_format: str = "f",
        check_block_number_len_format: str = "I",
    ) -> Optional[OnlinePacket]:
        if not from_multiple_files:
            packet_len = self.f.read(struct.calcsize("<" + packet_len_format))
            packet_len = struct.unpack("<" + packet_len_format, packet_len)[0]
            packet = self.f.read(int(packet_len))
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF: bool = True
            self.f.close()
            return None
        res = self.parse_raw_packet(
            packet,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            quality_len_format=quality_len_format,
            epsilon_len_format=epsilon_len_format,
            check_block_number_len_format=check_block_number_len_format,
        )
        if res is None:
            res = self.getNextValidPacket(
                from_multiple_files,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                quality_len_format=quality_len_format,
                epsilon_len_format=epsilon_len_format,
                check_block_number_len_format=check_block_number_len_format,
            )
        return res

    def parse_raw_packet(
        self,
        packet: bytes,
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        quality_len_format: str = "I",
        epsilon_len_format: str = "f",
        check_block_number_len_format: str = "I",
    ) -> Optional[OnlinePacket]:
        crc_len = -struct.calcsize("<" + crc_len_format)
        if self.error_correction.__code__.co_name == crc32.__code__.co_name:
            payload = packet[:crc_len]
            crc = struct.unpack("<" + crc_len_format, packet[crc_len:])[0]
            calced_crc = calc_crc(payload)
            if crc != calced_crc:  # If the Packet is corrupt, try next one
                logger.warning("CRC-Error - %s != %s", hex(crc), hex(calced_crc))
                self.corrupt += 1
                return None
        else:
            crc_len = None
            try:
                packet = self.error_correction(packet)
            except:
                return None  # if RS or other error correction cannot reconstruct this packet
        struct_str = (
            "<"
            + number_of_chunks_len_format
            + quality_len_format
            + epsilon_len_format
            + check_block_number_len_format
        )
        struct_len = struct.calcsize(struct_str)
        data = packet[struct_len:crc_len]
        len_data = struct.unpack(struct_str, packet[0:struct_len])
        if self.static_number_of_chunks is None:
            number_of_chunks, quality, epsilon_val, check_block_number = len_data
            self.number_of_chunks = xor_mask(number_of_chunks, number_of_chunks_len_format)
            self.epsilon = round(epsilon_val, 6)
        else:
            quality, epsilon_val, check_block_number = len_data
            self.epsilon = round(epsilon_val, 6)
        self.quality = xor_mask(quality, quality_len_format)
        if self.dist is None:
            self.dist = OnlineDistribution(self.epsilon)
        if self.correct == 0:
            self.createAuxBlocks()
            # Create MockUp AuxBlocks with the given Pseudo-Random Number -> we will know which Packets are Encoded in which AuxBlock
        self.correct += 1
        res = OnlinePacket(
            data,
            self.number_of_chunks,
            self.quality,
            self.epsilon,
            check_block_number,
            dist=self.dist,
            read_only=True,
            error_correction=self.error_correction,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            quality_len_format=quality_len_format,
            epsilon_len_format=epsilon_len_format,
            check_block_number_len_format=check_block_number_len_format,
            save_number_of_chunks_in_packet=self.static_number_of_chunks is None,
        )
        return res

    def decodeHeader(self, last_chunk_len_format: str = "I") -> None:
        if self.headerChunk is not None:
            return  # Header already set
        for decoded in self.degreeToPacket[1]:
            if decoded.get_used_packets().issubset({0}):
                self.headerChunk = HeaderChunk(
                    decoded,
                    last_chunk_len_format=last_chunk_len_format,
                    checksum_len_format=self.checksum_len_str,
                )

    # def saveDecodedFile(self, last_chunk_len_format: str = "I", null_is_terminator: bool = False,
    #                    print_to_output: bool = False, return_file_name:bool = False) -> None:
    def saveDecodedFile(
        self,
        last_chunk_len_format: str = "I",
        null_is_terminator: bool = False,
        print_to_output: bool = True,
        return_file_name: bool = False,
        partial_decoding: bool = True,
    ) -> Union[bytes, str]:
        assert self.is_decoded() or partial_decoding, "Can not save File: Unable to reconstruct."
        if partial_decoding:
            self.solve(partial=True)
        if self.use_headerchunk:
            self.headerChunk = HeaderChunk(
                OnlinePacket(
                    self.GEPP.b[0],
                    self.number_of_chunks,
                    self.quality,
                    self.epsilon,
                    0,
                    {0},
                    self.dist,
                    read_only=True,
                ),
                last_chunk_len_format=last_chunk_len_format,
                checksum_len_format=self.checksum_len_str,
            )
        file_name = "DEC_" + os.path.basename(self.file) if self.file is not None else "ONLINE.BIN"
        output_concat = b""
        if self.headerChunk is not None:
            file_name = self.headerChunk.get_file_name().decode("utf-8")
        file_name = file_name.split("\x00")[0]
        with open(file_name, "wb") as f:
            # for decoded in sorted(self.degreeToPacket[1]):
            for x in self.GEPP.result_mapping:
                if 0 != x or not self.use_headerchunk:
                    if x == self.number_of_chunks - 1 and self.use_headerchunk:
                        output = self.GEPP.b[x][0][0 : self.headerChunk.get_last_chunk_length()]
                        output_concat += output.tobytes()
                        f.write(output)
                    else:
                        if null_is_terminator:
                            splitter = self.GEPP.b[x].tostring().decode().split("\x00")
                            output = splitter[0].encode()
                            output_concat += output
                            f.write(output)
                            if len(splitter) > 1:
                                break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                        else:
                            output = self.GEPP.b[x]
                            try:
                                output_concat += output.tobytes()
                            except TypeError as te:
                                raise te
                            f.write(output)
        logger.info("Saved file as '%s'", file_name)
        if self.checksum_len_str is not None and self.checksum_len_str != "":
            decoded_crc = calc_file_crc(file_name, self.checksum_len_str)
            if self.headerChunk.checksum != decoded_crc:
                logger.warning("Decoded CRC: %s", decoded_crc)
                logger.warning("Header CRC: %s", self.headerChunk.checksum)
                raise ValueError(
                    "Checksum of decoded file does not match checksum in header chunk!"
                )
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))
        if return_file_name:
            return file_name
        return output_concat

    def getNumberOfAuxBlocks(self) -> int:
        return ceil(0.55 * self.quality * self.epsilon * self.number_of_chunks)


def main(
    file: str,
    number_of_chunks: int,
    error_correction: Callable[[bytes], bytes] = nocode,
    insertheader: bool = False,
    _header_crc_str: Optional[str] = None,
):
    decoder = OnlineDecoder(
        file,
        error_correction=error_correction,
        use_headerchunk=insertheader,
        static_number_of_chunks=number_of_chunks,
        checksum_len_str=_header_crc_str,
    )
    decoder.decode(
        quality_len_format="B", check_block_number_len_format="H", number_of_chunks_len_format="H"
    )
    decoder.saveDecodedFile(null_is_terminator=False)


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
                                    nocode, crc, reedsolomon (default=nocode)",
    )
    parser.add_argument(
        "--repair_symbols",
        metavar="repair_symbols",
        type=int,
        required=False,
        default=2,
        help="number of repair symbols for ReedSolomon (default=2)",
    )
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=True, type=int)
    parser.add_argument(
        "--header_crc_str", metavar="header_crc_str", required=False, type=str, default=""
    )
    parser.add_argument(
        "--insert_header", metavar="insert_header", required=False, type=bool, default=False
    )
    args = parser.parse_args()
    _file = args.filename
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _number_of_chunks = args.number_of_chunks
    _header_crc_str = args.header_crc_str
    _error_correction = get_error_correction_decode(args.error_correction, _repair_symbols)
    logger.info("File / Folder to decode: %s", _file)
    main(
        _file,
        _number_of_chunks,
        _error_correction,
        insertheader=_insert_header,
        _header_crc_str=_header_crc_str,
    )
    logger.info("Decoding finished.")
