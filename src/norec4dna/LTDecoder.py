#!/usr/bin/python
# -*- coding: latin-1 -*-
from __future__ import annotations

import argparse
import logging
import os
import struct
from io import BytesIO
from typing import Any, BinaryIO, Callable, Dict, Optional, Set, Tuple, Union

import numpy as np

from .Decoder import Decoder
from .distributions.Distribution import Distribution
from .distributions.ErlichZielinskiRobustSolitonDisribution import (
    ErlichZielinskiRobustSolitonDistribution,
)
from .ErrorCorrection import crc32, nocode, reed_solomon_decode
from .GEPP import GEPP, GEPP_intern
from .HeaderChunk import HeaderChunk
from .helper import calc_crc, calc_file_crc, xor_mask
from .helper.quaternary2Bin import quat_file_to_bin, tranlate_quat_to_byte
from .Packet import Packet

logger = logging.getLogger(__name__)

BoolArray = np.ndarray[Any, np.dtype[np.bool_]]


class LTDecoder(Decoder):
    def __init__(
        self,
        file: Optional[str] = None,
        error_correction: Callable[[bytes], bytes] = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: Optional[int] = None,
        implicit_mode: bool = True,
        dist: Optional[Distribution] = None,
        checksum_len_str: Optional[str] = None,
        config_map: Any = None,
    ):
        super().__init__(file)
        if checksum_len_str is None:
            self.checksum_len_str = ""
        if not use_headerchunk and (checksum_len_str != "" and checksum_len_str is not None):
            raise Exception("Header-checksums are only supported with headerchunks.")
        self.checksum_len_str = checksum_len_str
        self.use_headerchunk: bool = use_headerchunk
        self.isPseudo: bool = False
        self.file: Optional[str] = file
        self.degreeToPacket: Dict[int, Set[Packet]] = {}
        if self.file is not None:
            self.isFolder: bool = os.path.isdir(self.file)
            if not self.isFolder:
                self.f: Optional[BinaryIO] = open(self.file, "rb")
        self.correct: int = 0
        self.corrupt: int = 0
        self.number_of_chunks: int = 1000000
        self.headerChunk: Optional[HeaderChunk] = None
        self.GEPP: Optional[GEPP_intern] = None
        self.pseudoCount: int = 0
        self.read_all_before_decode: bool = True
        self.count: bool = True
        self.counter: Dict[int, int] = {}
        self.error_correction: Callable[[bytes], bytes] = error_correction
        self.static_number_of_chunks: Optional[int] = static_number_of_chunks
        if static_number_of_chunks is not None:
            self.number_of_chunks = static_number_of_chunks
        self.implicit_mode: bool = implicit_mode
        self.dist: Optional[Distribution] = dist
        self.EOF: bool = False
        self.config_map: Any = config_map

    def _update_number_of_chunks_format(self, number_of_chunks_len_format: str) -> str:
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            return ""
        return number_of_chunks_len_format

    def _open_lt_fasta(self) -> None:
        assert self.file is not None, "file must not be None!"
        if hasattr(self, "f") and self.f is not None:
            self.f.close()
        self.f = open(self.file, "rb")

    def _read_lt_fasta_entry(self) -> Optional[Tuple[bytes, bytes, bytes]]:
        assert self.f is not None, "fasta file handle must be open"
        line = self.f.readline()
        if not line:
            self.EOF = True
            return None
        try:
            error_prob, seed = line[1:].replace(b"\n", b"").split(b"_")
        except ValueError:
            error_prob, seed = b"0", b"0"
        line = self.f.readline()
        if not line:
            self.EOF = True
            return None
        return error_prob, seed, line.replace(b"\n", b"")

    def _decode_lt_fasta_file(
        self,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        degree_len_format: str,
        seed_len_format: str,
    ) -> bool:
        self._open_lt_fasta()
        decoded = False
        while not (decoded or self.EOF):
            entry = self._read_lt_fasta_entry()
            if entry is None:
                break
            _, seed, dna_str = entry
            new_pack = self.parse_raw_packet(
                BytesIO(tranlate_quat_to_byte(str(dna_str))).read(),
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=number_of_chunks_len_format,
                degree_len_format=degree_len_format,
                seed_len_format=seed_len_format,
            )
            if isinstance(new_pack, Packet):
                decoded = self.input_new_packet(new_pack)
            else:
                logger.warning(
                    "Could not add a packet to the decoder: Seed: %s - %s",
                    seed,
                    new_pack,
                )
            if self.progress_bar is not None:
                self.progress_bar.update(self.correct, Corrupt=self.corrupt)
        return decoded

    def _decode_lt_binary_file(
        self,
        packet_len_format: str,
        crc_len_format: str,
        number_of_chunks_len_format: str,
        degree_len_format: str,
        seed_len_format: str,
        last_chunk_len_format: str,
    ) -> bool:
        decoded = False
        while not (decoded or self.EOF):
            new_pack = self.getNextValidPacket(
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
        return decoded

    def _build_lt_header_chunk(self, last_chunk_len_format: str) -> None:
        if not self.use_headerchunk or self.GEPP is None:
            return
        self.headerChunk = HeaderChunk(
            Packet(self.GEPP.b[0], {0}, self.number_of_chunks, read_only=True),
            last_chunk_len_format=last_chunk_len_format,
            checksum_len_format=self.checksum_len_str,
        )

    def _default_output_file_name(self) -> str:
        return "DEC_" + os.path.basename(self.file) if self.file is not None else "LT.BIN"

    def _resolved_output_file_name(self) -> str:
        file_name = self._default_output_file_name()
        if self.headerChunk is None:
            return file_name.split("\x00")[0]
        header_file_name = self.headerChunk.get_file_name()
        resolved = (
            header_file_name.decode("utf-8")
            if isinstance(header_file_name, bytes)
            else header_file_name
        )
        return resolved.split("\x00")[0]

    def _write_lt_chunk(self, x: int, null_is_terminator: bool) -> Tuple[bytes, bool, bool]:
        if self.GEPP is None:
            raise RuntimeError("GEPP not initialized")
        if x < 0:
            return b"\x00" * len(self.GEPP.b[x][0]), False, True
        if self.number_of_chunks - 1 == x and self.use_headerchunk and self.headerChunk is not None:
            output = self.GEPP.b[x][0][0 : self.headerChunk.get_last_chunk_length()]
            return output.tobytes(), False, False
        if null_is_terminator:
            splitter = self.GEPP.b[x].tobytes().decode().split("\x00")
            return splitter[0].encode(), len(splitter) > 1, False
        return self.GEPP.b[x].tobytes(), False, False

    def _validate_decoded_checksum(self, file_name: str) -> None:
        if self.checksum_len_str is None or self.checksum_len_str == "":
            return
        decoded_crc = calc_file_crc(file_name, self.checksum_len_str)
        if self.headerChunk is not None and self.headerChunk.checksum != decoded_crc:
            logger.warning("Decoded CRC: %s", decoded_crc)
            logger.warning("Header CRC: %s", self.headerChunk.checksum)
            raise ValueError("Checksum of decoded file does not match checksum in header chunk!")

    def _convert_xor_value(self, value: Any) -> int:
        return int(xor_mask(value))

    def _decode_lt_crc(
        self, packet: bytes, crc_len_format: str
    ) -> Tuple[Optional[bytes], Optional[int]]:
        crc_len = -struct.calcsize("<" + crc_len_format)
        if self.error_correction.__name__ == crc32.__name__:
            payload: bytes = packet[:crc_len]
            crc: int = struct.unpack("<" + crc_len_format, packet[crc_len:])[0]
            calced_crc: int = calc_crc(payload)
            if crc == calced_crc:
                return packet, crc_len
            logger.warning("CRC-Error - %s != %s", hex(crc), hex(calced_crc))
            self.corrupt += 1
            return None, None
        try:
            return self.error_correction(packet), None
        except Exception:
            self.corrupt += 1
            return None, None

    def _decode_lt_header(
        self,
        len_data: Tuple[Any, ...],
        number_of_chunks_len_format: str,
        degree_len_format: str,
        seed_len_format: str,
    ) -> Tuple[int, int]:
        degree: Optional[int] = None
        if self.static_number_of_chunks is None:
            if self.implicit_mode:
                number_of_chunks, seed = len_data
            else:
                number_of_chunks, degree, seed = len_data
            self.number_of_chunks = int(xor_mask(number_of_chunks, number_of_chunks_len_format))
        else:
            if self.implicit_mode:
                seed = len_data[0]
            else:
                degree, seed = len_data
        seed = int(xor_mask(seed, seed_len_format))
        if degree is None:
            if self.dist is not None:
                self.dist.set_seed(seed)
                degree = self.dist.getNumber()
            else:
                degree = 1
        else:
            degree = int(xor_mask(degree, degree_len_format))
        return int(seed), int(degree)

    def decodeFolder(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        degree_len_format: str = "I",
        seed_len_format: str = "I",
        last_chunk_len_format: str = "I",
    ) -> Optional[int]:
        decoded: bool = False
        self.EOF: bool = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        assert self.file is not None, ""
        for filename in os.listdir(self.file):
            if filename.endswith(".LT") or filename.endswith("DNA"):
                self.EOF = False
                if filename.endswith("DNA"):
                    self.f = quat_file_to_bin(self.file + "/" + filename)
                else:
                    self.f = open(self.file + "/" + filename, "rb")
                new_pack = self.getNextValidPacket(
                    True,
                    packet_len_format=packet_len_format,
                    crc_len_format=crc_len_format,
                    number_of_chunks_len_format=number_of_chunks_len_format,
                    degree_len_format=degree_len_format,
                    seed_len_format=seed_len_format,
                    last_chunk_len_format=last_chunk_len_format,
                )
                if new_pack is not None:
                    decoded = self.input_new_packet(new_pack)
                if decoded:
                    break
        self.EOF = True
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets : %s", self.corrupt)
        if self.GEPP is not None and self.GEPP.isPotentionallySolvable():
            decoded = self.GEPP.solve()
        else:
            decoded = False
        if hasattr(self, "f") and self.f is not None:
            self.f.close()
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve file from chunks. Too many errors?")
            return -1
        return None

    def decodeFile(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        degree_len_format: str = "I",
        seed_len_format: str = "I",
        last_chunk_len_format: str = "I",
    ) -> Optional[int]:
        assert self.file is not None, "file must not be None!"
        self.EOF: bool = False
        decoded: bool = False
        number_of_chunks_len_format = self._update_number_of_chunks_format(
            number_of_chunks_len_format
        )
        if self.file.lower().endswith("fasta"):
            decoded = self._decode_lt_fasta_file(
                crc_len_format,
                number_of_chunks_len_format,
                degree_len_format,
                seed_len_format,
            )
        else:
            decoded = self._decode_lt_binary_file(
                packet_len_format,
                crc_len_format,
                number_of_chunks_len_format,
                degree_len_format,
                seed_len_format,
                last_chunk_len_format,
            )
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets : %s", self.corrupt)
        if hasattr(self, "f") and self.f is not None:
            self.f.close()
        if self.GEPP is not None and self.GEPP.isPotentionallySolvable():
            return self.GEPP.solve()
        if not decoded and self.EOF:
            logger.warning("Unable to retrieve file from chunks. Too many errors?")
            return -1
        return None

    def decodeZip(self, *args: Any, **kwargs: Any) -> Optional[Any]:
        raise NotImplementedError("Not implemented in current version!")

    def input_new_packet(self, packet: Packet) -> bool:
        self.pseudoCount += 1
        packets: BoolArray = packet.get_bool_array_used_packets()
        if self.count:
            for i in range(len(packets)):
                if i in self.counter.keys():
                    if packets[i]:
                        self.counter[i] += 1
                else:
                    self.counter[i] = 1
        if self.GEPP is None:
            self.GEPP = GEPP(
                np.array([packet.get_bool_array_used_packets()], dtype=bool),
                np.frombuffer(packet.get_data(), dtype="uint8"),
            )
        else:
            self.GEPP.addRow(
                packet.get_bool_array_used_packets(),
                np.frombuffer(packet.get_data(), dtype="uint8"),
            )
        if (
            self.isPseudo
            and not self.read_all_before_decode
            and self.GEPP.isPotentionallySolvable()
        ):
            return self.GEPP.solve(partial=False)
        return False

    def solve(self, partial: bool = False) -> bool:
        if self.GEPP is None:
            return False
        return self.GEPP.solve(partial)

    def getSolvedCount(self) -> int:
        if self.GEPP is None:
            return 0
        return self.GEPP.getSolvedCount()

    def choose_packet_numbers(self, degree: int, seed: int = 0) -> Set[int]:
        assert degree <= self.number_of_chunks
        res: Set[int] = set()
        rng = np.random
        rng.seed(seed)
        for _ in range(0, degree):
            tmp: int = rng.choice(range(0, self.number_of_chunks))
            while tmp in res:
                tmp = rng.choice(range(0, self.number_of_chunks))
            res.add(tmp)
        return res

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
        degree_len_format: str = "I",
        seed_len_format: str = "I",
        last_chunk_len_format: str = "I",
    ) -> Optional[Packet]:
        if not hasattr(self, "f") or self.f is None:
            return None
        if not from_multiple_files:
            packet_len_bytes = self.f.read(struct.calcsize("<" + packet_len_format))
            if not packet_len_bytes:
                self.EOF = True
                self.f.close()
                return None
            packet_len = struct.unpack("<" + packet_len_format, packet_len_bytes)[0]
            packet = self.f.read(int(packet_len))
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF = True
            self.f.close()
            return None
        res = self.parse_raw_packet(
            packet,
            number_of_chunks_len_format=number_of_chunks_len_format,
            degree_len_format=degree_len_format,
            seed_len_format=seed_len_format,
        )
        if res is None:
            res = self.getNextValidPacket(
                from_multiple_files=from_multiple_files,
                number_of_chunks_len_format=number_of_chunks_len_format,
                degree_len_format=degree_len_format,
                seed_len_format=seed_len_format,
            )
        return res

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
        dirty = False
        self._build_lt_header_chunk(last_chunk_len_format)
        file_name = self._resolved_output_file_name()
        output_concat: bytes = b""
        if self.GEPP is None:
            raise RuntimeError("GEPP not initialized")
        with open(file_name, "wb") as f:
            for x in self.GEPP.result_mapping:
                output, stop_writing, marked_dirty = self._write_lt_chunk(x, null_is_terminator)
                if x >= 0 and x == 0 and self.use_headerchunk:
                    continue
                output_concat += output
                f.write(output)
                dirty = dirty or marked_dirty
                if stop_writing:
                    break
        logger.info("Saved file as '%s'", file_name)
        self._validate_decoded_checksum(file_name)
        if dirty:
            logger.warning(
                "Some parts could not be restored, file WILL contain sections with \\x00 !"
            )
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))
        if return_file_name:
            return file_name
        return output_concat

    def removeAndXorAuxPackets(self, packet: Packet) -> BoolArray:
        """
        For LT this is an identity function (makes writing code for all three Coders easier)
        :param packet:
        :return:
        """
        packet_row = np.asarray(packet.get_bool_array_used_packets(), dtype=bool)
        return packet_row

    def parse_raw_packet(
        self,
        packet: bytes,
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        degree_len_format: str = "I",
        seed_len_format: str = "I",
    ) -> Optional[Packet]:
        decoded_packet, crc_len = self._decode_lt_crc(packet, crc_len_format)
        if decoded_packet is None:
            return None
        packet = decoded_packet
        if self.implicit_mode:
            degree_len_format = ""
        struct_str: str = "<" + number_of_chunks_len_format + degree_len_format + seed_len_format
        struct_len: int = struct.calcsize(struct_str)
        len_data = struct.unpack(struct_str, packet[0:struct_len])
        seed, degree = self._decode_lt_header(
            len_data, number_of_chunks_len_format, degree_len_format, seed_len_format
        )
        used_packets = self.choose_packet_numbers(degree, seed)
        data = packet[struct_len:crc_len] if crc_len is not None else packet[struct_len:]
        self.correct += 1

        return Packet(
            data,
            used_packets,
            self.number_of_chunks,
            read_only=True,
            error_correction=self.error_correction,
            save_number_of_chunks_in_packet=self.static_number_of_chunks is None,
        )


def main(
    file: str,
    number_of_chunks: int,
    error_correction: Callable[[bytes], bytes],
    insertheader: bool,
    _header_crc_str: Optional[str] = None,
) -> None:
    dist = ErlichZielinskiRobustSolitonDistribution(number_of_chunks, seed=2)

    decoder = LTDecoder(
        file,
        error_correction=error_correction,
        use_headerchunk=insertheader,
        static_number_of_chunks=number_of_chunks,
        implicit_mode=False,
        dist=dist,
        checksum_len_str=_header_crc_str,
    )
    decoder.decode(number_of_chunks_len_format="H", seed_len_format="I", degree_len_format="H")
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
    parser.add_argument(
        "--insert_header", metavar="insert_header", required=False, type=bool, default=False
    )
    parser.add_argument(
        "--header_crc_str", metavar="header_crc_str", required=False, type=str, default=""
    )
    parser.add_argument("--number_of_chunks", metavar="number_of_chunks", required=True, type=int)
    args = parser.parse_args()
    filename = args.filename
    e_correction_str = args.error_correction
    _repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _header_crc_str = args.header_crc_str
    _number_of_chunks = args.number_of_chunks
    e_correction_fn: Callable[..., bytes]
    if e_correction_str == "nocode":
        e_correction_fn = nocode
    elif e_correction_str == "crc":
        e_correction_fn = crc32
    elif e_correction_str == "reedsolomon":
        if _repair_symbols != 2:

            def custom_e_correction(data: bytes) -> bytes:
                return reed_solomon_decode(data, _repair_symbols)

            e_correction_fn = custom_e_correction
        else:
            e_correction_fn = reed_solomon_decode
    else:
        logger.error(
            "Selected Error Correction not supported, choose: 'nocode', 'crc' or 'reedsolomon'"
        )
        raise SystemExit(1)
    logger.info("File / Folder to decode: %s", filename)
    main(filename, _number_of_chunks, e_correction_fn, _insert_header, _header_crc_str)
    logger.info("Decoding finished.")
