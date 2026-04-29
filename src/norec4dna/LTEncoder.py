#!/usr/bin/python
# -*- coding: latin-1 -*-
from __future__ import annotations

import argparse
import configparser
import datetime
import glob
import os
import struct
import time
import typing
from math import ceil

import numpy as np

from .Decoder import Decoder
from .distributions.Distribution import Distribution
from .distributions.ErlichZielinskiRobustSolitonDisribution import (
    ErlichZielinskiRobustSolitonDistribution,
)
from .Encoder import ChunkData, Encoder
from .ErrorCorrection import crc32, nocode, reed_solomon_encode
from .helper import listXOR, should_drop_packet
from .Packet import Packet
from .rules.DNARules import DNARules
from .rules.DNARules2 import DNARules2
from .rules.DNARules_ErlichZielinski import DNARules_ErlichZielinski
from .rules.FastDNARules import FastDNARules

UInt8Array = np.ndarray[typing.Any, np.dtype[np.uint8]]


class LTEncoder(Encoder):
    def __init__(
        self,
        file: str,
        number_of_chunks: int,
        distribution: Distribution,
        insert_header: bool = True,
        pseudo_decoder: typing.Optional[Decoder] = None,
        prioritized_packets=None,
        chunk_size: int = 0,
        error_correction: typing.Callable = nocode,
        rules: typing.Optional[
            typing.Union[DNARules, DNARules2, FastDNARules, DNARules_ErlichZielinski]
        ] = None,
        implicit_mode: bool = True,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        used_packets_len_format: str = "I",
        id_len_format: str = "I",
        last_chunk_len_format: str = "I",
        save_number_of_chunks_in_packet: bool = True,
        drop_upper_bound=1.0,
        sequential_seed=True,
        checksum_len_str=None,
    ):
        super().__init__(
            file, number_of_chunks, distribution, insert_header, pseudo_decoder, chunk_size
        )
        if prioritized_packets is None:
            prioritized_packets = []
        assert number_of_chunks == distribution.get_size()
        self.out_file: typing.Optional[str] = None
        self.dist: Distribution = distribution
        self.rng: np.random.RandomState = np.random.RandomState()
        self.insert_header: bool = insert_header
        self.chunk_size: int = chunk_size
        self.file: str = file
        self.rules: typing.Optional[
            typing.Union[DNARules, DNARules2, FastDNARules, DNARules_ErlichZielinski]
        ] = rules
        self.upper_bound: float = drop_upper_bound
        if self.chunk_size == 0:
            self.number_of_chunks: int = number_of_chunks
        else:
            self.set_no_chunks_from_chunk_size()
        self.chunks: typing.List[ChunkData] = []
        self.overhead_limit: float = 0.20
        self.encodedPackets: typing.Set[Packet] = set()
        self.setOfEncodedPackets: typing.Set[int] = set()
        self.pseudo_decoder: typing.Optional[Decoder] = pseudo_decoder
        self.prioritized_packets: typing.List = prioritized_packets
        self.error_correction: typing.Callable[[bytes], bytes] = error_correction
        self.implicit_mode: bool = implicit_mode
        # Struct-Strings:
        self.packet_len_format: str = packet_len_format
        self.crc_len_format: str = crc_len_format
        self.number_of_chunks_len_format: str = number_of_chunks_len_format
        self.used_packets_len_format: str = used_packets_len_format
        self.id_len_format: str = id_len_format
        self.last_chunk_len_format: str = last_chunk_len_format
        self.save_number_of_chunks_in_packet: bool = save_number_of_chunks_in_packet
        self.ruleDrop: int = 0
        self.next_checkblock_id = -1
        self.sequential_seed = sequential_seed
        if checksum_len_str is not None:
            print("WARNING: Checksum not implemented yet, checksum_len_str will be ignored!")
        self.progress_bar = self.create_progress_bar(
            int(self.number_of_chunks + 0.02 * self.number_of_chunks)
        )

    def encode_file(self, split_to_multiple_files: bool = False):
        self.encode_to_packets()
        self.save_packets(split_to_multiple_files)

    def prepareEncoder(self):
        self.ruleDrop = 0
        file_size: int = self.get_file_size(self.file)
        if self.insert_header:
            self.chunk_size = ceil(1.0 * file_size / (self.number_of_chunks - 1))
            self.chunks = self.create_chunks(self.chunk_size)
            self.number_of_chunks += 1
            # First Chunk is a Header - convert NDArray to bytes
            header_info: UInt8Array = self.encode_header_info()
            self.chunks.insert(0, header_info.tobytes())
        else:
            self.chunk_size = ceil(1.0 * file_size / self.number_of_chunks)
            self.chunks = self.create_chunks(self.chunk_size)
        self.dist.update_number_of_chunks(self.number_of_chunks)
        self.fill_last_chunk()

    def encode_to_packets(self) -> float:
        self.prepareEncoder()
        start: float = time.time()
        if self.pseudo_decoder is not None:
            while not self.pseudo_decoder.is_decoded():
                # This process continues until the receiver signals that the
                # message has been received and successfully decoded.
                new_pack: Packet = self.create_new_packet()
                if self.rules is not None:
                    while should_drop_packet(self.rules, new_pack, self.upper_bound):
                        self.update_progress_bar()
                        new_pack = self.create_new_packet()
                        self.ruleDrop += 1

                self.pseudo_decoder.input_new_packet(new_pack)
                self.encodedPackets.add(new_pack)
        else:
            while (
                len(self.encodedPackets)
                < (self.number_of_chunks + (self.number_of_chunks * self.overhead_limit))
                or self.number_of_packets_encoded_already() < self.number_of_chunks
            ):  # 20% Aufschlag
                # This process continues until the receiver signals that the
                # message has been received and successfully decoded.
                pack: Packet = self.create_new_packet()
                if self.rules is not None:
                    while should_drop_packet(self.rules, pack, self.upper_bound):
                        pack = self.create_new_packet()
                        self.ruleDrop += 1
                self.encodedPackets.add(pack)
        return start

    def set_overhead_limit(self, n: float):
        self.overhead_limit = n

    def encodePriotizedPackets(self, error_correction: typing.Callable = nocode):
        for num in self.prioritized_packets:
            new_pack = Packet(
                self.chunks[num],
                {num},
                self.number_of_chunks,
                read_only=False,
                implicit_mode=self.implicit_mode,
                error_correction=error_correction,
                packet_len_format=self.packet_len_format,
                crc_len_format=self.crc_len_format,
                number_of_chunks_len_format=self.number_of_chunks_len_format,
                used_packets_len_format=self.used_packets_len_format,
                id_len_format=self.id_len_format,
            )
            if self.pseudo_decoder is not None:
                self.pseudo_decoder.input_new_packet(new_pack)
            self.encodedPackets.add(new_pack)

    def encode_header_info(self) -> UInt8Array:
        # Size of last Chunk
        # Filename
        # PAD-Bytes
        last_chunk: ChunkData = self.chunks[-1]
        file_name_length: int = len(self.file)
        assert file_name_length + 4 < self.chunk_size, "Chunks too small for HeaderInfo"
        # -4 for bytes to store length of last_chunk (I)
        last_chunk_len_struct_size = struct.calcsize("<" + self.last_chunk_len_format)
        struct_str = (
            "<"
            + self.last_chunk_len_format
            + str(file_name_length)
            + "s"
            + str(self.chunk_size - file_name_length - last_chunk_len_struct_size)
            + "x"
        )
        packed_data: bytes = struct.pack(
            struct_str, len(last_chunk), bytes(self.file, encoding="utf-8")
        )
        # Convert bytes to NDArray
        return np.frombuffer(packed_data, dtype=np.uint8)

    def number_of_packets_encoded_already(self) -> int:
        return len(self.setOfEncodedPackets)

    def _packet_output_mode(self, save_as_dna: bool) -> str:
        return "w" if save_as_dna else "wb"

    def _packet_output_data(self, packet, split_to_multiple_files: bool, save_as_dna: bool):
        return (
            packet.get_dna_struct(split_to_multiple_files)
            if save_as_dna
            else packet.get_struct(split_to_multiple_files)
        )

    def _default_output_folder(self, prefix: str, clear_output: bool) -> str:
        fulldir, filename = os.path.split(os.path.realpath(self.file))
        out_file = os.path.join(fulldir, prefix + filename)
        files = glob.glob(out_file + ("*" if out_file.endswith("/") else "/*"))
        if clear_output:
            for file_name in files:
                os.remove(file_name)
        return out_file

    def _save_packet_stream(
        self, out_file: str, split_to_multiple_files: bool, save_as_dna: bool
    ) -> None:
        with open(out_file, self._packet_output_mode(save_as_dna)) as f:
            for packet in self.encodedPackets:
                f.write(self._packet_output_data(packet, split_to_multiple_files, save_as_dna))

    def _save_packet_folder(
        self,
        out_file: str,
        file_ending: str,
        split_to_multiple_files: bool,
        save_as_dna: bool,
        seed_is_filename: bool,
    ) -> None:
        packet_index = 0
        error_prefix = ""
        if not os.path.exists(out_file):
            os.makedirs(out_file)
        for packet in sorted(
            self.encodedPackets, key=lambda elem: (elem.error_prob, elem.__hash__())
        ):
            if seed_is_filename:
                packet_index = packet.id
                error_prefix = (
                    (str(ceil(packet.error_prob * 100)) + "_")
                    if packet.error_prob is not None
                    else ""
                )
            packet_path = out_file + "/" + error_prefix + str(packet_index) + file_ending
            with open(packet_path, self._packet_output_mode(save_as_dna)) as f:
                f.write(self._packet_output_data(packet, split_to_multiple_files, save_as_dna))
            packet_index += 1

    def save_packets(
        self,
        split_to_multiple_files: bool,
        out_file: typing.Optional[str] = None,
        save_as_dna: bool = False,
        clear_output: bool = True,
        seed_is_filename: bool = False,
    ) -> None:
        """
        Saves the generated packets either to multiple files or to a single one. It's possible to save the packets
        either as DNA or binary.
        :param split_to_multiple_files: True: Saves the packets in multiple files. False: Saves all packets in one file.
        :param out_file: The location of the output file
        :param save_as_dna: True: Saves the information in bases. False: Saves the information binary.
        :param clear_output: Clears the location of the output file.
        :param seed_is_filename: True: Sets the seed as filename.
        :return:
        """
        file_ending: str = ".LT" + ("_DNA" if save_as_dna else "")
        if not split_to_multiple_files:
            if out_file is None:
                out_file = self.file + file_ending
            self._save_packet_stream(out_file, split_to_multiple_files, save_as_dna)
        elif out_file is None:
            out_file = self._default_output_folder("LT_", clear_output)
            self._save_packet_folder(
                out_file, file_ending, split_to_multiple_files, save_as_dna, seed_is_filename
            )
        else:
            self._save_packet_folder(
                out_file, file_ending, split_to_multiple_files, save_as_dna, seed_is_filename
            )
        self.out_file = os.path.relpath(out_file)
        print("Config: " + self.getConfigStr(out_file))

    def get_file_size(self, file: str) -> int:
        return os.stat(file).st_size

    def generate_new_checkblock_id(self, sequential=True) -> int:
        max_num: float = Encoder.calc_max_size(struct.calcsize("<" + self.id_len_format))
        if sequential:
            self.next_checkblock_id += 1
            if self.next_checkblock_id > max_num:
                raise RuntimeError("sequential checkblock_id > max allowed number!")
            return self.next_checkblock_id
        random_generator: np.random.RandomState = np.random.RandomState()
        return int(random_generator.randint(0, int(max_num), dtype=np.uint32))

    def create_new_packet(self, seed: typing.Optional[int] = None) -> Packet:
        generated_seed: int
        if seed is None:
            generated_seed = self.generate_new_checkblock_id(self.sequential_seed)
        else:
            generated_seed = seed
        if (
            self.implicit_mode
        ):  # in implicit mode we want to be able derive the used chunks from having only the seed
            self.dist.set_seed(generated_seed)
        degree: int = self.dist.getNumber()

        packet_numbers: typing.Set[int] = self.choose_packet_numbers(degree, seed=generated_seed)
        chunks: typing.List[ChunkData] = [self.chunks[i] for i in packet_numbers]
        self.setOfEncodedPackets |= set(packet_numbers)
        return Packet(
            listXOR(chunks),
            packet_numbers,
            self.number_of_chunks,
            read_only=False,
            seed=generated_seed,
            error_correction=self.error_correction,
            implicit_mode=self.implicit_mode,
            packet_len_format=self.packet_len_format,
            crc_len_format=self.crc_len_format,
            number_of_chunks_len_format=self.number_of_chunks_len_format,
            used_packets_len_format=self.used_packets_len_format,
            id_len_format=self.id_len_format,
            save_number_of_chunks_in_packet=self.save_number_of_chunks_in_packet,
        )

    def create_and_add_new_packet(self, error_correction=nocode):
        packet = self.create_new_packet()
        self.encodedPackets.add(packet)
        return packet

    def choose_packet_numbers(self, degree: int, seed: int = 0) -> typing.Set[int]:
        assert degree <= len(self.chunks)
        res: typing.Set[int] = set()
        self.rng.seed(seed)
        for _ in range(0, degree):
            tmp = self.rng.choice(range(0, len(self.chunks)))
            while tmp in res:
                tmp = self.rng.choice(range(0, len(self.chunks)))
            res.add(tmp)
        return res

    def save_config_file(
        self,
        default_map: typing.Optional[typing.Dict[str, typing.Any]] = None,
        section_name: typing.Optional[str] = None,
    ):
        if default_map is None:
            default_map = {}
        if section_name is None:
            section_name = str(self.out_file)
        config = configparser.ConfigParser()
        config[section_name] = {
            "algorithm": "LT",
            "error_correction": self.error_correction.__code__.co_name,
            "insert_header": str(self.insert_header),
            "savenumberofchunks": str(self.save_number_of_chunks_in_packet),
            "mode_1_bmp": str(self.mode_1_bmp),
            "upper_bound": str(self.upper_bound),
            "number_of_chunks": str(self.number_of_chunks),
            "config_str": self.getConfigStr(),
            "id_len_format": self.id_len_format,
            "number_of_chunks_len_format": self.number_of_chunks_len_format,
            "packet_len_format": self.packet_len_format,
            "crc_len_format": self.crc_len_format,
            "master_seed": "0",
            "distribution": self.dist.get_config_string(),
            "rules": str(list(self.rules.active_rules)) if self.rules else "[]",
            "chunk_size": str(self.chunk_size),
            "dropped_packets": str(self.ruleDrop),
            "created_packets": str(len(self.encodedPackets)),
        }
        for key, val in default_map.items():
            config[section_name][str(key)] = str(val)
        config_file_name = "{}_{}.ini".format(
            self.file, datetime.datetime.now().ctime().replace(" ", "_").replace(":", "_")
        )
        with open(config_file_name, "w") as config_file:
            config.write(config_file)
        return config_file_name

    def getConfigStr(self, out_file=""):
        res = (
            "USE_HEADER_CHUNK: "
            + str(self.insert_header)
            + ", NUMBER_OF_CHUNKS: "
            + str(self.number_of_chunks)
            + " NUMBER_OF_CHUNKS_LEN_FORMAT: "
            + self.number_of_chunks_len_format
            + " ID_LEN_FORMAT: "
            + self.id_len_format
            + " ERROR_CORRECTION: "
            + self.error_correction.__code__.co_name
            + " CRC_LEN_FORMAT(Optional): "
            + self.crc_len_format
            + " FILE: "
            + self.file
            + " OUT_FILE: "
            + out_file
            + " Distribution: "
            + self.dist.get_config_string()
        )
        return res


def main(
    file,
    number_of_chunks: int = 0,
    chunk_size: int = 0,
    error_correction: typing.Callable = nocode,
    as_dna: bool = False,
    insert_header: bool = False,
    save_number_of_chunks=False,
):
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    if as_dna:
        rules = FastDNARules()
    else:
        rules = None
    dist = ErlichZielinskiRobustSolitonDistribution(number_of_chunks, seed=2)
    encoder = LTEncoder(
        file,
        number_of_chunks,
        dist,
        insert_header=insert_header,
        rules=rules,
        error_correction=error_correction,
        number_of_chunks_len_format="H",
        id_len_format="I",
        used_packets_len_format="H",
        save_number_of_chunks_in_packet=save_number_of_chunks,
        implicit_mode=False,
    )
    encoder.encode_to_packets()
    print("Number of Chunks=%s" % encoder.number_of_chunks)
    encoder.save_packets(split_to_multiple_files=True, save_as_dna=as_dna)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--as_dna",
        help="convert packets to dna and use dna rules",
        action="store_true",
        required=False,
    )
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument(
        "--chunk_size",
        metavar="chunk_size",
        required=False,
        type=int,
        default=0,
        help="size of chunks to split the file into",
    )
    parser.add_argument(
        "--number_of_chunks",
        metavar="number_of_chunks",
        required=False,
        type=int,
        default=0,
        help="number of chunks to split the file into,only used if no chunk_size is given",
    )
    parser.add_argument(
        "--error_correction",
        metavar="error_correction",
        required=False,
        type=str,
        default="nocode",
        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)",
    )
    parser.add_argument(
        "--repair_symbols",
        metavar="repair_symbols",
        required=False,
        type=int,
        default=2,
        help="number of repair symbols for ReedSolomon (default=2)",
    )
    parser.add_argument(
        "--insert_header", metavar="insert_header", required=False, type=bool, default=False
    )
    parser.add_argument(
        "--save_number_of_chunks",
        metavar="save_number_of_chunks",
        required=False,
        type=bool,
        default=False,
    )

    args = parser.parse_args()
    filename = args.filename
    e_correction_str = args.error_correction
    _repair_symbols = args.repair_symbols
    _chunk_size = args.chunk_size
    _number_of_chunks = args.number_of_chunks
    _insert_header = args.insert_header
    _as_dna = args.as_dna
    _save_number_of_chunks = args.save_number_of_chunks
    if _chunk_size == _number_of_chunks == 0:
        print("Please set either a chunk_size or a number_of_chunks")
        exit()
    e_correction_fn: typing.Callable[..., bytes]
    if e_correction_str == "nocode":
        e_correction_fn = nocode
    elif e_correction_str == "crc":
        e_correction_fn = crc32
    elif e_correction_str == "reedsolomon":
        if _repair_symbols != 2:

            def custom_e_correction(data: bytes) -> bytes:
                return reed_solomon_encode(data, _repair_symbols)

            e_correction_fn = custom_e_correction
        else:
            e_correction_fn = reed_solomon_encode
    else:
        print("Selected Error Correction not supported, choose: 'nocode', 'crc' or 'reedsolomon'")
        raise SystemExit(1)
    filename = args.filename
    print("File to encode: " + str(filename))
    main(
        filename,
        _number_of_chunks,
        _chunk_size,
        e_correction_fn,
        _as_dna,
        _insert_header,
        save_number_of_chunks=_save_number_of_chunks,
    )
    print("File encoded.")
