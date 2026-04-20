import math
import os
import struct
import typing
from abc import ABC
from math import ceil
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import progressbar
from norec4dna import Decoder
from norec4dna.ErrorCorrection import nocode
from numpy.typing import NDArray
from PIL import Image


class Encoder(ABC):
    def __init__(
        self,
        file: str,
        number_of_chunks: int,
        distribution: typing.Any,
        insert_header: bool = True,
        pseudo_decoder: typing.Optional[Decoder] = None,
        chunk_size: int = 0,
        mode_1_bmp: bool = False,
    ):
        self.crc_len_format: str = "I"
        self.id_len_format: str = "I"
        self.number_of_chunks_len_format: str = "I"
        self.number_of_chunks: int = number_of_chunks
        self.chunk_size: int = chunk_size
        self.pseudo_decoder: typing.Optional[typing.Any] = pseudo_decoder
        self.insert_header: bool = insert_header
        self.dist: typing.Any = distribution
        self.file: str = file
        self.mode_1_bmp: bool = mode_1_bmp
        self.setOfEncodedPackets: typing.Set[int] = set()
        self.encodedPackets: typing.Set[typing.Any] = set()
        self.overhead_limit: typing.Optional[float] = None
        self.chunks: typing.List[np.ndarray] = []
        self.error_correction: typing.Callable[[bytes], bytes] = nocode
        self.progress_bar: typing.Optional[progressbar.ProgressBar] = None
        self.ruleDrop: int = 0
        self.out_file: typing.Optional[str] = None

    @staticmethod
    def create_progress_bar(max_value: int) -> progressbar.ProgressBar:
        widgets: typing.List[typing.Any] = [
            progressbar.Percentage(),
            progressbar.Bar(),
            " Encoded: ",
            progressbar.Counter(),
            ", ",
            progressbar.Variable("Dropped"),
            ", ",
            progressbar.AdaptiveETA(),
            " ",
            progressbar.Timer(),
        ]
        return progressbar.ProgressBar(
            max_value=max_value, widgets=widgets, max_error=False, redirect_stdout=False
        ).start()

    def encode_to_packets(self) -> None:
        pass  # implemented in subclasses

    def set_overhead_limit(self, n: float) -> None:
        self.overhead_limit = n

    def encode_file(self, split_to_multiple_files: bool = False) -> None:
        pass  # implemented in subclasses

    def create_new_packet(self) -> typing.Optional[typing.Any]:
        pass  # implemented in subclasses

    def set_no_chunks_from_chunk_size(self) -> None:
        file_size = self.get_file_size(self.file)
        self.number_of_chunks = ceil(1.0 * file_size / self.chunk_size) + (
            1 if self.insert_header else 0
        )

        print("number_of_chunks from chunk_size:" + str(self.number_of_chunks))

    def update_progress_bar(self) -> None:
        if self.progress_bar is not None:
            self.progress_bar.update(len(self.get_encoded_packets()), Dropped=self.ruleDrop)

    @staticmethod
    def get_number_of_chunks_for_file_with_chunk_size(
        file: str, chunk_size: int, insert_header: bool = True
    ) -> int:
        file_size = Encoder.get_file_size(file)
        return ceil(1.0 * file_size / chunk_size) + (1 if insert_header else 0)

    def create_and_add_new_packet(self) -> typing.Any:
        packet = self.create_new_packet()
        self.encodedPackets.add(packet)
        return packet

    def encode_header_info(
        self,
        checksum: typing.Optional[bytes] = None,
        checksum_len_str: typing.Optional[str] = None,
        last_chunk_len_format: str = "I",
    ) -> np.ndarray:
        # Size of last Chunk
        # Filename
        # PAD-Bytes
        if checksum_len_str is None:
            checksum_len_str = ""
        last_chunk = self.chunks[-1]
        file_name_only = os.path.basename(self.file)
        file_name_length = len(bytes(file_name_only, encoding="utf-8"))
        assert (
            file_name_length + 4 + struct.calcsize("" + last_chunk_len_format + checksum_len_str)
            < self.chunk_size
        ), "Chunks too small for HeaderInfo"
        struct_string = (
            "<"
            + last_chunk_len_format
            + checksum_len_str
            + str(file_name_length)
            + "s"
            + str(
                self.chunk_size
                - file_name_length
                - struct.calcsize("<" + last_chunk_len_format + checksum_len_str)
            )
            + "x"
        )
        if checksum_len_str == "" or checksum is None:
            return np.frombuffer(
                struct.pack(
                    struct_string, len(last_chunk), bytes(file_name_only, encoding="utf-8")
                ),
                dtype=np.uint8,
            )
        return np.frombuffer(
            struct.pack(
                struct_string, len(last_chunk), checksum, bytes(file_name_only, encoding="utf-8")
            ),
            dtype=np.uint8,
        )

    def fill_last_chunk(self, force_fill_zero: bool = False) -> None:
        last = self.chunks[-1]
        assert len(last) <= self.chunk_size, "Error, last Chunk ist bigger than ChunkSize"
        if len(last) < self.chunk_size:
            if self.insert_header and not force_fill_zero:
                # only add random bytes if xor_by_seed is False (otherwise packet will be scrambled by xor)
                filler = b"\x00" + os.urandom(self.chunk_size - len(last) - 1)
            else:
                filler = (self.chunk_size - len(last)) * b"\x00"
            struct_str = "<" + str(len(last)) + "s" + str(self.chunk_size - len(last)) + "s"
            self.chunks[-1] = np.frombuffer(
                struct.pack(struct_str, bytes(last), filler), dtype=np.uint8
            )

    def number_of_packets_encoded_already(self) -> int:
        return len(self.encodedPackets)

    def save_packets(
        self, split_to_multiple_files: bool, out_file: typing.Optional[str] = None
    ) -> None:
        pass  # implemented in subclasses

    def create_chunks(self, chunk_size: int) -> typing.List[NDArray[np.uint8]]:
        if hasattr(self, "") and self.mode_1_bmp:
            data = self.image_to_mode_1_bmp()
        else:
            with open(self.file, "rb") as f:
                data = f.read()
        res = [
            np.frombuffer(data[i : i + chunk_size], dtype=np.uint8)
            for i in range(0, len(data), chunk_size)
        ]
        if self.number_of_chunks != len(res):
            print(
                "Number of Chunks does not work for given file. New Number of Chunks = "
                + str(len(res) + (1 if self.insert_header else 0))
            )
            self.number_of_chunks = len(res)
        self.progress_bar = self.create_progress_bar(self.number_of_chunks)
        return res

    def image_to_mode_1_bmp(self) -> bytes:
        return self.create_bit_arr(Image.open(self.file))

    def create_bit_arr(self, img: Image.Image) -> bytes:
        arr = np.frombuffer(img.tobytes(), dtype=np.uint8)
        remain = len(arr) % 8
        if remain != 0:
            return self.translate_to_bytes(
                np.append(arr, np.zeros(8 - remain, dtype=np.uint8)), img
            )
        else:
            return self.translate_to_bytes(arr, img)

    def save_packets_zip(
        self,
        save_as_dna: bool = False,
        out_file: typing.Optional[str] = None,
        file_ending: str = "",
        seed_is_filename: bool = True,
    ) -> None:
        if out_file is None:
            out_file = self.file + file_ending
            self.out_file = os.path.relpath(out_file)
        if not out_file.endswith(".zip"):
            out_file = out_file + ".zip"
        i = 0
        abs_dir = os.path.split(os.path.abspath("../" + out_file))[0]
        if not os.path.exists(abs_dir):
            os.makedirs(abs_dir)

        with ZipFile(out_file, "w") as f:
            for packet in self.encodedPackets:
                if seed_is_filename:
                    i = packet.id
                if save_as_dna:
                    e_prob = (
                        (str(ceil(packet.error_prob * 100)) + "_")
                        if packet.error_prob is not None
                        else ""
                    )
                    f.writestr(f"{i}_{e_prob}{file_ending}", packet.get_dna_struct(True))
                else:
                    f.writestr(f"{i}{file_ending}", packet.get_struct(True))
                i += 1
        print(f"Saved result at: %s" % out_file)

    def save_packets_fasta(
        self,
        out_file: typing.Optional[str] = None,
        file_ending: str = "",
        seed_is_filename: bool = True,
    ) -> None:
        if out_file is None:
            out_file = self.file + file_ending
            self.out_file = os.path.relpath(out_file)
        if not out_file.endswith(".fasta"):
            out_file = out_file + ".fasta"
        i = 0
        abs_dir = Path(out_file).absolute().parent
        if not os.path.exists(abs_dir):
            os.makedirs(abs_dir)

        with open(out_file, "w") as f:
            for packet in self.encodedPackets:
                if seed_is_filename:
                    i = packet.id
                e_prob = (
                    (str(ceil(packet.error_prob * 100)) + "_")
                    if packet.error_prob is not None
                    else ""
                )
                f.write(
                    ">" + e_prob + str(i) + file_ending + "\n" + packet.get_dna_struct(True) + "\n"
                )
                i += 1
        print(f"Saved result at: %s" % out_file)

    @staticmethod
    def translate_to_bytes(bit_arr: np.ndarray, img: Image.Image) -> bytes:
        width, height = img.size
        img_byt = bytes(
            [i for i in np.packbits([bit_arr[i : i + 8] for i in range(0, len(bit_arr), 8)])]
        )
        return struct.pack(">HH", width, height) + img_byt

    @staticmethod
    def calc_max_size(number_bytes: int) -> float:
        if number_bytes == 1:
            return 256
        elif number_bytes == 2:
            return 65536
        elif number_bytes == 4:
            return 4294967296
        elif number_bytes == 8:
            return 18446744073709551616
        else:
            return math.pow(2, 8 * number_bytes)

    @staticmethod
    def get_file_size(file: str) -> int:
        return os.stat(file).st_size

    def get_encoded_packets(self) -> typing.Set[typing.Any]:
        return self.encodedPackets
