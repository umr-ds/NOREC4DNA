import struct
from numpy.typing import NDArray
import numpy
import typing

import numpy as np
from norec4dna.Packet import Packet


class HeaderChunk:
    def __init__(self, packet: Packet, last_chunk_len_format: str = "I", checksum_len_format: typing.Optional[str] = None):
        assert packet.get_used_packets().issubset({0}), "only first packet can be HeaderPacket"
        if isinstance(packet.data, numpy.ndarray):
            self.data: bytes = packet.data.tobytes()
        else:
            self.data: bytes = packet.data
        self.last_chunk_len_format: str = last_chunk_len_format
        self.checksum_len_format: str = checksum_len_format
        self.checksum: typing.Optional[int] = None
        self.additional_payload: typing.Optional[bytes] = None
        self.last_chunk_length, self.file_name = self.decode_header_info()

    def get_last_chunk_length(self) -> int:
        return self.last_chunk_length

    def get_file_name(self) -> typing.Union[str, bytes]:
        return self.file_name

    def decode_header_info(self) -> typing.Tuple[int, typing.Union[bytes, str]]:
        # Size of last Chunk
        # Filename
        # PAD-Bytes
        last_chunk_struct_len: int = struct.calcsize("<" + self.last_chunk_len_format)
        last_chunk_length: int = struct.unpack("<" + self.last_chunk_len_format,
                                               bytes(self.data[0:last_chunk_struct_len]))[0]
        if self.checksum_len_format is not None and self.checksum_len_format != "":
            checksum_struct_len: int = struct.calcsize("<" + self.checksum_len_format)
        else:
            checksum_struct_len = 0
        data: bytes = self.data
        end_of_file_name: int = data.find(0x00, last_chunk_struct_len + 1)
        if end_of_file_name < 0:
            end_of_file_name = len(data)
        if self.checksum_len_format is not None and self.checksum_len_format != "":
            checksum_struct_len: int = struct.calcsize("<" + self.checksum_len_format)
            self.checksum = struct.unpack("<" + self.checksum_len_format,
                                               data[last_chunk_struct_len:last_chunk_struct_len + checksum_struct_len])[
                0]
        file_name: typing.Union[str, bytes] = \
            struct.unpack("<" + str(len(data[last_chunk_struct_len + checksum_struct_len:end_of_file_name])) + "s",
                          data[last_chunk_struct_len + checksum_struct_len:end_of_file_name])[0]
        self.additional_payload = data[end_of_file_name:]  # should be all zero bytes for legacy packets
        return last_chunk_length, file_name

    def update_header(self, filename: str, checksum: int, additional_payload: bytes):
        # We may reduce the filename to "" (indicating "use old filename") while grating us more space for any metadata
        if filename is None:
            filename = ""
        self.file_name = filename
        self.checksum = checksum
        file_name_length = len(filename)
        if len(additional_payload) > 0:
            # +1 as we MUST have a zero-terminated filename!
            padding_str = str(1 + len(self.data) - file_name_length - struct.calcsize(
                "" + self.last_chunk_len_format + self.checksum_len_format + str(len(additional_payload)) + "s")) + "x"
            self.data = struct.pack(
                f"<{self.last_chunk_len_format}{self.checksum_len_format}{file_name_length}s{padding_str}{len(additional_payload)}s",
                self.last_chunk_length, checksum, filename, additional_payload)
        else:
            padding_str = str(len(self.data) - file_name_length - struct.calcsize(
                "" + self.last_chunk_len_format + self.checksum_len_format + str(len(additional_payload)) + "s")) + "x"
            self.data = struct.pack(str("<" + self.last_chunk_len_format + self.checksum_len_format + str(
                file_name_length) + "s" + padding_str),
                                    self.last_chunk_length, checksum, filename)

    @staticmethod
    def from_raw_array(raw_array: NDArray[np.uint8], last_chunk_len_format: str = "I",
                       checksum_len_format: typing.Optional[str] = None) -> 'HeaderChunk':
        packet = Packet(raw_array, {0},
                        total_number_of_chunks=0)  # we use 0 for # chunks as the content as we just need a stub to initialize the HeaderChunk
        return HeaderChunk(packet, last_chunk_len_format, checksum_len_format)

    def get_numpy(self):
        return np.frombuffer(self.data, dtype=np.uint8)

    def __str__(self) -> str:
        return "< last_chunk_length: " + str(self.last_chunk_length) + " , file_name: " + str(self.file_name) + " >"

    def __repr__(self):
        return self.__str__()
