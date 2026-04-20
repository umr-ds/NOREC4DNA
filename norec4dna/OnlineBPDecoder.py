#!/usr/bin/python
# -*- coding: latin-1 -*-

import os
import struct
import typing
from math import ceil

import numpy
import numpy as np
from norec4dna.BPDecoder import BPDecoder
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.ErrorCorrection import crc32, nocode
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.helper import calc_crc, logical_xor, xor_mask
from norec4dna.helper.quaternary2Bin import quat_file_to_bin
from norec4dna.OnlineAuxPacket import OnlineAuxPacket
from norec4dna.OnlinePacket import OnlinePacket
from numpy.typing import NDArray


class OnlineBPDecoder(BPDecoder):
    def __init__(
        self,
        file: typing.Optional[str] = None,
        error_correction: typing.Callable[[bytes], bytes] = nocode,
        use_headerchunk: bool = True,
        static_number_of_chunks: typing.Optional[int] = None,
    ):
        super().__init__(file, error_correction, use_headerchunk, static_number_of_chunks)
        self.use_headerchunk: bool = use_headerchunk
        self.file: typing.Optional[str] = file
        if file is not None:
            self.isFolder: bool = os.path.isdir(file)
            if not self.isFolder:
                self.f = open(file, "rb")
        self.rng: numpy.random.RandomState = numpy.random.RandomState()
        self.auxBlockNumbers: typing.Dict[int, typing.Set[int]] = dict()
        self.error_correction: typing.Callable[[bytes], bytes] = error_correction
        self.static_number_of_chunks: typing.Optional[int] = static_number_of_chunks
        self.epsilon: typing.Optional[float] = None
        self.quality: typing.Optional[int] = None

    def decodeFolder(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        quality_len_format: str = "I",
        epsilon_len_format: str = "f",
        check_block_number_len_format: str = "I",
    ) -> int:
        effective_number_of_chunks_len_format = number_of_chunks_len_format
        decoded: bool = False
        self.EOF: bool = False
        if self.static_number_of_chunks is not None:
            self.number_of_chunks: int = self.static_number_of_chunks
            effective_number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        if self.file is None:
            return 0
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
                    number_of_chunks_len_format=effective_number_of_chunks_len_format,
                    quality_len_format=quality_len_format,
                    epsilon_len_format=epsilon_len_format,
                    check_block_number_len_format=check_block_number_len_format,
                )
                if new_pack is not None:
                    # koennte durch input_new_packet ersetzt werden:
                    self.addPacket(new_pack)
                    decoded = self.updatePackets(new_pack)
                if decoded:
                    break
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        if hasattr(self, "f"):
            self.f.close()
        if not decoded and self.EOF:
            print("Unable to retrieve file from chunks. Too many errors?")
            return -1
        return 0

    def decodeFile(
        self,
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        number_of_chunks_len_format: str = "I",
        quality_len_format: str = "I",
        epsilon_len_format: str = "f",
        check_block_number_len_format: str = "I",
    ) -> int:
        effective_number_of_chunks_len_format = number_of_chunks_len_format
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            effective_number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        decoded = False
        self.EOF: bool = False
        while not (decoded or self.EOF):
            new_pack = self.getNextValidPacket(
                False,
                packet_len_format=packet_len_format,
                crc_len_format=crc_len_format,
                number_of_chunks_len_format=effective_number_of_chunks_len_format,
                quality_len_format=quality_len_format,
                epsilon_len_format=epsilon_len_format,
                check_block_number_len_format=check_block_number_len_format,
            )
            if new_pack is None:
                break
            self.addPacket(new_pack)
            decoded = self.updatePackets(new_pack)
            ##
        print("Decoded Packets: " + str(self.correct))
        print("Corrupt Packets : " + str(self.corrupt))
        self.f.close()
        if not decoded and self.EOF:
            print("Unable to retrieve file from chunks. Too many errors?")
            return -1
        return 0

    def decodeZip(self, *args: typing.Any, **kwargs: typing.Any) -> typing.Optional[typing.Any]:
        raise RuntimeError("Not implemented for BP-based decoders in the current version!")

    def input_new_packet(self, packet: OnlinePacket, last_chunk_len_format: str = "I") -> bool:
        self.number_of_chunks = packet.total_number_of_chunks
        self.quality = packet.quality
        self.epsilon = packet.epsilon
        if self.headerChunk is None and self.use_headerchunk:
            self.decodeHeader(last_chunk_len_format=last_chunk_len_format)
        if self.correct == 0:
            self.createAuxBlocks()
        self.correct += 1
        self.addPacket(packet)
        return self.updatePackets(packet)

    def getAuxPacketListFromPacket(self, packet: OnlinePacket) -> typing.List[typing.List[bool]]:
        res: typing.List[typing.List[bool]] = []
        aux_used_packets = packet.getBoolArrayAuxPackets()
        i = 0
        for aux in aux_used_packets:
            if aux:
                bool_array = self.auxBlocks[i].get_bool_array_used_packets()
                # Convert numpy array to list if needed, handle None case
                if bool_array is None:
                    res.append([])
                elif hasattr(bool_array, "tolist"):
                    res.append(bool_array.tolist())
                else:
                    res.append(bool_array)
            i += 1
        return res

    def removeAndXorAuxPackets(self, packet: OnlinePacket) -> typing.List[bool]:
        aux_mapping = self.getAuxPacketListFromPacket(packet)
        packet_bool_array = packet.get_bool_array_used_packets()
        # Convert numpy array to list if needed
        if hasattr(packet_bool_array, "tolist"):
            aux_mapping.append(packet_bool_array.tolist())
        else:
            aux_mapping.append(packet_bool_array)
        result = logical_xor(aux_mapping)
        # Ensure result is a list
        if hasattr(result, "tolist"):
            return result.tolist()
        return result if isinstance(result, list) else list(result)

    def createAuxBlocks(self) -> None:
        assert (
            self.number_of_chunks is not None
        ), "createAuxBlocks can only be called AFTER first Packet"
        # self.dist.update_number_of_chunks(self.number_of_chunks)
        self.rng.seed(int(self.number_of_chunks))
        if self.debug:
            print(
                "We should have "
                + str(self.getNumberOfAuxBlocks())
                + " Aux-Blocks and "
                + str(self.number_of_chunks)
                + " normal Chunks (+ 1 HeaderChunk)"
            )
        for i in range(0, self.getNumberOfAuxBlocks()):
            self.auxBlockNumbers[i] = set()
        for chunk_no in range(
            0, self.number_of_chunks
        ):  # + (1 if self.use_headerchunk else 0)):  # + 1 for HeaderChunk
            # Insert this Chunk into quality different Aux-Packets
            for i in range(0, self.quality if self.quality is not None else 0):
                # uniform choose a number of aux blocks
                aux_no = int(self.rng.randint(0, self.getNumberOfAuxBlocks()))
                self.auxBlockNumbers[aux_no].add(chunk_no)

        # XOR all Chunks into the corresponding AUX-Block
        for aux_number in self.auxBlockNumbers.keys():
            self.auxBlocks[aux_number] = OnlineAuxPacket(
                b"",
                self.auxBlockNumbers[aux_number],
                aux_number=aux_number,
                total_number_of_chunks=self.number_of_chunks,
            )  # , numberOfAuxPackets=self.getNumberOfAuxBlocks()) # We will add the Data once we have it.

    def solve(self):
        if self.use_headerchunk and self.headerChunk is None:
            self.decodeHeader()
        # Decoder.solve(self)
        super(OnlineBPDecoder, self).solve()

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
        quality_len_format: str = "I",
        epsilon_len_format: str = "f",
        check_block_number_len_format: str = "I",
    ) -> typing.Optional[OnlinePacket]:
        if not from_multiple_files:
            packet_len = self.f.read(struct.calcsize("<" + packet_len_format))
            packet_len = struct.unpack("<" + packet_len_format, packet_len)[0]
            packet: bytes = self.f.read(int(packet_len))
        else:
            packet = self.f.read()
            packet_len = len(packet)
        if not packet or not packet_len:  # EOF
            self.EOF = True
            self.f.close()
            return None

        crc_len: typing.Optional[int] = struct.calcsize("<" + crc_len_format)
        if self.error_correction.__code__.co_name == crc32.__code__.co_name:
            payload = packet[:crc_len]
            crc = struct.unpack("<L", packet[crc_len:])[0]
            calced_crc = calc_crc(payload)

            if crc != calced_crc:  # If the Packet is corrupt, try next one
                print("[-] CRC-Error - " + str(hex(crc)) + " != " + str(hex(calced_crc)))
                self.corrupt += 1
                return self.getNextValidPacket(from_multiple_files)
        else:
            crc_len = None
            try:
                packet = self.error_correction(packet)
            except:
                self.corrupt += 1
                return self.getNextValidPacket(from_multiple_files)
        struct_str: str = (
            "<"
            + number_of_chunks_len_format
            + quality_len_format
            + epsilon_len_format
            + check_block_number_len_format
        )
        struct_len: int = struct.calcsize(struct_str)
        data = packet[struct_len:crc_len]
        len_data: typing.Tuple = struct.unpack(struct_str, packet[0:struct_len])
        quality: int
        check_block_number: int

        def convert_xor_result(result: typing.Any) -> int:
            """Helper to convert xor_mask result to int"""
            if isinstance(result, np.ndarray):
                return int(result.item())
            elif isinstance(result, (int, float)):
                return int(result)
            else:
                return int(result)

        if self.static_number_of_chunks is None:
            number_of_chunks_raw, quality_raw, self.epsilon, check_block_number_raw = len_data
            number_of_chunks_result = xor_mask(number_of_chunks_raw, number_of_chunks_len_format)
            self.number_of_chunks = convert_xor_result(number_of_chunks_result)
            quality_result = xor_mask(quality_raw, quality_len_format)
            quality = convert_xor_result(quality_result)
            check_block_number = int(check_block_number_raw)
        else:
            quality_raw, self.epsilon, check_block_number_raw = len_data
            quality_result = xor_mask(quality_raw, quality_len_format)
            quality = convert_xor_result(quality_result)
            check_block_number = int(check_block_number_raw)
        if self.dist is None:
            eps_value: float = self.epsilon if self.epsilon is not None else 0.0
            self.dist = OnlineDistribution(eps_value)
        if self.correct == 0:
            # Create MockUp AuxBlocks with the given Pseudo-Random Number -> we will know which Packets are Encoded in which AuxBlock
            self.createAuxBlocks()

        self.correct += 1
        data_array: NDArray[np.uint8] = np.frombuffer(data, dtype=np.uint8)
        epsilon_value: float = self.epsilon if self.epsilon is not None else 0.0
        res = OnlinePacket(
            data_array.tobytes(),
            self.number_of_chunks,
            quality,
            epsilon_value,
            check_block_number,
            read_only=True,
            crc_len_format=crc_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            quality_len_format=quality_len_format,
            epsilon_len_format=epsilon_len_format,
            check_block_number_len_format=check_block_number_len_format,
            save_number_of_chunks_in_packet=self.static_number_of_chunks is None,
        )
        return res

    def decodeHeader(self, last_chunk_len_format: str = "I") -> None:
        if self.headerChunk is not None or 1 not in self.degreeToPacket.keys():
            return  # Header already set or no chunks decoded so far
        for decoded in self.degreeToPacket[1]:
            if decoded.get_used_packets().issubset({0}):
                self.headerChunk = HeaderChunk(decoded, last_chunk_len_format=last_chunk_len_format)
                return

    def saveDecodedFile(
        self, null_is_terminator: bool = False, print_to_output: bool = True
    ) -> None:
        assert self.is_decoded(), "Can not save File: Unable to reconstruct."
        if self.use_headerchunk:
            self.decodeHeader()
        assert self.file is not None, "Cannot save file: file is None"
        file_name = (
            "DEC_" + self.file.split("\x00")[0]
        )  # split is needed for weird  MAC / Windows bugs...
        output_concat = b""
        if self.headerChunk is not None:
            file_name = self.headerChunk.get_file_name().decode("utf-8")
        with open(file_name, "wb") as f:
            a = []
            for decoded in sorted(self.decodedPackets):
                [num] = decoded.get_used_packets()
                if 0 != num or not self.use_headerchunk or self.number_of_chunks - 1 == 0:
                    if isinstance(decoded, OnlineAuxPacket):
                        a.append(num)
                    if self.number_of_chunks - 1 == num and self.use_headerchunk:
                        output = decoded.get_data()[0 : self.headerChunk.get_last_chunk_length()]
                        if isinstance(output, bytes):
                            output_concat += output
                        else:
                            output_concat += output.tobytes()
                        f.write(output)
                    else:
                        if null_is_terminator:
                            data = decoded.get_data()
                            if isinstance(data, bytes):
                                splitter = data.decode().split("\x00")
                            else:
                                splitter = data.tobytes().decode().split("\x00")
                            output = splitter[0].encode()
                            if isinstance(output, bytes):
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
                            if len(splitter) > 1:
                                break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                        else:
                            output = decoded.get_data()
                            if isinstance(output, bytes):
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
        print("Saved file as '" + str(file_name) + "'")
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))

    def getNumberOfAuxBlocks(self) -> int:
        quality = self.quality if self.quality is not None else 0
        epsilon = self.epsilon if self.epsilon is not None else 0.0
        return int(ceil(0.55 * quality * epsilon * self.number_of_chunks))


if __name__ == "__main__":
    # test()
    example_file = "../ONLINE_logo.jpg"
    x = OnlineBPDecoder(example_file)
    x.decode()
    x.decode()
    x.saveDecodedFile()
