import argparse
import logging
import os
import struct
from io import BytesIO
from math import ceil, floor
from typing import Any, Callable, List, Optional, Union

import numpy as np
from numpy.typing import NDArray
from PIL import Image

from .BPDecoder import BPDecoder
from .distributions.Distribution import Distribution
from .distributions.RaptorDistribution import RaptorDistribution
from .ErrorCorrection import get_error_correction_decode, nocode
from .HeaderChunk import HeaderChunk
from .helper import bitSet, buildGraySequence, logical_xor, xor_mask
from .helper.quaternary2Bin import quad_file_to_bytes, quat_file_to_bin, tranlate_quat_to_byte
from .helper.RU10Helper import choose_packet_numbers, from_true_false_list, intermediate_symbols
from .RU10IntermediatePacket import RU10IntermediatePacket
from .RU10Packet import RU10Packet

logger = logging.getLogger(__name__)


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
        if self.static_number_of_chunks is not None:
            self.number_of_chunks = self.static_number_of_chunks
            number_of_chunks_len_format = (
                ""  # if we got static number_of_chunks we do not need it in struct string
            )
        for filename in os.listdir(self.file):
            if filename.endswith(".RU10") or filename.endswith("DNA"):
                self.EOF = False
                filepath = self.file + "/" + filename
                if filename.endswith("DNA"):
                    if self.error_correction.__name__ == "dna_reed_solomon_decode":
                        try:
                            self.f = quad_file_to_bytes(filepath)
                        except TypeError:
                            logger.warning("skipping CORRUPT file - contains illegal character(s)")
                            self.corrupt += 1
                            continue
                    else:
                        try:
                            self.f = quat_file_to_bin(filepath)
                        except TypeError:
                            logger.warning("skipping CORRUPT file - contains illegal character(s)")
                            self.corrupt += 1
                            continue
                else:
                    self.f = open(filepath, "rb")
                new_pack = self.getNextValidPacket(
                    True,
                    packet_len_format=packet_len_format,
                    crc_len_format=crc_len_format,
                    number_of_chunks_len_format=number_of_chunks_len_format,
                    id_len_format=id_len_format,
                )
                if new_pack is not None:
                    # koennte durch input_new_packet ersetzt werden:
                    # self.addPacket(new_pack)
                    decoded = self.input_new_packet(new_pack)
                if decoded:
                    break
                ##
        logger.info("Decoded Packets: %s", self.correct)
        logger.info("Corrupt Packets: %s", self.corrupt)
        if hasattr(self, "f"):
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
        decoded = False
        self.EOF = False
        if self.file.lower().endswith("dna"):
            try:
                if self.f is not None:
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
            if self.f is not None:
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
                    packet_len_format=packet_len_format,
                    id_len_format=id_len_format,
                )
                if new_pack is not None:
                    decoded = self.input_new_packet(new_pack)
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
                decoded = self.input_new_packet(new_pack)
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
        return logical_xor(aux_mapping)

    def input_new_packet(self, packet: RU10Packet):
        """
        Removes auxpackets (LDPC and Half) and adds the remaining data to the GEPP matrix.
        :param packet: A Packet to add to the GEPP matrix
        :return: True: If solved. False: Else.
        """
        if self.auxBlocks == dict() and self.dist is None:  # self.isPseudo and
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
        assert self.number_of_chunks is not None, (
            "createAuxBlocks can only be called AFTER first Packet"
        )
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
        if not from_multiple_files:
            packet_len = self.f.read(struct.calcsize("<" + packet_len_format))
            try:
                packet_len = struct.unpack("<" + packet_len_format, packet_len)[0]
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
        except:
            self.corrupt += 1
            return None

        data = packet[struct_len:]
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
                raise RuntimeError("Invalid method_data: %s" % method_data)
        len_data = struct.unpack(struct_str, packet[0:struct_len])
        if self.static_number_of_chunks is None:
            self.number_of_chunks = xor_mask(len_data[0], number_of_chunks_len_format)
            unxored_id = xor_mask(len_data[1], id_len_format)
        else:
            unxored_id = xor_mask(len_data[0], id_len_format)
        if self.dist is None:
            self.dist = RaptorDistribution(self.number_of_chunks)
            _, self.s, self.h = intermediate_symbols(self.number_of_chunks, self.dist)

        if self.correct == 0:
            self.createAuxBlocks()
        self.correct += 1
        if self.use_method:
            numbers = choose_packet_numbers(
                len(chunk_lst), unxored_id, self.dist, systematic=False, max_l=len(chunk_lst)
            )
            used_packets = set([chunk_lst[i] for i in numbers])
        else:
            used_packets = set(
                choose_packet_numbers(
                    self.number_of_chunks, unxored_id, self.dist, systematic=False
                )
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
                if bitSet(np.uint32(m[j]), np.uint32(i)):
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
        file_name: str = (
            "DEC_" + os.path.basename(self.file) if self.file is not None else "RU10.BIN"
        )
        output_concat = b""
        if self.headerChunk is not None:
            file_name: str = self.headerChunk.get_file_name().decode("utf-8")
        file_name: str = file_name.split("\x00")[0]
        with open(file_name, "wb") as f:
            a = []
            for decoded in sorted(self.decodedPackets):
                [num] = decoded.get_used_packets()
                if 0 != num or not self.use_headerchunk or self.number_of_chunks - 1 == 0:
                    if isinstance(decoded, RU10IntermediatePacket):
                        a.append(num)
                    if self.number_of_chunks - 1 == num and self.use_headerchunk:
                        output: Union[bytes, np.ndarray] = decoded.get_data()[
                            0 : self.headerChunk.get_last_chunk_length()
                        ]
                        if type(output) == bytes:
                            output_concat += output
                        else:
                            output_concat += output.tobytes()
                        f.write(output)
                    else:
                        if null_is_terminator:
                            data = decoded.get_data()
                            if type(data) == bytes:
                                splitter = data.decode().split("\x00")
                            else:
                                splitter = data.tostring().decode().split("\x00")
                            output = splitter[0].encode()
                            if type(output) == bytes:
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
                            if len(splitter) > 1:
                                break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                        else:
                            output = decoded.get_data()
                            if type(output) == bytes:
                                output_concat += output
                            else:
                                output_concat += output.tobytes()
                            f.write(output)
        logger.info("Saved file as '%s'", file_name)
        if print_to_output:
            print("Result:")
            print(output_concat.decode("utf-8"))
        if return_file_name:
            return file_name

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
        flip_bits: NDArray[np.int_] = np.logical_not(unpack).astype(int)
        new_img = self.draw_img(flip_bits, width, height)
        assert self.file is not None, "filename must be known, not None!"
        tmp_file_name = os.path.basename(self.file) + ".bmp"
        file_name = "DEC_" + tmp_file_name if self.file is not None else "RU10.BIN.bmp"
        new_img.save(file_name)
        return file_name

    @staticmethod
    def draw_img(unpacked_flipped_bits: NDArray[np.int_], width: int, height: int) -> Image.Image:
        new_img = Image.new("1", (width, height))
        pixels = new_img.load()

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
