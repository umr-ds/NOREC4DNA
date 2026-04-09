"""
This tool should allow a user to:
1) Decode a file encoded with NOREC4DNA
2) if there are not enough packets to decode the file, the user should get:
    - a list of missing chunks
    - a partial result with \x00 for missing chunks
    - ideally a ranking of the missing chunks based on how many additional chunks could be retrived if it was present
3) view the file (either as hex, image or as a text) and manually select corrupt chunks
    - based on the selected chunks the tool will then suggest which packet(s) might have caused the corruption
    - the used can then request a new decoding with the detected packet removed

Automatic mode:
1) if there are multiple packets with the same packet-id (or very close hamming distance in total):
    - the tool should try each combination of these packets
    - if there are (multiple) checksums in the header chunks, the tool could automatically find the corrupt packets and either:
        - remove them from the decoding because there are still enough packets left to decode the file
        - bruteforce the corrupt chunks until the checksums match (this can be done in parallel and using believe propagation)
2) if there is only a single packet with this id:
    - the tool can only try to bruteforce the corrupt chunks / packets:
        IF WE BRUTEFORCE THE CHUNK WE MIGHT HAVE A PROBLEM IF THE PACKET HAD A MUTATION AT THE START (wrong ID!)
            we can avoid this pitfall by NOT using the chunk-mapping of the corrupt packet!
        IF WE BRUTEFORCE THE PACKET WE CANT DIRECTLY USE THE CRC (we must always perform a belief propagation / gauss elimination) - this is slower
"""
import os
import shutil
import typing
from functools import reduce
from io import BytesIO
from itertools import combinations
from pathlib import Path
from time import sleep
import numpy as np
import magic
import crcmod
from norec4dna.GEPP import GEPP
from norec4dna.helper import xor_numpy

import NOREC4DNA.norec4dna.helper as helper
from NOREC4DNA.ConfigWorker import ConfigReadAndExecute
from norec4dna.HeaderChunk import HeaderChunk
from norec4dna.Packet import Packet
from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.OnlineDecoder import OnlineDecoder
from norec4dna.LTDecoder import LTDecoder
from numpy.linalg import matrix_rank


class SemiAutomaticReconstructionToolkit:
    def __init__(self, decoder: typing.Union[RU10Decoder, LTDecoder, OnlineDecoder]):
        self.last_chunk_len_format = "I"
        self.checksum_len_format = None
        self.decoder: typing.Union[RU10Decoder, LTDecoder, OnlineDecoder] = decoder
        decoder.read_all_before_decode = True
        self.headerChunk: typing.Optional[HeaderChunk] = None
        self.decoder.GEPP.insert_tmp()
        self.initial_A = self.decoder.GEPP.A.copy()
        self.initial_b = self.decoder.GEPP.b.copy()
        self.initial_packet_mapping = None  # self.decoder.GEPP.packet_mapping.copy()
        self.multi_error_packets_mode = False

    def calculate_rank_A(self):
        return np.linalg.matrix_rank(self.decoder.GEPP.A)

    def calculate_rank_augmented_matrix(self):
        return np.linalg.matrix_rank(np.c_[self.decoder.GEPP.A, self.decoder.GEPP.b])

    def set_multi_error_mode(self, mode: bool):
        """
        Set the multi error mode
        @param mode: True if the multi error mode should be enabled
        """
        self.multi_error_packets_mode = mode

    def calculate_unused_packets(self):
        # return all packets that are still possible erroneous after all chunks are tagged as valid
        return self.decoder.GEPP.get_common_packets([], [i for i in range(
            self.decoder.number_of_chunks)])

    def manual_repair(self, chunk_id, corrupt_packet_id, repaired_content):
        """
        Repair all chunks by propagating the repaired content of chunk "chunk_id"
        @param chunk_id: id of the chunk to repair
        @param corrupt_packet_id: id of the corrupt packet
        @param repaired_content: content of the chunk to repair
        """
        # xor the corrupt chunk with the repaired chunk
        chunk_diff = helper.xor_numpy(self.decoder.GEPP.b[chunk_id], repaired_content)

        for i in range(self.decoder.GEPP.chunk_to_used_packets.shape[0]):
            if self.decoder.GEPP.chunk_to_used_packets[i, corrupt_packet_id]:
                # if the corrupt packet is used by this chunk
                # xor the chunk with the chunk_diff
                self.decoder.GEPP.b[i] = helper.xor_numpy(self.decoder.GEPP.b[i], chunk_diff)

    def constrained_repair(self, error_delta, possible_packets, working_dir="constrained_repaired"):
        """
        Try to repair the file in the case of multiple possible corrupt packets with a known error delta
        If there is a header chunk with a checksum, the checksum will be used to verify the repair, otherwise all
        solutions will be stored and the user can select the correct one.
        @param error_delta: the error delta
        @param possible_packets: list of possible corrupt packets
        """
        bkp_b = self.decoder.GEPP.b.copy()
        res = ""
        if Path(working_dir).exists():
            shutil.rmtree(working_dir)
        # create the folder working_dir:
        Path(working_dir).mkdir(parents=True, exist_ok=True)
        for corrupt_packet_id in possible_packets:
            for i in range(self.decoder.GEPP.chunk_to_used_packets.shape[0]):
                if self.decoder.GEPP.chunk_to_used_packets[i, corrupt_packet_id]:
                    # if the corrupt packet is used by this chunk
                    # xor the chunk with the chunk_diff
                    self.decoder.GEPP.b[i] = helper.xor_numpy(self.decoder.GEPP.b[i], error_delta)
            self.parse_header(self.last_chunk_len_format, checksum_len_format=self.checksum_len_format)
            if self.headerChunk is not None and self.headerChunk.checksum_len_format is not None:
                is_correct = self.is_checksum_correct()
            else:
                is_correct = False
            try:
                filename = self.decoder.saveDecodedFile(return_file_name=True, print_to_output=False)
            except ValueError as ve:
                filename = ve.args[1]
            _file = Path(filename)
            stem = ("CORRECT_" if is_correct else "") + _file.stem + f"_{corrupt_packet_id}"
            new_filename = Path(working_dir + "/" + stem + _file.suffix)
            new_filename.parents[0].mkdir(parents=True, exist_ok=True)
            _file = _file.rename(new_filename)
            res += f"{_file.name}, "
            self.decoder.GEPP.b = bkp_b.copy()
            if is_correct:
                return f"Found a correct CRC. Saved to folder {working_dir}: {filename}"
        return f"Saved {len(possible_packets)} results to folder {working_dir}: {res}"

    def get_possible_invalid_chunks_from_common_packets(self, _common_packets: typing.List[bool]) -> typing.List[bool]:
        """
        Get a list of possible invalid chunks from the given list of invalid chunks by calculating which chunks use
        the same potentially invalid packets
        @param _common_packets: list of invalid chunks
        @return: list of possible invalid chunks

        """
        # iterate over the inverse mapping and add each row number if it uses an invalid packet
        res = [False] * self.decoder.GEPP.A.shape[0]
        for i, row_content in enumerate(self.decoder.GEPP.chunk_to_used_packets):
            if any(np.logical_and(_common_packets, row_content)):
                # if there is any overlap between the invalid packets and the packets used by this chunk
                # then this chunk could also be invalid
                res[i] = True
        return res

    def predict_file_type(self):
        """
        Predict the file type/informations based on the magic package
        """
        return magic.from_buffer(self.get_file_as_bytes())

    def remove_packet_from_equation_system(self, packet_id: int):
        """
        Remove the given packet from the equation system
        @param packet_id: id of the packet to remove, corresponds to the column in the inverse matrix
        """
        pass

    def repair_by_exclusion(self, comm_packet):
        # speedup: check if the matrix is still solvable after we remove all packets that are invalid
        tmp_gepp: GEPP = GEPP(self.initial_A.copy(), self.initial_b.copy())
        for i, is_common in reversed(list(enumerate(comm_packet))):
            if is_common:
                tmp_gepp.remove_row(i)
        try:
            res = tmp_gepp.isPotentionallySolvable() and not any(tmp_gepp.find_missing_chunks()) and tmp_gepp.solve()
        except:
            res = False
        return res, tmp_gepp

    def all_solutions_by_reordering(self, comm_packet, only_possible_invalid_packets=False):
        # speedup: we might want to check if the matrix is still solvable after we remove all packets that are invalid
        mapping = dict()
        # if only_possible_invalid_packets is True: remove all packets _i_ where comm_packets[i] is True
        if only_possible_invalid_packets:
            # count the number of True in comm_packet:
            valid_packets = len([i for i in comm_packet if not i])
            mod_tmp_gepp: GEPP = GEPP(self.initial_A.copy(), self.initial_b.copy())
            prob_invalid_a = []
            prob_invalid_b = []
            # remove all possible invalid packets from the temporary GEPP
            for i in np.arange(len(mod_tmp_gepp.A) - 1, -1, -1):
                if comm_packet[i]:
                    prob_invalid_a.append(mod_tmp_gepp.A[i])
                    prob_invalid_b.append(mod_tmp_gepp.b[i])
                    mod_tmp_gepp.remove_row(i)
            # add the possibly invalid packets to the end of the matrix
            for i in np.arange(len(prob_invalid_a) - 1, -1, -1):
                mod_tmp_gepp.addRow(prob_invalid_a[i], prob_invalid_b[i])
            mod_tmp_gepp.insert_tmp()
        else:
            valid_packets = 0

        for i in np.arange(valid_packets, self.decoder.GEPP.A.shape[0]):
            if only_possible_invalid_packets:
                tmp_gepp: GEPP = GEPP(mod_tmp_gepp.A.copy(), mod_tmp_gepp.b.copy())
            else:
                tmp_gepp: GEPP = GEPP(self.initial_A.copy(), self.initial_b.copy())
            a_row = tmp_gepp.A[i].copy()
            b_row = tmp_gepp.b[i].copy()
            tmp_gepp.remove_row(i)
            tmp_gepp.addRow(a_row, b_row)
            try:
                res = tmp_gepp.isPotentionallySolvable() and tmp_gepp.solve()
            except:
                res = False
            if res:
                mapping[i] = tmp_gepp
        return mapping

    def view_file_with_chunkborders(self, as_hex: bool = False, null_is_terminator=False,
                                    last_chunk_len_format: str = "I", add_line_numbers=False, checksum_len_format=None):
        """
        shows the content of decoder.b with borders after every n-th symbol
        """
        self.checksum_len_format = checksum_len_format
        self.last_chunk_len_format = last_chunk_len_format
        if self.decoder.GEPP is not None:
            if self.initial_A is None:
                # create an inital backup of the GEPP
                self.initial_A = self.decoder.GEPP.A.copy()
                self.initial_b = self.decoder.GEPP.b.copy()
                self.initial_packet_mapping = self.decoder.GEPP.packet_mapping.copy()
            self.decoder.solve(partial=True)
        dirty = False
        self.parse_header(last_chunk_len_format, checksum_len_format=checksum_len_format)
        file_name = "DEC_" + os.path.basename(self.decoder.file) if self.decoder.file is not None else "RU10.BIN"
        if self.headerChunk is not None:
            try:
                try:
                    file_name = self.headerChunk.get_file_name().decode("utf-8")
                except UnicodeDecodeError:
                    raise RuntimeError("Filename in headerchunk is not utf-8 encoded!")
                    # file_name = self.headerChunk.get_file_name().decode("latin-1")
                if self.headerChunk.data[-1] != 0x00:
                    raise RuntimeError("Headerchunk is not null terminated!" +
                                       "Either the headerchunk is corrupt or no headerchunk was used!")
            except RuntimeError as ex:
                print("Warning:", ex)
        file_name = file_name.split("\x00")[0]
        res = []
        for x in self.decoder.GEPP.result_mapping:
            if x < 0:
                res.append(b"\x00" * len(self.decoder.GEPP.b[x][0]))
                dirty = True
                continue
            if self.decoder.number_of_chunks - 1 == x and self.decoder.use_headerchunk and self.headerChunk is not None:
                # to show the last chunk padding remove: " self.headerChunk.get_last_chunk_length()":
                output = self.decoder.GEPP.b[x][0][0: self.headerChunk.get_last_chunk_length()]
                res.append(output)
            else:
                if null_is_terminator:
                    splitter = self.decoder.GEPP.b[x].tostring().split("\x00")
                    output = splitter[0].encode()
                    res.append(output)
                    if len(splitter) > 1:
                        break  # since we are in null-terminator mode, we exit once we see the first 0-byte
                else:
                    output = self.decoder.GEPP.b[x]
                    res.append(output)
        if dirty:
            print("Some parts could not be restored, file WILL contain sections with \\x00 !")
        ret = []
        for j, line in enumerate(res):
            try:
                line = line.tobytes()
            except Exception:
                pass
            s1 = " ".join([f"{i:02x}" for i in line])
            try:
                width = res[0].shape[1]
            except Exception:  # if first row is not decoded, it will be of type bytes!
                width = len(res[0])
            s2 = "".join([chr(i) if 32 <= i < 127 else "." for i in line])
            ret.append((f"{j:08x} | " if add_line_numbers else "") + f"{s1: <{width * 3}}  |{s2: <{width}}|")
        return ret

    def get_corrupt_chunks_by_packets(self, packets, chunk_tag=None, tag_num=1):
        """
        returns a list of all chunks that are affected by the given packets
        """
        if chunk_tag is None:
            chunk_tag = np.zeros(self.decoder.GEPP.chunk_to_used_packets.shape[0], dtype=np.uint8)
        for packet in packets:
            for i, is_common in enumerate(self.decoder.GEPP.chunk_to_used_packets[:, int(packet)]):
                if is_common:
                    chunk_tag[i] = tag_num
        return chunk_tag

    def bruteforce_repair(self, error_matrix):
        """
        Bruteforce the corrupt chunks
        @param error_matrix: the error matrix

        """

        pass

    def is_checksum_correct(self) -> bool:
        crc_len_str = self.headerChunk.checksum_len_format
        if crc_len_str == "B":
            algo = crcmod.predefined.mkPredefinedCrcFun("crc-8")
        elif crc_len_str == "H":
            algo = crcmod.predefined.mkCrcFun('crc-16')
        elif crc_len_str == "I":
            algo = crcmod.predefined.mkCrcFun('crc-32')  # zlib.crc32
        else:
            raise ValueError("crc_len_str must be one of B, H, I")

        f = BytesIO(self.decoder.GEPP.b[1:].reshape(-1)[
                    :self.headerChunk.last_chunk_length])
        checksum = 0
        while chunk := f.read():
            checksum = algo(chunk, checksum)
        return checksum == self.headerChunk.checksum

    def get_file_as_bytes(self):
        """
        Return the chunk content as bytes
        """
        start = 1 if self.decoder.use_headerchunk else 0
        return self.decoder.GEPP.b[start:].reshape(-1).tobytes()

    def parse_header(self, last_chunk_len_format, checksum_len_format=None):
        if self.decoder.use_headerchunk:
            header_row = self.decoder.GEPP.result_mapping[0]
            if header_row >= 0:
                self.headerChunk = HeaderChunk(
                    Packet(self.decoder.GEPP.b[header_row], {0}, self.decoder.number_of_chunks, read_only=True),
                    last_chunk_len_format=last_chunk_len_format, checksum_len_format=checksum_len_format)

    @staticmethod
    def solve_lin_dep(a, b):
        """
        Calculates which rows in vector a can be used to create the target b
        @param a: a matrix , where each row is either used to create b or not
        @param b: the target vector
        @return: a list of rows in a that can be used to create b or None if no solution exists
        """
        combs = [[x for x in combinations(a, i)] for i in range(1, min(4, len(a) + 1))]
        for comb in combs:
            for elem in comb:
                if len(elem) > 1:
                    r = reduce(lambda x, y: xor_numpy(x.astype("uint8"), y.astype("uint8")), elem)
                else:
                    r = elem[0]
                if np.array_equal(r.astype('uint8'), b):
                    return [x.astype("uint8") for x in elem]
        return None



    def repair_and_store_by_packet(self, chunk_id, packet_id, hex_value, clear_working_dir=False, correctness_function=None):
        # this function will be used if we have multiple invalid packets (and corrected chunks) to save multiple version,
        # where each saved version used a different possible packet to repair the chunk.
        bkp_A = self.decoder.GEPP.A.copy()
        bkp_b = self.decoder.GEPP.b.copy()
        self.manual_repair(chunk_id, packet_id, hex_value)
        working_dir = "multi_file_repair"
        if clear_working_dir:
            # delete the folder working_dir if it exists:
            if Path(working_dir).exists():
                shutil.rmtree(working_dir)
            # create the folder working_dir:
            Path(working_dir).mkdir(parents=True, exist_ok=True)
        # we might have to check if header chunk is used!
        self.parse_header("I")
        if self.headerChunk is not None and self.headerChunk.checksum_len_format is not None:
            is_correct = self.is_checksum_correct()
        else:
            if correctness_function is not None:
                is_correct = correctness_function(self.decoder.GEPP.b)
            else:
                is_correct = False
        try:
            filename = self.decoder.saveDecodedFile(return_file_name=True, print_to_output=False)
        except ValueError as ve:
            filename = ve.args[1]
        _file = Path(filename)
        stem = ("CORRECT_" if is_correct else "") + _file.stem + f"_{chunk_id}_{packet_id}"
        _new_file = _file.rename(Path(working_dir + "/" + stem + _file.suffix))
        self.decoder.GEPP.A = bkp_A
        self.decoder.GEPP.b = bkp_b
        return f"{_new_file.name}"


if __name__ == "__main__":
    x = ConfigReadAndExecute("NOREC4DNA/logo.jpg_Fri_Jan__7_13_18_39_2022.ini").execute(return_decoder=True)[0]
    semi_automatic_solver = SemiAutomaticReconstructionToolkit(x)
    print(semi_automatic_solver.view_file_with_chunkborders(False, False, "I"), flush=True)

    sleep(1)
    print("Enter the rows that are INVALID (as hex; separated by space): ")
    invalid_rows = input().split(" ")
    invalid_rows = [int(i, 16) for i in invalid_rows]

    print("Enter the rows that are VALID (as hex; separated by space): ")
    valid_rows = input().split(" ")
    valid_rows = [int(i, 16) for i in valid_rows]

    common_packets = semi_automatic_solver.decoder.GEPP.get_common_packets(invalid_rows, valid_rows)
    print("potentially invalid Packets:")
    print(" ".join(map(lambda x: "1" if x else "0", common_packets)), flush=True)
    while np.count_nonzero(common_packets == True) > 1:
        rem_possible_chunks = semi_automatic_solver.get_possible_invalid_chunks_from_common_packets(common_packets)
        print("possible invalid chunks:")
        print(" ".join(map(lambda _x: f"{_x[0]:08x}" if _x[1] else "_", enumerate(rem_possible_chunks))), flush=True)

        print(
            "Result unambiguous, enter additional rows that are INVALID (as hex; separated by space), if there are none, just hit [ENTER]: ",
            flush=True)
        tmp_invalid_rows = input()
        if len(tmp_invalid_rows) != 0:
            for new_invalid_line in tmp_invalid_rows.split(" "):
                invalid_rows.append(int(new_invalid_line, 16))
        print(
            "Result unambiguous, enter additional rows that are VALID (as hex; separated by space), if there are none, just hit [ENTER]: ",
            flush=True)
        tmp_valid_rows = input()
        if len(tmp_valid_rows) != 0:
            for new_valid_line in tmp_valid_rows.split(" "):
                valid_rows.append(int(new_valid_line, 16))
        common_packets = semi_automatic_solver.decoder.GEPP.get_common_packets(invalid_rows, valid_rows)
        print(" ".join(map(lambda _X: "1" if _X else "0", common_packets)), flush=True)
        if len(tmp_valid_rows) == 0 and len(tmp_invalid_rows) == 0:
            break
    print("Missing chunks:")
    print(" ".join(map(lambda _x: "1" if _x else "0", semi_automatic_solver.decoder.GEPP.find_missing_chunks())),
          flush=True)
