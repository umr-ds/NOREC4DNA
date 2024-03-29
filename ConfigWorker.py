import sys
import configparser
import typing

from demo_raptor_decode import demo_decode as demo_raptor_decode
from demo_online_decode import demo_decode as demo_online_decode
from demo_decode import demo_decode as demo_lt_decode
from norec4dna.Decoder import Decoder
from norec4dna.ErrorCorrection import get_error_correction_decode
from norec4dna.helper import find_ceil_power_of_four, fasta_cluster_and_remove_index, cluster_and_remove_index, \
    merge_parts


class ConfigReadAndExecute:
    def __init__(self, filename):
        self.config: configparser.ConfigParser = configparser.ConfigParser()
        self.config.read(filename)
        self.coder = None

    def execute(self, return_decoder: bool = False, skip_solve=False) -> typing.Optional[
        typing.List[typing.Union[str, Decoder]]]:
        decoders = []
        if len(self.config.sections()) == 0:
            print("Empty or missing config file. Does the config file exist?")
        for section in self.config.sections():
            print("Decoding {}".format(section))
            decoder = self.__decode(section, self.config[section], return_decoder, skip_solve=skip_solve)
            if return_decoder:
                decoders.append(decoder)
        return decoders

    def warn_unknown_items(self, config):
        known = ["error_correction", "repair_symbols", "as_mode_1_bmp", "number_of_splits", "split_index_position",
                 "split_index_length", "last_split_smaller", "is_null_terminated", "insert_header", "id_len_format",
                 "number_of_chunks_len_format", "packet_len_format", "crc_len_format", "algorithm", "number_of_chunks",
                 "read_all", "epsilon", "quality", "rules", "quality_len_format", "epsilon_len_format", "config_str",
                 "savenumberofchunks", "dropped_packets", "created_packets", "upper_bound", "asdna", "master_seed",
                 "mode_1_bmp", "chunk_size", "distribution", "checksum_len_str","xor_by_seed", "id_spacing", "mask_id"]
        for cfg in config:
            if cfg not in known:
                print(f"[Warning] Config-entry '{cfg}' not known!")

    # @staticmethod
    def __decode(self, filename, decode_conf, return_decoder=False, skip_solve=False):
        self.warn_unknown_items(decode_conf)
        algorithm = decode_conf.get("algorithm")
        number_of_chunks = decode_conf.getint("number_of_chunks", None)
        e_correction = decode_conf.get("error_correction", "nocode")  # optional
        repair_symbols = decode_conf.getint("repair_symbols", 2)  # optional
        mode_1_bmp = decode_conf.getboolean("as_mode_1_bmp")  # optional, bool
        number_of_splits = decode_conf.getint("number_of_splits", 0)  # optional
        split_index_position = decode_conf.getint("split_index_position", "end")  # optional
        split_index_length = decode_conf.getint("split_index_length")  # optional
        last_split_smaller = decode_conf.getboolean("last_split_smaller")  # optional, bool
        is_null_terminated = decode_conf.getboolean("is_null_terminated")  # optional, bool
        use_header_chunk = decode_conf.getboolean("insert_header")  # optional, bool
        id_len_format = decode_conf.get("id_len_format")  # optional, str
        number_of_chunks_len_format = decode_conf.get("number_of_chunks_len_format", "I")  # optional, str
        packet_len_format = decode_conf.get("packet_len_format", "I")  # optional, str
        crc_len_format = decode_conf.get("crc_len_format", "L")  # optional, str
        read_all_packets = decode_conf.getboolean("read_all", False)
        distribution_cfg_str = decode_conf.get("distribution", "")
        checksum_len_str = decode_conf.get("checksum_len_str", "")
        xor_by_seed = decode_conf.getboolean("xor_by_seed", False)
        id_spacing = decode_conf.getint("id_spacing", 0)
        mask_id = decode_conf.getboolean("mask_id", True)
        # extract preconfig steps:
        if number_of_splits != 0:
            split_index_length = find_ceil_power_of_four(number_of_splits)
        last_split_folder = None
        if split_index_length != 0 and split_index_length is not None:
            if filename.lower().endswith("fasta"):
                folders, last_split_folder = fasta_cluster_and_remove_index(split_index_position, split_index_length,
                                                                            filename)
            else:
                folders, last_split_folder = cluster_and_remove_index(split_index_position, split_index_length,
                                                                      filename)
            # check if the number of folders is equal to the number_of_splits given by user ( if this is != 0 )
            if number_of_splits > 0 and number_of_splits != len(folders):
                print("[WARNING] Number of Splits given by user differs from number of splits found!")
        else:
            folders = [filename]
        error_correction = get_error_correction_decode(e_correction, repair_symbols)
        decoded_files = []
        for f_file in folders:
            print("File / Folder to decode: " + str(f_file))
            if algorithm.lower() == "ru10":
                demo = demo_raptor_decode()
            elif algorithm.lower() == "lt":
                demo = demo_lt_decode()
            elif algorithm.lower() == "online":
                demo = demo_online_decode()
            else:
                raise RuntimeError("unsupported algorithm, this version supports: \"RU10\", \"Online\" and \"LT\"")
            self.coder = demo
            try:
                decoded_files.append(
                    demo.decode(f_file, error_correction=error_correction, null_is_terminator=is_null_terminated,
                                mode_1_bmp=mode_1_bmp, id_len_format=id_len_format,
                                number_of_chunks_len_format=number_of_chunks_len_format,
                                packet_len_format=packet_len_format, crc_len_format=crc_len_format,
                                number_of_chunks=number_of_chunks + (
                                    -1 if f_file == last_split_folder and last_split_smaller else 0),
                                use_header_chunk=use_header_chunk, read_all=read_all_packets,
                                distribution_cfg_str=distribution_cfg_str, return_decoder=return_decoder,
                                checksum_len_str=checksum_len_str, skip_solve=skip_solve, xor_by_seed=xor_by_seed,
                                id_spacing=id_spacing, mask_id=mask_id))
            except Exception as ex:
                raise ex
        if len(folders) > 1:
            merge_parts(decoded_files, remove_tmp_on_success=True)
        else:
            return decoded_files[0]


if __name__ == "__main__":
    if len(sys.argv) > 1:
        file = sys.argv[1]
    else:
        file = "./Dorn_Mon_Oct_10_10_45_15_2022.ini"  # SpringBlossoms.txt_Thu_Nov_26_16_13_20_2020.ini"
    x = ConfigReadAndExecute(file)
    x.execute()
