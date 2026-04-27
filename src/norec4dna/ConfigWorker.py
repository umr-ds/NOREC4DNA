import configparser
import logging
import sys
import typing

from .demo_decode import demo_decode as demo_lt_decode
from .demo_online_decode import demo_decode as demo_online_decode
from .demo_raptor_decode import demo_decode as demo_raptor_decode
from norec4dna.Decoder import Decoder
from norec4dna.ErrorCorrection import get_error_correction_decode
from norec4dna.helper import (
    cluster_and_remove_index,
    fasta_cluster_and_remove_index,
    find_ceil_power_of_four,
    merge_parts,
)

DecoderResult = typing.Union[str, bytes, Decoder]
DecodeCallableResult = typing.Union[DecoderResult, Decoder]
logger = logging.getLogger(__name__)


class DecoderDemo(typing.Protocol):
    @typing.overload
    def decode(
        self,
        file: str,
        *,
        error_correction: typing.Callable[[bytes], bytes],
        null_is_terminator: bool = False,
        mode_1_bmp: bool = False,
        number_of_chunks: typing.Optional[int] = None,
        use_header_chunk: bool = False,
        id_len_format: str = "",
        number_of_chunks_len_format: str = "I",
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        read_all: bool = False,
        distribution_cfg_str: str = "",
        return_decoder: typing.Literal[True],
        checksum_len_str: str = "",
        skip_solve: bool = False,
        xor_by_seed: bool = False,
        id_spacing: int = 0,
        mask_id: bool = True,
        store_parsed_packets: bool = False,
        config_map: configparser.SectionProxy,
    ) -> Decoder: ...

    @typing.overload
    def decode(
        self,
        file: str,
        *,
        error_correction: typing.Callable[[bytes], bytes],
        null_is_terminator: bool = False,
        mode_1_bmp: bool = False,
        number_of_chunks: typing.Optional[int] = None,
        use_header_chunk: bool = False,
        id_len_format: str = "",
        number_of_chunks_len_format: str = "I",
        packet_len_format: str = "I",
        crc_len_format: str = "L",
        read_all: bool = False,
        distribution_cfg_str: str = "",
        return_decoder: bool = False,
        checksum_len_str: str = "",
        skip_solve: bool = False,
        xor_by_seed: bool = False,
        id_spacing: int = 0,
        mask_id: bool = True,
        store_parsed_packets: bool = False,
        config_map: configparser.SectionProxy,
    ) -> DecoderResult: ...


class ConfigReadAndExecute:
    def __init__(self, filename: str) -> None:
        self.config: configparser.ConfigParser = configparser.ConfigParser()
        self.config.read(filename)
        self.coder: typing.Any = None

    @typing.overload
    def execute(
        self,
        return_decoder: typing.Literal[True],
        skip_solve: bool = False,
        store_parsed_packets: bool = False,
    ) -> typing.Sequence[Decoder]: ...

    @typing.overload
    def execute(
        self,
        return_decoder: typing.Literal[False] = False,
        skip_solve: bool = False,
        store_parsed_packets: bool = False,
    ) -> typing.Sequence[DecoderResult]: ...

    def execute(
        self,
        return_decoder: bool = False,
        skip_solve: bool = False,
        store_parsed_packets: bool = False,
    ) -> typing.Sequence[DecoderResult]:
        decoders: typing.List[DecoderResult] = []
        if len(self.config.sections()) == 0:
            logger.warning("Empty or missing config file. Does the config file exist?")
        for section in self.config.sections():
            logger.info("Decoding %s", section)
            decoder = self.__decode(
                section,
                self.config[section],
                return_decoder,
                skip_solve=skip_solve,
                store_parsed_packets=store_parsed_packets,
            )
            if return_decoder and decoder is not None:
                decoders.append(decoder)
        return decoders

    def warn_unknown_items(self, config: configparser.SectionProxy) -> None:
        known = [
            "error_correction",
            "repair_symbols",
            "as_mode_1_bmp",
            "number_of_splits",
            "split_index_position",
            "split_index_length",
            "last_split_smaller",
            "is_null_terminated",
            "insert_header",
            "id_len_format",
            "number_of_chunks_len_format",
            "packet_len_format",
            "crc_len_format",
            "algorithm",
            "number_of_chunks",
            "read_all",
            "epsilon",
            "quality",
            "rules",
            "quality_len_format",
            "epsilon_len_format",
            "config_str",
            "savenumberofchunks",
            "dropped_packets",
            "created_packets",
            "upper_bound",
            "asdna",
            "master_seed",
            "mode_1_bmp",
            "chunk_size",
            "distribution",
            "checksum_len_str",
            "xor_by_seed",
            "id_spacing",
            "mask_id",
            "last_chunk_len_str",
        ]
        for cfg in config:
            if cfg not in known:
                logger.warning("Config-entry '%s' not known!", cfg)

    @typing.overload
    def __decode(
        self,
        filename: str,
        decode_conf: configparser.SectionProxy,
        return_decoder: typing.Literal[True],
        skip_solve: bool = False,
        store_parsed_packets: bool = False,
    ) -> typing.Optional[Decoder]: ...

    @typing.overload
    def __decode(
        self,
        filename: str,
        decode_conf: configparser.SectionProxy,
        return_decoder: typing.Literal[False] = False,
        skip_solve: bool = False,
        store_parsed_packets: bool = False,
    ) -> typing.Optional[DecoderResult]: ...

    def __decode(
        self,
        filename: str,
        decode_conf: configparser.SectionProxy,
        return_decoder: bool = False,
        skip_solve: bool = False,
        store_parsed_packets: bool = False,
    ) -> typing.Optional[DecoderResult]:
        self.warn_unknown_items(decode_conf)
        algorithm = decode_conf.get("algorithm", "")
        number_of_chunks = decode_conf.getint("number_of_chunks", fallback=None)
        e_correction = decode_conf.get("error_correction", "nocode")
        repair_symbols = decode_conf.getint("repair_symbols", fallback=2)
        mode_1_bmp = decode_conf.getboolean("as_mode_1_bmp", fallback=False)
        number_of_splits = decode_conf.getint("number_of_splits", fallback=0)
        split_index_position = decode_conf.get("split_index_position", "end")
        split_index_length = decode_conf.getint("split_index_length", fallback=0)
        last_split_smaller = decode_conf.getboolean("last_split_smaller", fallback=False)
        is_null_terminated = decode_conf.getboolean("is_null_terminated", fallback=False)
        use_header_chunk = decode_conf.getboolean("insert_header", fallback=False)
        id_len_format = decode_conf.get("id_len_format", "")
        number_of_chunks_len_format = decode_conf.get("number_of_chunks_len_format", "I")
        packet_len_format = decode_conf.get("packet_len_format", "I")
        crc_len_format = decode_conf.get("crc_len_format", "L")
        read_all_packets = decode_conf.getboolean("read_all", fallback=False)
        distribution_cfg_str = decode_conf.get("distribution", "")
        checksum_len_str = decode_conf.get("checksum_len_str", "")
        xor_by_seed = decode_conf.getboolean("xor_by_seed", fallback=False)
        id_spacing = decode_conf.getint("id_spacing", fallback=0)
        mask_id = decode_conf.getboolean("mask_id", fallback=True)
        # extract preconfig steps:
        if number_of_splits != 0:
            split_index_length = find_ceil_power_of_four(number_of_splits)
        last_split_folder: typing.Optional[str] = None
        folders: typing.List[str]
        if split_index_length != 0:
            if filename.lower().endswith("fasta"):
                folders, last_split_folder = fasta_cluster_and_remove_index(
                    split_index_position, split_index_length, filename
                )
            else:
                folders, last_split_folder = cluster_and_remove_index(
                    split_index_position, split_index_length, filename
                )
            # check if the number of folders is equal to the number_of_splits given by user ( if this is != 0 )
            if number_of_splits > 0 and number_of_splits != len(folders):
                logger.warning(
                    "Number of Splits given by user differs from number of splits found!"
                )
        else:
            folders = [filename]
        error_correction = get_error_correction_decode(e_correction, repair_symbols)
        decoded_files: typing.List[typing.Any] = []
        for f_file in folders:
            logger.info("File / Folder to decode: %s", f_file)
            demo: DecoderDemo
            if algorithm.lower() == "ru10":
                demo = typing.cast(DecoderDemo, demo_raptor_decode())
            elif algorithm.lower() == "lt":
                demo = typing.cast(DecoderDemo, demo_lt_decode())
            elif algorithm.lower() == "online":
                demo = typing.cast(DecoderDemo, demo_online_decode())
            else:
                raise RuntimeError(
                    'unsupported algorithm, this version supports: "RU10", "Online" and "LT"'
                )
            self.coder = demo
            try:
                decoded_result = demo.decode(
                    f_file,
                    error_correction=error_correction,
                    null_is_terminator=is_null_terminated,
                    mode_1_bmp=mode_1_bmp,
                    id_len_format=id_len_format,
                    number_of_chunks_len_format=number_of_chunks_len_format,
                    packet_len_format=packet_len_format,
                    crc_len_format=crc_len_format,
                    number_of_chunks=(
                        number_of_chunks
                        + (-1 if f_file == last_split_folder and last_split_smaller else 0)
                        if number_of_chunks is not None
                        else None
                    ),
                    use_header_chunk=use_header_chunk,
                    read_all=read_all_packets,
                    distribution_cfg_str=distribution_cfg_str,
                    return_decoder=return_decoder,
                    checksum_len_str=checksum_len_str,
                    skip_solve=skip_solve,
                    xor_by_seed=xor_by_seed,
                    id_spacing=id_spacing,
                    mask_id=mask_id,
                    store_parsed_packets=store_parsed_packets,
                    config_map=decode_conf,
                )
                decoded_files.append(decoded_result)
            except Exception as ex:
                raise ex
        if len(folders) > 1:
            merge_parts(typing.cast(typing.List[str], decoded_files), remove_tmp_on_success=True)
        else:
            return typing.cast(DecoderResult, decoded_files[0])


if __name__ == "__main__":
    if len(sys.argv) > 1:
        file = sys.argv[1]
    else:
        file = "./Dorn_Mon_Oct_10_10_45_15_2022.ini"  # SpringBlossoms.txt_Thu_Nov_26_16_13_20_2020.ini"
    x = ConfigReadAndExecute(file)
    x.execute()
