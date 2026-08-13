#!/usr/bin/python
import argparse
import os
import typing
from configparser import SectionProxy

import numpy as np
from norec4dna.ErrorCorrection import get_error_correction_decode, nocode
from norec4dna.helper import (
    cluster_and_remove_index,
    fasta_cluster_and_remove_index,
    find_ceil_power_of_four,
    merge_parts,
)
from norec4dna.RU10BPDecoder import RU10BPDecoder
from norec4dna.RU10Decoder import RU10Decoder

# Both RU10 decoders expose the same packet-parse / save interface; the BP
# decoder is selected when the config requests solver="bp".
RU10DecoderLike = typing.Union[RU10Decoder, RU10BPDecoder]

STATIC_NUM_CHUNKS = None  # 149
ID_LEN_FORMAT = "I"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "L"
READ_ALL_BEFORE_DECODER = True


class demo_decode:
    @staticmethod
    def _save_partial_decode(decoder: RU10DecoderLike, null_is_terminator: bool):
        return decoder.saveDecodedFile(
            null_is_terminator=null_is_terminator,
            print_to_output=False,
            return_file_name=True,
            partial_decoding=True,
        )

    @staticmethod
    def _cleanup_decoded_output(decoder: RU10DecoderLike) -> None:
        if decoder.headerChunk is None:
            return
        try:
            raw_file_name = decoder.headerChunk.get_file_name()
            file_name = (
                raw_file_name.decode("utf-8") if isinstance(raw_file_name, bytes) else raw_file_name
            )
            os.remove(file_name.split("\x00")[0])
        except FileNotFoundError:
            pass

    @staticmethod
    def _retry_partial_decode(
        decoder: RU10DecoderLike,
        tmp_A,
        tmp_B,
        null_is_terminator: bool,
        failed_repeats: int,
    ):
        res = False
        attempts = 0
        while not res and attempts < failed_repeats:
            print("failed " + str(attempts))
            assert decoder.GEPP is not None
            gepp = decoder.GEPP
            assert len(gepp.A) == len(gepp.b)
            permutation = np.random.permutation(len(gepp.A))
            gepp.A = np.copy(tmp_A[permutation])
            gepp.b = np.copy(tmp_B[permutation])
            decoder.solve(partial=True)
            attempts += 1
            try:
                res = demo_decode._save_partial_decode(decoder, null_is_terminator)
            except ValueError:
                demo_decode._cleanup_decoded_output(decoder)
                res = False
        if not res:
            raise ValueError
        return res

    @staticmethod
    def decode(
        file: str,
        error_correction: typing.Callable[[bytes], bytes] = nocode,
        null_is_terminator: bool = False,
        mode_1_bmp: bool = False,
        number_of_chunks: typing.Optional[int] = STATIC_NUM_CHUNKS,
        use_header_chunk: bool = False,
        id_len_format: str = ID_LEN_FORMAT,
        number_of_chunks_len_format: str = NUMBER_OF_CHUNKS_LEN_FORMAT,
        packet_len_format: str = PACKET_LEN_FORMAT,
        crc_len_format: str = CRC_LEN_FORMAT,
        read_all: bool = READ_ALL_BEFORE_DECODER,
        distribution_cfg_str: str = "",
        return_decoder: bool = False,
        checksum_len_str: typing.Optional[str] = None,
        skip_solve: bool = False,
        failed_repeats: int = 1000,
        xor_by_seed: bool = False,
        id_spacing: int = 0,
        mask_id: bool = True,
        store_parsed_packets: bool = False,
        config_map: typing.Optional[SectionProxy] = None,
    ) -> typing.Union[RU10DecoderLike, bool, bytes, str]:
        use_bp_decoder = False
        if config_map is not None:
            solver_val = str(config_map.get("solver", "")).lower()
            use_bp_val = config_map.getboolean("use_bp", fallback=False)
            algo_val = str(config_map.get("algorithm", "")).lower()
            if (
                solver_val in ("bp", "belief_propagation", "beliefpropagation")
                or use_bp_val
                or algo_val.endswith("bp")
            ):
                use_bp_decoder = True

        if use_bp_decoder:
            print("Belief Propagation Mode")
            x = RU10BPDecoder(
                file,
                use_headerchunk=use_header_chunk,
                error_correction=error_correction,
                static_number_of_chunks=number_of_chunks,
                checksum_len_str=checksum_len_str,
                xor_by_seed=xor_by_seed,
                mask_id=mask_id,
                id_spacing=id_spacing,
                config_map=config_map,
            )
        else:
            print("Pure Gauss-Mode")
            x = RU10Decoder(
                file,
                use_headerchunk=use_header_chunk,
                error_correction=error_correction,
                static_number_of_chunks=number_of_chunks,
                checksum_len_str=checksum_len_str,
                xor_by_seed=xor_by_seed,
                mask_id=mask_id,
                id_spacing=id_spacing,
                config_map=config_map,
            )
        x.read_all_before_decode = read_all
        x.decode(
            id_len_format=id_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            packet_len_format=packet_len_format,
            crc_len_format=crc_len_format,
            store_parsed_packets=store_parsed_packets,
        )
        assert x.GEPP is not None, "GEPP must not be None at this point"
        x.GEPP.insert_tmp()
        tmp_A = np.copy(x.GEPP.A)
        tmp_B = np.copy(x.GEPP.b)
        if not skip_solve:
            x.solve(partial=True)
        if mode_1_bmp:
            return x.mode_1_bmp_decode()
        if return_decoder:
            return x

        try:
            return demo_decode._save_partial_decode(x, null_is_terminator)
        except (FileNotFoundError, ValueError):
            demo_decode._cleanup_decoded_output(x)
            return demo_decode._retry_partial_decode(
                x, tmp_A, tmp_B, null_is_terminator, failed_repeats
            )


def _build_parser() -> argparse.ArgumentParser:
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
    parser.add_argument(
        "--repair_symbols",
        metavar="repair_symbols",
        type=int,
        required=False,
        default=2,
        help="number of repair_symbols for ReedSolomon (default=2)",
    )
    parser.add_argument("--as_mode_1_bmp", required=False, action="store_true")
    parser.add_argument(
        "--number_of_splits",
        metavar="number_of_splits",
        required=False,
        type=int,
        default=0,
        help="(optional) number of parts the file has bin split into",
    )
    parser.add_argument(
        "--split_index_position",
        metavar="split_index_position",
        required=False,
        type=str,
        default="end",
        help="position of the split index. can be 'start' or 'end",
    )
    parser.add_argument(
        "--split_index_length",
        metavar="split_index_length",
        required=False,
        type=int,
        default=0,
        help="number of bases storing the split index",
    )
    parser.add_argument(
        "--last_split_smaller",
        required=False,
        action="store_true",
        help="If set, the number of chunks for the last split will be reduced by 1",
    )
    parser.add_argument(
        "--number_of_chunks",
        metavar="number_of_chunks",
        required=False,
        type=int,
        default=STATIC_NUM_CHUNKS,
        help="static number of chunks (only set this if not stored in each packet)",
    )
    parser.add_argument("--is_null_terminated", required=False, action="store_true")
    parser.add_argument(
        "--header_crc_str", metavar="header_crc_str", required=False, type=str, default=""
    )
    parser.add_argument("--use_header_chunk", required=False, action="store_true")
    parser.add_argument(
        "--failed_repeats",
        metavar="failed_repeats",
        default=1000,
        help="Number of permutations to try if the decoding fails",
        required=False,
        type=int,
    )
    parser.add_argument("--xor_by_seed", required=False, action="store_true")
    parser.add_argument("--id_spacing", metavar="id_spacing", required=False, type=int, default=0)
    return parser


def _prepare_folders(args) -> typing.Tuple[typing.List[str], typing.Optional[str], int]:
    split_index_length = args.split_index_length
    if args.number_of_splits != 0:
        split_index_length = find_ceil_power_of_four(args.number_of_splits)
    last_split_folder = None
    if split_index_length == 0:
        return [args.filename], None, split_index_length

    if args.filename.lower().endswith("fasta"):
        folders, last_split_folder = fasta_cluster_and_remove_index(
            args.split_index_position, split_index_length, args.filename
        )
    else:
        folders, last_split_folder = cluster_and_remove_index(
            args.split_index_position, split_index_length, args.filename
        )
    if args.number_of_splits > 0 and args.number_of_splits != len(folders):
        print("[WARNING] Number of Splits given by user differs from number of splits found!")
    return folders, last_split_folder, split_index_length


def _number_of_chunks_for_folder(
    args, folder: str, last_split_folder: typing.Optional[str]
) -> typing.Optional[int]:
    if args.number_of_chunks is None:
        return None
    return args.number_of_chunks + (
        -1 if folder == last_split_folder and args.last_split_smaller else 0
    )


def _run_cli() -> None:
    args = _build_parser().parse_args()
    folders, last_split_folder, _ = _prepare_folders(args)
    error_correction = get_error_correction_decode(args.error_correction, args.repair_symbols)
    decoded_files = []
    demo = demo_decode()
    for folder in folders:
        print("File / Folder to decode: " + str(folder))
        try:
            decoded_files.append(
                demo.decode(
                    folder,
                    error_correction=error_correction,
                    null_is_terminator=args.is_null_terminated,
                    mode_1_bmp=args.as_mode_1_bmp,
                    number_of_chunks=_number_of_chunks_for_folder(args, folder, last_split_folder),
                    use_header_chunk=args.use_header_chunk,
                    checksum_len_str=args.header_crc_str,
                    failed_repeats=args.failed_repeats,
                    xor_by_seed=args.xor_by_seed,
                    id_spacing=args.id_spacing,
                )
            )
        except (FileNotFoundError, ValueError):
            continue
    if len(folders) > 1:
        merge_parts(decoded_files, remove_tmp_on_success=True)


if __name__ == "__main__":
    _run_cli()
