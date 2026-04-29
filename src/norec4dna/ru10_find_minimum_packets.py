from __future__ import annotations

import argparse
import bisect
import glob
import multiprocessing
import os
import random
from functools import partial
from math import ceil, floor

from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.Encoder import Encoder
from norec4dna.ErrorCorrection import get_error_correction_encode, nocode, reed_solomon_encode
from norec4dna.helper import (
    find_ceil_power_of_four,
    merge_folder_content,
    number_to_base_str,
    should_drop_packet,
    split_file,
)
from norec4dna.Packet import ParallelPacket
from norec4dna.RU10Encoder import RU10Encoder

from .rules.FastDNARules import FastDNARules

ID_LEN_FORMAT = "I"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
DEFAULT_CHUNK_SIZE = 34
DEFAULT_SAVE_AS_FASTA = True
NUMBER_OF_PACKETS_TO_CREATE = 655360


def old_main(
    file=".INFILES/logo.jpg",
    asdna=True,
    insert_header=True,
    error_correction=nocode,
    save_number_of_chunks_in_packet=False,
):
    _ = DEFAULT_CHUNK_SIZE  # Chunksize 50, mit reedsolomon auf 6 repair_symbols


def run(
    seq_seed=None,
    file=".INFILES/logo.jpg",
    asdna=True,
    insert_header=False,
    error_correction=reed_solomon_encode,
    save_number_of_chunks_in_packet=False,
    l_size=1000,
    while_count=1000,
    chunk_size=0,
    number_of_chunks=300,
    prepend="",
    append="",
    seed_len_format=ID_LEN_FORMAT,
    drop_above=1.0,
    e_correction_name="nocode",
    repair_symbols=2,
    number_of_splits=1,
):
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    dist = RaptorDistribution(number_of_chunks)
    dna_rules = FastDNARules()
    if asdna:
        rules = dna_rules
    else:
        rules = None
    x = RU10Encoder(
        file,
        number_of_chunks,
        dist,
        chunk_size=chunk_size,
        insert_header=insert_header,
        rules=rules,
        error_correction=error_correction,
        id_len_format=seed_len_format,
        number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
        save_number_of_chunks_in_packet=save_number_of_chunks_in_packet,
        prepend=prepend,
        append=append,
    )
    x.prepare()
    i = 0
    tmp_list = []
    while i < while_count:
        if seq_seed is not None:
            if seq_seed + i >= NUMBER_OF_PACKETS_TO_CREATE:
                break
            packet = x.create_new_packet(seed=seq_seed + i)
        else:
            packet = x.create_new_packet()
        should_drop_packet(rules, packet)
        assert packet.error_prob is not None
        if packet.error_prob <= drop_above and (
            len(tmp_list) < l_size or packet.error_prob < tmp_list[-1].error_prob
        ):
            if packet not in tmp_list:
                bisect.insort_left(tmp_list, packet)
            else:
                elem = next((x for x in tmp_list if x == packet), None)
                if elem is not None and packet < elem:
                    tmp_list.remove(elem)
                    bisect.insort_left(tmp_list, packet)
            if len(tmp_list) > l_size:
                tmp_list = tmp_list[:l_size]
        i += 1
    print([x.error_prob for x in tmp_list])
    conf = {
        "error_correction": e_correction_name,
        "repair_symbols": repair_symbols,
        "asdna": asdna,
        "number_of_splits": number_of_splits,
        "find_minimum_mode": True,
        "seq_seed": seq_seed,
    }
    x.save_config_file(conf, section_name="RU10_" + file)
    return [ParallelPacket.from_packet(p) for p in tmp_list]


def save_packets(packets, out_file, clear_output=True, seed_is_filename=True):
    if not out_file.endswith("/"):
        files = glob.glob(out_file + "/*")
    else:
        files = glob.glob(out_file + "*")
    if clear_output:
        for f in files:
            os.remove(f)
    i = 0
    e_prob = ""
    if not os.path.exists(out_file):
        os.makedirs(out_file)
    for packet in packets:  # sorted(packets, key=lambda elem: (elem.error_prob, elem.__hash__())):
        if seed_is_filename:
            i = packet.id
        e_prob = (str(ceil(packet.error_prob * 100)) + "_") if packet.error_prob is not None else ""
        with open(out_file + "/" + e_prob + str(i) + ".RU10_DNA", "w") as f:
            f.write(packet.get_dna_struct(True))
        i += 1


def save_packets_fasta(packets, out_file, file_ending, clear_output=True, seed_is_filename=True):
    if not out_file.endswith("/"):
        files = glob.glob(out_file + "/*")
    else:
        files = glob.glob(out_file + "*")
    if clear_output:
        for f in files:
            try:
                os.remove(f)
            except FileNotFoundError:
                print("Error while removing file: {}".format(f))
    i = 0
    e_prob = ""
    if not os.path.exists(out_file):
        os.makedirs(out_file)
    with open(
        out_file
        + "/"
        + str(random.randint(1, 80000000000))
        + "_"
        + str(random.randint(1, 800000000000))
        + ".fasta",
        "w",
    ) as f:
        for (
            packet
        ) in packets:  # sorted(packets, key=lambda elem: (elem.error_prob, elem.__hash__())):
            if seed_is_filename:
                i = packet.id
            e_prob = (
                (str(ceil(packet.error_prob * 100)) + "_") if packet.error_prob is not None else ""
            )
            f.write(">" + e_prob + str(i) + file_ending + "\n" + packet.get_dna_struct(True) + "\n")
            i += 1


def reduceLists(base, input_list=None, l_size=100):
    if input_list is None and len(base) == 2:
        input_list = base[1]
        base = base[0]
    base = base[:l_size]
    if input_list is None:
        return base
    for packet in input_list:
        if packet.error_prob < base[-1].error_prob:
            if packet not in base:
                bisect.insort_left(base, packet)
            else:
                elem = next((x for x in base if x == packet), None)
                if elem is not None and packet < elem:
                    # print("new packet is better")
                    base.remove(elem)
                    bisect.insort_left(base, packet)
            # print("Found duplicate!:" + str(packet))
            if len(base) > l_size:
                base = base[:l_size]
        else:
            # we can skip the rest of the current result-row since all remaining packets have a higher error_prob
            # than the worst packet of the current answer
            return base
    return base


# multiprocess and merge all created lists at the end
# @DeprecationWarning
def main(
    filename=".INFILES/logo.jpg",
    while_count=1000,
    l_size=1000,
    chunksize=0,
    number_of_chunks=300,
    sequential=False,
    spare1core=False,
    prepend="",
    append="",
    insert_header=False,
    seed_size_str=ID_LEN_FORMAT,
    drop_above=1.0,
    save_as_fasta=DEFAULT_SAVE_AS_FASTA,
    error_correction=nocode,
    e_correction_name="nocode",
    repair_symbols=2,
    number_of_splits=1,
):
    cores = multiprocessing.cpu_count()
    if spare1core:
        cores = cores - 1
    p = multiprocessing.Pool(cores)
    param = [None] * cores
    if sequential:
        stepsize = NUMBER_OF_PACKETS_TO_CREATE / cores
        param = [int(floor(i * stepsize)) for i in range(cores)]
        while_count = int(ceil(stepsize)) + 1
    a = p.map(
        partial(
            run,
            file=filename,
            l_size=l_size,
            while_count=while_count,
            chunk_size=chunksize,
            number_of_chunks=number_of_chunks,
            prepend=prepend,
            append=append,
            insert_header=insert_header,
            seed_len_format=seed_size_str,
            drop_above=drop_above,
            error_correction=error_correction,
            e_correction_name=e_correction_name,
            repair_symbols=repair_symbols,
            number_of_splits=number_of_splits,
        ),
        param,
    )

    while len(a) > 1:
        if len(a) % 2 != 0:
            a.append([])
        tmp = []
        for i in range(0, len(a) - 1, 2):
            tmp.append([a[i], a[i + 1]])
        a = p.map(partial(reduceLists, input_list=None, l_size=l_size), tmp)
    a = a[0]
    print("Parallel map/reduce:")
    print([x.error_prob for x in a])
    if save_as_fasta:
        save_packets_fasta(
            a,
            os.path.dirname(os.path.realpath(filename))
            + "/RU10_"
            + os.path.basename(filename)
            + "/",
            ".RU10_DNA",
        )
    else:
        save_packets(
            a,
            os.path.dirname(os.path.realpath(filename))
            + "/RU10_"
            + os.path.basename(filename)
            + "/",
        )
    return a


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument(
        "--repair_symbols",
        metavar="repair_symbols",
        required=False,
        type=int,
        default=2,
        help="number of repair symbols for ReedSolomon (default=2)",
    )
    parser.add_argument(
        "--list_size",
        metavar="list_size",
        required=False,
        type=int,
        default=1000,
        help="size of operational list per thread [inferred by #cores if sequential is set to true]",
    )
    parser.add_argument(
        "--out_size",
        metavar="out_size",
        required=False,
        type=int,
        default=1000,
        help="number of packets to save",
    )
    parser.add_argument(
        "--chunk_size",
        metavar="chunk_size",
        required=False,
        type=int,
        default=0,
        help="chunk size (default=0 [infer from number of chunks])",
    )
    parser.add_argument(
        "--number_of_chunks",
        metavar="number_of_chunks",
        required=False,
        type=int,
        default=300,
        help="number of chunks (default=300) [ignored if chunk_size is set to value != 0]",
    )
    parser.add_argument("--sequential", required=False, default=False, action="store_true")
    parser.add_argument("--spare1core", required=False, default=False, action="store_true")
    parser.add_argument(
        "--split_input",
        metavar="split_input",
        required=False,
        type=int,
        default=1,
        help="number of subcodings to split input file into",
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
    parser.add_argument("--insert_header", required=False, default=False, action="store_true")
    parser.add_argument(
        "--drop_above",
        required=False,
        type=float,
        default=1.0,
        help="Instantly drop packet if error-prob above set threshold",
    )
    parser.add_argument(
        "--seed_size_str",
        required=False,
        type=str,
        default="I",
        help=(
            "struct-string for seed - possible values: I,H,B (see struct). "
            "This impacts the amount of generated packets in sequential mode"
        ),
    )
    parser.add_argument(
        "--store_as_fasta",
        required=False,
        default=DEFAULT_SAVE_AS_FASTA,
        action="store_true",
        help="if set, store result as fasta file",
    )
    return parser


def _prepare_input_files(
    input_file: str, number_of_splits: int
) -> tuple[list[str], dict[str, str]]:
    if number_of_splits <= 1:
        return [input_file], {input_file: ""}

    input_files = split_file(input_file, number_of_splits)
    power_of_four = find_ceil_power_of_four(len(input_files))
    print(
        "Spltting input into {} sub files. We need to prepend/append {} base(s)".format(
            len(input_files), power_of_four
        )
    )
    prepend_matching = {
        input_files[i]: number_to_base_str(i, power_of_four) for i in range(len(input_files))
    }
    return input_files, prepend_matching


def _run_cli() -> None:
    args = _build_parser().parse_args()
    _input_file = args.filename
    _repair_symbols = args.repair_symbols
    _list_size = args.list_size
    _out_size = args.out_size
    _chunk_size = args.chunk_size
    _number_of_chunks = args.number_of_chunks
    _sequential = args.sequential
    _spare1core = args.spare1core
    _number_of_splits = args.split_input
    _insert_header = args.insert_header
    e_correction = args.error_correction
    _drop_above = args.drop_above
    seed_size_str = args.seed_size_str
    store_as_fasta = args.store_as_fasta
    _error_correction = get_error_correction_encode(e_correction, _repair_symbols)
    input_files, prepend_matching = _prepare_input_files(_input_file, _number_of_splits)
    for _file in input_files:
        print("File to encode: " + str(_file))
        main(
            _file,
            _list_size,
            _out_size,
            _chunk_size,
            _number_of_chunks,
            _sequential,
            _spare1core,
            append=prepend_matching[_file],
            insert_header=_insert_header,
            seed_size_str=seed_size_str,
            drop_above=_drop_above,
            save_as_fasta=store_as_fasta,
            error_correction=_error_correction,
            e_correction_name=e_correction,
            repair_symbols=_repair_symbols,
            number_of_splits=_number_of_splits,
        )
    if len(input_files) > 1:
        merge_folder_content(
            "split_" + os.path.basename(_input_file),
            _input_file + "_combined_split_output",
            append_folder_name=True,
            clear_dest_folder=True,
        )


if __name__ == "__main__":
    _run_cli()
