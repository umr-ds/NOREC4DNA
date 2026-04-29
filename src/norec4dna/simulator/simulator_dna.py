#!/usr/bin/python
# -*- coding: latin-1 -*-
import argparse
import importlib
import importlib.util
import math
import os
import time
from random import random

import colorama

from ..distributions.IdealSolitonDistribution import IdealSolitonDistribution
from ..distributions.OnlineDistribution import OnlineDistribution
from ..distributions.RaptorDistribution import RaptorDistribution
from ..distributions.RobustSolitonDistribution import RobustSolitonDistribution
from ..Encoder import Encoder
from ..LTBPDecoder import LTBPDecoder
from ..LTDecoder import LTDecoder
from ..LTEncoder import LTEncoder
from ..OnlineBPDecoder import OnlineBPDecoder
from ..OnlineDecoder import OnlineDecoder
from ..OnlineEncoder import OnlineEncoder
from ..RU10Decoder import RU10Decoder
from ..RU10Encoder import RU10Encoder
from ..rules.DNARules import DNARules
from . import bcolors

if os.name == "nt" and "PYCHARM_HOSTED" not in os.environ:
    colorama.init()
lines = [
    "Algorithm,ID,A_Permutation,T_Permutation,C_Permutation,G_Permutation,dinucleotid_Runs,Homopolymers,GC_Content,"
    "Trinucleotid_Runs,Random_Permutation,Overall_Dropchance,Random_Number,Did_Drop"
]

useDNARules = True
algo_type = []
packetToDropChance = {}


def _handle_bp_packet(decoder, packet, dec_input, invalid_drop, scale_first, scale_second):
    if should_drop_packet(packet, scale=scale_first):
        return None, dec_input, invalid_drop
    if should_drop_packet(packet, False, scale=scale_second):
        return None, dec_input, invalid_drop + 1
    if decoder.input_new_packet(packet):
        dec_input += 1
        decoder.solve()
        return decoder.is_decoded(), dec_input, invalid_drop
    return None, dec_input, invalid_drop


def _handle_fast_packet(decoder, packet, dec_input, invalid_drop, scale_first, scale_second):
    if should_drop_packet(packet, scale=scale_first):
        return None, dec_input, invalid_drop + 1
    if should_drop_packet(packet, False, scale=scale_second):
        return None, dec_input, invalid_drop
    dec_input += 1
    decoder.input_new_packet(packet)
    if packet.total_number_of_chunks <= dec_input - invalid_drop and decoder.solve():
        return decoder.is_decoded(), dec_input, invalid_drop
    return None, dec_input, invalid_drop


def _feed_replacement_packets(encoder, decoder, invalid_drop, dec_input, scale_first, scale_second):
    inserted_packets = 0
    while inserted_packets < invalid_drop:
        packet = encoder.create_and_add_new_packet()
        if should_drop_packet(packet, scale=scale_first):
            continue
        if not should_drop_packet(packet, False, scale=scale_second):
            decoder.input_new_packet(packet)
            dec_input += 1
        inserted_packets += 1
    return dec_input


def blackbox(encoder, decoder, scale_first=1.0, scale_second=1.0):
    encoder.encode_to_packets()
    encoded_packets = encoder.get_encoded_packets()
    dec_input = 0
    invalid_drop = 0
    print(bcolors.BOLD + "[+] Created " + str(len(encoded_packets)) + " Packets" + bcolors.ENDC)
    for packet in encoded_packets:
        if isinstance(decoder, LTBPDecoder) or isinstance(decoder, OnlineBPDecoder):
            result, dec_input, invalid_drop = _handle_bp_packet(
                decoder, packet, dec_input, invalid_drop, scale_first, scale_second
            )
        else:
            result, dec_input, invalid_drop = _handle_fast_packet(
                decoder, packet, dec_input, invalid_drop, scale_first, scale_second
            )
        if result is not None:
            print("DNA-Simulator dropped " + str(invalid_drop) + " Packets and finished successful")
            return result, dec_input, invalid_drop

    dec_input = _feed_replacement_packets(
        encoder, decoder, invalid_drop, dec_input, scale_first, scale_second
    )
    decoder.solve()
    print("[!] DNA-Simulator dropped " + str(invalid_drop) + " Packets.")
    return decoder.is_decoded(), dec_input, invalid_drop


def should_drop_packet(packet, add_line=True, scale=1.0, rules=DNARules):
    if not useDNARules:
        return False
    rand = random()  # create number from [0,1)
    if add_line:
        drop_chance, data, annotated_packet = rules.apply_all_rules_with_data(packet)
        drop_chance = drop_chance * scale  # scale according to parameter
        line = (
            algo_type[0]
            + ","
            + str(hex(packet.id))
            + ","
            + ",".join([str(round(x, 4)) for x in data])
            + ","
            + str(drop_chance)
            + ","
            + str(rand)
            + ","
            + str(drop_chance > rand)
        )
        lines.append(line)
        packetToDropChance[packet] = drop_chance
    elif packet in packetToDropChance:
        drop_chance = packetToDropChance[packet]
        packetToDropChance.clear()
    else:
        drop_chance, data, annotated_packet = rules.apply_all_rules_with_data(packet)
    # drop packet if rand bigger than the drop_chance for this Packet.
    return drop_chance > rand


def blackboxOnlineTest(
    file, number_of_chunks=800, seed=2, overhead=0.20, scale_first=1.0, scale_second=1.0
):
    start = time.time()
    epsilon = 0.024343  # pruefOrdnung: 0.007084 # according to filesize this should make a length of 200...
    quality = 5
    dist = OnlineDistribution(epsilon, seed)
    number_of_chunks = dist.get_size()
    assert number_of_chunks is not None
    algo_type.clear()
    algo_type.append("Online_" + str(number_of_chunks) + "_" + str(dist.get_config_string()))
    print(
        bcolors.OK
        + "Starting Blackbox Test with "
        + str(number_of_chunks)
        + " Chunks"
        + bcolors.ENDC
    )
    encoder = OnlineEncoder(file, number_of_chunks, dist, epsilon, quality)
    encoder.set_overhead_limit(overhead)
    decoder = OnlineDecoder.pseudo_decoder(number_of_chunks)
    decoder.set_read_all_before_decode(True)

    result, dec_input, invalid_drop = blackbox(
        encoder, decoder, scale_first=scale_first, scale_second=scale_second
    )
    end = time.time() - start
    print(
        bcolors.BLUE
        + "Blackbox-Decode "
        + (bcolors.OK + "successful" if result else bcolors.ERR + "NOT successful")
        + bcolors.END
        + bcolors.BLUE
        + " after "
        + str(round(end, 4))
        + " sec."
        + bcolors.ENDC
    )
    return [
        "Online_eps=" + str(epsilon) + "_quality=" + str(quality),
        result,
        number_of_chunks,
        dec_input,
        invalid_drop,
        round(end, 4),
    ]


def blackboxLTTest(
    file,
    number_of_chunks=800,
    seed=2,
    chunk_size=0,
    overhead=0.20,
    scale_first=1.0,
    scale_second=1.0,
):
    print(
        bcolors.OK
        + "Starting Blackbox Test with "
        + str(number_of_chunks)
        + " Chunks"
        + bcolors.ENDC
    )
    start = time.time()
    dist = RobustSolitonDistribution(S=number_of_chunks, delta=overhead, seed=seed)
    algo_type.clear()
    algo_type.append("LT_" + str(number_of_chunks) + "_" + str(dist.get_config_string()))
    encoder = LTEncoder(file, number_of_chunks, dist, chunk_size=chunk_size)
    encoder.set_overhead_limit(overhead)
    decoder = LTDecoder.pseudo_decoder(number_of_chunks)
    decoder.set_read_all_before_decode(True)

    result, dec_input, invalid_drop = blackbox(
        encoder, decoder, scale_first=scale_first, scale_second=scale_second
    )
    end = time.time() - start
    print(
        bcolors.BLUE
        + "Blackbox-Decode "
        + (bcolors.OK + "successful" if result else bcolors.ERR + "NOT successful")
        + bcolors.END
        + bcolors.BLUE
        + " after "
        + str(round(end, 4))
        + " sec."
        + bcolors.ENDC
    )
    return [
        dist.get_config_string(),
        result,
        number_of_chunks,
        dec_input,
        invalid_drop,
        round(end, 4),
    ]


def blackboxLTIdealTest(
    file,
    number_of_chunks=800,
    seed=2,
    chunk_size=0,
    overhead=0.20,
    scale_first=1.0,
    scale_second=1.0,
):
    print(
        bcolors.OK
        + "Starting Blackbox Test with "
        + str(number_of_chunks)
        + " Chunks"
        + bcolors.ENDC
    )
    start = time.time()
    dist = IdealSolitonDistribution(S=number_of_chunks, seed=seed)
    algo_type.clear()
    algo_type.append("LT_" + str(number_of_chunks) + "_" + str(dist.get_config_string()))
    encoder = LTEncoder(file, number_of_chunks, dist, chunk_size=chunk_size)
    encoder.set_overhead_limit(overhead)
    decoder = LTDecoder.pseudo_decoder(number_of_chunks)
    decoder.set_read_all_before_decode(True)
    result, dec_input, invalid_drop = blackbox(
        encoder, decoder, scale_first=scale_first, scale_second=scale_second
    )
    end = time.time() - start
    print(
        bcolors.BLUE
        + "Blackbox-Decode "
        + (bcolors.OK + "successful" if result else bcolors.ERR + "NOT successful")
        + bcolors.END
        + bcolors.BLUE
        + " after "
        + str(round(end, 4))
        + " sec."
        + bcolors.ENDC
    )
    return [
        dist.get_config_string(),
        result,
        number_of_chunks,
        dec_input,
        invalid_drop,
        round(end, 4),
    ]


def blackboxRU10Test(
    file,
    number_of_chunks=800,
    seed=2,
    chunk_size=200,
    overhead=0.20,
    scale_first=1.0,
    scale_second=1.0,
):
    print(
        bcolors.OK
        + "Starting Blackbox Test with "
        + str(number_of_chunks)
        + " Chunks"
        + bcolors.ENDC
    )
    start = time.time()
    dist = RaptorDistribution(number_of_chunks)
    algo_type.clear()
    algo_type.append("RU10_" + str(number_of_chunks) + "_" + str(dist.get_config_string()))
    encoder = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size)
    encoder.set_overhead_limit(overhead)
    decoder = RU10Decoder.pseudo_decoder(number_of_chunks)
    decoder.set_read_all_before_decode(True)

    result, dec_input, invalid_drop = blackbox(
        encoder, decoder, scale_first=scale_first, scale_second=scale_second
    )
    end = time.time() - start
    print(
        bcolors.BLUE
        + "Blackbox-Decode "
        + (bcolors.OK + "successful" if result else bcolors.ERR + "NOT successful")
        + bcolors.END
        + bcolors.BLUE
        + " after "
        + str(round(end, 4))
        + " sec."
        + bcolors.ENDC
    )
    return [
        dist.get_config_string(),
        result,
        number_of_chunks,
        dec_input,
        invalid_drop,
        round(end, 4),
    ]


def get_random_int(max_int):
    return int(random() * max_int)


def _get_number_of_chunks(file, chunk_size, insert_header):
    if chunk_size != 0:
        return Encoder.get_number_of_chunks_for_file_with_chunk_size(
            file, chunk_size=chunk_size, insert_header=insert_header
        )
    return 800


def _run_simulation_case(file, mode, number_of_chunks, rnd, overhead, scale_first, scale_second):
    if mode == "Online":
        return blackboxOnlineTest(
            file,
            number_of_chunks=number_of_chunks,
            seed=rnd,
            overhead=overhead,
            scale_first=scale_first,
            scale_second=scale_second,
        )
    if mode == "LT":
        return blackboxLTTest(
            file,
            number_of_chunks=number_of_chunks,
            seed=rnd,
            chunk_size=100,
            overhead=overhead,
            scale_first=scale_first,
            scale_second=scale_second,
        )
    if mode == "LTIdeal":
        return blackboxLTIdealTest(
            file,
            number_of_chunks=number_of_chunks,
            seed=rnd,
            chunk_size=100,
            overhead=overhead,
            scale_first=scale_first,
            scale_second=scale_second,
        )
    return blackboxRU10Test(
        file,
        number_of_chunks=number_of_chunks,
        seed=rnd,
        chunk_size=100,
        overhead=overhead,
        scale_first=scale_first,
        scale_second=scale_second,
    )


def _format_simulation_line(
    file, overhead, name, number_of_chunks, dec_input, invalid_drop, rnd, result, time_needed
):
    return (
        str(file)
        + ","
        + str(overhead)
        + ","
        + str(name)
        + ","
        + str(number_of_chunks)
        + ","
        + str(dec_input)
        + ","
        + str(invalid_drop)
        + ","
        + str(rnd)
        + ","
        + str(result)
        + ","
        + str(time_needed)
    )


def _write_simulation_outputs(mode, overhead, csv):
    dtimeno = (
        mode
        + "_"
        + str(overhead)
        + "_sim"
        + str(time.strftime("%Y-%m-%d_%H-%M", time.localtime()))
        + ".csv"
    )
    with open("DNA_" + dtimeno, "w") as f:
        for line in lines:
            f.write(line + "\n")
    lines.clear()
    lines.append(
        "Algorithm,ID,A_Permutation,T_Permutation,C_Permutation,G_Permutation,dinucleotid_Runs,Homopolymers,"
        "GC_Content,Trinucleotid_Runs,Random_Permutation,Overall_Dropchance,Random_Number,Did_Drop"
    )
    with open(dtimeno, "w") as f:
        for line in csv:
            f.write(line + "\n")


def main(file=".INFILES/logo.jpg", repeats=5):
    csv = [
        "filename, overhead, codecName, number_of_chunks, dec_input, invalid_drop, seed, result, time_needed"
    ]
    name = "ERROR"
    scale_first = 1.0
    scale_second = 0.5
    for mode in ["RU10", "LT", "LTIdeal", "Online"]:
        for overhead in [
            0.05,
            0.10,
            0.15,
            0.20,
            0.25,
            0.30,
            0.35,
            0.40,
            0.45,
            0.50,
            0.55,
            0.60,
        ]:
            for _ in range(repeats):
                rnd = get_random_int(math.pow(2, 31) - 1)
                number_of_chunks = _get_number_of_chunks(file, 100, True)
                name, result, number_of_chunks, dec_input, invalid_drop, time_needed = (
                    _run_simulation_case(
                        file, mode, number_of_chunks, rnd, overhead, scale_first, scale_second
                    )
                )
                line = _format_simulation_line(
                    file,
                    overhead,
                    name,
                    number_of_chunks,
                    dec_input,
                    invalid_drop,
                    rnd,
                    result,
                    time_needed,
                )
                print(line)
                csv.append(line)
            _write_simulation_outputs(mode, overhead, csv)
        csv = [
            "filename, codecName, number_of_chunks, dec_input, invalid_drop, seed, result, time_needed"
        ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze network traffic")
    parser.add_argument("-f", "--file", help="File to use for Simulation", required=True)
    parser.add_argument(
        "-p",
        "--profile",
        help='If set, a profiling Graph "pycallgraph.png" will be created.',
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "-r",
        "--repeats",
        type=int,
        help="Number of repeats the Simulator should run (default=5)",
        default=5,
    )
    args = parser.parse_args()
    profile = bool(args.profile)
    filename = str(args.file)
    repeats = int(args.repeats)
    if profile:
        if importlib.util.find_spec("pycallgraph") is None:
            raise ImportError("pycallgraph is required for profiling")
        pycallgraph = importlib.import_module("pycallgraph")
        pycallgraph_output = importlib.import_module("pycallgraph.output")

        print(
            bcolors.WARN
            + "[!] running with Profiler - this might decrease performance"
            + bcolors.ENDC
        )
        with pycallgraph.PyCallGraph(output=pycallgraph_output.GraphvizOutput()):
            main(filename, repeats)
        print(bcolors.BLUE + '[*] profiling Graph saved as "pycallgraph.png"' + bcolors.ENDC)
    else:
        main(filename, repeats)
