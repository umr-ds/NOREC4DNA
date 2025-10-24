import math
import multiprocessing
import struct
from functools import partial
from multiprocessing import freeze_support

import numpy as np

from norec4dna import RU10Encoder, RaptorDistribution, nocode, Encoder, reed_solomon_encode
from norec4dna.Packet import ParallelPacket
from norec4dna.helper import should_drop_packet, bin2Quaternary
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.rules.RuleParser import longestSequenceOfChar


def run(seq_seed=None, file='logo.jpg', repair_symbols=2, insert_header=False,
        error_correction=nocode, save_number_of_chunks_in_packet=False, l_size=1000, while_count=1000,
        chunk_size=0, number_of_chunks=300, prepend="", append="", seed_len_format="I",
        number_of_chunks_len_format="I", packets_to_create=None, xor_by_seed=False, id_spacing=0):
    # global counter
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    dna_rules = FastDNARules()
    if packets_to_create is None:
        packets_to_create = math.pow(2, 8 * struct.calcsize(seed_len_format))
    rules = dna_rules
    if repair_symbols != 0 and error_correction != nocode:
        error_correction = lambda x: reed_solomon_encode(x, repair_symbols)
    dist = RaptorDistribution(number_of_chunks)
    x = RU10Encoder(file, number_of_chunks, dist, chunk_size=chunk_size, insert_header=insert_header, rules=rules,
                    error_correction=error_correction, id_len_format=seed_len_format,
                    number_of_chunks_len_format=number_of_chunks_len_format,
                    save_number_of_chunks_in_packet=save_number_of_chunks_in_packet, mode_1_bmp=False,
                    prepend=prepend, append=append, xor_by_seed=xor_by_seed, id_spacing=id_spacing)
    x.prepare()
    i = 0
    tmp_list = []
    while i < while_count:
        if seq_seed is not None:
            if seq_seed + i >= packets_to_create:
                break
            packet = x.create_new_packet(seed=seq_seed + i)
        else:
            packet = x.create_new_packet()
        _ = should_drop_packet(rules, packet)
        tmp_list.append(ParallelPacket.from_packet(packet))
        i += 1
    return tmp_list


def analyse_spacing(file, num_chunks, min_spacing, max_spacing, xor_by_seed, seed_len_format="H"):
    for i in range(min_spacing, max_spacing):
        print("Spacing: " + str(i))
        packets_to_create = math.pow(2, 8 * struct.calcsize(seed_len_format))
        stepsize = packets_to_create / cores
        param = [int(math.floor(i * stepsize)) for i in range(cores)]
        while_count = int(math.ceil(stepsize)) + 1
        a = p.imap(
            partial(run, file=file, repair_symbols=0, number_of_chunks=num_chunks,
                    insert_header=True, seed_len_format=seed_len_format, xor_by_seed=xor_by_seed,
                    id_spacing=i, while_count=while_count), param)
        tmp = [y for x in a for y in x]
        print(f"Mean: {np.mean([x.error_prob for x in tmp])}")
        print(f"Std: {np.std([x.error_prob for x in tmp])}")
        print(f"# Packets < 1.0: {len([x for x in tmp if x.error_prob < 1.0])}")

def calculate_impossible_seeds(seed_len_format="H"):
    bad = 0
    total_packets = int(math.pow(2, 8 * struct.calcsize(seed_len_format)))
    for i in range(0, total_packets):
        raw = struct.pack(seed_len_format, i)
        dna_seq = "".join(bin2Quaternary.string2QUATS(raw))
        res = longestSequenceOfChar(dna_seq, '*')[1]
        if res > 3:
            bad+=1
    print(f"Seed format: {seed_len_format}, Bad: {bad}, Total: {total_packets}, Percentage: {1.0*bad/total_packets}")

def create_pandas(filename):
    import pandas as pd

    # Read the file
    with open(filename, 'r') as file:
        lines = file.readlines()[4:]

    # Initialize empty lists to store the data
    spacings = []
    means = []
    stds = []
    packets = []

    # Process the lines and extract the data
    for i in range(0, len(lines), 4):
        spacing = int(lines[i].split(': ')[1])
        mean = float(lines[i + 1].split(': ')[1])
        std = float(lines[i + 2].split(': ')[1])
        packet = int(lines[i + 3].split(': ')[1])

        spacings.append(spacing)
        means.append(mean)
        stds.append(std)
        packets.append(packet)

    # Create a DataFrame
    data = {
        'Spacing': spacings,
        'Mean': means,
        'Std': stds,
        'Packets < 1.0': packets
    }
    return pd.DataFrame(data)

if __name__ == '__main__':
    freeze_support()
    calculate_impossible_seeds("H")
    calculate_impossible_seeds("I")
    calculate_impossible_seeds("Q")

    cores = multiprocessing.cpu_count() - 1
    p = multiprocessing.Pool(cores)
    analyse_spacing("Dorn", 288, 10, 20, False)
    analyse_spacing("Dorn", 288, 10, 20, True)

    create_pandas("spacing_results")
    create_pandas("spacing_results_w_xor")