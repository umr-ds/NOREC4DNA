#!/usr/bin/python
import argparse

from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode
from norec4dna.distributions.OnlineDistribution import OnlineDistribution


class demo_online_encode:
    @staticmethod
    # def encode(file, error_correction=nocode, asdna=True):
    def encode(file, asdna=True,  error_correction=nocode, insert_header=False, save_number_of_chunks_in_packet=False,
               save_as_fasta=True, save_as_zip=True, overhead=0.40, epsilon=0.068, quality=7, upper_bound=1.0):
        dist = OnlineDistribution(epsilon)
        number_of_chunks = dist.get_size()
        dna_rules = FastDNARules()
        if asdna:
            rules = dna_rules
        else:
            rules = None
        encoder = OnlineEncoder(
            file, number_of_chunks, dist, epsilon, quality, error_correction=error_correction, quality_len_format="B",
            insert_header=insert_header, check_block_number_len_format="H", number_of_chunks_len_format="H", rules=rules,
            save_number_of_chunks_in_packet=save_number_of_chunks_in_packet, drop_upper_bound=upper_bound)  # , pseudo_decoder=pseudo)
        encoder.set_overhead_limit(overhead)
        #encoder.encode_file(split_to_multiple_files=True, save_as_dna=asdna)
        encoder.encode_to_packets()
        if save_as_fasta:
            encoder.save_packets_fasta(file_ending="_Online", seed_is_filename=True)
        elif save_as_zip:
            encoder.save_packets_zip(save_as_dna=True, file_ending="_Online", seed_is_filename=True)
        else:
            encoder.save_packets(True, save_as_dna=True, seed_is_filename=True, clear_output=True)

        encoder.save_packets(split_to_multiple_files=True, save_as_dna=True)
        print("Number of Chunks=%s" % encoder.number_of_chunks)
        return encoder


if __name__ == "__main__":
    # try:
    parser = argparse.ArgumentParser()
    #parser.add_argument("--chunk_size", metavar="chunk_size", required=False, type=int, default=DEFAULT_CHUNK_SIZE,
    #                    help="size of chunks to split the file into")
    parser.add_argument("filename", metavar="file", type=str, help="the file to Encode")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", action="store_true", required=False, default=False)
    parser.add_argument("--save_number_of_chunks", metavar="save_number_of_chunks", required=False, type=bool,
                        default=False)
    parser.add_argument("--save_as_fasta", action="store_true", required=False)
    parser.add_argument("--save_as_zip", action="store_true", required=False)
    parser.add_argument("--drop_upper_bound", metavar="drop_upper_bound", required=False, type=float, default=0.5,
                        help="upper bound for calculated error probability of packet before dropping")
    parser.add_argument("--overhead", metavar="overhead", required=False, type=float, default=0.40,
                        help="desired overhead of packets")
    parser.add_argument("--epsilon", metavar="epsilon", required=False, type=float, default=0.068,
                        help="epsilon value of the coding")
    parser.add_argument("--quality", metavar="quality", required=False, type=int, default=7,
                        help="quality value of the coding")
    args = parser.parse_args()
    _file = args.filename
    #_chunk_size = args.chunk_size
    _no_repair_symbols = args.repair_symbols
    _insert_header = args.insert_header
    _save_number_of_chunks = args.save_number_of_chunks
    _upper_bound = args.drop_upper_bound
    _error_correction = get_error_correction_encode(args.error_correction, _no_repair_symbols)
    _save_as_fasta = args.save_as_fasta
    _save_as_zip = args.save_as_zip
    _overhead = args.overhead
    _epsilon = args.epsilon
    _quality = args.quality
    # for _file in input_files:
    print("File to encode: " + str(_file))
    demo = demo_online_encode()
    encoder_instance = demo.encode(_file, error_correction=_error_correction,
                                   insert_header=_insert_header, save_number_of_chunks_in_packet=_save_number_of_chunks,
                                   save_as_fasta=_save_as_fasta, save_as_zip=_save_as_zip, overhead=_overhead,
                                   epsilon=_epsilon, quality=_quality, upper_bound=_upper_bound)
    conf = {'error_correction': args.error_correction, 'repair_symbols': _no_repair_symbols, 'asdna': True,
            'number_of_splits': 0, 'quality': encoder_instance.quality, 'epsilon': encoder_instance.epsilon,
            'find_minimum_mode': False, 'seq_seed': False}
    config_filename = encoder_instance.save_config_file(conf, section_name="Online_" + _file)
    print("Saved config file: %s" % config_filename)
    # input("Press Enter to continue ...")
