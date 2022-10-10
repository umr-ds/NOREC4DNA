#!/usr/bin/python
import argparse
from norec4dna import Encoder
from norec4dna.LTEncoder import LTEncoder as LTEncoder
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode
from norec4dna.rules.DNARules_ErlichZielinski import DNARules_ErlichZielinski
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution
from norec4dna.distributions.IdealSolitonDistribution import IdealSolitonDistribution
from norec4dna.distributions.RobustSolitonDistribution import RobustSolitonDistribution

INSERT_HEADER = True
IMPLICIT_MODE = True
NUMBER_OF_CHUNKS_IN_PACKET = False
CHUNK_SIZE = 30


class demo_encode:
    @staticmethod
    def encode(file, error_correction=nocode, insert_header=INSERT_HEADER,
               save_number_of_chunks=NUMBER_OF_CHUNKS_IN_PACKET, save_as_fasta=True, save_as_zip=True, overhead=5.0,
               upper_bound=1.0, checksum_len_str=None):
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size=CHUNK_SIZE,
                                                                                 insert_header=insert_header)
        print("Number of Chunks=%s" % number_of_chunks)
        #dist = ErlichZielinskiRobustSolitonDistribution(number_of_chunks, seed=2)
        #dist = IdealSolitonDistribution(number_of_chunks, seed=2)
        dist = RobustSolitonDistribution(number_of_chunks, seed=2)
        encoder = LTEncoder(file, number_of_chunks, dist, insert_header=insert_header, rules=DNARules_ErlichZielinski(),
                            error_correction=error_correction, number_of_chunks_len_format="H", id_len_format="H",
                            used_packets_len_format="H", save_number_of_chunks_in_packet=save_number_of_chunks,
                            implicit_mode=IMPLICIT_MODE, drop_upper_bound=upper_bound,
                            checksum_len_str=checksum_len_str)

        encoder.set_overhead_limit(overhead)
        encoder.encode_to_packets()
        if save_as_fasta:
            encoder.save_packets_fasta(file_ending="_LT", seed_is_filename=True)
        elif save_as_zip:
            encoder.save_packets_zip(save_as_dna=True, file_ending="_LT", seed_is_filename=True)
        else:
            encoder.save_packets(True, save_as_dna=True, seed_is_filename=True, clear_output=True)

        encoder.save_packets(split_to_multiple_files=True, save_as_dna=True)
        print("Number of Chunks=%s" % encoder.number_of_chunks)
        return encoder


if __name__ == "__main__":
    # try:
    parser = argparse.ArgumentParser()
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
    parser.add_argument("--header_crc_str", metavar="header_crc_str", required=False, type=str, default="")
    parser.add_argument("--drop_upper_bound", metavar="drop_upper_bound", required=False, type=float, default=0.5,
                        help="upper bound for calculated error probability of packet before dropping")
    parser.add_argument("--overhead", metavar="overhead", required=False, type=float, default=0.40,
                        help="desired overhead of packets")

    args = parser.parse_args()
    filename = args.filename
    _insert_header = args.insert_header
    _repair_symbols = args.repair_symbols
    _save_number_of_chunks = args.save_number_of_chunks
    _upper_bound = args.drop_upper_bound
    _error_correction = get_error_correction_encode(args.error_correction, _repair_symbols)
    _save_as_fasta = args.save_as_fasta
    _save_as_zip = args.save_as_zip
    _overhead = args.overhead
    _header_crc_str = args.header_crc_str
    print("File to encode: " + str(filename))
    demo = demo_encode()
    encoder_instance = demo.encode(filename, error_correction=_error_correction,
                                   insert_header=_insert_header, save_number_of_chunks=_save_number_of_chunks,
                                   save_as_fasta=_save_as_fasta, save_as_zip=_save_as_zip, overhead=_overhead,
                                   upper_bound=_upper_bound, checksum_len_str=_header_crc_str)
    conf = {'error_correction': args.error_correction, 'repair_symbols': _repair_symbols, 'asdna': True,
            'number_of_splits': 0}
    config_filename = encoder_instance.save_config_file(conf, section_name="LT" + filename)
    print("Saved config file: %s" % config_filename)
    # input("Press Enter to continue ...")
