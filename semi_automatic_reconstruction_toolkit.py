"""
This tool should allow a user to:
1) Decoder a file encoded with NOREC4DNA
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
    - the tool can only try to brutefoce the corrupt chunks / packets:
        IF WE BRUTEFORCE THE CHUNK WE MIGHT HAVE A PROBLEM IF THE PACKET HAD A MUTATION AT THE START (wrong ID!)
            we can avoid this pitfall by NOT using the chunk-mapping of the corrupt packet!
        IF WE BRUTEFORCE THE PACKET WE CANT DIRECTLY USE THE CRC (we must always perform a belief propagation / gauss elimination) - this is slower
"""
import sys

from norec4dna import Decoder


class SemiAutomaticReconstructionToolkit:
    def __init__(self, input_file, output_file):
        # TODO
        self.decoder: Decoder = None
        self.input_file = input_file
        self.output_file = output_file

    def decode(self):
        # TODO
        pass

    def view_file_with_chunkborders(self):
        """
        shows the content of decoder.b with borders after every n-th symbol
        """

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    SemiAutomaticReconstructionToolkit(input_file, output_file).decode()