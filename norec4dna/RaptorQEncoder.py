from norec4dna import Encoder
from raptorq import Encoder as RQEncoder

from norec4dna.helper.bin2Quaternary import string2QUATS


class RaptorQEncoder(Encoder):
    def __init__(self, file, number_of_chunks: int = 0, distribution=None):
        super().__init__(file, number_of_chunks, distribution)
        data = b""
        with open(file, "rb") as in_file:
            data = in_file.read()
        encoder = RQEncoder.with_defaults(data, 80)
        x = encoder.get_encoded_packets(200000)
        print(["".join(string2QUATS(i)) for i in x][-50:])


if __name__ == "__main__":
    RaptorQEncoder("logo.jpg")
