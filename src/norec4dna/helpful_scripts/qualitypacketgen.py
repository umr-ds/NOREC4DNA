import importlib
import importlib.util
from typing import Any

from numpy import mean

from ..distributions.RaptorDistribution import RaptorDistribution
from ..Encoder import Encoder
from ..ErrorCorrection import reed_solomon_encode
from ..helper import should_drop_packet
from ..RU10Encoder import RU10Encoder
from ..rules.FastDNARules import FastDNARules

plt: Any = (
    importlib.import_module("matplotlib.pyplot")
    if importlib.util.find_spec("matplotlib.pyplot") is not None
    else None
)


class QualityPacketGen:
    def __init__(self, encoder):
        encoder.prepare()
        i = 0
        tmp_list = []

        while i < 1000:
            packet = encoder.create_new_packet()
            should_drop_packet(rules, packet)
            tmp_list.append(packet)
            i += 1

        QualityAnalyzer(tmp_list)


class QualityAnalyzer:
    def __init__(self, tmp_list):
        tmp_dict = {}
        for packet in tmp_list:
            err_prob = packet.error_prob
            for chunk_no in packet.get_used_packets():
                if chunk_no not in tmp_dict:
                    tmp_dict[chunk_no] = []
                tmp_dict[chunk_no].append(err_prob)

        # for key, data_list in tmp_dict.items():
        num_list = []
        mean_list = []
        for key in tmp_dict:
            mean_list.append(len(tmp_dict[key]) / mean(tmp_dict[key]))
            num_list.append(len(tmp_dict[key]))

        index = []
        data = []
        for _, (key, val) in enumerate(sorted(tmp_dict.items())):
            index.append(key)
            data.append(val)

        fig, (ax) = plt.subplots(ncols=1)
        ax.boxplot(data)
        ax.set_xticklabels(index)
        plt.show()

        plt.plot(num_list)
        plt.plot(mean_list)
        plt.show()


if __name__ == "__main__":
    file = ".INFILES/Dorn"
    chunk_size = 100
    norepairsymbols = 6
    save_number_of_chunks_in_packet = False
    insert_header = False
    rules = FastDNARules()

    def error_correction(data: bytes) -> bytes:
        return reed_solomon_encode(data, norepairsymbols)

    number_of_chunks = 50
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)

    dist = RaptorDistribution(number_of_chunks)
    x = RU10Encoder(
        file,
        number_of_chunks,
        dist,
        chunk_size=chunk_size,
        insert_header=insert_header,
        rules=rules,
        error_correction=error_correction,
        id_len_format="H",
        number_of_chunks_len_format="B",
        save_number_of_chunks_in_packet=save_number_of_chunks_in_packet,
    )
    aa = QualityPacketGen(x)
