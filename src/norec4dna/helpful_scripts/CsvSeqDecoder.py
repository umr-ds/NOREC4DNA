import importlib
import importlib.util


def main() -> None:
    if importlib.util.find_spec("pandas") is None:
        raise ImportError("pandas is required to run CsvSeqDecoder")
    pandas = importlib.import_module("pandas")
    csv = pandas.read_csv(
        "/home/michael/Code/norec4dna/unilogo/results/finalData/filtered_table.csv"
    )
    is_bigger_th = csv
    print(is_bigger_th)
    res = is_bigger_th[is_bigger_th["sequences"].str.len() == 164]
    print(res["sequences"])
    for index, sequence in enumerate(res["sequences"], start=1):
        with open("out/" + str(index) + ".RU10_DNA", "w") as f:
            f.write(sequence.upper())


if __name__ == "__main__":
    main()
