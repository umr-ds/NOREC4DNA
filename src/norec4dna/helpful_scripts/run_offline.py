import os

COMMAND = (
    "python3 -m norec4dna.optimizer.optimize_dist_de "
    "Dorn --pop_size=10 --pop_size=15 --pop_size=20 --pop_size=25 "
    "--max_gen=250 --max_gen=200 --max_gen=150 --max_gen=100 "
    "--cr=0.4 --cr=0.4 --cr=0.4 --cr=0.4 "
    "--f=0.5 --f=0.5 --f=0.5 --f=0.5 --log=True --cores=4"
)

if __name__ == "__main__":
    os.system("nohup sh -c '" + COMMAND + " >res.txt' &")
