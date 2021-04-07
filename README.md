# NOREC4DNA [![Build Status](https://travis-ci.com/thejanky/fountaindna.svg?token=ChxYCfqxwUJpgKstks1n&branch=master)](https://travis-ci.com/thejanky/fountaindna) [![codecov](https://codecov.io/gh/thejanky/fountaindna/branch/master/graph/badge.svg?token=WglFq3uNTI)](https://codecov.io/gh/thejanky/fountaindna)

NOREC4DNA is an all-in-one Suite for analyzing, testing and converting Data into DNA-Chunks to use for a
DNA-Storage-System using integrated DNA-Rules as well as the MOSLA DNA-Simulation-API.

NOREC4DNA implements LT, Online, and Raptor (RU10) Fountain Codes.

## Overview:

- [Install](#Install)
    * [Docker](#Using-docker)
    * [Source](#From-source)
- [Usage](#Usage)
    * [Docker](#Docker)
- [Tools](#Tools)
- [TODO](#TODO)

---

## Install

### Using docker

+ Building the docker container from source:
    - ````git clone git@github.com:umr-ds/NOREC4DNA.git````
    - ```docker build . --tag norec4dna```
+ pulling the container from Dockerhub:
    - TBA once NOREC4DNA is available on Dockerhub
    
### From source

+ Clone the repository:
    - ````git clone git@github.com:umr-ds/NOREC4DNA.git````


+ OPTIONAL create a virtual environment (recommended):
    - ````python3 -m venv <name_of_virtualenv>````
    - activate/source the newly created venv


+ Installing the dependencies:
    - depending on your distro you might need to manually install `LLVM` as well as `gcc` and `build-essential`
    - ```pip3 install -r requirements.txt```


+ Install NOREC4DNA:
    - ```python3 setup.py install```

---

## Usage

### Docker
To get the en- and decoded files from NOREC4DNA using docker you might need to map a volume into the container.

Alternatively you could use the `docker cp` command to transfer the desired files.


## Find minimum

### Building and running

Build the docker container:

`docker build --tag norec4dna_gd`

Run:

`docker run --name norec4dna_gd_multiple_files -d -t -v /tmp/norec4dna/:/norec4dna/tmp norec4dna_gd (Parameter...)`

Alternatively you can run the script directly:
`python3 ru10_find_minimum_packets.py <Parameters>`

### Parameters

First enter the filename of the file to generate the packets for.

`FILE (--parameters)`

The following parameters can be set:

`--repair_symbols=[no_symbols]`

The number of repair_symbols for ReedSolomon (default=2). This does only apply if --error_correction is set to reedsolomon

`--list_size=[size]`

Size of operational list per thread, inferred by the number cores if sequential is set to true. The list size should
always be greater than the out_size to ensure optimal results (default=1000).

`--out_size=[size]`

Number of packets to save after combining the lists and sorting them by the packets error_prob (default=1000).

`--chunk_size=[size]`

Size of chunks to split the file into, inferred from number of chunks and the filesize if not set (default=0).

`--number_of_chunks=[no_chunks]`

Number of chunks to split the file into, ignored if chunksize is set to value != 0 (default=300).

`--sequential`

If set, all seed will be generated in a sequential matter. (Recommended!)

`--spare1core`

If activated, one core is not used for the calculation of the lists.

`--method=[RU10/Online/LT]`

Sets the method to generate the packets with. Available are RU10, Online and LT.

`--seed_size_str=[I,H,...]`

Set the struct-string for the seed field. See [https://docs.python.org/3/library/struct.html#format-characters](https://docs.python.org/3/library/struct.html#format-characters) for more information

`--drop_above`

Sets an upper-limit for the error probability. WARNING: This might reduce the total number of sequences returned!

## With optimization

`--optimization`

Activates the automated optimization of the chunk distribution in the packets with different options.

`--overhead=[overhead (0.1=10%)]`

Overhead to use for the optimization, where 0.1 means 10% additional packets based on the number of packets needed to
decode the file (default=0.1).

`--overhead_factor=[factor (0.1=10%)]`

If the overhead is not enough to optimize the packets, the overhead factor is a factor that allows exceeding the given
overhead to try to optimize the chunk distribution (default=0.0).

`--errorprob_factor=[factor (0.1=10%)]`

A factor for the maximum allowed error_prob of the additional packets based on the average error_prob of the packets
needed to decode (default=0.1). If set to 0.0 no more packets may be added since the packets with the lowest error_probs
were already used to decode the file.

`--plot`
Generates and saves different plots to show the results.

---

## Tools

### demo_*.py

TODO

### ru10_find_minimum_packets.py
`--error_correction [nocode, crc, reedsolomon]`

Defines the error detection / correction algorithm to use per packet. (Default: nocode = no error-detection/correction)

`--split_input`

Sets the number of pre-splits to perform Default: 1 (= do not split the input file into multiple NOREC rounds)
WARNING: If set, this value _should_ be known during decoding (thus using a bruteforce approach this value might be reconstructed) 

`--store_as_fasta`

If set, stores the result in a .fasta file instead of one file per sequence

`--insert_header`

If set, besides the created chunks an additional header chunk will be added. This chunk stores the filename and the correct padding for the last chunk.
(Recommended!) WARNING: If not set, the reconstructed file will most likely be longer due to the \00-padding at the end.
 

### ConfigWorker.py
Allows easy en- and decoding used .ini files.
Since the supplied encoder can create such .ini files, this is especially useful for easy decoding.

### helpful_scripts

there are various more or less useful scripts inside `helpful_scripts/`

TODO


---

## TODO:

- generate _**all**_ possible packets for each Method, log the error_prob of all of them
- check error-resilience: on non-successful decode: how much can be read..
- check required overhead to decode
- check absolute overhead (header, additional packets,...)
- select a licence
- merge ru10_find_minimum into find_minimum!