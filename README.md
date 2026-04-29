# NOREC4DNA

NOREC4DNA is a fountain-code-based DNA storage toolkit with LT, Online, and
Raptor (RU10) implementations plus supporting analysis, simulation, and
decoding utilities.

## Overview

- [Install](#install)
  - [Docker](#docker)
  - [From source](#from-source)
- [Development workflow](#development-workflow)
- [Compatibility notes](#compatibility-notes)
- [Usage](#usage)
  - [Find minimum](#find-minimum)
- [Tools](#tools)
- [Example](#example)

---

## Install

### Docker

:warning: Docker builds are currently not recommended.

- Build locally:
  - `git clone git@github.com:umr-ds/NOREC4DNA.git`
  - `docker build . --tag norec4dna`
- Docker Hub image:
  - not published yet

### From source

1. Clone the repository:

   `git clone git@github.com:umr-ds/NOREC4DNA.git`

2. Optionally create and activate a virtual environment:

   - `python3 -m venv .venv`
   - `source .venv/bin/activate`

3. Install build/runtime prerequisites if your platform requires them:

   - some environments need `gcc`, `build-essential`, `llvm`, and
     Python development headers

4. Install the currently required Python dependencies:

   - `python -m pip install --upgrade pip`
   - `python -m pip install -r requirements.txt`

5. Install NOREC4DNA in editable mode:

   - `python -m pip install -e .`

6. Optional: build a wheel/sdist explicitly:

   - `python -m build`

### Project layout

The package now uses a standard `src/` layout. The Python sources live in
`src/norec4dna/`, and packaging metadata is defined in `pyproject.toml`.

**If you plan to build NOREC4DNA from source under Windows, Anaconda is still
recommended.**

---

## Development workflow

Useful commands for local development:

- install dependencies: `python -m pip install -r requirements.txt`
- install editable package: `python -m pip install -e .`
- run tests: `python -m pytest`
- run pre-commit hooks: `pre-commit run --all-files`

The repository includes a maintained `.pre-commit-config.yaml`, so running the
hooks locally should match the current CI-style validation workflow.

---

## Compatibility notes

Some modules that historically lived in NOREC4DNA are now compatibility layers
in this checkout:

- `norec4dna.metadata_coding` forwards to
  `norec4dna_multiversion.metadata_coding`
- `norec4dna.semi_automatic_reconstruction_toolkit.SemiAutomaticReconstructionToolkit`
  is deprecated and forwards to
  `norec4dna_multiversion.reconstruction.SemiAutomaticReconstructionToolkit`
- `norec4dna.file_update_coding` is deprecated; use
  `norec4dna_multiversion.coder` instead

If you need metadata embedding or multi-version coding, use the DR4DNA
multiversion package directly instead of the legacy NOREC4DNA entry points.

---

## Usage

### Find minimum

#### Building and running

Build the docker container:

`docker build --tag norec4dna_gd`

Run:

`docker run --name norec4dna_gd_multiple_files -d -t -v /tmp/norec4dna/:/norec4dna/tmp norec4dna_gd (Parameter...)`

Alternatively, run the script directly:

`python -m norec4dna.find_minimum_packets <Parameters>`

#### Parameters

First enter the filename of the file to generate the packets for.

`FILE (--parameters)`

The following parameters can be set:

`--repair_symbols=[no_symbols]`

The number of repair symbols for Reed-Solomon (default: `2`). This only applies
if `--error_correction=reedsolomon` is set.

`--list_size=[size]`

Size of the operational list per thread, inferred by the number of cores if
sequential mode is enabled. The list size should always be greater than the
output size to ensure optimal results (default: `1000`).

`--out_size=[size]`

Number of packets to save after combining the lists and sorting them by packet
error probability (default: `1000`).

`--chunk_size=[size]`

Size of chunks to split the file into, inferred from number of chunks and file
size if not set (default: `0`).

`--number_of_chunks=[no_chunks]`

Number of chunks to split the file into; ignored if `--chunk_size` is set to a
non-zero value (default: `300`).

`--sequential`

If set, all seeds are generated sequentially (recommended).

`--spare1core`

If activated, one CPU core is not used for list generation.

`--method=[RU10/Online/LT]`

Sets the encoding method. Available values are `RU10`, `Online`, and `LT`.

`--seed_size_str=[I,H,...]`

Set the `struct` format string for the seed field. See the Python `struct`
format documentation for details.

`--drop_above`

Sets an upper limit for the error probability. Warning: this may reduce the
total number of sequences returned.

### With optimization

`--optimization`

Activates automated optimization of the chunk distribution in the packets.

`--overhead=[overhead (0.1=10%)]`

Overhead to use for the optimization, where `0.1` means 10% additional packets
based on the number needed to decode the file (default: `0.1`).

`--overhead_factor=[factor (0.1=10%)]`

If the overhead is not enough to optimize the packets, this factor allows
exceeding the given overhead to keep searching (default: `0.0`).

`--errorprob_factor=[factor (0.1=10%)]`

A factor for the maximum allowed error probability of additional packets based
on the average packet error probability needed to decode (default: `0.1`).

`--plot`

Generates and saves plots showing the results.

---

## Tools

### `demo_*.py`

Demo applications for fast encoding and decoding of sequences.

### `ru10_find_minimum_packets.py` (deprecated)

`--error_correction [nocode, crc, reedsolomon]`

Defines the error detection/correction algorithm to use per packet (default:
`nocode`).

`--split_input`

Sets the number of pre-splits to perform (default: `1`, meaning no split into
multiple NOREC rounds).

`--store_as_fasta`

If set, stores the result in a `.fasta` file instead of one file per sequence.

`--insert_header`

If set, an additional header chunk storing the filename and the correct padding
for the last chunk is added. Recommended.

### `ConfigWorker.py`

Allows easy encoding and decoding via `.ini` files. Since the supplied encoder
can create such `.ini` files, this is especially useful for repeatable decode
workflows.

### `helpful_scripts`

There are additional helper scripts in `helpful_scripts/`.

---

## Example

To try out NOREC4DNA you can use the `demo_*_encode.py` scripts:

`python -m norec4dna.demo_raptor_encode .INFILES/Dorn --error_correction=reedsolomon --repair_symbols=3 --as_dna --insert_header`

This should create a new folder `RU10_Dorn` as well as a `Dorn_*.ini` file.

To decode the file from DNA, either use `demo_*_decode.py`:

`python -m norec4dna.demo_raptor_decode RU10_Dorn --use_header_chunk --error_correction=reedsolomon --repair_symbols=3 --number_of_chunks=145`

Or use the `ConfigWorker` module:

`python -m norec4dna.ConfigWorker <name of the .ini-file>`

The decoded file will be saved as `DEC_RU10_Dorn` if no header chunk was added;
otherwise it will be saved under the original filename.
