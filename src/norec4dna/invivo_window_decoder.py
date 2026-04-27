"""Sliding-window helpers for extracting RU10 packets from in-vivo FASTA input."""

import copy
import os
import random
import typing
from io import BytesIO
from pathlib import Path

from norec4dna import RU10Decoder
from norec4dna.ErrorCorrection import get_error_correction_decode
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte
from norec4dna.rules.FastDNARules import FastDNARules

REED_SOLOMON_PARITY_LENGTH = 3
_error_correction = get_error_correction_decode("reedsolomon", REED_SOLOMON_PARITY_LENGTH)
RULES_DROP_LIMIT = 2.0
RULES = FastDNARules()
NUM_CHUNKS = 70
PACKET_SEQ_LENGTH = 300
USE_HEADER_CHUNK = False
ID_LEN_FORMAT = "H"
CRC_LEN_FORMAT = "I"
# this file should include a SINGLE fasta-entry:
INPUT_FILE = "assembly.fasta"
# OPTIONAL - only required for creating a match-html file for whitebox analysis:
# this file should include all original sequences encoded into DNA (in any order)
GROUND_TRUTH = ".INFILES/merged_file.fasta"


def create_decoder() -> RU10Decoder:
    """
    create a fresh decoder instance
    """
    decoder = RU10Decoder(
        INPUT_FILE,
        use_headerchunk=USE_HEADER_CHUNK,
        error_correction=_error_correction,
        static_number_of_chunks=NUM_CHUNKS,
    )
    decoder.number_of_chunks = NUM_CHUNKS
    decoder.read_all_before_decode = True
    return decoder


def load_fasta_list(fasta_file: typing.Union[Path, str]) -> typing.List[typing.Tuple[str, str]]:
    """
    Load a FASTA file and return all entries as a list of ``(header, sequence)`` tuples.

    Unlike :func:`load_fasta`, this function preserves *every* entry in the file,
    including multiple entries that share the same header key.  The order of
    entries matches the order in the file.

    Args:
        fasta_file: Path to the FASTA file to read.

    Returns:
        A list of ``(header, sequence)`` tuples where *header* is the sequence
        name (the part after ``>`` up to the first whitespace) and *sequence*
        is the concatenated nucleotide string for that entry.

    Example:
        >>> entries = load_fasta_list("pool.fasta")
        >>> for header, seq in entries:
        ...     print(header, seq[:10])
    """
    entries: typing.List[typing.Tuple[str, str]] = []
    current_name = None
    current_seq: typing.List[str] = []
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    entries.append((current_name, "".join(current_seq)))
                current_name = line.split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name is not None:
        entries.append((current_name, "".join(current_seq)))
    return entries


def load_fasta(fasta_file: str) -> typing.Dict[str, str]:
    """
    Load a FASTA file and return a dictionary mapping header → sequence.

    .. deprecated::
        Use :func:`load_fasta_list` instead.  ``load_fasta`` silently drops all
        but the *last* entry when multiple entries share the same header key,
        which can cause data loss in pools that contain multiple packets with the
        same seed.  :func:`load_fasta_list` returns every entry as a
        ``(header, sequence)`` tuple list and is the preferred API going forward.

    Args:
        fasta_file: Path to the FASTA file to read.

    Returns:
        A ``dict`` mapping sequence name → sequence string.  When duplicate
        headers with differing content are present only the last entry for each
        key is kept and a ``UserWarning`` is emitted.

    Example:
        >>> sequences = load_fasta("pool.fasta")  # deprecated – prefer load_fasta_list
    """
    import warnings

    entries = load_fasta_list(fasta_file)

    # Build dict, detecting headers whose entries have differing content.
    fasta_dict: typing.Dict[str, str] = {}
    duplicates: typing.Set[str] = set()
    for name, seq in entries:
        if name in fasta_dict and fasta_dict[name] != seq:
            duplicates.add(name)
        fasta_dict[name] = seq  # keep the last entry for each key

    if duplicates:
        affected = sorted(duplicates)
        warnings.warn(
            f"load_fasta: {len(duplicates)} FASTA header(s) have multiple entries with "
            f"differing content – only the last entry per key is kept. "
            f"Affected headers: {affected[:10]}"
            f"{'...' if len(duplicates) > 10 else ''}. "
            "Use load_fasta_list() to retrieve all entries without data loss.",
            UserWarning,
            stacklevel=2,
        )

    return fasta_dict


def window_parse_packets():
    """
    This function uses a sliding window to find (correct) packets in the input sequnce
    """
    decoder = create_decoder()
    fasta = load_fasta(INPUT_FILE)
    window_start = 0
    packets = []
    correct = 0
    correct_seqs = []
    for k in fasta.keys():
        while (window_start + PACKET_SEQ_LENGTH) <= len(fasta[k]):
            line = fasta.get(k)[window_start : (window_start + PACKET_SEQ_LENGTH)]
            # ensure the packet adheres to the rules
            rule_err = RULES.apply_all_rules(line)
            if rule_err < RULES_DROP_LIMIT:
                new_pack = decoder.parse_raw_packet(
                    BytesIO(tranlate_quat_to_byte(line)).read(),
                    crc_len_format=CRC_LEN_FORMAT,
                    packet_len_format="",
                    number_of_chunks_len_format="",
                    id_len_format=ID_LEN_FORMAT,
                )
                # TODO: one could convert the RS-repaired data back to dna and check if it adheres to all rules - this
                # would be less restrictive since packets that violate constraints might still be repairable
                if new_pack is not None and new_pack != "CORRUPT":
                    # packet correct: add to correct_seqs and move the window by PACKET_SEQ_LENGTH
                    correct += 1
                    packets.append(new_pack)
                    correct_seqs.append(line)
                    window_start += PACKET_SEQ_LENGTH
                else:
                    # RS was unable to repair the packet: move the window by one base
                    window_start += 1
            else:
                # packet does not adhere to the rules: move the window by one base
                window_start += 1
    # (optional) write all correctly parsed sequences to a file
    # this will allow using the usual means of decoding the input (e.g. ConfigWorker)
    with open("correct_seqs.fasta", "w") as f:
        for i, seq in enumerate(correct_seqs):
            f.write(f">{i}\n{seq}\n")
    print(f"Correct sequences: {correct}")
    return packets


def main():
    packets = window_parse_packets()
    for i in range(1000):
        # start with a clean decoder
        _decoder = create_decoder()
        try:
            # shuffle the packets to reduce impact of single wrong packet
            random.shuffle(packets)
            for index, pack in enumerate(packets):
                if index != i:
                    _decoder.input_new_packet(pack)
            # if the decoder solved the matrix, the result should be correct
            # optionally, one could sanity check the result and reject the result if it is not correct
            # such a sanity check could be a check-sum at the ent of the file, or a check of the content
            # (e.g. only ASCII characters...) same for the filename stored in the header chunk (if available).
            if _decoder.solve(partial=True):
                _decoder.saveDecodedFile(
                    null_is_terminator=True,
                    print_to_output=True,
                    return_file_name=True,
                    partial_decoding=True,
                )
                print(f"Finished after {i} runs.")
                return
        except Exception:
            # The decoder will raise various exceptions if the decoding fails or yields an invalid result
            pass


def merge_files(folder: str):
    # write content of each file in the folder into a single file
    with open(folder + "/merged_file.fasta", "w") as outfile:
        for filename in os.listdir(folder):
            if filename.endswith(".txt"):
                with open(folder + "/" + filename, "r") as infile:
                    outfile.write(f"\n>{filename}\n")
                    for line in infile:
                        outfile.write(line)


def create_match_html():
    correct_seqs = 0
    ground_truth = load_fasta(GROUND_TRUTH)
    ground_truth = set(x for x in ground_truth.values())
    in_vivoed = load_fasta(INPUT_FILE)
    in_vivoed_seq = list(in_vivoed.values())[0]
    vivo_html = copy.deepcopy(in_vivoed_seq)
    for k in copy.deepcopy(ground_truth):
        if k in in_vivoed_seq:
            correct_seqs += 1
            vivo_html = vivo_html.replace(k, f"<br/><font color='green'>{k}</font><br/>")
            ground_truth.remove(k)
    print("Correct sequences:", correct_seqs)
    with open("match.html", "w") as f:
        f.write('<font color="red">' + vivo_html.replace("<br/><br/>", "<br/>") + "</font>")
    print(ground_truth)
    return "match.html"


if __name__ == "__main__":
    # main()
    create_match_html()
