import argparse
import logging
import struct
import typing
from typing import Any, Generator

import numpy as np
from PIL import Image

from NOREC4DNA.norec4dna.HeaderChunk import HeaderChunk
from NOREC4DNA.norec4dna.Packet import Packet
from NOREC4DNA.norec4dna.helper.helper_cpu_single_core import should_drop_packet, xor_numpy
from NOREC4DNA.norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.helper.RU10Helper import int31, choose_packet_numbers, from_true_false_list
from NOREC4DNA.ConfigWorker import ConfigReadAndExecute
from NOREC4DNA.metadata_coding import encoder_from_decoder
from NOREC4DNA.norec4dna import RU10Encoder
from NOREC4DNA.norec4dna.RU10Packet import RU10Packet
from NOREC4DNA.norec4dna.helper.bin2Quaternary import quads2dna, byte2QUATS, string2QUATS
from NOREC4DNA.norec4dna.helper.quaternary2Bin import quats_to_bytes, tranlate_quat_to_byte
from repair_algorithms.utils.select_numbers import select_numbers
from semi_automatic_reconstruction_toolkit import SemiAutomaticReconstructionToolkit
import imagehash

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def find_affected_chunks(
        semiautomatic_solver: SemiAutomaticReconstructionToolkit, new_file: bytes
) -> typing.Tuple[np.ndarray, np.ndarray]:
    """
    Calculate diff between existing chunks and new file
    - for now requires the new file to have equal length or be shorter!
    - index is defined WITHOUT header row as the header must be changed for any changed content!
    """
    semiautomatic_solver.get_file_as_bytes()

    # for now use this as a new_file:
    new_file = bytearray(semiautomatic_solver.get_file_as_bytes())
    new_file[500] = 0x00
    new_file[200] = 0x00
    assert (semiautomatic_solver.decoder.GEPP.b.shape[0] - 1) * semiautomatic_solver.decoder.GEPP.b.shape[1] == len(
        new_file)
    # create compatible numpy array with correct chunk-sizes from bytes:
    # TODO: allow growing by defining a binary diff and referencing a new "location" (e.g. new fountain coded file containing the diff) -> this can be stored in the new header version of the
    new_array = np.frombuffer(new_file, dtype=np.uint8).reshape(
        -1, semiautomatic_solver.decoder.GEPP.b.shape[1]
    )  # .copy()
    # new_array.flags.writeable = True
    # new_array[5,4] = 0
    diff = semiautomatic_solver.decoder.GEPP.b[1:] - new_array

    differing_rows = np.where(
        np.any(diff != 0, axis=1)
    )  # for correct mapping with a header: + np.array(1)

    return diff, differing_rows[0]


"""
def generate_packets(semiautomatic_solver: SemiAutomaticReconstructionToolkit, ru10_encoder: RU10Encoder,
                     number_of_packets_to_generate: int, changed_header: typing.Optional[HeaderChunk] = None) -> \
        typing.List[Packet]:
    # generate an encoder with the exact values used for decoding - including the header chunk!
    # if changed_header is None, we must avoid using the header chunk in ANY new packet!
    # for i in range(number_of_packets_to_generate):
    #    semiautomatic_solver.decoder.pack
    #    # pass

    seed = 0
    for seed in range(10000):
        packet_numbers = choose_packet_numbers(ru10_encoder.number_of_chunks, seed, ru10_encoder.dist, systematic=False)
        # if any([x in ])
        # break

    ru10_encoder.create_new_packet(systematic=False, seed=seed)
"""


def update_file(
        semiautomatic_solver: SemiAutomaticReconstructionToolkit,
        encoded_packets,
        new_file,
        packets_to_generate,
        updated_header,
):
    """
    Returns a list of encoded packets representing the changed content.
    The function calculates the changed chunks, the corresponding packets to update and how often these chunks should
    be updated in the final packets.

    @param semiautomatic_solver: SemiAutomaticReconstructionToolkit
    @param encoded_packets: list of encoded packets
    @param new_file: path to new file
    @param packets_to_generate: number of packets to generate
    @param updated_header: updated header2
    """

    pass


def create_perceptual_hash(image: Image.Image):
    hash = imagehash.phash(image)
    return str(hash)


"""
def main():
    # Encode:
    # Load the ini file
    # Decode all packets and find packets that use chunk 0 (header chunk)
    # For the last Calculate: Wanted DNA Sequence
    # Decode:
    # Load the ini file
    # decode using semi_automatic_reconstruction_toolkit
    # extract and parse all version (using inconsistent
    # save all versions + metadata into the filesystem
    pass

    header = semiautomatic_solver.decoder.headerChunk
    ru10_encoder = encoder_from_decoder(semiautomatic_solver, cfg_worker)
    #update_header(ru10_encoder, semiautomatic_solver, "", "", "")
    
    # generate_packets(semiautomatic_solver, ru10_encoder, 10, None)
"""


def generate_dna_version_string(fileversion: int = 0, base_length: int = 8) -> str:
    if fileversion >= 64:
        raise RuntimeError("fileversion must be in range [0, 63]!")
    return byte2QUATS(fileversion)[1:] + "GAGCCAGTGAGTCGTA"


def init_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ini",
        metavar="ini",
        type=str,
        help="config file (ini)",
        default="/home/michael/Code/DR4DNA/eval/sleeping_beauty_no_error.ini",
    )
    parser.add_argument(
        "--new_file",
        metavar="new_file",
        type=str,
        help="updated file",
    )
    parser.add_argument(
        "--packet_add_limit",
        metavar="packet_add_limit",
        type=int,
        default=5,
        help="maximum number of packets to add for each changed chunk",
    )
    return parser.parse_args()
    """
    parser.add_argument(
        "--ini",
        metavar="ini",
        type=str,
        help="config file (ini)",
        default="/home/michael/Code/DR4DNA/eval/sleeping_beauty_RU10.ini",
    )
    wanted_arg_group = parser.add_mutually_exclusive_group(required=True)
    wanted_arg_group.add_argument(
        "--wanted_metadata_file",
        metavar="wantedmetafile",
        type=str,
        help="file containing unwanted metadata in the fasta format, description should contain a multiplier and the corresponding metadata ",
    )
    wanted_arg_group.add_argument(
        "--wanted_metadata",
        metavar="wantedmeta",
        type=str,
        help="comma-separated list of metadata DNA sequences",
    )
    unwanted_arg_group = parser.add_mutually_exclusive_group(required=False)
    unwanted_arg_group.add_argument(
        "--unwanted_metadata_file",
        metavar="unwantedmetafile",
        type=str,
        help="file containing unwanted metadata in the fasta format, description should contain a multiplier and the corresponding metadata ",
    )
    unwanted_arg_group.add_argument(
        "--unwanted_metadata",
        metavar="wantedmeta",
        type=str,
        help="comma-separated list of metadata DNA sequences",
    )
    return parser.parse_args()
    """


def get_current_file_version(
        semi_automatic_solver: SemiAutomaticReconstructionToolkit,
) -> int:
    magic_string = "GAGCCAGTGAGTCGTA"
    from NOREC4DNA.invivo_window_decoder import load_fasta

    fasta_content = load_fasta(semi_automatic_solver.decoder.file)

    # normalize to an iterable of sequence strings:
    if hasattr(fasta_content, "values"):
        sequences = fasta_content.values()
    else:
        # assume list/iterable of tuples like (header, seq)
        try:
            sequences = (entry[1] for entry in fasta_content)
        except Exception:
            sequences = (str(entry) for entry in fasta_content)

    # check if magic_string is at the end of any sequence
    found_strings = set()
    for seq in sequences:
        if seq is None:
            continue
        if isinstance(seq, bytes):
            seq_str = seq.decode(errors="ignore")
        else:
            seq_str = str(seq)
        if seq_str.endswith(magic_string):
            found_strings.add(seq_str)
    highest_version = 0
    if len(found_strings) == 0:
        return highest_version

    for found_string in found_strings:
        found_string = found_string[len(found_string) - 3 - len(magic_string):-len(magic_string)]
        tmp = ord(quats_to_bytes(f"A{found_string}"))
        if tmp > highest_version:
            highest_version = tmp
    return highest_version


def insert_dna_version_string(packet, dna_version_string, insertion_position=30, diff=None):
    # TODO: this approach might break seed spacing if the spaced seed overlaps with the version string!
    # TODO: check seed-spacing for overlap!
    bkp_dna_data = packet.dna_data
    if packet.dna_data is None or packet.dna_data == "":
        packet.prepare_and_pack()
    assert not diff[insertion_position:insertion_position + int((len(dna_version_string) + 1) / 2) + 1].any(), \
        "Cannot insert version string into non-padding region!"
    packet.dna_data = packet.dna_data[:insertion_position * 2] + dna_version_string + packet.dna_data[
        insertion_position * 2 + len(dna_version_string):]
    deinterleaved_dna_data = Packet.deinterleave_spacing(packet.dna_data, packet.id_spacing, struct.calcsize(packet.id_len_format))
    bin_data = tranlate_quat_to_byte(deinterleaved_dna_data)
    tmp = packet.dna_data
    packet.dna_data = None
    dna_data = packet.get_dna_struct(False, packet.id_spacing, struct.calcsize(packet.id_len_format))
    np_packet = np.frombuffer(bin_data, dtype=np.uint8).copy()
    # FIXME: np_packet != packet.data for id (THIS IS TRUE FOR ANY PACKET DNA_DATA VS packed_data!!?!
    packet.packed = np_packet


# TO+DO: add the version_string (replace end of each generated packet with version string - use metadata_coding)
#  - convert version string to binary and replace end of potential_packet.data with this.
def reduce_packet_to_chunk(packet: RU10Packet,
                           semiautomatic_solver: SemiAutomaticReconstructionToolkit, chunk_to_reduce_to=0):
    # must be done such that the decoding is performed with all OLD chunks!
    # get normalized used chunks and remote chunk
    used_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(packet)
    packet.used_packets = set(i for i, x in enumerate(used_chunks) if x)
    for i, must_remove_chunk in enumerate(used_chunks):
        if must_remove_chunk and i != chunk_to_reduce_to:
            packet.data = xor_numpy(packet.data, semiautomatic_solver.decoder.GEPP.b[i])
            packet.used_packets.remove(i)
    assert len(packet.used_packets) == 1 and chunk_to_reduce_to in packet.used_packets
    return packet
    # packet contains only the header chunk
    pass
    # TODO:
    #  chunk_0 = forall id in used_chunks [except chunk 0] of new_packet: new_packet = new_packet XOR id
    #  for this, we must know the diff between the original changed chunk and the new version! (apply the diff!)


def insert_id_string(packet: RU10Packet,
                     insertion_position: int,
                     changed_chunk_id: int,
                     semiautomatic_solver: SemiAutomaticReconstructionToolkit):
    # TODO: for DECODING, we would have to make sure that the parse_header uses only the ORIGINAL packets (remove any sequence with the version flag in it!)
    # semiautomatic_solver.parse_header(checksum_len_format="H",
    #                                  last_chunk_len_format="I")  # ensure original header is present
    #
    # header = reduce_packet_to_chunk(packet, semiautomatic_solver)
    # header = HeaderChunk(header, last_chunk_len_format=semiautomatic_solver.last_chunk_len_format,
    #                     checksum_len_format=semiautomatic_solver.decoder.checksum_len_str)
    # TODO we replace the first n-bytes of the filename with the id of the changed chunk (n determined by the number of chunks len format):
    to_xor = np.zeros_like(packet.data)
    used_chunks = from_true_false_list(semiautomatic_solver.decoder.removeAndXorAuxPackets(packet))
    assert changed_chunk_id in used_chunks, "Changed chunk id must be in used chunks!"  # TODO: handle this case!
    index_of_chunk = used_chunks.index(changed_chunk_id)
    assert index_of_chunk < 256, "Cannot encode chunk id index >= 256 in this version!"  # TODO: handle this case! (use more bytes or different base packet...)
    to_xor[insertion_position] = struct.pack("<B", index_of_chunk)
    packet.data = xor_numpy(packet.data, to_xor)
    return packet


def find_insertion_position(version_string_length: int, chunk_diff: np.ndarray,
                            id_spacing: int, id_len: int) -> Generator[int, int, int]:
    """
    Find a suitable insertion position for the version string in the chunk diff.
    The position must be in a region of zeros of at least version_string_length.
    Returns the starting index of the insertion position, or -1 if not found.
    @param version_string_length: Length of the version string to insert. Length must be defined in stored bytes,
        not quats!
    @param chunk_diff: Numpy array representing the chunk diff.
    @param id_spacing: Spacing between chunk ID DNA-bases in the packet.
    @param id_len: Length of the chunk ID in bytes.
    @return: Starting index of the insertion position, or -1 if not found.
    """

    zero_regions = np.where(chunk_diff == 0)[0]
    for start_idx in reversed([x for x in zero_regions]):
        if start_idx + version_string_length <= len(chunk_diff) and id_spacing * id_len * 4 - 1 < start_idx:
            # TODO: check that the insertion position does not overlap with any seed positions (seed spacing) (NO off-by-one!)
            if np.all(chunk_diff[start_idx:start_idx + version_string_length] == 0):
                yield start_idx
    # Could not find suitable insertion position for version string!
    return -1


def generate_new_packets(semiautomatic_solver: SemiAutomaticReconstructionToolkit,
                         encoder: RU10Encoder,
                         diff: np.ndarray,
                         changed_chunk_ids: np.ndarray,
                         new_file_version: int) -> typing.Dict[int, typing.List[RU10Packet]]:
    """
    Scan ALL seeds and record ALL candidate seeds for each changed chunk.

    Optimizations:
    - Precompute changed set for O(1) membership checks.
    - Quick-skip seeds that don't mention any changed chunk.-
    - Quick-skip seeds that don't include header (0) or last chunk index before performing heavy removeAndXorAuxPackets.
    - Scan through ALL possible seeds to collect every valid candidate.

    Returns mapping {chunk_id: set(all_valid_seeds)} for downstream packet construction.
    """
    possible_seeds: typing.Dict[int, typing.Any] = {}
    max_num = min(encoder.calc_max_size(struct.calcsize("<" + encoder.id_len_format)), int31)

    last_chunk_idx = encoder.number_of_chunks - 1
    changed_set = set(int(x) for x in changed_chunk_ids)

    # Initialize mapping for results - one entry per changed chunk
    packet_to_seed_mapping: typing.Dict[int, typing.Set[int]] = {cid: set() for cid in changed_set}

    logger.info(f"Scanning {max_num} seeds for {len(changed_set)} changed chunks...")

    for seed in range(max_num):
        packet_numbers = choose_packet_numbers(encoder.number_of_chunks, seed, encoder.dist, systematic=False)

        # OPTIMIZATION 1: Quick filter using set intersection - skip if no changed chunks touched
        pn_set = set(packet_numbers)
        if not (pn_set & changed_set):
            continue

        # OPTIMIZATION 2: Skip if doesn't initially include header or last chunk
        if 0 not in pn_set and last_chunk_idx not in pn_set:
            continue

        # Now perform the expensive decoder call
        packet = RU10Packet(b"", packet_numbers, encoder.number_of_chunks, seed, encoder.dist, read_only=True)
        used_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(packet)

        # Verify header or last chunk remains after aux removal
        try:
            if not (bool(used_chunks[0]) or bool(used_chunks[-1])):
                continue
        except Exception:
            continue

        # Find which chunks this seed actually uses (after aux removal)
        try:
            true_indices = set(int(x) for x in np.nonzero(np.asarray(used_chunks))[0])
        except Exception:
            true_indices = {i for i, v in enumerate(used_chunks) if v}

        # Find intersection with changed chunks
        touched_changed = true_indices & changed_set
        if touched_changed:
            for cid in touched_changed:
                packet_to_seed_mapping[cid].add(seed)
            possible_seeds[seed] = used_chunks

        # Log progress every 10000 seeds
        if seed % 10000 == 0 and seed > 0:
            logger.info(
                f"Scanned {seed}/{max_num} seeds... Found {sum(len(v) for v in packet_to_seed_mapping.values())} total candidates")

    # Verify we found at least one seed for every changed chunk
    missing = [cid for cid, seeds in packet_to_seed_mapping.items() if not seeds]
    if missing:
        raise RuntimeError(f"Could not find ANY seed for chunks {missing} that contains either header or last chunk")

    logger.info(f"Scan complete! Total candidates found: {sum(len(v) for v in packet_to_seed_mapping.values())}")

    # Return the mapping for downstream processing (you can implement packet selection next)
    chunk_to_potential_seed_mapping = select_numbers(packet_to_seed_mapping,
                                                     n=10)  # TODO: mapping between packet numbers and results!
    generated_packets: typing.Dict[int, typing.List[RU10Packet]] = {}
    for tpl in chunk_to_potential_seed_mapping:
        chunk_id, seed_set = tpl
        for seed in seed_set:
            # create packets with the  a packet with the given seed and store it
            # check error probability: if it is 0, we can use it directly, otherwise we use the best packet
            packet = encoder.create_new_packet(False, seed)
            should_drop_packet(encoder.rules, packet, 1.0)  # calculate the errorprobability
            if chunk_id not in generated_packets:
                generated_packets[chunk_id] = []
            generated_packets[chunk_id].append(packet)
        generated_packets[chunk_id] = sorted(generated_packets[chunk_id], key=lambda x: x.error_prob)

    semiautomatic_solver.decoder.populate_header_chunk()

    # TODO: XOR the diff of old and current version of the chunk content for each changed chunk included in the packet
    for changed_chunk, potential_packets in generated_packets.items():
        for potential_packet in potential_packets:
            plain_used_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(potential_packet)
            version_string = generate_dna_version_string(new_file_version)
            if plain_used_chunks[0]:
                insertion_position = next(find_insertion_position(int((len(version_string) + 1) / 2), diff[
                    changed_chunk],
                                                                  packet))  # we need room for the version string and one bytes for the chunk id (as an offset!
                # we must update the header chunk
                assert insertion_position != -1, "Could not find insertion position for version string!"  # TODO: handle this case!
                insert_dna_version_string(potential_packet,
                                          version_string, insertion_position,
                                          diff[
                                              changed_chunk])  # pad with one A as we stripped it during version-string encoding!
                # TODO: add the version_string (replace end of each generated packet with version string - use metadata_coding)
                #  - convert version string to binary and replace end of potential_packet.data with this.
                insert_id_string(potential_packet, insertion_position + (len(version_string) + 1) / 2,
                                 changed_chunk,
                                 semiautomatic_solver)  # must be done such that the decoding is performed with all OLD chunks!
                # TODO:
                #  chunk_0 = forall id in used_chunks [except chunk 0] of new_packet: new_packet = new_packet XOR id
                #  for this, we must know the diff between the original changed chunk and the new version! (apply the diff!)
            elif plain_used_chunks[-1]:
                # FIXME: we must update the last chunk. Skipping for now...
                raise NotImplementedError("not handled yet!")
                new_packet = insert_version_string_last_chunk(potential_packet,
                                                              tranlate_quat_to_byte(version_string))
            # apply the diff for the changed chunk
            potential_packet.data = xor_numpy(potential_packet.data, diff[changed_chunk])

            org_headerchunk: HeaderChunk = semiautomatic_solver.decoder.headerChunk
            org_headerchunk.data

    # TODO: change a padding region in the header to the the id of the packet to change
    #  (decoding of that header should be performed with knowledge of the OLD file version, not the new one!)
    #  this can be achieved by using the new packet to solve chunk 0 using the already solved chunks (except chunk 0)
    #  the packet then has to be repaired by XORing the diff between the expected header and the actual header
    # TODO: during encoding the DNA value for the chunk id-field has to be calculated using the packet generated using
    #  the OLD chunk!

    # TODO: during decoding of the new chunk (id found in previous step) we reduce the same packet towards the chunk id
    #  previously retrieved while also using the clean header / last chunk for reduction.
    #  this can be done on the fly without performing a full decode

    # add the version-string to
    return generated_packets


def add_packets(encoder: RU10Encoder, new_packets: typing.Set[RU10Packet]):
    encoder.encodedPackets += new_packets


if __name__ == "__main__":
    print(generate_dna_version_string())

    parsed_args = init_args()

    ini_file = parsed_args.ini
    new_file = parsed_args.new_file
    packet_add_limit = parsed_args.packet_add_limit
    """
    if parsed_args.wanted_metadata_file is not None:
        wanted_metadata = parse_metadata_file(parsed_args.wanted_metadata_file)
    else:
        wanted_metadata = parsed_args.wanted_metadata.split(",")

    if parsed_args.unwanted_metadata_file is not None:
        unwanted_metadata = parse_metadata_file(parsed_args.unwanted_metadata_file)
    else:
        unwanted_metadata = parsed_args.unwanted_metadata.split(",")
    """

    cfg_worker = ConfigReadAndExecute(ini_file)
    x = cfg_worker.execute(return_decoder=True, skip_solve=True)[0]
    semiautomatic_solver = SemiAutomaticReconstructionToolkit(x)

    diff, changed_chunks = find_affected_chunks(semiautomatic_solver, new_file)

    new_file_version = get_current_file_version(semiautomatic_solver) + 1

    encoder = encoder_from_decoder(semiautomatic_solver, cfg_worker, rules=FastDNARules())

    packet_candidates = generate_new_packets(semiautomatic_solver, encoder, diff, changed_chunks, new_file_version)

    logger.info(f"Found candidate seeds per chunk: {dict((k, len(v)) for k, v in packet_candidates.items())}")

    for changed_chunk_packet_group, packets in packet_candidates.items():
        added_packets_for_chunk = 0
        for packet in packets:
            if added_packets_for_chunk >= packet_add_limit:
                break
            encoder.encodedPackets.add(packet)

    # TODO: implement packet selection and encoding logic here
    # For now, just log what was found
    logger.warning("Packet encoding not yet implemented - returning candidate mapping")

    # encoder.save_packets_fasta("1test_out.fasta", "", False)

    # add_packets(encoder, new_packets)

    # new_packets = create_new_metadata_packets(semiautomatic_solver, wanted_metadata, unwanted_metadata, True, True)

    # replace_packets(encoder, new_packets, 0.0, True, True)

    logger.error(encoder.encodedPackets)

    encoder.save_packets_fasta(f"file_out_v{new_file_version}.fasta", "", False)
    encoder.save_config_file(section_name=f"file_out_v{new_file_version}.fasta")
    logger.warning(semiautomatic_solver.decoder.packets)
