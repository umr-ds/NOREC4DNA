import argparse
import logging
import struct
import typing
from typing import Generator

import numpy as np
from PIL import Image

from NOREC4DNA.norec4dna.helper.helper_cpu_single_core import should_drop_packet, xor_numpy
from NOREC4DNA.norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.helper.RU10Helper import int31, choose_packet_numbers, from_true_false_list
from NOREC4DNA.ConfigWorker import ConfigReadAndExecute
from NOREC4DNA.metadata_coding import encoder_from_decoder
from NOREC4DNA.norec4dna import RU10Encoder
from NOREC4DNA.norec4dna.RU10Packet import RU10Packet
from norec4dna.helper.helper import xor_with_seed
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
    """
    TODO: convert the phash to DNA sequence
    """
    hash = imagehash.phash(image)
    return str(hash)


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
    """
    Adds a DNA version string at a specific position in an RU10Packet.
    :param packet: The RU10Packet to operate on.
    :param dna_version_string: The DNA string to insert (e.g., "ACGT").
    :param insertion_position: Byte position in the data section (not the DNA representation!).
    :param diff: array containing the diff between original and new packet data; used only for asserting correct insertion position.
    """
    if packet.dna_data is None or packet.dna_data == "":
        packet.calculate_packed_data()
    assert not diff[insertion_position:insertion_position + int((len(dna_version_string) + 1) / 2) + 1].any(), \
        "Cannot insert version string into non-padding region!"

    pre_padding = "A" * ((4 - len(dna_version_string) % 4) % 4)
    binary_dna_version_string = tranlate_quat_to_byte(pre_padding + dna_version_string)
    num_padding_bits = len(pre_padding) * 2
    # calculate the diff between the existing DNA at the selected position and the wanted DNA version string to insert
    # convert this to binary and XOR the existing data with this diff
    if num_padding_bits > 0:
        original_byte = packet.data[insertion_position]
        mask_keep = (0xFF << (8 - num_padding_bits)) & 0xFF
        mask_new = (0xFF >> num_padding_bits) & 0xFF
        combined_first_byte = (original_byte & mask_keep) | (binary_dna_version_string[0] & mask_new)
        binary_dna_version_string = bytes([combined_first_byte]) + binary_dna_version_string[1:]

    # get the current dna sequence at the insertion-position:
    """
    current_dna_seq = packet.dna_data[(packet.get_packet_header_size())*4:][(insertion_position) * 4:(insertion_position) * 4 + len(dna_version_string)]
    # as we may have xor_by_seed enabled, we apply the new content as a diff between the existing data and the wanted data:
    existing_binary = tranlate_quat_to_byte(pre_padding + current_dna_seq)
    wanted_binary = tranlate_quat_to_byte(pre_padding + dna_version_string)
    binary_diff = bytes(a ^ b for a, b in zip(existing_binary, wanted_binary))
    # xor the diff into packet.data at the correct position:
    #packet.data[insertion_position:insertion_position + len(binary_diff)] ^= np.frombuffer(binary_diff, dtype=np.uint8)
    """

    stub = np.zeros_like(packet.data)
    # TODO: use lenth of DNA-version-string to determine length here as length of binary string may be different due to omitted padding if starting with 0x00!
    stub[insertion_position:insertion_position + len(binary_dna_version_string)] = np.frombuffer(
        binary_dna_version_string, dtype=np.uint8)
    stub = xor_with_seed(stub, packet.id)
    # packet.data = xor_numpy(packet.data, stub)
    packet.data[insertion_position:insertion_position + len(binary_dna_version_string)] = np.frombuffer(
        stub[insertion_position:insertion_position + len(binary_dna_version_string)], dtype=np.uint8)
    packet.get_dna_struct(True, packet.id_spacing, packet.id_spacing_length, True)


#  add the version_string (replace end of each generated packet with version string - use metadata_coding)
#  - convert version string to binary and replace end of potential_packet.data with this.
def reduce_packet_to_chunk(packet: RU10Packet,
                           semiautomatic_solver: SemiAutomaticReconstructionToolkit, chunk_to_reduce_to=0):
    # must be done such that the decoding is performed with all OLD chunks!
    # get normalized used chunks and remote chunk
    used_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(packet)
    packet.used_packets = {i for i, x in enumerate(used_chunks) if x}
    for i, must_remove_chunk in enumerate(used_chunks):
        if must_remove_chunk and i != chunk_to_reduce_to:
            packet.data = xor_numpy(packet.data, semiautomatic_solver.decoder.GEPP.b[i])
            packet.used_packets.remove(i)
    assert len(packet.used_packets) == 1 and chunk_to_reduce_to in packet.used_packets
    return packet
    # packet contains only the header chunk
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
    to_xor = np.zeros_like(packet.data, dtype=np.uint8)
    used_chunks = from_true_false_list(semiautomatic_solver.decoder.removeAndXorAuxPackets(packet))
    assert changed_chunk_id in used_chunks, "Changed chunk id must be in used chunks!"  # TODO: handle this case!
    index_of_chunk = used_chunks.index(changed_chunk_id)
    assert index_of_chunk < 256, "Cannot encode chunk id index >= 256 in this version!"  # TODO: handle this case! (use more bytes or different base packet...)
    to_xor[insertion_position] = index_of_chunk & 0xFF
    packet.data = xor_numpy(packet.data, to_xor)
    # recalculate dna_data:
    packet.get_dna_struct(True, packet.id_spacing, packet.id_spacing_length, True)
    """
    res = RU10Packet(packet.data.copy(), packet.used_packets, packet.total_number_of_chunks, packet.id, dist=packet.dist,
                   read_only=False, error_correction=packet.error_correction,
                   packet_len_format=packet.packet_len_format, crc_len_format=packet.crc_len_format,
                   number_of_chunks_len_format=packet.number_of_chunks_len_format,
                   id_len_format=packet.id_len_format,
                   save_number_of_chunks_in_packet=packet.save_number_of_chunks_in_packet, prepend=packet.prepend,
                   append=packet.append, xor_by_seed=packet.xor_by_seed, mask_id=packet.mask_id,
                   id_spacing=packet.id_spacing)
    #assert dna_version_string in res.dna_data, "Version string insertion failed!"
    packet.dna_data = res.dna_data
    packet.data = res.data
    packet.packed = res.packed
    packet.packed_used_packets = res.packed_used_packets
    """


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


def find_insertion_position_with_seed(version_string_length: int, chunk_diff: np.ndarray,
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
    # Create a mask for valid positions (True = usable, False = blocked by seed)
    num_bytes = len(chunk_diff)
    num_quats = num_bytes * 4  # Each byte = 4 DNA bases (quaternary)

    # Create quaternary-level mask (True = usable)
    quat_mask = np.ones(num_quats, dtype=bool)

    if id_spacing > 0 and id_len > 0:
        # Number of seed bases (each byte = 4 DNA bases)
        num_seed_bases = id_len * 4
        # Seed bases are placed at positions 0, (id_spacing+1), 2*(id_spacing+1), ...
        # until all num_seed_bases are placed
        seed_step = id_spacing + 1
        for i in range(num_seed_bases):
            seed_quat_pos = i * seed_step
            if seed_quat_pos < num_quats:
                quat_mask[seed_quat_pos] = False

    # Convert quaternary mask to byte mask
    # A byte is only usable if ALL 4 of its quaternary positions are usable
    byte_mask = np.ones(num_bytes, dtype=bool)
    for byte_idx in range(num_bytes):
        quat_start = byte_idx * 4
        quat_end = min(quat_start + 4, num_quats)
        if not np.all(quat_mask[quat_start:quat_end]):
            byte_mask[byte_idx] = False

    # Combine with chunk_diff == 0 condition
    # Valid positions: chunk_diff is zero AND byte is not blocked by seed spacing
    valid_positions = (chunk_diff == 0) & byte_mask

    zero_regions = np.where(valid_positions)[0]
    for start_idx in reversed([x for x in zero_regions]):
        if start_idx + version_string_length <= len(chunk_diff):
            if np.all(valid_positions[start_idx:start_idx + version_string_length]):
                yield start_idx
    # Could not find suitable insertion position for version string!
    return -1


def generate_new_packets(semiautomatic_solver: SemiAutomaticReconstructionToolkit,
                         encoder: RU10Encoder,
                         diff: np.ndarray,
                         changed_chunk_ids: np.ndarray,
                         new_file_version: int) -> typing.Dict[
    int, typing.Union[typing.List[RU10Packet], typing.List[typing.Tuple[RU10Packet, RU10Packet]]]]:
    """
    Scan ALL seeds and record ALL candidate seeds for each changed chunk.

    Optimizations:
    - Precompute changed set for O(1) membership checks.
    - Quick-skip seeds that don't mention any changed chunk.-
    - Quick-skip seeds that don't include header (0) or last chunk index before performing heavy removeAndXorAuxPackets.
    - Scan through ALL possible seeds to collect every valid candidate.

    Returns mapping {chunk_id: set(all_valid_seeds)} for downstream packet construction.
    """
    res: typing.Dict[int, typing.List[typing.Union[RU10Packet, typing.Tuple[RU10Packet, RU10Packet]]]] = {}
    possible_seeds: typing.Dict[int, typing.Any] = {}
    max_num = min(encoder.calc_max_size(struct.calcsize("<" + encoder.id_len_format)), int31)

    last_chunk_idx = encoder.number_of_chunks - 1
    changed_set = {int(x) for x in changed_chunk_ids}

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
            true_indices = {int(x) for x in np.nonzero(np.asarray(used_chunks))[0]}
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
                                                     n=10)  # mapping between packet numbers and results!
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

    changed_chunk_to_new_packets: typing.Dict[int, typing.List[RU10Packet]] = {}
    changed_chunk_to_packet_pair_list: typing.Dict[int, typing.List[typing.Tuple[
        RU10Packet, RU10Packet]]] = {}  # list of (packet1, packet2) tuples, each applying half of the diff with the same version string and chunk id insertion at
    version_string = generate_dna_version_string(new_file_version)
    # XOR the diff of old and current version of the chunk content for each changed chunk included in the packet
    for changed_chunk, potential_packets in generated_packets.items():
        packet_added = False
        # generate all potential packets for this changed chunk and calculate the error probability; use the best one!
        for potential_packet in potential_packets:
            modified_packet = potential_packet.copy()
            plain_used_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(modified_packet)
            if plain_used_chunks[0]:
                insertion_position = next(find_insertion_position(int((len(version_string) + 1) / 4 + 1), diff[
                    changed_chunk], modified_packet.id_spacing, struct.calcsize(
                    modified_packet.id_len_format)))  # we need room for the version string and one bytes for the chunk id (as an offset!
                # we must update the header chunk
                if insertion_position == -1:
                    # skip this potential packet as we cannot insert the version string!
                    continue
                assert insertion_position != -1, "Could not find insertion position for version string!"  # TODO: handle this case!
                # add the version_string (replace end of each generated packet with version string - use metadata_coding)
                #  - convert version string to binary and replace end of modified_packet.data with this:
                insert_dna_version_string(modified_packet, version_string, insertion_position, diff[changed_chunk])

                # insert the id (offset from the chunk list) of the changed chunk
                # must be done such that the decoding is performed with all OLD chunks!
                insert_id_string(modified_packet, int(insertion_position + (len(version_string) + 1) / 4),
                                 changed_chunk,
                                 semiautomatic_solver)
            elif plain_used_chunks[-1]:
                # TODO: to increase the chance of low overhead insertion, we can update the padding of the last chunk or use ANY other chunk (e.g. always the "lowest" used chunk).
                continue
                # new_packet = insert_version_string_last_chunk(modified_packet,
                #                                              tranlate_quat_to_byte(version_string))
            # apply the diff for the changed chunk
            modified_packet.data = xor_numpy(modified_packet.data, diff[changed_chunk])
            # recalculate dna_data:
            modified_packet.get_dna_struct(True, modified_packet.id_spacing, modified_packet.id_spacing_length, True)

            should_drop_packet(encoder.rules, modified_packet)  # populate .error_prob
            if changed_chunk not in changed_chunk_to_new_packets:
                changed_chunk_to_new_packets[changed_chunk] = []
            changed_chunk_to_new_packets[changed_chunk].append(modified_packet)
            packet_added = True

        # we did not find a packet with a suitable insertion position for the version string (or the error prob is too high)!
        if not packet_added or sorted(changed_chunk_to_new_packets[changed_chunk], key=lambda x: x.error_prob)[
            0].error_prob >= 1.0:
            logger.warning(
                f"Could not find suitable packet for chunk {changed_chunk} to insert version string. Generating two packets with same seed with split diff (50/50)!")
            for potential_packet in potential_packets:
                modified_packet_first = potential_packet.copy()
                modified_packet_second = potential_packet.copy()
                # split the diff in half using a mask and & operation:
                diff_mask = np.zeros_like(diff[changed_chunk], dtype=bool)
                half_point = len(diff[changed_chunk]) // 2
                diff_mask[:half_point] = True
                first_diff = np.where(diff_mask, diff[changed_chunk], 0)
                second_diff = np.where(~diff_mask, diff[changed_chunk], 0)

                # apply the diff for the changed chunk
                if plain_used_chunks[0]:
                    insertion_position_first = next(
                        find_insertion_position(int((len(version_string) + 1) / 4 + 1), diff[
                            changed_chunk], modified_packet_first.id_spacing, struct.calcsize(
                            modified_packet_first.id_len_format)))

                    insertion_position_second = next(
                        find_insertion_position(int((len(version_string) + 1) / 4 + 1), diff[
                            changed_chunk], modified_packet_second.id_spacing, struct.calcsize(
                            modified_packet_second.id_len_format)))
                    # we must update the header chunk
                    if insertion_position_first == -1 or insertion_position_second == -1:
                        # skip this potential packet as we cannot insert the version string!
                        # more likely to happen for insertion_position_first as this should only happen if we have a large seed-spacing value!
                        # TODO: we may try to insert the version string into the gap between two seed positions if seed-spacing is large enough!
                        continue
                    assert insertion_position_first != -1 and insertion_position_second != -1, "Could not find insertion position for version string! This is due to insufficient space free space in the first modified diff version. (Most likely because the seed-spacing value is too big)!"
                    # operate on the first packet:
                    insert_dna_version_string(modified_packet_first, version_string, insertion_position_first,
                                              first_diff)
                    insert_id_string(modified_packet_second,
                                     int(insertion_position_second + (len(version_string) + 1) / 4),
                                     changed_chunk, semiautomatic_solver)

                    # same for second packet:
                    insert_dna_version_string(modified_packet_first, version_string, insertion_position_first,
                                              second_diff)
                    insert_id_string(modified_packet_second,
                                     int(insertion_position_second + (len(version_string) + 1) / 4),
                                     changed_chunk, semiautomatic_solver)
                elif plain_used_chunks[-1]:
                    # TODO: to increase the chance of low overhead insertion, we can update the padding of the last chunk or use ANY other chunk (e.g. always the "lowest" used chunk).
                    continue
                modified_packet_first.data = xor_numpy(modified_packet_first.data, first_diff)
                modified_packet_first.get_dna_struct(True, modified_packet_first.id_spacing,
                                                     modified_packet_first.id_spacing_length,
                                                     True)

                modified_packet_second.data = xor_numpy(modified_packet_second.data, second_diff)
                modified_packet_first.get_dna_struct(True, modified_packet_second.id_spacing,
                                                     modified_packet_second.id_spacing_length,
                                                     True)

                should_drop_packet(encoder.rules, modified_packet_first)  # populate .error_prob
                should_drop_packet(encoder.rules, modified_packet_second)

                changed_chunk_to_new_packets[changed_chunk].append(modified_packet_first)
                changed_chunk_to_new_packets[changed_chunk].append(modified_packet_second)
                if changed_chunk not in changed_chunk_to_packet_pair_list:
                    changed_chunk_to_packet_pair_list[changed_chunk] = []
                changed_chunk_to_packet_pair_list[changed_chunk].append((modified_packet_first, modified_packet_second))
                packet_added = True
        if not packet_added:
            logger.error(f"Could not create any packet for changed chunk {changed_chunk}!")
            # recreate the dna_data:
            """
            packet = RU10Packet(modified_packet.data.copy(), modified_packet.used_packets,
                            modified_packet.total_number_of_chunks, modified_packet.id, dist=modified_packet.dist,
                            read_only=False, error_correction=modified_packet.error_correction,
                            packet_len_format=modified_packet.packet_len_format,
                            crc_len_format=modified_packet.crc_len_format,
                            number_of_chunks_len_format=modified_packet.number_of_chunks_len_format,
                            id_len_format=modified_packet.id_len_format,
                            save_number_of_chunks_in_packet=modified_packet.save_number_of_chunks_in_packet,
                            prepend=modified_packet.prepend,
                            append=modified_packet.append, xor_by_seed=modified_packet.xor_by_seed,
                            mask_id=modified_packet.mask_id,
                            id_spacing=modified_packet.id_spacing)
            modified_packet.dna_data = packet.dna_data
            #modified_packet.packed = packet.packed
            #modified_packet.packed_used_packets = packet.packed_used_packets

            #org_headerchunk: HeaderChunk = semiautomatic_solver.decoder.headerChunk
            #org_headerchunk.data
            """
    # select the best packet(s) for each changed chunk (based on error probability):
    for changed_chunk in changed_chunk_to_new_packets.keys():

        if changed_chunk in changed_chunk_to_new_packets:
            # try using best single packet first:
            changed_chunk_to_new_packets[changed_chunk] = sorted(changed_chunk_to_new_packets[changed_chunk],
                                                                 key=lambda x: x.error_prob)
            for new_pack in changed_chunk_to_new_packets[changed_chunk]:
                if new_pack.error_prob < 1.0:
                    if changed_chunk not in res:
                        res[changed_chunk] = []
                    res[changed_chunk].append(new_pack)
                    logger.info(
                        f"Generated packet with error probability {new_pack.error_prob} for chunk {changed_chunk}.")
                else:
                    logger.warning(
                        f"Skipping packet with error probability {new_pack.error_prob} for chunk {changed_chunk}!")
        elif changed_chunk in changed_chunk_to_packet_pair_list:
            # there was no valid single packet, use best packet pair:
            changed_chunk_to_packet_pair_list[changed_chunk] = sorted(changed_chunk_to_packet_pair_list[changed_chunk],
                                                                      key=lambda x: max(x[0].error_prob,
                                                                                        x[1].error_prob))
            for new_pack_pair in changed_chunk_to_packet_pair_list[changed_chunk]:
                if new_pack_pair[0].error_prob < 1.0 and new_pack_pair[1].error_prob < 1.0:
                    if changed_chunk not in res:
                        res[changed_chunk] = []
                    res[changed_chunk].append(new_pack_pair)
                    logger.info(
                        f"Generated packet pair with error probability {new_pack_pair[0].error_prob, new_pack_pair[1].error_prob} for chunk {changed_chunk}!")
                else:
                    logger.warning(
                        f"Skipping packet pair with error probability {new_pack_pair[0].error_prob, new_pack_pair[1].error_prob} for chunk {changed_chunk}!")
        # logger.error("Could not find any valid packet (single or pair) for chunk {changed_chunk}!")

    # TODO: set a padding region in the header to the the id of the packet to change
    #  (decoding of that header should be performed with knowledge of the OLD file version, not the new one!)
    #  this can be achieved by using the new packet to solve chunk 0 using the already solved chunks (except chunk 0)
    #  the packet then has to be repaired by XORing the diff between the expected header and the actual header
    # TODO: during encoding the DNA value for the chunk id-field has to be calculated using the packet generated using
    #  the OLD chunk!

    # TODO: during decoding of the new chunk (id found in previous step) we reduce the same packet towards the chunk id
    #  previously retrieved while also using the clean header / last chunk for reduction.
    #  this can be done on the fly without performing a full decode

    # TODO: IMPORTANT: Theoretically, it is possible to use ANY chunk and not only the header / last chunk for the changed chunk insertion as long as we already know the prior version of the chunk!
    return res


def decode_versions(semiautomatic_solver: SemiAutomaticReconstructionToolkit, dna_version_string_prefix="") -> \
typing.Dict[int, str]:
    """
    Decode all versions of an encoded file from the packets known to the semiautomatic_solver.
    Returns a mapping {file_version: filename} for all versions found.
    param semiautomatic_solver: The SemiAutomaticReconstructionToolkit with the decoder containing the packets.
    param dna_version_string_prefix: The prefix used for the version strings in the DNA sequences.
    res: typing.Dict[int, str] = {}
    """
    # TODO: for base -version: decode using only packets that DO NOT contain the version string!
    #  then for each changed version:
    #  then decode all changed chunks using packets that DO contain the version string with the correct version number
    #  safe the new file version as v<version>_<base_filename>.<extension>
    #  then use the new version of the file as the base for the next version
    # semiautomatic_solver.
    pass


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

    logger.info(f"Found candidate seeds per chunk: {dict((k, v) for k, v in packet_candidates.items())}")

    packets_added = 0
    for changed_chunk_packet_group, packets in packet_candidates.items():
        added_packets_for_chunk = 0
        for packet in packets:
            if packet.error_prob >= 1.0:
                continue
            if added_packets_for_chunk >= packet_add_limit:
                break
            encoder.encodedPackets.add(packet)
            added_packets_for_chunk += 1
            packets_added += 1
        logger.info(f"Added {added_packets_for_chunk} packets for changed chunk {changed_chunk_packet_group}.")

    logger.info(f"Added total of {packets_added} packets for {len(packet_candidates)} changed chunks.")
    outfile = f"{semiautomatic_solver.decoder.file.split(".fasta")[0]}_v{new_file_version}"
    file_bkp = encoder.file
    encoder.file = outfile
    encoder.save_packets_fasta(None, "", False)
    encoder.save_config_file(add_dot_fasta=True)
    encoder.file = file_bkp
    logger.warning(semiautomatic_solver.decoder.packets)
