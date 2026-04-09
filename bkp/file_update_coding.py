r"""
File Update Coding for NOREC4DNA Multi-Version Support.

This module provides functionality for encoding file updates into DNA sequences,
enabling efficient multi-version support in NOREC4DNA-encoded files.

Key Features:
    1. Differential Encoding: Only encode changes between file versions
    2. Version String Insertion: Embed version identifiers in DNA packets
    3. Chunk-Level Diff: Calculate differences at the chunk level
    4. Optimal Packet Selection: Scan seeds to find best packets for changes
    5. Error Probability Calculation: Evaluate and select low-error packets
    6. Packet Pair Support: Split changes across multiple packets when needed

Example Usage:
    >>> from NOREC4DNA.file_update_coding import (
    ...     find_affected_chunks,
    ...     generate_new_packets,
    ...     generate_dna_version_string,
    ...     get_current_file_version
    ... )
    >>> from NOREC4DNA.ConfigWorker import ConfigReadAndExecute
    >>> from semi_automatic_reconstruction_toolkit import SemiAutomaticReconstructionToolkit
    >>>
    >>> # Load configuration
    >>> config = ConfigReadAndExecute("config.ini")
    >>> decoder = config.execute(return_decoder=True, skip_solve=True)[0]
    >>> solver = SemiAutomaticReconstructionToolkit(decoder)
    >>>
    >>> # Read new file version
    >>> with open("updated_file.pdf", "rb") as f:
    ...     new_file_data = f.read()
    >>>
    >>> # Calculate diff and find changed chunks
    >>> diff, changed_chunks = find_affected_chunks(solver, new_file_data)
    >>> print(f"Changed chunks: {changed_chunks}")
    >>>
    >>> # Get current version and generate version string
    >>> current_version = get_current_file_version(solver)
    >>> new_version = current_version + 1
    >>> version_string = generate_dna_version_string(new_version)
"""

import argparse
import logging
import struct
import typing
from pathlib import Path
from typing import Generator, Dict, List, Tuple, Set, Optional, Union

import numpy as np
from PIL import Image

import imagehash

from ConfigWorker import ConfigReadAndExecute
from metadata_coding import encoder_from_decoder
from norec4dna import RU10Encoder
from norec4dna.RU10Packet import RU10Packet
from norec4dna.helper.bin2Quaternary import byte2QUATS
from norec4dna.helper.helper import xor_with_seed
from norec4dna.helper.helper_cpu_single_core import (
    should_drop_packet,
    xor_numpy,
)
from norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.helper.RU10Helper import choose_packet_numbers, from_true_false_list, int31
from repair_algorithms.utils.select_numbers import select_numbers
from semi_automatic_reconstruction_toolkit import SemiAutomaticReconstructionToolkit

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def find_affected_chunks(
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    new_file: typing.Optional[bytes] = None,
) -> typing.Tuple[np.ndarray, np.ndarray]:
    """
    Calculate diff between existing chunks and new file.

    This function compares the current decoded file state with a new file version
    to identify which chunks have changed. It returns both the raw difference array
    and the indices of differing chunks.

    Args:
        semiautomatic_solver: SemiAutomaticReconstructionToolkit instance with
            decoder state.
        new_file: New file bytes to compare against. If None, creates artificial
            changes at positions 200 and 500 for testing.

    Returns:
        A tuple containing:
        - diff: Numpy array of differences between old and new chunks
        - differing_rows: Array of chunk indices that have changed

    Raises:
        ValueError: If file size doesn't match expected chunk layout.

    Note:
        Currently requires the new file to have equal length or be shorter than
        the original. Index is defined WITHOUT header row as the header must be
        changed for any changed content.

    Example:
        >>> diff, changed = find_affected_chunks(solver, new_file_data)
        >>> print(f"Chunks {changed} have changed")
    """
    # Get current file content from decoder
    current_file_bytes = semiautomatic_solver.get_file_as_bytes()

    # For testing/development: create artificial changes if no new_file provided
    if new_file is None:
        new_file = bytearray(current_file_bytes)
        if len(new_file) > 500:
            new_file[500] = 0x00
            new_file[200] = 0x00
            logger.debug("Created artificial diff for testing")

    # Validate file size - new file must not be larger than current
    if len(new_file) > len(current_file_bytes):
        raise ValueError(
            f"New file ({len(new_file)} bytes) is larger than current file "
            f"({len(current_file_bytes)} bytes). File growth not yet supported."
        )

    # Pad new_file to match current file size if smaller
    if len(new_file) < len(current_file_bytes):
        logger.debug(
            f"Padding new file from {len(new_file)} to {len(current_file_bytes)} bytes"
        )
        new_file_padded = bytearray(new_file)
        new_file_padded.extend(b'\x00' * (len(current_file_bytes) - len(new_file)))
        new_file = bytes(new_file_padded)

    # Get chunk size from GEPP
    chunk_size = semiautomatic_solver.decoder.GEPP.b.shape[1]

    # Create compatible numpy array with correct chunk sizes from bytes
    # Skip header chunk (index 0) as we compare data chunks only
    new_array = np.frombuffer(new_file, dtype=np.uint8).reshape(
        -1, chunk_size
    )

    # Calculate difference between old and new chunks (skip header chunk at index 0)
    # GEPP.b[1:] contains data chunks 0 to n-1 (where n is number_of_chunks - 1)
    diff = semiautomatic_solver.decoder.GEPP.b[1:] - new_array

    # Find indices of differing rows (chunks)
    differing_rows = np.where(np.any(diff != 0, axis=1))[0]

    logger.info(f"Found {len(differing_rows)} changed chunks")

    return diff, differing_rows


def update_file(
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    encoded_packets: Set[RU10Packet],
    new_file: bytes,
    packets_to_generate: int,
    updated_header: bytes,
) -> List[RU10Packet]:
    """
    Generate encoded packets representing changed content.

    This function calculates the changed chunks between the current file version
    and a new version, then generates DNA packets that encode only the differences.

    Args:
        semiautomatic_solver: SemiAutomaticReconstructionToolkit instance.
        encoded_packets: Set of existing encoded packets.
        new_file: Path to new file or bytes of new file content.
        packets_to_generate: Number of packets to generate per changed chunk.
        updated_header: Updated header bytes.

    Returns:
        A list of generated RU10Packet objects representing the changes.

    Raises:
        ValueError: If file update cannot be performed.

    Note:
        This is the main entry point for creating version updates. It orchestrates
        the entire process of diff calculation, packet generation, and version
        string insertion.

    Example:
        >>> packets = update_file(
        ...     solver,
        ...     existing_packets,
        ...     new_file_data,
        ...     packets_to_generate=5,
        ...     updated_header=header_bytes
        ... )
    """
    # Handle both file path and bytes input
    if isinstance(new_file, str):
        with open(new_file, "rb") as f:
            new_file_data = f.read()
    else:
        new_file_data = new_file

    # Calculate diff and find changed chunks
    diff, changed_chunks = find_affected_chunks(semiautomatic_solver, new_file_data)

    # Get current version and increment
    current_version = get_current_file_version(semiautomatic_solver)
    new_version = current_version + 1

    # Create encoder from decoder configuration
    encoder = encoder_from_decoder(
        semiautomatic_solver,
        ConfigReadAndExecute(semiautomatic_solver.decoder.file),
        rules=FastDNARules(),
    )

    # Generate packets for changed content
    packet_candidates = generate_new_packets(
        semiautomatic_solver, encoder, diff, changed_chunks, new_version
    )

    # Select best packets based on error probability
    generated_packets: List[RU10Packet] = []
    packets_added = 0

    for chunk_id, packets in packet_candidates.items():
        added_for_chunk = 0
        for packet in packets:
            if packet.error_prob >= 1.0:
                logger.warning(
                    f"Skipping packet for chunk {chunk_id} with error prob {packet.error_prob}"
                )
                continue

            if added_for_chunk >= packets_to_generate:
                break

            generated_packets.append(packet)
            added_for_chunk += 1
            packets_added += 1

        logger.info(
            f"Added {added_for_chunk} packets for changed chunk {chunk_id}"
        )

    logger.info(
        f"Generated {packets_added} total packets for {len(packet_candidates)} changed chunks"
    )

    return generated_packets


def create_perceptual_hash(image: Image.Image) -> str:
    """
    Create a perceptual hash for an image.

    This function generates a perceptual hash (phash) that can be used to
    compare images for similarity. The hash is returned as a string representation.

    Args:
        image: PIL Image object to hash.

    Returns:
        String representation of the perceptual hash.

    Note:
        Perceptual hashes are useful for detecting similar images even after
        minor modifications like compression or resizing.

    Example:
        >>> from PIL import Image
        >>> img = Image.open("image.png")
        >>> hash_value = create_perceptual_hash(img)
        >>> print(f"Image hash: {hash_value}")
    """
    hash_value = imagehash.phash(image)
    return str(hash_value)


def generate_dna_version_string(
    fileversion: int = 0, base_length: int = 8
) -> str:
    """
    Generate a DNA version string for embedding in packets.

    Creates a DNA sequence that encodes the file version number followed by
    a magic string marker. The version is encoded in the 3 bases immediately
    preceding the magic string.

    Args:
        fileversion: Version number to encode (must be in range [0, 63]).
        base_length: Length of version encoding in bases (default: 8).

    Returns:
        DNA version string (e.g., "AAAGAGCCAGTGAGTCGTA" for version 0).

    Raises:
        RuntimeError: If fileversion is >= 64 (exceeds 6-bit encoding capacity).
        ValueError: If fileversion is negative.

    Example:
        >>> v0 = generate_dna_version_string(0)
        >>> print(v0)  # "AAAGAGCCAGTGAGTCGTA"
        >>> v1 = generate_dna_version_string(1)
        >>> print(v1)  # "AACGAGCCAGTGAGTCGTA"
    """
    if fileversion < 0:
        raise ValueError("fileversion must be non-negative")
    if fileversion >= 64:
        raise RuntimeError("fileversion must be in range [0, 63]!")

    # Encode version in 3 bases (6 bits = 64 possible values)
    version_bases = byte2QUATS(fileversion)[1:]

    # Magic string marker
    magic_string = "GAGCCAGTGAGTCGTA"

    return version_bases + magic_string


def init_args() -> argparse.Namespace:
    """
    Parse command-line arguments for file_update_coding.

    Returns:
        Parsed arguments namespace.

    Example:
        >>> args = init_args()
        >>> print(f"Processing {args.new_file}")
    """
    parser = argparse.ArgumentParser(
        description="File Update Coding for NOREC4DNA multi-version support"
    )
    parser.add_argument(
        "--ini",
        metavar="ini",
        type=str,
        help="Configuration file (INI format)",
        default="/home/michael/Code/DR4DNA/eval/sleeping_beauty_Mon_Feb_16_13_45_57_2026.ini",
    )
    parser.add_argument(
        "--new_file",
        metavar="new_file",
        type=str,
        help="Updated file path",
        required=True,
    )
    parser.add_argument(
        "--packet_add_limit",
        metavar="packet_add_limit",
        type=int,
        default=5,
        help="Maximum number of packets to add for each changed chunk",
    )

    return parser.parse_args()


def get_current_file_version(
    semi_automatic_solver: SemiAutomaticReconstructionToolkit,
    magic_string: str = "GAGCCAGTGAGTCGTA",
) -> int:
    """
    Determine the current version number from encoded DNA sequences.

    Scans all DNA sequences in the pool and extracts the highest version number
    found encoded before the magic_string marker.

    Args:
        semi_automatic_solver: SemiAutomaticReconstructionToolkit instance.
        magic_string: Magic DNA string marking version sequences
            (default: "GAGCCAGTGAGTCGTA").

    Returns:
        Highest version number found. Returns 0 if no version strings are present.

    Example:
        >>> version = get_current_file_version(solver)
        >>> print(f"Current version: {version}")
    """
    from invivo_window_decoder import load_fasta

    fasta_content = load_fasta(semi_automatic_solver.decoder.file)

    # Normalize to an iterable of sequence strings
    if hasattr(fasta_content, "values"):
        sequences = fasta_content.values()
    else:
        # Assume list/iterable of tuples like (header, seq)
        try:
            sequences = (entry[1] for entry in fasta_content)
        except Exception:
            sequences = (str(entry) for entry in fasta_content)

    # Check if magic_string is at the end of any sequence
    found_strings: Set[str] = set()
    for seq in sequences:
        if seq is None:
            continue
        if isinstance(seq, bytes):
            seq_str = seq.decode(errors="ignore")
        else:
            seq_str = str(seq)

        if seq_str.endswith(magic_string):
            found_strings.add(seq_str)

    if len(found_strings) == 0:
        return 0

    # Extract version numbers from found sequences
    highest_version = 0
    for found_string in found_strings:
        try:
            # Extract 3 bases before magic string
            version_bases = found_string[
                len(found_string) - 3 - len(magic_string) : -len(magic_string)
            ]
            version_num = ord(tranlate_quat_to_byte(f"A{version_bases}"))
            highest_version = max(highest_version, version_num)
        except Exception as e:
            logger.warning(f"Failed to parse version from sequence: {e}")

    return highest_version


def insert_dna_version_string(
    packet: RU10Packet,
    dna_version_string: str,
    insertion_position: int = 30,
    diff: typing.Optional[np.ndarray] = None,
) -> None:
    """
    Add a DNA version string at a specific position in an RU10Packet.

    This function inserts a version string into the packet data at the specified
    position, handling padding and seed spacing as needed.

    Args:
        packet: The RU10Packet to modify.
        dna_version_string: The DNA string to insert (e.g., "ACGT").
        insertion_position: Byte position in the data section (not DNA representation).
        diff: Array containing the diff between original and new packet data.
            Used for asserting correct insertion position. If None, skips validation.

    Raises:
        AssertionError: If insertion position overlaps with non-padding region.

    Note:
        The function handles padding automatically if the version string length
        is not a multiple of 4 (since each byte encodes to 4 DNA bases).

    Example:
        >>> insert_dna_version_string(
        ...     packet,
        ...     "ACGTACGT",
        ...     insertion_position=50,
        ...     diff=chunk_diff
        ... )
    """
    # Ensure packet has DNA data
    if packet.dna_data is None or packet.dna_data == "":
        packet.calculate_packed_data()

    # Validate insertion position if diff provided
    if diff is not None:
        insert_len = int((len(dna_version_string) + 1) / 2) + 1
        assert not diff[insertion_position : insertion_position + insert_len].any(), (
            "Cannot insert version string into non-padding region!"
        )

    # Calculate padding needed to align to byte boundary
    pre_padding = "A" * ((4 - len(dna_version_string) % 4) % 4)
    binary_dna_version_string = tranlate_quat_to_byte(pre_padding + dna_version_string)
    num_padding_bits = len(pre_padding) * 2

    # Handle first byte if padding was added
    if num_padding_bits > 0:
        original_byte = packet.data[insertion_position]
        mask_keep = (0xFF << (8 - num_padding_bits)) & 0xFF
        mask_new = (0xFF >> num_padding_bits) & 0xFF
        combined_first_byte = (original_byte & mask_keep) | (
            binary_dna_version_string[0] & mask_new
        )
        binary_dna_version_string = (
            bytes([combined_first_byte]) + binary_dna_version_string[1:]
        )

    # Create stub array with version string at insertion position
    stub = np.zeros_like(packet.data)
    end_pos = insertion_position + len(binary_dna_version_string)
    stub[insertion_position:end_pos] = np.frombuffer(
        binary_dna_version_string, dtype=np.uint8
    )

    # Apply seed-based XOR if enabled
    stub = xor_with_seed(stub, packet.id)

    # Update packet data
    packet.data[insertion_position:end_pos] = stub[insertion_position:end_pos]

    # Recalculate DNA structure
    packet.get_dna_struct(True, packet.id_spacing, packet.id_spacing_length, True)

    logger.debug(
        f"Inserted version string at position {insertion_position} in packet {packet.id}"
    )


def reduce_packet_to_chunk(
    packet: RU10Packet,
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    chunk_to_reduce_to: int = 0,
) -> RU10Packet:
    """
    Reduce a packet to contain only a single target chunk.

    This function removes all chunks from a packet except the specified target
    chunk by XORing with the known chunk data from the decoder state.

    Args:
        packet: RU10Packet to reduce.
        semiautomatic_solver: SemiAutomaticReconstructionToolkit instance.
        chunk_to_reduce_to: Chunk index to keep (default: 0 for header chunk).

    Returns:
        Modified packet containing only the target chunk.

    Raises:
        AssertionError: If packet doesn't contain exactly the target chunk after reduction.

    Note:
        This operation is essential for isolating individual chunks during
        version decoding and repair operations.

    Example:
        >>> reduced = reduce_packet_to_chunk(packet, solver, chunk_to_reduce_to=0)
        >>> assert len(reduced.used_packets) == 1
    """
    # Get normalized used chunks and remove auxiliary packets
    used_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(packet)
    packet.used_packets = {i for i, x in enumerate(used_chunks) if x}

    # XOR out all chunks except the target
    for i, must_remove_chunk in enumerate(used_chunks):
        if must_remove_chunk and i != chunk_to_reduce_to:
            packet.data = xor_numpy(packet.data, semiautomatic_solver.decoder.GEPP.b[i])
            packet.used_packets.remove(i)

    # Verify only target chunk remains
    assert len(packet.used_packets) == 1 and chunk_to_reduce_to in packet.used_packets, (
        f"Packet reduction failed: expected only chunk {chunk_to_reduce_to}, "
        f"got {packet.used_packets}"
    )

    return packet


def insert_id_string(
    packet: RU10Packet,
    insertion_position: int,
    changed_chunk_id: int,
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
) -> None:
    """
    Insert the ID of a changed chunk into a packet.

    This function encodes the index of a changed chunk into the packet data,
    allowing the decoder to identify which chunk was modified.

    Args:
        packet: RU10Packet to modify.
        insertion_position: Byte position for ID insertion.
        changed_chunk_id: ID of the changed chunk to encode.
        semiautomatic_solver: SemiAutomaticReconstructionToolkit instance.

    Raises:
        AssertionError: If changed_chunk_id is not in used chunks or >= 256.

    Note:
        The chunk ID is encoded as a single byte (0-255 range). For larger
        chunk indices, this function will need to be extended.

    Example:
        >>> insert_id_string(packet, insertion_position=100, changed_chunk_id=5, solver=solver)
    """
    # Get used chunks list
    used_chunks = from_true_false_list(
        semiautomatic_solver.decoder.removeAndXorAuxPackets(packet)
    )

    # Validate chunk ID is in used chunks
    assert changed_chunk_id in used_chunks, (
        f"Changed chunk id {changed_chunk_id} must be in used chunks {used_chunks}!"
    )

    # Get index of chunk in used chunks list
    index_of_chunk = used_chunks.index(changed_chunk_id)

    # Validate chunk index fits in single byte
    assert index_of_chunk < 256, (
        f"Cannot encode chunk id index >= 256 in this version! Got {index_of_chunk}"
    )

    # Create XOR mask with chunk ID at insertion position
    to_xor = np.zeros_like(packet.data, dtype=np.uint8)
    to_xor[insertion_position] = index_of_chunk & 0xFF

    # Apply XOR to packet data
    packet.data = xor_numpy(packet.data, to_xor)

    # Recalculate DNA structure
    packet.get_dna_struct(True, packet.id_spacing, packet.id_spacing_length, True)

    logger.debug(
        f"Inserted chunk ID {changed_chunk_id} (index {index_of_chunk}) at position {insertion_position}"
    )


def find_insertion_position(
    version_string_length: int,
    chunk_diff: np.ndarray,
    id_spacing: int,
    id_len: int,
) -> Generator[int, int, int]:
    """
    Find a suitable insertion position for the version string in the chunk diff.

    Scans the chunk diff array for regions of zeros (padding) that are large
    enough to accommodate the version string without overlapping seed positions.

    Args:
        version_string_length: Length of version string in bytes (not quats).
        chunk_diff: Numpy array representing the chunk diff.
        id_spacing: Spacing between chunk ID DNA bases in the packet.
        id_len: Length of the chunk ID in bytes.

    Yields:
        Starting indices of valid insertion positions.

    Returns:
        -1 if no suitable position is found.

    Note:
        The function searches from the end of the array backwards to find
        the last suitable position, which often provides better stability.

    Example:
        >>> positions = list(find_insertion_position(10, chunk_diff, id_spacing=2, id_len=4))
        >>> if positions:
        ...     print(f"Found {len(positions)} insertion positions")
    """
    # Find all zero regions in chunk diff
    zero_regions = np.where(chunk_diff == 0)[0]

    # Search for suitable positions from end backwards
    for start_idx in reversed([x for x in zero_regions]):
        # Check if region is large enough and doesn't overlap with ID spacing
        if (
            start_idx + version_string_length <= len(chunk_diff)
            and id_spacing * id_len * 4 - 1 < start_idx
        ):
            # Verify entire region is zero
            if np.all(chunk_diff[start_idx : start_idx + version_string_length] == 0):
                yield start_idx

    # No suitable position found
    return -1


def find_insertion_position_with_seed(
    version_string_length: int,
    chunk_diff: np.ndarray,
    id_spacing: int,
    id_len: int,
) -> Generator[int, int, int]:
    """
    Find insertion position accounting for seed spacing constraints.

    Similar to find_insertion_position but creates a detailed mask that accounts
    for seed spacing at the quaternary (DNA base) level, ensuring the version
    string doesn't overlap with seed positions.

    Args:
        version_string_length: Length of version string in bytes (not quats).
        chunk_diff: Numpy array representing the chunk diff.
        id_spacing: Spacing between chunk ID DNA bases in the packet.
        id_len: Length of the chunk ID in bytes.

    Yields:
        Starting indices of valid insertion positions that don't overlap with seeds.

    Returns:
        -1 if no suitable position is found.

    Note:
        This function is more precise than find_insertion_position when seed
        spacing is used, as it checks at the DNA base level rather than byte level.

    Example:
        >>> positions = list(
        ...     find_insertion_position_with_seed(10, chunk_diff, id_spacing=2, id_len=4)
        ... )
    """
    # Create mask for valid positions (True = usable, False = blocked by seed)
    num_bytes = len(chunk_diff)
    num_quats = num_bytes * 4  # Each byte = 4 DNA bases (quaternary)

    # Create quaternary-level mask (True = usable)
    quat_mask = np.ones(num_quats, dtype=bool)

    if id_spacing > 0 and id_len > 0:
        # Number of seed bases (each byte = 4 DNA bases)
        num_seed_bases = id_len * 4
        # Seed bases are placed at positions 0, (id_spacing+1), 2*(id_spacing+1), ...
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

    # Find zero regions in valid positions
    zero_regions = np.where(valid_positions)[0]

    # Search for suitable positions from end backwards
    for start_idx in reversed([x for x in zero_regions]):
        if start_idx + version_string_length <= len(chunk_diff):
            # Verify entire region is valid
            if np.all(valid_positions[start_idx : start_idx + version_string_length]):
                yield start_idx

    # No suitable position found
    return -1


def generate_new_packets(
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    encoder: RU10Encoder,
    diff: np.ndarray,
    changed_chunk_ids: np.ndarray,
    new_file_version: int,
) -> Dict[int, Union[List[RU10Packet], List[Tuple[RU10Packet, RU10Packet]]]]:
    """
    Scan ALL seeds and record ALL candidate seeds for each changed chunk.

    This function performs an exhaustive search through all possible seeds to
    find packets that:
    1. Include at least one changed chunk
    2. Include the header chunk (0) or last chunk
    3. Have suitable insertion positions for version strings

    Args:
        semiautomatic_solver: SemiAutomaticReconstructionToolkit instance.
        encoder: RU10Encoder instance for packet generation.
        diff: Diff array from find_affected_chunks.
        changed_chunk_ids: Array of changed chunk indices.
        new_file_version: New version number to encode.

    Returns:
        Dictionary mapping chunk IDs to lists of generated packets (or packet pairs).

    Raises:
        RuntimeError: If no valid seeds found for any changed chunk.

    Optimizations:
        - Precompute changed set for O(1) membership checks
        - Quick-skip seeds that don't mention any changed chunk
        - Quick-skip seeds that don't include header or last chunk

    Example:
        >>> packets = generate_new_packets(solver, encoder, diff, changed_chunks, version)
        >>> for chunk_id, chunk_packets in packets.items():
        ...     print(f"Chunk {chunk_id}: {len(chunk_packets)} packets")
    """
    res: Dict[
        int, Union[List[RU10Packet], List[Tuple[RU10Packet, RU10Packet]]]
    ] = {}

    max_num = min(
        encoder.calc_max_size(struct.calcsize("<" + encoder.id_len_format)), int31
    )

    last_chunk_idx = encoder.number_of_chunks - 1
    changed_set = {int(x) for x in changed_chunk_ids}

    # Initialize mapping for results - one entry per changed chunk
    packet_to_seed_mapping: Dict[int, Set[int]] = {cid: set() for cid in changed_set}

    logger.info(f"Scanning {max_num} seeds for {len(changed_set)} changed chunks...")

    # Scan all seeds
    for seed in range(max_num):
        packet_numbers = choose_packet_numbers(
            encoder.number_of_chunks, seed, encoder.dist, systematic=False
        )

        # OPTIMIZATION 1: Quick filter using set intersection
        pn_set = set(packet_numbers)
        if not (pn_set & changed_set):
            continue

        # OPTIMIZATION 2: Skip if doesn't include header or last chunk
        if 0 not in pn_set and last_chunk_idx not in pn_set:
            continue

        # Create packet and get used chunks
        packet = RU10Packet(
            b"", packet_numbers, encoder.number_of_chunks, seed, encoder.dist, read_only=True
        )
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

        # Log progress every 10000 seeds
        if seed % 10000 == 0 and seed > 0:
            logger.info(
                f"Scanned {seed}/{max_num} seeds... "
                f"Found {sum(len(v) for v in packet_to_seed_mapping.values())} total candidates"
            )

    # Verify we found at least one seed for every changed chunk
    missing = [cid for cid, seeds in packet_to_seed_mapping.items() if not seeds]
    if missing:
        raise RuntimeError(
            f"Could not find ANY seed for chunks {missing} that contains either header or last chunk"
        )

    logger.info(
        f"Scan complete! Total candidates found: "
        f"{sum(len(v) for v in packet_to_seed_mapping.values())}"
    )

    # Select best seeds for each chunk
    chunk_to_potential_seed_mapping = select_numbers(packet_to_seed_mapping, n=10)

    # Generate packets for selected seeds
    generated_packets: Dict[int, List[RU10Packet]] = {}
    for chunk_id, seed_set in chunk_to_potential_seed_mapping:
        for seed in seed_set:
            packet = encoder.create_new_packet(False, seed)
            should_drop_packet(
                encoder.rules, packet, 1.0
            )  # Calculate error probability

            if chunk_id not in generated_packets:
                generated_packets[chunk_id] = []
            generated_packets[chunk_id].append(packet)

        # Sort by error probability (best first)
        generated_packets[chunk_id] = sorted(
            generated_packets[chunk_id], key=lambda x: x.error_prob
        )

    # Populate header chunk for version string insertion
    semiautomatic_solver.decoder.populate_header_chunk(
        last_chunk_len_str=semiautomatic_solver.decoder.config_map.get(
            "last_chunk_len_str", "I"
        )
    )

    # Process each changed chunk
    changed_chunk_to_new_packets: Dict[int, List[RU10Packet]] = {}
    changed_chunk_to_packet_pair_list: Dict[
        int, List[Tuple[RU10Packet, RU10Packet]]
    ] = {}

    version_string = generate_dna_version_string(new_file_version)

    for changed_chunk, potential_packets in generated_packets.items():
        packet_added = False

        # Try each potential packet
        for potential_packet in potential_packets:
            modified_packet = potential_packet.copy()
            plain_used_chunks = semiautomatic_solver.decoder.removeAndXorAuxPackets(
                modified_packet
            )

            if plain_used_chunks[0]:
                # Header chunk case - find insertion position
                try:
                    insertion_position = next(
                        find_insertion_position(
                            int((len(version_string) + 1) / 4 + 1),
                            diff[changed_chunk],
                            modified_packet.id_spacing,
                            struct.calcsize(modified_packet.id_len_format),
                        )
                    )
                except StopIteration:
                    # No suitable position found, try next packet
                    continue

                # Insert version string
                insert_dna_version_string(
                    modified_packet, version_string, insertion_position, diff[changed_chunk]
                )

                # Insert chunk ID
                insert_id_string(
                    modified_packet,
                    int(insertion_position + (len(version_string) + 1) / 4),
                    changed_chunk,
                    semiautomatic_solver,
                )

            elif plain_used_chunks[-1]:
                # Last chunk case - not yet implemented
                continue
            else:
                # Neither header nor last chunk - skip
                continue

            # Apply diff for changed chunk
            modified_packet.data = xor_numpy(modified_packet.data, diff[changed_chunk])

            # Recalculate DNA structure
            modified_packet.get_dna_struct(
                True, modified_packet.id_spacing, modified_packet.id_spacing_length, True
            )

            # Calculate error probability
            should_drop_packet(encoder.rules, modified_packet)

            # Store packet
            if changed_chunk not in changed_chunk_to_new_packets:
                changed_chunk_to_new_packets[changed_chunk] = []
            changed_chunk_to_new_packets[changed_chunk].append(modified_packet)
            packet_added = True

        # Handle case where no single packet worked - try packet pairs
        if not packet_added or (
            changed_chunk_to_new_packets.get(changed_chunk)
            and sorted(
                changed_chunk_to_new_packets[changed_chunk], key=lambda x: x.error_prob
            )[0].error_prob
            >= 1.0
        ):
            logger.warning(
                f"Could not find suitable packet for chunk {changed_chunk} to insert "
                f"version string. Generating two packets with same seed with split diff (50/50)!"
            )

            for potential_packet in potential_packets:
                modified_packet_first = potential_packet.copy()
                modified_packet_second = potential_packet.copy()

                # Split diff in half
                diff_mask = np.zeros_like(diff[changed_chunk], dtype=bool)
                half_point = len(diff[changed_chunk]) // 2
                diff_mask[:half_point] = True
                first_diff = np.where(diff_mask, diff[changed_chunk], 0)
                second_diff = np.where(~diff_mask, diff[changed_chunk], 0)

                # Try to find insertion positions for both halves
                if plain_used_chunks[0]:
                    try:
                        insertion_position_first = next(
                            find_insertion_position(
                                int((len(version_string) + 1) / 4 + 1),
                                diff[changed_chunk],
                                modified_packet_first.id_spacing,
                                struct.calcsize(modified_packet_first.id_len_format),
                            )
                        )
                        insertion_position_second = next(
                            find_insertion_position(
                                int((len(version_string) + 1) / 4 + 1),
                                diff[changed_chunk],
                                modified_packet_second.id_spacing,
                                struct.calcsize(modified_packet_second.id_len_format),
                            )
                        )
                    except StopIteration:
                        continue

                    # Insert version string and chunk ID in both packets
                    insert_dna_version_string(
                        modified_packet_first,
                        version_string,
                        insertion_position_first,
                        first_diff,
                    )
                    insert_id_string(
                        modified_packet_second,
                        int(insertion_position_second + (len(version_string) + 1) / 4),
                        changed_chunk,
                        semiautomatic_solver,
                    )

                    insert_dna_version_string(
                        modified_packet_first,
                        version_string,
                        insertion_position_first,
                        second_diff,
                    )
                    insert_id_string(
                        modified_packet_second,
                        int(insertion_position_second + (len(version_string) + 1) / 4),
                        changed_chunk,
                        semiautomatic_solver,
                    )

                elif plain_used_chunks[-1]:
                    # Last chunk case - not yet implemented
                    continue

                # Apply split diffs
                modified_packet_first.data = xor_numpy(
                    modified_packet_first.data, first_diff
                )
                modified_packet_first.get_dna_struct(
                    True,
                    modified_packet_first.id_spacing,
                    modified_packet_first.id_spacing_length,
                    True,
                )

                modified_packet_second.data = xor_numpy(
                    modified_packet_second.data, second_diff
                )
                modified_packet_second.get_dna_struct(
                    True,
                    modified_packet_second.id_spacing,
                    modified_packet_second.id_spacing_length,
                    True,
                )

                # Calculate error probabilities
                should_drop_packet(encoder.rules, modified_packet_first)
                should_drop_packet(encoder.rules, modified_packet_second)

                # Store packet pair
                changed_chunk_to_new_packets[changed_chunk].append(modified_packet_first)
                changed_chunk_to_new_packets[changed_chunk].append(modified_packet_second)

                if changed_chunk not in changed_chunk_to_packet_pair_list:
                    changed_chunk_to_packet_pair_list[changed_chunk] = []
                changed_chunk_to_packet_pair_list[changed_chunk].append(
                    (modified_packet_first, modified_packet_second)
                )
                packet_added = True

        if not packet_added:
            logger.error(f"Could not create any packet for changed chunk {changed_chunk}!")

    # Select best packets for each changed chunk
    for changed_chunk in changed_chunk_to_new_packets.keys():
        if changed_chunk in changed_chunk_to_new_packets:
            # Try using best single packet first
            changed_chunk_to_new_packets[changed_chunk] = sorted(
                changed_chunk_to_new_packets[changed_chunk], key=lambda x: x.error_prob
            )

            for new_pack in changed_chunk_to_new_packets[changed_chunk]:
                if new_pack.error_prob < 1.0:
                    if changed_chunk not in res:
                        res[changed_chunk] = []
                    res[changed_chunk].append(new_pack)
                    logger.info(
                        f"Generated packet with error probability {new_pack.error_prob} "
                        f"for chunk {changed_chunk}."
                    )
                else:
                    logger.warning(
                        f"Skipping packet with error probability {new_pack.error_prob} "
                        f"for chunk {changed_chunk}!"
                    )

        elif changed_chunk in changed_chunk_to_packet_pair_list:
            # Use best packet pair
            changed_chunk_to_packet_pair_list[changed_chunk] = sorted(
                changed_chunk_to_packet_pair_list[changed_chunk],
                key=lambda x: max(x[0].error_prob, x[1].error_prob),
            )

            for new_pack_pair in changed_chunk_to_packet_pair_list[changed_chunk]:
                if (
                    new_pack_pair[0].error_prob < 1.0
                    and new_pack_pair[1].error_prob < 1.0
                ):
                    if changed_chunk not in res:
                        res[changed_chunk] = []
                    res[changed_chunk].append(new_pack_pair)
                    logger.info(
                        f"Generated packet pair with error probability "
                        f"{(new_pack_pair[0].error_prob, new_pack_pair[1].error_prob)} "
                        f"for chunk {changed_chunk}!"
                    )
                else:
                    logger.warning(
                        f"Skipping packet pair with error probability "
                        f"{(new_pack_pair[0].error_prob, new_pack_pair[1].error_prob)} "
                        f"for chunk {changed_chunk}!"
                    )

    return res


def decode_versions(
    semiautomatic_solver: SemiAutomaticReconstructionToolkit,
    dna_version_string_prefix: str = "",
) -> typing.Dict[int, str]:
    """
    Decode all versions of an encoded file from the packets known to the solver.

    This function iterates through all available versions in the DNA pool and
    decodes each one, saving them with version-prefixed filenames.

    Args:
        semiautomatic_solver: The SemiAutomaticReconstructionToolkit with the
            decoder containing the packets.
        dna_version_string_prefix: The prefix used for version strings in the
            DNA sequences (default: "").

    Returns:
        A dictionary mapping file version numbers to filenames for all versions found.

    Note:
        The decoding process:
        1. Decodes base version using only packets WITHOUT version strings
        2. For each version: decodes changed chunks using packets WITH that version's string
        3. Saves each version as v<version>_<base_filename>.<extension>
        4. Uses each new version as the base for the next version

    Example:
        >>> versions = decode_versions(solver)
        >>> for version, filename in versions.items():
        ...     print(f"Version {version}: {filename}")
    """
    res: typing.Dict[int, str] = {}

    # Get base version string
    if not dna_version_string_prefix:
        dna_version_string_prefix = "GAGCCAGTGAGTCGTA"

    # Get maximum version in pool
    from MultiVersionDecoder import MultiVersionDecoder

    mv_decoder = MultiVersionDecoder(semiautomatic_solver.decoder)
    max_version = mv_decoder.get_versions_in_pool(dna_version_string_prefix)

    logger.info(f"Found versions 0 to {max_version} in pool")

    # Decode base version (version 0)
    logger.info("Decoding base version (v0)...")
    mv_decoder.decode_base_version(dna_version_string_prefix)
    res[0] = f"v0_{semiautomatic_solver.decoder.headerChunk.file_name.decode('utf-8')}"

    # Decode each subsequent version
    for version in range(1, max_version + 1):
        logger.info(f"Decoding version {version}...")
        try:
            decoded = mv_decoder.decode_to_version(dna_version_string_prefix, version)
            if version in decoded:
                res[version] = decoded[version]
                logger.info(f"Version {version} decoded successfully")
        except Exception as e:
            logger.error(f"Failed to decode version {version}: {e}")
            # Continue with next version

    return res


def add_packets(
    encoder: RU10Encoder, new_packets: typing.Union[Set[RU10Packet], List[RU10Packet]]
) -> None:
    """
    Add new packets to an encoder's packet set.

    This is a convenience function for adding generated packets to an encoder.

    Args:
        encoder: RU10Encoder instance to add packets to.
        new_packets: Set or list of RU10Packet objects to add.

    Example:
        >>> add_packets(encoder, generated_packets)
    """
    if isinstance(new_packets, list):
        encoder.encodedPackets.update(new_packets)
    else:
        encoder.encodedPackets |= new_packets


def main() -> None:
    """
    Main entry point for file_update_coding CLI.

    This function parses command-line arguments and performs file update encoding.
    """
    parsed_args = init_args()

    ini_file = parsed_args.ini
    new_file_path = parsed_args.new_file
    packet_add_limit = parsed_args.packet_add_limit

    # Load configuration
    cfg_worker = ConfigReadAndExecute(ini_file)
    x = cfg_worker.execute(return_decoder=True, skip_solve=True)[0]
    semiautomatic_solver = SemiAutomaticReconstructionToolkit(x)

    # Read new file content
    if not Path(new_file_path).exists():
        logger.error(f"New file not found: {new_file_path}")
        return
    
    with open(new_file_path, "rb") as f:
        new_file_content = f.read()
    
    logger.info(f"Loaded new file: {new_file_path} ({len(new_file_content)} bytes)")

    # Calculate diff and find changed chunks
    diff, changed_chunks = find_affected_chunks(semiautomatic_solver, new_file_content)

    # Get current version and increment
    new_file_version = get_current_file_version(semiautomatic_solver) + 1
    logger.info(f"Creating version {new_file_version}")

    # Create encoder
    encoder = encoder_from_decoder(
        semiautomatic_solver, cfg_worker, rules=FastDNARules()
    )

    # Generate packets for changed content
    packet_candidates = generate_new_packets(
        semiautomatic_solver, encoder, diff, changed_chunks, new_file_version
    )

    logger.info(
        f"Found candidate seeds per chunk: "
        f"{dict((k, len(v)) for k, v in packet_candidates.items())}"
    )

    # Add best packets to encoder
    packets_added = 0
    added_packets = []

    for changed_chunk_packet_group, packets in packet_candidates.items():
        added_packets_for_chunk = 0
        for packet in packets:
            if packet.error_prob >= 1.0:
                continue
            if added_packets_for_chunk >= packet_add_limit:
                break

            encoder.encodedPackets.add(packet)
            added_packets.append(packet)
            added_packets_for_chunk += 1
            packets_added += 1

        logger.info(
            f"Added {added_packets_for_chunk} packets for changed chunk "
            f"{changed_chunk_packet_group}."
        )

    logger.info(
        f"Added total of {packets_added} packets for {len(packet_candidates)} changed chunks."
    )

    # Save updated FASTA and config
    outfile = f"{semiautomatic_solver.decoder.file.split('.fasta')[0]}_v{new_file_version}"
    file_bkp = encoder.file
    encoder.file = outfile
    encoder.save_packets_fasta(None, "", False)
    encoder.save_config_file(add_dot_fasta=True)

    # Save added packets to debug file
    debug_outfile = (
        f"{semiautomatic_solver.decoder.file.split('.fasta')[0]}"
        f"_v{new_file_version}_added_packets.fasta"
    )
    with open(debug_outfile, "w") as f:
        for packet in added_packets:
            f.write(f">{packet.id}\n")
            f.write(f"{packet.dna_data}\n")

    encoder.file = file_bkp
    logger.info(f"Version {new_file_version} saved successfully")


if __name__ == "__main__":
    main()
