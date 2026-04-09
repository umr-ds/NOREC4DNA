r"""
Multi-Version Decoder for NOREC4DNA Encoded Files.

This module provides functionality to decode files encoded with NOREC4DNA that contain
multiple versions. It supports:

- Decoding files with multiple version packets
- Handling missing chunks with partial reconstruction
- Ranking missing chunks by importance
- Manual and automatic error correction
- Version-aware decoding with metadata filtering

Key Features:
    1. Version-Aware Decoding: Filter and decode specific versions from a pool of DNA sequences
    2. Metadata Filtering: Exclude metadata sequences during decoding
    3. Error Correction: Automatic and manual repair of corrupt packets
    4. Partial Reconstruction: Decode files even with missing chunks
    5. Checksum Validation: Verify decoded files using embedded checksums

Example Usage:
    >>> from MultiVersionDecoder import MultiVersionDecoder
    >>> from NOREC4DNA.ConfigWorker import ConfigReadAndExecute
    >>>
    >>> # Load decoder configuration
    >>> config = ConfigReadAndExecute("config.ini")
    >>> decoder = config.execute(return_decoder=True)[0]
    >>>
    >>> # Initialize MultiVersionDecoder
    >>> mv_decoder = MultiVersionDecoder(
    ...     decoder,
    ...     metadata_list=["ACGT", "TGCA"]  # Optional: metadata sequences to filter
    ... )
    >>>
    >>> # Decode base version (version 0)
    >>> base_version_string = "GAGCCAGTGAGTCGTA"
    >>> mv_decoder.decode_base_version(base_version_string)
    >>>
    >>> # Decode to version N
    >>> mv_decoder.decode_to_version(base_version_string, version=2)
"""

import argparse
import logging
import shutil
import struct
import typing
from functools import reduce
from io import BytesIO
from itertools import combinations
from pathlib import Path

import numpy as np

from NOREC4DNA.ConfigWorker import ConfigReadAndExecute
from MultiVersionCoder import reduce_packet_to_chunk
from NOREC4DNA.invivo_window_decoder import load_fasta
from NOREC4DNA.norec4dna.HeaderChunk import HeaderChunk
from NOREC4DNA.norec4dna.helper.helper_cpu_single_core import xor_numpy
from NOREC4DNA.norec4dna.helper.quaternary2Bin import tranlate_quat_to_byte
from NOREC4DNA.norec4dna.helper.RU10Helper import from_true_false_list
from NOREC4DNA.norec4dna.LTDecoder import LTDecoder
from NOREC4DNA.norec4dna.OnlineDecoder import OnlineDecoder
from NOREC4DNA.norec4dna.RU10Decoder import RU10Decoder
from NOREC4DNA.norec4dna.RU10Packet import RU10Packet
from semi_automatic_reconstruction_toolkit import SemiAutomaticReconstructionToolkit

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class MultiVersionDecoder(SemiAutomaticReconstructionToolkit):
    """
    Multi-version decoder for NOREC4DNA encoded files.

    This class extends SemiAutomaticReconstructionToolkit to provide version-aware
    decoding capabilities for DNA-encoded files containing multiple versions.

    Attributes:
        decoder: The underlying decoder instance (RU10Decoder, LTDecoder, or OnlineDecoder)
        metadata_list: List of DNA sequences to filter out as metadata
        headerChunk: Parsed header chunk from the encoded file
        multi_error_packets_mode: Flag for handling multiple error packets
    """

    def __init__(
        self,
        decoder: typing.Union[RU10Decoder, LTDecoder, OnlineDecoder],
        metadata_list: typing.Optional[typing.List[str]] = None,
    ) -> None:
        """
        Initialize the MultiVersionDecoder with a decoder instance.

        Args:
            decoder: The decoder instance to use for reconstruction.
                Can be RU10Decoder, LTDecoder, or OnlineDecoder.
            metadata_list: Optional list of metadata DNA sequences to filter out
                during decoding. If None, an empty list is used.

        Raises:
            ValueError: If decoder is None or invalid.
        """
        if decoder is None:
            raise ValueError("Decoder instance cannot be None")

        super().__init__(decoder)
        self.last_chunk_len_format = "I"
        self.checksum_len_format: typing.Optional[str] = None
        self.decoder: typing.Union[RU10Decoder, LTDecoder, OnlineDecoder] = decoder
        decoder.read_all_before_decode = True
        self.headerChunk: typing.Optional[HeaderChunk] = None
        self.decoder.GEPP.insert_tmp()
        self.initial_A = self.decoder.GEPP.A.copy()
        self.initial_b = self.decoder.GEPP.b.copy()
        self.initial_packet_mapping: typing.Optional[dict] = None
        self.multi_error_packets_mode = False

        if metadata_list is None:
            self.metadata_list: typing.List[str] = []
        else:
            self.metadata_list = metadata_list

        # Initialize version tracking
        self._version_cache: typing.Dict[int, typing.List[str]] = {}

    def get_versions_in_pool(self, base_dna_version_string: str) -> int:
        """
        Return the largest version number available in the pool.

        Scans all DNA sequences in the pool and extracts version numbers from
        sequences containing the base_dna_version_string marker.

        Args:
            base_dna_version_string: The magic DNA string marking version sequences
                (e.g., "GAGCCAGTGAGTCGTA").

        Returns:
            The highest version number found. Returns 0 if no versions are found
            (base version only).

        Note:
            Versions are indexed starting from 0, where version 0 is the base version
            and version 1 is the FIRST version after the base version.

        Example:
            >>> mv_decoder = MultiVersionDecoder(decoder)
            >>> max_version = mv_decoder.get_versions_in_pool("GAGCCAGTGAGTCGTA")
            >>> print(f"Available up to version: {max_version}")
        """
        res = 0
        fasta_entries = load_fasta(self.decoder.file)

        for seq in fasta_entries.values():
            idx = seq.find(base_dna_version_string)
            if idx != -1 and idx + len(base_dna_version_string) < len(seq):
                # Extract 3 bases before the magic string (version encoding)
                version_bases = seq[max(0, idx - 3) : idx]
                # Pad if necessary
                if len(version_bases) < 3:
                    version_bases = "A" * (3 - len(version_bases)) + version_bases
                version_num = tranlate_quat_to_byte(f"A{version_bases}")
                try:
                    version_value = struct.unpack("B", version_num)[0]
                    res = max(res, version_value)
                except struct.error:
                    logger.warning(f"Failed to parse version from sequence: {seq[:50]}...")

        return res

    def get_sequences_for_version(
        self, base_dna_version_string: str, version: int
    ) -> typing.List[str]:
        """
        Return a list of all sequences corresponding to the given version.

        Scans the DNA pool and returns all sequences that contain the specified
        version number encoded before the base_dna_version_string marker.

        Args:
            base_dna_version_string: The base DNA version string marker.
            version: The version number to retrieve (0-based, where 0 is base version).

        Returns:
            A list of DNA sequences for the specified version.
            Returns an empty list if the version is not found.

        Example:
            >>> sequences = mv_decoder.get_sequences_for_version("GAGCCAGTGAGTCGTA", version=1)
            >>> print(f"Found {len(sequences)} sequences for version 1")
        """
        res: typing.List[str] = []
        fasta_entries = load_fasta(self.decoder.file)

        for seq in fasta_entries.values():
            idx = seq.find(base_dna_version_string)
            if idx != -1 and idx + len(base_dna_version_string) < len(seq):
                # Extract 3 bases before the magic string
                version_bases = seq[max(0, idx - 3) : idx]
                # Pad if necessary
                if len(version_bases) < 3:
                    version_bases = "A" * (3 - len(version_bases)) + version_bases
                version_num = tranlate_quat_to_byte(f"A{version_bases}")
                try:
                    version_value = struct.unpack("B", version_num)[0]
                    if version_value == version:
                        res.append(seq)
                except struct.error:
                    logger.warning(f"Failed to parse version from sequence: {seq[:50]}...")

        return res

    def contains_metadata(
        self, seq: str, metadata_list: typing.Optional[typing.List[str]] = None
    ) -> bool:
        """
        Check if a DNA sequence contains any metadata sequences.

        Args:
            seq: DNA sequence to check.
            metadata_list: Optional list of metadata sequences to check against.
                If None, uses the instance's metadata_list.

        Returns:
            True if the sequence contains any metadata, False otherwise.

        Example:
            >>> if mv_decoder.contains_metadata(seq):
            ...     print("Sequence contains metadata - filtering out")
        """
        if metadata_list is None:
            metadata_list = self.metadata_list

        if not metadata_list:
            return False

        return any(metadata in seq for metadata in metadata_list)

    def decode_base_version(
        self, base_dna_version_string: str, known_base_file: typing.Optional[str] = None
    ) -> typing.Union[RU10Decoder, LTDecoder, OnlineDecoder]:
        """
        Decode the base version (version 0) and store it on disk.

        This method decodes the original file without any version modifications.
        It filters out packets containing version strings and metadata sequences.

        Args:
            base_dna_version_string: DNA version string to filter out during decoding.
            known_base_file: Optional path to a pre-decoded base file. If provided
                and exists, loads the base version from this file instead of decoding.

        Returns:
            The decoder instance after successful decoding.

        Raises:
            FileNotFoundError: If known_base_file is provided but doesn't exist.
            RuntimeError: If decoding fails to complete.

        Note:
            The decoded file is saved with a "v0_" prefix to distinguish it from
            other versions.

        Example:
            >>> # Decode from DNA sequences
            >>> mv_decoder.decode_base_version("GAGCCAGTGAGTCGTA")
            >>>
            >>> # Or load from existing file
            >>> mv_decoder.decode_base_version(
            ...     "GAGCCAGTGAGTCGTA",
            ...     known_base_file="v0_original_file.txt"
            ... )
        """
        # Check if we can load from a known base file
        if known_base_file is not None and Path(known_base_file).exists():
            logger.info(f"Base version already decoded, loading from {known_base_file}")
            try:
                with open(known_base_file, "rb") as f:
                    data = f.read()

                chunk_size = self.decoder.GEPP.chunk_size
                # Split data into chunks
                res = [
                    np.frombuffer(data[i : i + chunk_size], dtype=np.uint8)
                    for i in range(0, len(data), chunk_size)
                ]

                # Create identity matrix and update GEPP state
                iden = np.identity(self.decoder.number_of_chunks)
                self.decoder.GEPP.A = iden
                self.decoder.GEPP.b = np.array(res, dtype=np.uint8)
                self.decoder.input_new_packet()

                logger.info("Base version loaded successfully from file")
                return self.decoder
            except Exception as e:
                logger.error(f"Failed to load base version from file: {e}")
                # Fall through to decoding from DNA

        # Decode base version from DNA sequences
        logger.info("Decoding base version from DNA sequences...")

        # Reset decoder state
        self.decoder = type(self.decoder).from_config_map(self.decoder.config_map)

        # Load all FASTA entries
        fasta_entries = load_fasta(self.decoder.file)

        # Filter out sequences containing metadata or version strings
        fasta_seqs = [
            seq
            for seq in fasta_entries.values()
            if not self.contains_metadata(seq, self.metadata_list)
            and not self.contains_metadata(seq, [base_dna_version_string])
        ]

        logger.info(f"Found {len(fasta_seqs)} sequences for base version decoding")

        # Get configuration parameters
        id_len_format = self.decoder.config_map.get("id_len_format", "")
        crc_len_format = self.decoder.config_map.get("crc_len_format", "")
        packet_len_format = self.decoder.config_map.get("packet_len_format", "")

        # Process packets until we can solve
        for seq in fasta_seqs:
            # Revert seed spacing as it is DNA-based and not part of parse_raw_packet
            seq = self.decoder.revert_seed_spacing(seq, id_len_format)

            packet = self.decoder.parse_raw_packet(
                BytesIO(tranlate_quat_to_byte(seq)).read(),
                crc_len_format=crc_len_format,
                number_of_chunks_len_format="",
                packet_len_format=packet_len_format,
                id_len_format=id_len_format,
            )

            self.decoder.input_new_packet(packet)
            self.decoder.packets.append(packet)

            if len(self.decoder.packets) >= self.decoder.static_number_of_chunks:
                if self.decoder.solve():
                    logger.info("Base version decoded successfully")
                    break

        # Populate header chunk if enabled
        if self.decoder.use_headerchunk:
            self.decoder.populate_header_chunk()

        # Save decoded file with version prefix
        if self.decoder.headerChunk is not None and self.decoder.headerChunk.file_name is not None:
            try:
                file_name = self.decoder.headerChunk.file_name.decode("utf-8")
                Path(file_name).rename("v0_" + file_name)
            except FileNotFoundError:
                # File doesn't exist yet, which is fine
                pass

        # Save and rename the decoded file
        file_name = self.decoder.saveDecodedFile(
            last_chunk_len_format=self.decoder.config_map.get("last_chunk_len_str", "I"),
            return_file_name=True,
        )
        Path(file_name).rename("v0_" + file_name)

        logger.info(f"Base version saved as v0_{file_name}")
        return self.decoder

    def decode_to_version(
        self, base_dna_version_string: str, version: int
    ) -> typing.Dict[int, str]:
        """
        Decode the file up to the given version.

        Iteratively decodes each version from base (v0) up to the specified version,
        storing all intermediate versions on disk. Existing versions are loaded from
        disk and do not need to be decoded again.

        Args:
            base_dna_version_string: The base DNA version string marker.
            version: Target version number to decode to (must be >= 1).

        Returns:
            A dictionary mapping version numbers to file paths for all decoded versions.

        Raises:
            ValueError: If version is less than 1 or not available in the pool.
            RuntimeError: If decoding fails for any version.

        Note:
            Each version is based on the previous version, so all intermediate
            versions must be decoded in sequence.

        Example:
            >>> decoded = mv_decoder.decode_to_version("GAGCCAGTGAGTCGTA", version=3)
            >>> for ver, path in decoded.items():
            ...     print(f"Version {ver}: {path}")
        """
        if version < 1:
            raise ValueError("Version must be >= 1")

        # Check if requested version exists in pool
        max_version = self.get_versions_in_pool(base_dna_version_string)
        if version > max_version:
            raise ValueError(
                f"Version {version} is not available in the pool. "
                f"Maximum available version is {max_version}."
            )

        decoded_versions: typing.Dict[int, str] = {}
        id_len_format = self.decoder.config_map.get("id_len_format", "")

        # Process each version iteratively
        for i in range(1, version + 1):
            logger.info(f"Decoding version {i}/{version}...")

            version_seqs = self.get_sequences_for_version(base_dna_version_string, i)
            if not version_seqs:
                logger.warning(f"No sequences found for version {i}")
                continue

            solved_chunks: typing.Dict[int, typing.List[RU10Packet]] = {}
            bin_dna_version_str = tranlate_quat_to_byte(base_dna_version_string)

            for seq in version_seqs:
                # Revert seed spacing if it was used during encoding
                reverted_dna_str = self.decoder.revert_seed_spacing(seq, id_len_format)

                # Parse packet from DNA sequence
                packet = self.decoder.parse_raw_packet(
                    BytesIO(tranlate_quat_to_byte(reverted_dna_str)).read(),
                    crc_len_format=self.decoder.config_map.get("crc_len_format", ""),
                    number_of_chunks_len_format="",
                    packet_len_format=self.decoder.config_map.get("packet_len_format", ""),
                    id_len_format=id_len_format,
                )

                # Get used chunks list
                used_chunks_list = from_true_false_list(self.decoder.removeAndXorAuxPackets(packet))

                # Find version string position in packet
                find_result = packet.packed_used_packets.find(bin_dna_version_str)
                header_size = packet.get_packet_header_size()
                offset_pos = find_result - header_size + len(bin_dna_version_str)

                # Reduce packet to target chunk
                reduced = reduce_packet_to_chunk(packet.copy(), self, used_chunks_list[0])

                # Extract target chunk ID from reduced packet
                try:
                    (target_chunk,) = struct.unpack("<B", reduced.data[offset_pos : offset_pos + 1])
                except struct.error:
                    logger.warning(f"Failed to extract chunk ID from packet")
                    continue

                # Create mask to isolate version string region
                zeros_mask = np.zeros(len(packet.data), dtype=np.uint8)
                mask_start = offset_pos - len(bin_dna_version_str) - 1
                mask_end = offset_pos + 1
                if mask_start >= 0 and mask_end <= len(reduced.data):
                    zeros_mask[mask_start:mask_end] = np.frombuffer(
                        reduced.data[mask_start:mask_end], dtype=np.uint8
                    )

                # XOR to repair the packet data
                repaired_data = xor_numpy(packet.data, zeros_mask)

                # Create repaired packet
                res = RU10Packet(
                    repaired_data,
                    packet.used_packets,
                    packet.total_number_of_chunks,
                    packet.id,
                    read_only=True,
                    packet_len_format=packet.packet_len_format,
                    crc_len_format=packet.crc_len_format,
                    number_of_chunks_len_format=packet.number_of_chunks_len_format,
                    id_len_format=id_len_format,
                    save_number_of_chunks_in_packet=packet.total_number_of_chunks is None,
                )

                # Store packet for target chunk
                if used_chunks_list[target_chunk] not in solved_chunks:
                    solved_chunks[used_chunks_list[target_chunk]] = []

                # Reduce packet and store
                reduced_packet = reduce_packet_to_chunk(
                    res.copy(), self, used_chunks_list[target_chunk]
                )
                solved_chunks[used_chunks_list[target_chunk]].append(reduced_packet)

            # Combine data for chunks with multiple solutions
            for key, values in solved_chunks.items():
                if not values:
                    continue

                unique_data_parts = {bytes(v.data) for v in values}

                # XOR differences with original version
                tmp = np.zeros_like(self.decoder.GEPP.b[key], dtype=np.uint8)
                for part in unique_data_parts:
                    tmp = xor_numpy(tmp, xor_numpy(part, self.decoder.GEPP.b[key]))

                # Create insertion packet
                insertion_packet = values[0].copy()
                insertion_packet.data = xor_numpy(tmp, self.decoder.GEPP.b[key])
                self.decoder.packets.append(insertion_packet)
                self.decoder.GEPP.b[key] = insertion_packet.data

            # Save the new version
            try:
                file_name = self.decoder.saveDecodedFile(
                    last_chunk_len_format=self.decoder.config_map.get("last_chunk_len_str", "I"),
                    return_file_name=True,
                    ignore_crc=True,
                    print_to_output=False,
                )
                versioned_file = f"v{i}_{file_name}"
                Path(file_name).rename(versioned_file)
                decoded_versions[i] = versioned_file
                logger.info(f"Version {i} saved as {versioned_file}")
            except Exception as e:
                logger.error(f"Failed to save version {i}: {e}")
                raise RuntimeError(f"Failed to save version {i}: {e}")

        return decoded_versions

    @staticmethod
    def solve_lin_dep(
        a: typing.List[np.ndarray], b: np.ndarray
    ) -> typing.Optional[typing.List[np.ndarray]]:
        """
        Calculate which rows in vector `a` can be used to create the target `b`.

        This method tries combinations of up to 3 vectors from `a` to find which
        ones, when XORed together, produce the target vector `b`.

        Args:
            a: A list of numpy arrays where each array is either used to create b or not.
            b: The target numpy array to produce.

        Returns:
            A list of arrays from `a` that can be XORed to produce `b`,
            or None if no solution exists with up to 3 vectors.

        Note:
            This method checks combinations of 1, 2, and 3 vectors. For larger
            combinations, consider using a more efficient algorithm.

        Example:
            >>> vectors = [np.array([1, 0, 1, 0]), np.array([0, 1, 0, 1])]
            >>> target = np.array([1, 1, 1, 1])
            >>> solution = MultiVersionDecoder.solve_lin_dep(vectors, target)
        """
        # Try combinations of 1, 2, and 3 vectors
        max_combinations = min(3, len(a))
        for i in range(1, max_combinations + 1):
            for comb in combinations(a, i):
                if len(comb) > 1:
                    r = reduce(
                        lambda x, y: xor_numpy(x.astype("uint8"), y.astype("uint8")),
                        comb,
                    )
                else:
                    r = comb[0]

                if np.array_equal(r.astype("uint8"), b):
                    return [x.astype("uint8") for x in comb]

        return None

    def repair_and_store_by_packet(
        self,
        chunk_id: int,
        packet_id: int,
        hex_value: str,
        clear_working_dir: bool = False,
        correctness_function: typing.Optional[typing.Callable[[np.ndarray], bool]] = None,
    ) -> str:
        """
        Repair a chunk and store the result, trying different possible corrupt packets.

        This function is used when there are multiple invalid packets to save multiple
        versions, where each saved version uses a different possible packet to repair
        the chunk.

        Args:
            chunk_id: ID of the chunk to repair.
            packet_id: ID of the packet suspected to be corrupt.
            hex_value: Hexadecimal value to use for repair.
            clear_working_dir: If True, clear the working directory before saving.
            correctness_function: Optional function to verify repair correctness.
                Takes GEPP.b as input and returns True if correct.

        Returns:
            The name of the saved file.

        Raises:
            ValueError: If repair fails.

        Example:
            >>> filename = mv_decoder.repair_and_store_by_packet(
            ...     chunk_id=5,
            ...     packet_id=3,
            ...     hex_value="A1B2C3D4",
            ...     clear_working_dir=True
            ... )
        """
        # Backup current GEPP state
        bkp_A = self.decoder.GEPP.A.copy()
        bkp_b = self.decoder.GEPP.b.copy()

        # Perform manual repair
        self.manual_repair(chunk_id, packet_id, hex_value)

        # Setup working directory
        working_dir = "multi_file_repair"
        if clear_working_dir:
            if Path(working_dir).exists():
                shutil.rmtree(working_dir)
            Path(working_dir).mkdir(parents=True, exist_ok=True)

        # Parse header if using header chunk
        self.parse_header("I")
        is_correct = False

        if self.headerChunk is not None and self.headerChunk.checksum_len_format is not None:
            is_correct = self.is_checksum_correct()
        elif correctness_function is not None:
            is_correct = correctness_function(self.decoder.GEPP.b)

        # Save decoded file
        try:
            filename = self.decoder.saveDecodedFile(return_file_name=True, print_to_output=False)
        except ValueError as ve:
            filename = ve.args[1] if len(ve.args) > 1 else "unknown"

        # Rename with metadata
        _file = Path(filename)
        prefix = "CORRECT_" if is_correct else ""
        stem = f"{prefix}{_file.stem}_{chunk_id}_{packet_id}"
        _new_file = _file.rename(Path(working_dir) / f"{stem}{_file.suffix}")

        # Restore GEPP state
        self.decoder.GEPP.A = bkp_A
        self.decoder.GEPP.b = bkp_b

        return _new_file.name


def init_args() -> argparse.Namespace:
    """
    Parse command-line arguments for MultiVersionDecoder.

    Returns:
        Parsed arguments namespace.

    Example:
        >>> args = init_args()
        >>> print(f"Using config: {args.ini}")
    """
    parser = argparse.ArgumentParser(
        description="Multi-Version Decoder for NOREC4DNA encoded files"
    )
    parser.add_argument(
        "--ini",
        metavar="ini",
        type=str,
        help="Configuration file (INI format)",
        default="/home/michael/Code/DR4DNA/eval/sleeping_beauty_no_error.ini",
    )

    # Metadata argument group (mutually exclusive)
    metadata_arg_group = parser.add_mutually_exclusive_group(required=False)
    metadata_arg_group.add_argument(
        "--metadata_file",
        metavar="metafile",
        type=str,
        help="File containing metadata in FASTA format (comma-separated list)",
    )
    metadata_arg_group.add_argument(
        "--metadata",
        metavar="dmeta",
        type=str,
        help="Comma-separated list of metadata DNA sequences",
    )

    return parser.parse_args()


def main() -> None:
    """
    Main entry point for MultiVersionDecoder CLI.

    This function parses command-line arguments, initializes the decoder,
    and performs base version decoding followed by version 1 decoding.
    """
    parsed_args = init_args()

    # Load metadata
    metadata_list: typing.List[str] = []
    if parsed_args.metadata_file is not None:
        # Split and parse each metadata file
        for metadata_file in parsed_args.metadata_file.split(","):
            try:
                fasta_entries = load_fasta(metadata_file)
                metadata_list.extend(fasta_entries.values())
                logger.info(f"Loaded {len(fasta_entries)} metadata sequences from {metadata_file}")
            except Exception as e:
                logger.error(f"Failed to load metadata from {metadata_file}: {e}")
    elif parsed_args.metadata is not None:
        metadata_list = parsed_args.metadata.split(",")
        logger.info(f"Using {len(metadata_list)} inline metadata sequences")

    # Initialize decoder from config
    try:
        config = ConfigReadAndExecute(parsed_args.ini)
        decoder = config.execute(return_decoder=True)[0]
    except Exception as e:
        logger.error(f"Failed to initialize decoder: {e}")
        raise

    # Initialize toolkit and decoder
    semi_automatic_solver = SemiAutomaticReconstructionToolkit(decoder)

    # Display file with chunk borders
    print(semi_automatic_solver.view_file_with_chunkborders(False, False, "I"), flush=True)

    # Initialize multi-version decoder
    mv_decoder = MultiVersionDecoder(decoder, metadata_list)

    # Define base version string
    BASE_VERSION_STRING = "GAGCCAGTGAGTCGTA"

    # Decode base version
    logger.info("Decoding base version...")
    mv_decoder.decode_base_version(BASE_VERSION_STRING)

    # Decode to version 1
    logger.info("Decoding version 1...")
    mv_decoder.decode_to_version(BASE_VERSION_STRING, 1)

    logger.info("Decoding complete!")


if __name__ == "__main__":
    main()
