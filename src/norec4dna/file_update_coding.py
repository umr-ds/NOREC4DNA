r"""
DEPRECATED: File Update Coding for NOREC4DNA Multi-Version Support.

This module is DEPRECATED. Please use MultiVersionCoder instead.

This module provides backward compatibility by re-exporting all functions
from MultiVersionCoder.py. All new code should import from MultiVersionCoder.

Migration Guide:
    OLD: from norec4dna.file_update_coding import find_affected_chunks
    NEW: from norec4dna_multiversion.coder import find_affected_chunks

    OLD: from norec4dna.file_update_coding import get_current_file_version
    NEW: from norec4dna_multiversion.coder import get_current_file_version

    OLD: from norec4dna.file_update_coding import generate_new_packets
    NEW: from norec4dna_multiversion.coder import generate_new_packets

    OLD: from norec4dna.file_update_coding import MultiVersionCoder
    NEW: from norec4dna_multiversion.coder import MultiVersionCoder

Example Usage (NEW):
    >>> from norec4dna_multiversion.coder import MultiVersionCoder
    >>> coder = MultiVersionCoder("existing_pool.ini")
    >>> max_version = coder.get_max_version_in_pool()
    >>> version, packets = coder.encode_new_version(new_file_data)
    >>> coder.save_updated_pool("updated_pool.fasta")

For more information, see MULTIVERSION_INPUT_FIX.md
"""

import logging
import warnings

# Issue deprecation warning on import
warnings.warn(
    "norec4dna.file_update_coding is deprecated. "
    "Please use norec4dna_multiversion.coder instead. "
    "See MULTIVERSION_INPUT_FIX.md for migration guide.",
    DeprecationWarning,
    stacklevel=2,
)

logger = logging.getLogger(__name__)
logger.warning(
    "norec4dna.file_update_coding is deprecated. Please use norec4dna_multiversion.coder instead."
)

# Also export types for type checking

# Re-export all functions and classes from MultiVersionCoder
from norec4dna_multiversion.coder import (
    MultiVersionCoder,
    add_packets,
    create_perceptual_hash,
    decode_versions,
    encoder_from_decoder,
    find_affected_chunks,
    find_insertion_position,
    find_insertion_position_with_seed,
    generate_dna_version_string,
    generate_new_packets,
    get_current_file_version,
    init_args,
    insert_dna_version_string,
    insert_id_string,
    main,
    reduce_packet_to_chunk,
)

__all__ = [
    "MultiVersionCoder",
    "find_affected_chunks",
    "generate_dna_version_string",
    "get_current_file_version",
    "insert_dna_version_string",
    "insert_id_string",
    "reduce_packet_to_chunk",
    "find_insertion_position",
    "find_insertion_position_with_seed",
    "generate_new_packets",
    "create_perceptual_hash",
    "encoder_from_decoder",
    "add_packets",
    "decode_versions",
    "init_args",
    "main",
]

# Keep CLI functionality for backward compatibility
if __name__ == "__main__":
    warnings.warn(
        "Running file_update_coding.py directly is deprecated. "
        "Please use: python -m norec4dna_multiversion.coder",
        DeprecationWarning,
        stacklevel=2,
    )
    main()
