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
from importlib import import_module, util
from typing import Any

_coder = (
    import_module("norec4dna_multiversion.coder")
    if util.find_spec("norec4dna_multiversion.coder") is not None
    else None
)

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

MultiVersionCoder: Any
find_affected_chunks: Any
generate_dna_version_string: Any
get_current_file_version: Any
insert_dna_version_string: Any
insert_id_string: Any
reduce_packet_to_chunk: Any
find_insertion_position: Any
find_insertion_position_with_seed: Any
generate_new_packets: Any
create_perceptual_hash: Any
encoder_from_decoder: Any
add_packets: Any
decode_versions: Any
init_args: Any
main: Any

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

if _coder is not None:
    for exported_name in __all__:
        globals()[exported_name] = getattr(_coder, exported_name)
else:
    for exported_name in __all__:
        globals()[exported_name] = None


def __getattr__(name: str) -> Any:
    if _coder is not None and name in __all__:
        return getattr(_coder, name)
    if name in __all__:
        raise ImportError(
            "norec4dna_multiversion.coder is required for deprecated file_update_coding helpers"
        )
    raise AttributeError(name)


# Keep CLI functionality for backward compatibility
if __name__ == "__main__":
    warnings.warn(
        "Running file_update_coding.py directly is deprecated. "
        "Please use: python -m norec4dna_multiversion.coder",
        DeprecationWarning,
        stacklevel=2,
    )
    if _coder is None:
        raise ImportError("norec4dna_multiversion.coder is required to run this module")
    _coder.main()
