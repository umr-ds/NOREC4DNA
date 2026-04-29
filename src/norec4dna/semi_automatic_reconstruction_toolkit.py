"""Deprecated compatibility wrapper for the shared reconstruction toolkit."""

from __future__ import annotations

import warnings
from importlib import import_module
from typing import Any

from python_utils.types import deprecated

_reconstruction = import_module("norec4dna_multiversion.reconstruction")
_SharedReconstructionToolkit = _reconstruction.SemiAutomaticReconstructionToolkit

__all__ = ["SemiAutomaticReconstructionToolkit"]


@deprecated("Use norec4dna_multiversion.reconstruction.SemiAutomaticReconstructionToolkit instead.")
class SemiAutomaticReconstructionToolkit(_SharedReconstructionToolkit):
    """Backward-compatible alias for the shared reconstruction toolkit."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn(
            "norec4dna.semi_automatic_reconstruction_toolkit.SemiAutomaticReconstructionToolkit "
            "is deprecated; use norec4dna_multiversion.reconstruction.SemiAutomaticReconstructionToolkit instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)
