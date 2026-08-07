"""Compatibility shim for the shared metadata mapping module."""

from __future__ import annotations

from importlib import import_module

_metadata_mapping = import_module("norec4dna_multiversion.metadata_mapping")

main = _metadata_mapping.main


def __getattr__(name: str):
    return getattr(_metadata_mapping, name)


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(dir(_metadata_mapping)))


if __name__ == "__main__":
    import sys

    sys.exit(main())
