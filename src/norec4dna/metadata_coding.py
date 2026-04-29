"""Compatibility shim for the shared metadata coding module."""

from __future__ import annotations

from importlib import import_module

_metadata_coding = import_module("norec4dna_multiversion.metadata_coding")

main = _metadata_coding.main


def __getattr__(name: str):
    return getattr(_metadata_coding, name)


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(dir(_metadata_coding)))


if __name__ == "__main__":
    import sys

    sys.exit(main())
