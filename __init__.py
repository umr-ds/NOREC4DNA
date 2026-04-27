"""Compatibility shim for legacy ``NOREC4DNA.*`` imports."""

from pathlib import Path

_ROOT = Path(__file__).resolve().parent
_SRC_PACKAGE = _ROOT / "src" / "norec4dna"
__path__ = [str(_ROOT), str(_SRC_PACKAGE)]

if __spec__ is not None:
    __spec__.submodule_search_locations[:] = __path__

del Path
del _ROOT
del _SRC_PACKAGE
