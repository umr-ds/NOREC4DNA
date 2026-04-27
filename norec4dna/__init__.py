"""Compatibility shim that forwards in-repo imports to ``src/norec4dna``."""

from pathlib import Path

_SRC_PACKAGE = Path(__file__).resolve().parent.parent / "src" / "norec4dna"
__file__ = str(_SRC_PACKAGE / "__init__.py")
__path__ = [str(_SRC_PACKAGE)]

if __spec__ is not None:
    __spec__.submodule_search_locations[:] = __path__

with open(__file__, "rb") as _src_init:
    exec(compile(_src_init.read(), __file__, "exec"))

del Path
del _SRC_PACKAGE
del _src_init
