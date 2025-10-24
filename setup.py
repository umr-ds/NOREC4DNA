from setuptools import setup, Extension
import numpy, sys

is_windows = sys.platform.startswith("win")

extra_compile_args = []
define_macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]

if is_windows:
    extra_compile_args = ["/O2"]
else:
    extra_compile_args = ["-O3", "-D__BUILTIN_POPCOUNT"]

ext_modules = [
    Extension(
        "cdnarules",
        sources=["cdnarules.c"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=extra_compile_args,
        define_macros=define_macros,
    )
]
setup(ext_modules=ext_modules)
