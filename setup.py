import sys

import numpy
from setuptools import Extension, setup

is_windows = sys.platform.startswith("win")

extra_compile_args = []
define_macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
extra_link_args = []

if is_windows:
    extra_compile_args = ["/O2", "/arch:AVX2", "/fp:fast"]
else:
    extra_compile_args = [
        "-D__BUILTIN_POPCOUNT",
        "-O3",  # Maximum optimization
        "-march=native",  # Use all CPU features
        "-mtune=native",  # Tune for current CPU
        "-ffast-math",  # Fast floating point
        "-funroll-loops",  # Unroll loops
        "-mavx2",  # Enable AVX2 if available
        "-msse4.2",  # Enable SSE4.2
        "-ftree-vectorize",  # Auto vectorization
    ]

ext_modules = [
    Extension(
        "norec4dna.cdnarules",
        sources=["csrc/cdnarules.c"],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        include_dirs=[numpy.get_include()],
        define_macros=define_macros,
    )
]
setup(
    package_dir={"": "src"},
    py_modules=["cdnarules"],
    ext_modules=ext_modules,
)
