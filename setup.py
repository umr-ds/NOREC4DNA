from setuptools import setup, Extension
import numpy

ext_modules = [
    Extension(
        "cdnarules",
        sources=["cdnarules.c"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=["-O3", "-D__BUILTIN_POPCOUNT"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
    )
]

setup(
    ext_modules=ext_modules
)