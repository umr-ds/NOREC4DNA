import os
import sys
import platform
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

try:
    with open("README.md", encoding="utf-8") as f:
        long_description = f.read()
except IOError:
    long_description = ""

lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = os.path.join(lib_folder, "requirements.txt")
install_requires = []
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()


class BuildExtWithNumpy(build_ext):
    def finalize_options(self):
        super().finalize_options()
        # Delay numpy import until build time
        import numpy
        self.include_dirs.append(numpy.get_include())

    def build_extensions(self):
        # Get compiler type
        compiler_type = self.compiler.compiler_type

        # Platform-specific optimizations
        if compiler_type == 'msvc':
            # MSVC-specific flags
            extra_compile_args = [
                '/O2',  # Maximum optimization
                '/Ot',  # Favor speed over size
                '/GL',  # Whole program optimization
                '/arch:AVX2',  # Use AVX2 if available
            ]
            extra_link_args = ['/LTCG']  # Link-time code generation
        else:
            # GCC/Clang flags
            extra_compile_args = [
                '-O3',  # Aggressive optimization
                '-ffast-math',  # Fast floating point
                '-funroll-loops',  # Loop unrolling
                '-finline-functions',  # Aggressive inlining
                '-flto',  # Link-time optimization
                '-fomit-frame-pointer',  # Remove frame pointer
                '-D__BUILTIN_POPCOUNT',  # Enable popcount intrinsic
            ]

            # Architecture-specific optimizations
            machine = platform.machine().lower()
            if machine in ('x86_64', 'amd64'):
                extra_compile_args.extend([
                    '-march=native',  # CPU-specific optimizations
                    '-mtune=native',  # Tune for local CPU
                    '-msse4.2',  # SSE 4.2 instructions
                    '-mavx2',  # AVX2 vectorization
                ])
            elif machine.startswith('arm') or machine == 'aarch64':
                extra_compile_args.extend([
                    '-mcpu=native',  # ARM-specific optimizations
                    '-mfpu=neon',  # NEON SIMD for ARM
                ])

            extra_link_args = ['-flto', '-Wl,--strip-all']

        # Apply optimizations to all extensions
        for ext in self.extensions:
            if not hasattr(ext, 'extra_compile_args'):
                ext.extra_compile_args = []
            if not hasattr(ext, 'extra_link_args'):
                ext.extra_link_args = []

            ext.extra_compile_args.extend(extra_compile_args)
            ext.extra_link_args.extend(extra_link_args)

            # Add OpenMP support if available
            if compiler_type != 'msvc':
                ext.extra_compile_args.append('-fopenmp')
                ext.extra_link_args.append('-fopenmp')
            else:
                ext.extra_compile_args.append('/openmp')

        super().build_extensions()


# Define the extension
cdnarules_extension = Extension(
    "cdnarules",
    ["cdnarules.c"],  # Use your optimized C file
    language='c',
    # include_dirs will be set by BuildExtWithNumpy
)

setup(
    name="norec4dna",
    version="0.2.0",  # Bump version for optimized release
    description="NOREC4DNA - a Fountain Code based approach to DNA-Storage (Optimized)",
    author="Michael Schwarz",
    author_email="",
    url="",
    packages=find_packages(),
    python_requires=">=3.8",
    setup_requires=[
        "numpy>=1.26",
        "setuptools>=65.0",
        "wheel",
    ],
    install_requires=[
                         "numpy>=1.26",
                     ] + install_requires,
    zip_safe=False,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: C",
        "Operating System :: OS Independent",
    ],
    keywords="dna storage fountain codes bioinformatics",
    ext_modules=[cdnarules_extension],
    cmdclass={"build_ext": BuildExtWithNumpy},
    # Optional: Add entry points if your package provides command-line tools
    # entry_points={
    #     'console_scripts': [
    #         'norec4dna=norec4dna.cli:main',
    #     ],
    # },
)