# -*- coding: utf-8 -*-
from setuptools import find_packages
from distutils.core import setup, Extension
import os

try:
    long_description = open("README.md").read()
except IOError:
    long_description = ""
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + '/requirements.txt'
install_requires = []  # Here we'll get: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()


class get_numpy_include(object):
    def __str__(self):
        import numpy
        return numpy.get_include()


setup(
    name="norec4dna",
    version="0.1.2",
    description="NOREC4DNA - a Fountain Code based approach to DNA-Storage",
    author="Michael Schwarz",
    packages=find_packages(),
    setup_requires=["numpy~=1.23.5"],
    install_requires=install_requires,
    zip_safe=False,
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    ext_modules=[Extension('cdnarules', ['cdnarules.c'], include_dirs=[get_numpy_include()]), ],
    include_dirs=[get_numpy_include()]
)
