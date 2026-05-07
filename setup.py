"""Build configuration for the C++ extension only.

All package metadata (name, version, authors, classifiers, deps, …) lives
in pyproject.toml's [project] table per PEP 621. setuptools picks it up
automatically; this file just declares the C++ extension that PEP 621
can't express.
"""

from setuptools import Extension, setup

import pybind11

ext_modules = [
    Extension(
        "pycirclemedianfilter",
        sources=[
            "filters/CMF_python_bindings.cpp",
            "filters/CMF_library.cpp",
        ],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["-std=c++11"],
    ),
]

setup(ext_modules=ext_modules)
