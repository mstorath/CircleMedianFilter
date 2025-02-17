from setuptools import setup, find_packages, Extension
import pybind11

ext_modules = [
    Extension(
        "pycmf.pycmf", 
        ["filters/CMF_python_bindings.cpp", "filters/CMF_library.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["-std=c++11"],  # specify C++11 standard (or higher)
    ),
]

setup(
    name="pycmf",
    version="0.1.0",
    author="Martin Storath",
    author_email="martin.storath@thws.de",
    description="A Circle Median Filter package for Python",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mstorath/CircleMedianFilter", 
    packages=find_packages(),
    ext_modules=ext_modules,
    install_requires=["numpy", "pybind11"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)