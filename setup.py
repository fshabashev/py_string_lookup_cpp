import sys

from pybind11 import get_cmake_dir
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
from glob import glob


__version__ = "0.0.5"

ext_modules = [
    Pybind11Extension("cpp_string_lookup",
        ["src/main.cpp"] + sorted(glob("src/arithmetic_coding/*.cpp")),
        # Example: passing in the version to the compiled code
        define_macros = [('VERSION_INFO', __version__)],
        ),
]

setup(
    name="cpp_string_low_mem_lookup",
    version=__version__,
    author="Fedor Shabashev",
    author_email="fedor.shabashev@gmail.com",
    url="https://github.com/fshabashev/py_string_lookup_cpp",
    description="Low memory string lookup",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6",
)
