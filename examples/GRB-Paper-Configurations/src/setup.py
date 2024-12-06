from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "JET_ENERGY_FUNCTIONS",  # Module name
        ["example1.cpp"],       # C++ source file
        include_dirs=[pybind11.get_include()],
        language="c++"
    ),
]

setup(
    name="JET_ENERGY_FUNCTIONS",
    description="A C++ extension that evaluates the jet energy of a GRB using a GRB jet configuration structure.",
    ext_modules=ext_modules,
    install_requires=["pybind11"],
    zip_safe=False,
)