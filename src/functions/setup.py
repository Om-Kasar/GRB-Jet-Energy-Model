from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        'computationalfunctions',
        ['ComputationalFunctions.cpp'],
        include_dirs=[
            pybind11.get_include(),
            'C:/vcpkg/installed/x64-windows/include/'  # Path to yaml-cpp
        ],
        library_dirs=['C:/vcpkg/installed/x64-windows/lib/'],
        libraries=['yaml-cpp'],
        language='c++'
    )
]

setup(
    name='computationalfunctions',
    ext_modules=ext_modules,
    description = 'A Python wrapper for ComputationalFunctions.cpp.',
    zip_safe = False
)