cmake_minimum_required(VERSION 3.13.0)
cmake_policy(SET CMP0148 NEW)

project(JET_ENERGY_FUNCTIONS)

set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

set(CMAKE_PREFIX_PATH "C:/vcpkg/installed/x64-windows/share/pybind11/")
find_package(pybind11 REQUIRED)

pybind11_add_module(JET_ENERGY_FUNCTIONS example1.cpp)