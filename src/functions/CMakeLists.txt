cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0148 NEW)
project(computationalfunctions)

set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

set(CMAKE_PREFIX_PATH "C:/vcpkg/installed/x64-windows/share/pybind11/")
find_package(pybind11 REQUIRED)

set(CMAKE_PREFIX_PATH "C:/vcpkg/installed/x64-windows/share/yaml-cpp/")
find_package(yaml-cpp REQUIRED)

pybind11_add_module(computationalfunctions ComputationalFunctions.cpp)
target_link_libraries(computationalfunctions PRIVATE yaml-cpp::yaml-cpp)