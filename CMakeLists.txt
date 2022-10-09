cmake_minimum_required(VERSION 3.10)

project(cpp_arithmetic_encoding)
add_subdirectory(extern/pybind11)
add_subdirectory(src)

link_libraries(ArCompress)
pybind11_add_module(cpp_string_lookup src/main.cpp)
