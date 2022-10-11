cmake_minimum_required(VERSION 3.10)

project(cpp_arithmetic_encoding)
add_subdirectory(extern/pybind11)
add_subdirectory(src)

link_libraries(ArCompress)
pybind11_add_module(cpp_string_lookup src/main.cpp)


SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")


find_package(pybind11 REQUIRED)  # or `add_subdirectory(pybind11)`

add_executable(executable src/main.cpp)
target_link_libraries(executable PRIVATE pybind11::embed)
