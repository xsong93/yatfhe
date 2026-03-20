cmake_minimum_required(VERSION 3.10)
message(STATUS "Building base...")
project(base)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(yatfhe)
add_subdirectory(ya_benchmark)
