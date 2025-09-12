cmake_minimum_required(VERSION 3.22)
message(STATUS "Building base...")
project(base)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(yatfhe)
add_subdirectory(benchmark)
add_subdirectory(google_tests)
