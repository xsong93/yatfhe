cmake_minimum_required(VERSION 3.5...3.12)
message(STATUS "Building base...")
project(base)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(yatfhe)
#add_subdirectory(benchmark)
add_subdirectory(google_tests)
add_subdirectory(google_benchmark)
