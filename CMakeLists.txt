cmake_minimum_required(VERSION 3.22)
project(hls-yatfhe)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(yatfhe)
add_subdirectory(benchmark)
add_subdirectory(google_tests)


#add_subdirectory(ntl)

option(INSTALL_HEADERS "Install headers to system directories" OFF)