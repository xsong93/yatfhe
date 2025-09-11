cmake_minimum_required(VERSION 3.22)
project(hls-yatfhe)

set(CMAKE_CXX_STANDARD 17)

option(HEXL_REQUIRED "Set to true if you want to use HEXL library. Otherwise, set it to false." TRUE)

add_subdirectory(yatfhe)
add_subdirectory(benchmark)
add_subdirectory(google_tests)


#add_subdirectory(ntl)

option(INSTALL_HEADERS "Install headers to system directories" FALSE)