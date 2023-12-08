cmake_minimum_required(VERSION 3.25)
project(hls-yatfhe)

file(GLOB SOURCES
        *.cpp
        src/*.cpp
        )

set(CMAKE_CXX_STANDARD 17)

add_executable(hls-yatfhe ${SOURCES})

include_directories(
        .
        ./util
        ./include
        )
