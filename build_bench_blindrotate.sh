#! /bin/bash

cmake -B CMAKE_BUILD -DINSTALL="OFF" -DTORUS_TYPE="36" -DPRINTER_ON="ON" -DENABLE_TIMER="ON"
cmake --build CMAKE_BUILD -j$(nproc)