#! /bin/bash

cmake -B CMAKE_BUILD -DINSTALL="OFF" -DTORUS_TYPE="42" -DPRINTER_ON="ON" -DENABLE_TIMER="ON"
cmake --build CMAKE_BUILD -j$(nproc)

sudo cp CMAKE_BUILD/ya_benchmark/blindrotate_cache benchmark_run/server/.
sudo cp CMAKE_BUILD/ya_benchmark/gen_benchkeys benchmark_run/server/.