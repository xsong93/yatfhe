# yatfhe: Yet another TFHE library

## Repository Structure
```
yatfhe/
├── yatfhe/ # Core library source code (yatfhe_lib)
├── ya_benchmark/ # Benchmark programs
├── google_tests/ # GoogleTest suite
├── build_bench_blindrotate.sh # Configure + build everything into CMAKE_BUILD/
├── cp_bench.sh # Copies the built benchmarks to benchmark_run/server/
├── gtest_run.sh # Runs the test binary, optional --gtest_filter argument
├── cmake_run.sh # Configure + `make install` the library system-wide
├── clean_install.sh # Removes CMAKE_BUILD/ and the installed headers/libs
├── LICENSE.txt
└── README.md
```

## Dependencies

- **CMake** (>= 3.10)
- **C++17** compatible compiler
- **GMP** (GNU Multiple Precision Arithmetic Library)
  - Please refer to https://gmplib.org/ for the installation instructions on your system.
  - The build finds it through pkg-config, so `gmp.pc` has to be visible on the default search path.
- **Intel HEXL** (Homomorphic Encryption Acceleration Library)
  - You can install it manually from https://github.com/intel/hexl.
  - Override by `HEXL_ROOT` if it is not on the default CMake search path (see Build Instructions).
- **nlohmann/json**
  - Debian/Ubuntu: `sudo apt install nlohmann-json3-dev`. Otherwise, install the header-only library from https://github.com/nlohmann/json.

## Build Instructions

1. **Clone the repository**
2. **Configure and build**
    ```
    mkdir -p benchmark_run/server
    build_bench_blindrotate.sh
    ./build_bench_blindrotate.sh
    bash cp_bench.sh
    ```

    The script configures and builds into `CMAKE_BUILD/`, and reads three environment variables:

    | variable | default | meaning                                             |
    | --- | --- |-----------------------------------------------------|
    | `TORUS` | `32` | Torus type, `32` or `56` (`-DTORUS_TYPE`)           |
    | `TARGET_ARCH` | `native` | ISA baseline: `native`, `avx512`, `avx2`, `generic` |
    | `HEXL_ROOT` | `/usr/local/lib` | HEXL install path                                   |

    Use an AVX-512 HEXL for `TARGET_ARCH=avx512`. A HEXL built without it costs roughly 4x on the NTT.

3. **Running Benchmarks**
    ```
    cd benchmark_run/server
    ./gen_benchkeys      # first time only, writes the BSK_*.bin key files
    ./blindrotate_cache  # writes one json per method per cache capacity
    ```

    On a multi-CCD CPU, pin the run to one CCD. The pipeline hands off buffers between threads, so spreading those 
    threads across core complexes costs, making the latency bimodal. Find a CCD's CPUs from the CPUs sharing its L3, then pin:
    ```
    cat /sys/devices/system/cpu/cpu0/cache/index3/shared_cpu_list  # e.g. 0-7,16-23
    taskset -c 0-7 ./blindrotate_cache
    ```

## License
This project is licensed under the terms of the MIT License. See LICENSE.txt for details.