# yatfhe: Yet another TFHE library

## Repository Structure
```
yatfhe/
├── yatfhe/ # Core library source code (yatfhe_lib)
├── ya_benchmark/ # Benchmark programs, one executable per .cpp
├── google_tests/ # GoogleTest suite (google_tests_run)
├── CMakeLists.txt # Top-level CMake build configuration
├── build_bench_blindrotate.sh # Configure + build everything into CMAKE_BUILD/
├── cp_bench.sh # Copies the built benchmarks to benchmark_run/server/
├── gtest_run.sh # Runs the test binary, optional --gtest_filter argument
├── cmake_run.sh # Configure + `make install` the library system-wide
├── clean_install.sh # Removes CMAKE_BUILD/ and the installed headers/libs
├── .gitignore
├── LICENSE.txt # License file
└── README.md
```

## Dependencies

- **CMake** (>= 3.10)
- **C++17** compatible compiler
- **pkg-config** -- used to locate GMP
- **GMP** (GNU Multiple Precision Arithmetic Library)
  - Please refer to https://gmplib.org/ for the installation instructions on your system.
  - The build finds it through pkg-config, so `gmp.pc` has to be visible. Installing
    to a prefix outside the default search path means exporting it yourself.
- **Intel HEXL** (Homomorphic Encryption Acceleration Library)
  - You can install it manually from https://github.com/intel/hexl.
  - Point `HEXL_ROOT` at the install if it is not on the default CMake search path (see Build Instructions).
- **nlohmann/json** -- required by the `blindrotate_cache` benchmark, which writes its results as json.
  - Debian/Ubuntu: `sudo apt install nlohmann-json3-dev`. Otherwise, install the
    header-only library from https://github.com/nlohmann/json.

## Build Instructions

1. **Clone the repository**
2. **Configure and build**
    ```
    mkdir -p benchmark_run/server
    chmod +x build_bench_blindrotate.sh
    ./build_bench_blindrotate.sh
    bash cp_bench.sh
    ```

    The script configures and builds into `CMAKE_BUILD/`, and reads three
    environment variables:

    | variable | default | meaning |
    | --- | --- | --- |
    | `TORUS` | `32` | Torus type, `32` or `56` (`-DTORUS_TYPE`) |
    | `TARGET_ARCH` | `native` | ISA baseline: `native`, `avx512`, `avx2`, `generic` |
    | `HEXL_ROOT` | `/usr/local/lib` for `native`/`avx512`, else unset | HEXL install to put on `CMAKE_PREFIX_PATH` |

    `TARGET_ARCH` governs this project's own code only -- the HEXL you link has to
    match separately, which is what `HEXL_ROOT` selects. Use an AVX-512-tuned HEXL
    for `avx512`. A HEXL built without it costs roughly 4x on the NTT.
    Set `HEXL_ROOT=` (empty) to force the system install.

3. **Running Benchmarks**
    ```
    cd benchmark_run/server
    ./gen_benchkeys      # first time only, writes the BSK_*.bin key files
    ./blindrotate_cache  # writes one json per method per cache capacity
    ```

    On a multi-CCD CPU, pin the run to one CCD. The pipeline hands off buffers
    between threads once per LWE coefficient, so spreading those threads across core
    complexes costs ~45% and makes the latency bimodal. Find a CCD's CPUs from the
    CPUs sharing its L3, then pin:
    ```
    cat /sys/devices/system/cpu/cpu0/cache/index3/shared_cpu_list  # e.g. 0-7,16-23
    taskset -c 0-7 ./blindrotate_cache
    ```
    The thread pool sizes itself from the CPU affinity mask, so no other change is
    needed. Single-die CPUs need no pinning.

    On a NUMA machine (multi-socket, or one socket with Sub-NUMA Clustering enabled),
    bind CPU *and* memory to a single node instead. Crossing sockets costs more than
    crossing CCDs, since it adds remote memory to the remote cache traffic:
    ```
    numactl --hardware                                        # nodes and their CPUs
    numactl --cpunodebind=0 --membind=0 ./blindrotate_cache
    ```
    Binding costs nothing even on a large node: one blind rotation keeps at most
    `3 * lApprox` tasks in flight, and measured speedup saturates by 6-8 cores. Put
    the spare cores to work with concurrent bootstraps in separate processes rather
    than more threads per bootstrap -- the pool is a process-wide singleton behind a
    single queue, and oversized pools slow the per-coefficient barriers down.


## License
This project is licensed under the terms of the MIT License. See LICENSE.txt for details.