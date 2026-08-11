# yatfhe: Yet another TFHE library

## Repository Structure
```
yatfhe/
├── benchmark_run/ # Benchmarking scripts and configuration
├── ya_benchmark/ # Additional benchmarking helpers
├── yatfhe/ # Core library source code
├── .gitignore
├── CMakeLists.txt # Main CMake build configuration
├── LICENSE.txt # License file
├── build_bench_blindrotate.sh # Script to build benchmarks
└── cp_bench.sh # Utility to copy benchmarks
```

## Dependencies

- **CMake** (>= 3.10)
- **C++17** compatible compiler
- **GMP** (GNU Multiple Precision Arithmetic Library)
  - Please refer to https://gmplib.org/ for the installation instructions on your system.
- **Intel HEXL** (Homomorphic Encryption Acceleration Library)
  - You can install it manually from https://github.com/intel/hexl.

## Build Instructions

1. **Clone the repository**
2. **Configure and build**
    ```
    mkdir benchmark_run/server
    chmod +x build_bench_blindrotate.sh
    ./build_bench_blindrotate.sh
    bash cp_bench.sh
    ```

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