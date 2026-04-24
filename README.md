# yatfhe: A Pipelined TFHE Framework for Multi-Tenant FHE

This repository contains the source code for the paper:

> **Mitigating I/O Jitter in Multi-Tenant FHE with a Co-Designed Pipelined Architecture**  
> *Submission to IEEE Transactions on Computers (under review).*

The code is provided for **review purposes only**. It will be made public with the final paper upon acceptance.

## Overview

This project implements a **pipelined TFHE bootstrapping framework** that hides I/O latency by overlapping key loading, scheme switching, and external product operations. It includes:

- A restructured bootstrapping key format (replacing the first RGSW with an RLWE).
- A three-stage pipeline (NS', SS', EP) with ping-pong buffers to resolve hazards.
- A multi‑tenant key‑cache simulator to evaluate performance under pressure.
- A fully functional TFHE framework with multiple optimizations (key grouping, mcrt, etc.).
- Google Test suites for performance and correctness.

## Repository Structure
```
yatfhe/
├── ya_benchmark/ # Benchmarking helpers
├── yatfhe/ # Core library source code
├── .gitignore
├── CMakeLists.txt # Main CMake build configuration
├── LICENSE.txt # License file
├── build_bench_blindrotate.sh # Script to build benchmarks
├── clean_install.sh # Clean install script
├── cmake_run.sh # CMake configuration wrapper
└── cp_bench.sh # Utility to copy benchmarks
```

## Dependencies

- **CMake** (>= 3.10)
- **C++17** compatible compiler (GCC 7+, Clang 6+, or MSVC 2019+)
- **GMP** (GNU Multiple Precision Arithmetic Library)
  - Please refer to https://gmplib.org/ for the installation instructions on your system.
- **Intel HEXL** (Homomorphic Encryption Acceleration Library)
   - You can install it manually from https://github.com/intel/hexl.

## Build Instructions

1. **Clone the repository**
   ```bash
   git clone https://anonymous.4open.science/r/yatfhe-FEDF/
   cd yatfhe-FEDF
   ```
2. **Configure and build**
    ```
    mkdir -p benchmark_run/server
    chmod +x build_bench_blindrotate.sh
    ./build_bench_blindrotate.sh
    ```

3. **Running Benchmarks**
```
   cd benchmark_run/server
   sudo ./gen_benkeys (this binary will generate required key files; only run for the first time)
   sudo ./blindrotate_cache (multiple json data files will be generated upon finish)
```

## License
This project is licensed under the terms of the MIT License. See LICENSE.txt for details.