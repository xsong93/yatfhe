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
    sudo bash cp_bench.sh
    ```

3. **Running Benchmarks**
```
   cd benchmark_run/server
   sudo ./gen_benkeys (only run this command for the first time)
   sudo ./blindrotate_cache (will create multiple json data files upon finish)
```

## License
This project is licensed under the terms of the MIT License. See LICENSE.txt for details.