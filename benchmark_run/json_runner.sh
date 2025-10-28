#! /bin/bash
BENCHMARK_EXEC="bench_keygen"
OUTPUT_PREFIX="result_i5"

# Run the benchmark
"./$BENCHMARK_EXEC" \
    --benchmark_out="${OUTPUT_PREFIX}_${BENCHMARK_EXEC}.json" \
    --benchmark_out_format=json