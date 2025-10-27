#! /bin/bash
BENCHMARK_EXEC="bench_blindrotate_new"
OUTPUT_PREFIX="result_9950x3d"

# Run the benchmark
"./$BENCHMARK_EXEC" \
    --benchmark_out="${OUTPUT_PREFIX}_${BENCHMARK_EXEC}.json" \
    --benchmark_out_format=json