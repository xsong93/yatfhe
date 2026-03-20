#!/bin/bash

# Configuration variables
BENCHMARK_EXEC="bench_blindrotate_all"
OUTPUT_PREFIX="result_i5_binary"
TEST_PATTERNS=(
    "BlindRotateBenchmark/JP22_MULTITHREAD/batchSize:16/tasksPerThread:9"
    "BlindRotateBenchmark/OURS_MULTITHREAD/batchSize:16/tasksPerThread:9"
    "BlindRotateBenchmark/MP21_SINGLETHREAD"
    "BlindRotateBenchmark/JP22_SINGLETHREAD"
    "BlindRotateBenchmark/OURS_SINGLETHREAD"
)

# Generate the filter regex by joining patterns with |
FILTER_REGEX=$(IFS="|"; echo "${TEST_PATTERNS[*]}")

# Run the benchmark
"./$BENCHMARK_EXEC" \
    --benchmark_out="${OUTPUT_PREFIX}_${BENCHMARK_EXEC}.json" \
    --benchmark_out_format=json \
    --benchmark_filter="$FILTER_REGEX"

# Verify output
if [ -f "${OUTPUT_PREFIX}_${BENCHMARK_EXEC}.json" ]; then
    echo "Benchmark completed successfully. Results saved to ${OUTPUT_PREFIX}_${BENCHMARK_EXEC}.json"
else
    echo "Error: Benchmark output file was not created!"
    exit 1
fi