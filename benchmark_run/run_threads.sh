#!/bin/bash
# 5a: thread scaling and CCD placement
#
# The pool sizes from sched_getaffinity, so taskset sets both the thread count and core complex.
# Cap 100 with 50 tenants is pure compute and the numbers are not contaminated by I/O.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$(dirname "$0")/server" || exit 1
REQS="${REQS:-200}"
run() {  # run <label> <cpulist>
    local label="$1" cpus="$2"
    local n; n=$(taskset -c "$cpus" nproc)
    echo "=== $label  cpus=$cpus  threads=$n ==="
    taskset -c "$cpus" ./blindrotate_cache --out="$OUTDIR" --caps=100 --users=50 --reqs="$REQS" \
        --warm=1000 --seed=1 --tag="thr_${label}" 2>&1 | grep -E "^Saved"
}
# scaling within the V-Cache CCD
run t1  0
run t2  0-1
run t4  0-3
run t8  0-7
# placement at 8 threads
run t8b 8-15
run t8x 0-3,8-11
# wider
run t16smt 0-7,16-23
run t16phys 0-15
run t32 0-31
