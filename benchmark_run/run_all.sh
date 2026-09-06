#!/bin/bash

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
say(){ echo "===== $1  $(date -Is) ====="; }

# Deploy fresh binaries (build with TORUS=32 ./build_bench_blindrotate.sh first)
for f in blindrotate_cache blindrotate_comp blindrotate_serving gen_benchkeys \
         firststep_ablate stagea_ablate stagea_width_ablate stage_ratio_final; do
    cp "$B/../CMAKE_BUILD/ya_benchmark/$f" "$B/server/" 2>/dev/null || true
done

say "1/12 pressure sweep (5 seeds x 6 caps x 2000 req, pre-touch + settle warmup)"
REQS=2000 SEEDS="1 2 3 4 5" "$B/run_sweep.sh"

say "2/12 component ablation (100 warm / 100 cold, rows A C D E G)"
( cd "$B/server" && taskset -c 0-7 ./blindrotate_comp --out="$OUTDIR" --warmreps=100 --coldreps=100 )

say "3/12 thread scaling + CCD placement"
"$B/run_threads.sh"

say "4/12 tenant sensitivity"
taskset -c 0-7 "$B/run_tenants.sh"

say "5/12 SLO-capacity B=1  (8000 req/point, bracketed Poisson grids)"
"$B/run_slo.sh"

say "6/12 bursty-arrival SLO-capacity (4x/200ms-on, 0.25x/800ms-off)"
"$B/run_slo_burst.sh"

say "7/12 decision-tree depth sweep (B = 1/2/4/8)"
"$B/run_trees.sh"

say "8/12 concurrent crossover"
"$B/run_crossover.sh"

say "9/12 storage tiers (SATA then N2)"
"$B/run_sata.sh"
"$B/run_n2.sh"

say "10/12 sharded two-pool placement (one serving pool per CCD)"
"$B/run_shard.sh"

say "11/12 scheduler comparison (soft vs shared vs affinity)"
"$B/run_sched.sh"

say "12/12 paired micro-benchmarks (firststep, NTT-split, dispatch width, stage ratio)"
mkdir -p "$OUTDIR/paired"
for i in 1 2 3 4 5 6; do
    taskset -c 0-7 "$B/server/firststep_ablate" --reps=2000 --out="$OUTDIR/paired/fs$i" >/dev/null 2>&1
done
for i in 1 2 3 4; do
    taskset -c 0-7 "$B/server/stagea_ablate" --reps=2000 --out="$OUTDIR/paired/sa$i" >/dev/null 2>&1
done
for i in 1 2 3 4; do
    taskset -c 0-7 "$B/server/stagea_width_ablate" --reps=2000 --out="$OUTDIR/paired/sw$i" >/dev/null 2>&1
done
taskset -c 0-7 "$B/server/stage_ratio_final" --runs=10

echo "===== COMPLETE $(date -Is) ====="
