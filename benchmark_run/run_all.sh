#!/bin/bash

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
say(){ echo "===== $1  $(date -Is) ====="; }

say "1/9 pressure sweep (5 seeds x 6 caps x 2000 req)"
REQS=2000 SEEDS="1 2 3 4 5" "$B/run_sweep.sh"

say "2/9 component (100 warm / 100 cold)"
( cd "$B/server" && taskset -c 0-7 ./blindrotate_comp --out="$OUTDIR" --warmreps=100 --coldreps=100 )

say "3/9 thread scaling + CCD placement"
"$B/run_threads.sh"

say "4/9 tenant sensitivity"
taskset -c 0-7 "$B/run_tenants.sh"

say "5/9 SLO-capacity B=1  (8000 req/point, bracketed grids)"
"$B/run_slo.sh"

say "6/9 decision-tree depth sweep (B = 1/2/4/8)"
"$B/run_trees.sh"

say "7/9 concurrent crossover"
"$B/run_crossover.sh"

say "8/9 storage tiers (SATA then N2)"
"$B/run_sata.sh"
"$B/run_n2.sh"

say "9/9 sharded two-pool placement (one serving pool per CCD)"
"$B/run_shard.sh"

echo "===== COMPLETE $(date -Is) ====="
