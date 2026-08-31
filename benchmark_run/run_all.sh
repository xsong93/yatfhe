#!/bin/bash

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
say(){ echo "===== $1  $(date -Is) ====="; }

say "1/7 pressure sweep (5 seeds x 6 caps x 1000 req)"
REQS=1000 SEEDS="1 2 3 4 5" "$B/run_sweep.sh"

say "2/7 component (40 warm / 50 cold)"
( cd "$B/server" && taskset -c 0-7 ./blindrotate_comp --out="$OUTDIR" --warmreps=40 --coldreps=50 )

say "3/7 thread scaling + CCD placement"
"$B/run_threads.sh"

say "4/7 tenant sensitivity"
taskset -c 0-7 "$B/run_tenants.sh"

say "5/7 SLO-capacity B=1  (8000 req/point, bracketed grids)"
"$B/run_slo.sh"

say "6/7 concurrent crossover"
"$B/run_crossover.sh"

say "7/7 storage tiers (SATA then N2)"
"$B/run_sata.sh"
"$B/run_n2.sh"

echo "===== COMPLETE $(date -Is) ====="
