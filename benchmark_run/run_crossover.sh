#!/bin/bash
# Locate the CONCURRENT crossover in batch size B

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$B/server" || exit 1
declare -A CAP=( [tfhe]=10 [wwl24]=20 [ours]=40 )
echo "crossover sweep start $(date -Is)"
for b in 1 2 3 5 7; do
  rate=$(( 200 / b ))          # comfortably above capacity for every method
  for m in tfhe wwl24 ours; do
    taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m --workers=8 --cap=${CAP[$m]} \
        --users=50 --rate=$rate --requests=600 --warm=100 --batch=$b \
        --arrival=poisson --sched=soft --slo=300 --seed=1 \
        --tag="xo_b${b}" 2>&1 | grep -E "^  throughput" | sed "s/^/  [B=$b $m] /"
  done
done
echo "crossover sweep done $(date -Is)"
