#!/bin/bash
# Scheduler comparison: soft vs shared vs pure affinity, 30 and 120 rps.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
mkdir -p "$OUTDIR"
cd "$B/server" || exit 1
echo "sched sweep start $(date -Is)"
for r in 30 120; do
  for s in soft shared affinity; do
    taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=ours --workers=8 --cap=40 \
        --users=50 --rate=$r --requests=800 --warm=300 --batch=1 \
        --arrival=poisson --sched=$s --slo=300 --seed=1 --devstat=nvme0n1 --tag="sched_${s}_r${r}" \
        2>&1 | grep -E "^  (mean|p99|fairness|locality)" | sed "s/^/  [r=$r $s] /"
  done
done
echo "sched sweep done $(date -Is)"
