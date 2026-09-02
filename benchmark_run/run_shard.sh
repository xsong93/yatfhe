#!/bin/bash
# Config D: two independent serving pools, one per CCD (sharded by tenant).
# Pool 0 pinned to cores 0-7, pool 1 to cores 8-15.
# Each process has its own worker set, soft-affinity scheduler and in-memory key cache.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
HARNESS="$B/server"
OUTDIR="${OUTDIR:-$B/out}"
mkdir -p "$OUTDIR"
cd "$HARNESS" || exit 1

# wait for the tree-depth sweep to release cores 0-7
while pgrep -f "blindrotate_serving.*tree_b" > /dev/null; do sleep 60; done

run_cell() {  # cpus method rate requests tag cap
  taskset -c "$1" "$HARNESS/blindrotate_serving" --out="$OUTDIR" --method="$2"       --workers=8 --cap="$6" --users=50 --rate="$3" --requests="$4" --warm=300       --batch=1 --arrival=poisson --sched=soft --slo=300 --seed=1 --tag="$5"       > "$OUTDIR/$5.log" 2>&1
}

echo "shard sweep start $(date -Is) out=$OUTDIR"
# saturated-throughput point: both pools driven above capacity
for m in tfhe wwl24 ours; do
  case $m in tfhe) cap=10;; wwl24) cap=20;; ours) cap=40;; esac
  run_cell 0-7  $m 200 800 "shard0_${m}" $cap &
  p0=$!
  run_cell 8-15 $m 200 800 "shard1_${m}" $cap &
  p1=$!
  wait "$p0" "$p1"
done
# pressure point (PR 5 anchored): each pool at 30 rps
for m in tfhe wwl24 ours; do
  case $m in tfhe) cap=10;; wwl24) cap=20;; ours) cap=40;; esac
  run_cell 0-7  $m 30 8000 "shardp0_${m}" $cap &
  p0=$!
  run_cell 8-15 $m 30 8000 "shardp1_${m}" $cap &
  p1=$!
  wait "$p0" "$p1"
done
echo "shard sweep done $(date -Is)"
