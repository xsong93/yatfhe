#!/bin/bash
# Queue-depth instrumentation

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
mkdir -p "$OUTDIR"
cd "$B/server" || exit 1
declare -A CAP=( [tfhe]=10 [wwl24]=20 [ours]=40 )
declare -A FNAME=( [tfhe]=tfhe [wwl24]=wwl+24 [ours]=ours )

GRID="${GRID:-tfhe:8:8000 ours:8:8000 ours:62:8000 wwl24:16:8000 tfhe:200:600 ours:200:600 wwl24:200:600}"

failures=0
echo "queue-depth points start $(date -Is) out=$OUTDIR"
for cell in $GRID; do
    m="${cell%%:*}"; rest="${cell#*:}"; r="${rest%%:*}"; n="${rest##*:}"
    if [ "$r" = "200" ]; then tag="qd_${m}_sat"; else tag="qd_${m}_r${r}"; fi
    taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
        --workers=8 --cap=${CAP[$m]} --users=50 --rate=$r --requests=$n \
        --warm=500 --batch=1 --arrival=poisson --sched=soft \
        --slo=300 --seed=1 --devstat=nvme0n1 --tag="$tag" >/dev/null 2>&1
    rc=$?
    f="$OUTDIR/serving_${FNAME[$m]}_${tag}.json"
    if [ "$rc" -ne 0 ] || [ ! -f "$f" ]; then
        echo "  FAILED $m rate=$r requests=$n (exit $rc)"; failures=$((failures + 1)); continue
    fi
    python3 - "$f" "$m" "$r" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print("  %-7s rate=%-4s n=%-5s loads mean %5.2f peak %-3s  device reads mean %5.2f peak %-3s" % (
    sys.argv[2], sys.argv[3], d['config']['requests'], d['load_concurrency_mean'],
    d['load_concurrency_peak'], d['device_inflight_reads_mean'], d['device_inflight_reads_peak']))
PY
done
if [ "$failures" -ne 0 ]; then
    echo "queue-depth points INCOMPLETE: $failures point(s) missing $(date -Is)" >&2
    exit 1
fi
echo "queue-depth points done $(date -Is)"
