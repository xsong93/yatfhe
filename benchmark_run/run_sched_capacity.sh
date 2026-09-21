#!/bin/bash
# Capacity under a SINGLE SHARED READY QUEUE.
# Separates what the device costs from what the scheduler costs
#
set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
mkdir -p "$OUTDIR"
cd "$B/server" || exit 1
declare -A CAP=( [tfhe]=10 [wwl24]=20 [ours]=40 )
declare -A FNAME=( [tfhe]=tfhe [wwl24]=wwl+24 [ours]=ours )
GRID="${GRID:-tfhe:1 tfhe:2 tfhe:3 tfhe:4 wwl24:16 wwl24:17 wwl24:18 ours:60 ours:64}"
SLO="${SLO:-300}"
REQS="${REQS:-8000}"

failures=0
echo "shared-queue capacity start $(date -Is) out=$OUTDIR"
for cell in $GRID; do
    m="${cell%%:*}"; r="${cell##*:}"
    timeout 20000 taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
        --workers=8 --cap=${CAP[$m]} --users=50 --rate=$r --requests=$REQS \
        --warm=500 --batch=1 --arrival=poisson --sched=shared \
        --slo=$SLO --seed=1 --devstat=nvme0n1 --tag="shared_r${r}" >/dev/null 2>&1
    rc=$?
    f="$OUTDIR/serving_${FNAME[$m]}_shared_r${r}.json"
    if [ "$rc" -ne 0 ] || [ ! -f "$f" ]; then
        echo "  FAILED $m rate=$r (exit $rc, artifact $( [ -f "$f" ] && echo present || echo absent ))"
        failures=$((failures + 1))
        continue
    fi
    python3 - "$f" "$m" "$r" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); s = d['sojourn']; q = d['queue_delay']
ok = s['p99_9_ms'] <= d['config']['slo_ms'] and d.get('slo_violation_rate_300', 1) <= 0.001
print("  %-7s rate=%-4s p99.9 %7.1f ms  v300 %6.3f%%  q p99.9 %6.2f  %s" % (
    sys.argv[2], sys.argv[3], s['p99_9_ms'], 100 * d.get('slo_violation_rate_300', -1),
    q['p99_9_ms'], "inside" if ok else "OUTSIDE"))
PY
done
if [ "$failures" -ne 0 ]; then
    echo "shared-queue capacity INCOMPLETE: $failures point(s) missing $(date -Is)" >&2
    exit 1
fi
echo "shared-queue capacity done $(date -Is)"
