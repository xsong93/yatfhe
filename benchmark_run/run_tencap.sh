#!/bin/bash
# Tenant ceiling at a common service rate.
#
set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
mkdir -p "$OUTDIR"
cd "$B/server" || exit 1
RATE="${RATE:-30}"
SLO="${SLO:-300}"
REQS="${REQS:-8000}"
CAPUNITS="${CAPUNITS:-100}" # memory budget, in GINX-key bytes
METHODS="${METHODS:-tfhe wwl24 ours}"
declare -A FNAME=( [tfhe]=tfhe [wwl24]=wwl+24 [ours]=ours )
# slots per budget unit, from the 170.11 / 85.30 / 42.81 MB key sizes (1 : 2 : 4)
declare -A SLOTS=( [tfhe]=1 [wwl24]=2 [ours]=4 )
# population grid
grid() {
    case "$2" in
        tfhe)  echo "50 100 112 118 125 150 200" ;;          # holds 112, fails 118
        wwl24) echo "50 100 200 300 350 362 368 375 400" ;;  # holds 362, fails 368
        ours)  echo "50 100 200 400 600 700 800" ;;         # holds the whole pool
    esac
}

tag_for() { echo "tencap${1}_u${2}"; } # tag carries the memory budget so cells of different budgets cannot collide
failures=0
echo "tenant-ceiling sweep start $(date -Is) rate=$RATE req/s budgets='$CAPUNITS' out=$OUTDIR"
for units in $CAPUNITS; do
for m in $METHODS; do
    for u in $(grid "$units" "$m"); do
        tag="$(tag_for "$units" "$u")"
        f="$OUTDIR/serving_${FNAME[$m]}_${tag}.json"
        timeout 7200 taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
            --workers=8 --cap=$(( units * ${SLOTS[$m]} )) --users="$u" --rate="$RATE" \
            --requests="$REQS" --warm=500 --batch=1 --arrival=poisson --sched=soft \
            --slo="$SLO" --seed=1 --devstat=nvme0n1 --tag="$tag" >/dev/null 2>&1
        rc=$?
        if [ "$rc" -ne 0 ] || [ ! -f "$f" ]; then
            echo "  FAILED $m users=$u (exit $rc, artifact $( [ -f "$f" ] && echo present || echo absent ))"
            failures=$((failures + 1))
            continue
        fi
        python3 - "$f" "$m" "$u" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); s = d['sojourn']; q = d['queue_delay']
v = d.get('slo_violation_rate_300')
ok = s['p99_9_ms'] <= d['config']['slo_ms'] and v is not None and v <= 0.001
print("  %-7s users=%-4s slots=%-3s p99.9 %8.1f ms  v300 %6.3f%%  queue p99.9 %9.2f  "
      "hit %.3f  thr %6.2f  %s" % (
          sys.argv[2], sys.argv[3], d['config']['cache_capacity_keys'], s['p99_9_ms'],
          100 * (v if v is not None else float('nan')), q['p99_9_ms'],
          d['hit_rate'], d['throughput_rps'], "inside" if ok else "OUTSIDE"))
PY
    done
done
done
if [ "$failures" -ne 0 ]; then
    echo "tenant-ceiling sweep INCOMPLETE: $failures cell(s) missing $(date -Is)" >&2
    exit 1
fi
echo "tenant-ceiling sweep done $(date -Is)"
