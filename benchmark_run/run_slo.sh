#!/bin/bash
# SLO-capacity sweep: the highest offered rate P99.9 sojourn still meets the budget.
#
# 8 workers = request-level parallelism.
# ISO-BYTES: equal byte budget of 10 GINX keys is 10 / 20 / 40 slots for GINX / WWL+24 / OURS.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$B/server" || exit 1

BATCH="${BATCH:-1}"
REQS="${REQS:-8000}"
SLO="${SLO:-300}"

declare -A CAP=(   [tfhe]=10        [wwl24]=20        [ours]=40      )

declare -A RATES=( [ours]="40 50 60 62 64" [wwl24]="8 12 16 17 18" [tfhe]="1 2 3 4 5 8" )

for m in ours wwl24 tfhe; do
    v="RATES_$m"
    [ -n "${!v+set}" ] && RATES[$m]="${!v}"
done

TAGBASE="${TAGBASE:-slo8k}"
SKIP_GRID="${SKIP_GRID:-0}"
REPEAT_SPECS="${REPEAT_SPECS:-}"

if [ "$SKIP_GRID" != "1" ]; then
echo "SLO sweep start $(date -Is)  batch=$BATCH reqs=$REQS slo=${SLO}ms"
for m in ours wwl24 tfhe; do
    for r in ${RATES[$m]}; do
        taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
            --workers=8 --cap=${CAP[$m]} --users=50 --rate=$r --requests=$REQS \
            --warm=500 --batch=$BATCH --arrival=poisson --sched=soft \
            --slo=$SLO --seed=1 --devstat=nvme0n1 --tag="${TAGBASE}_r${r}" >/dev/null 2>&1
        f="$OUTDIR/serving_${m/wwl24/wwl+24}_${TAGBASE}_r${r}.json"
        p=$(python3 -c "import json,sys;print(round(json.load(open(sys.argv[1]))['sojourn']['p99_9_ms']))" "$f" 2>/dev/null)
        echo "  $m rate=$r  p99.9=${p:-?}ms  $(date +%H:%M)"
    done
done
echo "SLO sweep done $(date -Is)"
fi

# verify
for spec in $REPEAT_SPECS; do
    m="${spec%%:*}"; rest="${spec#*:}"; r="${rest%%:*}"; n="${rest##*:}"
    for i in $(seq 1 "$n"); do
        taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
            --workers=8 --cap=${CAP[$m]} --users=50 --rate=$r --requests=$REQS \
            --warm=500 --batch=$BATCH --arrival=poisson --sched=soft \
            --slo=$SLO --seed=1 --devstat=nvme0n1 --tag="${TAGBASE}_r${r}_rep${i}" >/dev/null 2>&1
        f="$OUTDIR/serving_${m/wwl24/wwl+24}_${TAGBASE}_r${r}_rep${i}.json"
        python3 - "$f" "$m" "$r" "$i" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); s = d['sojourn']
print("  repeat %s rate=%s #%s  p99.9=%.1f ms  violations=%.3f%%  thr=%.2f req/s" % (
    sys.argv[2], sys.argv[3], sys.argv[4], s['p99_9_ms'],
    100 * d['slo_violation_rate_300'], d['throughput_rps']))
PY
    done
done
