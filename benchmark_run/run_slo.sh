#!/bin/bash
# SLO-capacity sweep (R1-S2): the highest offered rate P99.9 sojourn still meets the budget.
# Shows latency, tail behaviour and throughput together in one number.
#
# *** 8000 REQUESTS PER RATE POINT. ***
# P99.9 of 1000 samples IS the single worst observation.
#
# *** PER-METHOD RATE GRIDS. ***
# Each method's crossing is bracketed rather than swept with one wide grid. Wall time per point is requests/rate,
# so low-rate points are the expensive ones.
#
# 8 workers = request-level parallelism (all methods peak there).
# Caps are ISO-BYTES: --cap is a literal slot count, so an equal byte budget of 10 GINX keys is 10 / 20 / 40 slots for GINX / WWL+24 / OURS.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$B/server" || exit 1

BATCH="${BATCH:-1}"
REQS="${REQS:-8000}"
SLO="${SLO:-300}"

declare -A CAP=(   [tfhe]=10        [wwl24]=20        [ours]=40      )

# override with e.g. RATES_ours="30 40 50"
declare -A RATES=( [ours]="40 50 60" [wwl24]="16 12 8" [tfhe]="8 5" )

for m in ours wwl24 tfhe; do
    v="RATES_$m"
    [ -n "${!v+set}" ] && RATES[$m]="${!v}"
done

TAGBASE="${TAGBASE:-slo8k}"

echo "SLO sweep start $(date -Is)  batch=$BATCH reqs=$REQS slo=${SLO}ms"
for m in ours wwl24 tfhe; do
    for r in ${RATES[$m]}; do
        taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
            --workers=8 --cap=${CAP[$m]} --users=50 --rate=$r --requests=$REQS \
            --warm=500 --batch=$BATCH --arrival=poisson --sched=soft \
            --slo=$SLO --seed=1 --tag="${TAGBASE}_r${r}" >/dev/null 2>&1
        f="$OUTDIR/serving_${m/wwl24/wwl+24}_${TAGBASE}_r${r}.json"
        p=$("$B/venv/bin/python" -c "import json,sys;print(round(json.load(open(sys.argv[1]))['sojourn']['p99_9_ms']))" "$f" 2>/dev/null)
        echo "  $m rate=$r  p99.9=${p:-?}ms  $(date +%H:%M)"
    done
done
echo "SLO sweep done $(date -Is)"
