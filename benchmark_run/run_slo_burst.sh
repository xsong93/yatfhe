#!/bin/bash
# Bursty-arrival SLO-capacity sweep (R1-S2): on-phase 4x rate for 200 ms, off-phase 0.25x for 800 ms

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$B/server" || exit 1

declare -A CAP=( [tfhe]=10 [wwl24]=20 [ours]=40 )
declare -A FNAME=( [tfhe]=tfhe [wwl24]=wwl+24 [ours]=ours )

# override with e.g. BRATES_ours="50 55 60 62"
declare -A RATES=( [ours]="25 30 35 40 50 55 60 62" [wwl24]="6 7 8 12 14 16" [tfhe]="" )
for m in ours wwl24 tfhe; do
    v="BRATES_$m"
    [ -n "${!v+set}" ] && RATES[$m]="${!v}"
done

echo "burst SLO sweep start $(date -Is)  burst=4x/200ms-on/800ms-off"
for m in ours wwl24 tfhe; do
    for r in ${RATES[$m]}; do
        taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
            --workers=8 --cap=${CAP[$m]} --users=50 --rate=$r --requests=8000 \
            --warm=500 --batch=1 --arrival=bursty --burston=200 --burstoff=800 \
            --burstfactor=4 --sched=soft --slo=300 --seed=1 --tag="slo_burst_r${r}" \
            >/dev/null 2>&1
        f="$OUTDIR/serving_${FNAME[$m]}_slo_burst_r${r}.json"
        p=$("$B/venv/bin/python" -c "import json,sys;print(round(json.load(open(sys.argv[1]))['sojourn']['p99_9_ms']))" "$f" 2>/dev/null)
        echo "  $m rate=$r  p99.9=${p:-?}ms  $(date +%H:%M)"
    done
done
echo "burst SLO sweep done $(date -Is)"
