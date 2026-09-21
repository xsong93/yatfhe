#!/bin/bash
# Bursty-arrival SLO-capacity sweep: on-phase 4x rate for 200 ms, off-phase 0.25x for 800 ms

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$B/server" || exit 1

declare -A CAP=( [tfhe]=10 [wwl24]=20 [ours]=40 )
declare -A FNAME=( [tfhe]=tfhe [wwl24]=wwl+24 [ours]=ours )

declare -A RATES=( [ours]="8 12 16 20 25 28 30 32" [wwl24]="3 4 6" [tfhe]="1" )
for m in ours wwl24 tfhe; do
    v="BRATES_$m"
    [ -n "${!v+set}" ] && RATES[$m]="${!v}"
done

failures=0
echo "burst SLO sweep start $(date -Is)  burst=4x/200ms-on/800ms-off"
for m in ours wwl24 tfhe; do
    for r in ${RATES[$m]}; do
        timeout 12000 taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
            --workers=8 --cap=${CAP[$m]} --users=50 --rate=$r --requests=8000 \
            --warm=500 --batch=1 --arrival=bursty --burston=200 --burstoff=800 \
            --burstfactor=4 --sched=soft --slo=300 --seed=1 --devstat=nvme0n1 --tag="slo_burst_r${r}" \
            >/dev/null 2>&1
        rc=$?
        f="$OUTDIR/serving_${FNAME[$m]}_slo_burst_r${r}.json"
        if [ "$rc" -ne 0 ]; then
            # 124 is the timeout, 137 the kill that follows it. Either way the cell has no result
            echo "  FAILED $m rate=$r (exit $rc: $( [ "$rc" = 124 ] && echo 'timed out after 12000s' || echo 'error' ))"
            echo "         no artifact written: $( [ -f "$f" ] && echo 'file exists from an earlier run, treat it as stale' || echo 'file absent' )"
            failures=$((failures + 1))
            continue
        fi
        p=$(python3 -c "import json,sys;print(round(json.load(open(sys.argv[1]))['sojourn']['p99_9_ms']))" "$f" 2>/dev/null)
        if [ -z "$p" ]; then
            echo "  FAILED $m rate=$r (exited 0 but wrote no readable $f)"
            failures=$((failures + 1))
            continue
        fi
        echo "  $m rate=$r  p99.9=${p}ms  $(date +%H:%M)"
    done
done
if [ "$failures" -ne 0 ]; then
    echo "burst SLO sweep INCOMPLETE: $failures point(s) produced no result $(date -Is)" >&2
    exit 1
fi
echo "burst SLO sweep done $(date -Is)"
