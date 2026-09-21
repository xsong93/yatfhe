#!/bin/bash
# Saturated throughput against tenant population, at one equal memory budget.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
mkdir -p "$OUTDIR"
cd "$B/server" || exit 1
RATE="${RATE:-1000}"
REQS="${REQS:-1000}"
STEP="${STEP:-5}"
STEP_HI_START="${STEP_HI_START:-200}" # coarser spacing beyond this population
STEP_HI="${STEP_HI:-25}"
STOP_BELOW="${STOP_BELOW:-30}" # the load the ceiling experiment offers
STOP_EXTRA="${STOP_EXTRA:-2}" # points kept after the first crossing
declare -A CAP=(   [tfhe]=100 [wwl24]=200 [ours]=400 )
declare -A FNAME=( [tfhe]=tfhe [wwl24]=wwl+24 [ours]=ours )
declare -A MAXU=( [tfhe]=800 [wwl24]=800 [ours]=800 ) # Population range per scheme
for m in tfhe wwl24 ours; do
    v="MAXU_$m"; [ -n "${!v+set}" ] && MAXU[$m]="${!v}"
done
# STARTU_<method> resumes a curve at a later population, for extending a curve
declare -A STARTU=( [tfhe]=1 [wwl24]=1 [ours]=1 )
for m in tfhe wwl24 ours; do
    v="STARTU_$m"; [ -n "${!v+set}" ] && STARTU[$m]="${!v}"
done
METHODS="${METHODS:-tfhe wwl24 ours}"

failures=0
echo "saturation curve start $(date -Is)  rate=$RATE req/s  step=$STEP out=$OUTDIR"
for m in $METHODS; do
    u=${STARTU[$m]}
    below=0
    while [ "$u" -le "${MAXU[$m]}" ]; do
        tag="sat_u${u}"
        f="$OUTDIR/serving_${FNAME[$m]}_${tag}.json"
        timeout 1800 taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m \
            --workers=8 --cap=${CAP[$m]} --users="$u" --rate="$RATE" \
            --requests="$REQS" --warm=500 --batch=1 --arrival=poisson --sched=soft \
            --slo=300 --seed=1 --devstat=nvme0n1 --tag="$tag" >/dev/null 2>&1
        rc=$?
        if [ "$rc" -ne 0 ] || [ ! -f "$f" ]; then
            echo "  FAILED $m users=$u (exit $rc)"
            failures=$((failures + 1))
        else
            read -r thr stop <<EOF
$(python3 - "$f" "$m" "$u" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print("%7.2f %.3f" % (d['throughput_rps'], d['hit_rate']))
PY
)
EOF
            python3 - "$f" "$m" "$u" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print("  %-7s users=%-4s thr %7.2f req/s  service %6.2f ms  hit %.3f  device inflight %.2f  wall %5.1fs"
      % (sys.argv[2], sys.argv[3], d['throughput_rps'], d['service']['mean_ms'],
         d['hit_rate'], d['device_inflight_reads_mean'], d['wall_sec']))
PY
            if python3 -c "import sys; sys.exit(0 if float('$thr') < float('$STOP_BELOW') else 1)"; then
                below=$((below + 1))
                if [ "$below" -gt "$STOP_EXTRA" ]; then
                    echo "  $m stays below $STOP_BELOW req/s; stopping at $u tenants"
                    break
                fi
            fi
        fi
        if [ "$u" -ge "$STEP_HI_START" ]; then step="$STEP_HI"; else step="$STEP"; fi
        next=$((u + step))
        # land the last point exactly on the limit
        if [ "$next" -gt "${MAXU[$m]}" ]; then
            if [ "$u" -ge "${MAXU[$m]}" ]; then break; fi
            next="${MAXU[$m]}"
        fi
        u=$next
    done
done
if [ "$failures" -ne 0 ]; then
    echo "saturation curve INCOMPLETE: $failures cell(s) missing $(date -Is)" >&2
    exit 1
fi
echo "saturation curve done $(date -Is)"
