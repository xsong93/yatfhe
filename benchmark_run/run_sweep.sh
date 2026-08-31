#!/bin/bash
# Full cache-pressure sweep

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"

CPUS="${CPUS:-0-7}"
REQS="${REQS:-1000}"
WARM="${WARM:-500}"
CAPS="${CAPS:-100,50,20,10,5,1}"
SEEDS="${SEEDS:-1 2 3 4 5}"
ISO="${ISO:-bytes}"

cd "$(dirname "$0")/server" || exit 1

echo "sweep start $(date -Is)  cpus=$CPUS reqs=$REQS caps=$CAPS seeds='$SEEDS' iso=$ISO"
for seed in $SEEDS; do
    echo "=== seed $seed  $(date -Is) ==="
    taskset -c "$CPUS" ./blindrotate_cache --out="$OUTDIR" \
        --caps="$CAPS" --reqs="$REQS" --warm="$WARM" \
        --seed="$seed" --iso="$ISO" --tag="s${seed}" \
        2>&1 | grep -E "^affinity|^===|^Saved|^warning|^error|parameter mismatch"
    status=${PIPESTATUS[0]}
    if [ "$status" -ne 0 ]; then
        echo "seed $seed failed with status $status, stopping"
        exit "$status"
    fi
done
echo "sweep done $(date -Is)"
