#!/bin/bash
# 4b (R3-4): sensitivity to the number of tenants
#
# Sweep A - cache fixed, tenants increase.
#           Pressure rises with tenants.
#           This is the deployment question: what happens as the service takes on more users without adding memory.
# Sweep B - pressure ratio held at 5.0.
#           This isolates scale effects that are NOT pressure: distinct key files, larger LRU structures, longer Zipf tail.
#           Flat curves mean the design scales; rising curves does not, independently of memory pressure.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$(dirname "$0")/server" || exit 1
REQS="${REQS:-1000}"
SEEDS="${SEEDS:-1 2}"

for seed in $SEEDS; do
    echo "########## seed $seed ##########"
    echo "===== sweep A: cap=10 fixed, tenants rising ====="
    for u in 25 50 100 200; do
        echo "--- users=$u (PR $(echo "scale=1; $u/10" | bc)) ---"
        ./blindrotate_cache --out="$OUTDIR" --caps=10 --users="$u" --reqs="$REQS" --warm=1200 \
            --seed="$seed" --tag="tenA_u${u}_s${seed}" 2>&1 | grep -E "^Saved"
    done
    echo "===== sweep B: PR=5.0 fixed, tenants and cache scale together ====="
    for pair in "25 5" "50 10" "100 20" "200 40"; do
        set -- $pair
        echo "--- users=$1 cap=$2 (PR 5.0) ---"
        ./blindrotate_cache --out="$OUTDIR" --caps="$2" --users="$1" --reqs="$REQS" --warm=1200 \
            --seed="$seed" --tag="tenB_u${1}_s${seed}" 2>&1 | grep -E "^Saved"
    done
done
echo "tenant sweep done $(date -Is)"
