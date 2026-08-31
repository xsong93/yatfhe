#!/bin/bash
# Re-run the N2 (USB 3.1 Gen 2) tier on the derived-first-key build.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
MNT=/media/spark/N2
DST="$MNT/exp1"
SRC="$B/server"
USERS=50
while [ -n "$(ps -C blindrotate_cache -o pid= --no-headers 2>/dev/null)" ]; do sleep 20; done
rm -rf "$DST"; mkdir -p "$DST" || exit 1
echo "copying keys $(date -Is)"
for i in $(seq 1 $USERS); do
    for p in GINX LAZY "WWL+24"; do cp "$SRC/BSK_${p}_${i}.bin" "$DST/" || exit 1; done
done
cp "$SRC/blindrotate_cache" "$DST/"; sync
echo "staged $(ls "$DST"/BSK_*.bin | wc -l) keys ($(du -sh --apparent-size "$DST" | cut -f1))"
cd "$DST" || exit 1
echo "benchmark start $(date -Is)"
taskset -c 0-7 ./blindrotate_cache --out="$OUTDIR" --caps=10 --users=$USERS --reqs=1000 \
    --warm=500 --seed=1 --tag=n2 2>&1 | grep -E "^Saved|pages still resident"
cd "$B" && rm -rf "$DST"
echo "benchmark done $(date -Is)"
