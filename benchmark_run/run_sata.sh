#!/bin/bash

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
MNT="${MNT:-/mnt/sata}"
DST="$MNT/exp1"
SRC="$B/server"
USERS=50

# 1. wait for the machine to be free
while [ -n "$(ps -C blindrotate_cache -o pid= --no-headers 2>/dev/null)" ] \
   || [ -n "$(ps -C blindrotate_comp  -o pid= --no-headers 2>/dev/null)" ]; do sleep 30; done
echo "machine free $(date -Is)"

# 2. wait for a writable mount
echo "waiting for a writable $MNT ..."
while true; do
    if mountpoint -q "$MNT" 2>/dev/null && touch "$MNT/.wtest" 2>/dev/null; then
        rm -f "$MNT/.wtest"; break
    fi
    sleep 20
done
echo "$MNT is mounted and writable $(date -Is)"
df -hT "$MNT" | tail -1

# 3. stage the 50-tenant key set
mkdir -p "$DST" || exit 1
need=$(( ( $(stat -c%s "$SRC/BSK_GINX_1.bin") + $(stat -c%s "$SRC/BSK_LAZY_1.bin") \
        + $(stat -c%s "$SRC/BSK_WWL+24_1.bin") ) * USERS / 1048576 ))
avail=$(df -m --output=avail "$MNT" | tail -1)
echo "need ~${need} MB, ${avail} MB available"
[ "$avail" -lt "$((need + 2048))" ] && { echo "not enough space"; exit 1; }

echo "copying keys $(date -Is)"
for i in $(seq 1 $USERS); do
    for p in GINX LAZY "WWL+24"; do cp "$SRC/BSK_${p}_${i}.bin" "$DST/" || exit 1; done
done
cp "$SRC/blindrotate_cache" "$DST/"
sync
echo "copied $(ls "$DST"/BSK_*.bin | wc -l) keys, $(du -sh --apparent-size "$DST" | cut -f1)"

# 4. run
cd "$DST" || exit 1
echo "benchmark start $(date -Is)"
taskset -c 0-7 ./blindrotate_cache --out="$OUTDIR" --caps=10 --users=$USERS --reqs=1000 \
    --warm=500 --seed=1 --tag=sata 2>&1 | grep -E "^Saved|pages still resident"
echo "benchmark done $(date -Is)"

# 5. leave the disk as we found it
cd "$B" && rm -rf "$DST"
echo "staged keys removed from $MNT"
