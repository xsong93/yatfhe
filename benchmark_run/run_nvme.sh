#!/bin/bash

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
cd "$B/server" || exit 1
while [ -n "$(ps -C blindrotate_cache -o pid= --no-headers 2>/dev/null)" ] \
   || [ -n "$(ps -C blindrotate_comp  -o pid= --no-headers 2>/dev/null)" ]; do sleep 30; done
echo "NVMe tier point start $(date -Is)  load $(cut -d' ' -f1 /proc/loadavg)"
taskset -c 0-7 ./blindrotate_cache --out="$OUTDIR" --caps=10 --users=50 --reqs=1000 \
    --warm=500 --seed=1 --tag=nvme 2>&1 | grep -E "^Saved|pages still resident" || true
rc=${PIPESTATUS[0]}
for m in tfhe wwl+24 ours; do
    f="$OUTDIR/${m}_benchmark_results_10_nvme.json"
    if [ ! -f "$f" ]; then echo "  MISSING $f (exit $rc)"; exit 1; fi
    python3 - "$f" "$m" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); s = d['statistics']
print("  %-7s mean %7.2f ms  hit %.3f  cold %6.2f ms  n=%d" % (
    sys.argv[2], s['mean_us'] / 1000, s['hit_rate'],
    d.get('miss_statistics', {}).get('mean_us', float('nan')) / 1000, s['count']))
PY
done
echo "NVMe tier point done $(date -Is)"
