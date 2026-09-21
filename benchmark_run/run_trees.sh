#!/bin/bash
# End-to-end decision-tree depth sweep.
# A depth-d tree query is modeled as B = 2^(d-1) sequential blind rotations on one tenant key.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
mkdir -p "$OUTDIR"
cd "$B/server" || exit 1
declare -A CAP=( [tfhe]=10 [wwl24]=20 [ours]=40 )
declare -A FNAME=( [tfhe]=tfhe [wwl24]=wwl+24 [ours]=ours )
echo "tree-depth sweep start $(date -Is) out=$OUTDIR"
for b in 1 2 4 8; do
  for m in tfhe wwl24 ours; do
    f="$OUTDIR/serving_${FNAME[$m]}_tree_b${b}.json"
    # Resume only when the existing file was produced by THIS protocol
    if python3 - "$f" "$b" "${CAP[$m]}" <<'PY' 2>/dev/null
import json, sys
try:
    c = json.load(open(sys.argv[1]))["config"]
except Exception:
    sys.exit(1)
want = dict(offered_rate_rps=8.0, requests=8000, warm=500, batch=int(sys.argv[2]),
            cache_capacity_keys=int(sys.argv[3]), users=50, workers=8,
            arrival="poisson", sched="soft", seed=1)
sys.exit(0 if all(c.get(k) == v for k, v in want.items()) else 1)
PY
    then
      echo "  [B=$b $m] already complete under this protocol, skipping"
      continue
    fi
    echo "  [B=$b $m] running (existing file absent or from another protocol)"
    taskset -c 0-7 ./blindrotate_serving --out="$OUTDIR" --method=$m --workers=8 --cap=${CAP[$m]} \
        --users=50 --rate=8 --requests=8000 --warm=500 --batch=$b \
        --arrival=poisson --sched=soft --slo=300 --seed=1 --devstat=nvme0n1 \
        --tag="tree_b${b}" 2>&1 | grep -E "^  (throughput|mean|p99|hit)" | sed "s/^/  [B=$b $m] /"
  done
done
echo "tree-depth sweep done $(date -Is)"
