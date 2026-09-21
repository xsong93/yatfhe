#!/bin/bash
# Age freshly written key files on the SSD before measuring cold reads.
#
# Usage: fold_device.sh [GB per round] [probe file] [target GB/s] [max rounds]

set -u
GB="${1:-32}"
DIR="$(cd "$(dirname "$0")" && pwd)"
PROBE="${2:-$DIR/server/COMP_A_GINX.bin}"
TARGET="${3:-1.6}"
ROUNDS="${4:-8}"

probe() {
    python3 - "$PROBE" <<'PY'
import os, sys, time
f = sys.argv[1]
fd = os.open(f, os.O_RDONLY)
os.fdatasync(fd)
os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
os.close(fd)
fd = os.open(f, os.O_RDONLY)
t0 = time.perf_counter()
n = 0
while True:
    b = os.read(fd, 1 << 22)
    if not b:
        break
    n += len(b)
ms = (time.perf_counter() - t0) * 1000
os.close(fd)
print("%.2f" % (n / ms / 1e6))
PY
}

if [ ! -f "$PROBE" ]; then
    echo "probe file $PROBE does not exist; nothing to fold"
    exit 1
fi
for round in $(seq 1 "$ROUNDS"); do
    gbps="$(probe)"
    echo "round $round: $(basename "$PROBE") reads at ${gbps} GB/s (target <= ${TARGET})"
    if python3 -c "import sys; sys.exit(0 if float('$gbps') <= float('$TARGET') else 1)"; then
        echo "device is in the aged regime; safe to measure cold reads"
        exit 0
    fi
    echo "  writing ${GB} GB through the device to evict the fast cache"
    dd if=/dev/zero of="$DIR/.fold.$$" bs=1M count=$((GB * 1024)) oflag=direct 2>&1 | tail -1
    sync
    rm -f "$DIR/.fold.$$"
    sync
    sleep 60
done
echo "still reading at ${gbps} GB/s after ${ROUNDS} rounds: the keys are not aged,"
echo "so a cold column measured now would be optimistic"
exit 1
