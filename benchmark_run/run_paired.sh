#!/bin/bash
# Paired micro-benchmarks on fixed inputs.
#
# The first-step A/B on one fixed key,
# the dispatch-width and NTT-split ablations,
# the stage ratio T_SS'/T_EP,
# and the caller-side derivation cost.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"

mkdir -p "$OUTDIR/paired" "$OUTDIR/paired_fk"
taskset -c 0-7 "$B/server/firststep_fixedkey" --reps=200 \
    --keyfile="$OUTDIR/paired_fk/fs_key.bin" --out="$OUTDIR/paired_fk/seed" >/dev/null 2>&1
for i in 1 2 3 4 5 6 7 8; do
    taskset -c 0-7 "$B/server/firststep_fixedkey" --reps=5000 \
        --keyfile="$OUTDIR/paired_fk/fs_key.bin" --out="$OUTDIR/paired_fk/run$i" >/dev/null 2>&1
done
for i in 1 2 3 4; do
    taskset -c 0-7 "$B/server/stagea_width_ablate" --reps=2000 --out="$OUTDIR/paired/sw$i" >/dev/null 2>&1
done
for i in 1 2; do
    taskset -c 0-7 "$B/server/ntt_split_ablate" --reps=1200 --out="$OUTDIR/paired/ns$i" >/dev/null 2>&1
done

taskset -c 0-7 "$B/server/stage_ratio_final" --runs=10 | tee "$OUTDIR/paired/ratio.txt"
taskset -c 0-7 "$B/server/firststep_derive_timing" --reps=20000 --out="$OUTDIR/paired" >/dev/null 2>&1
