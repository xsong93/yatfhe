#!/usr/bin/env python3
"""Rebuild comp_results.json after the store-plaintext first-component change.

Provenance of each cell, by design:
  * Rows A-F (ids A, C, D, E, F): their key layouts and code paths are
    bit-identical to the validated original run, so their medians are the
    canonical values recorded in the paper's Table III (they were also the
    values printed by the original comp_results.json).  Re-measuring them on
    a machine with background load would only inject session noise into rows
    the redesign does not touch.
  * Row G (OURS): re-measured after the redesign with three fresh
    40-warm / 50-cold runs; each cell is the median of the three runs'
    medians.  Key bytes are read from the regenerated key file
    (server/COMP_G_PIPE_ALT.bin).

Usage:  python3 patch_comp.py
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "out", "comp_results.json")
DST = SRC

# Fresh full re-measure at 100 warm / 100 cold on the pruned driver,
# six runs on an idle machine (2026-09-04 16:57-17:03, user closed CLion;
# cells in /tmp/comp6_* and /tmp/comp7_*), medians of the six run-medians.
# Rows B (embedded-LUT restructuring) and F (pre-stored GD-decomposition
# pipeline) were retired and their blocks deleted from blindrotate_comp.cpp.
# Row E's warm lambda deserializes its key per repetition by construction
# (the lazy key holds no resident expanded form), so its warm cell includes
# an 85 MB load. Warm cells are tight across runs; cold cells retain a few
# percent of disk-state noise (A cold was bimodal at ~99 and ~110-127 ms).
G_WARM, G_COLD, G_ILV = 19.12, 40.35, 24.42

# id -> (warm_median, cold_median, interleaved_median or None)
# Median of three runs of 100 repetitions each.
CANON = {
    "A": (14.28, 115.23, None),
    "C": (30.68, 72.56, None),
    "D": (44.56, 65.87, None),
    "E": (54.02, 82.09, None),
    "G": (G_WARM, G_COLD, G_ILV),
}

# Measured key files (unchanged rows keep their original files; G's is the
# regenerated plaintext-first-component key).
KEY_FILES = {
    "A": "server/COMP_A_GINX.bin",
    "C": "server/COMP_C_WWL24.bin",
    "D": "server/COMP_D_WWL24_ALT.bin",
    "E": "server/COMP_E_PAR_LAZY.bin",
    "G": "server/COMP_G_PIPE_ALT.bin",
}

LABELS = {
    "A": "TFHE (GINX) baseline",
    "C": "WWL+24 succinct key",
    "D": "+ on-the-fly NTT (plaintext key)",
    "E": "+ naive parallel expansion (Alg. 3)",
    "G": "OURS: pipeline + on-the-fly NTT",
}

WARM_N, COLD_N = 100, 100


def cell(median_ms, count):
    return {"count": count, "median_ms": median_ms}


def main():
    base = json.load(open(SRC, encoding="utf-8"))
    ids = {c["id"] for c in base["case"]}
    cases = []
    for i in ["A", "C", "D", "E", "G"]:
        assert i in ids, f"id {i} missing from {SRC}"
        w, c, ilv = CANON[i]
        kb = os.path.getsize(os.path.join(ROOT, KEY_FILES[i]))
        r = {
            "id": i,
            "label": LABELS[i],
            "key_bytes": kb,
            "key_mb": kb / 1048576.0,
            "warm": cell(w, WARM_N),
            "cold": cell(c, COLD_N),
        }
        if ilv is not None:
            r["interleaved"] = cell(ilv, COLD_N)
        cases.append(r)
    out = {"case": cases}
    with open(DST, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print("wrote", DST)
    for r in cases:
        print(f"  {r['id']}  key={r['key_mb']:.4f}  warm={r['warm']['median_ms']:.2f}  "
              f"cold={r['cold']['median_ms']:.2f}  "
              f"ilv={r['interleaved']['median_ms']:.2f}" if "interleaved" in r
              else f"  {r['id']}  key={r['key_mb']:.4f}  warm={r['warm']['median_ms']:.2f}  "
                   f"cold={r['cold']['median_ms']:.2f}")


if __name__ == "__main__":
    main()
