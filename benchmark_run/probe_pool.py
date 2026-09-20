#!/usr/bin/env python3
"""Probe the cold-read bandwidth of every key in the serving pool.

    probe_pool.py [--server DIR] [--out out/pool_readband.json]
"""
import argparse
import glob
import json
import os
import re
import statistics
import sys
import time

FAMILIES = ["BSK_GINX", "BSK_WWL+24", "BSK_LAZY"]
FRESH = 1.6 # GB/s: above this a file is still in the fast regime


def probe(path):
    """Sequential cold read of the whole file, in GB/s."""
    fd = os.open(path, os.O_RDONLY)
    os.fdatasync(fd)
    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    os.close(fd)
    fd = os.open(path, os.O_RDONLY)
    t0 = time.perf_counter()
    n = 0
    while True:
        b = os.read(fd, 1 << 22)
        if not b:
            break
        n += len(b)
    ms = (time.perf_counter() - t0) * 1000.0
    os.close(fd)
    return n / ms / 1e6


def drop(path):
    fd = os.open(path, os.O_RDONLY)
    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    os.close(fd)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--server", default=os.path.join(here, "server"))
    ap.add_argument("--out", default=os.path.join(here, "out", "pool_readband.json"))
    ap.add_argument("--families", default=",".join(FAMILIES))
    args = ap.parse_args()

    report = {"schema_version": 1, "fresh_threshold_gbps": FRESH,
              "method": "cold sequential read, page cache dropped before and after",
              "families": {}}
    for fam in args.families.split(","):
        files = glob.glob(os.path.join(args.server, fam + "_*.bin"))
        files.sort(key=lambda p: int(re.search(r"_(\d+)\.bin$", p).group(1)))
        if not files:
            print(f"  {fam}: no files", file=sys.stderr)
            continue
        rates = []
        for f in files:
            rates.append(round(probe(f), 3))
        for f in files:
            drop(f)                       # leave no pages behind
        fresh = [i for i, r in enumerate(rates) if r > FRESH]
        report["families"][fam] = {
            "files": len(files),
            "min_gbps": min(rates),
            "median_gbps": round(statistics.median(rates), 3),
            "max_gbps": max(rates),
            "fresh_files": len(fresh),
            "first_id": int(re.search(r"_(\d+)\.bin$", files[0]).group(1)),
            "last_id": int(re.search(r"_(\d+)\.bin$", files[-1]).group(1)),
        }
        r = report["families"][fam]
        print("  %-12s %3d files  ids %d-%d  cold read %.2f / %.2f / %.2f GB/s "
              "(min/median/max), %d above %.1f"
              % (fam, r["files"], r["first_id"], r["last_id"], r["min_gbps"],
                 r["median_gbps"], r["max_gbps"], r["fresh_files"], FRESH))
    allmin = min(v["min_gbps"] for v in report["families"].values())
    allmax = max(v["max_gbps"] for v in report["families"].values())
    report["band_gbps"] = [allmin, allmax]
    report["fresh_files_total"] = sum(v["fresh_files"] for v in report["families"].values())
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("  band %.2f-%.2f GB/s over %d files, %d still in the fast regime -> %s"
          % (allmin, allmax, sum(v["files"] for v in report["families"].values()),
             report["fresh_files_total"], args.out))
    return 1 if report["fresh_files_total"] else 0


if __name__ == "__main__":
    sys.exit(main())
