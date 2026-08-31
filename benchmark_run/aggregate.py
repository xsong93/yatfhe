#!/usr/bin/env python3
"""Aggregate blindrotate_cache results (schema_version 2) into a report.

Pools every seed found for a given (method, pressure ratio) and reports
percentiles with bootstrap confidence intervals, so the tail numbers carry an
uncertainty instead of being a single unqualified sample.

    python3 aggregate.py server/ --format markdown
    python3 aggregate.py server/ --format latex --baseline OURS
"""
import argparse
import glob
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

METHOD_ORDER = ["TFHE", "WWL+24", "OURS"]


def load(paths):
    runs = defaultdict(lambda: {"times": [], "hits": [], "seeds": set(), "cfg": None})
    for path in paths:
        with open(path) as fh:
            doc = json.load(fh)
        if doc.get("schema_version") != 2:
            print(f"skipping {path}: schema_version "
                  f"{doc.get('schema_version', 'missing')}, expected 2", file=sys.stderr)
            continue
        cfg = doc["config"]
        method = doc["benchmark_name"].split("/")[-1]
        key = (method, round(cfg["pressure_ratio"], 4), cfg["iso_mode"])
        entry = runs[key]
        entry["cfg"] = cfg
        entry["seeds"].add(cfg["seed"])
        for it in doc["iterations"]:
            entry["times"].append(it["time_us"] / 1000.0)
            entry["hits"].append(bool(it["hit"]))
    return runs


def bootstrap_ci(sample, stat, reps=2000, alpha=0.05, rng=None):
    """Percentile bootstrap CI. Returns (lo, hi)."""
    rng = rng or np.random.default_rng(12345)
    n = len(sample)
    if n < 2:
        return (float("nan"), float("nan"))
    idx = rng.integers(0, n, size=(reps, n))
    stats = stat(np.asarray(sample)[idx], axis=1)
    return (float(np.percentile(stats, 100 * alpha / 2)),
            float(np.percentile(stats, 100 * (1 - alpha / 2))))


def describe(times, hits, seeds):
    t = np.asarray(times, dtype=float)
    out = {
        "n": t.size,
        "seeds": len(seeds),
        "hit_rate": float(np.mean(hits)) if len(hits) else float("nan"),
        "mean": float(t.mean()),
        "sd": float(t.std(ddof=0)),
        "p50": float(np.percentile(t, 50)),
        "p99": float(np.percentile(t, 99)),
        "p999": float(np.percentile(t, 99.9)),
        "max": float(t.max()),
        # how many samples actually sit above p99.9 -- if this is <10 the
        # tail estimate is not resolvable and the CI will say so
        "tail_support": int((t >= np.percentile(t, 99.9)).sum()),
    }
    out["mean_ci"] = bootstrap_ci(t, lambda a, axis: a.mean(axis=axis))
    out["sd_ci"] = bootstrap_ci(t, lambda a, axis: a.std(axis=axis, ddof=0))
    out["p999_ci"] = bootstrap_ci(t, lambda a, axis: np.percentile(a, 99.9, axis=axis))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory", nargs="?", default="server")
    ap.add_argument("--format", choices=["markdown", "latex"], default="markdown")
    ap.add_argument("--baseline", default="OURS",
                    help="method used as the 1.0x reference in ratio columns")
    ap.add_argument("--iso", default=None, help="only report this cache budget model")
    args = ap.parse_args()

    # Match ONLY the seeded pressure sweep. A bare *_benchmark_results_*.json
    # also catches the thread runs (cap=100, tag thr_*) and the tenant runs
    # (caps 5/10/20/40, tags tenA_*/tenB_*), which share capacities with the
    # sweep and so land in the same pressure-ratio buckets: that silently
    # inflated PR 0.5 to n=6800, PR 5 to n=17000, and invented a PR 20 row.
    pat = re.compile(r"_benchmark_results_\d+_s\d+\.json$")
    paths = sorted(p for p in glob.glob(
        os.path.join(args.directory, "*_benchmark_results_*.json"))
        if pat.search(os.path.basename(p)))
    if not paths:
        sys.exit(f"no seeded sweep results under {args.directory} "
                 f"(expecting <method>_benchmark_results_<cap>_s<seed>.json)")
    runs = load(paths)
    if not runs:
        sys.exit("no schema_version 2 results found")

    stats = {}
    for (method, pr, iso), entry in runs.items():
        if args.iso and iso != args.iso:
            continue
        stats[(method, pr, iso)] = describe(entry["times"], entry["hits"], entry["seeds"])

    isos = sorted({k[2] for k in stats})
    for iso in isos:
        prs = sorted({k[1] for k in stats if k[2] == iso})
        print(f"\n## cache budget model: {iso}\n")
        if args.format == "markdown":
            print("| PR | method | n | seeds | HR | avg (ms) | sd (ms) | P50 | P99 | "
                  "P99.9 [95% CI] | tail n |")
            print("|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|")
        for pr in prs:
            base = stats.get((args.baseline, pr, iso))
            for method in METHOD_ORDER:
                s = stats.get((method, pr, iso))
                if not s:
                    continue
                if args.format == "markdown":
                    rel = ""
                    if base and method != args.baseline and base["sd"] > 0:
                        rel = f" ({s['sd']/base['sd']:.1f}x)"
                    print(f"| {pr:g} | {method} | {s['n']} | {s['seeds']} | {s['hit_rate']:.2f} "
                          f"| {s['mean']:.1f} | {s['sd']:.1f}{rel} | {s['p50']:.1f} "
                          f"| {s['p99']:.1f} | {s['p999']:.1f} "
                          f"[{s['p999_ci'][0]:.1f}, {s['p999_ci'][1]:.1f}] "
                          f"| {s['tail_support']} |")
                else:
                    print(f"{pr:g} & {method} & {s['hit_rate']:.2f} & {s['mean']:.1f} & "
                          f"{s['sd']:.1f} & {s['p50']:.1f} & {s['p99']:.1f} & "
                          f"{s['p999']:.1f} \\\\")

        # crossover: lowest pressure at which the baseline beats every other method
        print()
        for pr in prs:
            base = stats.get((args.baseline, pr, iso))
            if not base:
                continue
            others = {m: stats[(m, pr, iso)]["mean"] for m in METHOD_ORDER
                      if m != args.baseline and (m, pr, iso) in stats}
            if not others:
                continue
            verdict = "wins" if all(base["mean"] < v for v in others.values()) else "loses"
            detail = ", ".join(f"{m} {v/base['mean']:.2f}x" for m, v in others.items())
            print(f"  PR {pr:>5g}: {args.baseline} {verdict:>5}  ({detail})")

    under = [k for k, s in stats.items() if s["tail_support"] < 10]
    if under:
        print(f"\nwarning: {len(under)} series have fewer than 10 samples at or above "
              f"P99.9; that percentile is not resolvable from them.", file=sys.stderr)


if __name__ == "__main__":
    main()
