#!/usr/bin/env python3
"""Generate the decision-tree depth table (supplement Table S4) from
trees_out/serving_*_tree_b*.json. Run after run_trees.sh completes.

Usage:  python3 trees_table.py [--sync]
"""
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, "out")
OUT = "/home/spark/Desktop/Github/p1_tc/src/supplement/tab_trees.tex"

METHODS = [("tfhe", "GINX"), ("wwl24", "WWL+24"), ("ours", "OURS")]
BATCHES = [1, 2, 4, 8]
DEPTH = {1: 1, 2: 2, 4: 3, 8: 4}


def fmt(x, nd=1):
    return f"{x:.{nd}f}"


def load(method, b):
    stem = {"wwl24": "wwl+24"}.get(method, method)
    p = os.path.join(DATA, f"serving_{stem}_tree_b{b}.json")
    if not os.path.exists(p):
        return None
    d = json.load(open(p, encoding="utf-8"))
    return d


def main():
    rows = []
    for b in BATCHES:
        for stem, name in METHODS:
            d = load(stem, b)
            if d is None:
                print(f"missing serving_{stem}_tree_b{b}.json", file=sys.stderr)
                continue
            s, q = d["sojourn"], d["queue_delay"]
            rows.append((DEPTH[b], b, name, d["hit_rate"], s["mean_ms"], s["p50_ms"],
                         s["p99_ms"], s["p99_9_ms"], q["mean_ms"],
                         d["slo_violation_rate"]))
    tex = []
    tex.append(r"\begin{table}[!ht]")
    tex.append(r"    \centering")
    tex.append(r"    \caption{End-to-end query latency of a depth-$d$ decision-tree query, modeled as $B = 2^{d-1}$ sequential blind rotations on one tenant key, the bootstrap count of such a query. The homomorphic additions between bootstraps are omitted because they are cheap, key-independent and identical across schemes. Offered rate 8 req/s, soft affinity, pressure ratio 5.0, 8000 requests per configuration.}")
    tex.append(r"    \label{tab:trees}")
    tex.append(r"    \renewcommand{\arraystretch}{1.15}")
    tex.append(r"    \setlength{\tabcolsep}{3pt}")
    tex.append(r"    \begin{tabular}{lrrrrrrrr}")
    tex.append(r"        \toprule")
    tex.append(r"        Depth & $B$ & Method & Hit rate & Mean & P50 & P99 & P99.9 & SLO viol.")
    tex.append(r"        \\")
    tex.append(r"        \midrule")
    for d, b, name, hr, mean, p50, p99, p999, qm, viol in rows:
        tex.append(
            f"        {d} & {b} & {name} & {fmt(hr, 2)} & {fmt(mean, 0)} & {fmt(p50, 0)} & "
            f"{fmt(p99, 0)} & {fmt(p999, 0)} & {fmt(viol * 100, 1)}\\% \\\\")
    tex.append(r"        \bottomrule")
    tex.append(r"    \end{tabular}")
    tex.append(r"    \tabnote{Latency columns are sojourn (queueing included) in ms; SLO violation rate is the fraction of queries exceeding the 300\,ms budget.}")
    tex.append(r"\end{table}")
    content = "\n".join(tex) + "\n"

    if "--sync" in sys.argv:
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        with open(OUT, "w", encoding="utf-8") as f:
            f.write(content)
        print("wrote", OUT)
    else:
        print(content)
        # also print a compact summary for prose writing
        print("\nSUMMARY (depth, B, method, HR, mean, p50, p99, p99.9, viol%):")
        for d, b, name, hr, mean, p50, p99, p999, qm, viol in rows:
            print(f"  d={d} B={b} {name:7s} HR={hr:.2f} mean={mean:.1f} "
                  f"p50={p50:.1f} p99={p99:.1f} p99.9={p999:.1f} viol={viol*100:.1f}%")


if __name__ == "__main__":
    main()
