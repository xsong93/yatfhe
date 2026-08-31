#!/usr/bin/env python3
"""Regenerate every paper figure from the result JSONs.

    venv/bin/python make_figures.py                 # out/ -> fig/
    venv/bin/python make_figures.py --sync          # also copy into the paper
    venv/bin/python make_figures.py --out DIR --fig DIR

Figures are written to fig/. The paper consumes them from p1_tc/src/pics/;
--sync copies them across so the two never drift silently.
"""
import argparse
import json
import os
import shutil
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

PAPER_PICS = os.path.expanduser("~/Desktop/Github/p1_tc/src/pics")

# capacity -> pressure ratio, at the 50-tenant sweep configuration
CAPS = [(50, 1.0), (20, 2.5), (10, 5.0), (5, 10.0), (1, 50.0)]
SEEDS = range(1, 6)
METHODS = [
    ("tfhe", "GINX", "#c44e52", "--", "o"),
    ("wwl+24", "WWL+24", "#dd8452", "-.", "s"),
    ("ours", "OURS", "#4c72b0", "-", "^"),
]


def latencies(outdir, method, cap):
    """Pooled per-request latencies in ms across all seeds, sorted."""
    vals = []
    for s in SEEDS:
        path = os.path.join(outdir, f"{method}_benchmark_results_{cap}_s{s}.json")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            vals += [it["time_us"] / 1000.0 for it in json.load(fh)["iterations"]]
    return np.sort(np.array(vals))


def cdf_panels(outdir, figdir):
    """One CDF per pressure ratio. Log abscissa, minor tick labels suppressed --
    five narrow panels cannot carry '2x10^1 3x10^1 ...' legibly."""
    panels = [(50, "fig_cdf_1"), (20, "fig_cdf_2p5"), (10, "fig_cdf_5"),
              (5, "fig_cdf_10"), (1, "fig_cdf_50")]
    for cap, name in panels:
        fig, ax = plt.subplots(figsize=(2.35, 1.42))
        drew = False
        for m, label, colour, style, _ in METHODS:
            v = latencies(outdir, m, cap)
            if not len(v):
                continue
            ax.plot(v, np.arange(1, len(v) + 1) / len(v),
                    label=label, color=colour, ls=style, lw=1.5)
            drew = True
        if not drew:
            plt.close(fig)
            print(f"  skip {name}: no data", file=sys.stderr)
            continue
        ax.set_xscale("log")
        ax.set_ylim(0, 1.03)
        ax.set_xlim(10, 260)
        ax.xaxis.set_major_locator(FixedLocator([10, 30, 100, 200]))
        ax.xaxis.set_major_formatter(FixedFormatter(["10", "30", "100", "200"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel("Latency (ms)", fontsize=7, labelpad=1)
        ax.set_ylabel("CDF", fontsize=7, labelpad=1)
        ax.tick_params(labelsize=6.5, pad=1.5, length=2.5)
        ax.grid(alpha=0.3, lw=0.4)
        ax.set_axisbelow(True)
        ax.legend(fontsize=6, loc="lower right", frameon=False,
                  handlelength=1.6, borderpad=0.2, labelspacing=0.25)
        fig.tight_layout(pad=0.15)
        fig.savefig(os.path.join(figdir, name + ".pdf"))
        plt.close(fig)
        print(f"  {name}.pdf")


def mean_and_jitter(outdir, figdir):
    """Mean latency and standard deviation against pressure, with the crossover
    band shaded. Set the log scale BEFORE the locators -- set_yscale resets
    them, which silently reverts the axis to 10^n labels."""
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.05))
    pressures = [p for _, p in CAPS]
    for m, label, colour, style, marker in METHODS:
        means, sds = [], []
        for cap, _ in CAPS:
            v = latencies(outdir, m, cap)
            if not len(v):
                means, sds = [], []
                break
            means.append(v.mean())
            sds.append(v.std(ddof=1))
        if not means:
            print(f"  skip {label} in fig_avg_analysis: no data", file=sys.stderr)
            continue
        axs[0].plot(pressures, means, style, color=colour, marker=marker,
                    ms=4, lw=1.5, label=label)
        axs[1].plot(pressures, sds, style, color=colour, marker=marker,
                    ms=4, lw=1.5, label=label)

    for ax, ylabel in zip(axs, ["Mean latency (ms)", "Std. deviation (ms)"]):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_xlabel("Pressure ratio (tenants / cache capacity)", fontsize=8)
        ax.xaxis.set_major_locator(FixedLocator([1, 2.5, 5, 10, 50]))
        ax.xaxis.set_major_formatter(
            FixedFormatter(["1", "2.5", "5", "10", "50"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.tick_params(labelsize=7, pad=2)
        ax.grid(alpha=0.3, lw=0.4)
        ax.set_axisbelow(True)
        ax.legend(fontsize=7, frameon=False, loc="lower right", handlelength=2)

    axs[0].yaxis.set_major_locator(FixedLocator([15, 20, 30, 50, 80, 120, 180]))
    axs[0].yaxis.set_major_formatter(
        FixedFormatter(["15", "20", "30", "50", "80", "120", "180"]))
    axs[1].yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20, 50, 80]))
    axs[1].yaxis.set_major_formatter(
        FixedFormatter(["0.5", "1", "2", "5", "10", "20", "50", "80"]))
    for ax in axs:
        ax.yaxis.set_minor_locator(NullLocator())

    axs[0].axvspan(1.0, 2.5, color="0.86", zorder=0)
    axs[0].annotate("crossover", xy=(1.58, axs[0].get_ylim()[1] * 0.94),
                    ha="center", va="top", fontsize=6.5, color="0.3")
    fig.tight_layout(pad=0.3)
    fig.savefig(os.path.join(figdir, "fig_avg_analysis.pdf"))
    plt.close(fig)
    print("  fig_avg_analysis.pdf")


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--out", default=os.path.join(here, "out"),
                    help="directory holding result JSONs")
    ap.add_argument("--fig", default=os.path.join(here, "fig"),
                    help="directory to write figures into")
    ap.add_argument("--sync", action="store_true",
                    help=f"also copy figures into {PAPER_PICS}")
    args = ap.parse_args()

    os.makedirs(args.fig, exist_ok=True)
    print(f"reading {args.out}")
    cdf_panels(args.out, args.fig)
    mean_and_jitter(args.out, args.fig)

    if args.sync:
        os.makedirs(PAPER_PICS, exist_ok=True)
        n = 0
        for f in sorted(os.listdir(args.fig)):
            if f.endswith(".pdf"):
                shutil.copy2(os.path.join(args.fig, f), os.path.join(PAPER_PICS, f))
                n += 1
        print(f"synced {n} figures to {PAPER_PICS}")


if __name__ == "__main__":
    main()
