#!/usr/bin/env python3
"""The model behind the saturation curve: what a miss costs, and what bends it.

Inputs: the curve cells of run_satcurve.sh, the ceiling cells of run_tencap.sh and
the pool's measured cold-read rates from probe_pool.py. Output: out/satcurve_fit.json.

    satcurve_model.py [--out out] [--dest out/satcurve_fit.json]
"""
import argparse
import glob
import json
import os
import re
import sys

# key file sizes, as gen_benchkeys writes them (one file per tenant and scheme)
KEY_BYTES = {"tfhe": 178372849, "wwl+24": 89446336, "ours": 44890016}
LABEL = {"tfhe": "GINX", "wwl+24": "WWL+24", "ours": "OURS"}
FAMILY = {"tfhe": "BSK_GINX", "wwl+24": "BSK_WWL+24", "ours": "BSK_LAZY"}
CEILING_FROM = 300 # tenants past which GINX's demand has pinned
BUDGET_MS = 300.0 # the SLO the ceiling cells are measured against
SINGLE_STREAM_FILE = "tfhe_benchmark_results_10_nvme.json" # one cold key read
SINGLE_STREAM_PATH = "miss_statistics.mean_us"
KEY_GINX_BYTES = 178372849
PLATEAU_FROM = 300 # tenants past which the demand has settled (same as the ceiling)


def dig(obj, path):
    """Value at a dotted path, or None."""
    for part in path.split("."):
        if not isinstance(obj, dict) or part not in obj:
            return None
        obj = obj[part]
    return obj


def curve(outdir, method):
    d = {}
    for f in glob.glob(os.path.join(outdir, f"serving_{method}_sat_u*.json")):
        u = int(re.search(r"_sat_u(\d+)\.json$", f).group(1))
        j = json.load(open(f))
        n = j["config"]["requests"]
        misses = (1.0 - j["hit_rate"]) * n
        reads = max(misses - j.get("coalesced_loads", 0), 0.0)
        d[u] = (j["throughput_rps"], 1.0 - j["hit_rate"],
                j["device_inflight_reads_mean"], reads / n,
                j.get("coalesced_loads", 0) / max(misses, 1.0))
    return d


def ceiling_cell(outdir, method, u):
    f = os.path.join(outdir, f"serving_{method}_tencap100_u{u}.json")
    if not os.path.exists(f):
        return None
    j = json.load(open(f))
    n = j["config"]["requests"]
    misses = (1.0 - j["hit_rate"]) * n
    reads = max(misses - j.get("coalesced_loads", 0), 0.0)
    demand = j["throughput_rps"] * (reads / n) * KEY_BYTES[method] / 1e9
    return (round(demand / 1024 ** 3 * 1e9, 3), round(demand, 3),
            round(j["device_inflight_reads_mean"], 2))


def ceiling_fails(outdir, method):
    """Populations that do not hold the budget, from the stage-16 cells."""
    out = []
    for f in glob.glob(os.path.join(outdir, f"serving_{method}_tencap100_u*.json")):
        u = int(re.search(r"_u(\d+)\.json$", f).group(1))
        if json.load(open(f))["slo_violation_rate"] > 0.001:
            out.append(u)
    return sorted(out)


def ceiling_holds(outdir, method):
    """Populations whose 300 ms budget holds, from the stage-16 bisection."""
    holds = []
    for f in glob.glob(os.path.join(outdir, f"serving_{method}_tencap100_u*.json")):
        u = int(re.search(r"_u(\d+)\.json$", f).group(1))
        if json.load(open(f))["slo_violation_rate"] <= 0.001:
            holds.append(u)
    return sorted(holds)


def fit(ms, ts):
    """Least squares of worker time 8/T = W + S*m, in ms."""
    xs = list(ms)
    ys = [8000.0 / t for t in ts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    S = sxy / sxx
    W = my - S * mx
    ss = sum((y - my) ** 2 for y in ys)
    sr = sum((y - (W + S * x)) ** 2 for x, y in zip(xs, ys))
    return W, S, (1 - sr / ss if ss else 0.0)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(here, "out"))
    ap.add_argument("--dest", default=os.path.join(here, "out", "satcurve_fit.json"))
    args = ap.parse_args()

    band_path = os.path.join(args.out, "pool_readband.json")
    band = json.load(open(band_path)) if os.path.exists(band_path) else {}
    report = {"schema_version": 1,
              "budget_ms": BUDGET_MS,
              "plateau_from_tenants": PLATEAU_FROM,
              "ceiling_from_tenants": CEILING_FROM,
              "read_rate_source": os.path.basename(band_path),
              "schemes": {}}
    raw = {}
    for method in ("tfhe", "wwl+24", "ours"):
        c = curve(args.out, method)
        if not c:
            print(f"  no cells for {method}", file=sys.stderr)
            continue
        raw[LABEL[method]] = (method, c)
    infl = {}
    for name, (method, c) in raw.items():
        us = sorted(c)
        read_gbs = band.get("families", {}).get(FAMILY[method], {}).get("median_gbps", 1.2)
        qd1 = KEY_BYTES[method] / 1e6 / read_gbs
        W, S, _ = fit([c[u][1] for u in us], [c[u][0] for u in us])
        infl[name] = S / qd1
        raw[name] = (method, c, us, read_gbs, qd1, W, S)
    worst = max(infl.values())
    _, c_g, us_g = raw["GINX"][0], raw["GINX"][1], raw["GINX"][2]
    demand_g = {u: c_g[u][0] * c_g[u][3] * KEY_BYTES["tfhe"] / 1e9 for u in us_g}
    pinned = [demand_g[u] for u in us_g if u >= CEILING_FROM] or list(demand_g.values())
    CEILING_GBS = max(pinned)
    print(f"  device ceiling from GINX's demand asymptote, {CEILING_FROM}+ tenants: "
          f"{CEILING_GBS:.2f} GB/s = {CEILING_GBS / 1024 ** 3 * 1e9:.2f} GiB/s")
    for name, (method, c, us, read_gbs, qd1, W, S) in raw.items():
        key_mb = KEY_BYTES[method] / 1e6
        r2 = fit([c[u][1] for u in us], [c[u][0] for u in us])[2]
        demand = {u: c[u][0] * c[u][3] * KEY_BYTES[method] / 1e9 for u in us}
        hi = [demand[u] for u in us if u >= PLATEAU_FROM] or list(demand.values())
        holds = ceiling_holds(args.out, method)
        fails = [u for u in ceiling_fails(args.out, method) if u > (holds[-1] if holds else 0)]
        ceil_u = fails[0] if fails else (holds[-1] if holds else None)
        near = min(us, key=lambda u: abs(u - ceil_u)) if ceil_u else None
        to_gibs = lambda x: x / 1024 ** 3 * 1e9 # GB/s -> GiB/s
        worst_miss_ms = worst * qd1
        report["schemes"][name] = {
            "key_mb": round(key_mb, 2),
            "read_rate_gbs": round(read_gbs, 3),
            "qd1_read_ms": round(qd1, 1),
            "cells": len(us),
            "tenants_first": us[0], "tenants_last": us[-1],
            "quiet_ms": round(W, 1),
            "ms_per_miss": round(S, 1),
            "fit_r2": round(r2, 3),
            "miss_cost_over_read": round(S / qd1, 2),
            "throughput_first": round(c[us[0]][0], 2),
            "throughput_last": round(c[us[-1]][0], 2),
            "loss_first_to_last_pct": round((c[us[-1]][0] / c[us[0]][0] - 1) * 100, 1),
            "worker_time_first_ms": round(8000.0 / c[us[0]][0], 1),
            "worker_time_last_ms": round(8000.0 / c[us[-1]][0], 1),
            "worker_time_growth_ms": round(8000.0 / c[us[-1]][0] - 8000.0 / c[us[0]][0], 1),
            "miss_cost_at_last_ms": round(S * c[us[-1]][1], 1),
            "ms_per_miss_upper_ms": round(
                (8000.0 / c[us[-1]][0] - 8000.0 / c[us[0]][0]) / max(c[us[-1]][1], 1e-9), 1),
            "miss_cost_upper_over_read": round(
                (8000.0 / c[us[-1]][0] - 8000.0 / c[us[0]][0])
                / max(c[us[-1]][1], 1e-9) / qd1, 2),
            "demand_last_gbs": round(demand[us[-1]], 2),
            "coalesced_share_last": round(c[us[-1]][4], 3),

            "reads_per_request_last": round(c[us[-1]][3], 3),
            "demand_last_gibs": round(to_gibs(demand[us[-1]]), 2),
            "demand_200_gibs": (round(to_gibs(demand[200]), 2) if 200 in demand
                                 else round(to_gibs(demand[min(us, key=lambda u: abs(u-200))]), 2)),
            "demand_plateau_min_gbs": round(min(hi), 2),
            "demand_plateau_max_gbs": round(max(hi), 2),
            "demand_plateau_min_gibs": round(to_gibs(min(hi)), 2),
            "demand_plateau_max_gibs": round(to_gibs(max(hi)), 2),
            "reads_in_flight_last": round(c[us[-1]][2], 2),
            "miss_rate_at_last": round(c[us[-1]][1], 3),
            "hit_rate_at_last": round(1.0 - c[us[-1]][1], 3),
            "plateau_demand_gbs": round(max(hi), 2),
            "byte_bound": bool(len(us) >= 3
                               and max(demand[us[-1]], demand[us[-2]], demand[us[-3]])
                               / min(demand[us[-1]], demand[us[-2]], demand[us[-3]]) < 1.05
                               and demand[us[-1]] > 0.70 * CEILING_GBS),
            "crossing_miss_rate_at_30rps": (
                round(demand[us[-1]] * 1e9 / (30.0 * KEY_BYTES[method]), 3)
                if (len(us) >= 3
                    and max(demand[us[-1]], demand[us[-2]], demand[us[-3]])
                    / min(demand[us[-1]], demand[us[-2]], demand[us[-3]]) < 1.05
                    and demand[us[-1]] > 0.70 * CEILING_GBS) else None),
            "ceiling_tenants": ceil_u,
            "fail_demand_gibs": (ceiling_cell(args.out, method, ceil_u) or (None,))[0],
            "fail_demand_gbs": (ceiling_cell(args.out, method, ceil_u) or (None, None))[1],
            "fail_reads_in_flight": (ceiling_cell(args.out, method, ceil_u) or (None, None, None))[2],
            "floor_bytes_req_s": round(CEILING_GBS * 1e9 / KEY_BYTES[method], 1),
            "worst_miss_cost_ms": round(worst_miss_ms, 1),
            "floor_workers_req_s": round(8.0 / ((W + S) / 1000.0), 1),
            "floor_pessimistic_req_s": round(8.0 / ((W + worst_miss_ms) / 1000.0), 1),
        }
    ceiling_gbs = CEILING_GBS
    report["worst_miss_inflation"] = round(worst, 2)
    report["device_ceiling_basis"] = ("asymptote of GINX's miss-stream demand, "
                                      f"cells from {CEILING_FROM} tenants on, "
                                      "counting only reads that were issued")
    ss_path = os.path.join(args.out, SINGLE_STREAM_FILE)
    if ceiling_gbs and os.path.exists(ss_path):
        us = dig(json.load(open(ss_path)), SINGLE_STREAM_PATH)
        if us:
            single_gbs = KEY_GINX_BYTES / (us * 1e-6) / 1e9
            report["single_stream_gbs"] = round(single_gbs, 3)
            report["single_stream_gibs"] = round(single_gbs / 1024 ** 3 * 1e9, 3)
            report["single_stream_source"] = f"{SINGLE_STREAM_FILE}:{SINGLE_STREAM_PATH}"
            report["device_ceiling_gbs"] = round(ceiling_gbs, 3)
            report["device_ceiling_gibs"] = round(ceiling_gbs / 1024 ** 3 * 1e9, 3)
            report["bandwidth_range"] = round(ceiling_gbs / single_gbs, 2)
    os.makedirs(os.path.dirname(args.dest), exist_ok=True)
    with open(args.dest, "w") as fh:
        json.dump(report, fh, indent=1, sort_keys=True)
        fh.write("\n")
    for name, s in report["schemes"].items():
        print(f"  {name:<7} W {s['quiet_ms']:6.1f} ms + {s['ms_per_miss']:6.1f} ms/miss "
              f"= {s['miss_cost_over_read']:.2f}x its {s['qd1_read_ms']:.0f} ms read (R2 {s['fit_r2']:.3f})")
        print(f"          {s['throughput_first']:.2f} -> {s['throughput_last']:.2f} req/s over "
              f"{s['tenants_first']}-{s['tenants_last']} tenants ({s['loss_first_to_last_pct']:+.1f}%); "
              f"{s['worker_time_last_ms']:.1f} ms per request there, {s['miss_cost_at_last_ms']:.1f} of it misses")
        print(f"          demand {s['demand_plateau_min_gbs']:.2f}-{s['demand_plateau_max_gbs']:.2f} GB/s "
              f"({s['demand_plateau_min_gibs']:.2f}-{s['demand_plateau_max_gibs']:.2f} GiB/s), "
              f"{s['reads_in_flight_last']:.2f} reads in flight; breaks its budget at {s['ceiling_tenants']} "
              f"tenants, where the device carries {s['fail_demand_gibs']} GiB/s with "
              f"{s['fail_reads_in_flight']} reads in flight")
        print(f"          floor: bytes {s['floor_bytes_req_s']:.1f}, workers {s['floor_workers_req_s']:.1f}, "
              f"pessimistic {s['floor_pessimistic_req_s']:.1f} req/s")
    print(f"wrote {args.dest}")


if __name__ == "__main__":
    main()
