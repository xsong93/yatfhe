#!/bin/bash
# The paper's benchmark suite.
#
#   stage  what it produces                                              artifact
#   1      per-tenant serving keys, all three families to 800 tenants      server/BSK_*.bin
#   2      pressure sweep: 5 seeds x 5 caps x 2000 req                     *_results_<cap>_s<seed>.json
#   3      component ablation, part 1: generate the five key files         server/COMP_*.bin
#   4      thread scaling and core placement                               *_results_100_thr_*.json
#   5      SLO capacity, Poisson, 8000 req per point                       serving_*_slo8k_r*.json
#   6      capacity-cell repeats at the crossing, 60 and 62 req/s          serving_ours_slo8k_ref_r*_rep*.json
#   7      SLO capacity under burst arrivals                               serving_*_slo_burst_r*.json
#   8      decision-tree depth sweep, B = 1/2/4/8                          serving_*_tree_b*.json
#   9      concurrent crossover                                            serving_*_xo_b*.json
#   10     storage tiers: NVMe matched cell, then SATA and N2              *_results_10_{nvme,sata,n2}.json
#   11     sharded two-pool placement, one serving pool per CCD            serving_*_shard*.json
#   12     scheduler comparison: soft vs shared vs affinity                serving_ours_sched_*.json
#   13     paired micro-benchmarks: first step, NTT split, dispatch width  paired_fk/*, paired/ns*, paired/sw*
#   14     queue-depth points: in-flight key loads and device reads        serving_*_qd_*.json
#   15     shared-queue capacity controls                                  serving_*_shared_r*.json
#   16     component ablation, part 2: 600 warm / 600 cold on aged keys    comp_results.json
#   17     tenant ceiling at a common service rate of 30 req/s             serving_*_tencap100_u*.json
#   18     saturated throughput against tenant population                  serving_*_sat_u*.json
#   19     key-pool cold-read band, page cache dropped per file            out/pool_readband.json
#   20     saturation-curve fit behind the ceiling and saturation numbers  out/satcurve_fit.json
#
# No stage reads a result file from an earlier protocol.

set -u
B="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$B/out}"
say(){ echo "===== $1  $(date -Is) ====="; }

# Deploy fresh binaries (build with TORUS=32 ./build_bench_blindrotate.sh first)
missing=""
for f in blindrotate_cache blindrotate_comp blindrotate_serving \
         gen_benchkeys ntt_split_ablate stagea_width_ablate \
         stage_ratio_final; do
    if [ -f "$B/../CMAKE_BUILD/ya_benchmark/$f" ]; then
        cp "$B/../CMAKE_BUILD/ya_benchmark/$f" "$B/server/"
    else
        echo "  MISSING BINARY $f in CMAKE_BUILD/ya_benchmark/" >&2
        missing="$missing $f"
    fi
done
if [ -n "$missing" ]; then
    echo "aborting: the suite cannot run without:$missing" >&2
    echo "build them with TORUS=32 ./build_bench_blindrotate.sh" >&2
    exit 1
fi
# Prerequisite: ./build_bench_blindrotate.sh has been run and ./cp_bench.sh has populated $B/server/.

say "1/20 generate per-tenant serving keys (all three families to 800)"
( cd "$B/server" && ./gen_benchkeys --users=800 ) || exit 1

say "2/20 pressure sweep (5 seeds x 5 caps x 2000 req, pre-touch + settle warmup)"
CAPS="50,20,10,5,1" REQS=2000 SEEDS="1 2 3 4 5" "$B/run_sweep.sh"

# This stage only GENERATES the ablation's five key files.
# The 600 warm / 600 cold measurement is stage 15, after the rest of the suite has run.
# So keys are old enough.
say "3/20 component ablation, part 1: generate the five key files"
( cd "$B/server" && taskset -c 0-7 ./blindrotate_comp --warmreps=1 --coldreps=1 --out="$(mktemp -d)" >/dev/null 2>&1 )

say "4/20 thread scaling + CCD placement"
"$B/run_threads.sh"

say "5/20 SLO-capacity B=1  (8000 req/point, bracketed Poisson grids)"
SKIP_GRID=0 TAGBASE=slo8k "$B/run_slo.sh"

say "6/20 capacity-cell repeats at the crossing (60 and 62 req/s, three runs each)"
SKIP_GRID=1 REPEAT_SPECS="ours:60:3 ours:62:3" TAGBASE=slo8k_ref "$B/run_slo.sh"

say "7/20 burst arrival SLO-capacity (4x/200ms-on, 0.25x/800ms-off)"
"$B/run_slo_burst.sh"

say "8/20 decision-tree depth sweep (B = 1/2/4/8)"
"$B/run_trees.sh"

say "9/20 concurrent crossover"
"$B/run_crossover.sh"

# Update to your system before running.
say "10/20 storage tiers (NVMe matched cell, SATA then N2)"
"$B/run_nvme.sh"
"$B/run_sata.sh"
"$B/run_n2.sh"

# Update to your system before running.
say "11/20 sharded two-pool placement (one serving pool per CCD)"
"$B/run_shard.sh"

# Different schedulers test.
say "12/20 scheduler comparison (soft vs shared vs affinity)"
"$B/run_sched.sh"

# Every number this one reports are under $OUTDIR/paired/.
#   paired_fk/run*/firststep_fixedkey_results.json for the first-step difference,
#   derive_timing.json for the derivation it trades against,
#   sw*/stagea_width_results.json for the dispatch width,
#   ns*/ntt_split_results.json for the NTT split,
#   ratio.txt for T_SS'/T_EP.
say "13/20 paired micro-benchmarks (firststep, NTT split x2, dispatch width, stage ratio)"
OUTDIR="$OUTDIR" "$B/run_paired.sh"

say "14/20 queue-depth points (in-flight key loads and device reads, Section IV-A)"
"$B/run_qd.sh"

say "15/20 shared-queue capacity controls (Table VI shared row, ~6.5 h)"
"$B/run_sched_capacity.sh"

# The ablation's cold column, on aged keys.
say "16/20 component ablation, part 2: 600 warm / 600 cold on aged keys"
"$B/agedkeys.sh"
"$B/fold_device.sh" 32 "$B/server/agedkeys/COMP_A_GINX.bin" 1.6
( cd "$B/server" && taskset -c 0-7 ./blindrotate_comp --keyprefix=agedkeys --out="$OUTDIR" --warmreps=600 --coldreps=600 )

say "17/20 tenant ceiling at a common service rate (30 req/s, two memory budgets)"
"$B/run_tencap.sh"

say "18/20 saturated throughput against tenant population (30 req/s ceiling marked)"
"$B/run_satcurve.sh"

# Reads all 1800 files and fails the suite if any of them is still served from the fast cache.
say "19/20 key-pool cold-read band (1800 files, page cache dropped per file)"
python3 "$B/probe_pool.py" || exit 1

say "20/20 saturation-curve fit behind the ceiling and saturation numbers"
"$B/satcurve_model.py"

echo "===== COMPLETE $(date -Is) ====="
