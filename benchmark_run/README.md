# benchmark_run

## Layout

```
server/          binaries and tenant key files (BSK_*.bin)
out/             every benchmark JSON lands here
fig/             every generated figure lands here
run_*.sh         one experiment each; run_all.sh drives them in order
aggregate.py     pools seeds into the main results table
make_figures.py  out/ -> fig/, with --sync to push into the paper
make_tables.py   out/ -> fig/*.tex, with --sync to push into the paper
```

All three benchmarks take `--out=DIR` (or `$YA_OUT_DIR`, default the working
directory). The run scripts pass `--out=$OUTDIR`, which defaults to `out/`;
override with `OUTDIR=/some/path ./run_sweep.sh`. Results no longer land beside
the binaries.

Binaries are built in the parent project and copied into `server/`:
`blindrotate_cache` (closed-loop, cache pressure), `blindrotate_serving` (open-loop,
Poisson arrivals), `blindrotate_comp` (7-rung ablation).

## Experiments

| Script | What it measures | Feeds |
|---|---|---|
| `run_sweep.sh` | Latency vs cache pressure, 5 seeds × 6 capacities × 1000 req | Table II, Fig. 5, Fig. 6 |
| `blindrotate_comp` | 7-configuration ablation, warm/cold/interleaved | Table III |
| `run_sata.sh`, `run_n2.sh` | Storage tiers — SATA and USB3, vs local NVMe | Table IV (top) |
| `run_slo.sh` | SLO-capacity: offered rate vs P99.9, **8000 req/point**, per-method bracketed grids | Table IV (bottom) |
| `run_crossover.sh` | Saturated throughput vs batch size — locates the B≥2 crossover | §IV-B3 |
| `run_threads.sh` | Thread scaling and CCD placement (pure compute, HR=1.00) | §IV-B5 |
| `run_tenants.sh` | Tenant count: sweep A (cap fixed) and B (pressure pinned) | §IV-B6 |

## Running

Full campaign (~several hours, stages are serial):

```bash
./run_all.sh 2>&1 | tee campaign.log
```

Individual stages take env overrides:

```bash
./run_sweep.sh                              # defaults now match the paper
./run_slo.sh                               # defaults: 8000 req/point
RATES_ours="30 40 50" RATES_wwl24="" RATES_tfhe="" ./run_slo.sh   # one method
( cd server && taskset -c 0-7 ./blindrotate_comp --warmreps=40 --coldreps=50 )
```

**Pinning is mandatory.** All scripts pin to `taskset -c 0-7`, the eight physical cores of
the 9950X3D's V-Cache complex. An unpinned run migrates between complexes and reports the
resulting scheduling noise as latency jitter — the metric under study. Cross-complex
placement costs 1.8×, SMT 1.2×.

**Storage runs need the device present.** `run_sata.sh` waits for a writable `/mnt/sata`,
then copies 50 tenants' keys to it; `run_n2.sh` expects the USB3 disk at `/media/spark/N2`.
Both stage keys, run, copy the JSON back to `server/`, and clean up.

Keys are not regenerated between runs. Regenerate only if the key structs change — the
gadget-decomposition optimisation was bit-identical, so existing keys stayed valid.

## Diagnostics

Not part of `run_all.sh` — one-off tools, built by the CMake glob into
`CMAKE_BUILD/ya_benchmark/`, run by hand. Two of them back claims in the paper,
so keep them buildable.

| Tool | What it establishes | Used by |
|---|---|---|
| `gd_equiv` | The digit-major gadget decomposition is bit-identical to the original, 200 random trials at each level L=1…8 | §IV-A "verified bit-identical" |
| `gd_compare` | Per-level cost of the two decomposition forms; locates the register-spill cliff above L=4 | §IV-A, the 25% / 6% split |
| `ntt_split` | Splits `applyNtt` into conversion vs transform. Concluded ~3% available — the NTT relocation was **reverted** | nothing in the paper |
| `pipe_stages` | SS'/EP steady-state stage ratio (NS' 5.0 / SS' ~30 / EP ~17 us) | §III-B4, the 1.7 figure |
| `prim_timing` | NTT / pointwise / decomposition primitive costs | PART 0d analysis |
| `stage_timing` | Superseded by `pipe_stages` — measures WWL+24's SS, not our SS'. Do not use | nothing |

The two non-AVX-512 build directories referenced above are made with
`yatfhe/build_avx2.sh` and `yatfhe/build_generic.sh`, which are copies of
`build_bench_blindrotate.sh` with `BUILD_DIR` changed so the campaign binaries in
`CMAKE_BUILD` are not overwritten. They link the matching HEXL from
`~/.hexl/{avx2,generic}`.

## Producing results

Main table, pooling all seeds with bootstrap CIs on the tail percentiles:

```bash
venv/bin/python aggregate.py out/ --format latex --baseline OURS
```

Figures:

```bash
venv/bin/python make_figures.py            # out/ -> fig/
venv/bin/python make_figures.py --sync     # and copy into p1_tc/src/pics/
```

Tables:

```bash
venv/bin/python make_tables.py             # out/ -> fig/tab_*.tex
venv/bin/python make_tables.py --sync      # and copy into p1_tc/src/tables/
```

Three tables are generated -- `tab_performance` (Table II), `tab_ablation`
(Table III) and `tab_storage` (Table IV). Each is a complete float; `exp.tex`
pulls them in with `\input{tables/tab_*}`, so **editing them by hand in the paper
is pointless -- the next sync overwrites it**. Change the generator instead.

The parameter table and the blind-rotation variants table are *not* generated:
neither comes from measurement (one is configuration, the other is arithmetic over
published key structures), so both stay inline in `exp.tex`.

`fig/` is the generation target; the paper compiles against `src/pics/` and
`src/tables/`. Use `--sync` so they cannot drift apart unnoticed. The acceptance
test for both generators is that a rebuild leaves the rendered PDF byte-identical.

Ablation numbers come from `out/comp_results.json` (key `case`, one object per rung with
`key_mb`, `warm`, `cold`, and `interleaved` on rungs F and G). Use **medians** throughout —
mixing medians and means across rungs produces attribution percentages that do not
reconcile.

## Gotchas

**Result filenames encode cache capacity, not pressure ratio.**
`ours_benchmark_results_10_sata.json` is capacity 10, and its `config.pressure_ratio` is
5.0. Read the config field. This collision put a table of nine wrong numbers into a draft.

**`_s[1-5]` globs collide with tenant runs.** `*_s[1-5].json` also matches
`..._tenA_u100_s1.json`. Match `_benchmark_results_\d+_s[1-5]\.json` instead.

**Do not try to synthesise storage tiers with cgroup `blk-throttle`.** Under a 500 MB/s
read cap, OURS reported implied bandwidth above 1 GB/s — impossible. `blk-throttle`
accounts block-layer submissions and our interleaved reads are charged elsewhere, so it
measures the instrument rather than the system. Use real devices. (`run_bandwidth.sh`,
which did this, is in `archive/old_scripts/`.)

**Never time from `cmake-build-debug/`.** It is `CMAKE_BUILD_TYPE=Debug` (`-O0`)
and it looks like an ordinary build directory. HEXL is a separately installed,
already-optimised library, so `-O0` slows only *our* code — which hits each stage
in proportion to how little HEXL it contains, and inverts the stage ordering:

| build | type | NS' | SS' | EP | SS'/EP |
|---|---|---|---|---|---|
| `CMAKE_BUILD` | Release, AVX-512 | 4.9 | 29.3 | 17.1 | **1.73** |
| `CMAKE_BUILD_AVX2` | Release, AVX2 | 6.9 | 134.0 | 55.2 | 2.41 |
| `CMAKE_BUILD_GENERIC` | Release, scalar | 9.2 | 145.4 | 62.5 | 2.34 |
| `cmake-build-debug` | **Debug** | 102 | 80 | 164 | **0.49** |

NS' is pure hand-written loops and slows 21x; EP is our glue around HEXL calls and
slows 9.6x; SS' is dominated by HEXL NTTs and slows only 2.7x. The debug run
therefore prints `ratio SS'/EP = 0.49 -> CLAIM HOLDS`, confidently and wrongly.
The release ordering is the real one, and it is stable: removing vectorisation
widens the ratio rather than inverting it. Always run from `CMAKE_BUILD`, pinned.

**Ignore `pipe_stages`' "PROPOSED REBALANCE".** It predicts a 45% gain from moving
the L NTTs out of SS' into NS', computed as `max(stageA, stageB, EP)`. That model
assumes the stages run one at a time on dedicated cores. They do not — both dispatch
into a shared pool, so nine units contend for eight cores. The relocation was
implemented and measured at **-3.5% (a regression)** and reverted. The printed
suggestion is stale; the ratio and the per-stage times above are the useful output.

**P99.9 needs 8000 requests per point, and seeds cannot substitute.**
`blindrotate_serving` stores only summary percentiles -- no per-request records -- so
runs cannot be pooled afterwards. At 1000 requests the reported P99.9 *is* the single
worst observation, and WWL+24's capacity moved between 10 and 30 req/s across
repetitions on that basis. Only `--requests` fixes it. `run_slo.sh` defaults to 8000.
(`blindrotate_cache` *does* store per-request `time_us` and can be pooled -- do not
confuse the two.)

**Cold runs need repeats.** The first cold repetition runs ~4 µs slow. `--coldreps` below
~20 gives a standard error larger than the effects being measured.

**Archived result sets are from older builds.** See `archive/README.md` — `n640` (wrong LWE
dimension), `embedded_lut` (pre-LUT-independent key), `pre_gdopt` (before the shared
gadget-decomposition fix, which sped the baseline up more than ours). Never pool these with
`server/`. `archive/final_2026-08-26/` is a snapshot of the campaign the paper reports.
