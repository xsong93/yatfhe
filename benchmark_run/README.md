# Benchmark workflow and result files

Paths below are relative to this directory.

## Running this on other systems

The key pool is 234 GB on disk, so key generation needs that much free space.

Every measurement stage pins cores 0-7 and expects 8-15 to stay idle,
and the sharded placement stage runs one serving pool per CCD.
You may want to play around the params on different machines prior to running the full suit.

The storage stages name the device (nvme0n1) and depend on the drive's own behavior.

The whole measurement path is GNU/Linux on x86-64.
The cold path drops the file from the page cache with posix_fadvise(POSIX_FADV_DONTNEED).
The device counters come from the kernel's view of the named block device,
fold_device.sh uses dd with oflag=direct,
and cores are pinned with taskset.
Windows has none of these, and macOS offers neither POSIX_FADV_DONTNEED nor taskset.
Intel HEXL is x86-only as well, which puts ARM hosts out of luck.

The binaries are built with -march=native and link a HEXL build for one ISA.
Rebuild with build_bench_blindrotate.sh on the target machine, choosing TARGET_ARCH=avx2 or generic when the CPU differs.

## Result files and protocol labels

Every open-loop serving run writes `out/serving_<method>_<tag>.json`,
which records its own protocol in `config` (`requests`, `warm`, `arrival`, `batch`, `sched`, `users`, `seed`).

The tags

| tag | meaning                                                                                                       |
| --- |---------------------------------------------------------------------------------------------------------------|
| `slo8k_r<rate>` | capacity series, 8000 requests per point, warm 500, Poisson, soft affinity                                    |
| `shared_r<rate>` | the same series with `--sched=shared`                                                      |
| `slo_burst_r<rate>` | capacity series under burst arrivals                                                                          |
| `tree_b<B>` | depth-`d` decision-tree queries                                         |
| `qd_<method>_r<rate>`, `qd_<method>_sat` | queue-depth instrumentation runs (`--devstat`), one per point                                                 |
| `xo_b<B>` | saturated-throughput crossover, 8000 requests, warm 500, batch `B`                                            |
| `sched_{soft,affinity,shared}_r<rate>` | scheduler comparison, 800 requests, warm 300                                                                  |
| `shard{0,1}_*`, `shardp{0,1}_*` | two serving pools, one per core complex                                                                       |
| `refpool_r<rate>` | single-pool reference for the sharding comparison                                                             |
| `tencap100_u<N>` | tenant ceiling: 8000 requests at 30 req/s on 17 GB, cap 100 slots, `N` tenants                                |
| `sat_u<N>` | saturation curve: 1000 requests offered at 1000 req/s, cap at the scheme's iso-byte share (100/200/400), `N` tenants |

Storage-tier runs use a different naming scheme, with `sata` for the SATA SSD.
`n2` for the USB 3 SSD and `nvme` for the local NVMe.
Change those to match the disk lables on your system.

### The one-pass suite

`run_all.sh` regenerates every artifact below, and the stage numbers are the ones its banners print.

| stage | experiment | artifact |
| --- | --- | --- |
| 1 | per-tenant key pool: 800 tenants of each scheme | `server/BSK_*.bin` |
| 2 | pressure sweep, 5 reps x 5 caps x 2000 req | `out/{method}_benchmark_results_<cap>_s<seed>.json` |
| 3 | component ablation, part 1: generate the key files | `server/COMP_*.bin` |
| 4 | thread scaling and core placement | `out/*_100_thr_*.json` |
| 5 | SLO capacity, Poisson, 8000 req per point | `out/serving_*_slo8k_r*.json` |
| 6 | capacity-cell repeats at the crossing | `out/serving_ours_slo8k_ref_r*_rep*.json` |
| 7 | SLO capacity, bursty arrivals | `out/serving_*_slo_burst_r*.json` |
| 8 | decision-tree depth | `out/serving_*_tree_b*.json` |
| 9 | concurrent crossover | `out/serving_*_xo_b*.json` |
| 10 | storage tiers: NVMe, SATA, N2 | `out/*_benchmark_results_10_{nvme,sata,n2}.json` |
| 11 | sharded two-pool placement | `out/serving_*_shard*_*.json` |
| 12 | scheduler comparison | `out/serving_ours_sched_*.json` |
| 13 | paired micro-benchmarks | `out/paired*/**`, `out/paired/ratio.txt` |
| 14 | queue-depth points | `out/serving_*_qd_*.json` |
| 15 | shared-queue capacity controls | `out/serving_*_shared_r*.json` |
| 16 | component ablation, part 2: the 600/600 cold measurement on aged keys | `out/comp_results.json` |
| 17 | tenant ceiling at a common load, 30 req/s | `out/serving_*_tencap100_u*.json` |
| 18 | saturated throughput against tenant population | `out/serving_*_sat_u*.json` |
| 19 | cold-read band of the pool keys | `out/pool_readband.json` |
| 20 | the fit behind the ceiling and saturation numbers | `out/satcurve_fit.json` |

Stage 13 holds the first-step A/B on one fixed key (`paired_fk/run1`-`run8`),
the NTT-split sweep (`paired/ns1`, `paired/ns2`),
the dispatch-width A/B (`paired/sw1`-`paired/sw4`),
the stage ratio and the caller-side derivation.
Stage 20 needs `numpy`.

`out/` holds one protocol generation per included experiment.
Re-running the suite rewrites the pool.

### Key age

A cold column is a device read.
Some devices might serve data written seconds ago about twice as fast as the same *aged* data.

**A measurement must never generate the keys it reads.** 
Folding the device removes the SLC state but does not age a file,
so generation and measurement have to be separated in time as well.

Identical bytes read faster on freshly written blocks,
so the cold column follows the age of the *blocks*.
The system is responsible for ensuring data fairness.

How the suite implements it.

* The serving pool (`BSK_*.bin`, one file per tenant and scheme) is written by `gen_benchkeys` as a separate step (stage 1) and is never regenerated inside a run.
Regenerating it makes every cold read in every serving stage optimistic.

* The component ablation needs five key types,
three of which also exist in the pool.
Stage 2 only *generates* them.
The measurement is stage 15,
after the rest of the suite has run,
so those two remaining files have aged for the length of the run.
`agedkeys.sh` points the five ablation names at the most aged file of each type through **symlinks**, and never copies.

* A key file rewritten after the pool was built is a fresh file.
The fix is simple:
regenerate deliberately,
write 64 GB through the device,
let it settle, then *check*.

* Stage 16c is that check for the whole pool.
`probe_pool.py` reads all 2000 key files cold with the page cache dropped per file and drops it again afterward.
It fails the suite if any file is still served from the drive's fast cache.

* Across sessions, a cold column is a snapshot of one session and must never be compared with a cold number measured in another, even on the same file.
