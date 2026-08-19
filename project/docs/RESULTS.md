# RESULTS.md — the measured findings

This is the tracked home for the numbers. `results/` is gitignored, so
`results/bench_summary.md` (which `make bench-plot` generates from the CSVs)
does not survive a fresh clone. This file does, and it carries the provenance
alongside each number so any of them can be challenged.

**Every timing here came from a rangpur compute node via `sbatch`.** The
provenance stamp is written into each CSV by `write_header()` in
`slurm/bench_common.sh` and reads `host=r730-2.compute.eait.uq.edu.au`, which
is a compute node and not the login node. `scripts/bench_plot.py` and
`scripts/roofline.py` refuse to render a figure from any CSV whose `host=`
header is not a `*.compute.eait.uq.edu.au` name.

Benchmark node: Intel Xeon E5-2670 v3 (Haswell, 2014), 32 KiB L1d, 256 KiB L2,
30 MiB L3 (verbatim from `lscpu` in every CSV header).
Compiler: `g++ (GCC) 8.5.0 20210514 (Red Hat 8.5.0-28)`.
Flags: `-std=c++17 -O2 -Wall -Wextra -Iinclude`.
Jobs: 556505 (A), 556506 (B, superseded), 556507 (D), 556508 (C),
556513 (peak), 556515 (alloc), 556517 (E) on 2026-08-09; 556833 (F),
556834 (B, pinned re-run) and 556839 (G) on 2026-08-11.

---

## 0. Nomenclature, units and measurement discipline

Every number in this file carries a unit, and the units are defined once,
here. Conventions first: memory is binary (KiB, MiB — a huge page is 2 MiB,
which is not 2 MB); money is US dollars, written `$`; throughput is
`Mcell-updates/s` as defined below.

**The instrument.** Everything priced here is one contract, from
`config/reference.cfg`: a European call on the S&P 500 index, strike
K = $5250, spot S = $5200, maturity T = 0.25 years, risk-free rate
r = 4.5 %/yr, dividend yield q = 1.3 %/yr, under Heston dynamics with
v0 = θ = 0.04 (a 20 %/yr starting volatility). A "price" is the value of
that option in dollars.

**The grid.** `ns×nv` means ns stock-price nodes by nv variance nodes; `nt`
is the number of timesteps. The solver keeps two solution buffers of
`ns·nv` doubles each, so the **working set** is `2 · ns · nv · 8` bytes —
16 MiB at the reference footprint of 2048×512. A **cell update** is one
application of the nine-point finite-difference stencil to one (S, v) node
for one timestep; a full solve performs exactly `ns · nv · nt` of them.

**Throughput** is reported in **Mcell-updates/s**: millions of cell updates
per second of wall-clock time, `ns · nv · nt / seconds / 10⁶`. This is the
CSV column `cell_updates_per_sec` that the harness computes from a
`std::chrono::steady_clock` timer wrapped around the stepping loop only
(allocation, payoff setup and readout excluded; snapshot I/O disabled in
every benchmark job). It is a *rate*, which is why a 400-step benchmark run
and a 980 000-step production run measure the same quantity. A **speedup
(×)** is always a ratio of two throughputs, higher meaning faster, taken on
the same grid, the same node and the same page regime — the regime is named
wherever a speedup is quoted, because §5 shows it is worth 1.66× on its own.

**Uncertainty.** Every throughput is the **median of 5 repetitions** run
inside one process in one slurm job; the quoted `min–max` is the observed
range of those 5, which is a spread, not a standard error. Each sweep is a
single job on a single node (n = 1 at job level), so job-to-job variability
is not characterised — §10 says what that does and does not permit. The
Monte Carlo `±` in §1 and §7 is different: it is the **1σ standard error of
the mean** over 400 000 paths.

**Errors.** A "relative gap" between a and b is `|a − b| / |b|`. Parity
residuals are absolute amounts in $. "ulp" means units in the last place of
an IEEE double. The reported vega is ∂V/∂v0 per unit of **variance**, not
of volatility — the market convention divides by 2√v0, and the code says so
where it computes it (`src/solver.cpp`).

**The sweeps.** Letters name benchmark jobs; operational detail and
reproduction commands live in `RUNNING.md` §3.

| sweep | what it measures | script |
|---|---|---|
| A | the 9-rung optimisation ladder, allocator left alone (mixed page regime) | `slurm/bench_ladder.sh` |
| B | scaling with problem size, baseline vs level 6 (re-run pinned, job 556834) | `slurm/bench_scaling.sh` |
| C | `-fno-tree-vectorize` control — vacuous as run, see §6 | `slurm/bench_novec.sh` |
| D | negative controls at their own nt, paired in-job | `slurm/bench_controls.sh` |
| E | the ladder with huge pages pinned — the headline regime | `slurm/bench_pages.sh` |
| F | compiler-flag control: `-O3` and `-ffast-math` legs (job 556833) | `slurm/bench_flags.sh` |
| G | the valid vectorisation control: `-O3` vs `-O3 -fno-tree-vectorize` (job 556839) | `slurm/bench_novec_o3.sh` |

---

## 1. Headline

**Table 1** — the findings in one screen. Grids 2048×512 unless stated;
speedups are throughput ratios (§0).

| claim | value | source |
|---|---|---|
| Reference price | $196.1683699 (discretisation-limited; §7 puts the trustworthy digits at ~$196.17) | `config/reference.cfg` |
| vs semi-analytic Fourier | $196.1692, a relative gap of 4.2e-6 | `make mc-check` |
| vs Monte Carlo, 400k paths | $196.1345 ± $0.5524 (1σ standard error), brackets both | `make mc-check` |
| Serial speedup, baseline to level 6 | **2.56× throughput** huge-page regime / 1.81× default allocator | sweeps E / A |
| of which strength reduction alone | **+148.7 % throughput** | sweep E |
| accuracy cost of that speedup | 3.3e-15 relative on the price, about 15 ulp | `test_opt_matches` |
| levels bit-identical to the baseline | 0 and 1, max_abs_diff exactly $0.0 | `test_opt_matches` |
| page size alone, 2 MiB vs 4 KiB pages | **1.66× throughput** | sweep E vs A, `diag_alloc.sh` |
| traversal order, `ctl-order` control | 2.1× slower (throughput) | sweep E |
| fusing the v=0 row, `ctl-branch` control | free, 1.0× | sweep E |
| auto-vectoriser's contribution to the ladder | none — L6 measures 135.0 with the vectoriser off at `-O3`, vs 134.8 in the headline sweep (sweep C had claimed this but was vacuous, §6) | sweeps F, G |
| what the vectoriser does to the *baseline* | ×1.92 throughput (52.6 → 100.7 Mcell-updates/s at `-O3`), still 31 % short of hand-optimised L2; falls back to 52.6 with `-fno-tree-vectorize` | sweep G |
| level 6 vs measured scalar peak | **82 %** of 5.117 GFLOP/s at 2048×512, 95 % at 4096×1024 | `roofline.py` + `bench_peak.sh` |

---

## 2. The optimisation ladder (sweep E, huge-page regime)

A single production solve is a fresh process that allocates once, so it gets
huge pages. This is therefore the regime a real user experiences, and it is the
headline table.

**Table 2** — the ablation ladder. Sweep E, job 556517, grid 2048×512,
nt = 400, `config/bench.cfg`, huge-page regime, single job on r730-2.
Throughput in Mcell-updates/s (§0); min–max is the range of the 5 reps.

| level | technique added | throughput (Mcell-updates/s, median of 5) | min–max (Mcell-updates/s) | speedup vs baseline (× throughput) | step gain (% throughput) |
|---|---|---|---|---|---|
| baseline | reference solver | 52.6 | 52.5–52.6 | 1.00× | — |
| opt-L0 | none (anchor) | 52.6 | 52.6–52.6 | 1.00× | −0.1 % |
| opt-L1 | hoisting, lookup, CSE | 52.7 | 52.7–52.8 | 1.00× | +0.3 % |
| opt-L2 | **strength reduction** | **131.2** | 130.8–131.4 | **2.49×** | **+148.7 %** |
| opt-L3 | traversal order | 131.5 | 131.4–131.7 | 2.50× | +0.3 % |
| opt-L4 | loop splitting | 131.3 | 131.0–132.0 | 2.50× | −0.2 % |
| opt-L5 | induction variable | 134.2 | 133.5–134.3 | 2.55× | +2.2 % |
| opt-L6 | unrolling ×4 | 134.8 | 134.6–135.2 | 2.56× | +0.5 % |

A note on the price column of the bench CSVs: `config/bench.cfg` shrinks the
maturity to 0.00051 years (≈ 4.5 hours) so a rep fits the QOS window, which
makes the printed price (~$0.0003) financially meaningless — it exists so the
result is consumed and the compiler cannot delete the loop (PLAN §6). The
grid, the working set and the stencil traversal are byte-identical to the
reference footprint, which is what a throughput number depends on.

Sweep A (job 556505) is the same ladder with the allocator left alone, where
reps 2 to 5 are served from the reused glibc heap on 4 KiB pages: baseline
44.8, level 2 80.0, level 6 81.1 Mcell-updates/s, so 1.81× overall.

### The finding

**Only one rung out of six did anything.** GCC 8.5 at `-O2` had already
performed hoisting, common subexpression elimination, induction-variable
simplification and unrolling, so every rung except level 2 came back inside the
run-to-run noise.

Level 2 is the only transformation on the list the compiler is **forbidden** to
make, because replacing `x / (2*ds)` with `x * (1/(2*ds))` changes the answer in
the last bits and `-ffast-math` is off.

> At `-O2`, the hand optimisations worth doing are the ones the compiler is not
> allowed to do for you.

That is only visible because each rung was measured separately instead of
shipping one before-and-after pair.

### The claim, tested rather than asserted (sweep F, job 556833)

The paragraph above is a causal claim, so a job was built to attack it:
rebuild at `-O3` (licensing everything legal that `-O2` withholds) and at
`-O2 -ffast-math` (licensing exactly the reciprocal trick), each in the same
job as an `-O2` reference, huge-page regime throughout. What happened:

- **The `-O2` leg replicated the headline rungs independently** — baseline
  52.6, L2 131.5, L6 135.5 Mcell-updates/s in a fresh job two days after
  sweep E measured 52.6 / 131.2 / 134.8. The headline table reproduces
  across jobs to three digits.
- **At `-O3` the baseline nearly doubled** — 52.6 → 100.7 Mcell-updates/s
  (×1.92) — **while L2 (131.9) and L6 (134.4) did not move.** So the
  compiler, given its full legal licence, closes about half the gap to the
  hand-optimised kernel and stops 31 % short of it; and once the divisions
  are gone by hand, `-O3` has nothing left to add. Sweep G then attributed
  the whole lift to the vectoriser (§6, Table 6b).
- **The `-ffast-math` leg was refused by the iron rule, exactly as
  pre-worded in the script.** Under `-ffast-math`, level 1 — which must be
  bit-identical to the baseline — drifted by 9.1e-13, `test_opt_matches`
  failed, and `build_and_validate()` exited rather than benchmark. The
  licence that permits the reciprocal substitution changes answers
  elsewhere too, which is the claim's own premise arriving as a test
  failure. The sharp prediction (a fast-math baseline jumping to L2's
  throughput) therefore remains unmeasured on this machine, and the honest
  wording of the finding stays "the transformation `-O2` is forbidden to
  make", not "the transformation the compiler cannot ever make".

Levels 3 and 4 were predicted null **by construction**, because
`BaselineSolver` was already written with the outer-variance, inner-stock loop
order and already peeled the `v = 0` row into its own loop. The prediction was
confirmed, and what those rungs would have been worth was measured instead by
the negative controls.

### Negative controls (sweep E)

Deliberately worse code doing arithmetically identical work.
`test_opt_matches` holds both to the same tolerance as the ladder, so a
different answer would mean a bug rather than a technique.

**Table 3** — negative controls. Sweep E, job 556517, same grid and regime
as Table 2; each control is compared against opt-L5 measured in the same job.

| control | what it ablates | throughput (Mcell-updates/s, median of 5) | slowdown vs opt-L5 in the same job (× throughput) |
|---|---|---|---|
| opt-ctl-order | loops swapped | 63.6 | 2.1× |
| opt-ctl-branch | v=0 row fused back in | 135.3 | 1.0× |

Both start from level 5 and change exactly one structural thing, and both are
compared against level 5 measured **in the same job at the same nt**, because a
cross-job comparison would measure run length rather than technique.

---

## 3. Scaling with problem size (sweep B)

**Table 4** — scaling, pinned re-run (job 556834). Per-grid nt keeps every
row at ~4e8 cell updates (equal work); `MALLOC_MMAP_THRESHOLD_` pinned so
every rep of every row is served the same way; the 4096×1024 row reads out
deep in the money (`config/bench_itm.cfg` — identical grid, readout point
moved) so its price carries digits instead of underflowing. The page-regime
column is **measured**, not assumed: the two small grids' throughput matches
sweep A's 4 KiB numbers, the two large grids match the pinned huge-page
sweeps (E, F). Transparent huge pages need a 2 MiB-aligned 2 MiB stretch to
promote, and only the 8 MiB and 32 MiB buffers got that; the pin therefore
buys within-row consistency everywhere, and huge pages only where the
buffers are big enough.

| grid (ns×nv) | working set (MiB) | page regime (measured) | baseline (Mcell-updates/s) | opt-L6 (Mcell-updates/s) | speedup (× throughput) |
|---|---|---|---|---|---|
| 512×128 | 1 | 4 KiB | 44.4 | 155.4 | 3.50× |
| 1024×256 | 4 | 4 KiB | 41.8 | 104.4 | 2.49× |
| 2048×512 | 16 | 2 MiB (THP) | 52.6 | 135.1 | 2.57× |
| 4096×1024 | 64 | 2 MiB (THP) | 57.8 | 157.5 | 2.72× |

Two readings, one honest caveat. The **baseline** is nearly flat from 1 MiB
to 64 MiB of working set — a division-bound kernel does not care where its
operands live, which is §4's roofline finding arriving from another
direction. The **level-6 curve** does feel the memory system: it dips where
the working set first spills past L2 into shared L3 (104.4 at 4 MiB) and
recovers as huge pages arrive with the bigger buffers. The caveat: the rows
still straddle two page regimes — that is now a physical constraint (a
buffer smaller than one huge page cannot be THP-backed), not an allocator
accident, and the within-regime comparisons are the clean ones.

The first run of this sweep (job 556506, retired to `results/superseded/`)
did confound the two: nothing pinned the allocator, so reps 1 and 2–5
straddled regimes within single rows, its 2048×512 row landed in the 4 KiB
regime (44.8 / 81.1 — the sweep A numbers), and its largest row's
out-of-the-money price underflowed to ~$2e-16. Its curve read backwards for
a cache story and was published with a warning; the re-run above replaces
it.

---

## 4. Two claims that did not survive being measured

Both belong in the video. A corrected assumption is a better story than an
assumption that was never tested.

### "The solver is memory-bound." It is not.

`scripts/roofline.py` counts 31 flops and 24 bytes of compulsory traffic per
cell update straight off the level-6 kernel — an arithmetic intensity of
31/24 ≈ 1.29 flop/byte, which is the x-axis of `results/roofline.png` — and
`slurm/bench_peak.sh` measures this node's achievable **single-core scalar
double-precision** peak at **5.117 GFLOP/s** rather than quoting a datasheet.

- Level 6 achieves 4.18 GFLOP/s at the reference grid, **82 % of peak**, and
  4.86 GFLOP/s at 4096×1024, **95 %**, while demanding about 3 GB/s from a
  four-channel DDR4 machine.
- The baseline sits at **27 % of peak** because it is *division*-bound, at
  five `divsd` per cell update.

So optimising the arithmetic moved this kernel from division-bound to
compute-bound, and the next win has to come from doing more work per
instruction. That is exactly what SIMD is, which makes it the natural bridge
into Milestone 2.

### "`ctl-order` costs about 30×." It costs 2.1×.

The 30× was a local measurement taken while the reference solve was running
concurrently, and a bandwidth-hungry kernel is precisely what a co-tenant hurts
most. On an exclusively allocated node it is 2.1×.

The lesson, now encoded as `--exclusive` in every benchmark job, is that it is
easy to measure something real that answers a different question than the one
you asked.

---

## 5. The page-size effect

The same binary on the same node runs at 134.8 or 81.1 Mcell-updates/s
depending on **nothing but whether the two 8 MiB grid buffers are backed by
2 MiB transparent huge pages or 4 KiB pages.** The arithmetic: 16 MiB is
8 pages of 2 MiB, or 4096 pages of 4 KiB — eight TLB entries one way, four
thousand the other, against a 1024-entry second-level TLB on this core.

glibc's dynamic `mmap` threshold means rep 1 of a `--bench` run gets huge pages
and reps 2 to 5 do not, so a plain median-of-5 was silently averaging two
different memory regimes.

Established by elimination rather than by guessing. `bench_peak.sh` sampled the
clock during an actual level-6 solve and found it pinned at 3100 MHz
throughout, ruling out throttling. `diag_alloc.sh` then pinned
`MALLOC_MMAP_THRESHOLD_` by hand and made the effect appear and disappear on
demand.

**Table 5** — the page-size effect in isolation. Same binary, same grid
(2048×512), same node; 4 KiB column from sweep A (job 556505), huge-page
column from sweep E (job 556517) with `MALLOC_MMAP_THRESHOLD_` pinned to
1 MiB. Cross-job, but same-day, same node, and reproduced on demand by
`diag_alloc.sh` (job 556515).

| solver | 4 KiB pages (Mcell-updates/s) | 2 MiB huge pages (Mcell-updates/s) | gain from pages alone (× throughput) |
|---|---|---|---|
| baseline | 44.8 | 52.6 | 1.17× |
| opt-L2 | 80.0 | 131.2 | 1.64× |
| opt-L6 | 81.1 | 134.8 | 1.66× |

---

## 6. Auto-vectorisation control (sweep C) — a control that turned out to be vacuous

The question was fair: does the ladder's speedup survive with the compiler's
vectoriser switched off, or were the gains hidden auto-SIMD? The control was
not: **GCC 8.5 does not enable the tree vectoriser at `-O2` in the first
place** (that default arrived in GCC 12), so `-O2 -fno-tree-vectorize`
disabled a pass that was already off. The 1.00× below is a tautology, not a
finding, and an earlier version of this file presented it as evidence. The
table stays as a record of what was actually run.

**Table 6** — sweep C (job 556508), compared against sweep A because both
ran with the default allocator; pairing it against the huge-page sweep would
credit the vectoriser with the page effect.

| solver | -O2 (Mcell-updates/s) | -O2 -fno-tree-vectorize (Mcell-updates/s) | ratio (× throughput) |
|---|---|---|---|
| baseline | 44.8 | 44.8 | 1.00× |
| opt-L6 | 81.1 | 81.1 | 1.00× |

What actually bounds the vectoriser's contribution is sweep F (job 556833),
which built at `-O3` — where the vectoriser IS on — in the same job as an
`-O2` reference, and it answers the original question cleanly in two halves:

- **On the shipped kernels, nothing.** L2 and L6 at `-O3` measure 131.9 and
  134.4 Mcell-updates/s against 131.5 and 135.5 at `-O2` — inside the
  spread, slightly negative for L6. The ladder's gains are the techniques,
  not hidden auto-SIMD, which is what sweep C was meant to establish.
- **On the baseline, a lot — and sweep G proves it is all the vectoriser.**
  `-O3` lifts the baseline 52.6 → 100.7 Mcell-updates/s (×1.92). Sweep G
  (job 556839) paired `-O3` against `-O3 -fno-tree-vectorize` in one job:
  with the vectoriser off, the baseline falls back to **52.6** — the `-O2`
  number, to three digits. The whole `-O3` lift is the vectoriser packing
  the five divisions two-per-`divpd`; every other `-O3` pass is worth
  nothing here.

**Table 6b** — sweep G (job 556839), the control sweep C should have been.
Grid 2048×512, huge-page regime, both legs in one job. Median of 5, in
Mcell-updates/s.

| solver | -O3 (Mcell-updates/s) | -O3 -fno-tree-vectorize (Mcell-updates/s) | vectoriser's share (× throughput) |
|---|---|---|---|
| baseline | 100.7 | 52.6 | 1.92× |
| opt-L2 | 131.0 | 131.7 | 0.99× |
| opt-L6 | 133.7 | 135.0 | 0.99× |

So the question sweep C was built to answer now has a measured answer with a
valid control: **the ladder's gains survive the vectoriser being genuinely
off** (L6 at 135.0 with `-fno-tree-vectorize` at `-O3`, vs 134.8 in the
headline sweep E) — the speedup is the techniques, not hidden auto-SIMD. The
flip side is a finding sweep C could never have seen: on the *unoptimised*
kernel the vectoriser is worth ×1.92, because division-heavy code gives it
real work. For Milestone 1 that is a footnote rather than a temptation —
compiler vectorisation is still SIMD, which the course counts as
parallelism, and the scalar hand-optimised kernel beats `-O3`'s best legal
effort by 31 % anyway.

---

## 7. Correctness

**Table 7** — the correctness gates. "Relative" always means relative gap on
the price (§0).

| gate | what it proves | result |
|---|---|---|
| `tests/test_bs_collapse` | with ξ = 0 and v0 = θ, Heston collapses to Black-Scholes, whose closed form the solver never sees | pass; 1.16e-3 relative against a 2.5e-3 tolerance |
| `tests/test_parity` | put-call parity holds, which is a linearity check on the payoff, the S boundaries and the discounting | pass; residual $5.4e-6 against a $1e-4 tolerance — fractions of a cent |
| `tests/test_opt_matches` | every ladder level and both controls reproduce the baseline | pass; L0/L1 bit-identical, L2+ within 3.3e-15 relative against a 1e-11 tolerance |
| `make check-heston` | **the ξ and ρ terms**, and the interpolated readout, against an external analytic price | pass; 3.3e-4 relative node-aligned, 4.2e-4 relative off-node, tolerance 2e-3 |
| `make convergence` | the answer approaches the analytic price as the grid refines | observed order ≈ 2.2 in dt between the two coarser rungs, then the error flattens at ≈ $0.014 (see below and `results/convergence.png`) |
| `make mc-check` | two methods sharing no code with the grid agree | Fourier $196.1692 (4.2e-6 relative), MC $196.1345 ± $0.5524 (1σ SE) |

On the convergence row: the ladder refines all three axes together (ds and
dv halve, dt quarters per rung, node-alignment preserved), so the observed
order belongs to that combined refinement path, not to dt alone. Between the
two coarser rungs the error against the Fourier anchor drops from $0.534 to
$0.026 — order ≈ 2.2 in dt. Between the two finer rungs it only reaches
$0.014 — order ≈ 0.4 — consistent with a floor that refinement does not
touch (domain truncation at s_max and v_max, and the anchor's own quadrature)
taking over once the discretisation error drops beneath it. The honest
statement of the reference price is therefore $196.17 ± a few cents, and the
ten digits printed by the solver are a reproducibility statement, not an
accuracy claim.

### What put-call parity does NOT prove, and what was added because of it

Parity says `C - P = S e^{-qT} - K e^{-rT}`, an identity that holds in **every**
model. It is therefore blind to the model by construction, and this was measured
rather than assumed. Holding the grid fixed at 421×51×56000 (the `test_parity`
grid: ns×nv×nt) and sweeping the Heston parameters:

**Table 8** — parity is model-blind. ξ is vol-of-vol, ρ the spot–variance
correlation; prices in $; the residual is `C − P − (S e^{-qT} − K e^{-rT})`
in $.

| ξ | ρ | call ($) | put ($) | parity residual ($) |
|---|---|---|---|---|
| 0.01 | +0.00 | 202.56839 | 210.70944 | +5.377e-06 |
| 0.35 | −0.70 | 195.63552 | 203.77657 | +5.377e-06 |
| 0.90 | +0.90 | 171.59651 | 179.73755 | +5.377e-06 |

The price moves by more than $30 and the residual does not change in four
significant figures. Combined with `test_bs_collapse` setting xi = 0, which
zeroes both the cross term `rho*xi*v*S*V_Sv` and the vol-of-vol diffusion
`0.5*xi^2*v*V_vv`, and with `test_opt_matches` only ever comparing the solver
against itself, **the two terms that make this a Heston solver rather than a
Black-Scholes solver had no automated gate at all.**

`make check-heston` closes that. It prices the full parameter set and fails if
the PDE disagrees with Fourier inversion beyond a relative tolerance. Its second
leg deliberately uses an off-node grid (527×100), so the reported price is the
bilinearly interpolated one, which gives the interpolated readout its first
absolute check against an external reference rather than against another copy of
the same solver.

It lives outside `make test` because it needs numpy, and `make test` has to run
on rangpur inside every benchmark job where there is no virtual environment.

`build_and_validate()` in `slurm/bench_common.sh` runs `make test` and **exits
rather than benchmark** a build that has not just re-proved itself. No number in
this file came from an unvalidated binary.

### Two defects the validation work exposed

**The readout was reading the nearest cell.** Today's market point almost never
lands exactly on a grid node, and reading the nearest one cost delta times the
offset, measured at **+4.67 %** on the smoke grid ($205.32 against a true
$195.25). `extract_result` now blends the four surrounding cells, and every Greek
stencil arm is blended the same way so the reported delta really is the
derivative of the reported price. On node-aligned grids the weight is exactly
0.0 and every previously recorded number is unchanged bit for bit, which was
verified rather than assumed.

**The refinement ladder had gone stale against its own config.** The ladder
chose `nv` so that v0 landed on a node when `v_max` was 1.0. P7 moved
`config/reference.cfg` to `v_max = 0.64` for stability (v is a *variance*,
so that ceiling is √0.64 = 80 %/yr volatility — generous for an index but
finite), which put every rung
off-node in the variance direction, and the solver's warning about it was being
written to a stderr stream the script discarded. Separately, the strike was
sitting mid-cell on the two coarsest rungs, so the payoff kink was discretised
differently on each and left a sawtooth on top of the trend.
`scripts/convergence_plot.py` now derives the alignment requirement from the
config at run time and refuses to start if any rung misses, re-raises the
solver's off-node warning instead of swallowing it, and anchors the error
against the semi-analytic Fourier price rather than against its own finest
rung.

---

## 8. Runtime of the full reference solve

**This has no rangpur measurement.** The reference grid is 2048×512×980000, or
1.028e12 cell updates, and a single solve exceeds the 15-minute QOS wall limit
by a wide margin, so it has never been run under `sbatch`.

The figure of 2295 seconds that appears in the project history is a **laptop**
measurement and must not be quoted alongside cluster numbers. Derived from the
measured rangpur throughput in section 2 it would be:

**Table 9** — derived, not measured. The arithmetic, shown once:
1.028e12 cell updates ÷ 52.6e6 cell-updates/s = 19 544 s ≈ 5.4 h.

| solver | measured throughput, huge pages (Mcell-updates/s) | derived runtime (hours) — a derivation, not a measurement |
|---|---|---|
| baseline | 52.6 | ~5.4 |
| opt-L6 | 134.8 | ~2.1 |

If these are used anywhere, they are to be described as derivations.

---

## 9. What had to be learned about rangpur

All four are now encoded in the scripts, and `RUNNING.md` is the operational
home for them.

1. **The `cosc3500` QOS caps wall time at 00:15:00.** A longer job is not
   rejected, it queues forever with reason `QOSMaxWallDurationPerJobLimit`,
   which is a far more confusing failure. The suite is therefore five small
   sweeps rather than one large job.
2. **`sbatch` is only on `PATH` for a login shell.** `ssh host 'sbatch ...'`
   fails with "command not found", so the Makefile prepends `/opt/slurm/bin`.
3. **The partition is heterogeneous**, mixing 24-core Xeon E5-2670 v3 machines
   with 8-core VMs. Every benchmark job pins `--nodelist=r730-2 --exclusive`,
   or the sweeps cannot be compared with each other.
4. **That node is slow, at about a tenth of a laptop core**, and this was
   tested rather than assumed. `slurm/diag_denormal.sh` ruled out denormal
   stalls by showing identical throughput with FTZ/DAZ set, with the grid full
   of NaNs, and with values around 1e33.

The two `diag_*.sh` scripts are kept deliberately. One documents a theory that
was wrong (denormals) and one a theory that was right (page size), and both are
evidence of how the surprising numbers were resolved.

---

## 10. Limitations and threats to validity

Ordered by how much they could move a conclusion. None of these is hidden in
the sections above, but a reader auditing this file should not have to
collect them.

1. **The auto-vectorisation control was vacuous as run** (§6). GCC 8.5 does
   not vectorise at `-O2`, so sweep C measured nothing, and an earlier
   version of this file quoted its 1.00× as a finding. Sweep F now carries
   the real bound (nothing on L2/L6, ×1.92 on the baseline), with sweep G
   attributing the baseline lift.
2. **The `-ffast-math` half of the causal claim remains unmeasured** (§2).
   The iron rule refused the build — under `-ffast-math` level 1 loses the
   bit-identity the test demands — so "the baseline would jump to L2's
   throughput if licensed" is still a prediction, not a measurement. The
   refusal itself is evidence the licence changes answers, but it leaves
   that specific number unknown on this machine.
3. **The published scaling curve still spans two page regimes** (§3), now
   as a physical constraint rather than an allocator accident: a buffer
   smaller than one 2 MiB huge page cannot be THP-backed, so the two small
   grids measure on 4 KiB pages whatever the allocator is told. The
   first run (job 556506), which additionally mixed regimes *within* rows
   and underflowed its largest readout, is retired to `results/superseded/`.
4. **n = 1 at the job level for most sweeps.** Five repetitions run inside
   one process, in one job, on one node, on one day; min–max ranges are
   within-job spreads, not error bars. Two mitigations are now measured
   rather than hoped for: sweep F's `-O2` leg replicated sweep E's headline
   rungs in an independent job to three digits (52.6 / 131.5 / 135.5 vs
   52.6 / 131.2 / 134.8 Mcell-updates/s), and sweep B's pinned re-run
   reproduced the 2048×512 rung again (52.6 / 135.1). The remaining gap:
   everything is still one node and one two-day window. In the unpinned
   sweep A the five reps additionally straddle a **bimodal allocator
   regime** — rep 1 huge-page-backed, reps 2–5 not — so its ranges describe
   a regime gap, not noise; sweep E removes exactly that, and its ranges
   (≲0.5 %) are the ones that behave like noise.
5. **No CPU or NUMA binding.** `--exclusive` keeps co-tenants off the node,
   but nothing pins the process to a core or a memory node (`taskset`,
   `numactl`), and r730-2 is dual-socket. A migration mid-solve, or memory
   on the far socket, is a real confound for a memory-sensitive kernel; the
   tight sweep-E spreads suggest it did not happen here, but that is an
   observation, not a control.
6. **No hardware counters.** The TLB mechanism behind §5 is inferred from
   timing plus the `diag_alloc.sh` hypothesis test (the effect appears and
   disappears with the threshold, on demand). Strong inference, but
   `perf stat -e dTLB-load-misses` would have been direct measurement, and
   it was not available in the job environment.
7. **The clock evidence samples cpu0.** `bench_peak.sh` reads the first
   `cpu MHz` line of `/proc/cpuinfo`, which is not necessarily the core the
   solve ran on. "Pinned at 3100 MHz" is therefore about the node's state,
   not a per-core guarantee.
8. **`dt_stable` as printed omits the advective and cross terms.** The
   stability margin (dt/dt_stable ≈ 0.7) held in every run, so the bound was
   conservative in practice, but the printed number is not the full CFL
   condition.
9. **The convergence ladder co-varies ns, nv and nt.** Legitimate for a
   convergence study along a fixed refinement path, but it means no single
   axis's order is isolated, and the reported ≈ 2.2 belongs to the path.
