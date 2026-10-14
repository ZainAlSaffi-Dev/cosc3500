# RESULTS.md, the measured findings

This file is the tracked record of the numbers. `results/` is gitignored, so
`results/bench_summary.md` (which `make bench-plot` generates from the CSVs)
is missing from a fresh clone. This file is in git, and it gives the
provenance of each number so that any of them can be checked.

**Every timing here came from a rangpur compute node via `sbatch`.** The
provenance stamp is written into each CSV by `write_header()` in
`slurm/bench_common.sh` and reads `host=r730-2.compute.eait.uq.edu.au`, which
is the name of a compute node. `scripts/bench_plot.py` and
`scripts/roofline.py` refuse to render a figure from any CSV whose `host=`
header is not a `*.compute.eait.uq.edu.au` name.

The benchmark node is an Intel Xeon E5-2670 v3 (Haswell, 2014) with 32 KiB
L1d, 256 KiB L2 and 30 MiB L3 (verbatim from `lscpu` in every CSV header).
The compiler is `g++ (GCC) 8.5.0 20210514 (Red Hat 8.5.0-28)`, and the flags
are `-std=c++17 -O2 -Wall -Wextra -Iinclude`. The jobs are 556505 (A),
556506 (B, superseded), 556507 (D), 556508 (C), 556513 (peak), 556515
(alloc) and 556517 (E) on 2026-08-09, and 556833 (F), 556834 (B, pinned
re-run) and 556839 (G) on 2026-08-11.

---

## 0. Nomenclature, units and measurement discipline

Every number in this file carries a unit, and the units are defined once,
here. Memory sizes are binary (KiB and MiB), so a huge page is 2 MiB, which
is slightly more than 2 MB. Money is in US dollars, written `$`. Throughput
is in `Mcell-updates/s`, defined below.

**The instrument.** Everything priced here is one contract from
`config/reference.cfg`, a European call on the S&P 500 index with strike
K = $5250, spot S = $5200, maturity T = 0.25 years, risk-free rate
r = 4.5 %/yr and dividend yield q = 1.3 %/yr, under Heston dynamics with
v0 = θ = 0.04 (a 20 %/yr starting volatility). A "price" is the value of
that option in dollars.

**The grid.** `ns×nv` means ns stock-price nodes by nv variance nodes, and
`nt` is the number of timesteps. The solver keeps two solution buffers of
`ns·nv` doubles each, so the **working set** is `2 · ns · nv · 8` bytes,
which is 16 MiB at the reference footprint of 2048×512. A **cell update** is
one application of the nine-point finite-difference stencil to one (S, v)
node for one timestep, and a full solve performs exactly `ns · nv · nt` of
them.

**Throughput** is reported in **Mcell-updates/s**, the number of millions of
cell updates per second of wall-clock time, `ns · nv · nt / seconds / 10⁶`.
This is the CSV column `cell_updates_per_sec`, which the harness computes
from a `std::chrono::steady_clock` timer around the stepping loop only.
Allocation, payoff setup and readout are excluded, and snapshot I/O is
disabled in every benchmark job. Since Phase 9g the code reads the same
monotonic clock through `omp_get_wtime()`, and what is timed has not
changed. Throughput is a *rate*, so a 400-step benchmark run and a
980 000-step production run measure the same quantity. A **speedup (×)** is
always a ratio of two throughputs, higher meaning faster, taken on the same
grid, the same node and the same page regime. The regime is named wherever
a speedup is quoted, because §5 shows that the page regime alone changes
throughput by 1.66×.

**Uncertainty.** Every throughput is the **median of 5 repetitions** run
inside one process in one slurm job. The quoted `min–max` is the observed
range of those 5 repetitions. It shows their spread within one job and
should not be read as a standard error. Each sweep is a single job on a single node (n = 1
at job level), so the variation between jobs is not characterised, and §10
says what that does and does not permit. The Monte Carlo `±` in §1 and §7
has a different meaning. It is the **1σ standard error of the mean** over
400 000 paths.

**Errors.** A "relative gap" between a and b is `|a − b| / |b|`. Parity
residuals are absolute amounts in $. "ulp" means units in the last place of
an IEEE double. The reported vega is ∂V/∂v0 per unit of **variance**. The
market convention quotes vega per unit of volatility, which differs from it
by a factor of 2√v0, and the code states the unit where it computes vega
(`src/solver.cpp`).

**The sweeps.** Each benchmark job is named by a letter. The operational
detail and the commands to reproduce each sweep are in `RUNNING.md` §3.

| sweep | what it measures | script |
|---|---|---|
| A | the 9-rung optimisation ladder, allocator left alone (mixed page regime) | `slurm/bench_ladder.sh` |
| B | scaling with problem size, baseline vs level 6 (re-run pinned, job 556834) | `slurm/bench_scaling.sh` |
| C | `-fno-tree-vectorize` control, vacuous as run, see §6 | `slurm/bench_novec.sh` |
| D | negative controls at their own nt, paired in-job | `slurm/bench_controls.sh` |
| E | the ladder with huge pages pinned, the headline regime | `slurm/bench_pages.sh` |
| F | compiler-flag control with `-O3` and `-ffast-math` legs (job 556833) | `slurm/bench_flags.sh` |
| G | the valid vectorisation control, `-O3` vs `-O3 -fno-tree-vectorize` (job 556839) | `slurm/bench_novec_o3.sh` |

---

## 1. Headline

**Table 1** gives the findings in one screen. Grids are 2048×512 unless
stated, and speedups are throughput ratios (§0).

| claim | value | source |
|---|---|---|
| Reference price | $196.1683699 (limited by the discretisation, and §7 puts the trustworthy digits at ~$196.17) | `config/reference.cfg` |
| vs semi-analytic Fourier | $196.1692, a relative gap of 4.2e-6 | `make mc-check` |
| vs Monte Carlo, 400k paths | $196.1345 ± $0.5524 (1σ standard error), brackets both | `make mc-check` |
| Serial speedup, baseline to level 6 | **2.56× throughput** huge-page regime / 1.81× default allocator | sweeps E / A |
| of which strength reduction alone | **+148.7 % throughput** | sweep E |
| accuracy cost of that speedup | 3.3e-15 relative on the price, about 15 ulp | `test_opt_matches` |
| levels bit-identical to the baseline | 0 and 1, max_abs_diff exactly $0.0 | `test_opt_matches` |
| page size alone, 2 MiB vs 4 KiB pages | **1.66× throughput** | sweep E vs A, `diag_alloc.sh` |
| traversal order, `ctl-order` control | 2.1× slower (throughput) | sweep E |
| fusing the v=0 row, `ctl-branch` control | no cost, 1.0× | sweep E |
| auto-vectoriser's contribution to the ladder | none, since L6 measures 135.0 with the vectoriser off at `-O3`, vs 134.8 in the headline sweep (sweep C had claimed this but was vacuous, §6) | sweeps F, G |
| what the vectoriser does to the *baseline* | ×1.92 throughput (52.6 → 100.7 Mcell-updates/s at `-O3`), still 31 % short of hand-optimised L2, and back to 52.6 with `-fno-tree-vectorize` | sweep G |
| level 6 vs measured scalar peak | **82 %** of 5.117 GFLOP/s at 2048×512, 95 % at 4096×1024 | `roofline.py` + `bench_peak.sh` |

---

## 2. The optimisation ladder (sweep E, huge-page regime)

A single production solve is a fresh process that allocates once, so it gets
huge pages. This is therefore the regime a real user experiences, and it is the
headline table.

**Table 2** is the ablation ladder from sweep E (job 556517), on grid
2048×512 with nt = 400 and `config/bench.cfg`, in the huge-page regime, in a
single job on r730-2. Throughput is in Mcell-updates/s (§0), and `min–max` is
the range of the 5 reps.

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

The price column of the bench CSVs needs one note. `config/bench.cfg`
shrinks the maturity to 0.00051 years (≈ 4.5 hours) so that a rep fits
within the QOS time limit, which makes the printed price (~$0.0003)
financially meaningless. The price is printed so that the result is used and
the compiler cannot delete the loop (PLAN §6). The grid, the working set and
the stencil traversal are byte-identical to the reference footprint, and a
throughput number depends only on those.

Sweep A (job 556505) is the same ladder with the allocator left alone, where
reps 2 to 5 get their memory from the reused glibc heap on 4 KiB pages.
There the baseline runs at 44.8, level 2 at 80.0 and level 6 at
81.1 Mcell-updates/s, so the overall gain is 1.81×.

### The finding

**Only one rung out of six changed the throughput.** GCC 8.5 at `-O2` had
already performed hoisting, common subexpression elimination,
induction-variable simplification and unrolling, so every rung except level 2
measured within the run-to-run noise.

Level 2 is the only transformation on the list the compiler is **forbidden** to
make, because replacing `x / (2*ds)` with `x * (1/(2*ds))` changes the answer in
the last bits and `-ffast-math` is off.

> At `-O2`, the hand optimisations worth doing are the ones the compiler is not
> allowed to do for you.

This result is visible only because each rung was measured separately. A
single before-and-after pair would not have shown which rung gave the gain.

### Testing the claim (sweep F, job 556833)

The paragraph above is a causal claim, so a job was built to test it. It
rebuilt the code at `-O3`, which enables every standard-conforming
optimisation that `-O2` leaves out, and at `-O2 -ffast-math`, which permits
replacing a division with a multiplication by the reciprocal. Each build ran
in the same job as an `-O2` reference, in the huge-page regime throughout.
The three results follow.

- **The `-O2` leg replicated the headline rungs independently.** It
  measured the baseline at 52.6, L2 at 131.5 and L6 at 135.5 Mcell-updates/s
  in a fresh job two days after sweep E measured 52.6 / 131.2 / 134.8. The
  headline table reproduces across jobs to three digits.
- **At `-O3` the baseline nearly doubled, from 52.6 to 100.7
  Mcell-updates/s (×1.92), while L2 (131.9) and L6 (134.4) did not change.**
  With every standard-conforming optimisation enabled, the compiler closes
  about half the gap to the hand-optimised kernel and stops 31 % short of it.
  Once the divisions have been removed by hand, `-O3` adds nothing. Sweep G
  then attributed the whole gain to the vectoriser (§6, Table 6b).
- **The `-ffast-math` leg was stopped by the validation rule, as the script
  was written to do.** Under `-ffast-math`, level 1, which must be
  bit-identical to the baseline, drifted by 9.1e-13. `test_opt_matches`
  failed, and `build_and_validate()` exited before any benchmark ran. The
  flag that permits the reciprocal substitution also changes answers
  elsewhere, which is the premise of the claim showing up as a test failure.
  The prediction that a fast-math baseline would reach L2's throughput
  therefore remains unmeasured on this machine. The finding is worded as
  "the transformation `-O2` is forbidden to make", because this job did not
  show what the compiler does when it is allowed to make it.

Levels 3 and 4 were predicted null **by construction**, because
`BaselineSolver` was already written with the outer-variance, inner-stock loop
order and already peeled the `v = 0` row into its own loop. The prediction was
confirmed, and the negative controls measured what those rungs would have
been worth.

### Negative controls (sweep E)

The negative controls are deliberately worse code that does arithmetically
identical work. `test_opt_matches` checks both against the same tolerance as
the ladder, so a different answer would point to a bug in the control.

**Table 3** gives the negative controls from sweep E (job 556517), on the
same grid and in the same regime as Table 2. Each control is compared
against opt-L5 measured in the same job.

| control | what it ablates | throughput (Mcell-updates/s, median of 5) | slowdown vs opt-L5 in the same job (× throughput) |
|---|---|---|---|
| opt-ctl-order | loops swapped | 63.6 | 2.1× |
| opt-ctl-branch | v=0 row fused back in | 135.3 | 1.0× |

Both start from level 5 and change exactly one structural thing, and both are
compared against level 5 measured **in the same job at the same nt**, because
a comparison across jobs would also include the effect of a different run
length.

---

## 3. Scaling with problem size (sweep B)

**Table 4** is the scaling sweep, from the pinned re-run (job 556834). The
nt of each grid keeps every row at ~4e8 cell updates, so every row does
equal work. `MALLOC_MMAP_THRESHOLD_` is pinned so that every rep of every row
gets its memory the same way. The 4096×1024 row reads the price deep in the
money (`config/bench_itm.cfg`, the identical grid with the readout point
moved) so that its price keeps significant digits and does not underflow.
The page-regime column was **measured**. The two small grids' throughput
matches sweep A's 4 KiB numbers, and the two large grids match the pinned
huge-page sweeps (E, F). Transparent huge pages need a 2 MiB-aligned 2 MiB
stretch of memory before the kernel promotes it, and only the 8 MiB and
32 MiB buffers had that. The pin therefore gives consistent reps within each
row on every grid, and huge pages only where the buffers are big enough.

| grid (ns×nv) | working set (MiB) | page regime (measured) | baseline (Mcell-updates/s) | opt-L6 (Mcell-updates/s) | speedup (× throughput) |
|---|---|---|---|---|---|
| 512×128 | 1 | 4 KiB | 44.4 | 155.4 | 3.50× |
| 1024×256 | 4 | 4 KiB | 41.8 | 104.4 | 2.49× |
| 2048×512 | 16 | 2 MiB (THP) | 52.6 | 135.1 | 2.57× |
| 4096×1024 | 64 | 2 MiB (THP) | 57.8 | 157.5 | 2.72× |

The table shows two results and has one caveat. The **baseline** is nearly
flat from 1 MiB to 64 MiB of working set, because a division-bound kernel
runs at the same speed whether its operands come from cache or from main
memory. This agrees with the roofline finding in §4. The **level-6 curve**
does depend on the memory system. It dips where the working set first no
longer fits in L2 and is read from the shared L3 (104.4 at 4 MiB), and it
recovers as the bigger buffers get huge pages. The caveat is that the rows
still cover two page regimes. That is now caused by a physical constraint,
since a buffer smaller than one huge page cannot be backed by a transparent
huge page, and the comparisons within one regime are the reliable ones.

The first run of this sweep (job 556506, retired to `results/superseded/`)
did mix the two effects. Nothing pinned the allocator, so rep 1 and reps 2 to
5 ran in different page regimes within single rows. Its 2048×512 row ran in
the 4 KiB regime (44.8 / 81.1, the sweep A numbers), and the out-of-the-money
price of its largest row underflowed to ~$2e-16. Its curve pointed the
opposite way to what the cache sizes predict, and it was published with a
warning. The re-run above replaces it.

---

## 4. Two claims that did not survive being measured

Both belong in the video, because an assumption that was tested and
corrected tells the viewer more than one that was never tested.

### The claim that the solver is memory-bound did not hold

`scripts/roofline.py` counts 31 flops and 24 bytes of compulsory traffic per
cell update directly from the level-6 kernel. That is an arithmetic
intensity of 31/24 ≈ 1.29 flop/byte, which sets the kernel's position on the
x-axis of `results/roofline.png`. `slurm/bench_peak.sh` measures this node's
achievable **single-core scalar double-precision** peak at **5.117 GFLOP/s**
on the node itself.

- Level 6 achieves 4.18 GFLOP/s at the reference grid, **82 % of peak**, and
  4.86 GFLOP/s at 4096×1024, **95 %**, while demanding about 3 GB/s from a
  four-channel DDR4 machine.
- The baseline runs at **27 % of peak** because it is *division*-bound, with
  five `divsd` per cell update.

Optimising the arithmetic therefore moved this kernel from division-bound to
compute-bound, and the next gain has to come from doing more work per
instruction. SIMD does more work per instruction, so it was the natural next
step for Milestone 2.

### `ctl-order` costs 2.1×, and the earlier 30× came from a busy laptop

The 30× was a local measurement taken while the reference solve was running
at the same time, and a kernel that needs a lot of memory bandwidth is the
kind of code that another program on the same machine slows the most. On an
exclusively allocated node it is 2.1×.

Every benchmark job now requests `--exclusive` for this reason. A
measurement can be real and still answer a different question from the one
that was asked.

---

## 5. The page-size effect

The same binary on the same node runs at 134.8 or 81.1 Mcell-updates/s
depending **only on whether the two 8 MiB grid buffers are backed by
2 MiB transparent huge pages or 4 KiB pages.** 16 MiB is 8 pages of 2 MiB,
or 4096 pages of 4 KiB. That needs eight TLB entries with huge pages and
about four thousand with small pages, and the second-level TLB on this core
has 1024 entries.

glibc's dynamic `mmap` threshold means rep 1 of a `--bench` run gets huge pages
and reps 2 to 5 do not, so a plain median-of-5 was silently averaging two
different memory regimes.

The cause was found by eliminating the alternatives. `bench_peak.sh` sampled
the clock during an actual level-6 solve and found it constant at 3100 MHz
throughout, which ruled out throttling. `diag_alloc.sh` then set
`MALLOC_MMAP_THRESHOLD_` by hand and made the effect appear and disappear on
demand.

**Table 5** isolates the page-size effect, with the same binary, the same
grid (2048×512) and the same node. The 4 KiB column is from sweep A (job
556505), and the huge-page column is from sweep E (job 556517) with
`MALLOC_MMAP_THRESHOLD_` pinned to 1 MiB. The two columns come from
different jobs, but on the same day and the same node, and `diag_alloc.sh`
(job 556515) reproduced the effect on demand.

| solver | 4 KiB pages (Mcell-updates/s) | 2 MiB huge pages (Mcell-updates/s) | gain from pages alone (× throughput) |
|---|---|---|---|
| baseline | 44.8 | 52.6 | 1.17× |
| opt-L2 | 80.0 | 131.2 | 1.64× |
| opt-L6 | 81.1 | 134.8 | 1.66× |

---

## 6. Auto-vectorisation control (sweep C), which turned out to be vacuous

The question was whether the ladder's speedup survives with the compiler's
vectoriser switched off, or whether the gains came from SIMD code that the
compiler generated by itself. The question was sound, and the control was
flawed. **GCC 8.5 does not enable the tree vectoriser at `-O2` at all**
(that default arrived in GCC 12), so `-O2 -fno-tree-vectorize` disabled a
pass that was already off. The 1.00× below is therefore true by construction
and measures nothing, and an earlier version of this file presented it as
evidence. The table stays as a record of what was actually run.

**Table 6** is sweep C (job 556508), compared against sweep A because both
ran with the default allocator. Pairing it with the huge-page sweep would
have counted the page effect as a vectoriser effect.

| solver | -O2 (Mcell-updates/s) | -O2 -fno-tree-vectorize (Mcell-updates/s) | ratio (× throughput) |
|---|---|---|---|
| baseline | 44.8 | 44.8 | 1.00× |
| opt-L6 | 81.1 | 81.1 | 1.00× |

Sweep F (job 556833) measures the vectoriser's contribution. It built at
`-O3`, where the vectoriser is on, in the same job as an `-O2` reference,
and it answers the original question in two parts.

- **On the shipped kernels the vectoriser adds nothing.** L2 and L6 at `-O3`
  measure 131.9 and 134.4 Mcell-updates/s against 131.5 and 135.5 at `-O2`,
  which is inside the spread and slightly lower for L6. The ladder's gains
  therefore come from the techniques themselves, and sweep C was meant to
  establish exactly this.
- **On the baseline the gain is large, and sweep G shows that all of it
  comes from the vectoriser.** `-O3` raises the baseline from 52.6 to
  100.7 Mcell-updates/s (×1.92). Sweep G (job 556839) paired `-O3` against
  `-O3 -fno-tree-vectorize` in one job. With the vectoriser off, the baseline
  falls back to **52.6**, the `-O2` number to three digits. The whole `-O3`
  gain comes from the vectoriser computing the five divisions two at a time
  with `divpd`, and no other `-O3` pass changes the rate here.

**Table 6b** is sweep G (job 556839), the control that sweep C should have
been, on grid 2048×512 in the huge-page regime with both builds in one job.
Each value is the median of 5, in Mcell-updates/s.

| solver | -O3 (Mcell-updates/s) | -O3 -fno-tree-vectorize (Mcell-updates/s) | vectoriser's share (× throughput) |
|---|---|---|---|
| baseline | 100.7 | 52.6 | 1.92× |
| opt-L2 | 131.0 | 131.7 | 0.99× |
| opt-L6 | 133.7 | 135.0 | 0.99× |

The question sweep C was built to answer now has a measured answer with a
valid control. **The ladder's gains remain with the vectoriser switched
off** (L6 at 135.0 with `-fno-tree-vectorize` at `-O3`, against 134.8 in the
headline sweep E), so the speedup comes from the techniques. Sweep G also
shows something sweep C could not. On the *unoptimised* kernel the
vectoriser is worth ×1.92, because the division-heavy code gives it real
work to vectorise. For Milestone 1 this result stays a side note and was not
adopted. Compiler vectorisation is still SIMD, which the course counts as
parallelism, and the scalar hand-optimised kernel is 31 % faster than the
best that `-O3` does with the baseline in any case.

---

## 7. Correctness

**Table 7** lists the correctness gates. "Relative" always means the
relative gap on the price (§0).

| gate | what it proves | result |
|---|---|---|
| `tests/test_bs_collapse` | with ξ = 0 and v0 = θ, Heston collapses to Black-Scholes, whose closed form the solver never sees | pass, 1.16e-3 relative against a 2.5e-3 tolerance |
| `tests/test_parity` | put-call parity holds, which is a linearity check on the payoff, the S boundaries and the discounting | pass, a residual of $5.4e-6 against a $1e-4 tolerance, a small fraction of a cent |
| `tests/test_opt_matches` | every ladder level and both controls reproduce the baseline | pass, L0/L1 bit-identical, L2+ within 3.3e-15 relative against a 1e-11 tolerance |
| `make check-heston` | **the ξ and ρ terms**, and the interpolated readout, against an external analytic price | pass, 3.3e-4 relative node-aligned, 4.2e-4 relative off-node, tolerance 2e-3 |
| `make convergence` | the answer approaches the analytic price as the grid refines | observed order ≈ 2.2 in dt between the two coarser rungs, then the error flattens at ≈ $0.014 (see below and `results/convergence.png`) |
| `make mc-check` | two methods sharing no code with the grid agree | Fourier $196.1692 (4.2e-6 relative), MC $196.1345 ± $0.5524 (1σ SE) |

The convergence row needs some explanation. The refinement ladder refines
all three axes together (ds and dv halve and dt quarters per rung, with node
alignment kept), so the observed order describes that combined refinement
path and does not isolate dt. Between the two coarser rungs the error
against the Fourier anchor drops from $0.534 to $0.026, an order of ≈ 2.2 in
dt. Between the two finer rungs it only falls to $0.014, an order of ≈ 0.4.
That is consistent with an error floor that refinement cannot reduce
(domain truncation at s_max and v_max, and the anchor's own quadrature),
which becomes the largest error once the discretisation error drops below
it. The reference price is therefore $196.17 ± a few cents. The ten digits
the solver prints show that a run reproduces exactly, and they do not claim
ten digits of accuracy.

### What put-call parity cannot prove, and what was added because of it

Parity says `C - P = S e^{-qT} - K e^{-rT}`, an identity that holds in **every**
model. It is therefore insensitive to the model by construction, and this
was checked by measurement. Table 8 holds the grid fixed at 421×51×56000
(the `test_parity` grid, ns×nv×nt) and sweeps the Heston parameters.

**Table 8** shows that parity does not depend on the model. ξ is
vol-of-vol, ρ is the correlation between spot and variance, prices are in
$, and the residual is `C − P − (S e^{-qT} − K e^{-rT})` in $.

| ξ | ρ | call ($) | put ($) | parity residual ($) |
|---|---|---|---|---|
| 0.01 | +0.00 | 202.56839 | 210.70944 | +5.377e-06 |
| 0.35 | −0.70 | 195.63552 | 203.77657 | +5.377e-06 |
| 0.90 | +0.90 | 171.59651 | 179.73755 | +5.377e-06 |

The price moves by more than $30 and the residual does not change in four
significant figures. `test_bs_collapse` sets xi = 0, which zeroes both the
cross term `rho*xi*v*S*V_Sv` and the vol-of-vol diffusion
`0.5*xi^2*v*V_vv`, and `test_opt_matches` only ever compares the solver
against itself. So **the two terms that distinguish a Heston solver from a
Black-Scholes solver had no automated test at all.**

`make check-heston` fills that gap. It prices the full parameter set and
fails if the PDE disagrees with Fourier inversion beyond a relative
tolerance. Its second leg deliberately uses an off-node grid (527×100), so
the reported price is the bilinearly interpolated one. That gives the
interpolated readout its first absolute check against an external
reference, where every earlier check compared the solver with another copy
of itself.

It is kept outside `make test` because it needs numpy, and `make test` has
to run on rangpur inside every benchmark job where there is no virtual
environment.

`build_and_validate()` in `slurm/bench_common.sh` runs `make test` and
**exits without benchmarking** any build that has not passed the tests in
the same job. No number in this file came from an unvalidated binary.

### Two defects the validation work exposed

**The readout was reading the nearest cell.** Today's market point almost
never falls exactly on a grid node, and reading the nearest node caused an
error of delta times the offset, measured at **+4.67 %** on the smoke grid
($205.32 against a true $195.25). `extract_result` now interpolates
bilinearly between the four surrounding cells, and every point of the Greek
stencils is interpolated the same way, so the reported delta is the
derivative of the reported price. On node-aligned grids the interpolation
weight is exactly 0.0, and a check confirmed that every previously recorded
number is unchanged bit for bit.

**The refinement ladder no longer matched its own config.** The ladder chose
`nv` so that v0 fell on a node when `v_max` was 1.0. P7 moved
`config/reference.cfg` to `v_max = 0.64` for stability (v is a *variance*,
so that ceiling is √0.64 = 80 %/yr volatility, which is generous for an
index and still finite). That put every rung off-node in the variance
direction, and the solver's warning about it was written to a stderr stream
the script discarded. Separately, the strike sat in the middle of a cell on
the two coarsest rungs, so the payoff kink was discretised differently on
each and added a sawtooth to the trend. `scripts/convergence_plot.py` now
derives the alignment requirement from the config at run time and refuses
to start if any rung is off-node. It also raises the solver's off-node
warning again where it used to discard it, and it now measures the error
against the semi-analytic Fourier price, where it used to use its own finest
rung.

---

## 8. Runtime of the full reference solve

**The full reference solve has no rangpur measurement.** The reference grid
is 2048×512×980000, or 1.028e12 cell updates, and a single solve exceeds the
15-minute QOS wall limit by a wide margin, so it has never been run under
`sbatch`.

The figure of 2295 seconds that appears in the project history is a **laptop**
measurement and must not be quoted alongside cluster numbers. Table 9 gives
the runtime derived from the measured rangpur throughput in section 2.

**Table 9** is derived from the measured throughput. For the baseline the
arithmetic is 1.028e12 cell updates ÷ 52.6e6 cell-updates/s = 19 544 s
≈ 5.4 h.

| solver | measured throughput, huge pages (Mcell-updates/s) | derived runtime (hours), derived from the measured throughput |
|---|---|---|
| baseline | 52.6 | ~5.4 |
| opt-L6 | 134.8 | ~2.1 |

If these are used anywhere, they are to be described as derivations.

---

## 9. What had to be learned about rangpur

All four are now built into the scripts, and `RUNNING.md` documents them for
day-to-day use.

1. **The `cosc3500` QOS caps wall time at 00:15:00.** A longer job is never
   rejected. It stays in the queue with reason
   `QOSMaxWallDurationPerJobLimit`, which makes the failure far harder to
   diagnose. The suite is therefore five small sweeps.
2. **`sbatch` is only on `PATH` for a login shell.** `ssh host 'sbatch ...'`
   fails with "command not found", so the Makefile prepends `/opt/slurm/bin`.
3. **The partition is heterogeneous**, mixing 24-core Xeon E5-2670 v3 machines
   with 8-core VMs. Every benchmark job requests
   `--nodelist=r730-2 --exclusive` so that the sweeps can be compared with
   each other.
4. **That node is slow, at about a tenth of a laptop core**, and this was
   tested. `slurm/diag_denormal.sh` ruled out denormal stalls by showing
   identical throughput with FTZ/DAZ set, with the grid full of NaNs, and
   with values around 1e33.

The two `diag_*.sh` scripts are kept deliberately. One documents a theory that
was wrong (denormals) and one a theory that was right (page size), and both are
evidence of how the surprising numbers were resolved.

---

## 10. Limitations and threats to validity

The limitations are ordered by how much they could change a conclusion. Each
one is also stated in the sections above, and this list collects them in one
place for a reader checking this file.

1. **The auto-vectorisation control was vacuous as run** (§6). GCC 8.5 does
   not vectorise at `-O2`, so sweep C measured nothing, and an earlier
   version of this file quoted its 1.00× as a finding. Sweep F now gives
   the real bound (nothing on L2/L6, ×1.92 on the baseline), and sweep G
   attributes the baseline's gain to the vectoriser.
2. **The `-ffast-math` part of the causal claim remains unmeasured** (§2).
   The validation rule stopped the build, because under `-ffast-math`
   level 1 loses the bit-identity the test demands. The statement that the
   baseline would reach L2's throughput under `-ffast-math` therefore
   remains a prediction. The refusal itself shows that the flag changes
   answers, but it leaves that specific number unknown on this machine.
3. **The published scaling curve still spans two page regimes** (§3). This
   now comes from a physical constraint. A buffer smaller than one 2 MiB
   huge page cannot be backed by a transparent huge page, so the two small
   grids run on 4 KiB pages whatever the allocator is told. The first run
   (job 556506), which also mixed regimes *within* rows and underflowed its
   largest readout, is retired to `results/superseded/`.
4. **n = 1 at the job level for most sweeps.** Five repetitions run inside
   one process, in one job, on one node, on one day, so the `min–max` ranges
   show the spread within one job and say nothing about the variation
   between jobs. Two checks across jobs have now been measured. Sweep F's
   `-O2` leg replicated sweep E's headline rungs in an independent job to
   three digits (52.6 / 131.5 / 135.5 against 52.6 / 131.2 / 134.8
   Mcell-updates/s), and sweep B's pinned re-run reproduced the 2048×512
   rung again (52.6 / 135.1). Everything still comes from one node within
   two days. In the unpinned sweep A the five reps also fall into a
   **bimodal allocator regime**, with rep 1 backed by huge pages and reps
   2 to 5 on small pages, so its ranges show the gap between the two
   regimes. Sweep E removes that, and its ranges (≲0.5 %)
   behave like noise.
5. **No CPU or NUMA binding.** `--exclusive` keeps other jobs off the node,
   but nothing binds the process to a core or a memory node (`taskset`,
   `numactl`), and r730-2 has two sockets. A migration to another core
   during a solve, or memory on the far socket, could change the rate of a
   memory-sensitive kernel. The tight sweep E spreads suggest it did not
   happen here, but that is an observation, and no control was run for it.
6. **No hardware counters.** The TLB mechanism behind §5 is inferred from
   timing plus the `diag_alloc.sh` hypothesis test (the effect appears and
   disappears with the threshold, on demand). The inference is strong, but
   `perf stat -e dTLB-load-misses` would have measured the misses directly,
   and it was not available in the job environment.
7. **The clock evidence samples cpu0.** `bench_peak.sh` reads the first
   `cpu MHz` line of `/proc/cpuinfo`, which may belong to a different core
   from the one the solve ran on. The constant 3100 MHz reading therefore
   describes the node's state and does not guarantee the clock of each core.
8. **`dt_stable` as printed omits the advective and cross terms.** The
   stability margin (dt/dt_stable ≈ 0.7) held in every run, so the bound was
   conservative in practice, but the printed number leaves out part of the
   CFL condition.
9. **The convergence ladder co-varies ns, nv and nt.** This is legitimate
   for a convergence study along a fixed refinement path, but no single
   axis's order is isolated, and the reported ≈ 2.2 describes the path.

---

# Milestone 2 results

Every number below comes from an exclusive r730 node (dual Xeon E5-2670 v3,
2 x 12 cores, no hyper-threading, two NUMA nodes, 30 MiB L3 per socket,
g++ 8.5.0, Open MPI 4.1.1) through `sbatch`. Each is the median of 5 reps,
with huge pages pinned (`MALLOC_MMAP_THRESHOLD_=1048576`), denormals
flushed to zero (`slurm/flush_denormals.c`, see M2.7), and
`OMP_PLACES=cores OMP_PROC_BIND=close` unless a row says otherwise. The CSVs
carry the full provenance in their headers. `make m2-figures` rebuilds
every figure in `results/m2/` from exactly these files. Rates are in
Mcell-updates/s (million cell updates per second) on the 2048 by 512 grid
unless stated.

The first round of sweeps (jobs 651147 to 651151 and 651154) ran without
flushing denormals and is kept in `results/superseded/`. It is not quoted
here, for the reason in M2.7.

From Phase 13 on, every rung that uses the AVX2 row was measured with the
split turned off. Given only `-mavx2`, g++ 8.5 compiles each unaligned
256 bit load in that row as two 128 bit loads joined by `vinsertf128`, and
each unaligned store as two 128 bit stores (E2 in M2.6). The shipped build
now adds `-mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store`
to the AVX2 file, so simd, omp, mpi and every experiment built on them run
with unsplit 256 bit loads and stores. The tables quote the Phase 13 reruns
(jobs 659238 to 659249 and 659436), which give every rung the same answers
as before to every printed digit. Numbers measured with the split on appear
only where a comparison needs them, mainly in E2, E2b and E1. The baseline,
serial and omp_scalar rungs do not use the AVX2 row. The omp_scalar sweep (job
651297), the optimisation levels (job 652696), the peak measurement (job
651302) and the profile (job 652697) were therefore not rerun, and only the
profile's simd row in M2.11 still comes from the build with the split on.

## M2.1 The ladder (job 659238, figure F1)

| Rung | Rate | Gain over the rung below | Total over baseline |
|---|---|---|---|
| baseline | 57 | | 1.0x |
| serial | 159 | 2.80x | 2.8x |
| simd (AVX2 and FMA) | 456 | 2.87x | 8.0x |
| omp, 24 threads | 6367 | 14.0x | 112x |
| mpi, 2 ranks x 12 threads, one node | 6984 | 1.10x | 123x |

Every rung matches the baseline to 3.3e-15 or better, and omp and mpi match
simd bit for bit (test log `slurm-659238.out`, figure F9).

## M2.2 Strong scaling (jobs 659239 and 651297, figure F2)

| Threads | omp (AVX2 row) | speedup | omp_scalar | speedup |
|---|---|---|---|---|
| 1 | 454 | 1.0 | 158 | 1.0 |
| 12 | 4354 | 9.6 | 1603 | 10.1 |
| 24 | 6311 | 13.9 | 2922 | 18.5 |

Amdahl fits give a serial fraction of 2.42% for the AVX2 row and 1.57% for
the scalar row. The boundary work alone is 0.44% of a one-thread step, so
most of the fitted fraction comes from other costs. The measured serial boundary
time grows from 0.010 ms per step on one thread to 0.015 ms on 24, because
the main thread then has to fetch rows that other cores have just written.
At 24 threads it is 9.1% of the step. A faster row shortens every step while
the boundaries and the fork and barrier cost the same, so the fitted
fraction grows as the row gets faster.

The AVX2 row scales worse than the scalar row at every thread count above 2
(efficiency 58% against 77% at 24). At 24 threads the AVX2 kernel moves about
6.3e9 x 24 bytes = 151 GB/s through the cache hierarchy, 2.2 times the
measured DRAM triad bandwidth, so it is running out of shared cache and
memory bandwidth before the scalar row does. The SIMD gain still carries
over to 24 threads, since the AVX2 rung is 2.16 times the scalar rung there.

On this node `OMP_PROC_BIND=close` alternates sockets from the second
thread, because the CPUs are numbered even on socket 0 and odd on socket 1
(job 651500 printed each thread's CPU). So every even thread count in this
sweep is split evenly over the two sockets.

## M2.3 Weak scaling (job 659240, figure F3)

Each thread gets 21 interior rows, and `vmax` grows with the rows so the
spacing stays the same. The step time rose from 0.1124 ms on 1 thread to
0.1658 ms on 24, so the scaled speedup on 24 threads is
24 x 0.1124 / 0.1658 = 16.3 against an ideal 24. At 24 threads the serial boundaries take 9.2% of the step, and
Gustafson's law with that serial share predicts 24 - 0.092 x 23 = 21.9, so
the measured 16.3 also falls short of Gustafson's line. The one-thread grid
(23 rows, 370 KiB per buffer) is larger than a core's 256 KiB L2, so both
ends of the sweep work from L3. The one thread has a whole socket's L3 and
memory bandwidth to itself, while the 24 threads share them, which favours
the one-thread point.

## M2.4 Grid size (job 659241, figure F4)

| Grid | Working set | serial | simd | omp 24 | omp over serial |
|---|---|---|---|---|---|
| 512 x 128 | 1 MiB | 145 | 422 | 1228 | 8x |
| 1024 x 256 | 4 MiB | 153 | 451 | 4264 | 28x |
| 2048 x 512 | 16 MiB | 157 | 454 | 6278 | 40x |
| 2896 x 724 | 32 MiB | 58 | 169 | 3288 | 56x, all three slow |
| 4096 x 1024 | 64 MiB | 156 | 392 | 2990 | 19x |
| 8192 x 1024 | 128 MiB | 153 | 348 | 2259 | 15x |

Job 659241 reran this sweep with the split turned off. simd gained 8% to 29%
at every grid except 2896 x 724, where it stayed at 169, and serial did not
change. omp on 24 threads gained 4% and 10% on the 4 and 16 MiB grids,
stayed within the spread of its reps on 512 x 128 and 4096 x 1024, and rose
from 1709 to 2259 on 8192 x 1024, inside the 1534 to 2260 range of the first
run (job 651299). At 2896 x 724 it fell from 6499 to 3288, with all five reps
between 3281 and 3328 against 6360 to 6748 in both earlier runs. A control
job (659437) ran the whole sweep again on the same day with the split turned
back on. It gave 3267 at 2896 x 724 and matched the earlier runs everywhere
else (omp 5712 against 5711 at 2048 x 512), so the drop does not come from
the flags. The prices match the earlier runs to every printed digit.

The serial rate changes little with grid size. simd loses 14% at
4096 x 1024 and 23% at 8192 x 1024. At 454 Mcell-updates/s it streams about
10.9 GB/s, close to the 12.4 GB/s triad bandwidth of one core, so one core's
memory bandwidth starts to limit it once the grid no longer fits in L3. The
24-thread rung peaks while the grid fits in the two L3 caches (60 MiB) and
drops to about half that rate once it does not. On small grids the fork and
barrier of every step take a large share of the time, because they cost the
same whatever the grid size. Every rung is slow at 2896 x 724 in this run.
That grid puts its two buffers exactly 2^24 bytes apart (M2.8), so its row
does not measure the effect of size, and its 56x does not reflect a real
speedup.

## M2.5 Roofline (jobs 651302, 659238, 659241, figure F5)

`bench_peak.sh` measured, on r730-2, a scalar peak of 5.12 GFLOP/s, an AVX2
FMA peak of 43.4 GFLOP/s on one core and 971 on 24, and a triad bandwidth of
12.4 GB/s on one thread and 67.9 GB/s on 24. The kernel does 31 flops per
cell update over 24 bytes, so the DRAM roof at its intensity of 1.29 is
about 88 GFLOP/s, and about 16 GFLOP/s for one core. The serial rung
reaches 159 x 31 = 4.9 GFLOP/s, 96% of the scalar peak, so the next gain
had to come from wider instructions. omp on 24 threads reaches 195 GFLOP/s
at 2048 by 512, 2.2 times the DRAM roof, because the grid is in L3. At 4096
by 1024 it is 93 GFLOP/s, about the DRAM roof. PLAN §5 hypothesis 3 held.
The kernel reads its data from L3 while the grid fits there, and it becomes
memory-bound once the grid no longer fits in L3. simd on one core
reaches 14.1 GFLOP/s at 2048 by 512 and 12.1 at 4096 by 1024, at most a
third of the one-core AVX2 FMA peak.

The first peak job (651154) reported only 24.8 GB/s on 24 threads, because
its arrays were `std::vector`, zero-filled on the main thread. Under the
first-touch policy Linux places each page in the memory of the socket whose
thread first writes it, so every page was placed in socket 0's memory.
Allocating with `new double[n]` and writing each array first from all the
threads in parallel gave 67.9 GB/s.

## M2.6 Experiments (jobs 659246, 659244, 659245, 659247, figures F1 side panel, F7 and F9)

Every row is at 24 threads on 2048 by 512 unless stated, relative to the omp
rung in the same job. The rows come from the Phase 13 reruns with the split turned
off, except E2 and E2b, which compare the two builds. The same sweeps
measured with the split on (jobs 651300, 651301, 651474 and 652647, and
their reruns 656960 to 656963) stay in `results/`.

| Id | Variant | Result |
|---|---|---|
| E4 | one parallel region for the whole solve | 1.02x at 2048 x 512, inside the spread of the reps, 1.13x at 512 x 128 |
| E5 | `schedule(dynamic,1)`, `(dynamic,8)`, `guided` | 0.57x, 0.63x, 0.78x, static is fastest |
| E8 | v = 0 row and stock edges inside the parallel loop | 1.06x at 2048 x 512, 1.10x at 512 x 128, serial share 9.4% to 2.3% |
| E6 | 12 threads `spread` (socket 0 only) against `close` (both sockets), 4096 x 1024 | 1894 against 3413 |
| E6 | unbound threads (`OMP_PROC_BIND=false`), 24 threads | 1.11x at 2048 x 512, 1.13x at 4096 x 1024, not explained |
| E7 | `numactl --interleave=all`, 24 threads | 1.10x at 2048 x 512, 1.74x at 4096 x 1024, and far less spread between reps |
| E11 | glibc default allocator against huge pages pinned | omp 24: 4154 against 6556 (0.63x), serial 89 against 159 (0.56x). Only the first rep of each run gets huge pages |
| E12 | denormals not flushed (job 659436) | see M2.7 |
| E3 | float buffers, `omp_float` (job 659247) | 1 thread 871 against simd 453 (1.92x), 24 threads 11043 against 6375 at 2048 x 512 (1.73x), 14353 against 2975 at 4096 x 1024 (4.82x). Rejected for accuracy, see below |
| E2 | alignment: column 1 moved 0 to 7 doubles along its cache line (`omp_align_loadu`), aligned load instructions (`omp_align_load`), east and west built by register shuffles (`omp_align_shuffle`) (jobs 657421, 657422) | 1 thread against simd 387: centre aligned 368, nothing aligned 383, aligned loads 425 (1.10x), shuffles 461 (1.19x). 24 threads against omp 5744: aligned loads 6062 (1.06x), shuffles 6161 (1.07x). The cause is the compiler, see below |
| E2b | E2 rebuilt with `-mno-avx256-split-unaligned-load` and `-store` on the AVX2 file (jobs 657498, 657499) | simd 460 against 387 in E2 (1.19x), omp 6431 against 5744 (1.12x). 1 thread: centre aligned 477, west or east aligned 439, nothing aligned 412, aligned loads 477, shuffles 460. 24 threads: centre aligned 6579 (1.02x over omp) |
| E1 | the whole program rebuilt under other compiler flags, one thread (jobs 659248, 659249) | `-O3 -mavx2 -mfma` takes the scalar serial rung to 346, and to 363 once its own loads are also left unsplit, against the hand-written simd 450 (0.81x). Adding `-ffast-math` takes serial to 468 (1.04x over simd) and the baseline from 58 to 380 (0.85x of simd). Answers move by at most 1.1e-14, see below |
| E9 | `mpi_overlap`, edge rows sent first and the other rows updated while the messages are in transit (job 659243) | 2 nodes, 1 rank x 24 each, 2048 x 512: 8009 against 5693 (1.41x), exchange share 39% down to 19%. Every other layout within 4%, see below |

E5 runs the omp solver with `schedule(runtime)`, so `OMP_SCHEDULE` picks the
schedule. Its own control, the same solver with `OMP_SCHEDULE=static`, ran at
6092 against the omp rung's 6374 in the same job. Against that control the
three schedules are 0.60x, 0.66x and 0.82x, so the static schedule is faster
against both references.

E4 and E8 give small gains, between 2% and 13% in this run. The persistent
region at 512 x 128 gained 12.5% here, against 15.5% and 5.4% in the two runs
with the split on, so the size of that gain changes from run to run. The
24 thread rows on 4096 x 1024 in E6 range from 2699 to 4009 Mcell/s within
one job, so the unbound ratio there is noise. E7 confirms the same page
placement problem in the solver. The main thread writes every page of the
grid first, so all of the grid is placed in socket 0's memory, and
interleaving the pages across both sockets gives the big grid 74% more in
this run. The interleaved runs are steady at 5022 to 5080 Mcell/s. Without
interleaving, the runs on the big grid reach either about 2800 or about 4000
(2768 to 4064 here), as they did in the two runs with the split on, which
gave 1.28x. The size of the gain therefore depends on which of the two rates
those runs reach.

E2 found that the compiler was splitting the loads. With only
`-mavx2`, g++ 8.5 tunes for a generic x86 and compiles every unaligned
256 bit load (`_mm256_loadu_pd`) as two 128 bit loads joined by
`vinsertf128`, and every unaligned store as two 128 bit stores. Aligned loads
are never split. The shipped row has 11 split loads per four cells. The
aligned load row has 6 and the shuffle row none, and the fewer split loads
account for all of their gain in E2 (disassembly in
`results/e2_disasm_657421.txt`). E2b turns the split off, and the shipped
row then runs at 460 Mcell/s on one core against 387, with no change to the
source. The aligned instruction gives exactly the same speed as the
unaligned one at the same address (477 and 477), as L04 says. Where the data
sits in its cache lines still changes the rate. A row whose centre load is
32 byte aligned runs at 477, one whose west or east load is aligned runs at
439, and one with no aligned load runs at 412, because each unaligned
32 byte load crosses a cache line boundary on every other iteration. The
shuffle row (460) is no faster once the loads are unsplit. At 24 threads
the gaps shrink to a few percent, and on 4096 x 1024 the reps spread more
than the variants differ. Phase 13 adopted the two flags in the shipped
build, and every rung and experiment above was measured again with them.

E3 is the one experiment that is much faster and was still not kept. On the
reference grid (2048 x 512, 980000 steps, 24 threads, one run each) the
float solver took 88.0 s against 149.9 s, but gave a price of 196.1256
against the double 196.1684 and the Fourier price 196.1692. The float error
of $0.044 is larger than the double error of the 841 x 65 grid in §7
($0.026). Gamma moved by 11.9%, vega by 2.3% and delta by 0.46%. On the
128 x 32 test grid with 5200 steps the worst float gap from the baseline is
8.2e-6 (put) and 2.1e-6 (call), within its own 1e-4 limit in
`test_ladder_matches`. The 4096 x 1024 gain is more than 2x because the
float buffers (32 MB) fit in the two L3 caches (60 MB) and the double
buffers (64 MB) do not.

E1 asks whether the compiler could have written the AVX2 rung by itself.
The whole program was rebuilt with each set of flags in place of `-O2`, and
each rung was timed on one thread at 2048 x 512 (median of five reps,
Mcell-updates/s). The last column is the largest relative change in the
price, delta, gamma or vega against the `-O2` build, printed to 17 digits.
In every build the AVX2 file also gets the two shipped flags that keep its
unaligned loads unsplit, so the simd column measures the shipped row. The two
rows marked without split loads add the same flags to the whole program,
which gives the compiler's own AVX2 code the same advantage.

| Build | baseline | serial | simd | largest change |
|---|---|---|---|---|
| `-O2` (shipped) | 57.6 | 157.0 | 449.6 | 0 |
| `-O3` | 112.5 | 241.3 | 446.9 | 0 |
| `-O3 -mavx2 -mfma` | 107.4 | 345.6 | 451.1 | 1.1e-14 |
| `-O3 -mavx2 -mfma` without split loads | 106.4 | 362.7 | 452.0 | 1.1e-14 |
| `-O3 -march=native` | 105.9 | 266.3 | 455.8 | 1.1e-14 |
| `-O2 -ffast-math` | 58.6 | 175.1 | 476.5 | 1.1e-14 |
| `-O3 -ffast-math` | 189.5 | 271.5 | 472.7 | 1.1e-14 |
| `-O3 -mavx2 -mfma -ffast-math` | 380.1 | 468.2 | 478.0 | 1.1e-14 |
| `-O3 -mavx2 -mfma -ffast-math` without split loads | 376.8 | 399.3 | 478.2 | 1.1e-14 |

Given AVX2 and FMA, the compiler vectorises the scalar row four lanes wide
and reaches 0.77 times the hand-written row, or 0.81 times once its own
unaligned loads are also left unsplit. It cannot prove that the input
and output rows do not overlap, so it also emits a run-time aliasing check
and a scalar fallback. The baseline gains nothing from AVX2 because it is
bound by its divisions, and Haswell executes a 256 bit divide as two 128 bit
divides. `-ffast-math` lets the compiler treat floating point arithmetic as
exact. It then computes the reciprocal of each grid spacing once and
multiplies inside the loop. The baseline's vectorised inner loop goes from 10
divide instructions to none, and the naive code reaches 380, 0.85 times the
hand-written rung. The serial rung already multiplies by precomputed
reciprocals, and it reaches 468 by regrouping its arithmetic, 1.04 times the
hand-written rung. Built without split loads it runs at 399 instead, and the
reason has not been found. Under `-ffast-math` the hand-written row itself
gains 5% to 6% (473 to 478 against 450). Its machine code is three 256 bit
instructions shorter (45 against 48), so the compiler regrouped some of the
intrinsics' arithmetic as well. The answers move by at most 1.1e-14.
`-ffast-math` relaxes the IEEE rules for the whole program, however, and lets
the compiler assume no value is ever NaN or infinite, while the tests detect
an unstable run by exactly such a value. The shipped build therefore does
not use it. `-march=native` tunes for Haswell, which also keeps each
unaligned load as one 256 bit load, as the two flags do, so simd stays at
456. Why it slows the auto-vectorised serial rung (266 against 363 with the
two flags) has not been found.

E9 asks whether the edge row exchange can be overlapped with the update.
`mpi_overlap` updates the two edge rows first, and the master thread then
posts their non-blocking sends and receives while the other threads update
the remaining rows. It calls `MPI_Waitall` only at the end of the step.
Both solvers ran back to back in every layout of M2.9 (median
Mcell-updates/s, with the exchange's share of a step in brackets).

| Layout | 2048 x 512, mpi | overlap | 4096 x 1024, mpi | overlap |
|---|---|---|---|---|
| 1 node, 1 rank x 24 | 6288 | 6068 | 3045 | 2988 |
| 1 node, 2 ranks x 12 | 7025 (22%) | 6983 (20%) | 6664 (11%) | 6583 (14%) |
| 2 nodes, 1 rank x 24 each | 5693 (39%) | 8009 (19%) | 9262 (31%) | 9332 (30%) |
| 2 nodes, 2 ranks x 12 each | 9364 (39%) | 9259 (38%) | 10848 (31%) | 10878 (29%) |

Only the layout with one rank on each of two nodes gains, by 1.41 times on
the small grid, and the ranges of its reps do not overlap (5480 to 6133
against 7688 to 8282). Everywhere else the rate moves by less than 4% and
the exchange share changes by at most 3 percentage points. The likely cause
is Open MPI's message protocol. Between nodes the TCP transport has an eager
limit of 64 KiB, so any message up to that size is sent as soon as
`MPI_Isend` is called. Between ranks on the same node the shared memory
transport has an eager limit of 4 KiB. A 16 or 32 KiB row is above that
limit and uses the rendezvous protocol, where the data is transferred only
once the receiving process is inside an MPI call, which here is
`MPI_Waitall` at the end of the step. This Open MPI build has no progress
thread, so no data is transferred while the threads compute, and no layout
with two ranks on one node can overlap the exchange with the update. That
explanation fits every small grid row. It does not explain the big grid's
1.01x on the same layout, and it has not been tested by changing the eager
limits.

## M2.7 Denormals (jobs 651284 and 659436)

| Rung | not flushed | flushed | gain |
|---|---|---|---|
| baseline | 53 | 57 | 1.07x |
| serial | 124 | 159 | 1.28x |
| simd | 320 | 454 | 1.42x |
| omp, 24 threads | 4865 | 6484 | 1.33x |

Short benchmark runs leave many cells holding numbers below 2.2e-308, which
this CPU handles through a slow microcode path. The number of such cells
depends on the grid and `nt`, so without flushing, rungs and thread counts
run with different `nt` were being compared in different regimes. Every
price in the table is identical to ten digits with and without flushing.

## M2.8 The 2^24 byte address effect (jobs 651439, 651445, 651450)

simd on one thread runs at 169 to 183 on any grid whose two buffers start
exactly 2^24 bytes apart (2896 x 724, 1448 x 1448, 2044 x 1026), and at about
405 on neighbouring grids. omp on 2 threads shows the same (336 against
814). Serial ran normally at 2896 x 724 in job 651299 (167). On the current
code it runs at 58 there (job 656959), so the serial rung is slow there as
well. omp on 24
threads ran at 6499 there in jobs 651299 and 656959. In jobs 659241 and
659437, run on the same day with the split turned off and on, it ran at 3288
and 3267, so the 24-thread rung is now slow there as well, whichever build
runs. The trigger is established, and the hardware mechanism is still
unknown (CHANGELOG Phase 5, Finding 1).

## M2.9 Hybrid MPI and OpenMP (job 659242, figures F6 and F8)

| Layout | 2048 x 512 | exchange share | 4096 x 1024 | exchange share |
|---|---|---|---|---|
| 1 node, 1 rank x 24 threads | 6233 | 0 | 2907 | 0 |
| 1 node, 2 ranks x 12 (one per socket) | 6843 | 20% | 6624 | 17% |
| 2 nodes, 1 rank x 24 each | 7064 | 28% | 9002 | 30% |
| 2 nodes, 2 ranks x 12 each | 8310 | 40% | 10863 | 28% |

One rank per socket is faster than one rank of 24 threads by 10% on the
small grid and by 2.28 times on the big one. Each rank allocates its own
copy of the grid and its own threads write it first, so Linux places those
pages in the memory of the rank's own socket (PLAN §5 hypothesis 5, held).
Across two nodes the best layout gains 1.21 times at 2048 x 512 (61%
efficiency) and 1.64 times at 4096 x 1024 (82%). Measured against the
unrounded ladder baseline of 56.7 Mcell/s (job 659238, M2.1), the 8310 of
2 ranks x 12 on each of two nodes is 147 times the baseline, and the 6984
of the one node ladder rung is 123 times. The edge row exchange is 40%
of a 0.126 ms step on the small grid, about 50 microseconds. That is far
more than 64 KiB needs at any network speed, so per-message latency and
waiting for the slower neighbour dominate it. PLAN §5 hypothesis 6 (above
80% efficiency) did not hold on the small grid and held on the big one. The
run excludes Open MPI's fast transports (`OMPI_MCA_mtl`, `OMPI_MCA_btl`), as
the matrixass MPI runs needed, so messages go over TCP.

The two node layouts on the small grid vary more between jobs than the reps
of one job do. The `mpi` column of E9 (job 659243), run minutes later with
the same build, gave 9364 for 2 ranks x 12 on each node and 5693 for 1 rank
x 24 on each node, against 8310 and 7064 here. With the E9 numbers the best
two node gain on the small grid is 1.33 times (67%). Every other layout
agrees between the two jobs within 5%.

## M2.10 Hypotheses from PLAN §5 against the data

| # | Hypothesis | Outcome |
|---|---|---|
| 1 | SIMD gives 2 to 3.5x over serial | held, 2.87x |
| 2 | OpenMP scales close to linearly to 12 threads | mostly held, 9.6x at 12 (80%) |
| 3 | the large grid becomes memory-bound | held, 2990 at 4096 x 1024 against 6278, on the DRAM roof |
| 4 | the serial boundary is about 0.5% and does not limit scaling | overturned, it grows to 9.1% at 24 threads |
| 5 | one rank per socket beats one rank of 24 threads | held, 1.10x and 2.28x |
| 6 | two nodes above 80% efficient | overturned at 2048 x 512 (61%), held at 4096 x 1024 (82%) |

E2 (aligned loads), E3 (float32) and E9 (overlapping the edge row exchange
with the update) ran later and are reported in M2.6.

## M2.11 Profile and optimisation levels (jobs 652696 and 652697, figures F11 and F12)

This section was added after the Milestone 1 feedback asked for a flat
profile of the baseline and a four-bar -O0 to -O3 chart. Both were measured
on one thread at 2048 x 512.

The gprof flat profile was built with `-O2 -pg` (job 652697).

| Rung | Hot function | Share | Everything else |
|---|---|---|---|
| baseline | `BaselineSolver::solve` | 99.9% | 0.1%, filling the new buffers with zeros (`vector::assign`) |
| serial | `update_cells_scalar` | 99.4% | `update_zero_row` 0.3%, `copy_top_row` 0.15%, filling the buffers 0.15% |
| simd | `update_row_avx` | 98.9% | the scalar remainder 0.4%, `update_zero_row` 0.4%, filling the buffers 0.4% |

The baseline is one function because the whole algorithm is one loop nest,
so the flat profile can only say the time is in that nest. The later rungs
split the work into functions, and the profile shows the interior update
holds 99% of the time and the boundaries under half a percent on one
thread. That is the serial fraction a parallel version starts from.

The optimisation levels (job 652696) are in Mcell-updates/s. Each pair was
built with its own level in place of `-O2`, and every level gave exactly
the `-O2` price, delta, gamma and vega.

| Level | baseline | serial |
|---|---|---|
| `-O0` | 10.7 | 44.8 |
| `-O1` | 57.6 | 156.1 |
| `-O2` | 57.7 | 158.4 |
| `-O3` | 112.5 | 238.6 |

`-O0` keeps every variable in memory, so the baseline is 5.4 times slower
than `-O1`. `-O1` and `-O2` are the same for this loop. `-O3` turns on
GCC 8's auto-vectoriser, which uses two-lane SSE2 because no `-mavx2` is
given, and gains 1.95 times on the baseline and 1.51 times on the serial
rung. The hand-written AVX2 rung, built at `-O2`, reaches 456 (job 659238),
which is 1.91 times the `-O3` serial rung. Experiment E1 in M2.6 repeats these
`-O3` bars (112.5 and 241.3) and adds AVX2, `-march=native` and
`-ffast-math` builds.

