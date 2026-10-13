# Milestone 2 Build Plan. Parallel Heston PDE Option Pricer

Milestone 1 (serial) is finished and scored 20/24. Most of the lost marks were
for presentation, not the engineering. The Milestone 1 code, ladder and
numbers are frozen at the git tag `milestone1` (`git checkout milestone1`
reproduces every M1 figure). This file replaces the M1 plan and is the single
source of truth for Milestone 2.

Spec: `COSC3500_Milestone.pdf` (one rubric shared by both milestones, page 2).
Model background: `docs/InitialSpecDesign.md`. M1 findings: `docs/RESULTS.md`.

---

## 0. The one-paragraph mental model

We price one European option under Heston by filling a 2-D grid of what-if
scenarios (stock price by variance) backwards in time from expiry to today.
Every timestep, each interior cell is rewritten from the nine cells around it
in the previous sheet. Two buffers, swap per step. Because each step only
reads the old sheet and only writes the new one, every cell in a step is
independent of every other cell. That independence is the whole Milestone 2
story. It lets four neighbouring cells share one AVX instruction (SIMD), lets
24 cores each own a band of rows (OpenMP), and lets two machines each own
half the rows and swap one edge row per step (MPI).

---

## 1. Task sheet cross-check (what earns the marks)

The same rubric as Milestone 1, with three Milestone 2 additions marked M2.
Every row says what we deliver and which figure or artefact proves it.

| Criterion | Weight | Spec asks for | Our evidence |
|---|---|---|---|
| Introduction and Background | 20% | Self-contained. Assessors will not remember M1. | Retell the model in under 2:30 with the M1 weather-map and BS-collapse animations (`results/demo_timevalue.mp4`, `results/bs_collapse.mp4`). One slide on why the grid is parallel (the independence argument in §0). |
| Optimisation | 25% | Optimisations and their performance. Techniques that failed or were considered and rejected. How the problem was structured and parallelised (M2). Trade-offs. | The M2 ladder (F1), one technique per rung, each measured. The tried-and-rejected panel on F1 and the experiments in §6. The decomposition diagram (F10). Trade-offs written into CHANGELOG (FMA rounding, full grid per MPI rank, float32 accuracy). |
| Benchmarking | 25% | Setup that avoids misleading results. Amdahl's and Gustafson's laws (M2). Hotspots and bottlenecks. Scaling with problem size. A UQ cluster. | Protocol in §7 (exclusive node, pinned page regime, thread binding, median of 5, provenance header). Strong scaling with Amdahl fit (F2). Weak scaling with Gustafson (F3). Size scaling with cache bands (F4). Roofline v2 (F5). Time breakdown hotspot chart (F6). All on rangpur r730 nodes. |
| Presentation | 20% | Topic choice and time per topic, delivery, slide design, data visualisation. | §11 time budget, M1-style figures with provenance footers, code shots per rung. Fix list from M1 in §11. |
| Reflection and Conclusion | 10% | Main conclusions, what surprised us, what we would do differently, how to extend. | §11 reflection beat. Candidate surprises are listed in §6 as hypotheses so we can say which ones the data overturned. |
| Identity verification | pass/fail | Face visible throughout the video. | Webcam in every take (unchanged from M1). |
| Code submission | pass/fail | Well-commented source, Makefile, slurm scripts, enough to compile and check. | §3 comment style, `make` and `make test` from a fresh clone, every slurm script in the zip. |
| Submission format | pass/fail | H.264, under 100 MB, `Milestone2_48008361.mp4`, code and video in one `Milestone2_48008361.zip`. | `make package` updated from `Milestone1_` to `Milestone2_`, same secret checks. |
| In-person interview (M2) | pass/fail | About 7 minutes, identity-verified, second week of exams. Must pass to pass the course. | STUDY_GUIDE M2 parts, CHANGELOG, COMMANDS.md demo runbook, rehearsal in P10. |

Path choice. The spec accepts one technique analysed very thoroughly, a
comparison of technologies, or a combination such as OpenMP+MPI. It also says
using every technique is not required. We take the combination path (AVX2,
then OpenMP, then hybrid MPI+OpenMP) and analyse each rung with the same
rigour as the M1 ladder. CUDA is a bonus only after every other phase is done.

The spec PDF is password protected. This table was built from the course
rubric as already extracted for M1 (`presentation/PRESENTATION.md` §1) and the
knowledge graph of the PDF. Before scripting the video, open the PDF once and
check two things it may state that this table does not, the M2 due date and
any M2-only wording in the FAQ section.

---

## 2. The M2 ladder (locked)

The M1 lesson was that markers liked a ladder where each rung adds one thing
and is measured. M2 keeps that shape but starts from a short ladder of the
techniques that actually mattered. The video retells M1 in two rungs, then
spends its time on the parallel rungs.

| Rung | `--solver` | Adds exactly this | File |
|---|---|---|---|
| 0 | `baseline` | The naive reference. Weights recomputed per cell, divisions in the hot loop. Every other rung must match its answer. | `src/solver_baseline.cpp` (kept, restyled) |
| 1 | `serial` | The M1 optimisations that moved the number. Per-row constants and S, S² tables hoisted out of the time loop, reciprocals instead of divisions, row pointers. | `src/solver_serial.cpp` |
| 2 | `simd` | AVX2. One `__m256d` holds 4 neighbouring cells, the stencil blend uses FMA, a scalar loop finishes the remainder. | `src/solver_simd.cpp` |
| 3 | `omp` | OpenMP. `#pragma omp parallel for schedule(static)` over variance rows. Each thread owns a band of rows. | `src/solver_omp.cpp` |
| 4 | `mpi` | Hybrid MPI+OpenMP. Rows split into one slab per rank, one halo row exchanged with each neighbour per step, OpenMP threads inside each rank. | `src/solver_mpi.cpp` |

Why rung 1 drops M1 levels 3 to 6. In M1 only hoisting and strength reduction
(levels 1 and 2) changed the throughput (53 to 131 Mcell-updates/s, +149%).
Levels 3 to 6 were inside run-to-run spread, and unrolling by four made the
kernel 180 lines longer. The video can say this in one sentence while the M1
ladder figure is on screen, which also covers the failed-techniques criterion
for the serial half.

### 2a. How the code stays short

One shared file holds the maths, so each rung's file is only its time loop.

- `include/stencil.h`, `src/stencil.cpp`
  - `struct StencilConstants` holds the reciprocals, the per-row factors and
    the S and S² tables. Built once before the time loop. This also removes
    the M1 per-timestep allocation of the tables.
  - `update_row_scalar(k, cur, next, var_j)` updates one variance row (rung 1 maths).
  - `update_row_avx(k, cur, next, var_j)` updates one row, 4 cells per
    instruction, then the remainder (rung 2 maths).
  - `update_zero_row`, `apply_stock_edges(row_begin, row_end, tau)`,
    `copy_top_row` for the boundaries.
- Each solver's `solve()` is the time loop plus a row loop that calls one of
  those. The step from `simd` to `omp` is visibly one pragma line. The step
  from `omp` to `mpi` is a row range plus the halo exchange. Those two diffs
  are the code shots for the video.
- A shared protected `Solver::report_stability()` replaces the copy of the
  stability and Feller printout in every solver.

### 2b. Rung details and the decisions behind them

Each decision below becomes one `### Decision N.` entry in CHANGELOG.md and
one interview answer in STUDY_GUIDE.md.

**Rung 2, SIMD.**
- Intrinsics from `<immintrin.h>`, not compiler auto-vectorisation and not
  assembly (L04 says auto-vectorisation gives subpar results and assembly is
  not worth it).
- Unaligned loads `_mm256_loadu_pd`. The nine-point stencil reads i-1, i and
  i+1, so at most one of those three loads can ever be 32-byte aligned. L04
  also says the unaligned penalty is gone on modern CPUs. Aligned buffers are
  kept as a measured experiment (E2).
- FMA via `_mm256_fmadd_pd`. FMA rounds once instead of twice, so results move
  in the last bits. The tolerance against baseline is re-measured and recorded.
- Build with `-mavx2 -mfma`, not `-march=native`, so the binary needs only
  Haswell features. `objdump -d heston | grep zmm` must be empty, because L04
  warns that some AVX-512 intrinsics also start with `_mm256_`.
- On the arm64 Mac the AVX row falls back to the scalar row behind one
  `#ifdef __AVX2__`, with a stderr warning. Local AVX correctness is checked
  with an x86_64 build under Rosetta (verified working on this Mac), and every
  rangpur job runs the full test suite with the real AVX path.

**Rung 3, OpenMP.**
- `parallel for` on the outer variance loop (L05, parallelise outermost loops
  first). `schedule(static)` because every row costs the same (T04).
- `default(none)` with an explicit `shared(...)` list (L05 recommendation).
  GCC 8.5 on rangpur and newer compilers disagree on whether `const` variables
  may appear in `shared()`, so the shared variables are kept non-const and the
  code is compiled on both before it counts.
- No synchronisation inside a step. Threads read the shared old sheet and
  write disjoint rows of the new one. The only barrier is the implicit one at
  the end of the `parallel for`. Never `critical` in the hot loop (L05 demo
  measured a 40x slowdown).
- The boundaries are O(ns + nv) work and stay serial at first. That serial
  part is the Amdahl term we measure in F6 and fit in F2.
- Thread count comes from `OMP_NUM_THREADS`. The CSV records the count the
  program actually got from `omp_get_max_threads()`, because L05 warns you may
  not get what you asked for.

**Rung 4, hybrid MPI+OpenMP.**
- `MPI_Init`, `MPI_Comm_rank`, `MPI_Comm_size`. Rows split evenly with the
  remainder given to the first ranks (T05).
- Each step, before the update, each rank swaps its edge rows with its
  neighbours using `MPI_Irecv`, `MPI_Isend`, `MPI_Waitall` (L06 non-blocking
  pattern). Rank 0 also owns the v=0 row and the last rank owns the top
  Neumann row.
- Every rank keeps the full 2 x 8 MiB grid but only updates its own rows. This
  keeps indexing identical to the other rungs and lets one `MPI_Gatherv` at the
  end hand rank 0 a complete sheet for the existing `extract_result`. The
  trade-off is memory, which is tiny here. A large-grid version would allocate
  only the local slab plus two halo rows (T05 convention). Recorded as a
  trade-off, not hidden.
- One binary. `heston` is built with `mpicxx` and always calls `MPI_Init` and
  `MPI_Finalize`, so the non-MPI rungs run as one rank. Only rank 0 prints.
- Launch with `mpiexec`, one rank per socket or per node, threads inside
  (L06 hybrid slide). On rangpur Open MPI 4.1.1 is
  `mpiexec -n R --map-by ppr:1:socket:PE=12 --bind-to core`, plus the
  `OMPI_MCA_mtl` and `OMPI_MCA_btl` exclusions that the matrixass MPI script
  needed.
- Communication time is accumulated separately so F6 can show compute vs
  halo exchange per rank count.

### 2c. What changes in the shared code

- `src/main.cpp` gains `MPI_Init`/`MPI_Finalize`, rank 0 printing, and drops
  `--opt-level`.
- CSV line gains two trailing columns, `threads` and `ranks`.
- `make_solver` learns `serial`, `simd`, `omp`, `mpi`. `opt` is removed.
- Removed after the tag. `src/solver_opt.cpp`, `src/solver_opt_kernels.cpp`,
  `include/solver_opt_kernels.h`, `tests/test_opt_matches.cpp`, and the
  M1-only slurm sweeps (ladder, scaling, controls, novec, novec_o3, pages,
  flags, diag_alloc, diag_denormal). `docs/RUNNING.md` is folded into
  `COMMANDS.md`.
- Makefile. `CXX = mpicxx`. OpenMP flags picked by compiler (Apple clang needs
  `-Xpreprocessor -fopenmp -lomp` with Homebrew libomp, GCC needs `-fopenmp`).
  `-mavx2 -mfma` only when the machine is x86_64. New object files and header
  dependencies listed by hand as before. `make test-mpi` runs the ladder test
  under `mpiexec -n 2` and `-n 3`. `make avx-check` builds the AVX row test
  for x86_64 under Rosetta.

---

## 3. Code style for M2 (replaces M1 §1b)

The user will walk an interviewer through this code, so it has to read like a
script. The model is `../matrixass/matrixMultiply.cpp`.

- A `//` comment on its own line directly above the line it describes, same
  indent. One or two straightforward sentences or noun phrases (CLAUDE.md, "How
  the user writes"), each ending in a full stop. No first person, no "This
  is the X that" openers, no "not X but Y" contrasts. No padding clauses tacked on the
  end. A reason only when the reader needs that fact to follow the code. This replaced the earlier 3 to 12 word
  fragments on 2026-10-13.
- A block opens with one sentence saying what the whole block does, then the
  per-line comments.
- Every function, class and struct has a two-line plain English overview
  directly above it saying what it is for and why it exists.
- Beginner layout. One statement per line, long lines are fine, `if`/`else`
  rather than `?:`, fields set in the constructor body rather than an
  initialiser list, the PDE split into named terms.
- Pragmas get their comment above the pragma. Else branches get theirs inside
  the branch.
- No lecture citations, plan references, phase names or essays in code. Those
  live in STUDY_GUIDE.md and CHANGELOG.md.
- No semicolons inside comments, no em dashes.
- Descriptive names (`row_above`, `owned_rows`, `rows_per_rank`), never cryptic
  abbreviations. A name says what the thing is, not a metaphor for it
  (`buffer` over sheet, `owned_rows` over slab, `row_below` over halo). The
  textbook term goes once in a comment.
- Each file opens with the same three-line AI declaration as matrixass.
- The kept M1 files on the demo path are restyled to match (`main.cpp`,
  `grid.*`, `params.*`, `solver.*`, `solver_baseline.cpp`, `black_scholes.*`,
  `io.*`, tests, slurm scripts). Behaviour must not change, and the tests prove
  it.
- The C++ shape from M1 stays (classes with private fields, `Solver` abstract
  base with `override`, plain structs for records, RAII, `unique_ptr` factory).

CLAUDE.md's code style paragraph and its Milestone 1 scope section are updated
to say this in P0.

---

## 4. Correctness contract (the iron rule carries over)

No rung is benchmarked until it matches the baseline. Every benchmark job
builds, runs the tests, and refuses to time anything if a test fails.

| Gate | Checks | Tolerance |
|---|---|---|
| `test_bs_collapse` | baseline vs closed-form Black-Scholes at xi=0 | 2.5e-3 relative (unchanged) |
| `test_parity` | put-call parity on baseline | 1e-4 absolute (unchanged) |
| `test_ladder_matches` (new, replaces `test_opt_matches`) | serial, simd, omp, mpi vs baseline, call and put, price, delta, gamma, vega | relative, measured first then set with headroom. M1 measured 3.3e-15 for reciprocals. FMA and thread-split are re-measured. |
| `make test-mpi` | the same test under `mpiexec -n 2` and `-n 3` (3 forces an uneven row split) | same as above |
| OpenMP split check | tests at `OMP_NUM_THREADS=1, 2, 7` (7 does not divide the row count) | same as above |
| `make check-heston` | Fourier reference price on the full parameter set, local only | 2e-3 (unchanged) |
| `make avx-check` | AVX row vs scalar row on random data, x86_64 build | measured, about 1e-15 relative |

Note on reproducibility. With `schedule(static)` and a fixed thread count each
cell is computed by exactly one thread in the same order, so the OpenMP rung
should be bit-identical to the SIMD rung. MPI likewise. If they are not, that
is a bug, not floating-point reordering, because no reduction happens inside a
step. This is a stronger claim than L05's warning about non-associativity and
is worth saying in the interview.

---

## 5. Measured hypotheses going in

Stated now so the video can say which ones the data confirmed and which it
overturned (the Reflection criterion rewards surprises).

1. SIMD gives 2 to 3.5x over rung 1. M1 showed the kernel is compute-bound
   (82 to 95% of scalar peak at 1.29 flop per byte), so the vector width should
   mostly carry through. L04's rule of thumb is about 4x.
2. At 2048x512 the two buffers (16 MiB) fit inside one socket's 30 MiB L3, so
   OpenMP should scale close to linearly up to 12 threads and keep scaling
   across the second socket.
3. At 4096x1024 (64 MiB) the working set spills to DRAM, and the 24-thread
   AVX kernel becomes memory-bound. The size sweep (F4) should show a
   throughput cliff at the L3 capacity line, and the roofline (F5) should show
   the operating point move from the compute roof to the memory roof.
4. The serial boundary work is about 0.5% of a step, so Amdahl's ceiling is far
   above 24 cores. The measured efficiency loss will instead come from
   fork-join overhead per step, NUMA traffic across sockets, and memory
   bandwidth. F6 separates these.
5. Two ranks per node (one per socket) beat one rank with 24 threads, because
   each rank's memory stays on its own socket.
6. Across two nodes the halo is 2 rows of 2048 doubles (32 KiB) per step, small
   next to about 1 ms of compute per step, so 2-node strong scaling should be
   above 80% efficient at 2048x512 and better at larger grids.

---

## 6. Experiments (tried, measured, kept or rejected)

Each one gets a CHANGELOG entry, kept or rejected, with the number. Rejected
ones go on the F1 side panel, which plays the role the negative controls
played in M1.

| Id | Experiment | Rung | What we measure | Why it is interesting |
|---|---|---|---|---|
| E1 | `-O3` and `-ffast-math` vs `-O2` on the serial rung | 1 | throughput, answer drift | Answers whether the compiler could have done the SIMD for us. M1 showed `-O3` auto-vectorised the baseline but not level 6. |
| E2 | aligned buffers (`_mm256_load_pd` on a 32-byte aligned grid) vs unaligned | 2 | throughput | Tests the L04 claim that the alignment penalty is gone. |
| E3 | float32 with `__m256` (8 lanes) | 2 | throughput and error vs double | Doubles the lane count. Trade-off of speed for accuracy. If time allows. |
| E4 | one `parallel for` per step vs one persistent `parallel` region with `omp for` and `omp single` | 3 | throughput at 24 threads, small and large grids | Fork-join cost per step. L05's boundary-exchange pattern. |
| E5 | `schedule(static)` vs `schedule(dynamic, chunk)` | 3 | throughput | T04 says static suits equal work. Confirms it. |
| E6 | `OMP_PROC_BIND` none / close / spread with `OMP_PLACES=cores` | 3 | throughput at 6, 12, 24 threads | Two sockets. Placement decides which L3 and which memory a thread uses. |
| E7 | first-touch initialisation (each thread touches the rows it will update) | 3 | throughput at 4096x1024, 24 threads | NUMA page placement. Matters only when the grid spills L3. |
| E8 | parallel boundaries vs serial boundaries | 3 | throughput, serial fraction | Shrinks the Amdahl serial term. Probably too small to see, which is itself a result. |
| E9 | overlapping the halo exchange with interior rows (post Irecv/Isend, update interior, wait, update edge rows) | 4 | throughput on 2 nodes | L06 overlap pattern. |
| E10 | ranks per node 1 vs 2 (per socket) | 4 | throughput on 1 and 2 nodes | Hypothesis 5. |
| E11 | `MALLOC_MMAP_THRESHOLD_` page regime re-checked at 24 threads | all | throughput | M1 found page size worth 1.66x serially. Checks whether it still matters when threaded. |

---

## 7. Benchmark protocol (Benchmarking 25%)

Carried over from M1 and extended for threads and ranks.

- Every number comes from a rangpur compute node through `sbatch`. Never the
  login node. The plotting scripts refuse any CSV whose `host=` header is not a
  `*.compute.eait.uq.edu.au` name (kept from M1).
- Single-node runs pin `--nodelist=r730-2 --exclusive`, the same machine as M1
  (dual-socket Xeon E5-2670 v3, 2 x 12 cores, 30 MiB L3 per socket, AVX2+FMA,
  GCC 8.5). Two-node runs use `r730-2,r730-3` exclusive. As of 2026-10-11
  `sinfo` shows r730-0 mixed, r730-1 down, r730-2 and r730-3 usable.
- `--cpus-per-task` set to the thread count, `--ntasks` and
  `--ntasks-per-node` to the rank layout. QOS caps wall time at 15 minutes, so
  each sweep is its own job sized under 14 minutes.
- Page regime pinned with `MALLOC_MMAP_THRESHOLD_` exactly as the M1 headline,
  passed to every rank with `mpiexec -x`. Recorded in the header.
- Thread binding recorded in the header (`OMP_PROC_BIND`, `OMP_PLACES`), and
  `--report-bindings` output kept in the `.err` file for MPI runs.
- Median of 5 repetitions, min to max range reported, as in M1.
- Each repetition runs for at least 1 second. Fast rungs get more timesteps,
  and the maturity is scaled with nt so dt stays at the same fraction of the
  stability bound. The rate (Mcell-updates/s) is unaffected.
- Header records host, CPU, caches, compiler, flags, threads, ranks, binding,
  page regime, job id, date.
- Weak scaling grows `nv` with the thread count and scales `v_max` with it so
  the spacing dv, and therefore dt, stays the same. Every step costs the same
  per cell.
- Peak and bandwidth ceilings are measured on the same node, not taken from a
  datasheet (scalar FMA peak from M1, new AVX2 FMA peak per core and for 24
  cores, triad bandwidth for 1 and 24 threads).

### 7a. Sweeps (one slurm script each)

| Script | What it runs | Feeds |
|---|---|---|
| `slurm/debug.sh` | build and `make test` on a compute node, 5 minutes, not exclusive | gate before any sweep |
| `slurm/bench_m2_ladder.sh` | baseline, serial, simd, omp(24), mpi(1 node, 2 ranks x 12) at 2048x512 | F1 |
| `slurm/bench_m2_strong.sh` | omp at 1, 2, 4, 6, 8, 10, 12, 16, 20, 24 threads, fixed 2048x512 | F2 |
| `slurm/bench_m2_weak.sh` | omp at the same thread counts, rows per thread fixed | F3 |
| `slurm/bench_m2_size.sh` | serial, simd, omp(24) at 512x128 up to 4096x1024 and one grid past both L3s | F4 |
| `slurm/bench_peak.sh` (extended) | scalar FMA peak, AVX2 FMA peak (1 core and 24), triad bandwidth (1 and 24 threads) | F5 |
| `slurm/bench_m2_experiments.sh` | E1, E2, E4, E5, E6, E7, E8, E11 | F1 side panel, F7 |
| `slurm/bench_m2_mpi.sh` | 1 node x {1, 2} ranks, 2 nodes x {1, 2} ranks per node, threads inside, E9 | F8, F6 |
| `slurm/bench_m2_float.sh` | E3 | F1 side panel, F9 |

---

## 8. Figures (Benchmarking and Presentation, 45% together)

Every figure follows the M1 house style from `scripts/figstyle.py`. Light
background, one highlight colour for real gains, grey for nulls, red for the
rejected or control panel, value and speedup printed on each bar, and a
provenance footer (host, CPU, compiler, flags, threads or ranks, date) along
the bottom. Each figure exists to highlight one HPC change and its impact.

| Id | Figure | Built from | The HPC change it highlights | The sentence it lets us say |
|---|---|---|---|---|
| F1 | **The M2 ladder.** Bars for baseline, serial, simd, omp(24), mpi, labelled with Mcell-updates/s and cumulative speedup, per-rung gain printed inside the bar. Side panel of tried-and-rejected variants (E2, E3, E4, E5, the losing binding), styled like the M1 negative-controls panel. Log y-axis if the range demands it, with the speedup labels carrying the story. | ladder, experiments | every rung at once | Each technique on top of the last, and what it was worth. |
| F2 | **Strong scaling.** Left, speedup vs threads with the ideal line and the fitted Amdahl curve, serial fraction f printed. Right, parallel efficiency vs threads. A vertical line at 12 threads marks the socket boundary. Both simd-based omp and a scalar-based omp line so SIMD and threads are shown to compound. | strong | OpenMP | Speedup flattens where the fit says, and here is the serial fraction behind it. |
| F3 | **Weak scaling.** Time per step vs threads with work per thread fixed, and scaled speedup against the Gustafson line. | weak | OpenMP | With more cores we can solve a proportionally bigger grid in the same time. |
| F4 | **Throughput against grid size.** serial, simd and omp(24) lines across grid sizes, with shaded bands for "fits in L2", "fits in one L3", "fits in both L3s", "DRAM", the same way M1 shaded the page regime. Right panel, speedup of omp(24) over serial per size. | size | memory hierarchy vs parallelism | Small grids lose to overhead, large grids hit memory bandwidth. |
| F5 | **Roofline v2.** M1 roofline extended with the AVX2 per-core peak, the 24-core AVX2 peak, and measured 1-thread and 24-thread bandwidth roofs. Points for baseline, serial, simd, omp(24), mpi, coloured by rung. | peak, ladder, size | SIMD lifts the compute roof, threads move the kernel toward the memory roof | Shows whether hypothesis 3 held, the kernel changed from compute-bound to memory-bound. |
| F6 | **Where the time goes.** Stacked bars of time per step split into stencil update, boundaries (the serial part), halo exchange, and idle wait, for omp at 1, 12, 24 threads and mpi at 1, 2, 4 ranks. | strong, mpi | hotspot and bottleneck identification | The serial part is X%, the halo exchange is Y%, here is what limits scaling. |
| F7 | **Thread placement and NUMA.** Grouped bars for binding none, close, spread, and first-touch on and off, at 12 and 24 threads, at 2048x512 and 4096x1024. | experiments | OpenMP on a two-socket node | Where threads and pages live matters as much as how many threads there are. |
| F8 | **Hybrid layouts.** Bars for 1 node x 1 rank x 24 threads, 1 node x 2 ranks x 12, 2 nodes x 1 x 24, 2 nodes x 2 x 12, with the communication fraction printed. Second panel, 2-node strong-scaling efficiency at two grid sizes. | mpi | MPI, then hybrid | Splitting by socket, then by node, and what the network costs. |
| F9 | **Every rung gives the same answer.** Max relative difference vs baseline per rung on a log axis with the test tolerance as a line. float32 shown if E3 runs, far above the line. | test output | correctness of each parallel change | Faster and still right, and the one variant that was faster but not right enough. |
| F11 | **What the compiler gives.** baseline and serial built at -O0, -O1, -O2, -O3, one thread. Added from the M1 feedback. | optlevels | the compiler before our changes | -O3's own vectoriser gives 1.5x on serial, the hand AVX2 rung 2.4x. |
| F12 | **Where the time goes.** gprof flat profile of baseline, serial and simd on one thread. Added from the M1 feedback. | profile | hotspot identification | 99% of the time is the interior update, the boundaries under half a percent. |
| F10 | **Who updates what.** Static diagram of the grid with rows coloured by thread (OpenMP), then slabs by rank with halo rows highlighted (MPI), and four cells in one AVX register. Optional animation in the `memory_anim.py` style. | none (diagram) | how the problem was structured and parallelised | Optimisation criterion asks for the structure explicitly. |

Kept from M1 for the Introduction: weather map, BS collapse, payoff diagram.
The M1 ladder figure (`results/bench_ladder_hugepage.png`) appears for one
slide as the serial backstory.

Code shots (`scripts/code_shots.py`, as in M1). One per rung: the scalar row
update, the four FMA lines of the AVX row, the one-line pragma, and the halo
exchange. Plus the CSV provenance header.

Scripts. `scripts/bench_plot.py` is extended for F1 to F4 and F6 to F9 and
keeps the compute-node refusal. `scripts/roofline.py` gains the extra roofs
for F5. `scripts/memory_anim.py` gains a decomposition style for F10.

---

## 9. Documents

| File | What it holds |
|---|---|
| `PLAN.md` | This plan. The phase checklist in §10 is ticked as work lands. |
| `CHANGELOG.md` (new) | Matrixass format. `## Phase N. Title`, then `### Decision N.`, `### Change N.` or `### Tried and rejected N.`, then `### Validation` with job ids and numbers, then `### Where this leaves the project`. Every code change appears here. |
| `COMMANDS.md` (new, replaces `docs/RUNNING.md`) | Every command, local and rangpur, with what each argument and environment variable does. Build, test, run each rung, `OMP_NUM_THREADS`, `OMP_PROC_BIND`, `OMP_PLACES`, `mpiexec` flags, sync, sbatch, squeue, fetch, plots. Ends with a five-minute live demo runbook for the interview. |
| `STUDY_GUIDE.md` | The maths and ideas. Updated with every code change, not only in P8 (see §10). M1 parts kept. New parts for SIMD vs MIMD, AVX2 and FMA, OpenMP, MPI, hybrid, NUMA, Amdahl and Gustafson, strong vs weak scaling, roofline. Interview practice questions with answers. Glossary. Links to the Intel Intrinsics Guide, the OpenMP 4.5 specification (what GCC 8.5 implements) and the Open MPI 4.1 `mpirun` manual. |
| `CODE_GUIDE.md` | The code, file by file in run order, each snippet copied from the real file with the comments stripped and two plain sentences under it, what it does and why (matrixass format). Each new rung file gets its own section. Practice questions on the code. |
| `docs/RESULTS.md` | New M2 section, every number with its job id, same discipline as M1. |
| `HPC_JUSTIFICATIONS.md` | Why each HPC technique was used, as interview answers in full sentences, each with the number that backs it. Covers the problem structure, memory layout, every rung, the measurement method, the scaling models and the rejected experiments. Updated when a decision changes. |
| `CLAUDE.md` | Milestone 2 scope and the new comment style. |

---

## 10. Phases

Each phase lists what it needs, the work, when it is done, and what validates
it (matrixass MASTERDOC format).

Every phase that changes code also updates the two guides in the same
session, alongside its CHANGELOG entry. New or changed code gets a block in
CODE_GUIDE.md, copied from the real file, followed by one sentence on what
the block does and one on why we do it. New ideas behind the code (SIMD,
OpenMP, MPI, scaling) go in STUDY_GUIDE.md. A diagram is added or updated
wherever a picture explains it faster than text (data layout, which worker
owns which rows, halo exchange, timing breakdowns). Snippets that no longer
match the code are fixed or removed, never left stale.

**P0. Setup** (done this session unless marked)
- [x] Tag `milestone1` on the M1 commit.
- [x] Open MPI 5.0 installed locally with Homebrew. Rangpur has Open MPI 4.1.1 (`module load mpi/openmpi-x86_64`, `mpicxx` wraps g++ 8.5).
- [x] Checked an x86_64 AVX2+FMA build runs on this Mac under Rosetta.
- [x] Checked r730 availability with `sinfo`.
- [x] Create CHANGELOG.md and COMMANDS.md, update CLAUDE.md and memory with the comment style.
- Done when the docs skeletons exist and `make` still builds the M1 tree.

**P1. Shared stencil, `serial` rung, restyle** (done 2026-10-11, see CHANGELOG Phase 1)
- Needs P0.
- `stencil.h/.cpp`, `SerialSolver`, `report_stability`, factory entry, `test_ladder_matches`, Makefile with `mpicxx` and OpenMP flags, restyle the kept files, remove the M1 opt ladder.
- Done when `make test` passes locally and serial matches baseline.
- Validated by `test_ladder_matches` and `make check-heston`.

**P2. `simd` rung** (done 2026-10-11, see CHANGELOG Phase 2)
- `update_row_avx`, `SimdSolver`, `make avx-check`.
- Done when the Rosetta AVX check and a rangpur `debug.sh` job both pass and `objdump` shows no zmm registers.

**P3. `omp` rung** (rung done 2026-10-11, see CHANGELOG Phase 3. E4 to E8 run with the P5 sweeps)
- `OmpSolver`. Then E4 to E8 and keep the winner.
- Done when tests pass at 1, 2 and 7 threads locally and on rangpur.

**P4. `mpi` rung** (rung done 2026-10-11, see CHANGELOG Phase 4. E9 runs with the P5 sweeps)
- `MpiSolver`, halo exchange, Gatherv, comm timer, `make test-mpi`. Then E9.
- Done when tests pass at `-n 1, 2, 3` locally and on rangpur across two nodes.

**P5. Benchmarks on rangpur** (done 2026-10-11, see CHANGELOG Phase 5. Every sweep re-run with denormals flushed, jobs 651294 to 651303, 651474. E2 and E9 not run)
- Rewrite `bench_common.sh`, write the §7a sweeps, extend `bench_peak.sh`.
- Done when every sweep CSV is fetched with a compute-node header.

**P6. Figures** (done 2026-10-12. F1 to F10 drawn by `scripts/m2_figures.py` and `scripts/decomposition_diagram.py` into `results/m2/`, code shots by `scripts/code_shots.py` into `results/m2/code_shots/`, see CHANGELOG Phase 6)
- F1 to F10 and the code shots.
- Done when every figure renders from fetched CSVs and its numbers match RESULTS.md.

**P7. float32 (E3), only if P1 to P6 are done** (done 2026-10-12, job 652647, see CHANGELOG Phase 7. Rejected for accuracy, 1.75x to 4.46x faster)

**P8. Documents** (done 2026-10-12, see CHANGELOG Phase 8. COMMANDS.md tested, demo runbook in §5, three missing STUDY_GUIDE answers added)
- STUDY_GUIDE M2 parts, CHANGELOG complete, COMMANDS.md tested command by command, RESULTS.md M2 section.
- Done when every CHANGELOG decision has a matching STUDY_GUIDE answer.

**P9. Presentation** (plan, script and deck done 2026-10-12, see CHANGELOG Phase 9 and 9b, `presentation/M2_SCRIPT.md` and `presentation/deck_m2.pdf`. Recording not done)
- New `presentation/` plan for M2 using §11. Video recorded in a later session.

**P10. Interview rehearsal**
- Final pass of STUDY_GUIDE.md and CODE_GUIDE.md before the interview. Check every snippet against the current code, check every snippet has its what and why sentences, check every diagram matches the final ladder, bring the older Milestone 1 sections in line with the writing style (no em dashes, no big adjectives), and trim history that does not help answer a question.
- Run the COMMANDS.md demo on rangpur end to end. Answer every STUDY_GUIDE practice question out loud.
- Done when each guide can be read start to finish in one sitting and every practice question has an answer in it.

**P11. CUDA bonus, only if everything above is done**

---

## 11. Presentation plan (Presentation 20%, and fixing M1)

### 11a. Milestone 1 feedback (82.6/100, one marker)

| Criterion | Mark | Comment | Change for M2 |
|---|---|---|---|
| Introduction | 15.9 | Excellent domain framing, maths and background. Two-buffer grid design but no data structures. Nearly half the video passed before any HPC content. Add the `index()` mapping, row-major stride, bytes per grid and the loop nest. | Finance and maths cut to one minute. A data structure beat at 1:20 shows `index()`, the row stride, 8 MiB per sheet at 2048 x 512 and the loop nest (code shot `shot_data_layout.png`). HPC content starts by 1:20. |
| Optimisation | 22.9 | Great methodology, nice strength-reduction insight. Most changes turned out null but handled well, push insights more. | Every rung beat states its insight as one sentence with its number, for example why SIMD gives 2.4x and not 4x, why `close` alternates sockets, why float is 4.5x faster on a big grid and still rejected. |
| Benchmarking | 18.9 | Full hardware and compiler details. Add a gprof or perf flat profile of the baseline and a four-bar -O0/-O1/-O2/-O3 chart for baseline. | `slurm/profile_m2.sh` (gprof flat profiles of baseline, serial and simd) and `slurm/bench_m2_optlevels.sh` (baseline and serial at -O0 to -O3), drawn as F11 and F12, shown in a where the time goes beat before the first rung. |
| Presentation | 15.9 | Professional slides with graphics and clean code diffs. A streaming overlay partly hid slide headers. Delivery fast and occasionally slurred on named theorems, use a cursor. | No overlay over the slide, webcam only in a reserved bottom-right box that no slide content uses. Script budgeted at 140 words per minute instead of 150. A pronunciation list for named results. A visible cursor or pointer on every figure, pointing at the number being said. |
| Conclusion | 9.0 | Excellent structured reflection, main insight and the page-size effect. Show speedup figures and the SIMD ceiling to quantify the target. | The conclusion slide shows the ladder numbers and the measured ceilings (AVX2 peak, DRAM bandwidth, the roofline position of each rung), then the hypotheses that held and the two that did not. |

### 11b. Rules carried over

- Target 9 minutes 20, ceiling 10. Word count per take checked at 140 words
  per minute.
- Self-contained introduction that never says "as in Milestone 1".
- One idea per slide. The figure is the slide, and the spoken line says the
  one sentence from the §8 table.
- Fewer, larger numbers on screen. The full tables live in RESULTS.md.

### 11c. Beats

| Beat | Time | Rubric criterion | Assets |
|---|---|---|---|
| Cold open, the final speedup | 0:20 | Presentation | F1 final bar |
| The problem, option, Heston PDE, grid | 1:00 | Introduction | payoff diagram, weather map |
| The data structure, `index()`, stride, bytes, loop nest | 0:40 | Introduction | `shot_data_layout.png`, F10 |
| Where the time goes, profile and -O levels | 0:40 | Benchmarking | F11, F12 |
| Serial rung, strength reduction | 0:30 | Optimisation | `shot_serial_row.png` |
| SIMD, and why not 4x | 1:00 | Optimisation | `shot_simd_fma.png`, F5 |
| OpenMP, Amdahl and placement | 1:15 | Optimisation, Benchmarking | `shot_omp_pragma.png`, F2, F7 |
| Hybrid MPI | 0:50 | Optimisation | `shot_mpi_exchange.png`, F8 |
| Same answer on every rung | 0:25 | Benchmarking | F9 |
| Weak scaling and grid size | 0:45 | Benchmarking | F3, F4 |
| Bottlenecks | 0:35 | Benchmarking | F6, F5 |
| Tried and rejected | 0:30 | Optimisation | F1 side panel |
| Conclusion, numbers and ceilings | 0:40 | Conclusion | conclusion slide |
| Close | 0:10 | | |

## 12. Interview preparation

The interview is about 7 minutes and pass/fail. The user should be able to
open any file, explain any line, and run the code live.

- STUDY_GUIDE practice questions cover at least these. Why rows and not
  columns for threads. Why no locks are needed. What `default(none)` does and
  why. What happens if `OMP_NUM_THREADS` exceeds the cores. Why unaligned
  loads. What FMA changes. Why the halo is one row. Why `Irecv` before
  `Isend`. Why every rank holds the full grid and what you would change for a
  huge grid. What Amdahl's fit says about this code. Strong vs weak scaling.
  Why the kernel became memory-bound. What the roofline shows.
- `HPC_JUSTIFICATIONS.md` holds the why answers for every HPC choice, ready to say out loud.
- COMMANDS.md demo runbook. Sync, submit `debug.sh`, watch `squeue`, show the
  CSV, run one rung interactively on a small grid with `salloc` if allowed, or
  show a saved job output.

---

## 13. Risks and cut order

| Risk | Mitigation |
|---|---|
| r730-3 busy or down, so no 2-node runs | 1-node 2-rank MPI still demonstrates decomposition and halo exchange. Two-node runs submitted early and retried. |
| GCC 8.5 rejects the `default(none)` shared list | Non-const shared variables, compile on rangpur in P3 before benchmarking. |
| AVX path only testable on rangpur | Rosetta x86_64 build locally, `debug.sh` on rangpur for every AVX change. |
| Queue congestion near the deadline | Sweeps written in P5 but `debug.sh` submitted from P2 onwards. |
| Scope creep | Cut order is CUDA first, then float32, then the F10 animation (keep the static diagram), then E9 overlap. Never cut correctness tests, F1, F2, F3 or F6. |

---

## 14. Definition of done (Milestone 2)

- [ ] `make clean && make && make test` passes locally and on a rangpur compute node, and `make test-mpi` passes at 1, 2 and 3 ranks.
- [ ] Every rung matches baseline within the measured tolerance (F9).
- [ ] All §7a sweeps run on r730 compute nodes with provenance headers.
- [x] F1 to F9 rendered, F10 diagram drawn, code shots made.
- [ ] Amdahl serial fraction fitted and Gustafson weak-scaling shown.
- [ ] Every experiment in §6 that ran has a kept-or-rejected CHANGELOG entry.
- [ ] STUDY_GUIDE M2 parts, COMMANDS.md and RESULTS.md complete.
- [ ] STUDY_GUIDE and CODE_GUIDE final pass done (P10), every snippet matches the code.
- [ ] 10-minute H.264 video under 100 MB, `Milestone2_48008361.zip` built by `make package`.
- [ ] Demo runbook rehearsed on rangpur.
