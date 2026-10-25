# Changelog

This file records every Milestone 2 change, in order, with what was done and
why. It is written to be read top to bottom so each decision can be explained
in the interview. The Milestone 1 code and numbers are frozen at the git tag
`milestone1`.

## Phase 0. Setup

### What exists and why

The Milestone 1 commit is tagged `milestone1`, so `git checkout milestone1`
reproduces every M1 figure. Open MPI 5.0 is installed locally with Homebrew.
Rangpur has Open MPI 4.1.1, loaded with `module add mpi/openmpi-x86_64`, and
its `mpicxx` wraps g++ 8.5. An x86_64 build with `-mavx2 -mfma` runs on this
Mac under Rosetta, which is how the AVX code will be checked locally before it
goes to rangpur. `sinfo` on 2026-10-11 showed r730-2 and r730-3 usable.

This file and `COMMANDS.md` were created, and `CLAUDE.md` now describes the
Milestone 2 scope and comment style.

## Phase 1. Shared stencil, serial rung, restyle

### What changed

The seven-level M1 optimisation ladder (`solver_opt.cpp`,
`solver_opt_kernels.cpp`, `solver_opt_kernels.h`, about 1100 lines) is gone.
In its place are two small files, `include/stencil.h` and `src/stencil.cpp`,
that hold the maths every rung shares, and `src/solver_serial.cpp`, which is
only the time loop. The `--solver` flag now takes `baseline` or `serial`, and
`--opt-level` no longer exists. `tests/test_ladder_matches.cpp` replaces
`tests/test_opt_matches.cpp`. The M1-only slurm sweeps were removed and
`debug_serial.sh` became `debug.sh`. Every kept source file was restyled to
the new comment style.

### Decision 1. Rung 1 keeps only the M1 techniques that moved the number

In M1 only hoisting (level 1) and replacing divisions with reciprocals
(level 2) changed the throughput, from 53 to 131 Mcell-updates/s. Levels 3 to
6 were inside run-to-run spread, and unrolling by four made the kernel 180
lines longer. So the serial rung has the hoisted constants, the S and S²
tables, the reciprocals and row pointers, and nothing else.

### Decision 2. One shared stencil file, one row at a time

`StencilConstants` holds every value that does not change during a solve.
That means the six reciprocals, the S and S² tables, and the four PDE
coefficients for each variance row (`half_v`, `rho_xi_v`, `half_xi2_v`,
`mean_reversion`). It is built once before the time loop. In M1 the S tables
were rebuilt every timestep, so this also removes one allocation per step.

`update_row_scalar(k, cur, next, var_j)` updates one interior row. Each later
rung only changes how rows are handed out. The SIMD rung adds an AVX version
of the row function, OpenMP puts a pragma on the row loop, and MPI gives each
rank a range of rows. That keeps each solver file to its time loop and makes
the step from one rung to the next a small diff.

The boundaries are split the same way. `update_zero_row` does the v = 0
transport row, `apply_stock_edges(k, next, row_begin, row_end, tau)` writes
the known S = 0 and S = s_max values for a range of rows, and `copy_top_row`
does the zero-slope top row. The row range on the stock edges is there
because an MPI rank will only own some rows.

### Decision 3. Copy the constants into locals inside the row function

`update_row_scalar` reads its constants out of the `StencilConstants`
struct, then writes through a `double*` into the next sheet. The compiler
cannot prove that pointer does not point into the struct, so without care it
would reload every constant from memory after every store. Copying them into
local variables at the top of the function tells the compiler they cannot
change, so they stay in registers.

### Decision 4. Shared stability report

The baseline and the M1 opt solver each had a copy of the stability estimate
and Feller printout. It is now one protected method, `Solver::report_stability`,
which labels the line with the solver's `name()`. The CSV label and the
off-node warning also use the solver's own name now, which removed
`solver_label` from `io.cpp`.

### Decision 5. Build with mpicxx and OpenMP flags from the start

The Makefile now sets `CXX = mpicxx`, because the final program is one binary
for every rung. OpenMP flags are chosen by compiler. Apple clang needs
`-Xpreprocessor -fopenmp` plus Homebrew's libomp, and GCC only needs
`-fopenmp`. The OpenMP flags sit in their own `OMPFLAGS` variable so a
benchmark script that overrides `CXXFLAGS` does not drop them.

`-mavx2 -mfma` is not added yet. GCC contracts `a*b + c` into an FMA by
default in C++ mode when FMA is enabled, so adding `-mfma` now would quietly
change the baseline and serial answers and timings. That flag lands with the
SIMD rung, where the effect can be measured on its own.

### Validation

The baseline gives bit-identical price, delta, gamma and vega to the M1
baseline on the smoke and demo grids, call and put. The serial rung matches
the baseline to 3.6e-16 relative on the call and 2.2e-15 on the put
(`test_ladder_matches`, tolerance 1e-11). `test_bs_collapse` (1.16e-3 and
1.12e-3) and `test_parity` (gap 5.4e-6) are unchanged from M1.
`make check-heston` passes both legs, with a relative gap of 4.2e-4 against
Fourier inversion on the off-node grid.

On rangpur, debug job 650550 (r730-0) built the tree with `mpicxx` over g++
8.5 and `-fopenmp` with no warnings and passed `make test`. There the serial
rung differs from the baseline by 3.3e-15 on the call and exactly 0 on the
put, the same as the M1 level 2 kernel under the same compiler.

### Change 1. Names that read like the PDE

Every single-letter or shorthand name in the solvers was replaced. The
constants struct is `constants` rather than `k`, the grid is `grid` rather
than `g`, the read sheet is `current`, the loop counter is `step` and the
discount time is `time_to_expiry`. In both the baseline and the stencil the
updated cell is `centre`, its neighbours are spelled out (`north_east` rather
than `ne`), and the derivatives are `dV_dS`, `d2V_dS2`, `dV_dv`, `d2V_dv2` and
`d2V_dS_dv`. The per-row coefficients are named after the PDE term they
multiply (`stock_diffusion`, `cross_diffusion`, `variance_diffusion`,
`variance_drift`, `stock_drift`). `kappa_theta` was dropped because it equals
`variance_drift[0]` exactly.

No arithmetic was reordered. The baseline and serial outputs on the smoke and
demo grids are bit-identical to before the rename and to the M1 baseline.

`compile_flags.txt` was added so the editor's language server can find
`include/` and jump from a name to its definition.

### Change 2. Beginner-style layout

Every source file was rewritten so it reads top to bottom without C++
shorthand. Each function, class and struct has a two-line overview above it.
Statements sit on one line even when long, `?:` became `if`/`else`,
constructors set their fields in the body instead of an initialiser list, and
the clock types use `auto`. In the serial stencil the PDE update is split into
five named terms (`stock_diffusion_term`, `cross_term`,
`variance_diffusion_term`, `stock_drift_term`, `variance_drift_term`) that are
then summed.

The baseline formula was left as one expression, laid out one term per line,
so its outputs stay bit-identical to M1. The serial rung now differs from the
baseline by at most 2.4e-15 on the Mac, inside the 1e-11 tolerance.

Splitting the PDE into named terms made the serial rung about 13% slower on
the Mac (1.33 to 1.16 Gcell/s, local, not a benchmark). Apple clang fuses a
multiply and an add into one FMA instruction only inside a single expression,
and the named terms are separate statements. Rebuilding with
`-ffp-contract=fast`, which fuses across statements, brought it back to 1.31
Gcell/s, which confirms the cause. Rangpur is unaffected because the build
there has no `-mfma`, so g++ cannot emit FMA at all. Debug job 650639 on
rangpur shows the serial rung still differs from the baseline by 3.3e-15 on
the call and exactly 0 on the put, the same as before the rewrite. This is
worth remembering for the SIMD rung, where FMA is used on purpose.

### Where this leaves the project

The serial rung is in place and matches the baseline. The SIMD rung adds
`update_row_avx` to `stencil.cpp` and a `SimdSolver` whose time loop calls
it.

## Phase 2. The SIMD rung

### What changed

There is a new rung, `--solver simd`. It runs the serial time loop but
updates each interior row with `update_row_avx`, which works on four
neighbouring cells per instruction using AVX2 intrinsics from
`<immintrin.h>`. The new files are `src/stencil_avx.cpp` (the AVX2 row),
`src/solver_simd.cpp` (the time loop) and `tests/test_avx_row.cpp` (the AVX2
row against the scalar row). `SimdSolver` is declared in `solver.h` and
`make_solver` knows the name `simd`. `test_ladder_matches` now checks
`serial` and `simd`. `slurm/debug.sh` also runs `make avx-check` and counts
AVX-512 instructions.

### Decision 1. The AVX2 row lives in its own file, the only one built with `-mavx2 -mfma`

PLAN §2a put `update_row_avx` inside `stencil.cpp`. It moved to its own file,
`stencil_avx.cpp`, because of the risk recorded in Phase 1 Decision 5. With
`-mfma` set, g++ is allowed to fuse any `a * b + c` into one FMA instruction,
in any file it compiles. Building `stencil.cpp` with `-mfma` would have
changed the serial rung's answers and timing as a side effect, and the
ladder would no longer show one change per rung.

The Makefile now has an explicit rule for `src/stencil_avx.o` that adds
`$(AVXFLAGS)`. `AVXFLAGS` is `-mavx2 -mfma` when `uname -m` is `x86_64` and
empty otherwise. Every other object file is compiled exactly as before.

`-mavx2 -mfma` was chosen over `-march=native` so the binary uses only the
Haswell features every r730 node has. `objdump -d heston | grep -c zmm`
printed 0 on rangpur, so no AVX-512 instruction slipped in.

### Decision 2. A scalar fallback behind one `#if`

The whole of `stencil_avx.cpp` sits inside
`#if defined(__AVX2__) && defined(__FMA__)`. The compiler defines those two
names only when it was given the flags. In the `#else` branch,
`update_row_avx` calls `update_row_scalar` and `avx_row_is_vectorised()`
returns false. `SimdSolver` prints a warning on stderr in that case. This
lets the arm64 Mac build and run every test, and nobody can mistake a Mac
run of the fallback for an AVX2 result.

### Decision 3. The scalar row takes a column range

The AVX2 loop handles cells in groups of four. A row has `ns - 2` interior
cells, which is often not a multiple of four, so the last zero to three cells
need the scalar maths. Rather than copy the scalar loop body into
`stencil_avx.cpp`, the body of `update_row_scalar` moved into a new function
`update_cells_scalar(constants, current, next, var_j, first_stock, end_stock)`
that updates only the columns it is given. `update_row_scalar` is now one
line that asks for columns 1 to `ns - 2`. The leftover cells of the AVX2 row
are therefore computed by exactly the serial rung's code. No arithmetic
moved, and the serial rung's answers on the smoke and demo grids and in
`test_ladder_matches` are unchanged.

### Decision 4. Unaligned loads everywhere

Every load and store uses `_mm256_loadu_pd` and `_mm256_storeu_pd`. The
stencil reads windows starting at `stock_i - 1`, `stock_i` and `stock_i + 1`,
which are 8 bytes apart, so at most one of them could be 32-byte aligned and
the aligned load would crash on the other two. Aligned buffers remain
experiment E2 for P5.

### Decision 5. FMA in the sum of the PDE terms

The five PDE terms are added into one running register `change` with
`_mm256_fmadd_pd`, in the same order as the scalar row's sum, and the
discount term is taken off with `_mm256_fnmadd_pd`. The final step
`centre + dt * change` is one more `fmadd`. The two second differences use
`fnmadd` for `east + west - 2 * centre`. That makes eight FMA instructions
out of 24 arithmetic instructions per four cells. Each intermediate value has
its own name so every line does one thing and can be read next to the scalar
version.

FMA rounds once where a multiply then an add rounds twice, so the SIMD rung is
not bit-identical to the serial rung. The test tolerance stays at 1e-11
because the measured gap is four orders of magnitude below it.

### Change 1. A row test for the AVX2 code

`tests/test_avx_row.cpp` fills a sheet with a call payoff plus random noise
(fixed seed), runs `update_row_scalar` and `update_row_avx` on every interior
row, and compares every cell, edges included. Both outputs start filled with
-7 so a stray write outside the interior would show up. The row lengths 130,
131, 132 and 133 leave a remainder of zero, one, two and three cells, and
2048 is the benchmark width. The test fails outright if the AVX2 row was not
built, because otherwise it would compare the scalar row with itself.

It runs with `make avx-check`. On the arm64 Mac that target uses
`clang++ -arch x86_64`, which makes an Intel program that runs under Rosetta.
On rangpur it builds natively with `mpicxx`. It is not part of `make test`
because on the Mac it needs that separate x86_64 build.

As a check of the check, changing the discount line from `fnmadd` to
`fmadd` made the test fail.

### Change 2. debug.sh checks the AVX2 build

`slurm/debug.sh` now runs a `simd` smoke solve, `make avx-check`, and prints
the number of `zmm` instructions in `heston` (must be 0) and `ymm`
instructions in `stencil_avx.o` (must not be 0).

### Validation

On the Mac (arm64, fallback row), `make test` passes and `simd` gives the
same answers as `serial`. Under Rosetta, `make avx-check` passes with a worst
gap of 1.2e-16 of the largest value for the short rows and 2.3e-15 for the
2048-wide row (tolerance 1e-13). An x86_64 build of `test_ladder_matches`
under Rosetta gave simd against baseline 1.06e-15 on the call and 0 on the
put.

Rangpur debug job 651064 (ran on `vcpu-2`, g++ 8.5.0 through `mpicxx`)
built with no warnings and passed everything.

| Check | Result |
|---|---|
| `test_ladder_matches` simd vs baseline | 1.059e-15 call, 0 put (tol 1e-11) |
| `test_ladder_matches` serial vs baseline | 3.296e-15 call, 0 put, unchanged from Phase 1 |
| `test_avx_row` | 1.167e-16 for ns 130 to 133, 2.305e-15 for ns 2048 (tol 1e-13) |
| `test_bs_collapse`, `test_parity` | 1.161e-3 and 1.116e-3, gap 5.377e-6, unchanged |
| `zmm` instructions in `heston` | 0 |
| `ymm` instructions in `stencil_avx.o` | 57 |

Rosetta and rangpur printed identical numbers in every line.

The smoke solve in the same job (256 by 64 grid, 20000 steps) took 2.43 s on
the baseline, 1.28 s on serial and 0.53 s on simd, so simd was about 2.4 times
serial. That node is a shared virtual machine and the job was not exclusive,
so this is a sanity check and not a result. The benchmark number comes from
r730-2 in P5.

### Where this leaves the project

The SIMD rung is in place and matches the baseline on rangpur. The OpenMP
rung adds a `#pragma omp parallel for` above the row loop of a copy of the
SIMD time loop, in `src/solver_omp.cpp`.

## Phase 3. The OpenMP rung

### What changed

There is a new rung, `--solver omp`, in `src/solver_omp.cpp`. It is a copy of
the SIMD time loop with one `#pragma omp parallel for` above the interior row
loop, and it records the thread count it got. `SolveResult` gained two
fields, `threads` and `ranks`, both defaulting to 1, and the CSV line gained
the same two columns at the end. `test_ladder_matches` checks `omp` against
the baseline and checks that it matches `simd` bit for bit, and `make test`
runs the ladder test again at 1, 2 and 7 threads. `slurm/debug.sh` now asks
for 8 cores and runs an 8-thread `omp` smoke solve.

### Decision 1. One pragma over the variance rows

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
```

The pragma goes on the outer `var_j` loop, so each thread gets a band of
whole rows. That is one fork and one barrier per step, each thread works on
one contiguous block of memory, and the inner loop is still the AVX2 row, so
the SIMD and thread speedups multiply. Threads read only `current` and each
writes only its own rows of `next`, so there is no data race and no lock,
`critical` or `atomic` anywhere.

`schedule(static)` was chosen because every row costs the same. When the
rows do not divide evenly the first threads get one extra row. A small test
program printed the same split (128, 128, 127, 127 for 510 rows on 4 threads)
under Apple clang's libomp and g++ 8.5's libgomp.

`default(none)` makes a forgotten variable a compile error instead of a
silent shared variable. PLAN §2b recorded a risk that GCC 8.5 and newer
compilers disagree about `const` variables in `shared()`. The list here holds
`current`, which is a reference to a const vector. GCC 8.5 on the rangpur
login node compiled it with `-fsyntax-only` and no warnings, and Apple clang
21 also accepts it, so no change was needed.

### Decision 2. The boundaries stay serial

The v = 0 row, the two stock edges and the top row run on the main thread
after the loop's implicit barrier. They are about 5100 cells per step
against about a million, roughly 0.5%. Leaving them serial keeps the code
shot to one pragma, and that serial part is the Amdahl term the strong
scaling sweep measures. Experiment E8 parallelises them later to see whether
the term is visible at all.

### Decision 3. The thread count in the result is what the program got

`result.threads = omp_get_max_threads()` is set before the time loop and
printed in the CSV. The team size can differ from what a script asked for,
for example when `OMP_NUM_THREADS` is unset and slurm limits the CPUs, so the
CSV records what the runtime actually used.

### Decision 4. The OpenMP rung must match the SIMD rung exactly

`test_ladder_matches` has a second list, `kSameBitsAsSimd`, holding `omp`.
Each rung in it must give a price, delta, gamma and vega exactly equal to the
SIMD rung's, with a tolerance of zero. This is a stronger check than the
1e-11 tolerance against the baseline. No value is ever added up across
threads, and each cell runs the same `update_row_avx` instructions it runs in
the SIMD rung, so any difference at all means a row was skipped, done twice,
or raced. `make test` runs the ladder test at 1, 2 and 7 threads. Seven
threads do not divide the 30 interior rows of the test grid, so the bands are
uneven, which is where an off-by-one at a band edge would show.

### Change 1. Two new CSV columns

`print_result_csv` now prints `threads` and `ranks` after
`cell_updates_per_sec`. They go at the end so any script that reads the
first twelve columns by position still works. `scripts/bench_plot.py` still
expects the M1 layout and is rewritten for M2 in P6. `make check-heston`
reads only the first column and is unaffected.

### Not done yet

Experiments E4 to E8 (persistent parallel region, dynamic schedule, thread
binding, first touch, parallel boundaries) need timings from an exclusive
r730 node, so they are run in P5 with the other sweeps, and the winner is
kept then.

### Validation

On the Mac, `make test` passes. `omp` against the baseline is 2.35e-15 on
the call and 2.24e-15 on the put, and `omp` against `simd` is exactly 0, at
1, 2 and 7 threads.

Rangpur debug job 651077 ran on r730-0 with 8 cores, built with no warnings,
and passed everything. With the real AVX2 row, `omp` against the baseline is
1.06e-15 on the call and 0 on the put, and `omp` against `simd` is exactly 0
at the default 8 threads and at 1, 2 and 7. `test_avx_row` and the zmm and
ymm counts (0 and 57) are unchanged from Phase 2.

The smoke solves in that job (256 by 64 grid, not an exclusive node, not a
benchmark) gave serial 171, simd 409 and omp on 8 threads 1838
Mcell-updates/s.

### Where this leaves the project

The OpenMP rung is in place and matches. The MPI rung splits the rows into
one slab per process, swaps one edge row with each neighbour per step, and
runs the OpenMP loop inside each slab.

## Phase 4. The MPI rung

### What changed

There is a new rung, `--solver mpi`, in `src/solver_mpi.cpp`. It is hybrid.
Each MPI rank owns one slab of variance rows, runs the OpenMP loop over its
own rows, and swaps one edge row with each neighbour per step. `main.cpp` now
starts and stops MPI, and only rank 0 prints. `SolveResult` gained
`exchange_seconds`, printed as a new last CSV column. `test_ladder_matches`
starts MPI, checks `mpi` against the baseline and bit for bit against
`simd`, and only rank 0 judges. `make test-mpi` runs it as 2 and 3 ranks.
`slurm/debug.sh` runs an `mpi` smoke solve and `make test-mpi`, and the new
`slurm/debug_mpi.sh` runs the ladder test across two nodes. The Makefile
adds `-DOMPI_SKIP_MPICXX`.

### Decision 1. Slabs of whole rows, every row shared out

All `nv` rows, boundaries included, are shared out with the first ranks
taking one extra row when they do not divide evenly (`slab_for_rank`). So
rank 0 owns the v = 0 row and the last rank owns the top row, and each
boundary is updated by exactly one rank. Each rank's interior rows are its
slab with row 0 and the top row cut off (`first_interior`,
`end_interior`). Rows were chosen over columns because a row is contiguous
in memory, so an edge row goes as one message straight out of the array with
no packing.

The rung refuses to run with fewer than two rows per rank, because the last
rank's `copy_top_row` reads row `nv - 2`, which must be its own. Every rank
makes the same check, so they all stop together.

### Decision 2. The edge row exchange runs on `next` at the end of each step

PLAN §2b put the exchange before the update. It is at the end instead, after
the boundaries and before `swap_buffers`, and it works on `next`. The effect
is the same, because the next step finds fresh copied rows in `current` either
way. Doing it on `next` means the exchange writes into a sheet the solver is
allowed to write, since `grid.current()` hands back a read-only reference,
and no new non-const accessor was needed on `Grid`. At step 1 every rank has
the whole payoff, so no exchange is needed before the first update.

`exchange_edge_rows` posts `MPI_Irecv` for each neighbour before
`MPI_Isend`, then one `MPI_Waitall`. Non-blocking calls make deadlock
impossible whatever order the ranks arrive in, and posting the receives first
lets the arriving row land straight in place. Both directions between a pair
use tag 0, which is safe because MPI matches on source as well as tag.

The exchange is timed on its own and added to `exchange_seconds`, which is also
included in `seconds`.

### Decision 3. Every rank keeps the whole grid

As PLAN §2b planned, each rank allocates both full sheets and only updates
its own rows. Row indices are then the same on every rank, so
`update_row_avx`, `update_zero_row`, `apply_stock_edges` (which already took
a row range for this) and `copy_top_row` are reused unchanged. The cost is
16 MiB per rank at 2048 by 512. A large-grid version would allocate only the
slab plus two copied rows and shift every index by `slab.begin - 1`.

### Decision 4. One gather at the end into rank 0's spare sheet

After the last step `current` holds each rank's finished rows and `next`
holds the stale previous sheet. `MPI_Gatherv` copies every slab into rank
0's `next` (each rank sends a different amount, hence `Gatherv`), then the
grid swaps so the gathered sheet is `current`, and rank 0 calls the shared
`extract_result`. This also avoided a writable `current` accessor. The gather
is after the timer, so it is not timed.

### Decision 5. One binary, MPI always started

`main` is now four lines. It calls `MPI_Init`, gets the rank, calls a new
`run(argc, argv, rank)` that holds the old body of `main`, then
`MPI_Finalize`. Moving the body into `run` means `MPI_Finalize` is called on
every exit path, errors included. Without `mpiexec` the program is a
singleton, so every other rung runs exactly as before. Under `mpiexec` with
a non-MPI rung each rank does the whole solve and only rank 0 prints.

`--dump-every` is refused by the `mpi` rung, because a snapshot on rank 0
would hold only its own rows.

### Decision 6. The MPI rung must match the SIMD rung exactly

`mpi` joined `omp` in `kSameBitsAsSimd`. Every cell runs the same
instructions, a copied row is an exact copy of the 64-bit values, and nothing
is summed across ranks, so the answer must be bit-identical. As a check of
the check, commenting out the edge row exchange made the two-rank answer move by
2.0e-9 and the test fail on both legs.

In the test, every rank runs every solver because the MPI rung needs all
ranks, but only rank 0 holds the gathered answer, so only rank 0 compares
and prints, and the other ranks return 0.

### Change 1. `make test-mpi`

Runs the ladder test with `mpiexec --oversubscribe --bind-to none` as 2 and
then 3 ranks, two threads each. Three ranks split the test grid's 32 rows
into 11, 11 and 10, which tests the uneven case. `--oversubscribe` is needed
on the Mac and inside the one-task debug job. `--bind-to none` stops each
rank being pinned to a single core, which would squeeze its threads onto it.

### Change 2. `-DOMPI_SKIP_MPICXX`

The first rangpur build printed `-Wcast-function-type` warnings from Open MPI
4.1's deprecated C++ bindings header, which `mpi.h` pulls in under g++ 8.5.
The code uses only MPI's C functions, so the Makefile now passes
`-DOMPI_SKIP_MPICXX`, which tells `mpi.h` to leave that header out. It sits
in its own `MPIFLAGS` variable so a script overriding `CXXFLAGS` keeps it.

### Change 3. Two-node debug job

`slurm/debug_mpi.sh` asks for 2 nodes with 2 tasks of 2 cores each, builds
once on the shared home directory, exports the two `OMPI_MCA_` exclusions
from the matrixass MPI script, and runs `test_ladder_matches` as 2 ranks
(`--map-by node`, one per node), 3 ranks and 4 ranks. It ends by printing
each rank's hostname to prove two machines were used.

### Validation

On the Mac, `make test` passes with `mpi` as a singleton at the default 14
threads and at 1, 2 and 7, and `make test-mpi` passes at 2 and 3 ranks. In
every case `mpi` against `simd` is exactly 0. The demo grid gives
196.1054098 at 1, 2 and 4 ranks.

Rangpur job 651092 ran on a100-a and a100-b (two ranks per node, hostnames
printed). `test_ladder_matches` passed at 2, 3 and 4 ranks with `mpi`
against `simd` exactly 0 and against the baseline 1.06e-15 on the call and
0 on the put.

Rangpur job 651091 ran on r730-0 and passed `make test` (singleton `mpi` at
8, 1, 2 and 7 threads), `make test-mpi` and `make avx-check`, with 0 zmm
instructions. That build still had the C++ binding warnings. Job 651116
rebuilt with `-DOMPI_SKIP_MPICXX`.

### Where this leaves the project

All five rungs exist and match the baseline, and the parallel ones match the
SIMD rung exactly. P5 rewrites `bench_common.sh` for threads and ranks,
writes the sweeps, and runs experiments E1 to E11 on exclusive r730 nodes.

## Phase 5. Benchmarks on rangpur

### What changed

`slurm/bench_common.sh` was rewritten for threads and ranks, and there are
new sweep scripts, `bench_m2_probe.sh`, `bench_m2_ladder.sh`,
`bench_m2_strong.sh`, `bench_m2_weak.sh`, `bench_m2_size.sh`,
`bench_m2_exp_a.sh`, `bench_m2_exp_b.sh` and `bench_m2_mpi.sh`.
`bench_peak.sh` also measures the AVX2 FMA peak and the triad memory
bandwidth. In the code, four experiment solvers were added in
`src/solver_experiments.cpp`, the omp and mpi rungs time their serial
boundary work into a new `boundary_seconds` CSV column, and `--vmax` sets the
top of the variance axis.

### Decision 1. Experiments are separate solvers, not switches in the rungs

Each experiment that needs different code is its own solver class, a copy of
the `omp` time loop with one thing changed (`omp_scalar`, `omp_persistent`,
`omp_runtime`, `omp_edges`). The alternative was a `--variant` flag with
`if` branches inside `OmpSolver`, which would have made the one-pragma code
shot for the video harder to read. As separate solvers, each experiment is
also checked by `test_ladder_matches` like any rung, and all except
`omp_scalar` must match `simd` bit for bit. That passed on the Mac at 1, 2,
7 and 14 threads, and `omp_runtime` also with `OMP_SCHEDULE=dynamic,1` at 7
threads.

`omp_runtime` uses `schedule(runtime)`, so one build can try `static`,
`dynamic,1`, `dynamic,8` and `guided` through `OMP_SCHEDULE`. That replaces
one solver per schedule.

### Decision 2. E7 uses numactl, not a first-touch rewrite

The plan for E7 was to have each thread first touch the rows it will
update, so the pages land on that thread's socket. That cannot work with the
current `Grid`. A `std::vector<double>(n)` writes zeros into every element
when it is made, so the main thread has already touched every page before
any parallel loop runs, and Linux does not move a page once it is placed.
Doing first touch properly would mean replacing the vectors in `Grid` with
uninitialised storage, which touches every rung. Instead E7 runs the same
binary under `numactl --interleave=all`, which spreads the pages evenly over
both sockets' memory. That tests the same question, whether page placement
matters, without changing any code.

### Decision 3. E2 aligned loads and E9 exchange overlap are not run yet

E2 needs 32-byte aligned buffers. `std::vector` only guarantees the
alignment of a double, so it needs an aligned allocator in `Grid`, and even
then at most one of the three loads per row (§15.4 of the study guide) can be
aligned. E9 needs the edge row exchange split into a start and a finish around
the interior update. Both are left for after the main sweeps. E9 is on the
PLAN §13 cut list.

### Decision 4. Timesteps chosen per run, maturity from the stability limit

Every rep should run for at least a second, so fast rungs get more
timesteps. Changing `nt` or the grid changes `dt`, so `bench_common.sh` has
`maturity_for ns nv vmax nt`, which prints the maturity that puts `dt` at 70%
of the same stability estimate `report_stability` prints. The rate in
Mcell-updates/s does not depend on the maturity, only on how many cells are
updated. The weak sweep grows `nv` by 21 rows per thread and grows `vmax`
with it so the variance spacing stays at 0.64/511.

### Decision 5. Pinned page regime, binding and two new columns

Every sweep exports `MALLOC_MMAP_THRESHOLD_=1048576`, the Milestone 1
headline regime, so every rep gets transparent huge pages, and passes it to
every MPI rank with `mpiexec -x`. `OMP_PLACES=cores` and `OMP_PROC_BIND=close`
are the default binding. Both go into the CSV header with the host, node
list, compiler, MPI version, flags, CPU, socket and NUMA layout and cache
sizes.

### Decision 6. Every job builds in its own copy of the tree

The first round of sweeps lost two jobs. `make sync` runs `rsync --delete`,
and I ran it twice while job 651152 (experiments E4, E5, E8, E11) was running
in the same directory. `--delete` removed the test binaries mid test and the
job's `slurm-651152.out`, so the job stopped at its test gate and wrote no
CSV. At the same time a diagnostic job on r730-3 ran `make clean` in that
same shared directory, which is what killed job 651153 (E6, E7) during its
build.

Two changes stop this happening again.

- `bench_common.sh`, `debug.sh` and `debug_mpi.sh` copy the tree into
  `~/cosc3500/work/job_<jobid>` and build and run there. The home directory is
  shared, so both nodes of an MPI job see the same copy. A `trap` removes it
  when the job ends. Results still go to the submitted tree's `results/`.
- `make sync` now also excludes `slurm-*.out` and the test binaries, so
  `--delete` cannot remove them. A dry run confirmed it deletes nothing.

### Decision 7. Denormals are flushed to zero in every timed run

The first size sweep (job 651151) showed the 2896 by 724 grid running 2.4
times slower than its neighbours on all three solvers, and the prices it
printed were as small as 1e-43. Short benchmark runs leave much of the grid
holding values far below 1e-300. Below about 2.2e-308 a double becomes
denormal (subnormal), and on this Xeon each arithmetic operation on one goes
through a slow microcode assist. How many such cells there are depends on
the grid and on `nt`, so the sweeps were comparing different numerical
regimes. Milestone 1 ruled denormals out for its own runs, but not for these
much shorter runs on bigger grids.

A counter-based check with `perf stat` did not count correctly inside a slurm
job. A direct check did work. Job 651284 ran the serial rung on 2048 by 512
at 381 steps with and without flush to zero, and got 146 and 171
Mcell-updates/s, 17% apart, with the same price to ten digits. So denormals
affect even the standard grid. Because the strong sweep uses
`nt = 1000 * threads` and the ladder gives each rung a different `nt`, they
would also bias those comparisons.

The fix pins the floating point regime, as M1 pinned the page regime.
`slurm/flush_denormals.c` is a five-line library whose constructor sets the
FTZ (flush to zero) and DAZ (denormals are zero) bits of the MXCSR register
before `main` runs. `bench_common.sh` builds it and loads it into every timed
run with `LD_PRELOAD`, and into every MPI rank with `mpiexec -x`. OpenMP
threads copy the setting from the thread that creates them. `heston` itself
is unchanged, the tests still run without it, and the CSV header records
`flush_denormals=1`. `FLUSH_DENORMALS=0` turns it off for one run.

New experiment E12 (`slurm/bench_m2_denormal.sh`) runs every ladder rung at
its ladder `nt` with and without it, so what denormals cost each rung is a
measured result in its own right.

The first-round CSVs (jobs 651147 to 651151 and 651154) were moved to
`results/superseded/`. Every sweep was resubmitted with the fix as jobs
651294 to 651303.

### Decision 8. The triad bandwidth benchmark had the first-touch problem

`bench_peak.sh` job 651154 measured 23.8 GB/s on 12 threads and only 24.8
GB/s on 24, which looks like the memory of one socket. The triad arrays were
`std::vector<double>`, which writes zeros to every element on the main thread
when it is created. That placed every page on socket 0 before the parallel
initialisation loop ran, the same first-touch effect described for `Grid` in
Decision 2. The arrays are now `new double[n]` held in `std::unique_ptr`,
which leaves the memory untouched until the parallel loop writes it, so each
page lands on the socket of the thread that will use it. The AVX2 FMA peaks
from that job, 43.4 GFLOP/s on one core and 951 on 24, do not depend on
memory and stand.

### Decision 9. A new figure script for Milestone 2

`scripts/m2_figures.py` draws F1 to F9 from the `bench_m2_*` CSVs and
`scripts/decomposition_diagram.py` draws F10. The Milestone 1
`scripts/bench_plot.py` was left alone so the M1 figures can still be
rebuilt from the M1 CSVs. The new script keeps the M1 house style
(`figstyle.py`) and the same rule that refuses any CSV whose host is not a
rangpur compute node. It fits Amdahl's serial fraction with a closed-form
least squares slope, because `1/S - 1/p = f (1 - 1/p)` is a straight line in
`f`, so no fitting library is needed.

### Finding 1. One grid size is slow because of where its two sheets land in memory

With denormals flushed, the size sweep (job 651299) still had one bad point.
The simd rung on 2896 by 724 ran at 170 Mcell-updates/s against 340 to 380
at every other size, while serial and the 24-thread omp rung were normal
there. Three diagnostic jobs on r730-3 (`slurm/diag_size_m2.sh`) narrowed it
down.

| Job | Grid (simd, 1 thread unless stated) | Mcell-updates/s |
|---|---|---|
| 651439 | 2832, 2880, 2912, 2960 by 724 | 402 to 407 |
| 651439 | 2896 by 512 and 2896 by 900 | 406 and 400 |
| 651439 | 2896 by 724 | 169, and 183 with `numactl --membind=0` |
| 651445 | 4096 by 512, 2048 by 1024, 8192 by 256, sheets of exactly 16 MiB | 380, 403, 360 |
| 651450 | 1448 by 1448 and 2044 by 1026 | 182 and 180 |
| 651450 | 2044 by 1027 | 236 |
| 651450 | omp on 2 threads, 2896 by 724 against 2880 by 724 | 336 against 814 |

Only grids whose two sheets end up exactly 2^24 bytes apart are slow. Each
sheet is its own `malloc`, which for a large size calls `mmap` with a 16 byte
header and rounds up to whole 4 KiB pages. 2896 by 724, 1448 by 1448 and
2044 by 1026 all have a sheet plus header just under 16 MiB, so each
mapping is exactly 16 MiB and the second sheet starts exactly 2^24 bytes
from the first. A sheet of exactly 16 MiB needs one more page for its
header, so 4096 by 512 is not affected, which is why the first guess at the
cause tested the wrong grids. 2044 by 1027, 20 KiB past the boundary, is
partly affected.

So the trigger is established, and the mechanism is not. Any two page
aligned sheets already share L1 and L2 cache sets, so plain set aliasing at
those levels cannot explain why only this distance is slow. The likely place
is further out, in how the L3 slices or the memory controller map physical
addresses, but proving that needs hardware counters, which `perf` did not
give inside a slurm job (Decision 7). The usual remedy is to pad one
allocation by a cache line or a page so the two sheets cannot sit an exact
power of two apart. It was not applied, because the benchmark grid
(2048 by 512) is not affected and every other size measured normally. The
point stays in the size figure with a note, rather than being dropped.

### Finding 2. `close` alternates sockets on this node

Job 651500 (`slurm/diag_placement.sh`) printed the CPU each OpenMP thread ran
on. With `OMP_PLACES=cores`, `close` on 12 threads used CPUs 0 to 11 and
`spread` on 12 threads used CPUs 0, 2, ..., 22. The r730 nodes number socket
0 even and socket 1 odd, so `close` alternates sockets and `spread` below 24
threads stays on socket 0. The first draft of F2 drew a line at 12 threads
labelled as where the second socket starts. It was removed, and the E6 and
study guide text now say which socket each binding uses.

### Change 1. E11 moved to its own job

Job 651300 hit the 14 minute limit before reaching E11, so E11 now runs in
`slurm/bench_m2_exp_c.sh` (job 651474).

### Validation

Every sweep job built from clean and passed `make test`, `make test-mpi` and
`make avx-check` before timing anything. The final data set is jobs 651294
(E12), 651295 (ladder), 651296 and 651297 (strong), 651298 (weak), 651299
(size), 651300 (E4, E5, E8), 651301 (E6, E7), 651474 (E11), 651302 (peak)
and 651303 (MPI layouts). The numbers and what they mean are in
`docs/RESULTS.md`, Milestone 2, and the hypotheses from PLAN §5 are scored
there. Four of the six held, the serial boundary (4) and two-node
efficiency (6) did not.

### Where this leaves the project

All sweeps are measured. Not run are E2 (aligned loads, needs an aligned
allocator in `Grid`), E9 (exchange overlap) and E3 (float32, phase P7).

## Phase 6. Figures

### What changed

`make m2-figures` draws F1 to F10 into `results/m2/` from the final CSVs,
named by job id in the Makefile, and writes `results/m2/summary.txt` with
every median it plotted.

| Figure | File | Shows |
|---|---|---|
| F1 | `f1_ladder.png` | the five rungs on a log axis with per rung and total gain, and every experiment as a ratio to the omp rung |
| F2 | `f2_strong.png` | speedup and efficiency against threads for the AVX2 and scalar rows, with the Amdahl fit |
| F3 | `f3_weak.png` | time per step with fixed work per thread, and scaled speedup against Gustafson |
| F4 | `f4_size.png` | serial, simd and omp 24 against working set with the cache bands, and omp over serial per size |
| F5 | `f5_roofline.png` | the measured scalar, AVX2 and DRAM roofs with every rung and the omp size points |
| F6 | `f6_breakdown.png` | time per step split into parallel rows, serial boundaries and edge row exchange |
| F7 | `f7_placement.png` | binding and page placement at 6, 12 and 24 threads on two grids |
| F8 | `f8_layouts.png` | the four hybrid layouts on two grids with the exchange share |
| F9 | `f9_correctness.png` | worst relative difference from the baseline per rung against the test tolerance |
| F10 | `f10_decomposition.png` | which cells each SIMD lane, thread and rank owns |

### Decision 1. Null rungs and nulls in the side panel are drawn quiet

As in M1, a rung within 3% of the one below is drawn light grey, and in the
side panel a ratio within 2% of 1 is light grey, a loss red and a gain blue.

### Decision 2. F6 uses two linear panels

The first version drew all seven bars on one log axis so the one-thread bar
(2.9 ms) and the 2-node bar (0.12 ms) fitted together, but the log axis hid
the parallel part of the shortest bar. Two linear panels, one for OpenMP
and one for the MPI layouts, show every part.

### Decision 3. The 2896 by 724 simd point is annotated, not dropped

F4 points at it with a note sending the reader to RESULTS.md M2.8.

### Change 1. Code shots for the M2 rungs

`scripts/code_shots.py` still pointed at `src/solver_opt_kernels.cpp`,
which M2 deleted, so `make code-shots` crashed. It now draws seven shots
into `results/m2/code_shots/`. The M1 shots in `results/code_shots/` were
made from the code at tag `milestone1` and are left alone.

| Shot | Shows |
|---|---|
| `shot_serial_row.png` | the scalar derivatives, five terms and time step in `src/stencil.cpp` |
| `shot_simd_fma.png` | the five scalar terms beside the FMA chain in `src/stencil_avx.cpp` |
| `shot_omp_pragma.png` | the simd row loop beside the same loop with the `parallel for` pragma |
| `shot_mpi_exchange.png` | the `Irecv`, `Isend` and `Waitall` lines of `exchange_edge_rows` |
| `shot_iron_rule.png` | `build_and_validate` in `slurm/bench_common.sh` |
| `shot_csv_header.png` | the provenance header of the ladder CSV, job 651295 |
| `shot_ladder_test.png` | the call rows of `test_ladder_matches` at 24 threads on r730, from the job 651295 log |

The M1 script picked its lines by number, so any edit above a shot moved
it onto the wrong code. The new script finds each range by searching for a
line of text, such as `"The five terms of the Heston PDE"` to
`"out_row[stock_i] = centre + dt * change;"`, and highlights lines by text
too. If a search text disappears from the source the script stops with an
error naming the file, rather than drawing a shot of the wrong lines.

The test shot uses the rangpur log, not a local run, because on the arm64
Mac the simd rung falls back to the scalar row, so its bit for bit match
with omp proves less there.

## Phase 7. float32 (experiment E3)

### What changed

A new experiment solver, `--solver omp_float`, runs the `omp` loop with
both sheets stored as `float` instead of `double`. Run on one thread it is
compared with `simd`, and on 24 threads with `omp`.

| File | Change |
|---|---|
| `include/stencil.h`, `src/stencil.cpp` | new `FloatStencilConstants`, built by rounding every double constant once, and `update_cells_scalar_float`, the scalar row in float |
| `src/stencil_avx.cpp` | new `update_row_avx_float`, the AVX2 row with `__m256` (eight floats) in place of `__m256d` (four doubles), plus its non-AVX fallback |
| `src/solver_float.cpp` (new) | `OmpFloatSolver`, with its own float zero row, stock edges and top row copy |
| `include/solver_experiments.h`, `src/solver.cpp` | the class and its `make_solver` name |
| `tests/test_ladder_matches.cpp` | `omp_float` against the baseline with its own limit, 1e-4 |
| `tests/test_avx_row.cpp` | the float AVX2 row against the float scalar row, row lengths 130 to 137 and 2048 |
| `slurm/bench_m2_exp_float.sh` (new) | the E3 sweep |
| `scripts/m2_figures.py`, `Makefile` | E3 in the F1 side panel, `omp_float` in F9, F9 now drawn from the job 652647 log |

### Decision 1. Float is an experiment solver, not a rung

It can never meet the 1e-11 limit every rung is held to, so it sits with
the other experiments and has its own limit in the test.

### Decision 2. The float code is copied, not templated

The float row, zero row and edges are copies of the double ones with
`float` in place of `double`. A C++ template could make one function serve
both types, but it would change the shared double functions that five
solvers must match bit for bit, and templates are harder to read for a
beginner. The copies sit next to the originals so they are easy to compare.

### Decision 3. The constants and the edge values are worked out in double first

`FloatStencilConstants` is built from the finished `StencilConstants`, and
the stock edge values are computed in double and then rounded. So both
solvers start from the same numbers, and the gap between them comes only
from doing the time steps in float.

### Decision 4. The answer goes back through the grid's double sheet

At the end the float sheet is copied into the grid's spare double sheet,
which is then swapped in, and `extract_result` reads the price and Greeks
as for every other solver. The grid class did not have to change.

### Decision 5. Rejected for accuracy, not speed

Float is 1.75x faster at 2048 by 512 and 4.46x at 4096 by 1024 on 24
threads, but on the reference grid its price is $0.044 from the Fourier
price where double is $0.0008, and gamma moves by 11.9%. The error comes
from rounding every cell to float on every one of the 980000 steps.

### Validation

On the Mac, `make test` passes with `omp_float` at 2.1e-6 (call) and
8.3e-6 (put) from the baseline, and `make avx-check` passes the float row
at 6.3e-8 of the largest value for every remainder 0 to 7 and 1.2e-6 on
the 2048 row (limit 1e-5). The double results did not move, and every
rung still matches `simd` bit for bit.

Rangpur job 652647 (r730-2) built with no warnings, passed every test
including 2.1e-6 and 8.2e-6 for `omp_float`, and measured, in
Mcell-updates/s,

| Run | double | float | float / double |
|---|---|---|---|
| 1 thread, 2048 x 512 | 383 (simd) | 735 | 1.92x |
| 24 threads, 2048 x 512 | 5732 | 10048 | 1.75x |
| 24 threads, 4096 x 1024 | 2902 | 12945 | 4.46x |
| 24 threads, reference grid, one run | 6037, price 196.1684 | 10754, price 196.1256 | 1.78x |

### Where this leaves the project

P7 is done. E2 (aligned loads) and E9 (exchange overlap) still have not run.
The guides cover the float solver in CODE_GUIDE §9e and STUDY_GUIDE §15.9.

## Phase 8. Documents

### What changed

`COMMANDS.md` was run command by command on the Mac. Every local command
worked. Four things were wrong or missing and were fixed.

- The first example ran `config/reference.cfg` with the serial rung, which
  is about 1e12 cell updates and 16 minutes on the Mac. It now runs
  `config/bench.cfg`, about 2 seconds, with a note on the reference grid.
- The column header shown first was the old ten-column one.
- The `--solver` row did not list the experiment solvers.
- The sweep table still had E11 in `bench_m2_exp_a.sh` and did not list
  `bench_m2_exp_c.sh`, `bench_m2_exp_float.sh` or `bench_m2_denormal.sh`.

A new §5 is the five minute interview demo. It starts the compute node
check first, then builds and tests on the Mac, runs each rung on the same
grid, shows a CSV provenance header and fetches the debug job's output.

Every Milestone 2 decision in this file was checked for an answer in
`STUDY_GUIDE.md`. Three had none (Phase 1 Decisions 3 and 5, Phase 5
Decision 1), and they are now questions 17 to 19 in the Part V drill.

`CODE_GUIDE.md` §12.3 showed the old two-entry `kSameBitsAsSimd` list and
§1.1 did not list the experiment files. Both now match the code.

### Validation

Every local command in `COMMANDS.md` §1, §2, §2a, §2b and §5 was run on
the Mac. `make queue` and `make fetch` were run against rangpur. The §5
runbook was run end to end. Its debug job 652682 passed `make test`,
`make test-mpi` and `make avx-check` on a compute node, with 0 `zmm`
instructions in `heston` and 114 `ymm` instructions in `stencil_avx.o`, up
from 57 because the float row now sits in the same file.

## Phase 9. Presentation plan, and what the Milestone 1 feedback asked for

### What changed

The Milestone 1 feedback (82.6/100) is written into `PLAN.md` §11a, each
comment next to the change it leads to. Three of them asked for new work.

| Feedback | New work |
|---|---|
| No data structures, half the video before any HPC | `shot_data_layout.png`, `index()` beside the loop nest with the row and sheet sizes, and a 40 second data structure beat at 1:20 |
| A flat profile of the baseline | `slurm/profile_m2.sh`, gprof on baseline, serial and simd, drawn as F12 |
| A four-bar -O0 to -O3 chart for the baseline | `slurm/bench_m2_optlevels.sh`, baseline and serial at each level, drawn as F11 |

`presentation/M2_SCRIPT.md` is the word-for-word script, 14 takes, about
1215 words, 8 minutes 40 at 140 words per minute, each take checked against
its slot. It opens with the checklist the presentation feedback led to (no
overlay over slides, a large cursor, a pronunciation list, re-record a take
that runs long).

`scripts/m2_figures.py` gained F11 and F12 and a `--profile` option, and
`make m2-figures` passes the new CSV and the three profile files.

### Decision 1. The opt level job refuses a build whose answer moved

Each level is built, run on the smoke grid, and compared with the `-O2`
answer field by field before it is timed. All four matched exactly, as
they should, since no level is allowed to reorder floating point.

### Decision 2. Profile with gprof, not perf

`perf stat` did not count properly inside a slurm job on rangpur (Phase 5),
and gprof only needs `-pg` and the `gprof` binary, which the nodes have. It
samples every 10 ms, which is enough for runs of several seconds.

### Decision 3. Serial is in the opt level chart too

The feedback asked for the baseline. The serial rung is next to it because
its `-O3` bar answers experiment E1, how much of the SIMD rung the compiler
could have done for us.

### Validation

Job 652696 (r730-2), every level gave the `-O2` answer exactly, baseline
10.7, 57.6, 57.7, 112.5 and serial 44.8, 156.1, 158.4, 238.6
Mcell-updates/s from `-O0` to `-O3`. Job 652697 (r730-2), baseline 99.9% in
`BaselineSolver::solve`, serial 99.4% in `update_cells_scalar`, simd 98.9%
in `update_row_avx`. Both jobs passed `build_and_validate` first. Numbers
in RESULTS M2.11, concepts in STUDY_GUIDE §18.7.

### Where this leaves the project

The slides for M2 are not built yet. The script names the figure or code
shot for each take, so the deck is one slide per take in the same order.


## Phase 9b. The Milestone 2 deck and the HPC justifications

### What changed

`presentation/deck_m2.tex` is the Milestone 2 recording deck, built with
`make slides-m2` into `presentation/deck_m2.pdf`. It has 22 frames, a title
frame and then one frame for each on-screen item in `M2_SCRIPT.md`, in the
script's order, so a take that says "F12, then F11" has two frames. The
Milestone 1 deck (`deck.tex`) is left as it was so its PDF can still be
rebuilt.

`HPC_JUSTIFICATIONS.md` is new. It answers why each HPC technique was used,
one question per heading, in full sentences that can be said in the
interview. Each answer shows the code the work was done in, copied from
the real file with its file and line numbers, explains the lines that
matter one by one, then gives the reason and the number that backs it, with
the job behind it in `docs/RESULTS.md`. Every snippet line was checked
against its source file by a script. It is linked from `README.md` and PLAN §9 and §12.

### Decision 1. A new deck file, with the M1 preamble copied in

The new deck copies the palette, theme and `sig*` macros from `deck.tex`
so the two decks and the figures look like one set. Copying rather than
sharing a preamble file keeps the M1 deck byte for byte as it was
submitted.

### Decision 2. The bottom-right corner is kept for the webcam

The Milestone 1 marker found a streaming overlay covering slide headers.
The callout band at the foot of every frame is redefined at 74% of the
width, so the bottom-right corner never holds slide content and the webcam
box can sit there.

### Decision 3. One callout sentence per frame, taken from the script

Each frame carries the figure and one or two sentences with the number the
script says over it, so the slide and the voice agree. The conclusion frame
shows the ladder beside the measured ceilings, which the Milestone 1
feedback asked for.

### Validation

`make slides-m2` compiles with no overfull boxes. The one remaining
warning, a missing character in the title frame, also appears in the M1
deck and does not change the output. Every number on the slides was checked
against `docs/RESULTS.md` Milestone 2 and `results/m2/summary.txt`. The two
node total of 155x is 8768 / 56.6 from jobs 651303 and 651295.

### Where this leaves the project

P9 needs only the recording. P10 is the interview rehearsal and the final
pass of both guides.

## Phase 9c. The equation, the grid and the Black-Scholes check in the M2 video

### What changed

The M2 deck and script assumed the marker had seen the Milestone 1 video.
The task sheet says each presentation must stand on its own. Three frames
now follow the problem frame. The equation frame shows the Heston PDE with
its terms grouped into randomness and drift, the central differences for
V_S, V_SS and V_Sv in grid indices, and the explicit update. The grid frame
shows F13 and the check frame shows F14. `scripts/method_figures.py` draws
both and `make m2-figures` now runs it.

The cold open frame (F1) was removed. Its line is now spoken over the title
frame, and F1 still appears in the rejected and close takes.

A frame on how the numbers were measured was added before the profile. The
task sheet asks how the benchmarks avoid misleading results, and no M2 take
said it. It gives the node, `--exclusive`, the median of 5, pinned threads,
the tests that run first in each job, the printed price, and the first
round of sweeps that was thrown away over denormals. The weak scaling take
and its callout now name Gustafson's law, which predicts 22 at 24 threads
with the measured 8.4 percent serial share, against 19.1 measured. The
course does not teach Gustafson's law, so the slide gives the formula and
cites the source, J. L. Gustafson, Reevaluating Amdahl's law,
Communications of the ACM 31(5), 532 to 533, 1988.

The same rule now covers every idea the slides use from outside the lecture
material on disk (L01, L03 to L06, T04 to T06). A `\sigcite` macro puts a
small grey source line under the callout, 74 percent wide so it stays clear
of the webcam corner. The protocol slide cites Goldberg (1991) for denormals,
both roofline slides cite Williams, Waterman and Patterson (2009), strong
scaling cites lecture 3 for Amdahl and Karp and Flatt (1990) for fitting the
serial fraction, and the placement and hybrid MPI slides cite Drepper (2007)
for first touch. The strong scaling callout now gives the fitted formula. The
figure slots on those slides shrank to make room. M2_SCRIPT.md lists the
sources.

F13 also labels how the cells are calculated. The expiry sheet says the value
is known there, with the zero and S minus K payoff regions labelled. The
arrow reads one step back. A key between the sheets shows which derivative
each of the nine cells feeds, and the update formula sits under the next
sheet. The grid callout and take now say why the solve starts at expiry.

| Figure | File | What it shows |
|---|---|---|
| F13 | `f13_grid.png` | one step, the current sheet at T holding the payoff, the red block of nine giving one cell of the next sheet at T - dt, the next sheet filling row by row, and its edges |
| F14 | `f14_bs_collapse.png` | the solver at xi = 0 and rho = 0 drawn over the closed-form Black-Scholes price today, with the gap underneath |

### Decision 1. F13 shows one step, not the whole solve

The first version showed a grid with a stencil and four small sheets from
expiry to today. It was redrawn to show only the step from T to T - dt,
with the same i and j outlined in both sheets, so the setup is clear and the
speaker can say the step repeats back to t = 0.

### Decision 2. F14 is redrawn from the saved collapse run

The Milestone 1 still of this comparison is a frame of a dark video. F14
reads the last snapshot in `results/collapse_work/dump`, the same 526 by 101
run, and draws it in the light figure style. It reuses `bs_call` and
`read_variance_row` from `scripts/bs_collapse_anim.py`. At S = 5200 the
solver gives 202.7705 and the closed form 202.8080, a gap of 3.75 cents or
1.85e-4 relative. The largest gap between 0.7 and 1.4 times the strike is
0.0375 dollars. The convergence plot was left out, because its error
flattens onto a floor and would take more time to explain than the slot has.

### Decision 3. Other takes were cut to pay for the new ones

The data structure take dropped the loop nest and swap sentences, which F13
now shows. Where the time goes dropped the boundary sentence, which is still
on its slide. SIMD dropped the O3 vectoriser sentence, which is on its
slide. OpenMP, MPI, weak scaling and the same-answer take each lost a
clause. The script is 1259 words, 9 minutes at 140 words a minute, and the
slots end at 9:40.

The M1 feedback asked for HPC content by 1:20. The grid take at 1:15 now
carries the two sheets and why every cell in a step can be done at once,
and the data structure take follows at 2:05.

### Validation

`make slides-m2` compiles with no overfull boxes. The equation, grid,
check and protocol frames were rendered and read for overflow. Every claim
on the protocol frame was checked against `slurm/bench_common.sh` and the
ladder job script. `make test` in `bench_common.sh` runs
`tests/test_bs_collapse` before any timing.

## Phase 9d. Before and after snippets in the code guide

### What changed

`CODE_GUIDE.md` has a new §15. It takes each code change in this file in
order and shows the code before it next to the code after it, with what
changed and why. It covers the renames, the hoisted constants, the named
PDE terms, the five-call time loop, the column range, the AVX2 row, the
simd, omp and mpi loops, the four experiment solvers, the denormal flush,
the triad arrays and the float row. The reading order note at the top of
the guide points to it.

§2.2 also explains `std::size_t`, what it is, which header it comes from
and how it differs from `int`. §4.4 and §5.2 say why `find` is compared
with `npos` and why the cell count is cast to `std::size_t`, and §14 has a
question on it. The Heston Code Walkthrough page was rebuilt from the guide.

### Change 1. The package secret check no longer matches itself

The `make package` guard greps the payload for the word pass word in
capitals, or in lower case followed by an equals sign. The Makefile is in the
payload, and its own grep line held that word, so the guard refused every
package. The pattern now wraps the first letter in square brackets, which
matches the same words but not the line that holds it. A copy of the payload
passed the new check, and a test file holding a fake setting was still
caught.

### Validation

A script checked every line of every snippet in §15 against the source
files and the `milestone1` baseline. Every line matched except the two
before snippets marked as rebuilt from this file (§15.5 and §15.12), whose
earlier versions were never committed on their own.

## Phase 9e. The study guide page

### What changed

The Heston Code Walkthrough page became the Heston Pricer Study Guide, at
the same link, so it can be studied from instead of STUDY_GUIDE.md. A new
file, `docs/STUDY_CARDS.md`, holds 40 short cards in seven lessons. Six
cover STUDY_GUIDE.md (finance, the PDE and grid, memory and serial speed,
benchmarking, the parallel rungs, results). The seventh, The C++ in the
code, lists every standard library call and C++ feature the code uses,
each with its real line, what it does and why it is used there. The page
puts these cards before the CODE_GUIDE.md lessons, and links like
STUDY_GUIDE §7.4 in the code pieces now open the matching card.

The drill has 94 questions, 55 on the ideas and the C++ with written
answers and 39 on the code that open the piece answering them. A search box
finds every piece that mentions a name, such as `std::clamp` or `npos`. In
test-yourself mode a C++ card hides its why column, a concept card hides
everything after its first paragraph, and a code piece hides its
explanation.

STUDY_GUIDE.md is unchanged and is still the long version. The cards are
about 10700 words against its 27700.

## Phase 9f. Simpler settings code

### What changed

The command line and config file code was rewritten with fewer C++ ideas,
so it can be explained line by line. Nothing inside the time loop changed.

`main.cpp` now requires `--config <file>` as the first two words, loads it,
then reads every other flag as a pair with `i = i + 2`. The first pass that
searched for `--config` anywhere and the `flag_value` helper with its `int&`
parameter are gone. Every script, slurm job and Makefile target already put
`--config` first.

`params.cpp` now splits each line into three words with
`std::istringstream`, like `line.split()` in Python, and converts values
with plain `std::stod` and `std::stoi`. The `trim`, `parse_double` and
`parse_int` helpers, and the `find`, `npos` and `substr` calls, are gone. One
`try` and `catch` around `apply` turns a value that is not a number into a
message with the file and line. Unknown keys and flags still stop the
program.

One check was given up. A value that starts with a number and has text
after it, such as `12x`, is now read as 12 instead of being rejected. Every
config line must have spaces around the `=`, which every config file and
every config the Python scripts write already has.

### Validation

Four runs on `smoke.cfg` and `demo.cfg` across four solvers and every flag
type gave the same price, delta, gamma and vega to every printed digit
before and after. `make test` and `make test-mpi` pass. A missing value, an
unknown flag, an unknown key, a missing `=`, a non-number and a missing file
each stop the program with exit code 1 and a message.

CODE_GUIDE.md §2.2, §3 and §4, docs/STUDY_CARDS.md §I7.2 and §I7.3,
HPC_JUSTIFICATIONS.md §6 and DESIGN_DECISIONS.md decisions 2 to 4 were
updated to match.

## Phase 9g. Beginner pass over every source file

### What changed

Every file in `src/`, `include/` and `tests/` was read through for C++ that a
reader who knows C and Python would have to look up, and for comments that
restate the code or use filler words. Nothing the solvers compute changed.

The step timer is now `omp_get_wtime()`, which returns seconds as a
`double`, in place of `std::chrono::steady_clock` with `auto` and
`duration<double>(...).count()`. Every file is already built with OpenMP, so
the baseline, serial and simd solvers include `<omp.h>` for it. On Linux both
timers read `CLOCK_MONOTONIC`, so the benchmark numbers taken with the old
timer still describe this code.

`locate` in `grid.cpp` keeps values in range with plain `if` statements
instead of three `std::clamp` calls. `grid.cpp` counts cells in an `int`
instead of casting to `std::size_t`, since `index()` already returns an `int`
and the largest grid run has 8.4 million cells.

The float constants and the float payoff sheet are filled with plain loops
and `(float)` casts instead of iterator ranges (`.begin()`, `.end()`), and
`explicit` was dropped from the `FloatStencilConstants` constructor.

`dump_snapshot` writes with C's `fopen`, `fprintf("%.10g")` and `fclose`
instead of `std::ofstream` with `<<` and `precision(10)`.

Constants named in Google style (`kNodeSnapCells`, `kRungs`,
`kRelativeTolerance` and the rest) now have plain names such as
`node_snap_tolerance` and `solver_names`. The tests print OK or MISMATCH with
`if`/`else` instead of `?:`, and still fail on NaN.

Comments were reworded where they used filler (simply, clean, fresh,
actually), plan words (rung), or a speed claim that was out of date (serial
is 2.8 times the baseline, not twice). New one-line comments explain the C++
that a C or Python reader meets for the first time, such as `Grid::`,
`: public Solver`, `override`, `&` on a return type, `make_unique` and the
range `for` loop.

### Validation

Every solver on `smoke.cfg`, call and put, the MPI solver on three
processes, and six off-node and edge quotes that reach every branch of
`locate` gave the same price, delta, gamma and vega to every printed digit
before and after. All 400 snapshot files from an omp run with
`--dump-every 50` are byte for byte the same, so `%.10g` writes what the
stream did. `make test`, `make test-mpi` and `make avx-check` (under
Rosetta) pass, and the build has no warnings.

A dry run of `make package` found that the Phase 9 entry on the package
guard quoted the words the guard searches for, so the guard would refuse the
zip. That entry now describes the pattern without quoting it.

CODE_GUIDE.md §2.2, §5.2, §5.6, §7.3, §7.10, §9b.3, §9c.6, §9e, §10.2 and §12,
docs/STUDY_CARDS.md §I7, DESIGN_DECISIONS.md decisions 9, 20, 21 and 23 and
the documentation table, HPC_JUSTIFICATIONS.md and docs/WALKTHROUGH.md line
references, and the timer note in docs/RESULTS.md were updated to match.

## Phase 9h. Comments in my own voice

### What changed

Every comment in `src/`, `include/` and `tests/` was rewritten to read the
way I prefer to read it, as short plain statements that end in a full
stop. Earlier drafts used first person and padded most lines with an added
clause such as ", so the time loop only has to read them". Both were taken
out, along with "This is the X that" openers, "not X but Y" contrasts and
trailing "which is" asides. A reason stays only where the reader needs that
fact, such as "so no locks are needed". The old rule of 3 to 12 word
fragments is retired. The voice is written down in `../CLAUDE.md` under "How
the user writes", and `PLAN.md` §3 points to it.

Only comment lines changed. No comment was added or removed, so every line
number stays where it was and the line references in `HPC_JUSTIFICATIONS.md`
and `docs/WALKTHROUGH.md` still hold.

### Validation

With every comment line removed, each file is identical to the version
before. After a forced rebuild, all 15 object files are byte for byte the
same as before, so the program the benchmarks ran is the program in this
tree. `make test` and `make test-mpi` pass. No comment has a semicolon or a
dash, and every one ends in a full stop.

The snippets in `HPC_JUSTIFICATIONS.md` and `DESIGN_DECISIONS.md` that quote
comments were updated to the new text.

## Phase 9i. Names that say what things are

Names in the code that were metaphors were replaced by names that say what
the thing is.

| Before | After | Where |
|---|---|---|
| `RowSlab`, `slab_for_rank`, `slab` | `RowRange`, `rows_for_rank`, `owned_rows` | `solver_mpi.cpp` |
| the exchange function | `exchange_edge_rows` | `solver_mpi.cpp` |
| its two received row names | `row_below`, `row_above` | `solver_mpi.cpp` |
| its message tag | `message_tag` | `solver_mpi.cpp` |
| its timer start, end and total | `exchange_start`, `exchange_end`, `exchange_seconds` | `solver_mpi.cpp`, `solver.h`, `io.cpp` |
| `my_rows` (the Gatherv send pointer) | `my_values` | `solver_mpi.cpp` |
| `sheet` parameter | `buffer` | `solver_mpi.cpp` |
| `blend_cell` | `bilinear_interpolate` | `solver.cpp` |

Comments that said sheet now say buffer, matching `swap_buffers`. The
comment above `exchange_edge_rows` names the standard HPC term for the
received rows once. The north, south, east and west
neighbour names in the stencil stayed, because they are the usual textbook
labels for a stencil's neighbours and a comment says which way each points.

The CSV column for the exchange time is now `exchange_seconds` in the header that
`slurm/bench_common.sh` writes and in `scripts/m2_figures.py`. The figure
script names the columns itself, so older CSVs still read. The F6 and F8
labels and the MPI code shot (now `shot_mpi_exchange.png`) follow the new
names. The study guides keep sheet and slab in their prose, with a note in
`STUDY_GUIDE.md` §6.1, §17.2, §17.3, `CODE_GUIDE.md`
§9c and `docs/STUDY_CARDS.md` §I5.6 giving the code name for each.

### Validation

No line was added or removed in any source file, so line references hold.
After a forced rebuild the machine code (`otool -t`) of all 15 object files
is identical to before the rename. `make test` and `make test-mpi` pass. The
code snippets in `CODE_GUIDE.md`, `HPC_JUSTIFICATIONS.md` and
`DESIGN_DECISIONS.md` match the new names.


## Phase 9j. The code guide follows the logic line by line

### What changed

No code changed. `CODE_GUIDE.md` has a new §1.5 that traces one run from
`main` to the printed CSV line. It has a call tree with the file and line
of every function, then four numbered walks, from the command line to a
solver, inside one solve, one timestep and one row. Each step links to the
exact line in the editor and says why the step is there. A table lists,
for each of the 10 solvers, where its `solve` starts, which row function it
calls and the lines that differ from the serial solver. A list gives where
each test's `main` starts.

Blocks were added for code the guide skipped or hid behind `...`.

| Section | New or filled in |
|---|---|
| §3.2, §3.2b | the `--type` check, `--vmax` and `--maturity`, `throw`, and `usage` |
| §4.1b | `MarketSpec`, `HestonParams` and `GridSpec`, and why every field has a default |
| §5.6, §5.7 | the four parts of `locate` with a worked example, and `stock_position` and `variance_position` |
| §6.5, §6.6 | the far corner locals, the Feller text, the off-node warning and the four neighbours read for the Greeks |
| §7.1 | the baseline's local names, `dt` and `max_stock_price` |
| §8.1b, §8.3 | the scalar part of the `StencilConstants` constructor, `resize` and the bracketed coefficients |
| §8.8 | the full `apply_stock_edges`, with the value at each end of the stock axis |
| §9a.4, §9a.5 | the row pointers and every derivative line of the AVX2 row |
| §9c.3 | the whole grid on every rank and the stability line on rank 0 |
| §9d.1, §9d.3 | `omp_scalar` and `omp_runtime`, with the E5 numbers |
| §9d.2 | why `omp_persistent` times the whole loop once |
| §9e.4, §10.2 | `update_zero_row_float`, and the `NULL` check and `snprintf` in `dump_snapshot` |
| §12.0, §12.3, §12.4 | how each test is built, why the test grids are 421 by 51, `max_relative_diff` and the random fill |

§9d was renumbered, so `omp_persistent` is now §9d.2 and `omp_edges` §9d.4.
The §1.1 map lists `OmpSolver`, `MpiSolver` and `solver_experiments.h`.
The §15.9 After snippet still had the MPI names from before Phase 9i and
now uses `RowRange`, `rows_for_rank`, `owned_rows` and
`exchange_edge_rows`. §14 has nine new questions. `STUDY_GUIDE.md` points
to §1.5.

### Validation

A script checked every snippet line in `CODE_GUIDE.md` against the source.
The only lines not found are the deliberate Before snippets in §15 and
lines that stand for code left out. A second script opened every line link
and every call tree line and printed the line it lands on, and each one is
the line named. The off-node claim was checked with a run on
`config/bench.cfg`, which does print the warning, and the guide says so.
`tests/test_ladder_matches` was timed at about 7 seconds on the Mac.

The Heston Pricer Study Guide page is built from `CODE_GUIDE.md` and has
not been rebuilt for this phase.

## Phase 9k. The rerun on the current code, and one point that moved

### What changed

The M2 sweeps were submitted again on the current code as jobs 656955 to
656968. When this entry was written, five had finished (ladder 656955,
strong 656956 and 656957, weak 656958, size 656959), `exp_a` (656960) was
running and the other eight were queued. Every finished job passed all its
tests before timing, and every price matches the first runs to every
printed digit.

| Sweep | First run | Rerun | Largest change |
|---|---|---|---|
| ladder | 651295 | 656955 | mpi 2 x 12, 6019 to 6272 (+4.2%) |
| strong omp | 651296 | 656956 | 24 threads, 5733 to 5973 (+4.2%) |
| strong omp_scalar | 651297 | 656957 | 8 threads, +3.6% |
| weak | 651298 | 656958 | 10 threads, +3.1% |
| size | 651299 | 656959 | serial at 2896 x 724, 167 to 58 (-65%) |

Everything except the size sweep moved by 4.2% or less, which is ordinary
run-to-run noise, so F1, F2, F3 and the numbers in the deck still stand and
were not redrawn.

Serial at 2896 x 724 is now hit by the same address effect as simd
(Phase 5, Finding 1), with all five reps between 58.3 and 58.5. Phase 9g
said the old numbers still describe the current code, and for this one
point they do not. The cause of the change was not found. F4 and F5 are now
drawn from job 656959 (`M2_CSVS` in the Makefile). F4's note points at the
serial point and names both rungs, and the 2896 bar in the speedup panel is
grey and labelled serial slow, because 111x there is the slow serial run
and not a real speedup. `docs/RESULTS.md` M2.4, M2.5 and M2.8 and
`HPC_JUSTIFICATIONS.md` (177 and 97 GFLOP/s on the roofline) were updated.
The deck was rebuilt for the new F4 and F5, and no slide text changed.

### The code shots and the deck were out of date

`presentation/deck_m2.pdf` had last been built on 12 October, before
Phases 9f to 9i, so its code shots showed the code before the comment
rewrite and the renames. `deck_m2.tex` already pointed at
`shot_mpi_exchange.png`, which had never been drawn, and the deck no longer
compiled. `scripts/code_shots.py` looked for the old comment
"Every interior row, four cells at a time", and now looks for "Update every
interior row, four cells at a time.". `make code-shots` drew all eight
shots from the current source, and the old MPI code shot was deleted.
`make slides-m2` now builds the 25 frames with no errors or overfull boxes.
The Grid size and Hybrid MPI frames were rendered and checked.

### The experiment reruns

`exp_a` (656960) hit its 14 minute limit at the same point as the first
run, after the same 85 rows. `exp_float` (656963) also ended on the limit,
after its CSV was complete. `exp_b` (656961) and `exp_c` (656962) finished.

| Sweep | First run | Rerun | Largest change |
|---|---|---|---|
| exp_a | 651300 | 656960 | E5 `guided`, 4831 to 4388 (-9.2%), E4 persistent at 512 x 128, 1424 to 1293 (-9.2%) |
| exp_b | 651301 | 656961 | E6 24 threads on 4096 x 1024, up to -19.9% |
| exp_c | 651474 | 656962 | 0.7% |
| exp_float | 652647 | 656963 | 2.7% |
| optlevels | 652696 | 656964 | 1.4% |
| mpi | 651303 | 656965 | 2 nodes, 1 rank x 24, 2048 x 512, 5617 to 6358 (+13.2%) |
| denormal | 651294 | 656966 | omp 24 with flushing, +3.7% |
| peak | 651302 | 656967 | 24 thread AVX2 peak 951 to 937 GFLOP/s (-1.5%), triad within 1.9% |
| profile | 652697 | 656968 | same top function for every rung, serial update_cells_scalar 99.4% to 99.7% |

All 14 rerun jobs have now finished.

The E6 rows on 4096 x 1024 at 24 threads range from 2699 to 5147 within one
job, so that change is noise. The E4 and E5 reps do not overlap between the
two runs, so the size of those effects changes from run to run. Static still
wins E5, and E4 and E8 are still small wins. The kernels did not change, so
the figures keep the first runs, and `docs/RESULTS.md` M2.6,
`HPC_JUSTIFICATIONS.md` and `CODE_GUIDE.md` §9d.3 now quote the rerun next
to them. The RESULTS line saying E4 and E8 matter most on small grids was
removed, because the rerun does not show it.

In the mpi rerun only the 2 node, 1 rank x 24 layout on 2048 x 512 moved,
and its reps do not overlap with the first run. The 2 x 12 layouts and
the two node gains are within 1.4%, so F6 and F8 keep job 651303 and
RESULTS M2.9 quotes the rerun.

## Phase 10. Experiment E2, alignment

### What changed

E2 is now built. The new solver `OmpAlignSolver` (`src/solver_align.cpp`) is
the omp solver running on its own two sheets. Each sheet starts part way into
an oversized `std::vector`, so the experiment picks where column 1 of every row
sits in its 64 byte cache line. `--align-shift K` puts it K doubles past the
boundary. The S and S² tables move with it. It needs `ns` to be a multiple of
8, so every row lands at the same place in its cache line. The vectors are
still big `std::vector`s from mmap, so the page regime matches the Grid's.

- `omp_align_loadu` is the shipped AVX2 row, unchanged. Only the address moves.
- `omp_align_load` is `update_row_avx_load`. It uses `_mm256_load_pd` for
  centre, north, south and the tables, and `_mm256_store_pd` for the store,
  so `--align-shift` must be 0 or 4.
- `omp_align_shuffle` is `update_row_avx_shuffle`. It does one aligned load
  per row per 4 cells and builds east and west from neighbouring blocks with
  `vpermpd` and `vblendpd`. The blocks before and after can read up to 3
  cells past the grid, so the solver keeps padding there.

To share code with the shipped rungs, `update_cells_scalar`, `update_row_avx`,
`update_zero_row`, `apply_stock_edges` and `copy_top_row` now have pointer
versions. The vector versions are one-line calls into them. No arithmetic
moved.

### Decision 1. The shipped rungs are untouched

The Grid keeps its plain `std::vector`. Every rung's maths and every existing
measurement stand. E2 only has to say whether alignment is worth adopting.

### Decision 2. Sweep the shift instead of comparing aligned against unaligned

The windows at stock_i - 1, stock_i and stock_i + 1 are 8 bytes apart. At
most one of them can be 32 byte aligned, so there are four cases, not two.
Shifts 0 and 4 align centre, north, south and the store. Shifts 1 and 5 align
west, shifts 3 and 7 align east, and shifts 2 and 6 align nothing. Every 32
byte load that is not aligned crosses a cache line on every other pass.
Shifts 4 to 7 repeat 0 to 3, so they show the noise.

### Validation

`test_ladder_matches` runs `omp_align_loadu` at every shift, and `load` and
`shuffle` at shifts 0 and 4. All of them on the call and the put must match
simd to the bit. Local run on WSL with g++ 13: all OK at 7 threads, and
`test_avx_row` passes. `make` itself was not run locally (no mpicxx here).
The build used g++ with the solver_mpi object stubbed out.

The local timing on a Ryzen 7 5700X at 1 thread, 2048 x 512, best of 3, is a
smoke test and not a result. Shifts 0 and 4 gave about 1220 Mcell/s. 1, 3, 5
and 7 gave about 1115. 2 and 6 gave about 1040. `load` was the same as
`loadu` at shift 0, and shuffle was about 760. The shipped Grid sat at byte
24 (shift 3) with the mmap threshold pinned, because glibc returns mmap
chunks 16 bytes past a page. simd matched shift 3.

g++ folds many of the loads into the arithmetic instruction (for example
`vaddpd -0x8(%rdx,%rax,1)`). VEX memory operands never check alignment, so
`_mm256_load_pd` only becomes `vmovapd` where it is not folded. The rangpur
job saves the g++ 8.5 disassembly as `results/e2_disasm_<job>.txt`.

### The rangpur run (jobs 657421 one thread, 657422 24 threads, r730-2)

The cosc3500 QOS caps a job at 15 minutes, so the script takes `one` or
`many`. All tests passed in both jobs, and every row prices the same. On
rangpur the shipped Grid has column 1 at byte 24 (shift 3).

One thread, 2048 x 512, median Mcell/s:

| Variant | Mcell/s |
|---|---|
| simd (start, end) | 387.3, 388.2 |
| loadu shift 0, 4 (centre aligned) | 367.4, 368.7 |
| loadu shift 1, 5 (west aligned) | 372.1, 370.6 |
| loadu shift 2, 6 (none aligned) | 382.6, 382.8 |
| loadu shift 3, 7 (east aligned, shipped) | 383.0, 382.0 |
| load shift 0, 4 | 425.1, 427.2 |
| shuffle shift 0, 4 | 461.4, 461.2 |

At 24 threads on 2048 x 512, omp gave 5744, loadu 5508 to 5734, load 6062
and shuffle 6161. On 4096 x 1024 the spread between reps is 2557 to 4369,
which is wider than any gap between variants.

### Decision 3. The cause is g++ 8.5 splitting unaligned loads, not alignment

The rangpur result is the reverse of the local one, so I read the
disassembly (`results/e2_disasm_657421.txt`). With `-mavx2` and no
`-mtune`, g++ 8.5 tunes for a generic x86. That turns every
`_mm256_loadu_pd` into `vmovupd xmm` plus `vinsertf128`, and the store
into `vmovupd xmm` plus `vextractf128`. The shipped row has 11 split loads.
`load` has 6, only the east and west windows. Shuffle has none. An aligned
load is never split. So load and shuffle win because they avoid the split,
and loadu at the aligned shifts gains nothing. g++ 13 (local) does not split,
which is why the local run showed the plain cache line effect instead.
The CPU is a Xeon E5-2670 v3 (Haswell).

E2b reruns the sweep with
`-mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store` on the
AVX2 file only, through a new opt-in Makefile variable `AVXEXTRA`. The
default build is unchanged. It answers whether the shipped row recovers with
no alignment work at all.

### E2b result (jobs 657498 one thread, 657499 24 threads, r730-2)

All tests passed. The disassembly (`results/e2b_disasm_657498.txt`) has no
`vinsertf128` or `vextractf128` left in any of the three rows.

| One thread, 2048 x 512 | E2 (split on) | E2b (split off) |
|---|---|---|
| simd | 387 | 460 |
| loadu centre aligned (0, 4) | 368 | 477 |
| loadu west or east aligned (1, 3, 5, 7) | 371 to 383 | 438 to 440 |
| loadu nothing aligned (2, 6) | 383 | 412 |
| load (0, 4) | 425 to 427 | 477 |
| shuffle (0, 4) | 461 | 460 to 462 |

24 threads on 2048 x 512: omp 6431 (E2 5744), loadu centre aligned 6579,
load 6489, shuffle 6428, nothing aligned 6025. 4096 x 1024: omp 3892,
loadu centre aligned 4026, shuffle 2916, with reps spread from 2683 to 4517.

With the split off the cache line effect is the one the local g++ 13 run
showed. `load` equals `loadu` at the same address. simd (460) runs faster than
loadu at its own shift 3 (440). That is probably because the shipped S and S²
tables are small heap vectors with their own alignment, while the experiment
moves them with the rows. This has not been checked.

### Decision 4. The shipped build is not changed yet

Adding the flag (or `-mtune=haswell`) to the shipped AVX2 build would make
simd about 1.19x and omp about 1.12x faster. Every M2 number, figure and the
deck were measured with the split on, though, so adopting it means rerunning
the sweeps. That is the user's call. RESULTS M2.6 reports E2 and E2b as
experiments and says the rungs stand as measured.

## Phase 11. Experiments E1 and E9, and the slow MPI test

### What changed

E1 rebuilds the whole program under other compiler flags, to ask whether the
compiler could have written the AVX2 rung by itself. The earlier `-O3` bars
(M2.11) used only two-lane SSE2, because no `-mavx2` was given, so they did
not answer that. `slurm/bench_m2_exp_compiler.sh` takes `vector` or
`fastmath`, because the two halves do not fit in one 15 minute job. Each
build times the baseline, serial and simd rungs on one thread at
2048 x 512. The compiler's vectorisation report and an instruction count of
the three hot functions go to `results/e1_vectorised_<part>_<job>.txt`. The
new `--digits D` flag prints the price and the Greeks to D significant
digits. The script prints all four to 17 digits on `config/smoke.cfg` for
every build and rung, which is enough to tell any two doubles apart
(`results/e1_answers_<part>_<job>.csv`).

E9 adds the `mpi_overlap` solver (`MpiOverlapSolver` in
`src/solver_mpi.cpp`). In each step the threads holding the two rows the
neighbours need update those rows first. After a barrier the master thread
applies their stock edges and posts the non-blocking sends and receives,
while the other threads carry on with the rest of their rows. The rows that
were not sent get their boundaries after the parallel region, and only then
does the rank wait for the messages. The master thread calls MPI inside a
parallel region, which MPI allows only from thread level FUNNELED up, so
`main` and the ladder test now start MPI with `MPI_Init_thread` at that
level. The solver refuses to run below it. `test_ladder_matches` holds
`mpi_overlap` to the baseline tolerance and to simd bit for bit, like `mpi`.
`slurm/bench_m2_exp_overlap.sh` runs `mpi` and `mpi_overlap` back to back in
the four layouts of the hybrid sweep (M2.9).

### The first jobs timed out (jobs 657626, 657627, 657628)

All three jobs ran out of their 14 minutes before a single benchmark
finished. The two E1 jobs passed every test and were killed during their
first timed run. The E9 job was killed inside `make test-mpi`. Each CSV
header is written the moment validation ends, and it shows that
`build_and_validate` alone took 13 min 54 s and 13 min 57 s in the E1 jobs.
In the E2 jobs a few hours earlier it took 8 min 20 s to 8 min 27 s.

A diagnostic job (659141) timed one solve of the 128 x 32 test grid with
5200 steps, run as 3 processes of 2 threads, under the environment the bench
scripts export and under two alternatives.

| Solver | `close` binding on cores (as exported) | `OMP_WAIT_POLICY=passive` | no binding |
|---|---|---|---|
| simd | 1.05 s | 1.04 s | 0.62 s |
| omp | 1.04 s | 1.12 s | 0.59 s |
| mpi | 14.24 s | 1.16 s | 0.58 s |
| mpi_overlap | 14.66 s | 1.17 s | 0.60 s |

The whole 3 process ladder test took 124.4 s with the binding and 4.9 s
without it. On 1 process with 7 threads every solver took 0.57 s, so
`make test` was never the problem.

`make test-mpi` passes `--bind-to none` so that each process's threads can
spread out (Phase 4, Change 1). But `slurm/bench_common.sh` exports
`OMP_PROC_BIND=close` and `OMP_PLACES=cores` before it validates, and with
those every process pins its two threads to the first two cores of the node.
All the threads of all the processes then share cores 0 and 1. A thread
that waits does not sleep. OpenMP threads spin at the end of a parallel
region, and MPI spins while it waits for a message. One rank's spinning
threads hold the cores its neighbour needs to finish the update that sends
the message, so every step waits for the scheduler to swap them. Letting
OpenMP threads sleep while they wait (`OMP_WAIT_POLICY=passive`) removes
nearly all of it, so the spinning OpenMP threads are most of the cost. The
solvers that never exchange messages only pay for sharing the cores. Each
test run solves with `mpi` four times (the
call and the put, against the baseline and against simd). Before E9,
`make test-mpi` took about six of the eight validation minutes, and the four
`mpi_overlap` solves doubled it.

### Change 1. `make test-mpi` sets `OMP_PROC_BIND=false`

The test line now sets `OMP_PROC_BIND=false` next to `OMP_NUM_THREADS=2`,
which overrides any binding the caller exported and makes the
`--bind-to none` intent hold inside the bench jobs. A second diagnostic job
(659152) timed `build_and_validate` with the change. The build took 9 s,
`make test` 127 s, `make test-mpi` 10 s and `make avx-check` 3 s, so
validation now takes about 2.5 minutes instead of 14. Only the test target
changed, so no benchmark is affected.

### Validation

The three jobs were resubmitted as 659153 (E1 `vector`), 659154 (E1
`fastmath`) and 659155 (E9), and finished in 6 min 23 s, 6 min 25 s and
5 min 55 s. Each passed `make test`, `make test-mpi` and `make avx-check`
before it benchmarked. In E9 both solvers gave the same price in every
layout.

### E1 result (jobs 659153 and 659154, r730-2, one thread, 2048 x 512)

Median of five reps in Mcell-updates/s. The last column is the largest
relative change of the price, delta, gamma or vega from the `-O2` build of
the same rung, over all three rungs. Each job ran its own `-O2` control,
and the two agree within 3% (simd 370 and 381).

| Build | baseline | serial | simd | largest change |
|---|---|---|---|---|
| `-O2` (shipped) | 57.6 | 158.4 | 370.2 | 0 |
| `-O3` | 112.5 | 241.5 | 377.0 | 0 |
| `-O3 -mavx2 -mfma` | 107.4 | 348.7 | 381.1 | 1.1e-14 |
| `-O3 -march=native` | 106.0 | 265.5 | 455.4 | 1.1e-14 |
| `-O2 -ffast-math` | 58.4 | 177.2 | 382.5 | 1.1e-14 |
| `-O3 -ffast-math` | 189.6 | 271.2 | 386.8 | 1.1e-14 |
| `-O3 -mavx2 -mfma -ffast-math` | 381.0 | 482.2 | 378.4 | 1.1e-14 |

`-O3` repeats M2.11 (112.5 and 238.6 there). Given AVX2 and FMA, the
compiler vectorises the scalar row four lanes wide and reaches 349, which is
0.94 times the hand-written row at the shipped flags. Its report says every
loop was "versioned for vectorization because of possible aliasing". The
compiler cannot prove the input and output rows do not overlap, so it emits
a run-time check and keeps a scalar copy of the loop for the case where they
do. The baseline gains nothing from AVX2 (112.5 to 107.4) because it is
bound by its divisions, and Haswell executes a 256 bit divide as two 128 bit
halves.

`-ffast-math` lets the compiler treat floating point arithmetic as if it
were exact. It may then reorder sums and replace a division by a
multiplication with the reciprocal. The baseline divides by the grid
spacings inside its inner loop. The disassembly of `src/solver_baseline.cpp`
built by rangpur's g++ 8.5 shows 10 divide instructions in the vectorised
inner loop under `-O3 -mavx2 -mfma` and none once `-ffast-math` is added,
because the reciprocals are computed once before the loop. That is the whole
jump from 107 to 381, which matches the hand-written rung. The serial rung
already multiplies by precomputed reciprocals, so its gain from 349 to 482
comes from regrouping the arithmetic. Its vectorised loop shrinks from 51
instructions to 46. The simd rung does not move under any flag, because its
intrinsics fix the instructions.

`-march=native` tunes for Haswell, which stops g++ splitting unaligned
256 bit loads (Phase 10). The hand-written row gains as it did in E2b (455
against 460). The auto-vectorised serial rung loses (265 against 349), even
though, built with `-march=haswell`, its loop is still four lanes wide, has
no split loads and is shorter (34 to 45 instructions against 51). Why it is
slower has not been found. A trap for anyone repeating the disassembly is
that `-march=native` on the rangpur login node targets its AMD EPYC 7542,
for which g++ 8.5 prefers 128 bit vectors. The check above therefore used
`-march=haswell` and not `native`.

### E9 result (job 659155, r730-2 and r730-3)

Median of five reps in Mcell-updates/s. The share is the exchange's part of
a step, the time to post plus the time spent waiting.

| Layout | grid | mpi | mpi_overlap | gain | exchange share |
|---|---|---|---|---|---|
| 1 node, 1 rank x 24 | 2048 x 512 | 5731 | 5482 | 0.96x | 0% |
| 1 node, 2 ranks x 12 | 2048 x 512 | 6069 | 6121 | 1.01x | 19% to 20% |
| 2 nodes, 1 rank x 24 each | 2048 x 512 | 5567 | 8135 | 1.46x | 38% to 15% |
| 2 nodes, 2 ranks x 12 each | 2048 x 512 | 8672 | 8309 | 0.96x | 37% to 35% |
| 1 node, 1 rank x 24 | 4096 x 1024 | 2821 | 2900 | 1.03x | 0% |
| 1 node, 2 ranks x 12 | 4096 x 1024 | 6298 | 6114 | 0.97x | 10% to 10% |
| 2 nodes, 1 rank x 24 each | 4096 x 1024 | 8519 | 8843 | 1.04x | 28% to 28% |
| 2 nodes, 2 ranks x 12 each | 4096 x 1024 | 9833 | 9868 | 1.00x | 28% to 30% |

The `mpi` column is within 4% of M2.9 in every layout. One layout gains, by
a lot. With one rank of 24 threads on each of two nodes, the small grid runs
1.46 times faster and the exchange falls from 38% of a step to 15%. The reps
do not overlap (5493 to 5801 against 7884 to 8383), and 8135 is above
anything `mpi` reached in that layout in M2.9 or its rerun. Every other
layout moves by less than 5%, and its exchange share stays where it was.
Without neighbours the reordered loop costs 4% on the small grid.

The likely reason is how Open MPI 4.1.1 moves a message, read from
`ompi_info` on rangpur. Between nodes it uses TCP, where a message up to
64 KiB is sent eagerly, so the data leaves as soon as `MPI_Isend` is called.
Between two ranks on one node it uses shared memory (the `vader`
transport), where the eager limit is 4 KiB. A row is 16 KiB on the small
grid and 32 KiB on the big one, so on one node it goes by the rendezvous
protocol. The sender only announces the message, and the data moves once the
receiver is inside an MPI call. This Open MPI has no progress thread, so
nothing moves while the threads compute, and the transfer still happens in
`MPI_Waitall`. The two node 1 x 24 layout is the only one with no pair of
ranks on the same node, which fits. The same layout on the big grid gains
only 1.04 times, which this does not explain. A 32 KiB row is larger than
the 16 KiB a Linux TCP socket buffer starts with (`net.ipv4.tcp_wmem`), so
part of it may wait for `MPI_Waitall` as well. Neither idea has been tested.
Raising `OMPI_MCA_btl_vader_eager_limit` above the row size, or turning on
`btl_tcp_progress_thread`, would test them.

### Decision 5. Nothing shipped changes

The shipped build keeps `-O2` and the hand-written rungs. `-ffast-math`
gives up the IEEE rules for the whole program, not one loop. It also lets
the compiler assume no value is ever NaN or infinite, and the tests catch an
unstable run by exactly such a value. The fast-math builds were checked on
one contract, not by the test suite. Whether to ship any of it is the
user's call. `mpi_overlap` stays an experiment, because the best layouts
(2 ranks x 12 per node) do not gain from it.

## Phase 12. Plainer words for the edge row exchange

The user asked for the textbook word for the rows a rank borrows from its
neighbours to be dropped from every document, in favour of direct terms. The
docs now use the names the code already had since Phase 9i.

| Thing | Word in the docs | Name in the code |
|---|---|---|
| sending each rank's first and last rows to its neighbours | edge row exchange | `exchange_edge_rows` |
| the neighbour's rows a rank receives and keeps a copy of | copied rows | `row_below`, `row_above` |
| the exchange's part of a step | exchange share | `exchange_seconds` |

The change covers PLAN, STUDY_GUIDE, STUDY_CARDS, CODE_GUIDE,
DESIGN_DECISIONS, HPC_JUSTIFICATIONS, COMMANDS, RESULTS, HANDOFF,
InitialSpecDesign, this changelog, `presentation/M2_SCRIPT.md` and
`presentation/deck_m2.tex`. The ASCII diagrams in STUDY_GUIDE §6 and §17.3
keep their column widths. The rename table in Phase 9i now describes the old
names instead of spelling them out. `presentation/deck_m2.pdf` still has the
old wording on two frames, because tectonic is not installed on the WSL
machine. `make slides-m2` on the Mac rebuilds it. The local xelatex fails on
the unchanged deck as well (line 249), so it cannot stand in for tectonic.

The code was not touched. The comment above the received rows in
`src/solver_mpi.cpp` still names the textbook term once, and
`scripts/m2_figures.py` uses it in two variable names. No figure shows it.

## Phase 13. The shipped build keeps unaligned loads whole, and every sweep reruns

### The question

E2b (Phase 10, jobs 657498 and 657499) showed that g++ 8.5, given only
`-mavx2`, tunes for a generic x86 and splits every unaligned 256 bit load in
the AVX2 row into a 128 bit load and a `vinsertf128`, and every unaligned
store into two 128 bit stores. Turning the split off made simd 1.19 times
faster (387 to 460 Mcell-updates/s) and omp 1.12 times faster (5744 to
6431), with identical answers. Every sweep until then, including the
Phase 9k reruns (656955 to 656968) and the E1 and E9 jobs (659153 to
659155), was measured with the split on. The user decided to ship the two
flags and rerun the sweeps, so that the deck shows numbers from the final
build.

### What changed

The Makefile now gives the AVX2 file
`-mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store`
whenever the compiler is GCC. The flags sit in a new variable `NOSPLIT`,
which the existing clang test sets in its GCC branch, because they are GCC
options. `AVXFLAGS` became a simply expanded variable (`:=`).
`make -p` prints a recursive variable without expanding it, so the CSV
headers had recorded `$(AVXEXTRA)` literally. They now record every flag in
full. `AVXEXTRA` still comes last and stays the opt-in hook for experiments.
The `test-mpi` line with `OMP_PROC_BIND=false` is unchanged.

`slurm/bench_m2_exp_align.sh` (E2) now exports the positive forms
`-mavx256-split-unaligned-load -mavx256-split-unaligned-store` through
`AVXEXTRA`. They come after the shipped flags and win, so a rerun of E2
still measures the build that splits its loads. A local g++ 13 build
confirmed the order rule. The positive forms placed last bring back exactly
the 11 `vinsertf128` and the one `vextractf128` that E2 found in the shipped
row. `slurm/bench_m2_exp_nosplit.sh` (E2b) was left alone. It only repeats
the shipped flags, so it still measures the same build.

`slurm/bench_m2_exp_compiler.sh` (E1) gained one build in each part,
`O3avx2nosplit` and `O3avx2fastnosplit`, which append the two flags to the
whole-program `-mavx2` builds. Under generic tuning g++ also splits the
loads in the code it vectorises itself, so the comparison with the
hand-written row is fair only when those builds keep their loads whole as
well. The original builds stay, because they show what the compiler does by
default. The instruction count in `results/e1_vectorised_<part>_<job>.txt`
now also counts `vinsertf128` and `vextractf128` in each function.

`make avx-check` still compiles its own copy of the AVX2 file with plain
`-mavx2 -mfma`. It compares the arithmetic of the row with the scalar row,
and the split does not change the arithmetic. `make test` runs the shipped
object through `test_ladder_matches`.

### Validation

The AVX2 file compiled locally under g++ 13 with the two flags, and
`make -n` showed them on the `stencil_avx.o` line. On the rangpur login node
the exact Makefile command under g++ 8.5 gave an `update_row_avx` with no
`vinsertf128` and no `vextractf128`. That check only compiled, and nothing
was timed there. The tree was synced once, without `--delete` and with
`.venv` excluded, and the remote `.venv` survived.

Every job below passed `make test`, `make test-mpi` and `make avx-check`
before it timed anything, and validation took about 2.5 minutes in each, as
Phase 11 predicted. Every CSV header reads
`avx_flags=-mavx2 -mfma -mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store`,
and every job log shows the two flags on the `stencil_avx.o` line. Every
price, delta, gamma and vega matches the earlier run of the same sweep to
every printed digit, and the E1 answers move from the `-O2` build by the
same 1.1e-14 at most as before. exp_a (659246) finished in 9 minutes with
all 60 rows. Its earlier runs had hit the time limit when validation took 8
minutes. No job failed or timed out.

The E1 jobs checked the machine code built on the compute node. In all ten
builds `update_row_avx` has 48 ymm instructions and no `vinsertf128` or
`vextractf128`, where job 659153 had found 57 ymm and 12 packed xmm
instructions with the split on. The compiler's own AVX2 loop in
`update_cells_scalar` has 33 `vinsertf128` and 3 `vextractf128` in
`O3avx2` and none in `O3avx2nosplit`.

| Sweep | Split on | Split off | Elapsed |
|---|---|---|---|
| ladder | 651295 | 659238 | 4:53 |
| strong omp | 651296 | 659239 | 5:19 |
| weak | 651298 | 659240 | 5:05 |
| size | 656959 | 659241 | 6:48 |
| mpi | 651303 | 659242 | 4:20 |
| exp_overlap (E9) | 659155 | 659243 | 5:44 |
| exp_b (E6, E7) | 651301 | 659244 | 7:40 |
| exp_c (E11) | 651474 | 659245 | 4:56 |
| exp_a (E4, E5, E8) | 651300 | 659246 | 9:00 |
| exp_float (E3) | 652647 | 659247 | 8:05 |
| E1 vector | 659153 | 659248 | 7:02 |
| E1 fastmath | 659154 | 659249 | 6:45 |
| denormal (E12) | 651294 | 659436 | 7:26 |
| size with the split turned back on, a control | | 659437 | 7:04 |

The omp_scalar sweep (651297), the optimisation levels (652696) and the
peak (651302) never use the AVX2 row, so they were not rerun. The profile
(652697) was not rerun either, so its simd row still comes from the build
with the split on.

### What was measured

Median Mcell-updates/s unless a row says otherwise. The split-on column
holds the numbers RESULTS quoted before this phase.

| Measurement | Split on | Split off |
|---|---|---|
| ladder baseline, serial | 57, 159 | 57, 159 |
| ladder simd | 379 | 456 (1.20x) |
| ladder omp, 24 threads | 5729 | 6367 (1.11x) |
| ladder mpi, 2 ranks x 12 on one node | 6019 | 6984 (1.16x) |
| strong omp on 1, 12 and 24 threads | 363, 3781, 5733 | 454, 4354, 6311 |
| strong speedup on 12 and 24 threads | 10.4, 15.8 | 9.6, 13.9 |
| Amdahl serial fraction of the AVX2 row | 1.18% | 2.42% |
| boundary share of a 24 thread step | 8.2% | 9.1% |
| weak step time on 1 and 24 threads | 0.142, 0.179 ms | 0.112, 0.166 ms |
| weak scaled speedup on 24 threads | 19.1 | 16.3 |
| size simd, 512 x 128 up to 8192 x 1024 | 327, 377, 382, 171, 334, 323 | 422, 451, 454, 169, 392, 348 |
| size omp on 24 threads, same grids | 1250, 4104, 5711, 6499, 3141, 1709 | 1228, 4264, 6278, 3288, 2990, 2259 |
| roofline omp at 2048 x 512 and 4096 x 1024 | 177, 97 GFLOP/s | 195, 93 GFLOP/s |
| roofline simd at the same grids | 11.8, 10.4 GFLOP/s | 14.1, 12.1 GFLOP/s |
| mpi 1 node 1 x 24, small and big grid | 5714, 2926 | 6233, 2907 |
| mpi 1 node 2 x 12, with exchange share | 6232 (18%), 6314 (10%) | 6843 (20%), 6624 (17%) |
| mpi 2 nodes 1 x 24 | 5617 (36%), 8483 (29%) | 7064 (28%), 9002 (30%) |
| mpi 2 nodes 2 x 12 | 8769 (37%), 9900 (28%) | 8310 (40%), 10863 (28%) |
| one rank per socket over one rank of 24 | 1.09x, 2.16x | 1.10x, 2.28x |
| best two node gain over the best one node layout | 1.41x (70%), 1.57x (78%) | 1.21x (61%), 1.64x (82%) |
| E3 float over double, 1 thread, 24 on 2048, 24 on 4096 | 1.92x, 1.75x, 4.46x | 1.92x, 1.73x, 4.82x |
| E3 reference solve, double and float | 170.2 s, 95.6 s | 149.9 s, 88.0 s |
| E4 persistent region, 2048 x 512 and 512 x 128 | 1.035x, 1.155x | 1.02x, 1.13x |
| E5 dynamic 1, dynamic 8 and guided | 0.60x, 0.72x, 0.84x | 0.57x, 0.63x, 0.78x |
| E8 boundaries in parallel, 2048 x 512 and 512 x 128 | 1.05x, 1.07x | 1.06x, 1.10x |
| E8 serial share of a step | 8.2% to 1.9% | 9.4% to 2.3% |
| E6 12 threads spread against close, 4096 x 1024 | 1880, 3319 | 1894, 3413 |
| E6 unbound threads, 2048 x 512 and 4096 x 1024 | 1.165x, 1.141x | 1.11x, 1.13x |
| E7 interleaved pages, 2048 x 512 and 4096 x 1024 | 1.065x, 1.28x | 1.10x, 1.74x |
| E11 default allocator, omp and serial | 0.68x, 0.56x | 0.63x, 0.56x |
| E12 flushing gain, baseline, serial, simd, omp | 1.08x, 1.29x, 1.38x, 1.26x | 1.07x, 1.28x, 1.42x, 1.33x |
| E1 `-O2` baseline, serial, simd | 57.6, 158.4, 370.2 | 57.6, 157.0, 449.6 |
| E1 `-O3 -mavx2 -mfma` serial | 348.7 | 345.6, and 362.7 without split loads |
| E1 `-O3 -march=native` serial, simd | 265.5, 455.4 | 266.3, 455.8 |
| E1 `-O3 -mavx2 -mfma -ffast-math` baseline, serial, simd | 381.0, 482.2, 378.4 | 380.1, 468.2, 478.0 |
| the same build without split loads | | 376.8, 399.3, 478.2 |
| E9 on 2 nodes 1 x 24, 2048 x 512, mpi and mpi_overlap | 5567, 8135 (1.46x) | 5693, 8009 (1.41x) |
| E9 exchange share in that layout | 38% to 15% | 39% to 19% |
| E9 every other layout | within 5% | within 4% |

M2.9 had printed 38%, 38% and 28% for three of the old exchange shares.
The calculation F8 draws from gives 36%, 37% and 29%, and the table above
and RESULTS now use it throughout.

### What it shows

The flags did what E2b promised. simd gained 1.08 to 1.29 times on every
grid except 2896 x 724, and 1.18 to 1.25 times on the benchmark grid across
the jobs. omp and mpi on one node gained 1.09 to 1.16 times on the
benchmark grid. Rungs that do not use the AVX2 row did not move.

A faster row scales worse. The serial boundaries and the fork and barrier
cost the same per step, while each step got shorter. The Amdahl fraction of
the AVX2 row doubled to 2.42%, the speedup on 24 threads fell from 15.8 to
13.9 and the weak scaled speedup from 19.1 to 16.3. PLAN §5 hypothesis 2 is
weaker, at 80% efficiency on 12 threads against 87%, and RESULTS M2.10 now
says it mostly held.

PLAN §5 hypothesis 6 moved. On the small grid the best two node layout now
gains 1.21 times (61%), down from 1.41 times (70%), because 2 ranks x 12 on
each node did not speed up while every one node layout did. On the big grid
it gains 1.64 times (82%), up from 1.57 times (78%), which clears the 80%
line. The small grid's two node layouts also vary between jobs. The `mpi`
column of the E9 job, run minutes after the M2.9 job with the same build,
gave 9364 for 2 ranks x 12 on each node and 5693 for 1 rank x 24 on each
node, against 8310 and 7064. With the E9 numbers the gain would be 1.33
times (67%). RESULTS M2.9 quotes both.

E1's answer changed. With the hand-written row at 450, the compiler's own
AVX2 code reaches 0.77 times it, or 0.81 times with its own loads kept
whole, where it had reached 0.94 times the split row. With `-ffast-math` the
naive baseline reaches 0.85 times the hand-written row, where before it
matched it, and the serial row reaches 1.04 times, down from 1.27 times. Two
results are new and unexplained. The fast-math serial row slows from 468 to
399 when its loads are kept whole, while the plain AVX2 build gains from 346
to 363. And the hand-written row gains 5% to 6% under `-ffast-math` (449 to
between 473 and 478), where its machine code loses three ymm instructions
(48 to 45). The
Phase 11 remark that the simd rung does not move under any flag held only
while its loads were split.

omp on 24 threads at 2896 x 724 halved, from 6499 to 3288, with all five
reps between 3281 and 3328. The question was whether the flags caused it.
The control job 659437 reran the whole size sweep on the same day with
`AVXEXTRA` set to the positive forms, which turns the split back on. It gave
3267 at 2896 x 724 and matched the earlier split-on runs everywhere else
(omp 5712 against 5711 and simd 383 against 382 at 2048 x 512). So the drop
does not come from the flags. It is the same kind of change that hit serial
at this grid between jobs 651299 and 656959 (Phase 9k), on the grid whose
two buffers sit exactly 2^24 bytes apart (RESULTS M2.8), and its cause was
not found. RESULTS M2.4 and M2.8 now say so. F4 still draws this omp point on
the line, because `scripts/m2_figures.py` treats only serial and simd as
slow at that grid.

E7 on the big grid now reads 1.74 times, up from 1.28 times. The interleaved
runs sat at about 5050 in all three runs. The default runs land either near
2800 or near 4000 in every run, and this time their median fell on the low
side. The other experiments kept their conclusions. Static scheduling still
wins E5, E4 and E8 are still small wins, and E9 still gains only with one
rank on each of two nodes on the small grid. The same layout on the big
grid now gains 1.01 times, down from 1.04 times.

The faster simd row also sits closer to what one core can stream. At 454
Mcell-updates/s it moves about 10.9 GB/s against the 12.4 GB/s one-thread
triad, and it loses 14% at 4096 x 1024 and 23% at 8192 x 1024. RESULTS M2.4
now says that instead of saying simd barely cares about size.

### Documents

`docs/RESULTS.md` has a note at the top of the Milestone 2 section about the
split, and every number in M2.1 to M2.11 that changed. E2 and E2b keep their
split-on numbers. `HPC_JUSTIFICATIONS.md` quotes the new numbers, shows the
new Makefile lines, and no longer says E2 was never run. In the Makefile,
`M2_CSVS` points at the new CSVs and `m2-figures` takes its test log from
`slurm-659247.out`. omp_scalar 651297, optlevels 652696, profile 652697 and
peak 651302 stay. `make m2-figures` was run once into a scratch folder to
check the list, and it drew every figure without an error. STUDY_GUIDE,
CODE_GUIDE, STUDY_CARDS, DESIGN_DECISIONS, README, `presentation/` and
`scripts/` were left for the sessions that own them, and several of their
numbers are now out of date.

## Phase 14. The Milestone 2 deck as Quarto slides

### What changed

`presentation/deck_m2.qmd` is a Quarto revealjs version of the beamer deck
`presentation/deck_m2.tex`, with its theme in `presentation/deck_m2.scss`.
`quarto render presentation/deck_m2.qmd` writes one self-contained file,
`presentation/deck_m2.html`, at 1600 by 900. The recorded slides keep the
beamer titles, subtitles, callouts and source lines, and each one carries
its take from `presentation/M2_SCRIPT.md` as speaker notes (press S). The
callout band still stops at 74 percent of the width, so the bottom-right
corner stays empty for the webcam. Every benchmark number is defined once in
the front matter under `n`, from `docs/RESULTS.md`, the same way the beamer
deck defines its macros.

The code screenshots are replaced by the code itself, with the file and the
real line numbers above each block. The serial, SIMD, OpenMP and MPI code
slides are pairs. The first slide shows the rung below, and the next press
changes it into the new code with reveal.js auto-animate, then a bar grows
from the old rate to the new one. The grid step, the data layout, the SIMD
row, the OpenMP row bands and the MPI edge rows are drawn as SVG and step
forward one press at a time. The ladder and the two rooflines are redrawn
from the RESULTS numbers and draw themselves when their slide opens. The
other charts are the beamer figures, cut from `deck_m2.pdf` into
`presentation/assets_m2/` because `results/m2/` is not in git.

`scripts/deck_m2_svgs.py` draws every SVG and writes it into the deck between
`<!-- svg:name -->` markers, so it can be run again after a change.

The recorded part opens with four manim videos in `presentation/anim_m2/`,
drawn by `presentation/anim_m2/scenes.py`. `model_pde` goes from the two
Heston processes through the expected discounted payoff and Itô's lemma to
the PDE. `pde_code` takes one explicit Euler step with its stability limit,
the central differences, the reciprocals and lines 121 and 122 of
`src/stencil.cpp`. `differences` works one cell of the 421 by 51 test grid
with its real numbers, ending at 0.04959, the value the solver's own
snapshot holds. `full_run` is drawn from snapshots of a local run of the
same test grid with `--dump-every`, 56,000 steps to a price of 195.6355.
Each press plays one section and stops on its last frame, from the section
end times in the json file next to each video. These numbers are prices
from the test grid, so they are not benchmark numbers.

Detail from the backup slides now sits on the main slides as a short note
under a figure or an extra bullet, such as the Amdahl fit on strong scaling,
the eager limits on the E9 slide, the float price on the rejected slide and
the four tests on the Black-Scholes check. The profile and the optimisation
levels share one slide. The problem, the equation table, the compiler
alone, E1 and the step breakdown moved to the start of the appendix.

### Why

The beamer deck shows code as images and every figure as a still. In the
video and the interview the step from one rung's code to the next, and the
rate that step bought, is the thing to explain, so the slides now show that
change happening.

### What to know

The deck has 38 recorded slides against the beamer's 35 frames, so the take
numbers in `M2_SCRIPT.md` follow the Quarto slides now, and the notes under
each slide hold the same words as its take. The takes are timed at 150 words
per minute, a little faster than the 140 before because the animations carry
the working, and the 1339 words come to 9 minutes 34 seconds. The takes on where the time goes, the split loads and placement were cut by about 20 seconds, and the detail they dropped stays on screen, with the compiler's 77 percent moved onto the split loads slide. On the SIMD code slide the take says the
serial loop turns into the AVX2 loop, because the two loops are no longer
shown one above the other.

The reflection slide names each item by its technical term. Every rank
allocates the whole grid, so one node's memory caps the grid size, and the
next step is allocating only the rank's own rows plus one ghost row on each
side. The E9 test is setting `OMPI_MCA_btl_vader_eager_limit` above the
16 KiB row. The MPI rows slide gives the copied rows their textbook name,
ghost rows, once, so the conclusion and the reflection slide can use it.

The stepped videos froze in Safari 26.6 on slide 5, One cell with real
numbers, during the fourth section, the new value. The video clock kept
running but the picture stopped. Every key press set `currentTime` to the
section start and called `play()` straight away, while the video already sat
at that time from the end of the last section. In a test in Safari that drew
the video into a canvas every 50 ms, the fourth section showed 12 different
pictures this way, against 66 when the seek was skipped. The video file and
the Mac were checked too. One core decodes the 1600 by 900 video at about
1760 frames a second, and the video played through without a stall when
nothing seeked.

`playSection` now seeks only when the video is more than 0.05 s from the
section start, and it waits for the `seeked` event before playing when a seek
is running. It no longer sets `currentTime` back to the section end after the
pause, because that was a second seek. The section end check also accepts
`ended`, because the last stop of `differences.mp4`, 17.533 s, is a fraction
past the end of the video at 17.5327 s. Three passes over slide 5 in Safari
gave 66 to 67 pictures in the fourth section, finishing in 4.4 s.

`differences.mp4` also has frame times that drift 10 ticks of 1/15360 s off
the 30 frames a second grid from frame 174 on. It played through cleanly in
Safari, so the file was left as it is.

### The opening as stepped diagrams

The lecturers said the slides do not need to be heavy on notation from the
start. The opening can show the PDE once and rely on diagrams of how a
finite difference is worked out and why each part of the method is needed.
The four manim videos also needed the stepped video player and its Safari
fix. The opening now keeps the order of the videos, from the model to the
PDE, the differences, the C++ and the real numbers, and draws each part as a
stepped SVG or an auto-animate move, the same as the grid step and the code
slides. One key press shows one step.

Slides 3 to 8 now run as follows.

What an option is starts from today's price of 5200 and the strike line at
5250. One path ends at 5700 and the option pays 450. One path ends at 5000,
below the strike, so nobody buys and it pays 0. The payoff line max(S - K, 0)
then joins every ending price to what it pays, eight more paths show that
the path is unknown, and the price is the average payoff discounted back to
today. The paths come from a seeded Euler simulation in
`scripts/deck_m2_svgs.py` with the deck's Heston parameters, picked from 4000
so the two examples end at 5700 and 5000. They are an illustration and not
output of the solver.

From the model to the PDE shows the PDE as a row of six coloured terms. The
presses label the share price terms, the volatility terms, then the
correlation and the discounting, in plain words, and the last press shows
five grids solved one step back at a time from the payoff at expiry.

Finite differences on one cell moves the PDE to the top of the slide with
auto-animate, each term keeping its colour. It works the cell at S = 5250,
v = 0.04 of the test grid at expiry, where west, centre and east hold 0, 0
and 50. The slope (50 - 0) / 100 = 0.5 is drawn as the line from west to
east, the curvature 50 / 2500 = 0.02 as the change from slope 0 to slope 1,
and each press marks the stencil cells the difference reads with their
weights. The variance and cross differences are 0 because every row holds
the same payoff.

From the PDE to C++ puts the same cell beside the lines of
`src/stencil.cpp` that update it. Lines 107 to 111 hold the five
differences, each a multiplication by a reciprocal from lines 28 to 32.
Lines 114 to 118 hold the terms, 11025 and 84 with the others 0. Lines 121
and 122 give change = 11109 and the new value 4.464e-6 x 11109 = 0.04959,
the value the solver's own snapshot holds, 0.04959375. The last press gives
the stability limit, 5.66e-6 on the test grid and 3.66e-7 on the full grid,
so the full grid needs 980,000 steps.

One step over the whole grid is unchanged. A full run draws the v = 0.04 row
of the test grid at expiry, 1 and 9 days before, 37 days before and today,
from the solver snapshots, ending at 195.64 at S = 5200.

`scripts/deck_m2_svgs.py` gained `option_paths`, `backwards_in_time`,
`finite_differences`, `pde_to_code` and `full_run`, and lost the video step
code. The PDE row is HTML in `deck_m2.qmd`, styled by `.terms` in
`deck_m2.scss`, so MathJax draws each term and auto-animate can move it by
its `data-id`. The stepped video player left the title slide and its style
left `deck_m2.scss`. The deck file dropped from 15.5 MB to 6.4 MB because
the videos are no longer embedded. The recorded part has 39 slides, one
more than before, and the takes in `M2_SCRIPT.md` still add up to 9 minutes
34 seconds, 1338 words.

### Interview prep for the new opening, package size, AI declarations

`presentation/INTERVIEW.md` section 4b now lists the 39 recorded slides of
the Quarto deck, with what each press shows on slides 3 to 8. Section 8 has
a new group of questions on the opening, from what a call option is to the
stability limit and the reciprocals, each with its answer and where the full
answer is.

`make package` no longer copies the solver snapshots in
`presentation/anim_m2/snapshots` or the manim videos and their json stop
files, because the deck no longer uses them. The files it zips, before the
video, dropped from 92.6 MB to 45.6 MB, so the zip stays well under the
100 MB limit with the video added. `scenes.py` and `full_run.cfg` stay,
because they hold the commands that make the snapshots.

Every script in `scripts/` now opens with the three-line AI declaration, the
same as the C++ files.

### What a cell update is, on the measurement slide

The measurement slide now says what the rate counts. A cell update is one
new value for one cell in one step, and every rung divides cells times steps
by the seconds its time loop took, `src/solver_serial.cpp` lines 66 and 67.
On the 2048 by 512 grid with 2000 steps that is 2.1 billion updates, so the
serial rung's 159 million a second is about 13.2 s. The take gained one
sentence, and the script now comes to 9 minutes 39 seconds, 1349 words.

The one-line definition read as too dense for someone new to it, so the
slide now opens with five boxes, one per press, that build the serial rung's
rate. 2048 x 512 = 1,048,576 cells, times 2000 steps, is 2,097,152,000 cell
updates, and divided by the 13.2 s of the time loop that is 159 million a
second. `rate_steps` in `scripts/deck_m2_svgs.py` works the boxes out from
the grid size, the step count and the measured rate. The five bullets became
four to make room, with the same facts. The take says the same sum in one
sentence, and the script now comes to 9 minutes 46 seconds, 1368 words.

### A voice pass over the takes

The takes in `M2_SCRIPT.md` and the speaker notes in `deck_m2.qmd` were
reworded in the style of the voice profile, with a fact and its consequence
joined in one sentence where they belong together, "To conclude" dropped
from the conclusion, and the reflection said in the first person with what
comes next. Every number and technical term is unchanged and no
analogies were added. The script now comes to 9 minutes 49 seconds, 1375
words, so a take that runs long when spoken should be trimmed back.

### A fuller reflection slide

The reflection criterion asks for what surprised us, what we would do
differently and how to extend the work, and the slide only covered the
first and the last, so its take read as a list. Slide 38 now has three
blocks. What surprised me keeps the split loads, the 2^24 byte address
effect and the system settings. What I would do differently is new, setting
the floating point mode before the first sweep because the whole first round
was rerun, reading the compiler's assembly for every rung from the start
because the split loads only showed up after simd ran at 387, and giving
every slurm job its own copy of the tree because two first-round jobs were
killed by a sync and a `make clean`. How I would extend it keeps the per-rank
rows, the eager limit test, an implicit or ADI scheme and a CUDA kernel. The
take says each part in my own words.

To keep the video under 10 minutes the measurement take drops the spoken
sum, because the five boxes show it, and the strong scaling take drops the
9 percent boundary share at 24 threads, which stays on the slide. The
script now comes to 9 minutes 45 seconds, 1365 words.

### Denormals defined, the per-rank rows explained, and an AI declaration slide

The reflection slide and its take now say what a denormal is, a double below
2.2e-308, too small for the normal format, on which this CPU spends tens of
extra cycles for every operation that touches one. Short benchmark runs
leave values that small far from the strike, which is why the first round of
sweeps was rerun. The measurement take no longer tells that story, because
the reflection now does, and the amber block on the measurement slide still
shows the numbers.

The extension on the reflection slide now says where it applies and why. On
the mpi rung every rank allocates the whole 2048 by 512 grid, 16 MiB
(`src/solver_mpi.cpp` line 126), and updates only its own rows. The update
reads only the row above and below (`src/stencil.cpp` lines 89 and 90), so a
rank needs just its own rows plus one copied row from each neighbour. With 4
ranks on two r730 nodes that is 128 own rows and 2 copied rows, about 4 MiB,
so the grid could grow with the number of nodes.

Slide 40, the last recorded slide, is the AI declaration, held for three
seconds with nothing said. To pay for the longer reflection take, the ladder
take no longer repeats that every version must match the baseline, which the
same answer slide says, and the profile take no longer reads out the
baseline's 58, which is on the slide. The script comes to 9 minutes 47
seconds, 1367 words, plus the hold on the declaration.

### The video as 19 teaching takes

The takes were too short to explain anything. 34 spoken slides in 9:47 left
about 20 seconds a slide, so a take could read out numbers but not the idea
behind them, and the strong and weak scaling takes in particular did not say
what either kind of scaling measures. The task sheet marks the choice of
topics and the time given to each, and says a shorter talk done well beats a
padded one, so the video now covers fewer topics and teaches each one.

There are now 19 takes over 32 recorded slides. Each take says what the idea
is, shows what was measured, then says why. A take that runs over two or
three slides has an italic line in `presentation/M2_SCRIPT.md` saying when
to press next, and each slide's speaker notes hold the words said over it.

Eight slides moved to the start of the appendix, where they are still
there for the interview with their old takes as notes. These are the
project summary, the Black-Scholes check (now one line in the band of the
full run slide), the ladder (each rung slide grows its own bar, and the
conclusion has the table), the split loads (the reflection says it), the
one-core roofline (the roofline is taught once, on the node roofline), the
two MPI time loop slides (the MPI take keeps the exchange code) and the
overlap experiment E9.

Some slides moved so each take tells one story. The data structure follows
the one step slide, placement follows MPI so the first-touch explanation
leads into one process per socket, and the same answer slide comes just
before the conclusion.

The slides changed in a few places to match the takes.

- The strong and weak scaling subtitles now say in one sentence what each
  measures, and the bands state Amdahl's and Gustafson's laws.
- The OpenMP slide says the grid has 512 rows, one per variance value, and
  that thread 0 updates rows 1 to 22, so 21 rows per thread on the weak
  scaling slide has a meaning.
- The grid size and roofline subtitles say what each figure shows.

The script comes to 9 minutes 46 seconds, 1385 words, plus the hold on the
AI declaration. `presentation/INTERVIEW.md` sections 4b, 4c and 8 and the
README map follow the new order, and the appendix is numbered A1 to A45 with
the anchor for each slide.
