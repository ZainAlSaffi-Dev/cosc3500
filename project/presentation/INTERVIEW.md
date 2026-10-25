# INTERVIEW, the Milestone 2 interview and what the deck has to do for it

This file was written on 2026-10-14. It records what the interview is, what a
friend who sat it said about it, and what the video, the deck and the README must do
because of that. It also lists the C++ in the project most likely to make an
interviewer ask why, with an answer for each, and a tour of the code for the
big screen.

The numbers in this file come from the rerun with unsplit AVX2 loads
(CHANGELOG Phase 13, jobs 659238 to 659249 and 659436), as `docs/RESULTS.md`
quotes them. The deck, the script and the README use the same numbers.

---

## 1. What the interview is

### From the task sheet

The Milestone 2 interview is an in-person, identity-verified hurdle. It takes
about 7 minutes and happens in the second week of exams, which is the week
after the video is due. It is pass or fail, and failing it fails the course.
The task sheet sends you to the course FAQ and the course profile for more
detail. Those have not been checked yet, so the exact date, the room and any
rule about what you may bring are still open.

### From a friend who sat it

- The tone was relaxed. The tutor who marked the Milestone 1 and Milestone 2
  videos is expected to run the interview, so they already have a lot of
  context and will have watched the Milestone 2 video.
- It was in the AEB long room, the room Code Jam runs in. Several tutors sat
  at large TV screens with a keyboard and mouse, and the screen showed the
  student's whole project.
- The questions lean towards the HPC code, because the interview is about the
  optimisation. Where the video shows code, expect to be asked what happens
  there and why.
- You should be able to say in plain English what the project is.
- How deep they go into the C++ depends on how much code changed between
  Milestone 1 and Milestone 2. They may ask you to open another part of the
  code and explain it, depending on how your own explanation goes.
- A project you chose yourself, which this one is, draws more questions about
  the model than a preset project.
- Anything unusual in the C++ makes them curious about why it is there.
- Graphics draw questions about how they were made and what every part of
  them means. The friend was asked about their own visualisations.
- The friend's advice was to make the Milestone 2 video really good, because
  the interview builds on it. A very clear deck with everything in it helped
  them a lot.

---

## 2. What this means for us

| What happens in the room | What we prepare | Where it lives |
|---|---|---|
| They start from the video they watched | Every slide can be explained from what is on it, without outside notes | the deck |
| They ask what the project is | A 30 second plain-English answer, said the same way every time | section 3, deck slide 2 |
| They ask about the optimisation code | Each rung's code on its slide, large and readable, and the same lines easy to find in the editor | the deck, section 6 |
| They ask why something unusual is in the C++ | One short answer per unusual construct | section 5 |
| They ask how a figure was made and what it shows | The script, the data and the meaning of each mark for every figure | section 7, deck appendix |
| They open another file | Know where everything lives | section 6, the README map |
| They push past what the video covered | Backup slides in the same PDF, after the closing slide | deck appendix |

Two practical points follow from the screen showing the submitted project.

- The README is likely the first file opened on that screen. It was
  rewritten for Milestone 2 with the ladder, how to build and run, and a map
  of where everything is.
- `make package` zips `results/` and `presentation/` with the code, so the
  figures and the deck PDF with its appendix are included in the submission.
  The figures in `results/` were lost in the move from the Mac and were
  redrawn on this machine from the rangpur CSVs. `make package` must run
  after `make m2-figures` and `make code-shots`, or the zip ships without
  them.

---

## 3. The project in plain English

This is the answer to "what is your project", to be said in about 30 seconds.

> I price a financial option, which is the right to buy a share at a fixed
> price on a fixed day. In the Heston model the share's volatility moves
> randomly too, so there is no simple formula. The price comes from solving a
> partial differential equation on a grid, with the share price along one axis
> and the variance along the other. At expiry the value of every cell is known,
> so the solver steps backwards in time to today. Each step rewrites every cell
> from its nine neighbours in the previous step, and no cell in a step depends
> on another cell in the same step. I parallelised the solver using that
> independence. AVX2 updates four neighbouring cells with one instruction,
> OpenMP shares the rows of each step between the 24 cores of a node, and MPI
> gives each process a block of rows. Every step, each process sends its first
> and last rows to its neighbours and receives theirs into its copied rows.
> Every version gives the same price as the original. On one node the fastest is
> 123 times faster than the original, and 147 times on two nodes.

The follow-up they are most likely to ask is why the randomness does not
appear in the program. The randomness is averaged out analytically. The
Feynman-Kac theorem turns the expected payoff under the two random processes
into this deterministic equation, and the noise terms become the second
derivative terms. There is no random number generator anywhere in the solver.

### The numbers to know

Every row is on the 2048 x 512 grid on one r730 node unless it says
otherwise. The last column is the section of `docs/RESULTS.md` that holds
the number and its job id.

| What | Number | RESULTS |
|---|---|---|
| baseline, serial, simd, omp on 24 threads, mpi 2 x 12 | 57, 159, 456, 6367 and 6984 million cell updates per second | M2.1 |
| gain of each rung over the one below | 2.8x, 2.9x, 14.0x, 1.10x | M2.1 |
| total over the baseline | 123x on one node, 147x on two | M2.1, M2.9 |
| serial on the roofline | 4.9 GFLOP/s, 96 percent of the 5.12 scalar peak | M2.5 |
| simd on the roofline | 14.1 GFLOP/s, a third of the 43.4 AVX2 FMA peak, just under the one-core memory roof of about 16 | M2.5 |
| unsplit loads, E2b | 1.19x on simd, 1.12x on omp, same answers | M2.6 |
| strong scaling on 24 threads | 13.9x and Amdahl f = 2.4 percent. The scalar row 18.5x and f = 1.6 percent | M2.2 |
| efficiency on 12 threads | 80 percent, 9.6x | M2.2 |
| serial boundaries | 0.44 percent of a step on 1 thread, 9.1 percent on 24 | M2.2 |
| weak scaling | 0.11 to 0.17 ms per step, scaled speedup 16.3 on 24 threads, Gustafson predicts 21.9 | M2.3 |
| grid size | omp 40x serial while the grid fits in the L3 caches, half the rate at 4096 x 1024 | M2.4 |
| node roofline | 195 GFLOP/s at 2048 x 512, 93 at 4096 x 1024, memory roof about 88 | M2.5 |
| one process per socket | 1.10x on the small grid, 2.28x on 4096 x 1024 | M2.9 |
| two nodes | 1.21x one node on the small grid (61 percent), 1.64x on the big one (82 percent). The exchange is 40 percent of a 0.126 ms step | M2.9 |
| E9 overlap | 1.41x with one process on each of two nodes, every other layout within 4 percent | M2.6 |
| E1 | the compiler's AVX2 code reaches 0.77x of the hand-written row, 0.81x with unsplit loads. With fast-math the baseline reaches 380 | M2.6 |
| E3 float | 1.73x at 2048 x 512 and 4.82x at 4096 x 1024, price off by 4 cents | M2.6 |
| E5 schedules | dynamic 1, dynamic 8 and guided at 0.57x, 0.63x and 0.78x of the omp rung | M2.6 |
| E6 and E7 on 4096 x 1024 | 12 threads on one socket 1894, on both 3413. Interleaved pages 1.74x | M2.6 |
| denormals flushed | 1.07x to 1.42x per rung, the same prices | M2.7 |

---

## 4. The deck plan

### Decisions made on 2026-10-14

- One PDF, `presentation/deck_m2.pdf`. The recorded part is made as
  explanatory as possible, with intermediate slides added where a step needs
  one. An appendix after the closing slide holds the detail for the
  interview, and it is never recorded.
- Code appears as larger before-and-after snippets on the slides, with the
  changed lines highlighted. The video does not switch to the editor. The
  editor walk-through in section 6 is rehearsed separately for the interview.
- The numbers come from the rerun with the no-split flags (CHANGELOG Phase
  13). The deck, the script and the README were updated from it on
  2026-10-14.
- Every M2 figure is redrawn on this machine from the CSVs copied down from
  rangpur, with the fixes below.
- The video stays under 10 minutes. The time for new slides comes from
  trimming words from slides that are already well understood.

### 4a. What was wrong with the 25-frame deck of 13 October

| Frame | Problem | Fix |
|---|---|---|
| 2, the problem | The value map still has tiny text, a colon title and an em dash. The payoff title uses a colon and an "X, not Y" phrase. | Redraw both figures with plain titles and larger text. |
| 6, data structure | The code is too small to read, and the three numbers overlap the callout box. | A larger, shorter code shot and a fixed layout. |
| 10, serial rung | Only the new code is shown, so the change itself is not visible. | Before and after, the baseline's divisions next to the serial rung's multiplications. |
| 11, SIMD | Two code panels at a size nobody can read. | A diagram slide of four cells in one register, then a larger before-and-after code slide. |
| 12 and 22, roofline | Overlapping labels and too many points. The message is lost. | Two simple rooflines, one for a single core and one for the whole node. |
| 13, OpenMP | The code shot says 15.1x and the callout says 15.8x. Both are right but measured against different references, which invites a question. | One number with its reference stated, and a diagram of each thread's contiguous block of rows. |
| 15, placement | Four series and two grids in tiny type, and the numbers in the callout do not appear on the chart. | One chart that shows one socket against two and pages on socket 0 against interleaved pages, with the numbers on the bars. |
| 16, hybrid MPI | Code and diagram both tiny, and an analogy word for the copied rows in the subtitle. | A diagram slide of the blocks of rows and copied rows, then a larger code slide. |
| 20, grid size | The labels of the shaded cache regions are tiny and the grey 111x bar needs a paragraph to explain. | Larger labels and the 2896 point marked as the address effect, with its explanation in the appendix. |
| 23, tried and rejected | Sixteen rows of tiny labels with experiment ids. The ladder, the most important chart in the deck, only appears here at the end. | The ladder gets its own slide early on. The tried experiments get a clearer chart of their own. |
| 24, conclusion | The bullets collide with the callout. "Next is overlapping the exchange" is out of date because E9 has run. | Fixed layout, and the E9 result in place of the old next step. |
| none | No plain-English summary, and nothing on E1, E2 and E2b or E9. | New slides for each. |

### 4b. The recorded part, slides 1 to 39

The recorded part is `presentation/deck_m2.html`, the Quarto deck, and
`presentation/M2_SCRIPT.md` has one take per slide. The opening, slides 3 to
8, builds one step per key press, from what an option is to the update in
C++. Every rung that changed the code has a before and after made with
auto-animate, the serial rung (14 and 15), SIMD (17 and 18), OpenMP (22 and
23) and MPI (30 and 31). Slide 38 is the reflection the task sheet asks for,
what surprised me and what I would do next.

| Slide | Title | Slide | Title |
|---|---|---|---|
| 1 | Title | 21 | OpenMP, one block of rows per thread |
| 2 | The project in one slide | 22, 23 | OpenMP code |
| 3 | What an option is | 24 | Strong scaling |
| 4 | From the model to the PDE | 25 | Weak scaling |
| 5 | Finite differences on one cell | 26 | Grid size |
| 6 | From the PDE to C++ | 27 | After threads, memory bandwidth |
| 7 | One step over the whole grid | 28 | Placement |
| 8 | A full run | 29 | MPI, one block of rows per process |
| 9 | Checking the answer | 30, 31 | MPI code, the time loop |
| 10 | The data structure | 32 | MPI code, the edge row exchange |
| 11 | How every number was measured | 33 | Hybrid layouts |
| 12 | Where the time goes | 34 | Overlapping the exchange with the update (E9) |
| 13 | The ladder | 35 | Same answer on every rung |
| 14, 15 | The serial rung | 36 | Tried and rejected |
| 16 | SIMD, four cells per instruction | 37 | Conclusion |
| 17, 18 | SIMD code, before and after | 38 | What surprised me, and what I would do next |
| 19 | The compiler was splitting the loads (E2, E2b) | 39 | Close |
| 20 | SIMD on the roofline | | |

The opening slides and what each step shows.

| Slide | Presses | Steps |
|---|---|---|
| 3 | 5 | a path ending at 5700 pays 450, a path ending at 5000 pays 0, the payoff line max(S - K, 0), eight more paths, the price as the average discounted payoff |
| 4 | 4 | labels for the share price terms, the volatility terms, the correlation and discount terms, then five grids solved back from the payoff |
| 5 | 4 | the slope 0.5, the curvature 0.02, the variance differences, the cross difference, each with the stencil cells it reads |
| 6 | 3 | opens on `src/stencil.cpp` lines 107 to 111, then the terms at lines 114 to 118, the update at lines 121 and 122 giving 0.04959, the stability limit |
| 7 | 6 | the stencil moving along row 1, the rest of the grid, the swap |
| 8 | 3 | the v = 0.04 row 1 and 9 days, 37 days and 91 days before expiry, ending at 195.64 |

### 4c. The appendix, pages 37 to 68

Page 36 is the divider. In the room, open `presentation/deck_m2.pdf` and type
the page number into the viewer to go straight to a backup slide. The
footer of every backup slide shows its A number.

| Page | A | Slide | Page | A | Slide |
|---|---|---|---|---|---|
| 37 | A1 | Where everything is | 53 | A17 | Thread binding on the r730 nodes |
| 38 | A2 | The Heston model | 54 | A18 | Page placement, numactl and huge pages |
| 39 | A3 | The five differences | 55 | A19 | Which rows a process owns |
| 40 | A4 | The boundaries | 56 | A20 | E9 in code |
| 41 | A5 | The timestep limit | 57 | A21 | Why E9 only helped across nodes |
| 42 | A6 | Why an explicit scheme | 58 | A22 | The two node numbers between jobs |
| 43 | A7 | Reading the price and the Greeks | 59 | A23 | MPI started for threads |
| 44 | A8 | Flops and bytes per cell | 60 | A24 | The test that timed out |
| 45 | A9 | The AVX2 row next to the scalar row | 61 | A25 | Denormal numbers |
| 46 | A10 | Why SIMD is below 4x | 62 | A26 | The 2^24 byte address effect |
| 47 | A11 | Where the rows sit in their cache lines | 63 | A27 | Float instead of double |
| 48 | A12 | Experiment E1 in full | 64 | A28 | The six hypotheses |
| 49 | A13 | Why there are no locks | 65 | A29 | What the tests check |
| 50 | A14 | The swap | 66 | A30 | How the figures were made |
| 51 | A15 | Constants copied into local variables | 67 | A31 | Tests before timing |
| 52 | A16 | Amdahl and Gustafson | 68 | A32 | Where a number came from |

---

## 5. The unusual C++, and the answer for each

The code was kept plain on purpose. These are the places where a tutor is
still likely to ask why, roughly in the order they appear on screen.

**1. The AVX2 intrinsics.** `src/stencil_avx.cpp` lines 58 to 109. A
`__m256d` is a 256 bit register that holds four doubles, so every line
operates on four neighbouring cells. `_mm256_loadu_pd` loads four doubles
from any address, `_mm256_fmadd_pd(a, b, c)` computes a times b plus c with
one rounding, and `_mm256_fnmadd_pd(a, b, c)` computes c minus a times b.
The loads are unaligned because the west, centre and east loads start 8
bytes apart, at cells `stock_i - 1`, `stock_i` and `stock_i + 1`, so at most
one of them can start on a 32 byte boundary. Intrinsics were used because
they map one to one onto instructions, so the rung measures exactly what was
written.

**2. The `#if defined(__AVX2__) && defined(__FMA__)` guard.**
`src/stencil_avx.cpp` line 11, with a scalar fallback at the end of the
file. Only this one file is compiled with `-mavx2 -mfma`. The Makefile sets
`AVXFLAGS` at line 37, and only the rule at lines 66 and 67 uses it. That
keeps every other rung compiled exactly as before, so the ladder compares
like with like, and the program still builds and runs on a machine without
AVX2, where `simd` prints a warning and uses the scalar row.

**3. The two flags that stop the split loads.** The AVX2 file also gets
`-mno-avx256-split-unaligned-load` and `-mno-avx256-split-unaligned-store`
(`NOSPLIT`, Makefile line 23). With only `-mavx2`, g++ 8.5 tunes for a
generic x86 and compiles every unaligned 256 bit load as two 128 bit loads
joined by `vinsertf128`, and every unaligned store as two 128 bit stores.
The disassembly in E2 showed it, and E2b measured 1.19x on simd with the
source unchanged. The flags are GCC options, so the Makefile leaves them out
when the compiler is clang. The rerun with them took simd from 379 to 456
and omp from 5729 to 6367, with every answer the same to every printed
digit. `-march=native` would also keep each unaligned load as one 256 bit
load, because it tunes for Haswell, but it would tie the build to the
machine it was compiled on.

**4. The OpenMP pragma.** `src/solver_omp.cpp` line 45.
`schedule(static)` gives each thread one contiguous block of rows, which
suits rows that all cost the same and keeps each thread on the same rows
every step. `default(none)` makes every variable used in the loop be listed,
so forgetting one is a compile error. Without it, that variable would
silently become shared between the threads. The loop counter and everything
declared inside the body are private to each thread automatically.

**5. Why there are no locks.** Each thread reads only `current` and writes
only its own rows of `next`. No two threads ever write the same cell, and
no thread reads a cell of `next` during the loop. The implied barrier at the
end of the parallel loop makes every thread finish its rows before the
boundaries are updated and the two buffers swap pointers.

**6. The two references that stay right after every swap.**
`src/solver_omp.cpp` lines 36 and 37 bind `current` and `next` once, before
the time loop. `Grid::swap_buffers` (`src/grid.cpp` line 32) calls
`std::swap` on the two vectors, which swaps their internal pointers and
copies no numbers. The references name the two member vectors, so after the
swap `current` holds what was just written.

**7. Copying the constants into local variables.** `src/stencil.cpp` lines
66 to 76. The row writes through a `double*`, and the compiler cannot prove
that pointer never points into the constants. Without the copies it would
have to reload every constant from memory after every store. As locals they
stay in registers for the whole row.

**8. Two versions of each row function.** `include/stencil.h`. One takes
`std::vector` references and the other takes raw pointers. The vector version
is a one-line call into the pointer version, so the maths exists once. The
pointer versions were added for E2, whose buffers start part way into a
vector.

**9. The class hierarchy and `make_solver`.** `include/solver.h` and
`src/solver.cpp` line 95. Every rung is a subclass of `Solver` with its own
`solve`. `virtual ... = 0` means every subclass must write it, and
`override` makes the compiler check the signature matches. `make_solver`
turns the `--solver` name into an object, so one binary runs every rung and
the benchmark scripts only change one flag. `std::unique_ptr` deletes the
solver automatically when `run` returns.

**10. `static` functions in `src/solver_mpi.cpp`.** `rows_for_rank`,
`post_edge_rows`, `exchange_edge_rows` and `gather_rows` are `static`, which
keeps them private to that file, so their names cannot clash with anything
else in the program.

**11. Receives posted before sends.** `src/solver_mpi.cpp` lines 59 to 87.
`MPI_Irecv` and `MPI_Isend` return immediately, and `MPI_Waitall` blocks
until all four sends and receives have completed. Posting the receive first
means the library already knows the destination buffer (`row_below` or
`row_above`) when the matching message arrives, so it can copy the data
straight into it. Blocking sends in both directions could deadlock, because
each rank would block in its send while its neighbour is also blocked in a
send. `(int)requests.size()` is there because `MPI_Waitall` takes an `int`
count and `size()` returns an unsigned type.

**12. `MPI_Gatherv`.** `src/solver_mpi.cpp` lines 92 to 104. The ranks own
different numbers of rows when the rows do not divide evenly, and `Gatherv`
takes a count and an offset per rank, while `MPI_Gather` takes the same
count from every rank. It runs once after the time loop, so its cost is
outside the timed steps.

**13. `MPI_Init_thread` with `MPI_THREAD_FUNNELED`.** `src/main.cpp` line
117. It promises MPI that only the main thread will call it while OpenMP
threads exist. The plain `mpi` rung only calls MPI outside parallel regions.
The `mpi_overlap` rung calls it from the master thread inside one, which MPI
only allows from FUNNELED up, so that rung checks the level with
`MPI_Query_thread` and refuses to run below it.

**14. The E9 overlap region.** `src/solver_mpi.cpp` lines 284 to 337. This
is the most involved control flow in the project. It belongs to experiment
E9 and was never added to the ladder. Each thread works out its own
contiguous block of rows by hand with `rows_for_rank`, the same split
`schedule(static)` makes, so it knows before it starts whether its block
contains an edge row a neighbour needs. Those two rows are updated first.
`#pragma omp barrier` makes every thread wait until they are done.
`#pragma omp master` then lets only the master thread post the non-blocking
sends and receives. A `master` region has no barrier at its end, which
`single` would have, so the other threads continue with their remaining
rows while the sends and receives are outstanding.

**15. `slurm/flush_denormals.c`.** This is the most unusual file in the
project.
`__attribute__((constructor))` makes the function run before `main`. It sets
two bits in the SSE control register, flush to zero (0x8000) and denormals
are zero (0x0040), so numbers below about 2.2e-308 become zero and never
take the slow microcode path. It is compiled as a shared library and loaded
into `heston` with `LD_PRELOAD` by the bench scripts, so the program itself
is unchanged and the tests still run with full IEEE behaviour. Every price
is identical to ten digits with and without it (RESULTS M2.7).

**16. `MALLOC_MMAP_THRESHOLD_=1048576`.** `slurm/bench_common.sh`. It makes
glibc take every buffer over 1 MiB straight from `mmap`, so each repetition
gets fresh pages that the kernel backs with 2 MiB huge pages. Without it only
the first repetition gets them, which made rates vary by up to 1.8x in E11.

**17. The `(double)` casts.** For example `src/solver_omp.cpp` line 75. The
number of cell updates is computed in floating point because nodes times rows
times steps overflows a 32 bit `int` on the large runs.

**18. Makefile details.** `-DOMPI_SKIP_MPICXX` stops Open MPI 4.1 including
its old C++ bindings, which make g++ 8.5 print warnings under `-Wextra`.
`FSLIB` checks whether the compiler needs `-lstdc++fs` for
`std::filesystem`, which g++ 8 does and newer compilers reject.

**19. Experiment-only code.** `update_row_avx_load`,
`update_row_avx_shuffle` (with `_mm256_permute4x64_pd` and
`_mm256_blend_pd` masks such as `0b0001`), `OmpAlignSolver`, the float
solver and the other experiment solvers exist to answer one question each
and are never part of the ladder. `make_solver` lists them after the five
rungs.

---

## 6. A code tour for the big screen

Open the files in this order. Each step is one thing to point at.

1. `include/grid.h`. `index()` and the two buffers. A cell is
   `var_j * num_stock_nodes + stock_i`, so the stock axis is contiguous.
2. `src/solver_baseline.cpp`. The original loop nest with the divisions
   inside it. Every other rung must match its answer.
3. `src/stencil.cpp`. `StencilConstants`, where the reciprocals are worked
   out once, then `update_cells_scalar`, one cell per iteration.
4. `src/stencil_avx.cpp`. `update_row_avx`, the nine loads, the FMA chain and
   the scalar remainder.
5. `src/solver_omp.cpp`. The pragma at line 45, then the serial boundaries
   and the swap.
6. `src/solver_mpi.cpp`. `rows_for_rank`, `post_edge_rows` and the time loop
   of `MpiSolver::solve`.
7. `tests/test_ladder_matches.cpp`. How every rung is checked against the
   baseline and the parallel rungs against simd bit for bit.
8. `slurm/bench_common.sh`. Exclusive node, binding, huge pages, denormal
   flushing, tests before timing and the provenance header.
9. A results CSV. The header names the node, job, compiler, flags and
   binding.

---

## 7. How each figure was made

Every figure is drawn by a Python script in `scripts/` with matplotlib, from
files the rangpur jobs wrote. Nothing is drawn by hand, so `make m2-figures`
redraws them all.

| Figure | Script | Data | What the marks mean |
|---|---|---|---|
| Ladder, tried and rejected, strong, weak, grid size, both rooflines, step breakdown, placement, hybrid layouts, E9 overlap, split loads, E1 builds, correctness, -O levels, gprof | `scripts/m2_figures.py` | the sweep CSVs listed in `M2_CSVS` in the Makefile, the E1, E2, E2b and E9 CSVs, the peak file, one job log and the gprof outputs | Each bar or point is the median of 5 repetitions in one exclusive job. Rates are million cell updates per second. The grey line under each figure names the node, compiler, flags and job ids. |
| Decomposition diagrams | `scripts/decomposition_diagram.py` | none, drawn from the grid shape | Which cells one register, one thread and one rank own. |
| The grid step | `scripts/method_figures.py` | none, a schematic | The nine old cells used to compute one new cell, and which derivative each one is used in. |
| Black-Scholes check | `scripts/method_figures.py` | snapshot dumps from `heston --dump-every` on the collapse config | The solver's price curve drawn over the closed-form curve, and the gap between them in a panel below. |
| Payoff | `scripts/payoff_diagram.py` | none, the payoff formula | The call's payoff at expiry against the share price. |
| Value map | `scripts/method_figures.py` | snapshot dumps from `heston --dump-every`, in `results/demo` | Each pixel is one share price and variance, half way through the solve. The colour is the option value minus the payoff, which is the value of the time left. The circle marks today's share price and variance. |
| Code shots | `scripts/code_shots.py` | the source files | The real source lines, coloured by pygments, with the changed lines highlighted. The label above each pane names the file and lines. The SIMD before and after leaves out the comments, and the MPI time loop leaves out the comments and the step timers, which both labels say, so a slide can show fewer lines than the file. |

---

## 8. Questions to expect

Each question has a full answer in the place named after it. "HPC" means
`HPC_JUSTIFICATIONS.md`, and the experiments are in `docs/RESULTS.md` M2.6.

### The opening, slides 3 to 8

- What is a call option, and why is its price the discounted average payoff?
  It pays max(S - K, 0) at expiry. Its price today is the expected payoff
  under the risk-neutral measure, discounted at the rate r, because a hedged
  position earns the risk-free rate. (STUDY_GUIDE sections 1, 2 and 5b)
- Were the paths on slide 3 made by the solver? No. They illustrate the
  model and are drawn by `scripts/deck_m2_svgs.py` with the deck's Heston
  parameters. The solver has no random numbers. `make mc-check` is a separate
  Monte Carlo check, $196.13 plus or minus $0.55 against the grid's $196.17.
  (RESULTS section 7, Correctness)
- Where did the randomness go? The Feynman-Kac theorem turns the expected
  discounted payoff into the Heston PDE. The random parts of the model become
  the second derivative terms, so the solver steps a grid with no sampling.
  (STUDY_GUIDE section 5b, appendix The Heston model)
- Why does the solve run backwards in time? V is known only at expiry, where
  it is the payoff. Each step goes from t to t - dt, so 980,000 steps on the
  full grid reach today. (STUDY_GUIDE section 6.2)
- Why central differences? They take one cell on each side, so the error
  falls with the square of the spacing, second order. A one-sided difference
  is first order. They also need only the nine cells around each cell, so
  every cell can be updated on its own. Halving the spacing cut the
  error 3.78 times on the convergence run. (STUDY_GUIDE sections 6.6, 7.1
  and 7.5)
- Why are V_v, V_vv and V_Sv zero on slide 5? At expiry every variance row
  holds the same payoff, so north equals south and the four corners cancel.
  They become nonzero once the rows differ after the first step.
  (STUDY_GUIDE section 7.4)
- Why do the cells above and below get different values in the same step?
  The 1/2 v S^2 V_SS term grows with v. In the first step v = 0.06 gives
  0.07420, v = 0.04 gives 0.04959 and v = 0.02 gives 0.02498.
  (STUDY_GUIDE section 7.4)
- What sets the stability limit? The largest diffusion on the grid, at the
  top corner where v and S are largest. The limit is 1 / (v_max S_max^2 /
  dS^2 + xi^2 v_max / dv^2 + r), 5.66e-6 on the 421 by 51 test grid and
  3.66e-7 on the full grid, worked out in `src/solver.cpp` line 81. Above it
  the values grow without bound, which `config/unstable.cfg` shows.
  (STUDY_GUIDE section 8)
- Why not an implicit scheme, which has no limit? An implicit step solves a
  linear system across the whole grid, which couples every cell and is much
  harder to split between lanes, threads and processes. The explicit step is
  the same local stencil for every cell. (appendix Why an explicit scheme)
- Why multiply by reciprocals when the compiler could? Without
  `-ffast-math` the compiler may not, because 1/(2 dS) is itself rounded and
  the result can differ in the last bit. Doing it by hand gave 2.8x.
  (STUDY_GUIDE section 13, slides 14 and 15)
- How do the edge cells get their values? The boundary rules, applied after
  the interior cells in each step, by the main thread in the OpenMP rungs. (STUDY_GUIDE section 9, appendix The
  boundaries)

### The rungs and the measurements

- What is your project, in plain English? (section 3 of this file)
- Why is the problem parallel at all? (HPC section 1)
- Why do the threads share out rows instead of columns? (HPC section 5)
- Why is SIMD below 4x? (HPC section 4)
- Why write intrinsics when the compiler can vectorise? (E1 in RESULTS M2.6)
- What does `default(none)` do, and why are there no locks? (HPC section 5)
- Why 13.9x on 24 threads, and what does the Amdahl fit say? (HPC section 5)
- Why does the AVX2 row scale worse than the scalar row? Each AVX2 thread
  needs data at about three times the rate, so 24 of them run out of shared
  cache and memory bandwidth first. At 24 threads the AVX2 rung is still 2.16
  times the scalar rung. (RESULTS M2.2)
- Why did the Amdahl fraction double after the split loads were removed? A
  faster row shortens every step, while the boundaries and the fork and
  barrier cost the same, so the fixed part becomes a larger share. (RESULTS
  M2.2)
- Why does thread placement matter on this node? (HPC section 5)
- Why one MPI rank per socket on one node? (HPC section 6)
- Why exchange exactly one row, and why post the receives first? (HPC section 6)
- Why is two-node efficiency only 61 percent on the small grid when it is
  82 percent on the big one? The exchange is about 50 microseconds whatever
  the grid, which is 40 percent of the small grid's 0.126 ms step. (HPC
  section 6)
- How repeatable are the two-node numbers? The small grid's two-node layouts
  differ between jobs, so its gain lies between 1.21x and 1.33x. Every other
  layout agrees within 5 percent. (appendix A22, RESULTS M2.9)
- Did you try to overlap the communication with the computation? (E9 in
  RESULTS M2.6)
- Why must every rung match the baseline before it is timed? (HPC section 7)
- How did you keep the benchmarks from misleading you? (HPC section 8)
- Why was float rejected when it was faster? (E3 in RESULTS M2.6)
- Is comparing dynamic schedules with the omp rung fair, when they go
  through `schedule(runtime)`? The same solver with `OMP_SCHEDULE=static`
  ran at 6092 against the omp rung's 6374, and against that control dynamic
  1 is still 0.60x. The static schedule is faster against both references.
  (RESULTS M2.6)
- Why does the 2896 x 724 point drop out of the grid size figure? Its two
  buffers sit exactly 2^24 bytes apart and every rung runs slowly there. The
  trigger is known and the hardware cause is not. (appendix A26, RESULTS
  M2.8)

---

## 9. Status

- The no-split rerun finished on 2026-10-14 (CHANGELOG Phase 13, jobs
  659238 to 659249 and 659436, with the control job 659437). Every job
  passed the tests, and the deck, the script, the README and this file use
  its numbers.
- The two-node headline quotes job 659242, 8310 and 147x, because the hybrid
  figure is drawn from that job. The E9 job measured 9364 for the same layout
  minutes later. Appendix A22 shows both jobs side by side.
- Every figure and code image was redrawn on this machine with
  `make m2-figures` and `make code-shots`. F14 and F18 read local snapshots
  from `make collapse-anim` and `make demo-snapshots`, which need a local
  build of heston and so an MPI compiler.
- `make slides-m2` builds the deck with tectonic when it is installed and
  with pdflatex twice otherwise. The deck has 35 recorded frames, the
  divider and 32 backup slides, 68 pages in all.
- The script runs 9:43 at 140 words a minute.
- Still open are the date, the room and any rules for the interview from the
  course FAQ and profile, recording the video, and `make package` once the
  figures are final.
