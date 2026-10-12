# HPC Justifications. Why Each Technique Was Used

This file answers the question "why did you do X" for every high-performance
computing choice in the project, in full sentences that can be said out loud
in the interview. Each answer shows the code the work was done in, copied
from the real file with its file and line numbers and with the comments
taken out. Under the code, each line that matters is explained, and then the
answer says why we do it and gives the number that backs it. Every number
comes from an exclusive r730 node on rangpur and has its job id in
`docs/RESULTS.md`. The whole code base, file by file, is in `CODE_GUIDE.md`,
and the concepts behind it are in `STUDY_GUIDE.md`.

The machine, so the numbers make sense. An r730 node has two Xeon E5-2670 v3
sockets with 12 cores each, no hyper-threading, 30 MiB of L3 cache per
socket and one memory controller per socket. Each socket is one NUMA node.
The standard benchmark grid is 2048 stock nodes by 512 variance nodes, which
is 8 MiB per sheet and 16 MiB for the two sheets.

---

## 1. The problem and why it is parallel at all

### Why is this problem a good fit for parallel computing?

```cpp
// src/stencil.cpp:87-90 and 95-121 (excerpt)
const double* row_centre = current.data() + var_j * num_stock_nodes;
const double* row_above = row_centre + num_stock_nodes;
const double* row_below = row_centre - num_stock_nodes;
double* out_row = next.data() + var_j * num_stock_nodes;
...
double centre = row_centre[stock_i];
double east = row_centre[stock_i + 1];
...
out_row[stock_i] = centre + dt * change;
```

- `row_centre`, `row_above` and `row_below` all point into `current`, the
  old sheet, so every value the update reads comes from the previous step.
- `out_row` points into `next`, the new sheet, and the only write in the
  whole update is `out_row[stock_i]`.

Every timestep rewrites every interior cell from the nine cells around it
in the previous sheet. Because a step only reads the old sheet and only
writes the new one, no cell in a step depends on any other cell in the same
step. That means the cells can be updated in any order, by any number of
workers, and the answer does not change. This independence is what lets
four cells share one AVX2 instruction, lets 24 threads each own a band of
rows, and lets two machines each own half the rows.

### Why an explicit scheme and not an implicit one?

```cpp
// src/stencil.cpp:120-121
double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
out_row[stock_i] = centre + dt * change;
```

- `change` is the right-hand side of the PDE at this cell, built only from
  values that are already known.
- `centre + dt * change` gives the new value directly, with no equation to
  solve.

An explicit scheme computes each new cell directly from known values, so
the work is a stencil, which is regular, local and easy to split up. An
implicit scheme would need a linear solve across the whole grid every step,
which couples the cells together and is much harder to parallelise. The
price of the explicit choice is a stability limit on the timestep, so the
solver needs many small steps. That makes the problem compute-heavy and
repetitive, which is the kind of work parallel hardware is good at.

### Why two buffers instead of updating the grid in place?

```cpp
// src/grid.cpp:24-25 and 30-32
current_.assign(total_cells, 0.0);
next_.assign(total_cells, 0.0);

void Grid::swap_buffers() {
    std::swap(current_, next_);
}
```

- The two `assign` lines make two sheets of the same size, one to read and
  one to write.
- `std::swap(current_, next_)` exchanges the two vectors' internal pointers,
  so the sheet just written becomes the sheet read next step. No numbers are
  copied.

If the solver wrote new values into the same array it was reading, a cell
updated early in a step would be read as a neighbour by a later cell in the
same step, and the answer would depend on the order of the loop. Two
buffers keep every step a pure function of the previous sheet. That removes
the dependency between cells, so the double buffer is what makes the
parallel rungs possible. The swap costs nothing, because it only exchanges
pointers.

### Why was the hot spot measured before anything was optimised?

```cpp
// src/solver_serial.cpp:40-50 (excerpt)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_scalar(constants, current, next, var_j);
}
update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
grid.swap_buffers();
```

- The `for` loop is the interior update, which gprof found holds 99.4% of
  the serial rung's time inside `update_cells_scalar`.
- `update_zero_row`, `apply_stock_edges` and `copy_top_row` are the
  boundaries, which together held under half a percent.

Optimising code that does not hold the time is wasted effort, so the first
step was a gprof flat profile on one thread (`slurm/profile_m2.sh`). It put
99.9% of the baseline in its one solve function, 99.4% of the serial rung
and 98.9% of the SIMD rung in the interior row update. That told me the
interior update was the only thing worth attacking. Splitting the work into
small named functions is also what made the profile useful, because the
baseline's single function could only say that the time was somewhere in
one loop nest.

---

## 2. Memory layout

### Why is the grid one flat array and not a 2-D array of rows?

```cpp
// src/grid.cpp:23-25
std::size_t total_cells = (std::size_t)num_stock_nodes_ * num_variance_nodes_;
current_.assign(total_cells, 0.0);
next_.assign(total_cells, 0.0);
```

- `total_cells` is the number of cells in the whole grid, worked out in
  `std::size_t` so a large grid cannot overflow an `int`.
- Each `assign` makes one contiguous block of `total_cells` doubles.

A flat array is one contiguous block of memory, so the hardware prefetcher
can stream it and every cache line that is loaded is used in full. A vector
of vectors would put each row in a separate allocation somewhere on the
heap, and every access would go through an extra pointer. The flat layout
also means a whole row can be handed to MPI as one message with no packing.

### Why row-major with the stock axis contiguous?

```cpp
// include/grid.h:38
int index(int stock_i, int var_j) const { return var_j * num_stock_nodes_ + stock_i; }
```

- `var_j * num_stock_nodes_` jumps to the start of row `var_j`, so one step
  in variance is a jump of a whole row, 16 KiB at 2048 nodes.
- `+ stock_i` moves along the row, so one step in stock price is the next
  double in memory, 8 bytes away.

The inner loop runs along the stock axis, so it walks memory in order,
which is the access pattern caches and prefetchers are built for. The AVX2
rung depends on this as well, because the four cells in one register have
to sit next to each other in memory. In Milestone 1 a control that ran the
loops in the wrong order was 2.1 times slower on an exclusive node.

### Why do huge pages matter, and why are they pinned in every benchmark?

```bash
# slurm/bench_common.sh:22
export MALLOC_MMAP_THRESHOLD_=1048576
```

- This tells glibc that any allocation over 1 MiB must come straight from
  `mmap`, every time, instead of being reused from the heap. Fresh `mmap`
  memory of this size gets 2 MiB transparent huge pages from Linux.

The processor translates every virtual address to a physical one, and it
caches recent translations in the TLB. With normal 4 KiB pages the two
16 MiB sheets need 4096 translations, far more than the 1024-entry
second-level TLB holds. With 2 MiB huge pages they need only eight. In
Milestone 1 the same binary ran at 134.8 Mcell-updates/s with huge pages and
81.1 without, 1.66 times from page size alone. By default glibc only gives
huge pages to the first allocation of a run, so a median of five reps was
mixing two memory regimes. Pinning the threshold makes every rep get huge
pages. Experiment E11 removed the setting again with
`env -u MALLOC_MMAP_THRESHOLD_` (`slurm/bench_m2_exp_c.sh:31`), and the
OpenMP rung dropped to 0.68 times the pinned rate.

---

## 3. The serial rung

### Why replace divisions with multiplications by hand?

```cpp
// src/solver_baseline.cpp:66-67, before
double dV_dS = (east - west) / (2.0 * grid.stock_spacing());
double d2V_dS2 = (east - 2.0 * centre + west) / (grid.stock_spacing() * grid.stock_spacing());

// src/stencil.cpp:28-29, once per solve
inv_two_ds = 1.0 / (2.0 * stock_spacing);
inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);

// src/stencil.cpp:106-107, after, every cell
double dV_dS = (east - west) * inv_two_ds;
double d2V_dS2 = (east - 2.0 * centre + west) * inv_ds_sq;
```

- In the baseline, each derivative divides by a grid spacing, five divisions
  per cell per step.
- `inv_two_ds` and `inv_ds_sq` are the reciprocals, worked out once when the
  constants are built.
- In the serial rung, each derivative multiplies by the reciprocal instead.

A division costs 5 to 14 cycles on this CPU and does not pipeline, where a
multiplication costs about one. The compiler is not allowed to make this
change itself, because multiplying by a reciprocal can change the last bits
of the answer and `-ffast-math` is off. That one change took the rate from
57 to 159 Mcell-updates/s, 2.81 times the baseline. The rule I took from
Milestone 1 is that at `-O2` the hand optimisations worth doing are the ones
the compiler is forbidden to do.

### Why hoist the constants and the S and S² tables out of the time loop?

```cpp
// src/stencil.cpp:38-42 and 53-59
for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
    double price = stock_i * stock_spacing;
    stock_price[stock_i] = price;
    stock_price_sq[stock_i] = price * price;
}
for (int var_j = 0; var_j < num_variance_nodes; ++var_j) {
    double variance = var_j * variance_spacing;
    stock_diffusion[var_j] = 0.5 * variance;
    cross_diffusion[var_j] = (rho * xi) * variance;
    variance_diffusion[var_j] = ((0.5 * xi) * xi) * variance;
    variance_drift[var_j] = kappa * (theta - variance);
}
```

- The first loop fills a table of S and S² for every column.
- The second loop fills the four PDE coefficients for every variance row.
- Both run once, in the `StencilConstants` constructor, before the time
  loop starts.

Anything that does not change during the solve should be computed once, not
once per cell per step. These values depend only on the grid and the model,
so they are built once. In Milestone 1 the tables were rebuilt every step,
so this also removed one allocation per step.

### Why copy the constants into local variables inside the row function?

```cpp
// src/stencil.cpp:67-70 and 90
double dt = constants.dt;
double rate = constants.rate;
double stock_drift = constants.stock_drift;
double inv_two_ds = constants.inv_two_ds;
...
double* out_row = next.data() + var_j * num_stock_nodes;
```

- Each constant is read out of the `constants` struct once, at the top of
  the function, into a local variable.
- `out_row` is the pointer every new value is written through.

The compiler cannot prove that `out_row` does not point into the
`constants` struct. Without the local copies it would have to reload every
constant from memory after every store, in case the store changed it. A
local variable cannot be changed through `out_row`, so the compiler keeps
it in a register for the whole row.

### Why did the serial rung drop the other Milestone 1 optimisations?

In Milestone 1, loop reordering, loop splitting, induction variables and
unrolling by four each measured inside run-to-run noise, because GCC at
`-O2` already did them. Keeping code that does nothing would make the kernel
longer and harder to explain, so the serial rung keeps only the hoisting and
the reciprocals, the two changes that moved the number.

### Why did the serial rung leave SIMD as the only way forward?

The serial rung runs at 96% of the measured scalar peak of 5.12 GFLOP/s on
one core. A kernel that is already close to the arithmetic limit of one
instruction at a time can only go faster by doing more work per
instruction, which is what SIMD does.

---

## 4. SIMD with AVX2 and FMA

### Why does SIMD fit this stencil so well?

```cpp
// src/stencil_avx.cpp:56-61 and 110-112
while (stock_i + 3 <= last_interior) {
    __m256d centre = _mm256_loadu_pd(row_centre + stock_i);
    __m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
    __m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);
    ...
    _mm256_storeu_pd(out_row + stock_i, new_value);
    stock_i = stock_i + 4;
}
```

- `__m256d` is a 256-bit register that holds four doubles.
- `_mm256_loadu_pd(row_centre + stock_i)` loads four neighbouring cells,
  columns `stock_i` to `stock_i + 3`, in one instruction.
- `east` and `west` are the same four cells shifted one column right and
  left, so each lane gets its own neighbours.
- `_mm256_storeu_pd` writes four new values at once, and `stock_i` moves on
  by four.

Every cell in a row runs the same formula, and its neighbours in the stock
direction sit next to each other in memory. So four neighbouring cells can
be loaded into one register and updated by the same instruction. There are
no branches in the interior update, so all four lanes always do useful
work.

### Why vectorise by hand instead of leaving it to the compiler?

At `-O3` GCC 8's auto-vectoriser used two-lane SSE2, because the build does
not enable AVX2 by default, and it gained 1.51 times on the serial rung. The
hand-written AVX2 rung uses four lanes and FMA and reaches 379
Mcell-updates/s, 1.59 times the `-O3` serial rung and 2.39 times the `-O2`
serial rung. Writing it by hand also means I know exactly which instructions
run, and I can check the answer against the scalar row.

### Why use FMA?

```cpp
// src/stencil_avx.cpp:101-109
__m256d change = _mm256_mul_pd(stock_diffusion_factor, d2V_dS2);
change = _mm256_fmadd_pd(cross_factor, d2V_dS_dv, change);
change = _mm256_fmadd_pd(variance_diffusion, d2V_dv2, change);
change = _mm256_fmadd_pd(stock_drift_factor, dV_dS, change);
change = _mm256_fmadd_pd(variance_drift, dV_dv, change);
change = _mm256_fnmadd_pd(rate, centre, change);
__m256d new_value = _mm256_fmadd_pd(dt, change, centre);
```

- The first line starts the sum with the stock diffusion term.
- Each `_mm256_fmadd_pd(a, b, c)` works out `a * b + c` for four cells in one
  instruction, so each line adds one more PDE term to the running sum.
- `_mm256_fnmadd_pd(rate, centre, change)` works out `change - rate * centre`,
  the discounting.
- The last line is the time step, `centre + dt * change`, as one more FMA.

A fused multiply-add does a multiply and an add as one instruction with one
rounding, so it does twice the arithmetic per instruction and is slightly
more accurate. Eight of the 24 arithmetic instructions per four cells are
FMAs. Because FMA rounds once instead of twice, the SIMD rung is not bit
identical to the serial rung, but the gap is 1.1e-15, four orders of
magnitude inside the 1e-11 tolerance.

### Why AVX2 and not AVX-512, and why is only one file built with the flags?

```make
# Makefile:25-27 and 55-56
ifeq ($(shell uname -m),x86_64)
  AVXFLAGS = -mavx2 -mfma
endif

src/stencil_avx.o: src/stencil_avx.cpp
	$(CXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) $(AVXFLAGS) -c $< -o $@
```

- `AVXFLAGS` is set only on an x86_64 machine, so the arm64 Mac builds the
  scalar fallback.
- `-mavx2 -mfma` turn on exactly the Haswell features, nothing newer.
- The rule for `stencil_avx.o` is the only one that adds `$(AVXFLAGS)`.

The Xeon E5-2670 v3 is a Haswell chip, which supports AVX2 and FMA but not
AVX-512. Naming the two features, instead of `-march=native`, means the
binary uses only what every r730 node has, and `objdump` shows zero `zmm`
registers, which proves no AVX-512 instruction slipped in. Only one file
gets the flags because with `-mfma` set, GCC may fuse any multiply followed
by an add into an FMA in any file it compiles. If the whole program were
built that way, the serial rung's answers and timing would change too, and
the ladder would no longer show one change per rung.

### Why unaligned loads?

```cpp
// src/stencil_avx.cpp:59-61
__m256d centre = _mm256_loadu_pd(row_centre + stock_i);
__m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
__m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);
```

- The `u` in `loadu` means the address does not have to be a multiple of
  32 bytes.
- The three loads start 8 bytes apart, at `stock_i - 1`, `stock_i` and
  `stock_i + 1`.

Only one of three addresses 8 bytes apart can ever be 32-byte aligned, and
an aligned load (`_mm256_load_pd`) on the other two would crash. On Haswell
an unaligned load of data that happens to be aligned costs the same as an
aligned load, so there was little to gain. Experiment E2 would have tested
aligned buffers, and it was not run for that reason.

### Why is the SIMD speedup 2.4 times and not 4?

```cpp
// src/stencil_avx.cpp:62-67 and 116
__m256d north = _mm256_loadu_pd(row_above + stock_i);
__m256d south = _mm256_loadu_pd(row_below + stock_i);
__m256d north_east = _mm256_loadu_pd(row_above + stock_i + 1);
__m256d north_west = _mm256_loadu_pd(row_above + stock_i - 1);
__m256d south_east = _mm256_loadu_pd(row_below + stock_i + 1);
__m256d south_west = _mm256_loadu_pd(row_below + stock_i - 1);
...
update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
```

- Together with `centre`, `east` and `west`, these are nine overlapping
  loads for every four cells, plus two more for the S and S² tables.
- The last line hands the zero to three cells left at the end of the row to
  the scalar code, one at a time.

The loads grow faster than the useful arithmetic, and the remainder runs at
scalar speed. On the roofline the rung does 12 GFLOP/s against an AVX2
ceiling of 43.4 on one core, so it is limited by the loads and the
remainder, not by arithmetic.

---

## 5. OpenMP threads

### Why OpenMP, and why one pragma over the rows?

```cpp
// src/solver_omp.cpp:46-49
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
```

- `#pragma omp parallel for` starts a team of threads and shares the
  iterations of the next loop between them.
- The loop runs over variance rows, so each thread gets whole rows.
- Inside the loop is the same `update_row_avx` call as the SIMD rung, so
  each thread still does four cells per instruction.

The 24 cores of one node share memory, so threads can all read the same old
sheet without copying anything. OpenMP turns the row loop into a parallel
loop with one line, which keeps the change from the SIMD rung to one line.
The answer is bit for bit identical to the SIMD rung, and it gave 15.1 times
over SIMD on 24 threads, the largest single gain in the ladder.

### Why parallelise over rows and not columns?

A row is a contiguous block of memory, so a thread that owns a band of rows
streams through its own memory, and threads only meet at the edges of their
bands. Splitting by columns would give each thread a narrow strip of every
row, so neighbouring threads would write into the same cache lines and fight
over them, which is called false sharing. Putting the pragma on the outer
loop also means one fork and one barrier per step, and the inner loop stays
the AVX2 row, so SIMD and threads multiply.

### Why are no locks needed?

In the pragma above, `current` is only read and each thread writes only
`next` at its own rows `var_j`. No two threads ever write the same cell, and
nothing is summed across threads. So there is no data race, and there is no
`critical`, `atomic` or lock anywhere in the code.

### Why `schedule(static)`?

```cpp
// src/solver_experiments.cpp:132, the experiment E5 version
#pragma omp parallel for schedule(runtime) default(none) shared(constants, current, next, num_variance_nodes)
```

- `schedule(static)` in the rung splits the rows into one equal band per
  thread, once, with no bookkeeping while the loop runs.
- `schedule(runtime)` in the experiment reads the schedule from
  `OMP_SCHEDULE`, so one build could try `dynamic,1`, `dynamic,8` and
  `guided`.

Every row has the same number of cells and the same formula, so every row
costs the same. Experiment E5 measured the alternatives. `dynamic,1` ran at
0.60 times, `dynamic,8` at 0.72 times and `guided` at 0.84 times, because
dynamic scheduling hands out work at run time and breaks the contiguous
bands, for no benefit when the load is already balanced.

### Why `default(none)`?

`default(none)` in the pragma forces every variable used inside the parallel
loop to be listed. That is why the pragma names `constants`, `current`,
`next` and `num_variance_nodes` in `shared(...)`. A forgotten variable
becomes a compile error instead of silently becoming shared, which is how
most data races start.

### Why do the boundaries stay serial?

```cpp
// src/solver_omp.cpp:55-57
update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
```

- These three calls run after the parallel loop's implicit barrier, on the
  main thread only.
- `copy_top_row` goes last, because it copies the row below, which must be
  finished first.

The boundaries are about 5100 cells against about a million interior cells.
Leaving them on the main thread kept the OpenMP change to one pragma.
Experiment E8 (`OmpEdgesSolver`, `src/solver_experiments.cpp:177-185`) then
measured what that cost by moving the zero row and the stock edges inside
the parallel loop. It gained 1.05 times at 2048 by 512 and cut the serial
share of the step from 8.2% to 1.9%.

### Why does the OpenMP rung reach 15.8 times on 24 threads and not 24?

```python
# scripts/m2_figures.py:119-128
def amdahl_fit(threads, speedups):
    top = 0.0
    bottom = 0.0
    for p, s in zip(threads, speedups):
        x = 1.0 - 1.0 / p
        y = 1.0 / s - 1.0 / p
        top += x * y
        bottom += x * x
    return top / bottom
```

- Amdahl's law, `1/S = f + (1 - f)/p`, rearranges to `1/S - 1/p = f (1 - 1/p)`,
  which is a straight line through the origin with slope `f`.
- `x` and `y` are the two sides of that line for each thread count.
- `top / bottom` is the least squares slope, which is the serial fraction.

The fit gives a serial fraction of 1.18%. The boundaries alone are 0.36% of
a one-thread step, so most of the fraction is something else. One part is
the fork and barrier on every step. Another is that the serial boundary time
grows from 0.010 ms on one thread to 0.015 ms on 24, because the main thread
then has to fetch rows that other cores have just written. The third is
bandwidth. At 24 threads the AVX2 kernel moves about 138 GB/s through the
cache hierarchy, so it runs short of shared cache and memory bandwidth. The
scalar row moves data at half the rate and scales better, 18.5 times on 24
threads.

### Why one persistent parallel region was tried

```cpp
// src/solver_experiments.cpp:81-97 (excerpt)
#pragma omp parallel default(none) shared(constants, current, next, grid, num_variance_nodes, num_timesteps, dt)
{
    for (int step = 1; step <= num_timesteps; ++step) {
        #pragma omp for schedule(static)
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }
        #pragma omp single
        {
            update_zero_row(constants, current, next);
            ...
            grid.swap_buffers();
        }
    }
}
```

- `#pragma omp parallel` starts the team once, outside the time loop, so
  every thread runs every step.
- `#pragma omp for` shares the rows inside that existing team, with no new
  fork.
- `#pragma omp single` lets one thread do the boundaries and the swap while
  the others wait at its end.

Experiment E4 tested whether starting the team once is faster than forking
it every step. It gained 1.035 times at 2048 by 512 and 1.155 times at 512
by 128. The gain is real on small grids, where the per-step fork is a
larger share, but it is small on the benchmark grid, and keeping the
one-pragma rung made the step from SIMD to OpenMP easier to show.

### What happens if `OMP_NUM_THREADS` is larger than the number of cores?

```cpp
// src/solver_omp.cpp:20
result.threads = omp_get_max_threads();
```

- This records the thread count the runtime actually used, which goes into
  the CSV.

The node has 24 cores and no hyper-threading, so a 25th thread has no core
of its own. The operating system then time-slices two threads on one core,
and because every step ends in a barrier, the whole team waits for the
slowest thread every step. So asking for more threads than cores makes the
solver slower, not faster. That is why every sweep stops at 24 threads.

### Why does thread placement matter?

```bash
# slurm/bench_common.sh:25-26
export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
```

- `OMP_PLACES=cores` gives each thread a whole physical core.
- `OMP_PROC_BIND=close` pins threads to places in order and stops the
  operating system moving them, unless a sweep overrides it.

The node has two sockets, each with its own L3 cache and memory controller.
On the 4096 by 1024 grid, twelve threads on one socket ran at 1880
Mcell-updates/s and the same twelve split over both sockets ran at 3319,
because they had two caches and two memory controllers. On these nodes the
CPUs are numbered even on socket 0 and odd on socket 1, so `close`
alternates sockets and `spread` below 24 threads stays on one socket. I
printed the CPU each thread ran on (`slurm/diag_placement.sh`) to confirm
this.

### Why does page placement matter, and why was it tested with numactl?

```cpp
// src/grid.cpp:24
current_.assign(total_cells, 0.0);
```
```bash
# slurm/bench_m2_exp_b.sh:48
OMP_NUM_THREADS=24 RUN_PREFIX="numactl --interleave=all" run E7 "interleave_24_$1x$2" $(grid_args $1 $2 24) --solver omp
```

- `assign` writes a zero into every element on the main thread, so the main
  thread is the first to touch every page.
- `numactl --interleave=all` runs the same binary with its pages spread
  evenly over both sockets' memory.

Linux places a page in the memory of the socket whose thread first writes
it, and does not move it afterwards. So every page of the grid ends up on
socket 0, and the threads on socket 1 read across the link between sockets.
Fixing that in code would mean replacing the vectors in `Grid` with
uninitialised storage in every rung. Running under `numactl` tests the same
question without changing code. It gave 1.065 times at 2048 by 512 and
1.28 times at 4096 by 1024, and much less spread between reps.

---

## 6. Hybrid MPI and OpenMP

### Why add MPI at all, when OpenMP already uses the whole node?

```cpp
// src/main.cpp:119-123
MPI_Init(&argc, &argv);
...
int exit_code = run(argc, argv, rank);
MPI_Finalize();
```

- `MPI_Init` starts MPI in every process, and `MPI_Finalize` shuts it down.
- The program body is in `run`, so `MPI_Finalize` is reached on every exit
  path, errors included.

OpenMP only works inside one shared memory, so it cannot use a second node.
MPI runs separate processes that share nothing and send messages, so it can
span machines. Using MPI between processes and OpenMP inside each one keeps
the number of messages small while still using every core.

### Why split the grid into slabs of rows?

```cpp
// src/solver_mpi.cpp:26-40
static RowSlab slab_for_rank(int num_rows, int rank, int num_ranks) {
    int rows_per_rank = num_rows / num_ranks;
    int extra_rows = num_rows % num_ranks;
    RowSlab slab;
    if (rank < extra_rows) {
        slab.begin = rank * (rows_per_rank + 1);
        slab.end = slab.begin + rows_per_rank + 1;
    } else {
        slab.begin = extra_rows * (rows_per_rank + 1) + (rank - extra_rows) * rows_per_rank;
        slab.end = slab.begin + rows_per_rank;
    }
    return slab;
}
```

- `rows_per_rank` is the even share and `extra_rows` is what is left over.
- The first `extra_rows` ranks take one extra row each, so no rank has more
  than one row more than another.
- The `else` branch skips past the ranks that took an extra row.

Each rank owns a contiguous slab of variance rows. A row is contiguous in
memory, so the edge row a neighbour needs goes as one message straight out
of the array with no packing. Rows also match the OpenMP split, so each
rank runs the same parallel row loop over its own slab
(`src/solver_mpi.cpp:136-139`).

### Why is the halo exactly one row, and why non-blocking messages with the receives first?

```cpp
// src/solver_mpi.cpp:47-50 and 57-80 (excerpt)
double* first_row = sheet.data() + slab.begin * num_stock_nodes;
double* last_row = sheet.data() + (slab.end - 1) * num_stock_nodes;
double* halo_below = sheet.data() + (slab.begin - 1) * num_stock_nodes;
double* halo_above = sheet.data() + slab.end * num_stock_nodes;
...
if (rank > 0) {
    MPI_Irecv(halo_below, num_stock_nodes, MPI_DOUBLE, rank - 1, halo_tag, MPI_COMM_WORLD, &receive_below);
    MPI_Isend(first_row, num_stock_nodes, MPI_DOUBLE, rank - 1, halo_tag, MPI_COMM_WORLD, &send_below);
}
...
MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
```

- `first_row` and `last_row` are this rank's own edge rows, which the
  neighbours need.
- `halo_below` and `halo_above` are the rows just outside the slab, where
  the neighbours' edge rows land.
- `MPI_Irecv` posts the receive first and returns straight away, then
  `MPI_Isend` starts the send and also returns straight away.
- Each message is `num_stock_nodes` doubles, one whole row, 16 KiB at 2048
  nodes.
- `MPI_Waitall` waits until every send and receive has finished, before
  anything touches those rows.

The stencil reaches one cell in each direction, so to update its first and
last rows a rank needs exactly one row from each neighbour. If two
neighbours both did a blocking send first, each could wait for the other to
receive, and the program would deadlock. Non-blocking calls with one
`Waitall` make the order the ranks arrive in irrelevant. Posting the
receives before the sends means MPI already knows where an arriving row
goes, so it can land straight in the array instead of in a temporary
buffer.

### Why one rank per socket on a single node?

```bash
# slurm/bench_m2_mpi.sh:40
run_mpi M "1node_2x12_$SIZE" "-n 2 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12" $(grid_args $1 $2 24)
```

- `-n 2` starts two processes.
- `--map-by ppr:1:socket:PE=12` places one process per socket and gives
  each one 12 cores.
- `-x OMP_NUM_THREADS=12` gives each process 12 OpenMP threads.
- `--report-bindings` prints where each process landed, so the placement is
  in the job log.

Each rank creates its own `Grid`, so its main thread first-touches its own
copy and the pages land on its own socket. That fixes the page placement
problem from the OpenMP rung without changing any code. Two ranks of 12
threads beat one rank of 24 threads by 1.09 times on the 2048 grid and by
2.16 times on the 4096 grid, where memory placement matters most.

### Why does every rank keep the whole grid?

```cpp
// src/solver_mpi.cpp:102 and 124-125
Grid grid(cfg.grid, cfg.option);
...
int first_interior = std::max(slab.begin, 1);
int end_interior = std::min(slab.end, num_variance_nodes - 1);
```

- Every rank makes a full-size grid, then only updates the rows from
  `first_interior` to `end_interior`.
- `std::max(slab.begin, 1)` and `std::min(slab.end, num_variance_nodes - 1)`
  cut off row 0 and the top row, which are boundaries.

Keeping the full grid on every rank means row indices are the same
everywhere, so the shared stencil functions are reused unchanged. The cost
is 16 MiB per rank at 2048 by 512, which is small on a node with plenty of
memory. For a grid too big for one node I would allocate only the slab plus
two halo rows and shift every index by the slab's first row.

### Why `MPI_Gatherv` at the end?

```cpp
// src/solver_mpi.cpp:170-177
for (int other_rank = 0; other_rank < num_ranks; ++other_rank) {
    RowSlab other_slab = slab_for_rank(num_variance_nodes, other_rank, num_ranks);
    counts[other_rank] = (other_slab.end - other_slab.begin) * num_stock_nodes;
    offsets[other_rank] = other_slab.begin * num_stock_nodes;
}
...
MPI_Gatherv(my_rows, my_count, MPI_DOUBLE, next.data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
```

- `counts` holds how many doubles each rank sends, which differs when the
  rows do not divide evenly.
- `offsets` holds where each rank's rows go in rank 0's sheet.
- `MPI_Gatherv` copies every slab into rank 0's spare sheet in one call.

Only rank 0 needs the finished grid to read off the price and Greeks. The
`v` in `Gatherv` is what allows each rank to send a different amount. The
gather runs after the timer stops, so it is not part of any benchmark
number.

### Why is two-node efficiency only 70%?

```cpp
// src/solver_mpi.cpp:155-158
auto halo_start = std::chrono::steady_clock::now();
exchange_halo_rows(next, num_stock_nodes, slab, rank, num_ranks);
auto halo_end = std::chrono::steady_clock::now();
result.halo_seconds += std::chrono::duration<double>(halo_end - halo_start).count();
```

- The halo exchange is timed on its own and written to the `halo_seconds`
  column, so its share of each step can be measured.

Across two nodes the best layout is 1.41 times one node at 2048 by 512 and
1.57 times at 4096 by 1024. The halo exchange takes about 45 microseconds
per step, which is 38% of a 0.12 ms step on the small grid. Sending 64 KiB
should take far less than that at any network speed, so the cost is the
per-message latency and waiting for the slower neighbour, not bandwidth.
The runs went over TCP because the job excludes Open MPI's fast transports
(`slurm/bench_common.sh:17-18`), as the matrixass MPI runs also needed. The
bigger grid is more efficient because each step has more interior work to
hide the same halo cost.

### Why was overlapping the halo with computation not done?

Overlap would start the exchange, update the interior rows that do not need
the halo, then wait and finish the edge rows. It is the obvious next step
given that the halo is 38% of a short step. It was on the cut list because
the rungs, sweeps and figures came first, and it needs `exchange_halo_rows`
split into a start and a finish around the update.

---

## 7. Correctness before speed

### Why must every rung match the baseline before it is timed?

```bash
# slurm/bench_common.sh:62-68
build_and_validate() {
  make clean
  make all || { echo "FATAL: build failed, refusing to benchmark"; exit 1; }
  make test || { echo "FATAL: make test failed, refusing to benchmark"; exit 1; }
  make test-mpi || { echo "FATAL: make test-mpi failed, refusing to benchmark"; exit 1; }
  make avx-check || { echo "FATAL: make avx-check failed, refusing to benchmark"; exit 1; }
}
```

- `make clean` then `make all` builds from scratch on the compute node.
- Each `||` stops the whole job with `exit 1` if that step fails, so a
  failing test means no timings are written.

A faster answer that is wrong is worth nothing, and a bug in a parallel
rung often shows up as a small numerical change rather than a crash. Every
rung must match the baseline to 1e-11 on price, delta, gamma and vega
(`kRelativeTolerance`, `tests/test_ladder_matches.cpp:25`). The measured
worst gap is 3.3e-15.

### Why must OpenMP and MPI match the SIMD rung exactly, to zero difference?

```cpp
// tests/test_ladder_matches.cpp:21 and 118-120
static const std::vector<std::string> kSameBitsAsSimd = {"omp", "mpi", "omp_persistent", "omp_runtime", "omp_edges"};
...
double simd_diff = max_relative_diff(simd_reference, candidate);
bool ok = simd_diff == 0.0;
```

- `kSameBitsAsSimd` lists every solver that must give exactly the SIMD
  rung's answer.
- `simd_diff == 0.0` allows no difference at all, not even in the last bit.

Each cell runs the same instructions in the OpenMP and MPI rungs as in the
SIMD rung, a halo row is an exact copy of 64-bit values, and nothing is
summed across threads or ranks. So any difference at all means a row was
skipped, done twice, or raced. A tolerance of zero catches errors that a
looser tolerance would hide. As a check of the check, removing the halo
exchange moved the two-rank answer by 2.0e-9 and the test failed.

### Why test at 7 threads and 3 ranks?

```make
# Makefile:75 and 83
	@for n in 1 2 7; do echo "== tests/test_ladder_matches with OMP_NUM_THREADS=$$n"; OMP_NUM_THREADS=$$n ./tests/test_ladder_matches 2>/dev/null || exit 1; done
	@for n in 2 3; do echo "== tests/test_ladder_matches on $$n processes"; OMP_NUM_THREADS=2 $(MPIEXEC) -n $$n ./tests/test_ladder_matches 2>/dev/null || exit 1; done
```

- `make test` runs the ladder test at 1, 2 and 7 threads, and `|| exit 1`
  stops at the first failure.
- `make test-mpi` runs it as 2 and 3 processes with 2 threads each.

Seven threads do not divide the 30 interior rows of the test grid, and three
ranks split its 32 rows into 11, 11 and 10. Uneven splits are where an off
by one error at a band edge would show, so the tests deliberately use them.

---

## 8. Measurement method

### Why benchmark only on compute nodes through sbatch?

```python
# scripts/m2_figures.py:62-69
def check_provenance(provenance):
    bad = []
    for path, lines in provenance.items():
        hosts = [line.split("=", 1)[1] for line in lines if line.startswith("host=")]
        if not hosts or not hosts[0].endswith(COMPUTE_NODE_SUFFIX):
            bad.append(path)
    return bad
```

- Each CSV header has a `host=` line written by the job.
- Any file whose host is missing or is not a rangpur compute node is
  collected in `bad`, and the script stops without drawing a figure.

The login node is shared by every user, so its timings measure other
people's work as well as mine, and the course rubric does not accept them.
A compute node through `sbatch --exclusive` gives the job the whole node.
In Milestone 1 a control measured 30 times slower while another job was
running and 2.1 times on an exclusive node, which showed how easily a shared
machine gives the wrong answer.

### Why the median of five reps?

```python
# scripts/m2_figures.py:98-102
rates = [r["cups"] for r in group]
out[label] = {
    "median": statistics.median(rates),
    "min": min(rates),
    "max": max(rates),
```

- `REPS=5` (`slurm/bench_common.sh:29`) makes `heston --bench 5` run each
  solve five times.
- The figure script takes the median of the five rates and keeps the min and
  max beside it.

The median ignores a single slow rep caused by something outside the
program, such as an operating system task. Five reps are enough to see the
spread, and the min and max show it on every bar as an error bar.

### Why does every CSV carry a provenance header?

```bash
# slurm/bench_common.sh:75-89 (excerpt)
echo "# host=$(hostname)"
echo "# compiler=$(mpicxx --version | head -1)"
echo "# page_regime=MALLOC_MMAP_THRESHOLD_=${MALLOC_MMAP_THRESHOLD_}, THP $(cat /sys/kernel/mm/transparent_hugepage/enabled)"
echo "# omp_places=${OMP_PLACES} omp_proc_bind=${OMP_PROC_BIND}"
lscpu | grep -E '^(Socket|Core|Thread|NUMA node)' | sed 's/^/# /'
```

- Each line records one condition the numbers depend on, the machine, the
  compiler, the page regime, the thread binding and the socket and NUMA
  layout.
- The lines start with `#` so the plotting script skips them as data.

Every number can then be traced to the exact conditions it came from, and
anyone can check whether two CSVs were measured the same way before
comparing them.

### Why are denormals flushed to zero in every timed run?

```c
// slurm/flush_denormals.c:12-14
__attribute__((constructor)) static void flush_denormals(void) {
    _mm_setcsr(_mm_getcsr() | 0x8040);
}
```
```bash
# slurm/bench_common.sh:45 and 114 (excerpt)
gcc -O2 -shared -fPIC slurm/flush_denormals.c -o "$WORK_DIR/flush_denormals.so"
env "$(preload_setting)" ${RUN_PREFIX:-} ./heston "$@" --bench "$REPS"
```

- `__attribute__((constructor))` makes the function run before `main`.
- `_mm_getcsr()` reads the MXCSR register, which controls floating point
  behaviour, and `| 0x8040` sets bit 15 (FTZ, tiny results become zero) and
  bit 6 (DAZ, tiny inputs are read as zero).
- The script builds it as a shared library and loads it into `heston` with
  `LD_PRELOAD`, so `heston` itself is unchanged.

Short benchmark runs leave many cells holding numbers below 2.2e-308, which
this CPU handles through a slow microcode path. How many such cells there
are depends on the grid and the number of steps, so rungs run with
different step counts were being compared in different numerical regimes.
Flushing gained up to 1.38 times on the SIMD rung and left every price
identical to ten digits. It is the same idea as pinning the page regime.
The thing being measured is the code, so every run has to happen under the
same conditions.

### Why does each job build in its own copy of the tree?

```bash
# slurm/bench_common.sh:34-39
WORK_DIR=$HOME/cosc3500/work/job_${SLURM_JOB_ID}
mkdir -p "$WORK_DIR" "$SUBMIT_DIR/results"
rsync -a --exclude results --exclude '*.o' --exclude heston --exclude 'slurm-*.out' ./ "$WORK_DIR/"
cd "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT
```

- `WORK_DIR` is named after the job id, so no two jobs share it.
- `rsync` copies the source in, and the job builds and runs there.
- `trap ... EXIT` deletes the copy however the job ends.

The home directory is shared between jobs. Twice in the first round of
sweeps, a sync or a `make clean` from another job deleted files while a
benchmark was running and killed it. A private copy per job means nothing
outside the job can touch its files.

### Why measure the peaks instead of quoting the datasheet?

A datasheet gives the peak under ideal conditions, not what this node
actually delivers. `slurm/bench_peak.sh` measured a scalar peak of 5.12
GFLOP/s, an AVX2 FMA peak of 43.4 GFLOP/s on one core and 971 on 24, and a
triad memory bandwidth of 67.9 GB/s on 24 threads. Those are the ceilings
the roofline uses, so every rung is compared against what is really
achievable. The first triad run reported only 24.8 GB/s because its arrays
were `std::vector`, zero-filled on the main thread, the same first-touch
problem as the grid. Allocating them untouched and filling them in parallel
gave 67.9.

---

## 9. Scaling laws and models

### Why use Amdahl's law?

Amdahl's law says the speedup on p threads is limited by the part of the
work that stays serial. Fitting it to the strong scaling data (the
`amdahl_fit` code in section 5) gives one number, the serial fraction, that
summarises how well the code scales. Here it is 1.18%, which caps the
speedup at about 85 however many threads are added.

### Why also measure weak scaling with Gustafson's law?

```python
# scripts/m2_figures.py:252 and 256-258
scaled = [p * step_ms[0] / t for p, t in zip(threads, step_ms)]
...
widest = points[-1][1]
serial_fraction = widest["boundary_ms"] / widest["step_ms"]
gustafson = [p - serial_fraction * (p - 1) for p in threads]
```

- `scaled` is the scaled speedup, `p` times the one-thread step time over
  the `p`-thread step time, because each thread always has the same work.
- `serial_fraction` is the measured serial share of the widest run.
- `gustafson` is the line Gustafson's law predicts from that share.

Strong scaling keeps the problem fixed, so on many threads each thread has
very little work and the fixed costs dominate. Weak scaling grows the grid
with the threads so each thread always has 21 rows. The step time rose
from 0.142 ms on one thread to 0.179 ms on 24, a scaled speedup of 19.1
against an ideal 24. That shows the code would use more cores well if the
problem grew with them, which is how a real user would use more hardware.

### Why use a roofline?

```python
# scripts/m2_figures.py:517-519 and 559
FLOPS_PER_CELL = 31
BYTES_PER_CELL = 24
...
gflops = ladder[key]["median"] * FLOPS_PER_CELL / 1e9
```

- `FLOPS_PER_CELL` is the arithmetic in one cell update, counted off the
  kernel, with an FMA counted as two.
- `BYTES_PER_CELL` is the memory traffic that cannot be avoided. Each cell
  reads one double of `current` and writes one double of `next`, and the
  store first reads its cache line in (write-allocate), so 3 x 8 bytes. The
  neighbour rows are assumed to be reused from cache.
- `gflops` turns a measured rate in cell updates per second into GFLOP/s so
  it can be placed under the measured roofs.

The roofline compares a kernel's arithmetic intensity, flops per byte
moved, with the machine's measured compute and memory ceilings. It tells me
whether a rung is limited by arithmetic or by memory, and so what kind of
change could make it faster. The intensity here is 31/24 = 1.29. On 24
threads at 2048 by 512 the kernel reaches 178 GFLOP/s, twice the DRAM roof,
because the grid fits in the two L3 caches. At 4096 by 1024 it falls to 93
GFLOP/s, on the DRAM roof, because the grid no longer fits.

### Why does grid size change the result so much?

Small grids lose to the fork and barrier every step, which costs the same
whatever the size, so 512 by 128 only reaches 8 times serial. The 2048 grid
fits in the two 30 MiB L3 caches and reaches 36 times serial. At 4096 by
1024 the sheets spill to main memory and the rate halves. So the best
speedup comes from a grid large enough to amortise the per-step costs and
small enough to stay in cache.

---

## 10. Experiments that were rejected

### Why was float rejected when it was up to 4.46 times faster?

```cpp
// src/stencil_avx.cpp:121 and 125
void update_row_avx_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j) {
    __m256 dt = _mm256_set1_ps(constants.dt);
```

- `std::vector<float>` stores each cell in 4 bytes instead of 8.
- `__m256` holds eight floats in the same 256 bits that hold four doubles,
  and the `_ps` functions are the float versions of the `_pd` ones.

Float halves the bytes per cell and doubles the lanes per register, so it
was 1.75 times faster at 2048 by 512 and 4.46 times at 4096 by 1024, where
the float sheets fit in cache and the double sheets do not. But on the
reference grid it moved the price by $0.044 from the Fourier price, where
double is $0.0008 away, and it moved gamma by 11.9%. Rounding every cell to
float on every one of 980000 steps builds up error. A pricer that is fast
and wrong is not useful, so float stayed an experiment.

### Why not dynamic scheduling?

Every row costs the same, so there is no imbalance for dynamic scheduling to
fix, and it adds runtime overhead and breaks the contiguous bands. It lost
up to 40% (experiment E5, the `schedule(runtime)` code in section 5).

### Why was `-ffast-math` not used?

It lets the compiler reorder floating point operations, which changes the
answer in ways I cannot predict or test against the baseline. In Milestone 1
it broke the bit-identical match between the first two levels, so it failed
the correctness rule before it was timed.

---

## 11. What I would do next

The measurements point at three changes, in order of expected value. First,
overlap the halo exchange with the interior update, since the halo is 38% of
a short step across two nodes. Second, first-touch the grid in parallel by
replacing the zero-filled vectors with uninitialised storage, which E7
showed is worth up to 28% on a large grid. Third, tile the time loop so
several steps run on a block of rows while it is still in cache, which would
attack the memory roof that the large grid sits on.
