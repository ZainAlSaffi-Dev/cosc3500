# Study Cards. The Heston Pricer in Short

These cards are the short form of [STUDY_GUIDE.md](../STUDY_GUIDE.md). They
hold every idea needed to explain the code, and the numbers to quote with
them. The Heston Pricer Study Guide page is built from this file and from
[CODE_GUIDE.md](../CODE_GUIDE.md). References like §8.2 point at the code
guide, and references like §I2.4 point at a card here.

---

## I1. The problem

### I1.1 Options and payoffs

An option is the right, but not the duty, to trade at a fixed price on a
fixed date. A call is the right to buy at the strike K, and a put is the right
to sell at K. Ours are European, so they can only be used on the expiry date.
The contract is on the S&P 500 with K = 5250 and T = 0.25 years, with today's
index at 5200.

```text
call payoff at expiry:  max(S - K, 0)
put  payoff at expiry:  max(K - S, 0)
```

The payoff fills the first sheet of the grid (`Grid::init_payoff`, §5.5).
Its kink at the strike is where all the later action starts, because it is
the only place the payoff is curved.

### I1.2 The price and the Greeks

The payoff on the expiry date is easy. The question is what the contract is
worth today, while the final index level is still unknown. That fair value is
the number the program prints.

The Greeks say how the value moves when the market moves. Traders use them to
hedge.

| Greek | Meaning | On the grid |
|---|---|---|
| delta | change in value per $1 move in S | slope along a row, `(east - west) / (2 ds)` |
| gamma | how fast delta changes | curvature along a row, `(east - 2 centre + west) / ds²` |
| vega | change in value as variance moves | slope up a column, `(north - south) / (2 dv)` |

Once the last sheet is filled they cost almost nothing, because they are
differences between neighbouring cells (`extract_result`, §6.6).

### I1.3 The Heston model

Volatility is how much the price moves around, and variance v is volatility
squared. Black-Scholes treats volatility as a constant, which gives a formula
but does not match markets, where calm and panic alternate. Heston makes the
variance random too.

```text
dS = (r - q) S dt + sqrt(v) S dW1                 the stock
dv = kappa (theta - v) dt + xi sqrt(v) dW2        its variance
corr(dW1, dW2) = rho
```

| Parameter | Meaning | Value |
|---|---|---|
| `v0` | today's variance | 0.04, which is 20% volatility |
| `kappa` | how fast variance returns to normal | 1.5 |
| `theta` | the normal long-run variance | 0.04 |
| `xi` | vol of vol, how hard variance is kicked | 0.35 |
| `rho` | correlation of stock moves and variance moves | -0.70 |

rho is negative because crashes and volatility spikes arrive together. The
rate is r = 0.045 and the dividend yield is q = 0.013. Black-Scholes is the
special case xi = 0, which the collapse test uses (§I2.8).

---

## I2. From model to grid

### I2.1 The PDE and its terms

No-arbitrage theory turns the two random equations into one equation with no
randomness in it, which the option value V(S, v, t) must satisfy.

```text
V_t + (r-q) S V_S + kappa (theta - v) V_v
    + ½ v S² V_SS + ½ xi² v V_vv + rho xi v S V_Sv  =  r V
```

Each term is one named variable in `update_cells_scalar` (§8.6).

| PDE term | What it does | Code variable | Coefficient in `StencilConstants` |
|---|---|---|---|
| ½ v S² V_SS | spreads value along S, gamma lives here | `stock_diffusion_term` | `stock_diffusion` = 0.5v, `stock_price_sq` = S² |
| rho xi v S V_Sv | the correlation, needs the corner cells | `cross_term` | `cross_diffusion` = ρξv |
| ½ xi² v V_vv | spreads value along v | `variance_diffusion_term` | `variance_diffusion` = 0.5ξ²v |
| (r-q) S V_S | drift of the stock | `stock_drift_term` | `stock_drift` = r - q |
| kappa (theta - v) V_v | mean reversion of variance | `variance_drift_term` | `variance_drift` = κ(θ - v) |
| r V | discounting | `rate * centre` | `rate` |
| V_t | time passing | `centre + dt * change` | `dt` |

Every second-derivative term has a factor v, so at v = 0 all the diffusion
switches off. That shapes the bottom boundary (§I2.7).

### I2.2 Where the strike and the randomness went

The strike is not in the PDE, because pricing has two separate parts. The
model (the PDE) says how the index behaves and fits any contract. The contract
(the payoff and the edge values) says what you are owed. So K lives in
`init_payoff` and `apply_stock_edges`, never in the stencil. Pricing a put
instead of a call changes only those two.

There is no random number generator anywhere in the C++. The fair price is an
average over random futures, and the Feynman-Kac theorem says that average
satisfies the PDE. So we solve the equation the average obeys instead of
drawing paths.

| Route | Cost |
|---|---|
| Monte Carlo, millions of random paths | error shrinks like 1/sqrt(N), slow |
| PDE, ours | error shrinks with the grid spacing, and every (S, v) comes out at once |

Over one small step the noise moves S up or down by about sqrt(v) S. Averaging
the two cases cancels the slope, and what survives is ½ (sqrt(v) S)² V_SS,
which is ½ v S² V_SS. That is why only v appears in the code, never sqrt(v).
The cross term has no ½, because the mixed partial appears twice and the ½
cancels. An extra 0.5 there is the most common Heston bug.

`scripts/monte_carlo_check.py` checked this with three methods that share no
code. On `demo.cfg` the PDE gave 196.1054, a 400k path Monte Carlo 196.1345
± 0.55, and the Fourier formula 196.1692.

### I2.3 The stack of sheets

One sheet is a table of option values at one moment. Columns are stock prices
(2048 of them, from 0 to 4K), rows are variance levels (512 of them), and each
cell is the value in that scenario.

```text
  var_j
nv-1|  T  T  T  T  T  T  T  T     T  copy of the row below      copy_top_row
    |  E  I  I  I  I  I  I  E     I  full 9-point update        update_row_*
  2 |  E  I  I  I  *  I  I  E     *  today's market, S = 5200, v = 0.04
  1 |  E  I  I  I  I  I  I  E     E  known edge value           apply_stock_edges
  0 |  E  Z  Z  Z  Z  Z  Z  E     Z  drift-only update, v = 0   update_zero_row
    +---------------------------> stock_i
```

Only the expiry sheet is known, so the solver starts there and steps backwards
to today, one sheet per timestep. Each sheet is computed only from the one
before it, so only two sheets ever exist in memory (§I3.1). The price is read
at the * cell of the last sheet. Every parallel rung splits this rectangle
into rows.

### I2.4 Finite differences

A computer has cells, not derivatives. Finite differences turn each derivative
into arithmetic on a cell and its neighbours.

```text
                 stock_i - 1   stock_i      stock_i + 1
  var_j + 1      north_west    north        north_east     row_above
  var_j          west          centre       east           row_centre
  var_j - 1      south_west    south        south_east     row_below
```

| Derivative | Formula | Reads |
|---|---|---|
| `dV_dS` | (east - west) / (2 ds) | west, east |
| `d2V_dS2` | (east - 2 centre + west) / ds² | west, centre, east |
| `dV_dv` | (north - south) / (2 dv) | south, north |
| `d2V_dv2` | (north - 2 centre + south) / dv² | south, centre, north |
| `d2V_dS_dv` | (NE - NW - SE + SW) / (4 ds dv) | the four corners |

The cross derivative is the S slope on the row above minus the S slope on the
row below, divided by 2 dv. It is the only one that reads the corners, so the
correlation is why the stencil has nine points and not five. The `v = 0` row
uses the one-sided `(north - centre) / dv`, because there is no row below.

One explicit Euler step adds dt times the PDE's rate of change.

```cpp
double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
out_row[stock_i] = centre + dt * change;
```

It only reads the current sheet, so every cell can be done on its own, in any
order. Every parallel rung relies on that.

### I2.5 One step worked by hand

The test grid is 421 by 51 with 56000 steps, so ds = 50, dv = 0.02 and
dt = 4.464e-6. Take the cell at the strike, S = 5250 and v = 0.04, in the
first step. The nine values come from the payoff sheet.

```text
              S = 5200    S = 5250    S = 5300
  v = 0.06        0           0          50
  v = 0.04        0         [ 0 ]        50
  v = 0.02        0           0          50
```

```text
dV_dS   = (50 - 0) / 100        = 0.5      every v derivative is 0,
d2V_dS2 = (50 - 0 + 0) / 2500   = 0.02     because the payoff does not depend on v

stock_diffusion_term = 0.5 × 0.04 × 5250² × 0.02   = 11025
stock_drift_term     = 0.032 × 5250 × 0.5          = 84
change = 11109,  new value = 0 + 4.464e-6 × 11109  = 0.04959
```

The solver's own snapshot holds 0.04959375. The cell was worth nothing at
expiry and is worth 5 cents one step earlier, and almost all of it came from
the curvature at the kink. Rows above and below give 0.07420 and 0.02498,
growing with v. Today's spot at 5200 stays at 0 until step 2, because
information moves one cell per step.

### I2.6 The nine weights and the stability limit

Regroup the update by cell and each of the nine cells gets one weight. With
A = ½vS², B = (r-q)S, C = ½ξ²v, D = κ(θ-v), X = ρξvS:

```text
east    dt (A/ds² + B/(2ds))        west    dt (A/ds² - B/(2ds))
north   dt (C/dv² + D/(2dv))        south   dt (C/dv² - D/(2dv))
NE, SW  dt X/(4 ds dv)              NW, SE  -dt X/(4 ds dv)
centre  1 - dt (2A/ds² + 2C/dv² + r)
```

For the worked cell the centre weight is 0.997976 and the neighbours are
around 0.001. The weights add up to 1 - r dt, one step of discounting. So the
update is a lopsided weighted average, not a plain average.

The centre weight shrinks as dt grows. Once it goes negative, a cell is pushed
the wrong way from its old value, small errors flip sign every step and grow,
and the sheet turns into a checkerboard. It goes negative first at the high-S,
high-v corner, where A and C are largest. Setting it to zero there gives

```text
dt_stable = 1 / ( v_max S_max² / ds²  +  ξ² v_max / dv²  +  r )
```

which `report_stability` prints (§6.5). So dt must shrink like ds², and
doubling the resolution costs about four times the steps. On the smoke grid
the real cliff was 9% above the printed bound, so the bound is on the safe
side. Implicit or ADI schemes remove the limit but solve linear systems every
step.

### I2.7 Boundary conditions

| Edge | Condition | Code |
|---|---|---|
| t = T | payoff | `init_payoff` |
| S = 0 | call 0, put K e^(-rτ), exact | `apply_stock_edges` |
| S = S_max | call S e^(-qτ) - K e^(-rτ), put 0, approximate | `apply_stock_edges` |
| v = 0 | solve the reduced PDE on that row | `update_zero_row` |
| v = v_max | ∂V/∂v = 0, the top row copies the row below | `copy_top_row` |

The v = 0 row is the one interviewers ask about. All diffusion vanishes there,
leaving transport, `V_t + (r-q) S V_S + κθ V_v - r V = 0`. The variance drift
κθ is positive, so information flows into the grid and you solve that equation
rather than impose a value. The v derivative is one-sided upward, because a
central one needs a row at negative variance.

`copy_top_row` reads `next`, so it must run after the interior update has
finished. In every parallel rung that is why it comes after the barrier.

The Feller condition 2κθ ≥ ξ² says whether variance can reach zero. Ours is
0.12 against 0.1225, violated by a little, which is normal for equity
calibrations. The scheme stays valid, but the region near v = 0 is slightly
under-resolved on a uniform grid. A non-uniform grid is the fix, listed as
future work.

### I2.8 How we know the answer is right

Comparing with market prices would test the model, not the code. So the tests
compare with mathematical facts.

1. **Black-Scholes collapse** (`test_bs_collapse`). With xi = 0 and v0 = theta
   Heston becomes Black-Scholes, which has a formula. Error 1.16e-3 call,
   1.12e-3 put, tolerance 2.5e-3.
2. **Put-call parity** (`test_parity`). C - P = S e^(-qT) - K e^(-rT) holds in
   any model, and it exercises the cross term and v edges that test 1 turns
   off. Gap 5.4e-6 dollars, tolerance 1e-4.
3. **Fourier reference** (`make check-heston`). Heston has a semi-analytic
   price. On the reference grid, 2048 by 512 by 980000 steps, the PDE gives
   196.1684 against 196.1692, 4.2e-6 relative.
4. **Ladder match** (`test_ladder_matches`). Every rung against the baseline
   within 1e-11, and omp and mpi against simd with zero tolerance.

Tolerances were measured first and set just above the measurement, while a
real bug moves the price by dollars. No rung is benchmarked until it passes.

Tests 1 and 2 compare the solver with itself, so a shared offset passes both.
That happened once. The price was read at the nearest node, and on the smoke
grid v0 snapped 19% too high, giving +4.67%. The Fourier check caught it, and
`extract_result` now blends the four surrounding cells (§6.6, §6.7). On an
aligned grid the blend weight is exactly 0, so old answers stay bit-identical.

---

## I3. Memory and serial speed

### I3.0 Bytes, MiB and memory pages

Memory is counted in bytes. One byte is 8 bits, and one `double` takes 8
bytes. Bigger units go up in steps of 1024, because memory hardware is built
in powers of two, and 1024 is 2 to the power 10.

| Unit | Bytes | In this project |
|---|---|---|
| byte | 1 | one `double` is 8 bytes |
| cache line | 64 | 8 doubles, the amount the cache moves at once |
| KiB, kibibyte | 1024 | 4 KiB is one normal memory page |
| MiB, mebibyte | 1024 KiB = 1,048,576 | one sheet is 8 MiB |
| GiB, gibibyte | 1024 MiB | |
| TiB, tebibyte | 1024 GiB | keeping every sheet would take 7.5 TiB |

The i in KiB and MiB marks the 1024 version. KB and MB strictly mean 1000 and
1,000,000, so 1 MiB is about 5% more than 1 MB. Chip makers often write MB for
cache sizes when they mean MiB, so the 30 MB L3 cache on an r730 socket is 30
MiB.

The sizes to know by heart.

```text
one row     2048 doubles × 8 bytes          = 16,384 bytes    = 16 KiB
one sheet   2048 × 512 doubles × 8 bytes    = 8,388,608 bytes = 8 MiB
two sheets  current and next                                  = 16 MiB
L3 cache    one r730 socket                                   = 30 MiB, so both sheets fit
```

A memory page is the fixed-size chunk the operating system hands out memory
in. A program never sees real RAM addresses. Every address it uses is a
virtual address, and the operating system keeps a page table that says which
piece of real RAM each page of addresses lives in. On Linux on x86 a normal
page is 4 KiB. In Python terms, the page table is a dict from page number to a
place in RAM, and every load and every store has to look its page up.

```text
address the program uses  =  page number  +  position inside the page
                                  |              (last 12 bits, because 4 KiB = 2^12)
                                  v
              TLB, a small cache of recent lookups    hit   almost free
                                  | miss
                                  v
              page table in RAM, four levels deep     page walk, several extra memory reads
```

The page table itself sits in RAM, so a full lookup on every access would
make every load several loads. The core keeps the translations it used
recently in the TLB, the translation lookaside buffer. In Python terms it is
an `lru_cache` in front of the dict. On the r730 Haswell cores the
second-level TLB holds 1024 entries, one entry per page.

Because the TLB counts pages, not bytes, the page size decides how much memory
it covers. 1024 entries of 4 KiB cover only 4 MiB. 1024 entries of 2 MiB huge
pages cover 2 GiB.

| Page size | Pages needed for the 16 MiB of sheets | Fits in 1024 TLB entries? |
|---|---|---|
| 4 KiB, normal page | 4096 | no, so the core keeps doing page walks |
| 2 MiB, huge page | 8 | yes, with room to spare |

That is where the 1.66 times in §I3.5 comes from. The same binary on
the same grid ran at 134.8 Mcell/s on huge pages and 81.1 on 4 KiB pages.
Nothing in the code changed, only the page size the memory happened to get.
Linux gives a block huge pages by itself (transparent huge pages) when the
block comes fresh from the operating system through `mmap`.
`MALLOC_MMAP_THRESHOLD_=1048576` makes every block over 1 MiB come that way,
so every benchmark sets it.

Pages also decide where memory sits on a two-socket node. Each page is placed
on the socket of the thread that first writes to it, one page at a time. This
is first touch (§I5.5).

### I3.1 Two buffers

Every new cell must come from the old sheet only. Updating in place mixes
moments in time, because a cell's west and south neighbours are already new
while east and north are still old. On the test grid one step gave three
answers for one cell depending on loop order (50.0007471 with two buffers,
50.0007951 left to right, 50.0007465 right to left). With threads the answer
would depend on timing, which is a race.

So there are two sheets. During a step `current` is only read and `next` is
only written, then they swap. `std::swap` on two vectors swaps their internal
pointers and copies none of the million doubles (§5.4). Keeping every sheet
would take 980000 sheets of 8 MiB, about 7.5 TiB.

This rule is what makes the parallel rungs need no locks. Threads and ranks
read a sheet nobody is writing, and write only their own cells.

### I3.2 The flat array and the cache

A sheet is one flat `std::vector<double>`, not a list of lists.

```text
index(stock_i, var_j) = var_j * num_stock_nodes + stock_i

  row 2 |  8   9  10  11 |        memory  0 1 2 3 | 4 5 6 7 | 8 9 10 11
  row 1 |  4   5   6   7 |                 row 0  |  row 1  |  row 2
  row 0 |  0   1   2   3 |
```

A whole row is contiguous. East and west are one place apart, north and south
a whole row apart, which is why `row_above = row_centre + num_stock_nodes`
(§8.5). The inner loop walks along a row, so it walks memory in order.

RAM costs about 100 ns per access and cache about 1 to 10 ns. Data moves in
64-byte lines, eight doubles at a time, so walking in order makes seven of
every eight loads nearly free. SIMD also needs adjacent data. Swapping the
loops on purpose (the M1 `ctl-order` control) moved 4.8 times the bytes and
ran 2.1 times slower.

### I3.3 Compute-bound, not memory-bound

The cell update does 31 flops and, in the best case, moves 24 bytes (read one
new double of `current`, write one of `next`, plus the read of the line the
write lands in, called write-allocate). That is 1.29 flops per byte.

M1 measured the node's own limits with `slurm/bench_peak.sh`. The baseline
runs at the same speed whether its data is 1 MiB or 64 MiB, so memory does not
limit it. Its five divisions per cell do, at 13 to 20 cycles each, about the
whole budget of 66 cycles per cell. After the M1 optimisations the serial
kernel ran at 82 to 95% of the measured scalar peak while asking for only
about 3 GB/s. So it was compute-bound, and the next gain had to come from more
work per instruction, which is SIMD.

### I3.4 The serial optimisations

Milestone 1 tried six named techniques one at a time on rangpur, at -O2. The
M1 code is at git tag `milestone1`.

| Level | Technique | Mcell/s | Kept now as |
|---|---|---|---|
| baseline | | 52.6 | `solver_baseline.cpp` |
| L1 | hoisting, lookup tables | 52.7 | `StencilConstants`, S and S² tables |
| L2 | reciprocals instead of division | 131.2 | the `inv_*` fields |
| L3 | row-major loop order | 131.5 | `var_j` outside, `stock_i` inside |
| L4 | loop splitting | 131.3 | separate boundary functions |
| L5 | induction variables | 134.2 | row pointers |
| L6 | unrolling by four | 134.8 | dropped |

One rung did all the work. GCC at -O2 already did L1 and L3 to L6, because
they leave the answer unchanged. L2 changes the answer, since 1/(2ds) is
itself rounded, so the compiler is not allowed to do it without
`-ffast-math`. The trade is 2.5 times the speed for about 15 units in the last
place (3.3e-15 relative). The optimisations worth doing by hand are the ones
the compiler is not allowed to do.

The serial rung also copies every constant into a local variable (§8.4). The
loop writes through a `double*`, the compiler cannot prove that pointer never
points into the struct, and so it would reload each constant after every
store. This is pointer aliasing.

### I3.5 Page size and loop-order controls

Pages, huge pages and the TLB are explained in §I3.0.
The same binary ran 1.66 times faster in one memory regime than another. The
first 8 MiB allocation is served by `mmap` and gets 2 MiB huge pages, so 16
MiB needs about 8 TLB entries. Later allocations reuse the heap on 4 KiB
pages, needing 4096 entries and page walks. `MALLOC_MMAP_THRESHOLD_=1048576`
forces the fast regime for every rep, and every benchmark sets it. A real
solve is one fresh allocation, so the fast regime is the one reported.

Two controls showed what L3 and L4 were worth. Swapping the loops cost 2.1
times, smaller than its 4.8 times extra traffic, because 16 MiB fits in the
30 MB L3. Putting the v = 0 row back behind a per-cell `if` cost nothing,
because the branch predictor gets it right 510 times out of 511.

---

## I4. Benchmarking

### I4.1 The protocol

| Rule | Why |
|---|---|
| `sbatch` on an exclusive r730 node, never the login node | another job steals cores, cache and memory bandwidth |
| build clean and run every test first, refuse to time if one fails | a fast wrong answer is worth nothing |
| `MALLOC_MMAP_THRESHOLD_=1048576` | huge pages for every rep (§I3.5) |
| `OMP_PLACES=cores`, `OMP_PROC_BIND=close` | threads that move lose their cache |
| median of 5 reps, min and max kept | one rep can be unlucky |
| each rep runs for seconds, `maturity_for` sets dt at 70% of the limit | short runs are noise, unstable runs fill the grid with infinities |
| flush denormals with `LD_PRELOAD` | §I4.2 |
| each job builds in its own copy of the tree | two jobs once deleted each other's files |
| provenance header in every CSV | node, compiler, flags, binding and regime travel with the number |
| print the price | a compiler can delete a loop whose result is unused |

The rate is cell updates per second, `ns × nv × nt / seconds`, so grids of
different sizes compare. The node is a 2014 Haswell Xeon E5-2670 v3, two
sockets of 12 cores at 3.1 GHz, with AVX2 and FMA and no AVX-512. A job is
capped at 15 minutes, and a longer request waits in the queue forever rather
than being rejected.

### I4.2 Denormals

Below about 2.2e-308 a double becomes denormal, trading precision to reach
down to 4.9e-324. On this Xeon each operation on one takes a microcode path of
about a hundred cycles. Our grid makes them naturally, because value creeps
into the zero region of the payoff one cell per step, so far cells hold
1e-200, then 1e-310, then 0. How many cells are in that band depends on the
grid and the step count, so different runs measured different regimes.

The first size sweep found one grid 2.4 times slower than its neighbours, and
the standard grid lost 17%. The fix sets two bits of the MXCSR register before
`main` runs. FTZ writes tiny results as 0 and DAZ reads tiny inputs as 0.
`slurm/flush_denormals.c` does it and is loaded with `LD_PRELOAD`, so `heston`
and the tests are unchanged and the price is the same to ten digits.
`-ffast-math` sets the same bits. Experiment E12 measured the cost per rung,
1.08 to 1.38 times.

### I4.3 Profiling and -O levels

`-pg` makes the program note which function it is in every 10 ms, and
`gprof -b -p` prints a flat profile, the share of time per function. The
baseline spends 99.9% in `BaselineSolver::solve`, serial 99.4% in
`update_cells_scalar`, simd 98.9% in `update_row_avx`. The boundaries are
under half a percent, so the interior update is the only thing worth
speeding up.

| Level | baseline | serial |
|---|---|---|
| -O0 | 10.7 | 44.8 |
| -O1 | 57.6 | 156.1 |
| -O2 | 57.7 | 158.4 |
| -O3 | 112.5 | 238.6 |

Every level gave the -O2 answer exactly. -O3 adds the auto-vectoriser, here
two-lane SSE2. The project builds at -O2 so each rung measures one change of
ours, and the hand AVX2 rung (379) is still 1.59 times the -O3 serial rung.

---

## I5. The parallel rungs

### I5.1 Where the parallelism is

Inside one timestep every cell of the new sheet depends only on the old
sheet, so all million cells could be done at once. Across timesteps the work
is strictly in order, since sheet k+1 needs all of sheet k. So there is a lot
of parallelism inside a step and a hard barrier between steps. It is not
embarrassingly parallel like Monte Carlo.

```text
  SIMD, inside a row          OpenMP, rows to threads     MPI, slabs of rows to ranks
  [c c c c][c c c c]          row 7  thread 3             rank 1  rows 4..7, halo row 3
  4 cells per instruction     row 5  thread 2             rank 0  rows 0..3, halo row 4
                              row 1  thread 0
```

The ladder is one technique per rung, each tested against the rung below
before it is timed. baseline, serial, simd (AVX2 and FMA), omp (threads over
rows), mpi (ranks over slabs, OpenMP inside each).

### I5.2 SIMD basics

SIMD is single instruction, multiple data. One instruction works on several
numbers held side by side in a wide register. The r730s have AVX2, 256-bit
`ymm` registers holding four doubles, and FMA. They do not have AVX-512 and
its `zmm` registers.

| Intrinsic | Does in every lane | numpy on 4 elements |
|---|---|---|
| `_mm256_set1_pd(x)` | fill all lanes with x | `np.full(4, x)` |
| `_mm256_loadu_pd(p)` | load p[0..3] | `a = p[i:i+4]` |
| `_mm256_storeu_pd(p, a)` | store four lanes | `p[i:i+4] = a` |
| `_mm256_add_pd`, `_sub_pd`, `_mul_pd` | +, -, × | `a + b` |
| `_mm256_fmadd_pd(a, b, c)` | a × b + c, rounded once | `a * b + c` |
| `_mm256_fnmadd_pd(a, b, c)` | c - a × b, rounded once | `c - a * b` |

`__m256d` is the type of a four-double register. The names read `_mm256`,
then the operation, then `_pd` for packed doubles. The row vectorises
cleanly because every cell runs the same formula with no `if`, four
neighbouring cells are one 32-byte load, and no cell reads another's new
value. The per-row coefficients go into registers once with `set1`, and S and
S² are loaded four at a time from the tables.

SIMD happens inside one core, and threads multiply it, so 24 cores doing four
cells each is 96 cells in flight.

### I5.3 SIMD details

Four cells starting at column i need nine registers, the same four-wide
window shifted one column or one row. Lane 0 of `east` is the east neighbour
of lane 0 of `centre`, so after the loads the code is the scalar formula with
intrinsics (§9a.4).

The windows start at i-1, i and i+1, 8 bytes apart, so at most one can be
32-byte aligned and an aligned load would crash on the others. So every load
is `loadu`. Aligned buffers (E2) were not run.

FMA does a multiply and an add as one instruction with one rounding. Eight of
the 24 arithmetic instructions per four cells are FMAs. It changes the last
bits, so simd matches serial to about 1e-15, not bit for bit. With `-mfma` g++
may fuse any `a*b + c` in any file, so only `stencil_avx.cpp` gets `-mavx2
-mfma` and the other rungs compile as before (§13.5). Not `-march=native`,
because a build on a newer machine could use AVX-512 and crash on the r730s.
`objdump -d heston | grep -c zmm` printing 0 proves there is none.

The last zero to three cells of a row go through `update_cells_scalar`, the
serial code (§8.10). The Mac has no AVX2, so the file is inside
`#if defined(__AVX2__) && defined(__FMA__)` with a scalar fallback, and
`make avx-check` runs an x86_64 build under Rosetta. Rosetta and rangpur gave
the same numbers.

### I5.4 OpenMP threads

A thread is a stream of instructions inside a process, and all threads of a
process share its memory, while each has its own stack for locals. Python's
GIL lets one thread run Python at a time, C++ has no such lock. OpenMP is a
set of pragmas plus a runtime that creates and manages the threads (§9b).

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
```

| Part | Meaning |
|---|---|
| `parallel` | start a team of threads |
| `for` | split the loop's iterations across the team |
| `schedule(static)` | one equal contiguous band per thread, decided up front |
| `default(none)` | every variable's sharing must be stated, or it will not compile |
| `shared(...)` | the same object for every thread, `var_j` is private automatically |

Each step forks the team, each thread does its band of rows, and an implicit
barrier at the end waits for the last one. Then the main thread does the
boundaries and the swap. The runtime keeps threads alive between loops, so a
fork costs microseconds, not a new thread.

Rows suit it because it is the outer loop (one fork per step), each band is
contiguous memory, and the AVX2 inner loop is untouched. 510 rows on 4
threads split 128, 128, 127, 127.

### I5.5 Why OpenMP needs no locks, and where threads sit

A data race is two threads touching the same memory at once with at least one
writing. Here every thread reads `current`, which nobody writes in the step,
and writes only its own rows of `next`. So there is no `critical`, `atomic` or
lock. There is also no reduction, so each cell runs exactly the simd
instructions on a different core, and omp matches simd bit for bit at any
thread count. The test checks that with zero tolerance at 7 threads, where the
bands are uneven.

More threads than cores is oversubscription, and with a barrier every step
the team waits for whichever thread was paused, so it is slower.

Each r730 socket has its own L3 and memory. A page lands on the socket of the
thread that first writes it (first touch). `Grid` zero-fills both sheets on
the main thread, so every page starts on socket 0. On this node even CPUs are
socket 0 and odd are socket 1, so `close` alternates sockets, measured with
`slurm/diag_placement.sh`. The triad bandwidth benchmark in `bench_peak.sh`
had the same first-touch problem and was fixed with `new double[n]` in a
`unique_ptr` (§15.12).

### I5.6 MPI ranks and halos

MPI runs separate processes called ranks, each with its own memory, possibly
on another machine. They share data only by messages. `mpiexec -n 4 ./heston`
starts four copies. `MPI_Init`, `MPI_Comm_rank` (which am I), `MPI_Comm_size`
(how many) and `MPI_Finalize` frame every run. Started without `mpiexec` it is
a group of one, which is how the other rungs run in the same binary (§3.6).

Each rank owns a slab of whole rows, the first ranks one extra when they do
not divide (§9c.2). The code calls it `owned_rows` and the received halo rows
`row_below` and `row_above`. Rank 0 owns the v = 0 row and the last rank the top row.
A row is contiguous, so it goes as one message straight from the array. A
column would need copying into a buffer first.

The stencil reads one row above and below, so a slab's edge rows need one row
from each neighbour every step. Those copies are halo rows, one row deep. A
row of 2048 doubles is 16 KiB, so a middle rank moves 64 KiB per step.

```text
  each step on every rank
    1. own interior rows with OpenMP      reads current, writes next
    2. own boundaries
    3. exchange halo rows of next         MPI_Irecv, MPI_Isend, MPI_Waitall
    4. swap_buffers
```

Blocking `MPI_Send` from both neighbours at once can deadlock, and small
messages hide the bug. So the rung posts non-blocking receives first (the
data then lands straight in place), then non-blocking sends, then waits for
all four.

### I5.7 Hybrid MPI and the gather

Inside each rank the OpenMP loop runs over its own rows, so the rung is MPI
between ranks and threads within. One rank per socket (2 ranks of 12 threads
on a node) keeps each rank's memory on its own socket. `--bind-to none` stops
Open MPI pinning a rank to one core, which would squeeze all its threads onto
it.

Every rank keeps both full sheets and updates only its rows, so every row has
the same index everywhere and every row function is reused unchanged. The
cost is 16 MiB per rank. A grid too big for one machine would give each rank
only its slab plus two halo rows, with indices shifted.

After the loop `MPI_Gatherv` collects every slab into rank 0's spare sheet.
`Gatherv` takes a count and offset per rank because slabs can differ by a
row. Rank 0 then reads the price as usual. Halo rows are exact copies and
nothing is summed across ranks, so mpi matches simd bit for bit. With the
halo exchange removed on purpose the answer moved by 2e-9, and the zero
tolerance test failed as it should.

---

## I6. Results and what they mean

### I6.1 The ladder

On r730-2, 2048 by 512, Mcell-updates/s.

| Rung | Rate | vs baseline |
|---|---|---|
| baseline | 57 | 1 |
| serial | 159 | 2.8 |
| simd | 379 | 6.7 |
| omp, 24 threads | 5729 | 101 |
| mpi, 2 ranks of 12, one node | 6019 | 106 |
| mpi, two nodes | 8768 | 155 |

SIMD gave 2.4 times, not 4. Four lanes is the ceiling, the nine windows
overlap, the remainder cells are scalar, and the serial rung already kept the
arithmetic ports busy.

24 threads gave 15.8 times, not 24. The serial boundaries grow from 0.4% of a
one-thread step to 8.2% of a 24-thread step, because the main thread reads
rows other cores just wrote. The fork and barrier cost the same every step
however short the step. And at 24 threads the AVX2 kernel moves about 138
GB/s through the caches, so shared bandwidth starts to run out. The scalar
omp rung, at half the data rate, scales 18.5 times.

### I6.2 Amdahl and Gustafson

Strong scaling keeps the grid fixed and adds threads. Speedup is
S(p) = T(1)/T(p) and efficiency S(p)/p. Amdahl's law caps it with the serial
fraction f.

```text
S(p) = 1 / ( f + (1 - f)/p )        rearranged   1/S - 1/p = f (1 - 1/p)
```

The second form is a straight line through the origin in f, so f is the least
squares slope `sum(x*y) / sum(x*x)`, with no fitting library. The fit gives
f = 1.2%, three times the boundary share, so most of the loss is per-step
overhead and bandwidth, not the boundaries.

Weak scaling keeps 21 rows per thread and grows the grid with p, stretching
`--vmax` so the row spacing stays the same. Gustafson's law gives the scaled
speedup S = f + (1 - f) p, which has no cap, where f is the serial share of
the p-thread run (here `boundary_seconds / seconds`, 8.4%). That predicts 22
at 24 threads against 19.1 measured. The two laws use different f, so the
strong fit is not reused. The one-thread weak point has 23 rows, which fit in
L2, so it looks a little better than the rest. The course does not teach
Gustafson, so the slide cites Gustafson 1988.

### I6.3 Size, memory and nodes

512 by 128 on 24 threads is only 8 times serial, because each step is too
short to hide the fork and barrier. 2048 by 512 fits in the two L3 caches and
gets 36 times. 4096 by 1024 does not fit, sits on the measured DRAM roof and
drops to 19 times.

The roofline plots GFLOP/s against flops per byte (1.29 here). The flat roofs
are the measured scalar peak, the AVX2 FMA peak (43.4 GFLOP/s on one core,
951 on 24) and the sloped roofs are measured triad bandwidth. Every rung sits
on the same vertical line and rises as it gets faster. A point can sit above
the DRAM roof when the grid fits in cache.

Page placement mattered on the big grid. `numactl --interleave=all` gave 28%
more, and one MPI rank per socket, each first-touching its own grid, gave 2.16
times one 24-thread rank. Two nodes gave 1.41 to 1.57 times one node, not 2,
because a step is only 0.12 ms and each halo exchange costs about 45 µs of
latency and waiting.

Three things outside the code moved the numbers more than some
optimisations. Huge pages 1.5 to 1.8 times (E11), denormals 1.08 to 1.38
times (E12), page placement up to 2.16 times. That is why the protocol pins
all three.

### I6.4 Tried and rejected

| Experiment | Result | Decision |
|---|---|---|
| dynamic and guided schedules (E5) | 16 to 40% slower | rejected, every row costs the same |
| one parallel region per solve (E4) | 3.5% faster, 15% on small grids | not in the ladder, one technique per rung |
| boundaries in the parallel loop (E8) | 5% faster | not in the ladder, same reason |
| float instead of double (E3) | 1.92 times on 1 thread, 4.46 times at 4096 by 1024 | rejected for accuracy |
| aligned loads (E2), halo overlap (E9) | not run | E2 needs an aligned allocator, E9 is on the cut list |

Float is 4 bytes and about 7 digits, against 8 bytes and 16. A register
holds eight floats, and the 4096 by 1024 float sheets fit back in L3, hence
4.46 times. But every cell is rounded every step, and on the 980000 step
reference grid the price was $0.044 from the Fourier price (double is
$0.0008) and gamma moved 11.9%. Gamma suffers most because it subtracts nearly
equal numbers and divides by ds². That is worse than a far cheaper double
grid, so float throws away the accuracy the fine grid paid for.

### I6.5 Limitations and extensions

- The uniform grid under-resolves the region near v = 0, where Feller is
  violated. A non-uniform v grid would fix it.
- The explicit scheme needs dt ∝ ds², so fine grids are expensive. Implicit
  or ADI schemes remove that limit at the cost of linear solves.
- Cutting the grid off at S_max and v_max = 0.64 adds about 2e-4 dollars,
  three orders below the Fourier agreement.
- Every MPI rank holds the whole grid. Past one machine's memory each rank
  would hold only its slab and halos.
- The halo exchange is not overlapped with compute (E9), which would help the
  two-node runs most.
- CUDA was a bonus and was not built. It would use one GPU thread per cell
  with both sheets kept on the GPU.

---

## I7. The C++ in the code

Each card lists the standard library and C++ features of one kind. Every row
has a real line from the code, what it does and why it is used there. `std::`
means the name comes from the C++ standard library, which ships with every
compiler the way Python's built-ins and standard modules do. The `#include`
at the top of a file says which part of the library the file uses.

### I7.1 std::vector, the container for every sheet and table

`#include <vector>`. A `std::vector<double>` is a resizable array that owns
its memory and keeps every element side by side, like a Python list that
holds only doubles and frees itself.

| Line | What it does | Why here |
|---|---|---|
| `std::vector<double> current_;` (grid.h) | declares an empty vector member | the sheet must be one contiguous block for the cache and SIMD (§I3.2) |
| `current_.assign(total_cells, 0.0);` (grid.cpp) | sets the size to `total_cells` and fills every element with 0 | the sheet size is only known once the config is read, and zeros give a defined start |
| `stock_price.resize(num_stock_nodes);` (stencil.cpp) | grows the vector to that many elements | the S table is then filled in a loop, one entry per column |
| `const double* row_centre = current.data() + var_j * num_stock_nodes;` | `.data()` gives a raw pointer to the first element | the hot loop indexes rows through plain pointers, which the compiler handles well. The vector still owns the memory |
| `current[grid.index(stock_i, var_j)]` | `[]` reads one element with no bounds check | `.at()` would check every access and cost time in the hot loop. The tests catch index mistakes instead |
| `std::vector<MPI_Request> requests;` then `requests.push_back(receive_below);` | starts empty and appends one item, like `list.append` | how many messages there are depends on whether the rank has neighbours above and below |
| `MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);` | `.size()` is the count, `.data()` the pointer | MPI is a C library, so it takes a count and a pointer, not a vector |
| `current[cell] = (float)grid.current()[cell];` (solver_float.cpp) | copies one cell, rounding the double to a float | a plain loop fills the float copy of the payoff sheet |
| `stock_price[stock_i] = (float)constants.stock_price[stock_i];` (stencil.cpp) | copies one table entry, rounded to float | the float constants start from the double ones, rounded once |
| `std::swap(current_, next_);` (grid.cpp) | `#include <utility>`. Exchanges the two vectors' internal pointers | ends each step without copying a million doubles (§5.4) |

### I7.2 Ownership, references and pointers

| Line | What it does | Why here |
|---|---|---|
| `std::unique_ptr<Solver> solver = make_solver(cfg.solver);` (main.cpp) | `#include <memory>`. A pointer that owns its object and deletes it when it goes out of scope | `main` holds a solver chosen at run time without knowing its class, and nothing has to remember to free it |
| `return std::make_unique<SerialSolver>();` (solver.cpp) | creates the object and its owner in one step | this avoids writing `new`, which the project never does in its own C++ |
| `virtual ~Solver() = default;` (solver.h) | a virtual destructor with the compiler's default body | deleting a `SerialSolver` through a `Solver` pointer must run the subclass's destructor too |
| `std::unique_ptr<double[]> a(new double[n]), b(new double[n]), c(new double[n]);` (bench_peak.sh) | owns a plain array. `new double[n]` does not write the memory | a vector would zero every page on the main thread and put it all on socket 0. Here the parallel loop writes it first (§I5.5) |
| `const Config& cfg` (every `solve`) | a reference, another name for the caller's object. `const` means it cannot be changed | passing by reference copies nothing, and `const` says the solver only reads the settings |
| `const std::vector<double>& current = grid.current();` | a read-only name for the grid's sheet | the compiler refuses any write to `current` during a step, which is the two-buffer rule (§I3.1) |
| `double* out_row = next.data() + var_j * num_stock_nodes;` | a pointer that looks into memory it does not own | pointer arithmetic gives the start of a row once, then the loop indexes along it |

RAII means a resource is taken when an object is made and given back when it
is destroyed. Vectors, `unique_ptr` and `ifstream` all work this way, which
is why there is no `delete` anywhere. The one `fclose` is in io.cpp, because
`fopen` is plain C and C files are closed by hand.

### I7.3 Strings and parsing

`#include <string>`. `std::string` is a text value that owns its characters,
and `+` joins strings as in Python.

| Line | What it does | Why here |
|---|---|---|
| `std::istringstream words(line);` (params.cpp) | `#include <sstream>`. Treats one line of text as a stream you can read words from | lets the next line split the line into words |
| `words >> key >> equals >> value;` | reads the first three words, separated by spaces | like `key, equals, value = line.split()[:3]` in Python. A missing word stays `""` |
| `if (key.empty() \|\| key[0] == '#') {` | true for a blank line or a line whose first word starts with `#` | blank lines and comment lines are skipped |
| `cfg.option.strike = std::stod(value);` | turns text into a double, like `float("5250.0")` | config values arrive as text. Text that starts with a number, such as `12x`, is read as 12 |
| `} catch (const std::invalid_argument&) {` | catches the error `stod` and `stoi` throw for text that is not a number | turns it into a message that names the file and line (§4.3) |
| `cfg.grid.num_stock_nodes = std::stoi(value);` (main.cpp) | turns text into an int, like `int("1024")` | command line values arrive as text |
| `for (int i = 3; i < argc; i = i + 2) {` (main.cpp) | walks the command line two words at a time | every flag is followed by exactly one value (§3.2) |
| `std::string where = path + ":" + std::to_string(line_number);` | turns a number into text | error messages name the file and line, like `smoke.cfg:12` |
| `std::string(argv[1]) != "--config"` (main.cpp) | makes a `std::string` from a C string | `argv` holds C strings (`char*`), and comparing with `!=` needs a `std::string` on one side |
| `path.c_str()` (io.cpp) | gives the C string inside a `std::string` | `fopen` and `printf` are C functions and need a `char*` |

### I7.4 Files and printing

| Line | What it does | Why here |
|---|---|---|
| `std::ifstream file(path);` (params.cpp) | `#include <fstream>`. Opens a file for reading | the config file. It closes itself at the end of the function (RAII) |
| `while (std::getline(file, line)) {` | reads one line into `line`, and is false at the end of the file | the config is parsed line by line, like `for line in f` in Python |
| `FILE* out = std::fopen(path.c_str(), "w");` (io.cpp) | the C way to open a file for writing. It gives `NULL` if it fails | snapshots for the weather map figures |
| `std::fprintf(out, "%.10g", ...);` then `std::fclose(out);` | writes to the file with 10 significant digits, then closes it | `%g` with 6 digits, the default, would lose detail in the figures |
| `std::filesystem::create_directories(dump_dir);` | `#include <filesystem>`, C++17. Makes the folder and any parents, and does nothing if it exists | like `os.makedirs(exist_ok=True)`. GCC 8 needs `-lstdc++fs`, which the Makefile adds only when the compiler wants it |
| `std::snprintf(file_name, sizeof file_name, "snap_%06d.csv", step);` | `#include <cstdio>`. Formats into a fixed `char` buffer, never writing past `sizeof file_name` | gives zero-padded names like `snap_000100.csv` that sort in time order |
| `std::printf("%.10g,%.8g,%.8g,...\n", result.price, ...);` (io.cpp) | prints a formatted line to stdout | one exact CSV line per solve that the scripts parse. `%.10g` fixes the digits |
| `std::fprintf(stderr, "[%s] feller %s ...", ...);` (solver.cpp) | prints to stderr instead of stdout | stdout stays one clean CSV line, while messages go to the other stream |

`printf` is used rather than `std::cout` because a format string pins the
digits of every column in one place, and the CSV line reads as one line of
code.

### I7.5 Errors

| Line | What it does | Why here |
|---|---|---|
| `throw std::runtime_error(flag + " needs a value");` (main.cpp) | `#include <stdexcept>`. Stops the function and sends an error object up the call chain, like `raise` | a bad flag or config line must stop the run, not be ignored |
| `try {` ... `} catch (const std::invalid_argument&) {` (params.cpp) | catches the error `stod` and `stoi` throw in the `try` block | `stod` throws its own error for text like `abc`, and the code replaces it with one that names the file and line |
| `} catch (const std::exception& e) {` (main.cpp) | catches every error once, at the top | the program prints the message and returns 1 instead of crashing |
| `std::fprintf(stderr, "heston: %s\n", e.what());` | `.what()` is the error's message | shows the reason to the user and in the Slurm log |

The work sits in `run()` so that `main` still calls `MPI_Finalize` after an
error (§3.6). Catching by `const&` avoids copying the error and keeps the
subclass's message.

### I7.6 Maths functions

`#include <cmath>` for the maths, `#include <algorithm>` for `max` and
`min`.

| Line | What it does | Why here |
|---|---|---|
| `payoff = std::max(price - opt.strike, 0.0);` (grid.cpp) | the larger of two values | the payoff max(S - K, 0). Both arguments must be the same type, so it is `0.0` and not `0` |
| `int first_interior = std::max(owned_rows.begin, 1);` (solver_mpi.cpp) | the larger of two ints | rank 0's slab starts at row 0, which is not an interior row |
| `int end_interior = std::min(owned_rows.end, num_variance_nodes - 1);` | the smaller of two ints | the last rank's slab includes the top row, which is not interior either |
| `if (position < 0.0) { position = 0.0; }` (grid.cpp) | keeps a value inside a range with two `if`s, like `max(0, min(x, top))` in Python | a spot price off the grid still gives a valid position |
| `if (out.cell < 1) { out.cell = 1; }` and `if (out.cell > highest_cell) {...}` | the same for the cell | the blend and the gamma stencil read one cell further out on each side, so the cell stays one node in from each edge (§5.6) |
| `if (out.weight > 1.0) { out.weight = 1.0; }` | keeps the blend weight between 0 and 1 | once the cell is moved in from an edge, the raw fraction could fall outside one cell |
| `double nearest_node = std::round(position);` | rounds to the nearest whole number | to test whether the point sits on a node |
| `bool on_node = std::fabs(position - nearest_node) < node_snap_tolerance;` | absolute value of a double | a point within 1e-9 cells of a node counts as on it, so aligned grids give exactly the node value |
| `int lower_node = (int)std::floor(position);` | rounds down, then `(int)` makes it a whole number | the node just below the point |
| `double strike_today = constants.strike * std::exp(-constants.rate * time_to_expiry);` (stencil.cpp) | e to the power x | discounts the strike back from expiry for the edge values |
| `return 0.5 * std::erfc(-x / std::sqrt(2.0));` (black_scholes.cpp) | `erfc` is the complementary error function, `sqrt` the square root | N(x) = ½ erfc(-x/√2) is exact, so the normal CDF needs no approximation |
| `std::log(spot / strike)` (black_scholes.cpp) | natural log | part of d1 in the Black-Scholes formula |
| `worst = std::fmax(worst, std::fabs(next_scalar[cell] - next_avx[cell]));` (test_avx_row.cpp) | the larger of two doubles | keeps the largest difference seen so far |

### I7.7 Timing, and the random numbers in the tests

`#include <omp.h>` for the clock, `#include <random>` for the test's noise.

| Line | What it does | Why here |
|---|---|---|
| `double step_start = omp_get_wtime();` | OpenMP's wall clock, in seconds as a `double`, like `time.perf_counter()` | every solver uses it, the serial ones too, so all are timed the same way. On Linux it reads the clock that never jumps when the time is synced |
| `result.seconds += step_end - step_start;` | one subtraction gives the seconds for this step | adds up each step's compute time. The timer stops before snapshots, so file writing is never counted |
| `std::mt19937 random_engine(3500);` (test_avx_row.cpp) | a random number generator started from a fixed seed | the same seed gives the same numbers every run, so a failure can be repeated |
| `std::uniform_real_distribution<double> noise(-1.0, 1.0);` then `noise(random_engine)` | turns the generator's output into doubles between -1 and 1 | noise on the payoff makes every cell different from its neighbours, so a wrong neighbour in the AVX row shows up |

Random numbers appear only in this test. The solver itself has none (§I2.2).

### I7.8 Classes, structs and other C++ features

| Line | What it does | Why here |
|---|---|---|
| `struct AxisPosition {` `int cell = 0;` `double weight = 0.0;` `};` (grid.h) | a plain record with default values, like a dataclass | groups the two results `locate` returns |
| `int index(int stock_i, int var_j) const { return var_j * num_stock_nodes_ + stock_i; }` (grid.h) | a method written inside the class. `const` promises it does not change the grid | written in the header so the compiler can inline it, which makes the call free |
| `virtual SolveResult solve(const Config& cfg) = 0;` (solver.h) | `virtual` means a call through a `Solver` pointer runs the subclass's version. `= 0` means every subclass must write one | `main` calls `solve` without knowing which rung it has, like an abstract method in Python |
| `SolveResult solve(const Config& cfg) override;` | says this replaces a base class method | a typo in the signature becomes a compile error instead of a new method |
| `static const double node_snap_tolerance = 1e-9;` (grid.cpp) | a constant visible only in this file | a named value instead of a bare 1e-9 in the code |
| `static void apply(...)` (params.cpp) | `static` on a function means only this file can call it | a helper that is not part of the file's public interface |
| `double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;` | `(double)` converts the first int, so the whole product is done in double | 2048 × 512 × 980000 is about 1e12, far past the 2.1 billion an `int` holds |
| `lower_node = (int)nearest_node;` (grid.cpp) | `(int)` drops the fraction of a double | `floor` and `round` return doubles, and a cell index must be an int |
| `dt = (float)constants.dt;` (stencil.cpp) | rounds a double to the nearest float | the float solver starts from the same constants, rounded once |
| `for (int step = 1; step <= num_timesteps; ++step) {` | `++step` adds one | the C++ way to write `step += 1` in a loop |
| `#if defined(__AVX2__) && defined(__FMA__)` | the preprocessor keeps one branch of the file before compiling | the AVX2 code only exists when the compiler was given the AVX2 flags (§9a.1) |
| `#pragma omp parallel for ...` | an instruction to the compiler, ignored if OpenMP is off | turns the row loop into a threaded loop (§I5.4) |

`std::size_t` has its own explanation in §2.2.

### I7.9 OpenMP, MPI and intrinsic calls

These are not `std::`. They come from libraries the Makefile links in.

| Call | From | Why here |
|---|---|---|
| `omp_get_max_threads()` | `<omp.h>` | records how many threads the run really got, which goes in the CSV |
| `MPI_Init(&argc, &argv)`, `MPI_Finalize()` | `<mpi.h>` | start and stop MPI once in `main`, so every rung runs in the same binary |
| `MPI_Comm_rank`, `MPI_Comm_size` | `<mpi.h>` | which rank this is and how many there are, which decides the slab |
| `MPI_Irecv`, `MPI_Isend`, `MPI_Waitall` | `<mpi.h>` | the halo exchange without deadlock (§I5.6) |
| `MPI_Gatherv` | `<mpi.h>` | collects slabs of different sizes into rank 0 (§I5.7) |
| `_mm256_loadu_pd` and the other `_mm256_` calls | `<immintrin.h>` | each becomes one AVX2 instruction on four doubles (§I5.2) |
| `_mm_setcsr(_mm_getcsr() \| 0x8040)` | `<xmmintrin.h>` | sets the FTZ and DAZ bits in `flush_denormals.c` (§I4.2) |

---

## Drill

- **What is the payoff of a put at expiry?** max(K - S, 0). (§I1.1)
- **Why does Heston exist when Black-Scholes has a formula?** Volatility moves in real markets, and Heston makes the variance random with its own mean reversion. (§I1.3)
- **What does each Heston parameter control, and why is rho negative?** v0 today's variance, kappa the speed back to normal, theta the normal level, xi how hard variance is kicked, rho the correlation. Crashes and volatility spikes come together, so rho is negative. (§I1.3)
- **Where is the strike in the PDE?** Nowhere. It is in the payoff and the S edge values, so only those change between a call and a put. (§I2.2)
- **Where did dW1 go, and why is the coefficient v and not sqrt(v)?** The PDE is the equation the average over all paths satisfies (Feynman-Kac). Averaging squares the noise size, so sqrt(v) S becomes v S² in the diffusion term. (§I2.2)
- **Why does the cross term have no ½?** The mixed partial appears twice in Itô's lemma, so the ½ cancels. (§I2.2)
- **Why solve backwards from expiry?** That is the only sheet whose values are known. (§I2.3)
- **Why nine points and not five?** The correlation term needs the cross derivative, which reads the four corners. (§I2.4)
- **Work out the strike cell after one step on the test grid.** d2V_dS2 = 0.02 and dV_dS = 0.5 give change 11025 + 84 = 11109, times dt 4.464e-6 is 0.04959. (§I2.5)
- **Is the nine-point update an average of the nine cells?** No. It is a weighted sum with very uneven weights that add up to 1 - r dt. (§I2.6)
- **What makes the explicit scheme blow up, and where does it start?** A dt big enough to make the centre weight negative. It starts at the high-S, high-v corner, where the diffusion coefficients are largest. (§I2.6)
- **Why a forward difference in v on the v = 0 row?** A central difference needs a row at negative variance, and the drift κθ points into the grid, so information comes from above. (§I2.7)
- **What is the Feller condition, and do our parameters meet it?** 2κθ ≥ ξ². 0.12 against 0.1225, so it is violated by a little, which is normal and leaves v near 0 slightly under-resolved. (§I2.7)
- **How do you know the price is right without market data?** Black-Scholes collapse, put-call parity, and the Fourier price, which agrees to 4.2e-6 on the reference grid. (§I2.8)
- **Why was the external Fourier check needed when the other tests passed?** Self-consistency checks cannot see a shared offset. The nearest-node readout was 4.67% off and both tests passed. (§I2.8)
- **What goes wrong with one buffer?** Some neighbours are already this step's values, so the update mixes two moments in time and depends on loop order, or on thread timing. (§I3.1)
- **What does the buffer swap cost?** Two pointer swaps. No doubles are copied. (§I3.1)
- **Why is a flat array with a contiguous inner loop faster?** Cache lines bring in eight neighbouring doubles at once, so walking in order makes most loads nearly free, and SIMD needs adjacent data. (§I3.2)
- **Is the solver memory-bound or compute-bound?** The baseline is limited by division, and the optimised serial kernel ran at 82 to 95% of the measured scalar peak, so compute-bound. (§I3.3)
- **Which serial optimisation won, and why that one?** Reciprocals instead of division, 2.5 times. GCC already did the others at -O2, and it is not allowed to make this one because it changes the rounding. (§I3.4)
- **Why copy the constants into locals before the loop?** Pointer aliasing. The compiler cannot prove the output pointer does not point into the struct, so it would reload constants after every store. (§I3.4)
- **Why did reps 2 to 5 run 1.66 times slower than rep 1?** Rep 1 got huge pages through mmap and later reps got 4 KiB heap pages and TLB misses. The benchmarks pin the fast regime. (§I3.5)
- **Why benchmark on an exclusive compute node, with a median of 5?** The login node is shared, another job steals cache and bandwidth even on other cores, and one rep can be unlucky. (§I4.1)
- **What is a denormal and why did it matter?** A double below about 2.2e-308. Each operation on one is far slower, and the number of them depended on the run, so runs measured different regimes. FTZ and DAZ turn them into 0. (§I4.2)
- **Why build at -O2 when -O3 is faster?** So each rung measures one change of ours and not the auto-vectoriser. The hand AVX2 rung is still 1.6 times the -O3 serial rung. (§I4.3)
- **Where is the parallelism, and what limits it?** Every cell inside a step is independent, but steps are strictly in order, so there is a barrier every step and Amdahl's law applies. (§I5.1)
- **What is a __m256d?** A 256-bit register holding four doubles. One AVX2 instruction works on all four. (§I5.2)
- **Why loadu and not load?** The three windows start 8 bytes apart, so at most one can be 32-byte aligned and the aligned load crashes on the others. (§I5.3)
- **Why do only some files get -mfma?** With -mfma g++ may fuse multiply-adds in any file, which would change the baseline and serial answers. (§I5.3)
- **Why not -march=native?** It could enable AVX-512 on a newer build machine, which would crash on the r730s. objdump shows no zmm. (§I5.3)
- **What does each word of the OpenMP pragma do?** parallel starts a team, for splits the loop, schedule(static) gives equal bands up front, default(none) makes every variable's sharing explicit, shared lists the common ones. (§I5.4)
- **Why rows and not columns for threads?** It is the outer loop, so one fork per step, each thread gets contiguous memory, and the AVX2 inner loop is kept. (§I5.4)
- **Why are no locks needed, and why is omp bit-identical to simd?** Threads only read current and each writes its own rows of next. There is no reduction, so each cell runs the same instructions. (§I5.5)
- **What is first touch?** A page goes to the socket of the thread that first writes it. Grid fills both sheets on the main thread, so every page starts on socket 0. (§I5.5)
- **MPI or OpenMP, what is the difference?** Threads share one process's memory. Ranks are separate processes with separate memory and share data only by messages. (§I5.6)
- **Why is the halo one row deep, and why rows for slabs?** The stencil reaches one row each way, and a row is contiguous, so it is one message with no copying. (§I5.6)
- **Why Irecv before Isend, and why non-blocking?** Two blocking sends to each other can deadlock. Posting the receive first lets the data land in place. (§I5.6)
- **Why does every rank hold the whole grid?** Every row keeps its index, so the row functions and extract_result are reused. It costs 16 MiB per rank. (§I5.7)
- **Why Gatherv and not Gather?** Slabs can differ by a row, and Gatherv takes a count and offset per rank. (§I5.7)
- **Why did SIMD give 2.4 times and not 4?** Overlapping loads, the scalar remainder, and a serial rung that already kept the ports busy. (§I6.1)
- **Why 15.8 times on 24 threads?** The serial boundaries grow to 8% of a step, the fork and barrier cost the same every step, and the AVX2 kernel starts to share out the bandwidth. (§I6.1)
- **How is Amdahl's f fitted with no library?** 1/S - 1/p = f (1 - 1/p) is a line through the origin, so f is sum(x*y) / sum(x*x). (§I6.2)
- **Why are Amdahl's and Gustafson's f different numbers?** Amdahl's is the serial share of the one-thread run, Gustafson's the serial share of the many-thread run. (§I6.2)
- **Why does one rank per socket beat one rank of 24 threads?** Each rank first-touches its own grid on its own socket, while the single process has every page on socket 0. (§I6.3)
- **Why do two nodes not give twice the speed?** Each step ends with a halo exchange of about 45 µs, against a step of 0.12 ms. (§I6.3)
- **Why was float rejected when it was faster?** Rounding every cell every step put the reference price $0.044 off and moved gamma 12%, worse than a cheap double grid. (§I6.4)
- **What would you do with more time?** A non-uniform v grid, an implicit or ADI scheme, slab-only memory per rank, overlapping the halo exchange, and CUDA. (§I6.5)
- **Why does the code write std::max(price - opt.strike, 0.0) and not 0?** std::max needs both arguments to be the same type, and price - strike is a double. (§I7.6)
- **What do the three pairs of ifs in locate do?** Keep the position on the grid, keep the cell one node in from each edge so the Greek stencils fit, and keep the blend weight between 0 and 1. (§I7.6)
- **How is a step timed?** `omp_get_wtime()` before and after, and the difference is added to `seconds`. It returns seconds as a double, so there is nothing to convert. (§I7.7)
- **Why is there no delete anywhere?** RAII. vector, unique_ptr and ifstream free or close themselves when they go out of scope. The one `fclose` is for the C `fopen` in io.cpp. (§I7.2)
- **Why does stod take &used?** stod stops at the first bad character, so the code checks it used every character and rejects text like 1.5abc. (§I7.3)
- **Why printf and not cout for the result line?** The format string pins the digits of every column in one place. stdout carries only that line, and messages go to stderr. (§I7.4)
- **Why is the cell count cast to double before multiplying?** ns times nv times nt is about 1e12 on the reference grid, far past what an int can hold. (§I7.8)
- **Why is the test's random seed fixed at 3500?** The same seed gives the same numbers every run, so a failure can be repeated. (§I7.7)
- **Why does the halo code keep its requests in a std::vector?** The number of messages depends on whether the rank has neighbours, and MPI_Waitall takes the count and pointer the vector gives. (§I7.1)
