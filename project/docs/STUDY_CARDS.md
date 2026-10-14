# Study Cards. The Heston Pricer in Short

These cards are the short form of [STUDY_GUIDE.md](../STUDY_GUIDE.md). They
hold every idea needed to explain the code, and the numbers to quote with
them. Every number comes from [docs/RESULTS.md](RESULTS.md). The Heston
Pricer Study Guide page is built from this file and from
[CODE_GUIDE.md](../CODE_GUIDE.md). References like §8.2 point at the code
guide, and references like §I2.4 point at a card here.

---

## I1. The problem

### I1.1 Options and payoffs

An option gives its holder the right to trade at a fixed price on a fixed
date, with no obligation to do so. A call is the right to buy at the strike
K, and a put is the right to sell at K. Ours are European, so they can only
be used on the expiry date. The contract is on the S&P 500 with K = 5250 and
T = 0.25 years, with today's index at 5200.

```text
call payoff at expiry:  max(S - K, 0)
put  payoff at expiry:  max(K - S, 0)
```

The payoff fills the current buffer before the first step
(`Grid::init_payoff`, §5.5). The kink at the strike is the only place the
payoff is curved, so the values near the strike are the first to change.

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

Once the last step is done they cost almost nothing, because they are
differences between neighbouring cells (`extract_result`, §6.6).

### I1.3 The Heston model

Volatility is how much the price moves around, and variance v is volatility
squared. Black-Scholes treats volatility as a constant, which gives a closed
formula. Real markets move between calm and panic, so Heston makes the
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
| `xi` | vol of vol, how strongly the variance itself moves | 0.35 |
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
| ½ v S² V_SS | spreads value along S, and gamma enters here | `stock_diffusion_term` | `stock_diffusion` = 0.5v, `stock_price_sq` = S² |
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
(the payoff and the edge values) says what you are owed. So K appears in
`init_payoff` and `apply_stock_edges` and nowhere in the stencil. Pricing a
put instead of a call changes only those two.

There is no random number generator anywhere in the solver. The fair price
is an average over random futures, and the Feynman-Kac theorem says that
average satisfies the PDE. So we solve the equation the average obeys and
draw no paths.

| Route | Cost |
|---|---|
| Monte Carlo, millions of random paths | error shrinks like 1/sqrt(N), slow |
| PDE, ours | error shrinks with the grid spacing, and every (S, v) comes out at once |

Over one small step the noise moves S up or down by about sqrt(v) S.
Averaging the two cases cancels the slope, and what survives is
½ (sqrt(v) S)² V_SS, which is ½ v S² V_SS. That is why the code uses v and
never takes sqrt(v). The cross term has no ½, because the mixed partial
appears twice and the ½ cancels. An extra 0.5 there is the most common
Heston bug.

Three methods that share no code agree on the price. The PDE on the
reference grid gives 196.1684, the semi-analytic Fourier formula 196.1692,
and a Monte Carlo run of 400,000 paths 196.1345 with a standard error of
0.55 (`make mc-check`).

### I2.3 The grid of option values

The grid is a table of option values at one moment. Columns are stock prices
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

Only the values at expiry are known, so the solver starts there and steps
backwards to today, one timestep at a time. Each step is computed only from
the step before it, so only two buffers of values ever exist in memory
(§I3.1). The price is read at the * cell after the last step. Every parallel
rung splits this rectangle into rows.

### I2.4 Finite differences

A computer only holds values at cells. Finite differences turn each
derivative into arithmetic on a cell and its neighbours.

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
correlation term is the reason the stencil has nine points. Without it five
points would do. The `v = 0` row uses the one-sided `(north - centre) / dv`,
because there is no row below.

One explicit Euler step adds dt times the PDE's rate of change.

```cpp
double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
out_row[stock_i] = centre + dt * change;
```

It only reads the current buffer, so every cell can be done on its own, in
any order. Every parallel rung relies on that.

### I2.5 One step worked by hand

The test grid is 421 by 51 with 56000 steps, so ds = 50, dv = 0.02 and
dt = 4.464e-6. Take the cell at the strike, S = 5250 and v = 0.04, in the
first step. The nine values come from the payoff.

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
the curvature at the kink. The cells above and below grow with v. Today's
spot at 5200 stays at 0 until step 2, because information moves one cell per
step.

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
update is a weighted average with very uneven weights.

The centre weight shrinks as dt grows. Once it goes negative, a cell is pushed
the wrong way from its old value, small errors flip sign every step and grow,
and the grid turns into a checkerboard pattern. It goes negative first at the
high-S, high-v corner, where A and C are largest. Setting it to zero there
gives

```text
dt_stable = 1 / ( v_max S_max² / ds²  +  ξ² v_max / dv²  +  r )
```

which `report_stability` prints (§6.5). So dt must shrink like ds², and
doubling the resolution costs about four times the steps. The printed bound
leaves out the drift and cross terms, and every benchmark runs at about 0.7
of it, which held in every run (RESULTS §10). Implicit or ADI schemes remove
the limit but solve linear systems every step.

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
κθ is positive, so information flows into the grid from above, and the code
solves that equation on the row and imposes no fixed value. The v derivative
is one-sided upward, because a central one needs a row at negative variance.

`copy_top_row` reads `next`, so it must run after the interior update has
finished. In every parallel rung that is why it comes after the barrier.

The Feller condition 2κθ ≥ ξ² says whether variance can reach zero. Ours is
0.12 against 0.1225, violated by a little, which is normal for equity
calibrations. The scheme stays valid, but the region near v = 0 is slightly
under-resolved on a uniform grid. A non-uniform grid would fix that, and it
is listed as future work.

### I2.8 How we know the answer is right

Comparing with market prices would test the model. To test the code, the
tests compare with mathematical facts.

1. **Black-Scholes collapse** (`test_bs_collapse`). With xi = 0 and v0 = theta
   Heston becomes Black-Scholes, which has a formula. The relative error is
   1.16e-3 against a tolerance of 2.5e-3.
2. **Put-call parity** (`test_parity`). C - P = S e^(-qT) - K e^(-rT) holds in
   every model, so it checks the payoff, the S edges and the discounting. The
   gap is 5.4e-6 dollars against a tolerance of 1e-4. RESULTS §7 measured
   that parity cannot see the xi and rho terms at all.
3. **Fourier reference** (`make check-heston`). Heston has a semi-analytic
   price, and this check covers the xi and rho terms. On the reference grid,
   2048 by 512 by 980000 steps, the PDE gives 196.1684 against 196.1692,
   4.2e-6 relative.
4. **Ladder match** (`test_ladder_matches`). Every rung against the baseline
   within 1e-11, and omp, mpi, mpi_overlap and the alignment solvers against
   simd with zero tolerance.

Tolerances were measured first and set just above the measurement, while a
real bug moves the price by dollars. No rung is benchmarked until it passes.

Tests 1 and 2 run on grids where today's spot and v0 sit exactly on nodes,
so they could not see an error in reading the price between nodes. That
happened once. The price was read at the nearest node, which gave +4.67% on
the smoke grid. The Fourier check caught it, and `extract_result` now blends
the four surrounding cells (§6.6, §6.7). On an aligned grid the blend weight
is exactly 0, so old answers stay bit-identical.

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
| MiB, mebibyte | 1024 KiB = 1,048,576 | one buffer is 8 MiB |
| GiB, gibibyte | 1024 MiB | |
| TiB, tebibyte | 1024 GiB | keeping every time step would take 7.5 TiB |

The i in KiB and MiB marks the 1024 version. KB and MB strictly mean 1000 and
1,000,000, so 1 MiB is about 5% more than 1 MB. Chip makers often write MB for
cache sizes when they mean MiB, so the 30 MB L3 cache on an r730 socket is 30
MiB.

The sizes to know by heart.

```text
one row       2048 doubles × 8 bytes          = 16,384 bytes    = 16 KiB
one buffer    2048 × 512 doubles × 8 bytes    = 8,388,608 bytes = 8 MiB
two buffers   current and next                                  = 16 MiB
L3 cache      one r730 socket                                   = 30 MiB, so both buffers fit
```

A memory page is the fixed-size block in which the operating system hands out
memory. A program never sees real RAM addresses. Every address it uses is a
virtual address, and the operating system keeps a page table that says which
piece of real RAM each page of addresses maps to. On Linux on x86 a normal
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

The TLB holds one entry per page, so the page size decides how much memory it
covers. 1024 entries of 4 KiB cover only 4 MiB. 1024 entries of 2 MiB huge
pages cover 2 GiB.

| Page size | Pages needed for the 16 MiB of buffers | Fits in 1024 TLB entries? |
|---|---|---|
| 4 KiB, normal page | 4096 | no, so the core keeps doing page walks |
| 2 MiB, huge page | 8 | yes, with room to spare |

That is where the 1.66 times in §I3.5 comes from. The same binary on the same
grid ran at 134.8 Mcell-updates/s on huge pages and 81.1 on 4 KiB pages. The
code was the same in both runs, and only the page size of the memory changed.
Linux gives a block huge pages by itself (transparent huge pages) when the
block comes fresh from the operating system through `mmap`.
`MALLOC_MMAP_THRESHOLD_=1048576` makes every block over 1 MiB come that way,
so every benchmark sets it.

Pages also decide where memory sits on a two-socket node. Linux places each
page in the memory of the socket whose thread first writes it, one page at a
time. The official name is the first-touch policy (§I5.5).

### I3.1 Two buffers

Every new value must come from the previous time step only. Updating in place
mixes moments in time, because a cell's west and south neighbours are already
new while east and north are still old. On the test grid one step gave three
different answers for the same cell, one with two buffers, one looping left to
right and one looping right to left. With threads the answer would depend on
timing, which is a race.

So there are two buffers. During a step `current` is only read and `next` is
only written. Then the two buffers swap pointers (`std::swap` on the
vectors), so no values are copied (§5.4). Keeping every time step of the
reference solve would take 980000 buffers of 8 MiB, about 7.5 TiB.

Because of this rule the parallel rungs need no locks. Threads and ranks read
a buffer that nobody writes during the step, and each writes only its own
cells.

### I3.2 The flat array and the cache

Each buffer is one flat `std::vector<double>`, like a one-dimensional numpy
array, with the rows stored one after another.

```text
index(stock_i, var_j) = var_j * num_stock_nodes + stock_i

  row 2 |  8   9  10  11 |        memory  0 1 2 3 | 4 5 6 7 | 8 9 10 11
  row 1 |  4   5   6   7 |                 row 0  |  row 1  |  row 2
  row 0 |  0   1   2   3 |
```

A whole row is contiguous. East and west are one place apart, and north and
south a stride of one row apart, which is why
`row_above = row_centre + num_stock_nodes` (§8.5). The inner loop walks along
a row, so it walks memory in order.

RAM costs about 100 ns per access and cache about 1 to 10 ns. Data moves in
64-byte lines, eight doubles at a time, so walking in order makes seven of
every eight loads nearly free. SIMD also needs neighbouring cells side by side
in memory. Swapping the loops on purpose (the M1 `ctl-order` control) gave
the inner loop a stride of one row and ran 2.1 times slower.

### I3.3 The serial kernel is compute-bound

The cell update does 31 flops and, in the best case, moves 24 bytes (it
reads one new double of `current`, writes one of `next`, and reads the cache
line that the store goes to, which is called write-allocate). That is 1.29
flops per byte.

M1 measured the node's own limits with `slurm/bench_peak.sh`, a scalar peak
of 5.117 GFLOP/s. The baseline runs at nearly the same speed whether its data
is 1 MiB or 64 MiB, so memory does not limit it. Its five divisions per cell
do, and it reaches only 27% of that peak. After the M1 optimisations the
serial kernel ran at 82 to 95% of the peak while asking for only about
3 GB/s. So it was compute-bound, and the next gain had to come from more work
per instruction, which is SIMD. The Milestone 2 serial rung reaches
159 × 31 = 4.9 GFLOP/s, 96% of the scalar peak (RESULTS M2.5).

### I3.4 The serial optimisations

Milestone 1 tried six named techniques one at a time on rangpur, at -O2. The
M1 code is at git tag `milestone1`.

| Level | Technique | Mcell-updates/s | Kept now as |
|---|---|---|---|
| baseline | | 52.6 | `solver_baseline.cpp` |
| L1 | hoisting, lookup tables | 52.7 | `StencilConstants`, S and S² tables |
| L2 | reciprocals instead of division | 131.2 | the `inv_*` fields |
| L3 | row-major loop order | 131.5 | `var_j` outside, `stock_i` inside |
| L4 | loop splitting | 131.3 | separate boundary functions |
| L5 | induction variables | 134.2 | row pointers |
| L6 | unrolling by four | 134.8 | dropped |

Only level 2 changed the rate. GCC at -O2 already did L1 and L3 to L6,
because they leave the answer unchanged. L2 changes the answer, since
1/(2ds) is itself rounded, so the compiler is not allowed to make it without
`-ffast-math`. The trade is 2.5 times the speed for about 15 units in the
last place (3.3e-15 relative). The optimisations worth doing by hand are the
ones the compiler is not allowed to do.

The serial rung also copies every constant into a local variable (§8.4). The
loop writes through a `double*`, the compiler cannot prove that pointer never
points into the struct, and so it would reload each constant after every
store. This is pointer aliasing.

### I3.5 Page size and loop-order controls

Pages, huge pages and the TLB are explained in §I3.0. The same binary ran
1.66 times faster in one memory regime than in the other. The first 8 MiB
allocation is served by `mmap` and gets 2 MiB huge pages, so 16 MiB needs
about 8 TLB entries. Later allocations reuse the heap on 4 KiB pages, which
needs 4096 entries and page walks. `MALLOC_MMAP_THRESHOLD_=1048576` forces
the fast regime for every rep, and every benchmark sets it. A real solve is
one fresh allocation, so the fast regime is the one reported.

Two controls showed what L3 and L4 were worth. Swapping the loops gave the
inner loop a stride of one row and cost 2.1 times. The cost stays that small
because the 16 MiB of buffers fit in the 30 MiB L3. Putting the v = 0 row
back behind a per-cell `if` cost nothing, because the branch predictor gets
it right 510 times out of 511.

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
| provenance header in every CSV | node, compiler, flags, binding and regime are recorded with the number |
| print the price | a compiler can delete a loop whose result is unused |

The rate is cell updates per second, `ns × nv × nt / seconds`, so grids of
different sizes compare. The node is a 2014 Haswell Xeon E5-2670 v3, two
sockets of 12 cores at 3.1 GHz, with AVX2 and FMA and no AVX-512. A job is
capped at 15 minutes. A longer request is never rejected and waits in the
queue forever.

### I4.2 Denormals

Below about 2.2e-308 a double becomes denormal, trading precision to reach
down to 4.9e-324. On this Xeon each operation on one goes through a slow
microcode path. Our grid makes them naturally, because value creeps into the
zero region of the payoff one cell per step, so far cells hold 1e-200, then
1e-310, then 0. How many cells hold such values depends on the grid and the
step count, so different runs measured different regimes.

The fix sets two bits of the MXCSR register before `main` runs. FTZ writes
tiny results as 0 and DAZ reads tiny inputs as 0.
`slurm/flush_denormals.c` does it and is loaded with `LD_PRELOAD`, so `heston`
and the tests are unchanged and the price is the same to ten digits.
`-ffast-math` sets the same bits. Experiment E12 measured the gain per rung,
1.07 times on the baseline up to 1.42 times on simd. Without flushing simd ran
at 320 Mcell-updates/s against 454, and omp on 24 threads at 4865 against
6484 (RESULTS M2.7).

### I4.3 Profiling and -O levels

`-pg` makes the program note which function it is in every 10 ms, and
`gprof -b -p` prints a flat profile, the share of time per function. The
baseline spends 99.9% in `BaselineSolver::solve`, serial 99.4% in
`update_cells_scalar`, and simd 98.9% in `update_row_avx`. The boundaries are
under half a percent, so the interior update is the only part worth speeding
up.

| Level | baseline | serial |
|---|---|---|
| -O0 | 10.7 | 44.8 |
| -O1 | 57.6 | 156.1 |
| -O2 | 57.7 | 158.4 |
| -O3 | 112.5 | 238.6 |

Every level gave the -O2 answer exactly. -O3 adds the auto-vectoriser, here
two-lane SSE2. The project builds at -O2 so each rung measures one change of
ours, and the hand AVX2 rung (456) is 1.91 times the -O3 serial rung.
Experiment E1 later gave the compiler AVX2 as well (§I6.4).

---

## I5. The parallel rungs

### I5.1 Where the parallelism is

Inside one timestep every new cell depends only on the previous step, so all
million cells could be done at once. Across timesteps the work is strictly in
order, since step k+1 needs every cell of step k. So there is a lot of
parallelism inside a step and a hard barrier between steps. It is not
embarrassingly parallel like Monte Carlo.

```text
  SIMD, inside a row          OpenMP, rows to threads     MPI, blocks of rows to ranks
  [c c c c][c c c c]          row 7  thread 3             rank 1  rows 4..7, copy of row 3
  4 cells per instruction     row 5  thread 2             rank 0  rows 0..3, copy of row 4
                              row 1  thread 0
```

The ladder is one technique per rung, each tested against the rung below
before it is timed. The rungs are baseline, serial, simd (AVX2 and FMA), omp
(threads over rows) and mpi (ranks over blocks of rows, with OpenMP inside
each).

### I5.2 SIMD basics

SIMD is single instruction, multiple data. One instruction works on several
numbers held side by side in a wide register. The r730s have AVX2, with
256-bit `ymm` registers holding four doubles, and FMA. They do not have
AVX-512 and its `zmm` registers.

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

Four cells starting at column i need nine loads, the same four-cell load
shifted one column or one row. Lane 0 of `east` is the east neighbour of lane
0 of `centre`, so after the loads the code is the scalar formula written
with intrinsics (§9a.4).

The west, centre and east loads start at i-1, i and i+1, 8 bytes apart, so
at most one of them can be 32-byte aligned, and the aligned load
`_mm256_load_pd` crashes on an address that is not. So the shipped row uses
`loadu` for every load. Experiment E2 tested whether alignment matters, and
found a compiler problem in its place. Given only `-mavx2`, g++ 8.5 tunes for
a generic x86 and compiles each unaligned 256 bit load as two 128 bit loads
joined by `vinsertf128`, 11 of them per four cells, and each unaligned store
as two 128 bit stores. E2b added `-mno-avx256-split-unaligned-load` and
`-mno-avx256-split-unaligned-store` to the AVX2 file, and simd went from 387
to 460 Mcell-updates/s on one core with the same answers. The shipped build
uses those two flags since Phase 13 (§13.5). With unsplit loads the aligned
and unaligned instructions run at the same speed on the same address, and a
row whose centre load is 32-byte aligned runs at 477 against 412 with nothing
aligned, because an unaligned 32-byte load crosses a cache line on every
other iteration (§9f).

FMA does a multiply and an add as one instruction with one rounding. Eight of
the 24 arithmetic instructions per four cells are FMAs. It changes the last
bits, so simd matches the baseline within 3.3e-15 and the tests hold it to
the 1e-11 tolerance. With `-mfma` g++ may fuse any `a*b + c` in any file, so
only `stencil_avx.cpp` gets `-mavx2 -mfma` and the other rungs compile as
before (§13.5). The build names those two flags directly, because
`-march=native` on a newer build machine could use AVX-512 and crash on the
r730s. `objdump -d heston | grep -c zmm` printing 0 shows the binary has no
AVX-512 instruction.

The last zero to three cells of a row go through `update_cells_scalar`, the
serial code (§8.10). The Mac has no AVX2, so the file is inside
`#if defined(__AVX2__) && defined(__FMA__)` with a scalar fallback, and
`make avx-check` runs an x86_64 build under Rosetta. Rosetta and rangpur gave
the same test results.

### I5.4 OpenMP threads

A thread is a stream of instructions inside a process, and all threads of a
process share its memory, while each has its own stack for locals. Python's
GIL lets one thread run Python at a time, and C++ has no such lock. OpenMP is
a set of pragmas plus a runtime that creates and manages the threads (§9b).

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
```

| Part | Meaning |
|---|---|
| `parallel` | start a team of threads |
| `for` | split the loop's iterations across the team |
| `schedule(static)` | one contiguous block of rows per thread, decided up front |
| `default(none)` | every variable's sharing must be stated, or it will not compile |
| `shared(...)` | the same object for every thread, and `var_j` is private automatically |

Each step forks the team, each thread updates its block of rows, and an
implicit barrier at the end waits for the last one. Then the main thread does
the boundaries and the swap. The runtime keeps threads alive between loops, so
a fork costs microseconds and creates no new threads.

Rows suit it because the row loop is the outer loop (one fork per step), each
block of rows is contiguous memory, and the AVX2 inner loop is untouched. 510
rows on 4 threads split 128, 128, 127, 127.

### I5.5 Why OpenMP needs no locks, and where threads sit

A data race is two threads accessing the same memory at the same time with
at least one of them writing. Here every thread reads `current`, which nobody
writes in the step, and writes only its own rows of `next`. So there is no
`critical`, `atomic` or lock. There is also no reduction, so each cell runs
exactly the simd instructions on a different core, and omp matches simd bit
for bit at any thread count. The test checks that with zero tolerance at 7
threads, where the blocks of rows are uneven.

More threads than cores is oversubscription, and with a barrier every step
the team waits for whichever thread was paused, so it is slower.

Each r730 socket has its own L3 and memory. Linux places each page in the
memory of the socket whose thread first writes it. `Grid` fills both buffers
with zeros on the main thread, so every page is in socket 0's memory. On this
node even CPUs are socket 0 and odd are socket 1, so `close` alternates
sockets, measured with `slurm/diag_placement.sh`. Experiment E7 spread the
pages over both sockets with `numactl --interleave=all` (§I6.3). The triad
bandwidth benchmark in `bench_peak.sh` had the same placement problem and
measured 24.8 GB/s on 24 threads, and `new double[n]` in a `unique_ptr`
raised it to 67.9 GB/s (§15.12).

### I5.6 MPI ranks and copied rows

MPI runs separate processes called ranks, each with its own memory, possibly
on another machine. They share data only by messages. `mpiexec -n 4 ./heston`
starts four copies. `MPI_Init_thread`, `MPI_Comm_rank` (which am I),
`MPI_Comm_size` (how many) and `MPI_Finalize` frame every run.
`MPI_Init_thread` asks for thread level FUNNELED, which lets the master
thread call MPI from inside a parallel region, as the overlap experiment does
(§I5.8). Started without `mpiexec` the program is a group of one, which is
how the other rungs run in the same binary (§3.6).

Each rank owns a block of whole rows, and the first ranks get one extra row
when the rows do not divide evenly (§9c.2). The code calls the block
`owned_rows`, and the copied rows it receives from its neighbours `row_below`
and `row_above`. Rank 0 owns the v = 0 row and the last rank the top row. A
row is contiguous, so it goes as one message straight from the array. A
column would first need copying into a separate contiguous array.

The stencil reads one row above and one below, so the first and last rows of
a block need one row from each neighbour every step. A row of 2048 doubles is
16 KiB, so a middle rank sends and receives 64 KiB per step.

```text
  each step on every rank
    1. own interior rows with OpenMP      reads current, writes next
    2. own boundaries
    3. edge row exchange on next          post_edge_rows, then MPI_Waitall
    4. swap_buffers
```

Two ranks that each call a blocking `MPI_Send` to the other can deadlock once
the message is above the eager limit. Below it the send returns immediately,
so a test with small messages may never show the problem. So `post_edge_rows`
posts non-blocking receives first, which lets each row be received straight
into place, then non-blocking sends, and `exchange_edge_rows` calls
`MPI_Waitall`, which blocks until every send and receive has completed
(§9c.5, §9c.6).

### I5.7 Hybrid MPI and the gather

Inside each rank the OpenMP loop runs over its own rows, so the rung is MPI
between ranks and threads within each rank. One rank per socket (2 ranks of
12 threads on a node) puts each rank's pages in its own socket's memory,
because each rank writes its own copy of the grid first. The benchmark jobs
start the ranks with `--map-by ppr:1:socket:PE=12 --bind-to core`, which
gives each rank the 12 cores of one socket. `make test-mpi` uses
`--bind-to none` instead, so a test rank is never squeezed onto one core.

Every rank keeps both full buffers and updates only its rows, so every row
has the same index everywhere and every row function is reused unchanged. The
cost is 16 MiB per rank. A grid too big for one machine would give each rank
only its block of rows plus two copied rows, with the indices shifted.

After the loop `gather_rows` calls `MPI_Gatherv`, which collects every block
of rows into rank 0's next buffer. `Gatherv` takes a count and an offset per
rank because blocks of rows can differ by a row. Rank 0 then swaps the
buffers and reads the price as usual. Copied rows are exact copies and
nothing is summed across ranks, so mpi matches simd bit for bit. With the
edge row exchange removed on purpose, the zero tolerance test failed, as it
should.

### I5.8 Overlapping the exchange with the update

Experiment E9 asks whether the edge row exchange can overlap with the
interior update. `mpi_overlap` updates the two rows the neighbours need
first, posts their messages from the master thread with `post_edge_rows`,
updates every other row while the messages are in transit, and calls `MPI_Waitall`
only at the end of the step (§9g). The master thread calls MPI inside a
parallel region, which is why `main` starts MPI at thread level FUNNELED.

```text
  each step in mpi_overlap
    1. threads update the two rows to send     then #pragma omp barrier
    2. master posts the sends and receives     #pragma omp master, no barrier
    3. every thread updates its other rows     while the messages are in transit
    4. boundaries of the rows not sent
    5. MPI_Waitall, then swap_buffers
```

It helped in one layout. With one rank of 24 threads on each of two nodes the
2048 by 512 grid ran at 8009 Mcell-updates/s against 5693 (1.41 times), and
the exchange share fell from 39% to 19%. Every other layout moved by less
than 4%. Between nodes Open MPI sends a message of up to 64 KiB over TCP
under the eager limit, so the row leaves as soon as `MPI_Isend` is called.
Between two ranks on one node the eager limit is 4 KiB, so a 16 KiB row uses
the rendezvous protocol and moves only once the receiver is inside an MPI
call, which is the `MPI_Waitall` at the end of the step. This Open MPI has no
progress thread, so nothing moves while the threads compute. The limits were
not changed to test this (RESULTS M2.6).

---

## I6. Results and what they mean

### I6.1 The ladder

On r730-2, 2048 by 512, in Mcell-updates/s. The ladder is job 659238, and the
two node row is job 659242.

| Rung | Rate | Gain over the rung below | vs baseline |
|---|---|---|---|
| baseline | 57 | | 1 |
| serial | 159 | 2.80 | 2.8 |
| simd | 456 | 2.87 | 8.0 |
| omp, 24 threads | 6367 | 14.0 | 112 |
| mpi, 2 ranks of 12, one node | 6984 | 1.10 | 123 |
| mpi, 2 ranks of 12 on each of two nodes | 8310 | 1.21 over one node in the same job | 147 |

SIMD gave 2.87 times with four lanes, and four lanes set the upper limit at
4. The nine loads per four cells overlap, and the unaligned ones cross a
cache line on every other iteration. The last cells of each row run scalar.
At 454 Mcell-updates/s one core streams about 10.9 GB/s against its measured
triad bandwidth of 12.4 GB/s (RESULTS M2.4). Before Phase 13 the split loads
held simd at 387 (E2).

The strong scaling sweep gave 13.9 times on 24 threads, and the ladder job
measured 14.0. The serial boundaries grow from 0.44% of a one-thread step to
9.1% of a 24-thread step, because the main thread reads rows that other
cores have just written. The fork and barrier cost the same every step
however short the step is. At 24 threads the AVX2 kernel also moves about
151 GB/s through the caches, so shared bandwidth runs out. The scalar omp
rung, whose threads each need data at a third of the rate, scales 18.5 times.

### I6.2 Amdahl and Gustafson

Strong scaling keeps the grid fixed and adds threads. Speedup is
S(p) = T(1)/T(p) and efficiency S(p)/p. Amdahl's law caps it with the serial
fraction f.

```text
S(p) = 1 / ( f + (1 - f)/p )        rearranged   1/S - 1/p = f (1 - 1/p)
```

The second form is a straight line through the origin in f, so f is the least
squares slope `sum(x*y) / sum(x*x)`, with no fitting library. The fit gives
f = 2.42%, more than five times the boundary share, so most of the loss is
per-step overhead and bandwidth. The scalar row fits 1.57%, because a faster
row makes the fixed costs a larger share of each step.

Weak scaling keeps 21 rows per thread and grows the grid with p, stretching
`--vmax` so the row spacing stays the same. Gustafson's law gives the scaled
speedup S = f + (1 - f) p, which has no cap, where f is the serial share of
the p-thread run (here `boundary_seconds / seconds`, 9.2%). That predicts
21.9 at 24 threads against 16.3 measured. The two laws use different f, so the
strong fit is not reused. The one-thread weak point has 23 rows, 370 KiB per
buffer, which is larger than a core's 256 KiB L2. That one thread has a whole
socket's L3 and memory bandwidth to itself, so it looks a little better than
the rest. The course does not teach
Gustafson, so the slide cites Gustafson 1988.

### I6.3 Size, memory and nodes

512 by 128 on 24 threads is only 8 times serial, because each step is so
short that the fork and barrier take a large share of it. 2048 by 512 fits in
the two L3 caches and gets 40 times. 4096 by 1024 does not fit, sits on the
measured DRAM roof and drops to 19 times.

The roofline plots GFLOP/s against flops per byte (1.29 here). The flat roofs
are the measured scalar peak and the AVX2 FMA peak (43.4 GFLOP/s on one core,
971 on 24), and the sloped roofs are the measured triad bandwidth. Every rung
sits on the same vertical line and rises as it gets faster. A point can sit
above the DRAM roof when the grid fits in cache.

Page placement mattered on the big grid. `numactl --interleave=all` gave 1.74
times in the latest run and 1.28 times in the two before it, because the runs
without it run at either about 2800 or about 4000, and the median fell on a
different side. One MPI rank per socket, each writing its own copy of the
grid first so its pages sit in its own socket's memory, gave 2.28 times one
24-thread rank. Two nodes gave 1.21 times one node on the small grid and 1.64
on the big one, because a small grid step is only 0.126 ms and the edge row
exchange takes about 50 µs of it in latency and waiting.

Three things outside the code moved the numbers more than some
optimisations. Huge pages gave 1.6 to 1.8 times (E11), flushing denormals
1.07 to 1.42 times (E12), and page placement up to 2.28 times. That is why
the protocol pins all three.

### I6.4 Tried and rejected

| Experiment | Result | Decision |
|---|---|---|
| dynamic and guided schedules (E5) | 22 to 43% slower | rejected, every row costs the same |
| one parallel region per solve (E4) | 2% faster, 13% on small grids | not in the ladder, one technique per rung |
| boundaries in the parallel loop (E8) | 6% faster | not in the ladder, same reason |
| float instead of double (E3) | 1.92 times on 1 thread, 4.82 times at 4096 by 1024 | rejected for accuracy |
| aligned loads and register shuffles (E2, E2b) | g++ 8.5 split every unaligned load in two, and turning the split off gave simd 1.19 times and omp 1.12 times with no change to the source | the two flags were adopted in Phase 13, and the aligned rows stay experiments |
| compiler flags (E1) | `-O3 -mavx2 -mfma` takes the scalar row to 0.77x of simd (0.81x with unsplit loads), and adding `-ffast-math` takes the naive baseline to 380 (0.85x of simd) and serial to 1.04x of simd | rejected, fast-math drops the IEEE rules for the whole program |
| overlapping the edge row exchange (E9) | 1.41x with one rank on each of two nodes, under 4% in every other layout | kept as an experiment, because the best layouts do not gain |

Float is 4 bytes and about 7 digits, against 8 bytes and 16. A register holds
eight floats, and the float buffers of the 4096 by 1024 grid (32 MB) fit in
the two L3 caches while the double buffers (64 MB) do not, hence 4.82 times.
But every cell is rounded every step, and on the 980000 step reference grid
the price ended $0.044 from the Fourier price (double is $0.0008) and gamma
moved 11.9%. Gamma suffers most because it subtracts nearly equal numbers and
divides by ds². The double error on the far cheaper 841 by 65 grid is $0.026,
so float loses the accuracy the fine grid was chosen for.

### I6.5 Limitations and extensions

- The uniform grid under-resolves the region near v = 0, where Feller is
  violated. A non-uniform v grid would fix it.
- The explicit scheme needs dt ∝ ds², so fine grids are expensive. Implicit
  or ADI schemes remove that limit at the cost of linear solves.
- Cutting the grid off at S_max and v_max = 0.64 adds a truncation error, and
  RESULTS §7 counts it in the floor the convergence study reached.
- Every MPI rank holds the whole grid. Past one machine's memory each rank
  would hold only its block of rows and its copied rows.
- Overlapping the edge row exchange with the update (E9) helped only with one
  rank per node. Between two ranks on one node a row is above the 4 KiB eager
  limit, so it moves only once the receiver calls MPI, and the 2 x 12 layouts
  cannot overlap. Raising the eager limit or turning on a progress thread
  would test that.
- CUDA was a bonus and was not built. It would use one GPU thread per cell
  with both buffers kept on the GPU.

---

## I7. The C++ in the code

Each card lists the standard library and C++ features of one kind. Every row
has a real line from the code, what it does and why it is used there. `std::`
means the name comes from the C++ standard library, which ships with every
compiler the way Python's built-ins and standard modules do. The `#include`
at the top of a file says which part of the library the file uses.

### I7.1 std::vector, the container for every buffer and table

`#include <vector>`. A `std::vector<double>` is a resizable array that owns
its memory and keeps every element side by side, like a Python list that
holds only doubles and frees itself.

| Line | What it does | Why here |
|---|---|---|
| `std::vector<double> current_;` (grid.h) | declares an empty vector member | the buffer must be one contiguous block for the cache and SIMD (§I3.2) |
| `current_.assign(total_cells, 0.0);` (grid.cpp) | sets the size to `total_cells` and fills every element with 0 | the buffer size is only known once the config is read, and zeros give a defined start |
| `stock_price.resize(num_stock_nodes);` (stencil.cpp) | grows the vector to that many elements | the S table is then filled in a loop, one entry per column |
| `update_cells_scalar(constants, current.data(), next.data(), var_j, first_stock, end_stock);` (stencil.cpp) | `.data()` gives a raw pointer to the first element | the vector version passes the address of cell 0 to the pointer version, which holds the maths. The vector still owns the memory |
| `current[grid.index(stock_i, var_j)]` (solver_baseline.cpp) | `[]` reads one element with no bounds check | `.at()` would check every access and cost time in the hot loop. The tests catch index mistakes instead |
| `std::vector<MPI_Request> requests;` then `requests.push_back(receive_below);` (solver_mpi.cpp) | starts empty and appends one item, like `list.append` | how many messages there are depends on whether the rank has neighbours above and below |
| `MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);` | `.size()` is the count, `.data()` the pointer | MPI is a C library, so it takes a count and a pointer to the first element |
| `current[cell] = (float)grid.current()[cell];` (solver_float.cpp) | copies one cell, rounding the double to a float | a plain loop fills the float copy of the payoff |
| `stock_price[stock_i] = (float)constants.stock_price[stock_i];` (stencil.cpp) | copies one table entry, rounded to float | the float constants start from the double ones, rounded once |
| `std::swap(current_, next_);` (grid.cpp) | `#include <utility>`. Swaps the two vectors' internal pointers | ends each step with no values copied (§5.4) |

### I7.2 Ownership, references and pointers

| Line | What it does | Why here |
|---|---|---|
| `std::unique_ptr<Solver> solver = make_solver(cfg.solver);` (main.cpp) | `#include <memory>`. A pointer that owns its object and deletes it when it goes out of scope | `main` holds a solver chosen at run time without knowing its class, and nothing has to remember to free it |
| `return std::make_unique<SerialSolver>();` (solver.cpp) | creates the object and its owner in one step | it avoids writing `new`, which the solver code never does |
| `virtual ~Solver() = default;` (solver.h) | a virtual destructor with the compiler's default body | deleting a `SerialSolver` through a `Solver` pointer must run the subclass's destructor too |
| `std::unique_ptr<double[]> a(new double[n]), b(new double[n]), c(new double[n]);` (bench_peak.sh) | owns a plain array. `new double[n]` does not write the memory | a vector would write zeros to every page on the main thread and put every page in socket 0's memory. Here the parallel loop writes each page first (§I5.5) |
| `const Config& cfg` (every `solve`) | a reference, another name for the caller's object. `const` means it cannot be changed | passing by reference copies nothing, and `const` says the solver only reads the settings |
| `const std::vector<double>& current = grid.current();` | a read-only name for the grid's current buffer | the compiler refuses any write to `current` during a step, which is the two-buffer rule (§I3.1) |
| `double* out_row = next + var_j * num_stock_nodes;` (stencil.cpp) | a pointer that looks into memory it does not own | pointer arithmetic gives the start of a row once, then the loop indexes along it |
| `double* start = nullptr;` (solver_align.cpp) | a pointer that points at nothing yet, like `None` | `make_shifted_buffer` sets it to cell 0 of a buffer placed part way into a vector (§9f.1) |

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
| `FILE* out = std::fopen(path.c_str(), "w");` (io.cpp) | the C way to open a file for writing. It gives `NULL` if it fails | snapshots for the value map figures |
| `std::fprintf(out, "%.10g", ...);` then `std::fclose(out);` | writes to the file with 10 significant digits, then closes it | `%g` with 6 digits, the default, would lose detail in the figures |
| `std::filesystem::create_directories(dump_dir);` | `#include <filesystem>`, C++17. Makes the folder and any parents, and does nothing if it exists | like `os.makedirs(exist_ok=True)`. GCC 8 needs `-lstdc++fs`, which the Makefile adds only when the compiler wants it |
| `std::snprintf(file_name, sizeof file_name, "snap_%06d.csv", step);` | `#include <cstdio>`. Formats into a fixed `char` array, never writing past `sizeof file_name` | gives zero-padded names like `snap_000100.csv` that sort in time order |
| `std::printf("%.10g,%.8g,%.8g,...\n", result.price, ...);` (io.cpp) | prints a formatted line to stdout | one exact CSV line per solve that the scripts parse. `%.10g` fixes the digits |
| `std::printf("%.*g,%.*g,...", digits, result.price, ...);` (io.cpp) | `%.*g` reads the number of significant digits from the argument before the value | `--digits 17` prints enough digits to tell any two doubles apart, which experiment E1 needs (§10.1) |
| `std::fprintf(stderr, "[%s] feller %s ...", ...);` (solver.cpp) | prints to stderr | stdout stays one clean CSV line, while messages go to the other stream |

The code uses `printf` for the result line because a format string pins the
digits of every column in one place, and the CSV line reads as one line of
code.

### I7.5 Errors

| Line | What it does | Why here |
|---|---|---|
| `throw std::runtime_error(flag + " needs a value");` (main.cpp) | `#include <stdexcept>`. Stops the function and sends an error object up the call chain, like `raise` | a bad flag or config line must stop the run |
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
| `payoff = std::max(price - opt.strike, 0.0);` (grid.cpp) | the larger of two values | the payoff max(S - K, 0). Both arguments must be the same type, so the code writes `0.0` |
| `int first_interior = std::max(owned_rows.begin, 1);` (solver_mpi.cpp) | the larger of two ints | rank 0's block of rows starts at row 0, which is not an interior row |
| `int end_interior = std::min(owned_rows.end, num_variance_nodes - 1);` | the smaller of two ints | the last rank's block of rows includes the top row, which is not interior either |
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
| `enum class AlignKernel { loadu, load, shuffle };` (solver_experiments.h) | a type with three named values, like a Python `Enum` | the alignment solver stores which row function it runs |
| `explicit OmpAlignSolver(AlignKernel kernel) : kernel_(kernel) {}` | a constructor that sets `kernel_` from its argument. `explicit` stops a plain `AlignKernel` turning into a solver by accident | the three alignment solvers are one class with a different kernel each |
| `void (*row_update)(const StencilConstants&, const double*, double*, int, const double*, const double*) = update_row_avx;` (solver_align.cpp) | a function pointer, a variable that holds the address of a function | the solver picks its row function once, before the time loop (§9f.2) |
| `std::uintptr_t column_one = (std::uintptr_t)(cell_zero + 1);` (solver_align.cpp) | `#include <cstdint>`. Turns an address into an unsigned integer of the same width | `% 64` then gives the byte position inside a cache line (§9f.1) |
| `int early_rows[2] = {-1, -1};` (solver_mpi.cpp) | a fixed-size C array of two ints, both set to -1 | a thread in `mpi_overlap` updates at most two rows before the barrier (§9g.3) |
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
| `omp_get_thread_num()`, `omp_get_num_threads()` | `<omp.h>` | which thread this is and how many are in the team, so each thread in `mpi_overlap` works out its own rows (§9g.3) |
| `#pragma omp barrier`, `#pragma omp master` | OpenMP | wait for every thread, and run a block on the master thread only with no barrier at its end (§9g) |
| `MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level)`, `MPI_Finalize()` | `<mpi.h>` | start and stop MPI once in `main`. FUNNELED lets the master thread call MPI inside a parallel region (§3.6) |
| `MPI_Query_thread(&thread_level)` | `<mpi.h>` | reads back the thread level MPI gave, and `mpi_overlap` refuses to run below FUNNELED (§9g.1) |
| `MPI_Comm_rank`, `MPI_Comm_size` | `<mpi.h>` | which rank this is and how many there are, which decides the block of rows |
| `MPI_Irecv`, `MPI_Isend`, `MPI_Waitall` | `<mpi.h>` | the edge row exchange without deadlock (§I5.6) |
| `MPI_Gatherv` | `<mpi.h>` | collects blocks of rows of different sizes into rank 0 (§I5.7) |
| `_mm256_loadu_pd` and the other `_mm256_` calls | `<immintrin.h>` | each becomes one AVX2 instruction on four doubles (§I5.2) |
| `_mm256_load_pd`, `_mm256_permute4x64_pd`, `_mm256_blend_pd` | `<immintrin.h>` | the aligned load, a lane rotation and a lane blend, used only by the E2 rows (§9f) |
| `_mm_setcsr(_mm_getcsr() \| 0x8040)` | `<xmmintrin.h>` | sets the FTZ and DAZ bits in `flush_denormals.c` (§I4.2) |

---

## Drill

- **What is the payoff of a put at expiry?** The put pays max(K - S, 0) at expiry. (§I1.1)
- **Why does Heston exist when Black-Scholes has a formula?** Volatility moves in real markets, and Heston makes the variance random with its own mean reversion. (§I1.3)
- **What does each Heston parameter control, and why is rho negative?** v0 is today's variance, kappa the speed back to normal, theta the normal level, xi how strongly the variance moves and rho the correlation. Crashes and volatility spikes come together, so rho is negative. (§I1.3)
- **Where is the strike in the PDE?** The strike does not appear in the PDE. It is in the payoff and the S edge values, so only those change between a call and a put. (§I2.2)
- **Where did dW1 go, and why is the coefficient v instead of sqrt(v)?** The PDE is the equation the average over all paths satisfies (Feynman-Kac). Averaging squares the size of the noise, so sqrt(v) S becomes v S² in the diffusion term. (§I2.2)
- **Why does the cross term have no ½?** The mixed partial appears twice in Itô's lemma, so the ½ cancels. (§I2.2)
- **Why solve backwards from expiry?** The values at expiry are the only ones known, because they are the payoff. (§I2.3)
- **Why does the stencil have nine points?** The correlation term needs the cross derivative, which reads the four corner cells. (§I2.4)
- **Work out the strike cell after one step on the test grid.** d2V_dS2 = 0.02 and dV_dS = 0.5 give a change of 11025 + 84 = 11109, and times dt = 4.464e-6 that is 0.04959. (§I2.5)
- **Is the nine-point update an average of the nine cells?** It is a weighted sum whose weights are very uneven and add up to 1 - r dt. (§I2.6)
- **What makes the explicit scheme blow up, and where does it start?** A dt big enough to make the centre weight negative makes it blow up. It starts at the high-S, high-v corner, where the diffusion coefficients are largest. (§I2.6)
- **Why a forward difference in v on the v = 0 row?** A central difference needs a row at negative variance, and the drift κθ points into the grid, so information comes from above. (§I2.7)
- **What is the Feller condition, and do our parameters meet it?** The condition is 2κθ ≥ ξ². Ours is 0.12 against 0.1225, so it is violated by a little, which is normal and leaves v near 0 slightly under-resolved. (§I2.7)
- **How do you know the price is right without market data?** The tests are the Black-Scholes collapse, put-call parity, `make check-heston` against the Fourier price, which agrees to 4.2e-6 on the reference grid, and the ladder test. (§I2.8)
- **Why was the external Fourier check needed when the other tests passed?** The collapse and parity tests run where the quote sits exactly on a node, so they could not see the nearest-node readout error of 4.67%. Parity also cannot see the xi and rho terms. (§I2.8)
- **What goes wrong with one buffer?** Some neighbours already hold this step's values, so the update mixes two moments in time and depends on loop order, or on thread timing. (§I3.1)
- **What does the buffer swap cost?** The two buffers swap pointers (`std::swap` on the vectors), so no values are copied. (§I3.1)
- **Why is a flat array with a contiguous inner loop faster?** Cache lines bring in eight neighbouring doubles at once, so walking in order makes most loads nearly free, and SIMD needs neighbouring cells side by side. (§I3.2)
- **Is the solver memory-bound or compute-bound?** The baseline is limited by its divisions, and the optimised serial kernel ran at 82 to 95% of the measured scalar peak, so the serial code is compute-bound. The 24-thread rung on the 4096 by 1024 grid sits on the DRAM roof. (§I3.3, §I6.3)
- **Which serial optimisation won, and why that one?** Reciprocals in place of division won, at 2.5 times. GCC already did the others at -O2, and it is not allowed to make this one because it changes the rounding. (§I3.4)
- **Why copy the constants into locals before the loop?** The reason is pointer aliasing. The compiler cannot prove the output pointer does not point into the struct, so it would reload the constants after every store. (§I3.4)
- **Why did reps 2 to 5 run 1.66 times slower than rep 1?** Rep 1 got huge pages through mmap and later reps got 4 KiB heap pages and TLB misses. The benchmarks pin the fast regime. (§I3.5)
- **Why benchmark on an exclusive compute node, with a median of 5?** The login node is shared, another job takes cache and bandwidth even on other cores, and one rep can be unlucky. (§I4.1)
- **What is a denormal and why did it matter?** A denormal is a double below about 2.2e-308. Each operation on one is far slower, and how many there were depended on the run, so runs measured different regimes. FTZ and DAZ turn them into 0, which gave 1.07 to 1.42 times. (§I4.2)
- **Why build at -O2 when -O3 is faster?** So each rung measures one change of ours and leaves out the auto-vectoriser. The hand AVX2 rung is 1.91 times the -O3 serial rung. (§I4.3)
- **Where is the parallelism, and what limits it?** Every cell inside a step is independent, but steps are strictly in order, so there is a barrier every step and Amdahl's law applies. (§I5.1)
- **What is a __m256d?** It is a 256-bit register holding four doubles. One AVX2 instruction works on all four. (§I5.2)
- **Why does the shipped row use loadu for every load?** The west, centre and east loads start 8 bytes apart, so at most one can be 32-byte aligned, and the aligned load crashes on an address that is not. (§I5.3)
- **What did experiment E2 find?** Given only -mavx2, g++ 8.5 compiled each unaligned 256 bit load as two 128 bit loads. Two flags turn that off, which took simd from 387 to 460 with the same answers, and the shipped build now uses them. (§I5.3)
- **Why do only some files get -mfma?** With -mfma g++ may fuse multiply-adds in any file, which would change the baseline and serial answers. (§I5.3)
- **Why does the build name -mavx2 -mfma directly?** -march=native on a newer build machine could use AVX-512, which would crash on the r730s. objdump shows no zmm register in the binary. (§I5.3)
- **What does each word of the OpenMP pragma do?** parallel starts a team, for splits the loop, schedule(static) gives each thread one contiguous block of rows up front, default(none) makes every variable's sharing explicit, and shared lists the common ones. (§I5.4)
- **Why split rows between the threads?** The row loop is the outer loop, so there is one fork per step, each thread gets contiguous memory, and the AVX2 inner loop is kept. (§I5.4)
- **Why are no locks needed, and why is omp bit-identical to simd?** Threads only read current and each writes its own rows of next. There is no reduction, so each cell runs the same instructions. (§I5.5)
- **Where does Linux put a page of memory on a two-socket node?** Linux places each page in the memory of the socket whose thread first writes it. Grid fills both buffers with zeros on the main thread, so every page is in socket 0's memory. (§I5.5)
- **What is the difference between MPI and OpenMP?** Threads share one process's memory. Ranks are separate processes with separate memory and share data only by messages. (§I5.6)
- **Why is one copied row enough, and why blocks of rows?** The stencil reaches one row each way, and a row is contiguous, so it goes as one message with no packing. (§I5.6)
- **Why Irecv before Isend, and why non-blocking?** Two ranks that each make a blocking send to the other can deadlock. Posting the receive first lets each row be received straight into place. (§I5.6)
- **Why does every rank hold the whole grid?** Every row keeps its index, so the row functions and extract_result are reused. It costs 16 MiB per rank. (§I5.7)
- **Why Gatherv instead of Gather?** Blocks of rows can differ by a row, and Gatherv takes a count and an offset per rank. (§I5.7)
- **Why does main start MPI with MPI_Init_thread at FUNNELED?** The overlap solver posts its messages from the master thread inside a parallel region, which MPI allows only at thread level FUNNELED or higher. (§I5.6, §I5.8)
- **How does mpi_overlap overlap the exchange with the update?** It updates the two rows the neighbours need first, posts their messages from the master thread, updates the other rows while the messages are in transit, and waits only at the end of the step. (§I5.8)
- **Why did overlapping the exchange help only one layout?** Between nodes a 16 KiB row is under TCP's 64 KiB eager limit and leaves as soon as it is posted. Between two ranks on one node the eager limit is 4 KiB, so the row uses the rendezvous protocol and moves only when the receiver calls MPI_Waitall at the end of the step. (§I5.8)
- **Why did SIMD give 2.87 times with four lanes?** The nine loads per four cells overlap and the unaligned ones cross cache lines, the last cells of a row are scalar, and one core already streams about 10.9 GB/s of the 12.4 GB/s it can reach. (§I6.1)
- **Why only 13.9 times on 24 threads?** The serial boundaries grow to 9.1% of a step, the fork and barrier cost the same every step, and the AVX2 kernel runs out of shared bandwidth. (§I6.1)
- **How is Amdahl's f fitted with no library?** 1/S - 1/p = f (1 - 1/p) is a line through the origin, so f is sum(x*y) / sum(x*x). (§I6.2)
- **Why are Amdahl's and Gustafson's f different numbers?** Amdahl's is the serial share of the one-thread run, Gustafson's the serial share of the many-thread run. (§I6.2)
- **Why does one rank per socket beat one rank of 24 threads?** Each rank writes its own copy of the grid first, so its pages sit in its own socket's memory, while the single process has every page in socket 0's memory. (§I6.3)
- **Why do two nodes not give twice the speed?** The edge row exchange takes about 50 µs of a 0.126 ms step on the small grid. (§I6.3)
- **Why was float rejected when it was faster?** Rounding every cell every step put the reference price $0.044 from the Fourier price and moved gamma by 11.9%, worse than a far cheaper double grid. (§I6.4)
- **Could the compiler have written the SIMD rung for you?** With `-O3 -mavx2 -mfma` the scalar row reaches 0.77x of it, or 0.81x with unsplit loads. With `-ffast-math` the compiler also replaces the baseline's divisions with reciprocals, and the baseline reaches 0.85x and the serial row 1.04x of the hand-written row, but only by dropping the IEEE rules for the whole program. (§I6.4)
- **What would you do with more time?** I would add a non-uniform v grid, an implicit or ADI scheme, keeping only each rank's own rows in memory, testing the eager limit behind E9, and CUDA. (§I6.5)
- **Why does the code write std::max(price - opt.strike, 0.0) with 0.0?** std::max needs both arguments to be the same type, and price - strike is a double. (§I7.6)
- **What do the three pairs of ifs in locate do?** They keep the position on the grid, keep the cell one node in from each edge so the Greek stencils fit, and keep the blend weight between 0 and 1. (§I7.6)
- **How is a step timed?** The code reads `omp_get_wtime()` before and after, and the difference is added to `seconds`. It returns seconds as a double, so there is nothing to convert. (§I7.7)
- **Why is there no delete anywhere?** The code relies on RAII. vector, unique_ptr and ifstream free or close themselves when they go out of scope. The one `fclose` is for the C `fopen` in io.cpp. (§I7.2)
- **What happens when a config value is not a number?** std::stod and std::stoi throw std::invalid_argument, and load_config catches it and throws a message that names the file and line. Text that starts with a number, such as 12x, is read as 12. (§I7.3, §I7.5)
- **Why printf for the result line?** The format string pins the digits of every column in one place. stdout carries only that line, and messages go to stderr. (§I7.4)
- **Why is the cell count cast to double before multiplying?** ns times nv times nt is about 1e12 on the reference grid, far past what an int can hold. (§I7.8)
- **Why is the test's random seed fixed at 3500?** The same seed gives the same numbers every run, so a failure can be repeated. (§I7.7)
- **Why does the exchange code keep its requests in a std::vector?** The number of messages depends on whether the rank has neighbours, and MPI_Waitall takes the count and pointer the vector gives. (§I7.1)
- **What is row_update in the alignment solver?** It is a function pointer, set once before the time loop to the row function the solver name picked. (§I7.8)
