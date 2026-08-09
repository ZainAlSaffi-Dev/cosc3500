# Study Guide — Understanding the Heston Option Pricer From Scratch

This document explains every concept in the project in beginner's terms: the
finance, the maths, the numerics, the C++ and the HPC. PLAN.md says *what to
build*; this says *why it works*. Each section ends with self-test questions —
if you can answer those out loud without notes, you can survive the interview.

---

## Part I — The Finance (what problem are we solving?)

### 1. What is an option?

An **option** is a contract that gives you the *right, but not the obligation*,
to buy or sell something at a fixed price before a deadline.

- A **call option** = the right to *buy* at a fixed price `K` (the **strike**).
- A **put option** = the right to *sell* at strike `K`.
- The deadline is the **expiry** (for us: `T = 0.25` years = 3 months).
- "European" means you can only exercise *at* expiry, not before. (American =
  any time before expiry; we only do European — much simpler.)

Example: a call on the S&P 500 with strike 5,250. If at expiry the index is at
5,400, you exercise: buy at 5,250, worth 5,400 → payoff 150. If the index is
at 5,100, you walk away → payoff 0. In general:

```text
call payoff at expiry:  max(S − K, 0)
put  payoff at expiry:  max(K − S, 0)
```

That `max(...)` shape is the "hockey stick" you see on the expiry sheet.

### 2. What does "pricing" an option mean?

The payoff at expiry is easy. The hard question is: **what is this contract
worth today**, while we still don't know where the stock will end up? That
number is the option's *fair price* — the price at which neither buyer nor
seller has a built-in advantage. Computing it is the entire job of this
project.

### 3. What are the Greeks?

The Greeks measure how sensitive the option's value is to market changes.
We compute three:

- **Delta (Δ)** = how much the option value moves when the *stock* moves $1.
  (Slope of value vs stock price.)
- **Gamma (Γ)** = how fast delta itself changes as the stock moves.
  (Curvature of value vs stock price.)
- **Vega** = how much the value moves when *volatility* changes.
  (Slope of value vs variance, in our grid.)

Traders need these to hedge: delta tells you how many shares to hold against
the option so small stock moves cancel out. In our solver the Greeks are free:
once the grid is filled, they're just finite differences between neighbouring
cells.

### 4. Volatility, and why Black-Scholes isn't enough

**Volatility** measures how violently the stock price jiggles. **Variance**
`v` is volatility squared (we work in variance because the maths is cleaner;
volatility = `sqrt(v)`).

The classic **Black-Scholes** model assumes volatility is a *fixed constant*.
That gives a beautiful closed-form formula — but it's wrong about the real
world: markets alternate between calm and panic. Volatility *moves*.

The **Heston model** fixes this by making variance a random process of its
own, described by two coupled equations:

```text
dS = (r − q)·S·dt + sqrt(v)·S·dW1        (the stock)
dv = kappa·(theta − v)·dt + xi·sqrt(v)·dW2   (its variance)
corr(dW1, dW2) = rho
```

In words: the stock drifts and jiggles, and the *size of the jiggle* (`v`)
is itself being pulled toward a long-run level while being randomly kicked.

The five Heston parameters, one line each (memorise these):

| Parameter | Meaning | Our value |
|---|---|---|
| `v0` | today's variance (starting point) | 0.04 (= 20% vol) |
| `kappa` | speed variance is pulled back to normal | 1.5 |
| `theta` | the long-run "normal" variance level | 0.04 |
| `xi` | vol-of-vol: how hard variance gets kicked around | 0.35 |
| `rho` | correlation between stock moves and variance moves | −0.70 |

Why is `rho` negative? A famous market fact: **crashes and volatility spikes
arrive together**. When stocks fall hard, fear rises. `rho = −0.7` bakes that
in.

Other symbols: `r` = interest rate (0.045), `q` = dividend yield (0.013),
`S` = stock price, `K` = strike (5250), `T` = expiry (0.25y).

**Self-test:** What's the payoff of a put at expiry? Why does Black-Scholes
fall short? What does each of kappa/theta/xi/rho control? Why is rho negative
in practice?

---

## Part II — The Maths (how does pricing become a grid problem?)

### 5. From model to PDE

Standard financial maths (no-arbitrage arguments — you don't need to derive
this) turns the two random equations into one deterministic **partial
differential equation (PDE)** that the option value `V(S, v, t)` must satisfy
at every point:

```text
V_t + (r−q)·S·V_S + kappa·(theta−v)·V_v
    + (1/2)·v·S²·V_SS + (1/2)·xi²·v·V_vv + rho·xi·v·S·V_Sv
    = r·V
```

(Subscripts = partial derivatives: `V_S` = dV/dS, `V_SS` = d²V/dS², etc.)

You should be able to say what each term *does*:

- `V_t` — value decays as time passes (time decay).
- `(r−q)·S·V_S` — the stock's drift pushes value along the price axis.
- `kappa·(theta−v)·V_v` — mean reversion pushes value along the variance axis.
- `(1/2)·v·S²·V_SS` — diffusion (spreading) in the price direction; gamma
  lives here.
- `(1/2)·xi²·v·V_vv` — diffusion in the variance direction.
- `rho·xi·v·S·V_Sv` — the **cross-term**: price and variance moves are
  correlated, so diagonal neighbours matter. This term is why the stencil is
  9 points, not 5.
- `= r·V` — discounting: money later is worth less than money now.

Key observation for later: **every second-derivative term has a factor `v`**.
At `v = 0` all diffusion switches off. Remember this for boundaries.

### 6. The spreadsheet metaphor (how we actually solve it)

We never solve the PDE symbolically. Instead:

- Build a 2-D grid ("spreadsheet"): columns = candidate stock prices
  (0 → 4·K in `Ns = 2048` steps), rows = candidate variance levels
  (0 → 1.0 in `Nv = 512` steps).
- Each cell answers: "what would the option be worth if the market were in
  *that* state?"
- There's one sheet per moment in time, and we fill them **backwards**:
  1. **Expiry sheet is free** — it's just the payoff, `max(S−K, 0)`.
  2. **Step back** a small `dt`: each cell of the earlier sheet is a
     **weighted blend of 9 cells** of the later sheet (itself + 8
     neighbours). The weights come from the PDE terms.
  3. Repeat `nt = 2000` times until we reach today.
  4. **Read off the answer** at the cell nearest (today's spot 5200, today's
     variance 0.04). Neighbouring cells give the Greeks.

Why 9 cells? Left/right = slightly lower/higher stock price (S-derivatives).
Up/down = calmer/wilder variance (v-derivatives). **Diagonals = the
cross-term** (correlation). No cross-term → 5-point stencil like the heat
equation; Heston's correlation forces the diagonals.

### 7. Finite differences in one minute

How do you turn derivatives into cell arithmetic? Approximate them with
differences between neighbours (spacing `h`):

```text
slope:      V_S  ≈ (V[i+1] − V[i−1]) / (2h)        "central difference"
curvature:  V_SS ≈ (V[i+1] − 2·V[i] + V[i−1]) / h²
one-sided:  V_S  ≈ (V[i+1] − V[i]) / h              "forward difference"
```

Substitute these into the PDE, rearrange, and you get: *new cell value =
weighted sum of old neighbour values*. That's the entire numerical method.
The weights depend on the row's `v` and the column's `S` — which is why the
optimised solver precomputes them per row instead of re-deriving them per
cell.

### 8. Explicit scheme and the stability catch

Our scheme is **explicit**: the new sheet is computed *directly* from the old
sheet — no equation solving inside a step. Simple, fast per step, and
perfectly parallelisable later.

The price you pay: **conditional stability**. If the timestep `dt` is too
large relative to the cell spacing, errors don't die out — they *amplify*,
oscillating cell-to-cell, and the sheet explodes into a checkerboard of
garbage. Rule of thumb: `dt` must scale like the *square* of the spacing
(`dt ∝ ds²`), so doubling grid resolution means ~4× more timesteps.

Where does it blow up first? Where the diffusion coefficients are biggest —
the **high-S, high-v corner** (both coefficients grow with S and v). Our
instability animation shows exactly this.

This isn't a flaw to hide: *measuring* the stability boundary is one of our
deliverables. The alternative (implicit/ADI schemes) is unconditionally
stable but requires solving linear systems each step — out of scope, noted
as an extension.

### 9. Boundary conditions (know these cold)

The PDE needs values at the grid edges. Five conditions total:

| Edge | Story | Condition | Type |
|---|---|---|---|
| `t = T` | expiry: no uncertainty left | payoff `max(±(S−K), 0)` | terminal (exact) |
| `S = 0` | dead stock stays dead (drift & diffusion both ∝ S) | call: 0; put: `K·e^(−rτ)` | Dirichlet (exact) |
| `S = Smax` | deep in/out of the money; optionality gone | call: `S·e^(−qτ) − K·e^(−rτ)`; put: 0 | Dirichlet (approx) |
| `v = 0` | calm market, but mean reversion pushes v back up | **solve the degenerate PDE** on that row | degenerate (exact) |
| `v = vmax` | panic saturated: more vol barely adds value | `∂V/∂v = 0` | Neumann (approx) |

**The v = 0 boundary is the interview trap — be ready.** The full story:

1. Every diffusion term has a factor `v`, so at `v = 0` they all vanish.
   The PDE *degenerates* by itself into pure transport:
   `V_t + (r−q)·S·V_S + kappa·theta·V_v − r·V = 0`.
2. The drift there is `kappa·(theta − 0) = kappa·theta > 0` — pointing
   **into** the domain. Information flows inward from the boundary, so you
   don't *impose* a value — you *solve the degenerate equation* on that row.
3. Discretisation: 4-point stencil — central difference in S, **forward
   (one-sided) difference in v** (a central difference would need a cell at
   negative variance, which doesn't exist).
4. Two classic mistakes: imposing Dirichlet there (over-determines the
   problem) or using central-in-v (needs the impossible ghost point).

**The Feller condition** (bonus points): `2·kappa·theta ≥ xi²` determines
whether the variance *process* can actually reach zero. Ours:
`2(1.5)(0.04) = 0.12 < 0.1225 = 0.35²` — **violated, marginally**, which is
normal for real equity markets. The PDE treatment stays valid, but solution
detail piles up near `v = 0` which our uniform grid slightly under-resolves.
Fix would be a non-uniform grid — stated as a known limitation, not built.

**Self-test:** Why does S=0 give an *exact* condition but S=Smax only an
approximate one? What happens to the PDE at v=0 and what do we do about it?
Why forward difference, not central, at v=0? What is the Feller condition
and do our parameters satisfy it?

### 10. How we know the answer is right (validation)

We can't check against "the market" — that would test the *model*, not the
*code*. We check against mathematical ground truths:

1. **Black-Scholes collapse.** Set `xi = 0, v0 = theta` → variance never
   moves → Heston *becomes* Black-Scholes, which has an exact formula. Our
   solver must converge to it. Catches wrong weights, wrong time direction,
   most S-boundary bugs.
2. **Put-call parity.** For European options, `C − P = S·e^(−qT) − K·e^(−rT)`
   holds in *any* model, including full Heston. Exercises the cross-term and
   v-boundaries that test 1 can't reach (test 1 has `xi = 0`, killing them).
3. **Grid convergence.** Halve the spacings → the answer must settle toward
   a fixed value at the expected rate (first order in dt). Proves we
   converge to *the* solution, not a stable wrong number.
4. *(Optional but strong)* **Fourier reference.** Heston has a
   semi-analytical price via its characteristic function — a one-page Python
   integration gives the exact full-Heston price to compare against.

The iron rule: **no version is benchmarked until it matches the reference.**
A fast wrong answer is worth nothing.

---

## Part III — The Programming (what am I writing and why?)

### 11. The architecture in one picture

```text
main.cpp        parse CLI → build Config → make_solver(name) → solve → print CSV
params.h/.cpp   Config structs + config-file parser        (plain data)
grid.h/.cpp     class Grid: the two sheets + geometry      (owns all memory)
solver.h        class Solver (abstract) + SolveResult      (the interface)
solver_baseline.cpp  BaselineSolver: correct-first         (the reference)
solver_opt.cpp       OptSolver: same answers, faster       (the M1 story)
black_scholes.cpp    closed-form formula                   (validation only)
io.cpp          snapshot dumps + result CSV line           (never in timed loop)
```

Python analogies for the C++ constructs (also in code comments):

| C++ | Python equivalent | Watch out |
|---|---|---|
| `class Grid { private: ... }` | class with `_attrs` | privacy is compiler-*enforced* |
| `class Solver` + `virtual solve() = 0` | `ABC` + `@abstractmethod` | can't instantiate the base |
| `BaselineSolver : public Solver` + `override` | subclassing | `override` catches typos at compile time |
| `struct Config` | `@dataclass` | just data, all public |
| `std::vector<double>` | `list`, but typed & contiguous | THE workhorse; owns the memory |
| `const Config& cfg` param | passing an object | can't be None, can't be mutated |
| `std::unique_ptr<Solver>` | invisible GC ownership made explicit | sole owner; auto-deleted at scope exit |
| `const` method suffix | (nothing) | promise: doesn't mutate the object |

### 12. The two ideas Python never taught you

**Value vs reference semantics.** In Python, every variable is a reference
to an object on a heap somewhere. In C++, `Grid g;` **is** the object —
it lives right there (on the stack), and when the scope ends, it dies. If
you want reference behaviour you opt in with `Grid&`. Consequence: you
always know where your data physically is — which is exactly what
high-performance code needs.

**Ownership.** Python's garbage collector silently figures out when memory
can be freed. C++ makes ownership explicit: every byte belongs to exactly
one owner, and when the owner dies, the memory is freed — deterministically.
This is **RAII** (Resource Acquisition Is Initialization): the `Grid`
constructor allocates its two vectors; the destructor (run automatically at
scope exit) frees them. Our safety rules: no raw `new`/`delete` anywhere;
`std::vector` owns the sheets, `std::unique_ptr` owns the solver object;
raw pointers appear only as short-lived locals inside the hot kernel.

Why can't we hold a `Solver` by value? Because a variable of type `Solver`
must be *some concrete subclass* whose size isn't known from the interface —
same reason Java/Python interface variables are references under the hood.
C++ just makes you say it: `std::unique_ptr<Solver>`.

### 13. The double-buffer trick

Within one timestep, every new cell reads *only* the old sheet. So we keep
exactly two buffers:

- `cur` — the finished sheet, read-only during a step;
- `next` — the sheet being written, each cell exactly once;
- after the step: `std::swap(cur_, next_)` — vectors swap their internal
  pointers, O(1), no copying of 8 MB of doubles.

Nothing ever reads a half-written sheet. This discipline is what makes the
later parallel versions correct *by construction* — threads only read shared
immutable data and write disjoint cells.

### 14. Memory layout — why `idx = j*Ns + i`

The sheet is stored as ONE flat array, not a list-of-lists (Python's
`list[list]` scatters rows all over the heap — cache poison). Layout:

```text
index(i, j) = j * Ns + i      i = stock index, j = variance index
```

means an entire variance-row (`j` fixed, `i` running) is **contiguous in
memory**. The inner loop walks `i` → walks memory in order → the CPU's
caches and prefetcher love it. This single decision is the foundation of
every optimisation that follows, including SIMD in Milestone 2 (vector
instructions need adjacent data).

Cache mental model for the interview: RAM is slow (~100ns); caches are fast
(~1–10ns) but small, and data moves in 64-byte lines (8 doubles at a time).
Walk memory contiguously and 7 of every 8 loads are nearly free; jump around
randomly and every load pays full price. Cell-update arithmetic is only ~20
flops, so **this solver's speed is governed by memory access patterns, not
arithmetic** — that's why layout and loop order are the optimisation story.

### 15. The serial optimisation story (Milestone 1's 25%)

BaselineSolver is written deliberately naive: recompute all 9 stencil
weights per cell, every step. Correct, slow, honest "before". OptSolver then
applies measured steps (same answers each time, proven by `test_opt_matches`).
These are the course's own named techniques (L03) — use these names in the
video and interview:

1. **Loop-invariant code motion (hoisting) + lookup table** — the weights
   depend on the row (`v`) and column (`S`) only, not on time. Precompute
   them into a per-row table before the time loop: thousands of redundant
   per-cell multiplications become one table lookup.
2. **Strength reduction** — no division in the hot loop. Division costs
   5–14 cycles per instruction vs ~1 for multiplication, so precompute
   reciprocals (`1/ds²` etc.) once and multiply.
3. **Cache-friendly loop order** — outer `j`, inner `i`, contiguous walk
   matching the memory layout (the lecturer measured 8× from loop order
   alone on a large 2-D array).
4. **Hoisting conditionals out of loops / loop splitting** — interior loop
   runs branch-free (`i = 1..Ns−2`); edges handled in small separate loops.
   No `if` per cell.
5. **Induction-variable simplification** — replace per-cell `j*Ns + i`
   index arithmetic with a running index that just increments (the lecture
   frames this exact 2D→1D address-calculation case).
6. *(only if it measures faster)* **loop unrolling** of the inner loop.

Explicitly *not* doing: tiling/blocking (not course content), fast-math
(would compromise the correctness story), float32 (saved as an M2
experiment — half the memory traffic and twice the SIMD width).

The metric: **cell-updates per second** = `Ns·Nv·nt / seconds` — comparable
across grid sizes, which raw seconds are not.

### 16. Benchmarking properly (the other 25%)

- Numbers come **only from rangpur compute nodes** via `sbatch` (login node
  is shared and contended — numbers from it are garbage).
- ≥5 repetitions; report the **median** and the spread, never a single run.
- Record hostname, compiler version, flags, grid sizes with every result.
- Time only the solve loop, never file I/O.
- Flags fixed at `-O2`; one extra run with vectorisation disabled proves the
  serial-opt gains aren't hidden auto-SIMD (M1 forbids SIMD; `-O3` is the
  level that "enables vectorization" per the lectures — another reason we
  stay at `-O2`).
- **Profile before optimising** — with `gprof` (compile with `-pg`, run,
  then `gprof prog gmon.out`): find where time actually goes, optimise
  that. Don't polish code that's 1% of runtime.
- **Consume your results** — print the price. The lecturer demoed a
  benchmark whose runtime dropped to zero because the compiler noticed the
  result was never used and deleted the whole loop.
- **The problem must not fit in cache** — otherwise you're benchmarking
  cache, not the algorithm. Our reference grid is two 8 MB buffers; say so.
- **Verify correctness while benchmarking** — "it's very easy to make code
  run fast that doesn't work." Our bench job runs `make test` first.
- **Report honest null results** — if an optimisation doesn't help, the
  compiler may have already done it; saying so shows methodology.

**Self-test:** Why is the solver memory-bound rather than compute-bound?
Why is `j*Ns + i` the right index formula and not `i*Nv + j`? What does the
buffer swap cost? Why median of 5 runs on a compute node?

---

## Part IV — The Parallel Future (Milestone 2 preview — know the shape)

You'll be asked "where's the parallelism?" even at Milestone 1. The answer:

**Within one timestep: embarrassingly parallel.** Every cell of the new
sheet depends only on the old sheet — all `Ns·Nv ≈ 1M` cells could be
computed simultaneously in any order.

**Across timesteps: strictly sequential.** Sheet `n` can't start until
sheet `n+1` is complete. `nt = 2000` steps = 2000 mandatory
synchronisation points.

That combination — huge parallelism inside a step, a hard barrier between
steps — is what makes this a *genuine* parallel-computing problem rather
than an embarrassingly parallel one (like Monte Carlo, where paths never
talk to each other).

The M2 ladder exploits it rung by rung, same `Solver` interface each time,
using exactly what the course teaches (L04/L05) and nothing fancier:

- **SIMD = AVX intrinsics** (`<immintrin.h>`, `__m256d` = 4 doubles, FMA
  for the blend). NOT compiler auto-vectorisation — the lectures call its
  results "pretty subpar" — and not assembly. Adjacent cells in a row share
  the same formula and sit contiguously → one vector instruction computes 4
  cells. Memory aligned to 64-byte cache-line boundaries via
  `aligned_alloc`. Realistic expectation from the lectures: **~4×, not the
  full vector width** — the solver is memory-bound. A float32 variant
  (8-wide) is the natural follow-on experiment.
- **OpenMP**: `#pragma omp parallel for` on the **outer** j loop (lecture
  rule: parallelise outermost loops first), `default(none)` with explicit
  sharing lists (lecture recommendation), thread count via
  `OMP_NUM_THREADS`. Our design needs **zero synchronisation inside a
  step** — threads read the shared immutable `cur`, write disjoint rows of
  `next` — just the implicit barrier at the end of the parallel for. The
  lecture's cautionary demo: putting `critical` in a hot loop made code
  **40× slower** than serial; our no-locks-needed structure is the answer
  to "how did you avoid that?". Expect speedup ≈ physical cores (not
  hyper-threads), and accept tiny run-to-run FP differences from summation
  reordering (floating-point addition isn't associative) within test
  tolerance.
- **CUDA** (optional): one GPU thread per cell, both sheets resident on the
  GPU so nothing crosses the PCIe bus inside the time loop.

Vocabulary to have ready: **Amdahl's law** (speedup is capped by the serial
fraction — our per-step barrier is exactly such a fraction; lecture example:
90% parallel → max 10×) and **Gustafson's law** (scale the *problem* up
with the machine and efficiency returns — finer grids justify more cores).

---

## Part V — Rapid-fire interview drill

Answer these out loud. If any stumble, reread its section.

**Finance/model**
1. What's the difference between a call and a put? European vs American?
2. Why does Heston exist when Black-Scholes has a formula? (Vol moves.)
3. What does each Heston parameter do? Why is rho negative?
4. What are delta, gamma, vega, and who cares? (Hedgers.)

**Numerics**
5. Why solve backwards from expiry? (That's where the answer is known.)
6. Why 9 points in the stencil and not 5? (Cross-term/correlation.)
7. What makes the explicit scheme blow up, where first, and why there?
8. Explain the v=0 boundary to a child. Then the Feller condition.
9. How do you *know* the price is right without market data? (Three tests.)
10. What does the convergence plot prove that the other tests don't?

**Code/HPC**
11. Who owns the grid memory? When is it freed? (RAII, scope exit.)
12. Why two buffers and a swap instead of updating in place?
13. Why is the flat array + contiguous inner loop faster than list-of-lists?
14. What's the first serial optimisation and why does it win? (Hoisting
    redundant per-cell weight computation.)
15. Why benchmark on a compute node, and why median of several runs?
16. Where will the parallelism come from in M2, and what limits it?
    (Per-step barrier → Amdahl.)

**Honest-limitations questions (have these ready — they impress)**
17. Uniform grid under-resolves near v=0 (Feller violated) — would use a
    non-uniform grid given more time.
18. Explicit scheme's dt ∝ ds² restriction makes fine grids expensive —
    implicit/ADI schemes remove it at the cost of solving linear systems.
19. S_max/v_max truncation introduces (measurably small) boundary error —
    doubling the domain doesn't move the price.
