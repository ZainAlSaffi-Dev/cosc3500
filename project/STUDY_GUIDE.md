# Study Guide. The Maths and Ideas Behind the Heston Pricer

This document explains every concept in the project in beginner's terms. It
covers the finance, the maths, the numerics and the HPC ideas. PLAN.md says
*what to build*, and this guide says *why it works*. Each section ends with
self-test questions. If you can answer those out loud without notes, you are
ready for the interview.

**Start here.** If you only have one evening, read these sections in order.

1. §6 explains the grid as one table of values per time level, with pictures.
2. §7 explains finite differences, with the solver's first step worked by hand.
3. §11 explains why there are two buffers, and shows the bug you get with one.

Everything else is background, validation and Milestone 1 history.

**This guide is about the maths and the ideas.** The code, file by file and
snippet by snippet, is in [CODE_GUIDE.md](CODE_GUIDE.md). Read this guide
first, because the code guide points back here for the why behind each
formula. To follow the logic through the code itself, CODE_GUIDE §1.5
traces one run from `main` to the printed result, with a link to the line
of every step.

**How the code snippets work.** Every snippet is copied from the real file,
mostly with its comments removed so it stays short. The file and function are
named above it. Under each one, the first sentence says what the block does
and the second says why we do it that way.

---

## Part I. The Finance (what problem are we solving?)

### 1. What is an option?

An **option** is a contract that gives you the *right, but not the
obligation*, to buy or sell something at a fixed price before a deadline.

- A **call option** is the right to *buy* at a fixed price `K`, called the
  **strike**.
- A **put option** is the right to *sell* at the strike `K`.
- The deadline is the **expiry**. For us it is `T = 0.25` years, which is
  three months.
- A European option can only be exercised *at* expiry. An American option can
  be exercised any time before expiry. We only price European options, which
  are much simpler.

Here is an example. Take a call on the S&P 500 with strike 5,250. If the index
is at 5,400 at expiry, you exercise, buying at 5,250 something worth 5,400, so
the payoff is 150. If the index is at 5,100, you walk away and the payoff is 0.
In general the payoffs are these.

```text
call payoff at expiry:  max(S − K, 0)
put  payoff at expiry:  max(K − S, 0)
```

This block gives what each option pays on the expiry date. We need it because
it fills the grid at expiry (`Grid::init_payoff`), and the grid at every
earlier time is computed from that one.

The `max(...)` gives the payoff a hockey-stick shape, flat at zero and then
rising in a straight line. You can see it in the grid at expiry.

### 2. What does pricing an option mean?

The payoff at expiry is easy to work out. The hard question is what this
contract is worth **today**, while we still don't know where the stock will
end up. That number is the option's *fair price*, the price at which neither
the buyer nor the seller has a built-in advantage. Computing it is the whole
job of this project.

### 3. What are the Greeks?

The Greeks measure how sensitive the option's value is to changes in the
market. We compute three of them.

- **Delta (Δ)** is how much the option value moves when the stock moves by $1.
  It is the slope of value against stock price.
- **Gamma (Γ)** is how fast delta itself changes as the stock moves. It is the
  curvature of value against stock price.
- **Vega** is how much the value moves when volatility changes. In our grid it
  is the slope of value against variance.

Traders need these to hedge. Delta tells you how many shares to hold against
the option so that small stock moves cancel out. In our solver the Greeks come
almost for free. Once the grid is filled, they are just finite differences
between neighbouring cells.

### 4. Volatility, and why Black-Scholes isn't enough

**Volatility** measures how much the stock price moves around. **Variance**
`v` is volatility squared. We work in variance because the maths is cleaner,
and volatility is `sqrt(v)`.

The classic **Black-Scholes** model assumes volatility is a *fixed constant*.
That gives a closed-form formula, but real markets alternate between calm and
volatile periods, so volatility *changes over time*.

The **Heston model** fixes this by making variance a random process of its
own. It is described by two coupled equations.

```text
dS = (r − q)·S·dt + sqrt(v)·S·dW1        (the stock)
dv = kappa·(theta − v)·dt + xi·sqrt(v)·dW2   (its variance)
corr(dW1, dW2) = rho
```

This block is the Heston model itself, one random equation for the stock and
one for its variance. We write it down because the PDE in §5 comes from it,
but the C++ never uses these equations directly (§5b explains why).

In words, the stock price drifts and moves randomly, and the *size of those
random moves* (`v`) is pulled toward a long-run level while it also moves
randomly.

These are the five Heston parameters. Memorise one line for each.

| Parameter | Meaning | Our value |
|---|---|---|
| `v0` | today's variance, the starting point | 0.04 (20% volatility) |
| `kappa` | how fast variance is pulled back to normal | 1.5 |
| `theta` | the long-run normal variance level | 0.04 |
| `xi` | vol-of-vol, how large the random moves in the variance are | 0.35 |
| `rho` | correlation between stock moves and variance moves | −0.70 |

Why is `rho` negative? It is a well-known market fact that **crashes and
volatility spikes arrive together**. When stocks fall hard, fear rises, and
`rho = −0.7` builds that in.

The other symbols are the interest rate `r` (0.045), the dividend yield `q`
(0.013), the stock price `S`, the strike `K` (5250) and the expiry `T` (0.25
years).

**Self-test.** What's the payoff of a put at expiry? Why does Black-Scholes
fall short? What does each of kappa, theta, xi and rho control? Why is rho
negative in practice?

---

## Part II. The Maths (how does pricing become a grid problem?)

### 5. From model to PDE

Standard financial maths turns the two random equations into one
deterministic **partial differential equation (PDE)** that the option value
`V(S, v, t)` must satisfy at every point. The argument is a no-arbitrage one,
and you don't need to derive it.

```text
V_t + (r−q)·S·V_S + kappa·(theta−v)·V_v
    + (1/2)·v·S²·V_SS + (1/2)·xi²·v·V_vv + rho·xi·v·S·V_Sv
    = r·V
```

The subscripts are partial derivatives, so `V_S` is dV/dS, `V_SS` is d²V/dS²,
and so on.

This block is the equation the code actually solves. We solve this instead of
the random equations because it has no randomness left in it, so a grid and
plain arithmetic can handle it.

Each term becomes one named variable in `update_cells_scalar` in `stencil.cpp`,
which `update_row_scalar` calls for the whole row.
Keep this table next to you when reading the code.

| PDE term | Variable in the code | Coefficient built once in `StencilConstants` |
|---|---|---|
| `(1/2)·v·S²·V_SS` | `stock_diffusion_term` | `stock_diffusion` = 0.5v, `stock_price_sq` = S² |
| `rho·xi·v·S·V_Sv` | `cross_term` | `cross_diffusion` = ρξv |
| `(1/2)·xi²·v·V_vv` | `variance_diffusion_term` | `variance_diffusion` = 0.5ξ²v |
| `(r−q)·S·V_S` | `stock_drift_term` | `stock_drift` = r − q |
| `kappa·(theta−v)·V_v` | `variance_drift_term` | `variance_drift` = κ(θ − v) |
| `r·V` | `rate * centre` | `rate` |
| `V_t` | `centre + dt * change` | `dt` |

You should be able to say what each term *does*.

- `V_t` is time decay, the value changing as time passes.
- `(r−q)·S·V_S` is the stock's drift, which pushes value along the price axis.
- `kappa·(theta−v)·V_v` is mean reversion, which pushes value along the
  variance axis.
- `(1/2)·v·S²·V_SS` is diffusion, or spreading, in the price direction. Gamma
  lives here.
- `(1/2)·xi²·v·V_vv` is diffusion in the variance direction.
- `rho·xi·v·S·V_Sv` is the **cross term**. Price moves and variance moves are
  correlated, so the diagonal neighbours matter. This term is why the stencil
  has 9 points and not 5.
- `= r·V` is discounting, because money later is worth less than money now.

Keep one observation for later. **Every second-derivative term has a factor
`v`**, so at `v = 0` all diffusion switches off. This matters for the
boundaries.

### 5a. Where is the strike? (the model and the contract are separate)

Look at the two random equations, or at the PDE above, and you will see the
strike `K` is nowhere in them. If you are used to the Black-Scholes formula,
where `S` and `K` sit side by side, this is confusing. The answer is that
pricing always has two separate ingredients.

1. **The model** describes how the share price behaves. The two random
   equations, and the PDE they become, describe the model. They apply to *any*
   contract on this share and know nothing about yours.
2. **The contract** describes what you are owed at expiry. That is the payoff
   `max(S−K, 0)` for a call. It is the *terminal condition*, the values the
   grid is filled with at `t = T` before the solver steps backwards (the
   first row of the table in §9).

So the strike is in the program. It is used in `init_payoff` and in the
boundary values, and it does not appear in the update stencil. The
Black-Scholes *formula* feels different only because it is the special case
where both ingredients can be merged into one closed expression. Heston does
not build on that formula. It relaxes one of its assumptions (constant
volatility) and contains Black-Scholes as the `xi = 0` special case, and the
collapse test uses that special case.

Quiz yourself. If the same model prices a put instead of a call, what changes
in the program? Only the payoff and the `S` boundary values change. The PDE,
the stencil and every optimisation stay the same.

### 5b. Where did the randomness go? (the question every marker will ask)

Search the C++ for a random number generator and you will find **nothing**.
There is no `<random>`, no `rand()`, no `mt19937` and not one seed. The model
is defined by two *random* equations, and the code that solves it is
completely deterministic, because the PDE route solves directly for the
average over all random paths. You need to be able to say why in one
sentence.

**The one-sentence answer.** We never simulate the randomness because we
solved for its *average* in advance, since the PDE is the equation that the
average over all random paths satisfies.

**The theorem behind it is Feynman-Kac.** The fair price is, by definition,
an average over random futures.

```text
V(S, v, t) = E^Q[ e^(−r(T−t)) · payoff(S_T)  |  S_t = S, v_t = v ]
```

This block is the definition of the fair price, the average discounted payoff
over every random future that starts from today's state. It matters because
Feynman-Kac turns this average into the PDE, and that is the reason we can
solve an equation instead of simulating paths.

Feynman-Kac says that this average, viewed as a function of today's state
`(S, v, t)`, satisfies exactly the PDE in §5. So there are two valid ways to
get the number.

| Route | What it does | Cost |
|---|---|---|
| Monte Carlo | draws millions of random paths and averages the payoffs | the error shrinks like 1/√N, which is slow |
| **PDE (ours)** | solves the equation the average obeys | the error shrinks with the grid spacing, which is fast, and one solve gives the value at every `(S, v)` |

**The intuition behind the squaring, which is the version the video uses.**
Over one small step the noise shocks `S` up or down by an amount of size
`sqrt(v)·S`, with both directions equally likely. Average the value over the
two cases. The slope contributions cancel, because what the up shock gains
the down shock loses. The first term of the Taylor expansion that survives is
`½·(shock size)²·V_SS`, and squaring the shock size `sqrt(v)·S` gives exactly
`½·v·S²·V_SS`. This is the plain-English version of the quadratic variation
entries in the table below. `dW` itself never appears in the formula. It gets
averaged out, and the averaging leaves the curvature term behind.

**Where each piece of noise ended up.** This is the part worth memorising,
because you can point at the code and name the term.

| In the SDE | In the PDE (`update_cells_scalar` in `stencil.cpp`) | Why |
|---|---|---|
| `sqrt(v)·S·dW1` | `0.5 · v · S² · V_SS` | quadratic variation of S, where the diffusion coefficient gets **squared**, so `sqrt(v)² = v` |
| `xi·sqrt(v)·dW2` | `0.5 · xi² · v · V_vv` | the same, for variance |
| `corr(dW1, dW2) = rho` | `rho · xi · v · S · V_Sv` | cross variation `d⟨S,v⟩ = rho·(sqrt(v)S)·(xi·sqrt(v)) dt` |
| drift of S under Q | `(r − q) · S · V_S` | risk-neutral pricing, so the real-world drift `mu` does not matter |
| drift of v | `kappa · (theta − v) · V_v` | mean reversion |

Notice that **`sqrt(v)` appears nowhere in the numerics, only `v` does.**
That squaring is the sign that the noise has been averaged out. If you ever
see `sqrt(v)` in a diffusion coefficient in a PDE, something is wrong.

Two details are worth knowing before someone asks.

- **The cross term has no ½.** Itô's lemma contributes
  `½·(V_Sv + V_vS)·d⟨S,v⟩ = V_Sv·d⟨S,v⟩`. The mixed partial appears twice, so
  the ½ cancels. An extra `0.5` there is the most common Heston
  implementation bug, and ours does not have it.
- **`kappa` and `theta` are already risk-neutral.** The variance drift has no
  market price of volatility risk term `lambda`. That is standard and correct
  *for calibrated inputs*, but it is a hidden assumption. If you fed in
  `kappa` and `theta` estimated from historical data, you would be mixing the
  real-world and risk-neutral measures.

**How we checked all this.** We ran `scripts/monte_carlo_check.py` on
2026-08-09 on `config/demo.cfg`.

```text
PDE (deterministic grid)          196.1054
Monte Carlo (400k random paths)   196.1345  ± 0.5524   <- actually simulates dW1, dW2
semi-analytic (Fourier)           196.1692             <- no randomness, no mesh
PDE delta 0.566967   vs   MC delta 0.566829
sample corr(z1, z2) = −0.6991     (target rho = −0.7000)
```

This block is the output of `scripts/monte_carlo_check.py`. We ran it because
three methods that share no code and agree on the same number are the
evidence that the deterministic grid really computes the random average.

The Monte Carlo *does* use a random number generator. It correlates its two
random increments with the Cholesky factor `z2 = rho·z1 +
sqrt(1−rho²)·z_indep`, and its price agrees with the grid's answer within its
standard error. This agreement shows that the deterministic grid really
computes an average over random paths.

**Self-test.** Where is `dW1` in the code? Why is the coefficient `v` and not
`sqrt(v)`? Why does the cross term have no ½? What would break if you fed in
real-world `kappa` and `theta`?

### 6. The grid picture (how we actually solve it)

We never solve the PDE with algebra. Instead we fill a table of numbers once
for every time level. This section goes through it in order, using the small
test grid so the numbers are real.

#### 6.1 What the grid at one time level is

The grid at one time level is a table of option values at one moment in
time. The code stores it in a buffer. It keeps only two buffers (`current`
and `next`) and reuses them every step.

- Each column is a possible stock price. On the test grid the columns are
  50 apart, so column 0 is S = 0, column 104 is S = 5200 and column 105 is
  S = 5250. On the full grid there are 2048 columns from 0 to 4·K.
- Each row is a possible variance level. On the test grid the rows are 0.02
  apart, so row 0 is v = 0, row 1 is v = 0.02 and row 2 is v = 0.04. On the
  full grid there are 512 rows.
- Each cell holds one number. It is what the option would be worth if, at
  that moment, the stock price were that column's S and the variance were
  that row's v.

```text
  part of the grid         S = 5200    S = 5250    S = 5300
                         (col 104)   (col 105)   (col 106)
  v = 0.06  (row 3)         ...         ...         ...
  v = 0.04  (row 2)         ...         ...         ...
  v = 0.02  (row 1)         ...         ...         ...
```

This picture shows how to read a cell. The cell in row 2, column 104 is the
option's value in the scenario where the stock is at 5200 and the variance is
0.04, which happens to be today's market.

Here is the whole grid at one time level, with each cell marked by the
function that writes it.

```text
  var_j (variance v, rows)
    ^
nv-1|  T  T  T  T  T  T  T  T  T  T  T  T     top row, copy of the row below
    |  E  I  I  I  I  I  I  I  I  I  I  E
    |  E  I  I  I  I  I  I  I  I  I  I  E
  2 |  E  I  I  I  I  I  *  I  I  I  I  E     * = today's market, S = 5200, v = 0.04
  1 |  E  I  I  I  I  I  I  I  I  I  I  E
  0 |  E  Z  Z  Z  Z  Z  Z  Z  Z  Z  Z  E     v = 0 row
    +-------------------------------------->  stock_i (stock price S, columns)
       0                                ns-1
      S=0                              S=S_max

  I  full 9-point update            update_row_scalar
  Z  drift-only update at v = 0     update_zero_row
  E  known value at S = 0, S_max    apply_stock_edges
  T  copy of the row below          copy_top_row
```

This picture is the grid at one time level. We draw it because every rung of
the M2 ladder splits up exactly this rectangle. OpenMP hands rows to threads
and MPI hands blocks of rows to processes.

#### 6.2 One grid per time level, filled backwards

There is one grid of values per time level. We only know the values at
expiry, so we start there and work backwards towards today.

```text
   step 0            step 1             step 2                step nt
   t = T (expiry)    t = T - dt         t = T - 2dt           t = 0 (today)
  +-----------+     +-----------+      +-----------+         +-----------+
  |  payoff   | --> | computed  | -->  | computed  | --> ... | read the  |
  |  (known)  |     | from      |      | from      |         | price at *|
  |           |     | step 0    |      | step 1    |         |           |
  +-----------+     +-----------+      +-----------+         +-----------+
```

Each arrow is one timestep. The grid at each step is computed only from the
grid directly to its left. Step 1 reads only step 0, step 2 reads only step 1,
and so on. After step 1 is finished, step 0 is never read again, so only two
buffers exist in memory at any moment (§11).

#### 6.3 Building the grid at step 1 from the expiry grid

The interactive
[Stencil Stepper](https://claude.ai/artifact/AQ6V5ph6Uoyjs4hSszAR5F) runs
these same steps on the real test grid. You can click any cell to see its nine
inputs, weights and products, and animate a timestep cell by cell.

**Step 1. Fill the expiry grid (step 0) with the payoff.** At expiry there is
no uncertainty left, so each cell is just `max(S − K, 0)`. The payoff does not
depend on the variance, so every row of the expiry grid is the same.

```text
  step 0 (expiry)           S = 5200    S = 5250    S = 5300
  v = 0.06                     0           0           50
  v = 0.04                     0           0           50
  v = 0.02                     0           0           50
```

**Step 2. Pick a cell of the step 1 grid to fill.** Take row 2, column 105,
which is S = 5250 (the strike) and v = 0.04.

**Step 3. Look up the same position and its eight neighbours in the expiry
grid.** These are the nine numbers in the table above. In the first timestep
every value we read comes from the expiry grid, because that is the only
finished time level.

**Step 4. Multiply each of the nine numbers by its own weight and add them
up.** This is the nine-point update. Each of the nine cells gets its own
weight, so the result is a weighted sum, and for this cell the weights are
these.

```text
  weights for the cell at S = 5250, v = 0.04 on the test grid

                 S = 5200      S = 5250      S = 5300
  v = 0.06      +0.000057     +0.000027     −0.000057
  v = 0.04      +0.000977     +0.997976     +0.000992
  v = 0.02      −0.000057     +0.000027     +0.000057
```

So the new value of this one cell is a weighted sum of nine values from the
previous time level.

```text
  new value = 0.997976 × centre
            + 0.000992 × east       + 0.000977 × west
            + 0.000027 × north      + 0.000027 × south
            − 0.000057 × north_east + 0.000057 × north_west
            + 0.000057 × south_east − 0.000057 × south_west
```

Only the east column of the expiry grid is non-zero (50 on each of the three
rows), so only three terms survive. The two east corners cancel each other,
which leaves the east cell.

```text
  new value = 0.000991875 × 50 = 0.04959
```

The solver's own output for this cell after one step is `0.04959375`, so this
matches.

**Step 5. Do the same for every other cell of the step 1 grid.** Every cell
has its own set of nine weights, because the weights depend on that cell's S
and v. A cell with a higher stock price or higher variance gives more weight
to its neighbours.

**Step 6. Repeat for every timestep.** The grid at step 2 is built from step 1
in exactly the same way, then step 3 from step 2, and so on for all the
timesteps until we reach the grid for today.

**Step 7. Read off the answer.** The price is the value in today's grid at
S = 5200 and v = 0.04. The neighbouring cells give the Greeks.

#### 6.4 How to think about the weights

Look at the weights again and three things stand out.

- **The centre weight is almost 1.** Over one tiny timestep a cell mostly
  keeps its own old value.
- **The neighbours get very small weights.** Each step moves a cell's value a
  little towards its neighbours' values. That change is the spreading out
  (diffusion) and shifting (drift) that the PDE describes.
- **The weights add up to slightly less than 1.** They add up to exactly
  `1 − r·dt`, which is 0.9999998 here. If all nine cells held the same
  number, the new value would be that number times `1 − r·dt`, which is one
  timestep of discounting.

So the update is close to a weighted average with tiny, uneven weights,
followed by a little discounting, and some of the weights are negative. The
four corner weights come from the cross term, and their signs depend on the
correlation `rho`.

The code never stores these nine weights. It computes the five derivatives
and the five PDE terms instead, which is the same arithmetic grouped
differently (§7.4 does the same cell that way). §7.5 shows where the weights
come from, and why the solver is stable only while the centre weight stays
positive.

#### 6.5 Why nine cells?

The left and right neighbours are slightly lower and higher stock prices, and
they give the S derivatives. The cells above and below are lower and higher
variance levels, and they give the v derivatives. **The diagonal cells give
the cross term**, which comes from the correlation between the two. Without a
cross term this would be a five-point stencil like the heat equation's, and
the correlation in Heston adds the four diagonal cells.

#### 6.6 The Milestone 1 animations of the grid

These notes describe the four clips made for Milestone 1 (rendered
2026-08-09) and why they look the way they do.

- `weather_map.mp4` shows raw V(S, v) with an inferno colour ramp on a dark
  background. The first version used blue on white and was unreadable,
  because the low-value half of a light-to-dark ramp disappears into a white
  background. On a dark background with black, red and yellow, cold areas
  recede and hot areas stand out.
- `weather_map_timevalue.mp4` is the clip that shows the diffusion. Raw V
  spends the whole colour range on the straight deep in-the-money ramp (up
  to about $15,750), while the PDE changes values only within about $400 of
  the strike. Subtracting the expiry payoff (time value = V − payoff) removes
  the ramp. What remains is a region of time value centred on K = 5250. It
  widens as variance grows and narrows to a point at v = 0, where diffusion
  vanishes (every diffusion term has a factor v) and only drift remains.
  Time value is the quantity that diffuses, and this clip shows the grid
  picture of §6 directly.
- `surface.mp4` and `surface_timevalue.mp4` show the same data as a 3-D
  surface with a slow camera drift. You can see the kinked payoff ramp
  smoothing out, and the time-value surface growing out of the flat expiry
  grid, frame by frame.
- The clips use one fixed colour scale across all frames, so the animation
  does not appear to pulse. A 98th-percentile cap stops the large late values
  at high v from compressing the colour range of the early frames. NaN cells
  are magenta on dark backgrounds and black on white ones, so broken cells are
  easy to see.
- The framing was changed after feedback to suit a marker. Headlines count
  calendar time ("46 of 91 days before expiry"), with the raw solver step on
  a smaller second line. Every flat frame has a dashed line at the strike, a
  ring on today's market point (S = 5200, v = 0.04) labelled "the price we
  quote", and a one-line caption saying what the frame shows and which
  direction the solver moves. The smooth run dumps every 100 steps (200
  frames) so the motion looks continuous. One bug was fixed on the
  way. The step total must be read from the full dump list before
  `--max-frames` cuts it short, or every title shows the wrong fraction.

**The clips were measurably motionless, and three changes fixed that
(2026-08-09).** The feedback was that it was hard to see anything change.
That turned out to be literally true. These were measured on
`results/smooth` with the default render (raw V, linear scale).

- The colour bar spans $15,337 (the deep in-the-money ramp), and the entire
  change over time occupies **5.66%** of it.
- **None** of the 199 frame-to-frame transitions had a single pixel that
  moved by even one 8-bit colour level. The largest change in a frame was
  0.51 of one level, and the median frame changed by 0.05.
- 83.6% of the frame area never changed visibly (CIELAB ΔE < 2).

Three scale mismatches stacked up, and each one has a fix.

1. **The colour range was wrong.** The payoff ramp used all of it. The fix
   was `--time-value` (which already existed), plus `--gamma 0.5`, a
   power-law scale, because heat grows like `sqrt(tau)` and a linear ramp
   leaves the first half of the clip black. `--cap 99.9` replaced the old
   98th-percentile cap, which had clipped the final frame flat.
2. **The viewport was wrong.** With `s_max = 4K` and `v_max = 1.0`, the
   smoothed region around the payoff is about 5% of the S axis, and the
   reachable variance range is 12.8% of the v axis. `v_max = 1.0` sits **32.7
   standard deviations** above where the variance process actually goes. The
   fix is `--zoom`, which crops to ±3 CIR standard deviations of `v` and the
   matching lognormal spread in `S`. It uses the exact CIR variance
   formula.
3. **The time sampling was wrong.** Diffusion widths grow like `sqrt(tau)`,
   so evenly spaced frames put most of the motion at the start. 75% of an
   evenly spaced clip falls in the second half of the `sqrt(tau)`
   progression. The fix is `--warp sqrt` and fewer frames overall. Fewer
   frames was the largest single win. At 200 frames the median frame's
   visibly moving area is 0.00%, and at 128 frames it is 10.5%.

Two additions make the motion easy to see. `--envelope`
draws the ±1σ cone `S = K·exp(±sqrt(v·tau))`, the set of scenarios whose
outcome is still in doubt, and it visibly opens wider every frame.
`--price-trace` puts a live dollar readout and a growing curve of the quoted
price in the corner, so numbers move even where pixels don't.

On `results/demo` the result is that **about 10% of pixels change visibly
every frame, up from 0.00%.**

**The deeper fix was a change to the grid itself, which became
`config/demo.cfg`.** Cropping in the renderer showed that the smoke grid had
only **9 variance rows** inside the reachable range, so the zoomed frame was
blocky. The better fix is to stop spending nodes on variance levels the
market cannot reach. Holding `dv` fixed at 0.002 and `ns = 421`, `v_max = 0.4`
gives 195.8869 and `v_max = 0.2` gives 195.8876, a **0.0004%** difference.
Two-thirds of the frame, and two-thirds of the flops, were being spent on
variance levels that never matter. Combined with node alignment (spot on stock
node 130, `v0` on variance node 20), `demo.cfg` reaches **196.105 in
2.84 s**, against the P4 Richardson estimate of about 196.16. That is better
than the P4 ladder's 421×101×56000 rung, which needed 4.92 s to reach
195.925. Placing the nodes where the variance actually goes gave a better
price with fewer nodes. This is an HPC result about where to spend the
arithmetic, and it belongs in the video.

**The two Heston-becomes-Black-Scholes animations (2026-08-09).** Both answer
the same question visually, and they complement each other.

- `results/bs_collapse.mp4` comes from `scripts/bs_collapse_anim.py`. Set
  `xi = 0`, `rho = 0` and `v0 = theta`, and the Heston PDE reduces *exactly*
  to Black-Scholes. `xi = 0` removes both `0.5·xi²·v·V_vv` and the cross
  term, and `v = theta` makes `kappa(theta − v)` zero, so the `v = v0` row
  separates completely and runs the pure Black-Scholes stencil. The clip
  animates that row against the closed-form curve at the same `tau`, with a
  dollar-error panel underneath. The final-frame maximum error is
  **$0.0375**, which is 1.85 basis points of the at-the-money price. Halving
  `ds` cuts it to $0.0099, a 3.78× ratio, which is clean second order. The
  refined quoted price is **202.808075 against the closed-form 202.8079942,
  a relative error of 4.0e−07**. This is the best single validation number
  the project has. `v0` must sit exactly on a variance node or the collapse
  is only approximate, so the script stops with an error when it does not.
- `results/xi_smile.mp4` comes from `scripts/xi_sweep_smile.py`. It makes 30
  solver runs, sweeping `xi` from 0.35 down to 0, and converts each final
  grid into **Black-Scholes implied volatility per strike**. At `xi = 0.35`
  you get a skewed smile spanning **12.58 vol points**. At `xi = 0` it
  collapses to a flat line at 20.0%, with a remaining spread of **0.076 vol
  points**. That spread is the numerical noise floor of the whole pipeline,
  since the exact answer is a flat 20%. This is the picture of the sentence
  "Heston contains Black-Scholes as a special case". The other strikes come
  from the scaling property of a European call (rescaling `K` maps the grid
  onto itself because `ds` scales with `K`). This was checked against
  re-solving 5 strikes from scratch, and the worst disagreement was
  **3.9e−10**, which is the solver's print precision.

### 7. Finite differences, step by step

The PDE is written in derivatives, and a computer only has numbers stored in
cells. Finite differences turn each derivative into arithmetic on a cell and
its neighbours. This section builds that up in steps and works the solver's
real first timestep number by number. It is the same cell as §6.3, done the
way the code does it.

#### 7.1 One direction first

Take three neighbouring cells along the stock axis, spacing `h` apart.

```text
slope:      V_S  ≈ (V[i+1] − V[i−1]) / (2h)        central difference
curvature:  V_SS ≈ (V[i+1] − 2·V[i] + V[i−1]) / h²
one-sided:  V_S  ≈ (V[i+1] − V[i]) / h              forward difference
```

This block is the whole toolkit. The slope is rise over run between the two
outer cells. The curvature measures how far the middle cell sits below the
average of its two neighbours, and it is zero on a straight line. The
one-sided version is only used on the `v = 0` row, where there is no cell
below.

Here are the three expiry values around the strike on the test grid, where
`ds = 50` and the strike is 5250.

```text
  V
 50 |                          o  east    S = 5300, V = 50
    |                        /
    |                      /
    |                    /
  0 |  o---------------o          west    S = 5200, V = 0
    |  west          centre               centre  S = 5250, V = 0
    +--------------------------------> S
      5200            5250        5300
```

```text
slope      = (east − west) / (2 × 50)             = (50 − 0) / 100   = 0.5
curvature  = (east − 2 × centre + west) / 50²     = (50 − 0 + 0) / 2500 = 0.02
```

This block works both formulas on real numbers, and it shows why the
payoff's kink matters. The payoff is a straight line everywhere except at the
strike, so the curvature is zero everywhere except at the strike. The
curvature term spreads value out, so time value appears at the strike
first.

#### 7.2 The nine cells in two directions

The grid has two directions, so each cell has eight neighbours. These are the
names used in both `solver_baseline.cpp` and `stencil.cpp`.

```text
                 stock_i - 1   stock_i       stock_i + 1
               +-------------+-------------+-------------+
  var_j + 1    | north_west  |   north     | north_east  |   row_above
               +-------------+-------------+-------------+
  var_j        |   west      |   centre    |   east      |   row_centre
               +-------------+-------------+-------------+
  var_j - 1    | south_west  |   south     | south_east  |   row_below
               +-------------+-------------+-------------+
               <-- lower S                      higher S -->
```

East and west are one stock node apart. North and south are one variance row
apart. Each derivative uses a different subset of the nine cells.

| Derivative | Formula | Cells it reads | Meaning |
|---|---|---|---|
| `dV_dS` | (east − west) / (2·ds) | west, east | slope in S, which is delta |
| `d2V_dS2` | (east − 2·centre + west) / ds² | west, centre, east | curvature in S, which is gamma |
| `dV_dv` | (north − south) / (2·dv) | south, north | slope in v, which is vega |
| `d2V_dv2` | (north − 2·centre + south) / dv² | south, centre, north | curvature in v |
| `d2V_dS_dv` | (NE − NW − SE + SW) / (4·ds·dv) | the four corners | the cross term |

The cross derivative is the only one that reads the corners. It is the S
slope on the row above minus the S slope on the row below, divided by the
distance between those rows.

```text
d2V_dS_dv = [ (NE − NW) / (2·ds)  −  (SE − SW) / (2·ds) ] / (2·dv)
          = (NE − NW − SE + SW) / (4·ds·dv)
```

This block derives the corner formula. We need it because the correlation
term `rho·xi·v·S·V_Sv` asks how the S slope changes as v changes, and only
the corners carry that information. Without correlation the stencil would be
the five-point cross shape.

#### 7.3 From derivatives to the new value

The PDE says how fast V changes as we move away from expiry. Written in time
to expiry, it is `dV/dτ = (the five terms) − r·V`. One explicit Euler step
adds `dt` times that rate to the old value.

```cpp
double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
out_row[stock_i] = centre + dt * change;
```

This block is the update for one cell, copied from `update_cells_scalar`. We do
it this way because it uses only values from the current buffer, so every cell
can be computed on its own, in any order. The OpenMP and MPI rungs rely on
that.

#### 7.4 Worked example, the solver's real first step

The test grid is 421 × 51 cells with 56000 timesteps. That gives `ds = 50`,
`dv = 0.02` and `dt = 0.25 / 56000 = 4.464e-6`. We update the cell at
`stock_i = 105`, `var_j = 2`, which is S = 5250 (the strike) and v = 0.04.
The current buffer holds the expiry grid, filled with the payoff, so these
are the nine values the stencil reads.

```text
              S = 5200    S = 5250    S = 5300
  v = 0.06        0           0          50       row_above
  v = 0.04        0         [ 0 ]        50       row_centre, centre in brackets
  v = 0.02        0           0          50       row_below
```

**Step A is the five derivatives.** The payoff does not depend on v, so every
north value equals its south value and all the v derivatives are zero.

```text
dV_dS     = (50 − 0) / (2 × 50)                = 0.5
d2V_dS2   = (50 − 2 × 0 + 0) / 50²             = 0.02
dV_dv     = (0 − 0) / (2 × 0.02)               = 0
d2V_dv2   = (0 − 2 × 0 + 0) / 0.02²            = 0
d2V_dS_dv = (50 − 0 − 50 + 0) / (4 × 50 × 0.02) = 0
```

**Step B is the named terms.**

```text
stock_diffusion_term    = 0.5 × 0.04 × 5250² × 0.02    = 11025
cross_term              = ρξ × 0.04 × 5250 × 0         = 0
variance_diffusion_term = 0.5ξ² × 0.04 × 0             = 0
stock_drift_term        = (0.045 − 0.013) × 5250 × 0.5 = 84
variance_drift_term     = 1.5 × (0.04 − 0.04) × 0      = 0
rate * centre           = 0.045 × 0                    = 0
change                  = 11025 + 84                   = 11109
```

**Step C is the new value.**

```text
new value = centre + dt × change = 0 + 4.464e-6 × 11109 = 0.04959
```

The solver's own snapshot after step 1 holds `0.04959375` in that cell, so
the hand calculation and the code agree to every printed digit. It is also
the same number §6.3 got from the nine weights. The cell was worth nothing at
expiry, and one step earlier it is worth about 5 cents. That is time value,
the chance that the stock still finishes above the strike. Almost all of it
came from the curvature term, because the kink is where the payoff bends.

The same column on the neighbouring rows after step 1 shows how v enters.

```text
  v = 0.06    0.07420
  v = 0.04    0.04959
  v = 0.02    0.02498
```

The values grow in proportion to v because `stock_diffusion = 0.5 × v`. More
variance means more chance of finishing in the money. After step 1 the rows
differ, so from step 2 onwards the v derivatives are no longer zero and the
variance terms start contributing.

#### 7.5 The same step as nine weights

§6.3 described the update as nine weights times nine cells. The weights come
from the same update the code does, with the arithmetic regrouped by cell.
Each derivative formula is already a mix of cells, the PDE mixes the
derivatives together, and adding up how much of each cell ends up in the
total gives that cell's weight.

**Step 1. Write down the update.** To keep it short, call the five PDE
coefficients A = 0.5·v·S², B = (r − q)·S, C = 0.5·ξ²·v, D = κ(θ − v) and
X = ρ·ξ·v·S. The update from §7.3 is then

```text
new = centre + dt × ( A·d2V_dS2 + B·dV_dS + C·d2V_dv2 + D·dV_dv + X·d2V_dS_dv − r·centre )
```

**Step 2. Treat each derivative as a recipe of cells.** For example
`d2V_dS2 = (east − 2·centre + west) / ds²` takes 1/ds² of east, −2/ds² of
centre and 1/ds² of west. The §7.2 table gives the recipe for each of the
five.

**Step 3. Multiply each recipe by its coefficient and write down the amount
each cell receives.** These are the numbers for the worked example cell (S =
5250, v = 0.04, ds = 50, dv = 0.02). D is 0 here because v equals θ.

| Term | west | centre | east | north | south | NE and SW | NW and SE |
|---|---|---|---|---|---|---|---|
| A·d2V_dS2, which is A/ds² × (1, −2, 1) | 220.5 | −441 | 220.5 | | | | |
| B·dV_dS, which is B/(2·ds) × (−1, +1) | −1.68 | | +1.68 | | | | |
| C·d2V_dv2, which is C/dv² × (1, −2, 1) | | −12.25 | | 6.125 | 6.125 | | |
| D·dV_dv, which is D/(2·dv) × (+1, −1) | | | | 0 | 0 | | |
| X·d2V_dS_dv, which is X/(4·ds·dv) × (+1, −1) | | | | | | −12.86 | +12.86 |
| −r·centre | | −0.045 | | | | | |
| **total** | 218.82 | −453.295 | 222.18 | 6.125 | 6.125 | −12.86 | +12.86 |

This table shows where each cell's share comes from. We build it because it
makes visible that east and west get most of their weight from the curvature
term (220.5 each) and only a little from the drift term (±1.68), while the
centre loses something in every column.

**Step 4. Multiply every total by dt, then add 1 to the centre.** The extra 1
is the `centre +` at the front of the update, and dt is 4.464e-6. This gives
the weights in §6.3.

```text
east     222.18 × dt          = 0.000991875
west     218.82 × dt          = 0.000976875
centre   1 − 453.295 × dt     = 0.997976
north    6.125 × dt           = 0.000027
south    6.125 × dt           = 0.000027
NE, SW   −12.86 × dt          = −0.000057
NW, SE   +12.86 × dt          = +0.000057
```

The same steps with letters instead of numbers give the general formula for
any cell.

```text
weight on east        = dt × ( A / ds²  +  B / (2·ds) )
weight on west        = dt × ( A / ds²  −  B / (2·ds) )
weight on north       = dt × ( C / dv²  +  D / (2·dv) )
weight on south       = dt × ( C / dv²  −  D / (2·dv) )
weight on NE and SW   = dt × X / (4·ds·dv)
weight on NW and SE   = −dt × X / (4·ds·dv)
weight on centre      = 1 − dt × ( 2·A / ds²  +  2·C / dv²  +  r )
```

This block is the update rearranged so each cell has one weight. We show it
because it explains the weights in §6.3, and because the centre weight
decides whether the scheme is stable.

**Why the weights add up to `1 − r·dt`.** The drift terms (B and D) and the
cross term (X) each add the same amount to one cell as they take from
another, so they cancel. The two curvature terms put +1, −2 and +1 on three
cells, which also cancels. Only the leading 1 and the −r·dt are left.

**Why the centre weight must stay positive.** Look at the centre weight. It
starts at 1 and loses a bit for every unit of `dt`. If `dt` is too large, it
goes below zero. A negative centre weight means a cell's new value is pushed
the opposite way from its old value. A tiny error then flips sign every step
and grows, which is the checkerboard pattern in §8. The centre weight is
smallest at the far corner, where S and v are largest. Setting it to zero
there gives

```text
dt_stable = 1 / ( v_max · S_max² / ds²  +  ξ² · v_max / dv²  +  r )
```

which is exactly the formula `Solver::report_stability` prints as
`dt_stable`. This block is the solver's stability estimate. We print it at
the start of every run because a `dt` above it makes the run blow up.

#### 7.6 Information moves one cell per step

This is the row v = 0.04 around the strike, taken from the solver's
snapshots.

```text
             S = 5150    S = 5200     S = 5250    S = 5300    S = 5350
  step 0       0           0            0           50          100
  step 1       0           0            0.04959     50.00075    100.00074
  step 2       0           0.0000455    0.09909     50.00155    100.00149
```

S = 5200 (today's spot) only becomes non-zero at step 2. The stencil reads
one cell in each direction, so in one step nothing can travel further than
one cell. That is another way to picture why the explicit scheme needs many
small steps.

#### 7.7 Edge cells

The nine-point stencil needs a neighbour on every side, so the outer ring of
cells is handled separately, as in the picture in §6.1. The `S = 0` and
`S = S_max` columns get known values. The `v = 0` row has no row below, so it
uses the one-sided difference for `dV_dv`, and the diffusion terms drop out
because they all have a factor of v. The top row copies the row below, so the
value stops changing with v there. §9 explains each one.

### 8. Explicit scheme and its stability limit

Our scheme is **explicit**, which means the grid at the new time level is
computed *directly* from the old one, with no equation solving inside a step.
It is simple, fast per step, and easy to parallelise.

The cost is **conditional stability**. If the timestep `dt` is too large
compared with the cell spacing, errors don't die out. They *grow*,
oscillating from cell to cell, and the grid turns into a checkerboard of
meaningless numbers. As a rule of thumb, `dt` must scale like the *square* of
the spacing (`dt ∝ ds²`), so doubling the grid resolution needs about 4× more
timesteps. §7.5 shows where this comes from.

It blows up first where the diffusion coefficients are biggest, which is the
**high-S, high-v corner**, because both coefficients grow with S and v. Our
instability animation shows exactly this.

*Measuring* the stability boundary is one of our deliverables. The
alternative, implicit or ADI schemes, is stable for any `dt` but has to solve
linear systems every step. That is out of scope and noted as an extension.

**The measured boundary (P5, 2026-08-09, smoke grid 256×64, T = 0.25).** The
solver prints dt_stable ≈ 1.526e-5. We bracketed the real limit by sweeping
nt.

- nt = 15000 (dt = 1.092× the bound) still finishes with a finite price of
  205.3226. The 4th decimal is already drifting, which is a borderline mode
  decaying slowly.
- nt = 14900 (1.100×) ends at 2.8e+69, caught in the middle of blowing up.
- nt = 14600 (1.122×) ends in NaN.

So the real stability limit is about 9% above the printed estimate. The
estimate is cautious because it uses the single worst cell's coefficients (the
far corner), while an unstable mode spreads over a *region*, so its growth
rate averages in neighbouring cells with smaller coefficients. A bound that
errs on the safe side is the right kind of estimate, and it is worth saying so
in the video. The `unstable.cfg` config runs at 4× the bound. The checkerboard
starts at the high-S, high-v corner and overflows to NaN by about step 450 of
4096. `weather_map.py --diverging` films it, with red and blue showing the
plus and minus oscillation and black showing overflowed cells.

### 9. Boundary conditions (know these cold)

The PDE needs values at the edges of the grid. There are five conditions in
total.

| Edge | Story | Condition | Type |
|---|---|---|---|
| `t = T` | at expiry there is no uncertainty left | payoff `max(±(S−K), 0)` | terminal (exact) |
| `S = 0` | a stock worth zero stays at zero, since drift and diffusion are both ∝ S | call 0, put `K·e^(−rτ)` | Dirichlet (exact) |
| `S = Smax` | deep in or out of the money, so the option behaves like the stock or like nothing | call `S·e^(−qτ) − K·e^(−rτ)`, put 0 | Dirichlet (approximate) |
| `v = 0` | a calm market, but mean reversion pushes v back up | **solve the reduced PDE** on that row | degenerate (exact) |
| `v = vmax` | at very high variance, more volatility barely adds value | `∂V/∂v = 0` | Neumann (approximate) |

**The v = 0 boundary is the most likely interview question, so be ready.**
This is the full story.

1. Every diffusion term has a factor `v`, so at `v = 0` they all vanish. The
   PDE reduces by itself to pure transport, `V_t + (r−q)·S·V_S +
   kappa·theta·V_v − r·V = 0`.
2. The variance drift there is `kappa·(theta − 0) = kappa·theta > 0`, which
   points **into** the grid. Information flows inward from the boundary, so
   the code *solves the reduced equation* on that row, and imposing a value
   there would over-determine the problem.
3. The discretisation is a 4-point stencil. It uses a central difference in S
   and a **forward (one-sided) difference in v**, because a central
   difference would need a cell at negative variance, which doesn't exist.
4. There are two classic mistakes. One is imposing a fixed value there, which
   over-determines the problem. The other is using a central difference in v,
   which needs the impossible cell below zero.

**The Feller condition is worth bonus points.** The condition
`2·kappa·theta ≥ xi²` decides whether the variance *process* can actually
reach zero. Ours gives `2(1.5)(0.04) = 0.12 < 0.1225 = 0.35²`, so it is
**violated, but only just**, which is normal for real equity markets. The PDE
treatment stays valid, but detail in the solution piles up near `v = 0`,
which our uniform grid slightly under-resolves. The fix would be a
non-uniform grid. We state it as a known limitation and did not build it.

**Self-test.** Why does S = 0 give an *exact* condition but S = Smax only an
approximate one? What happens to the PDE at v = 0, and what do we do about
it? Why use a forward difference at v = 0? What is the
Feller condition, and do our parameters satisfy it?

### 10. How we know the answer is right (validation)

A check against market prices would test the *model* and say nothing about
the *code*, so we check against mathematical ground truths.

1. **Black-Scholes collapse.** Set `xi = 0` and `v0 = theta`. Then the
   variance never moves, Heston *becomes* Black-Scholes, and Black-Scholes
   has an exact formula. Our solver must converge to it. This catches wrong
   weights, a wrong time direction and most S-boundary bugs.
2. **Put-call parity.** For European options, `C − P = S·e^(−qT) − K·e^(−rT)`
   holds in *any* model, including full Heston. This exercises the cross term
   and the v boundaries, which test 1 can't reach because `xi = 0` switches
   them off.
3. **Grid convergence.** When the spacings are halved, the answer must settle
   toward a fixed value at the expected rate (first order in dt). This shows
   the answer converges to the solution of the PDE itself.
4. **Fourier reference (done 2026-08-09).** Heston has a semi-analytical
   price through its characteristic function. `scripts/monte_carlo_check.py`
   computes it (numpy only, using Albrecher's little Heston trap form) *and*
   runs a Monte Carlo that simulates the random equations directly. Any
   disagreement between two methods can then be settled by the third.

The validation rule is that **no version is benchmarked until it matches the
reference.** A fast wrong answer is worth nothing.

**Tests 1 to 3 share a blind spot, which is why test 4 mattered.** Tests 1
and 2 compare the solver *against itself* at the same grid cell. At the time,
`extract_result` read the price at the cell **nearest** to (spot, v0) with no
interpolation. A shared offset in the quote point therefore cancels exactly
in put-call parity, and appears identically on both sides of the
Black-Scholes comparison, so both tests pass while the absolute price is
wrong. That happened on `config/smoke.cfg`. There `dv = 1/63 = 0.015873`, so
`v0 = 0.04` snapped to node 3, which is **0.047619**. That is a 19% increase
in variance, worth about +$16 at the reported vega. The stock snap was worth
about −$6.6, so the net error was +$9.2, and smoke's quoted 205.32 was indeed
**+4.67%** above the true 196.1692. The cause was the readout, because the
same run matched the semi-analytic price *evaluated at its own snapped node*
to within $0.29, which rules out instability and discretisation error.
Node-aligned configs (the P4 ladder and `config/demo.cfg`) never had this
problem, and on `demo.cfg` the three routes agree to 0.03%. **A validation
suite made only of self-consistency checks cannot see a systematic offset,
and one comparison against an external absolute reference can.** That rule
carried into Milestone 2.

**The fix (P7) was bilinear interpolation at the quote point.**
`extract_result` no longer snaps. `Grid` gained `stock_position()` and
`variance_position()`, which give the cell just below the quote point plus
how far across that cell the point sits. The price is a blend of the four
surrounding cells, and each arm of the Greek stencils (east, west, north and
south) is blended the *same* way. That way the reported delta really is the
derivative of the reported price.

Three details are worth being able to defend on camera.

- **Node alignment still had to give bit-identical answers.** The whole P3
  and P4 record is quoted from aligned grids. That includes `demo.cfg`'s
  196.1054098, the four convergence rungs and both test grids. On an aligned
  grid the weight is exactly 0, and the blend is written so that
  `weight = 0` reduces to `1.0*V + 0.0 + 0.0 + 0.0`, which returns the node's
  value unchanged. We checked this, and every one of those numbers is still
  identical to the last digit after the change. Anything else would
  have meant the interpolation was wrong.
- **The weight is snapped to 0 within 1e-9 of a node.** `0.04 / (0.2/100)`
  happens to come out as exactly 20.0 in IEEE double on these configs, but
  relying on that would be relying on luck. A grid whose position came out
  as 19.999999999999996 would quietly interpolate between nodes 19 and 20 and
  break the bit-identical guarantee. The snap makes the guarantee part of the
  design.
- **The cell index is clamped to `[1, n-3]`.** The blend reads
  `cell` and `cell+1`, and the gamma stencil reads one further out on each
  side. So `cell-1` and `cell+2` must both exist, which costs one more node of
  margin at each edge than the old nearest-node readout needed.

The solver also **says when a config is off-node**. It prints one stderr line
naming the offset in cells, with the same `[baseline]` or `[serial]` prefix as
the Feller line. Node alignment is still the rule. The warning exists so that
breaking it is always visible.

The measured effect is that `smoke.cfg` went from 205.3242902 to
**195.247466**. The remaining gap to 196.1692 is discretisation error on a
256×64 grid (`dv = 1/63` is very coarse), and the +4.67% snap error is
gone.

**Notes on how the two test programs work (P3).**

- **The normal CDF is computed without an approximation.** Black-Scholes needs
  N(x), the standard normal CDF. `<cmath>` has no `norm_cdf`, but it has
  `erfc`, the *complementary error function*, and N(x) = ½·erfc(−x/√2) is an
  exact identity. It takes one line and has no polynomial fit to defend.
- **The test grid is node-aligned.** When the tests were written,
  `extract_result` read the *nearest* cell with no interpolation. On the smoke
  grid the nearest stock node is $11.80 away from spot 5200, and at delta
  ≈ 0.56 that is about $6.60 of price error from the snap alone, bigger than
  everything we were trying to measure. So the tests pick node counts that put
  spot and v0 *exactly* on nodes. `ns = 421` gives a spacing of 50, so 5200 is
  node 104 and the strike 5250 is node 105. `nv = 51` gives a spacing of 0.02,
  so v0 = 0.04 is node 2.
- **Parity is near-exact on the grid as well as on paper.** The parity
  portfolio (long a call, short a put) has payoff S − K, which is *linear* in
  S. Central differences are exact on linear functions, every boundary
  formula treats it exactly, and it does not depend on v. So the scheme
  carries it forward with only compounding error. After nt steps the grid
  holds (1 − r·dt)^nt where the formula says e^(−rT), a gap of about 1e-5
  dollars here. A tight dollar tolerance is therefore fair for parity, while
  the Black-Scholes test carries real discretisation error and needs a
  measured tolerance.
- **Tolerances come from measurements.** Each test was run once to record its
  actual error, and the tolerance was then set just above it. On the
  node-aligned 421×51 grid with nt = 56000 (2026-08-09), the Black-Scholes
  relative error was 1.16e-3 for the call and 1.12e-3 for the put, so the
  tolerance is 2.5e-3, about twice the measurement. The parity gap was 5.4e-6
  dollars, so the tolerance is 1e-4. That is 20× the measurement, which leaves
  room for compiler floating-point differences on rangpur and is still about 5
  orders of magnitude below the size of a bug. A real bug moves the answer by
  dollars, so these margins catch regressions and leave room for
  floating-point noise. The plan's target of a relative error below 1e-3 is
  for the *reference-sized* grid, with a spacing of about 10 against 50 here,
  and the P4 convergence study demonstrates it.

**The measured convergence (P4, redone 2026-08-10 because the earlier study
was invalid, see below).** The ladder is node-aligned, dt drops by a quarter
per rung, and the error is measured against the **semi-analytic Fourier price
196.1692090** (point 3 below says why).

| Grid | dt | Price | Error vs analytic |
|---|---|---|---|
| 421×33 | 4.464e-6 | 195.635524 | 0.533685 |
| 841×65 | 1.116e-6 | 196.143695 | 0.025514 |
| 1681×129 | 2.790e-7 | 196.155137 | 0.014072 |

The successive error ratios are 20.9 and 1.81, so the observed order in dt is
**2.19, then 0.43**. The error falls fast and then *flattens onto a floor*,
so the data has no single clean rate. Quoting the mean (1.31) would describe
a straight line that the data does not follow. RESULTS.md §7 reads the floor
as error that this refinement does not reduce, from truncating the domain at
`s_max` and `v_max` and from the Fourier anchor's own quadrature, which takes
over once the discretisation error drops beneath it. Confirming that needs
the fourth rung (3361×257, about 1.3 hours), which is in `LADDER` and has not
been run.

**Three things were wrong with the earlier version of this study, and all
three are worth understanding.**

1. **The ladder had gone out of date against its own config.** It picked `nv`
   from multiples of 25 so that v0 = 0.04 sits on a node. That was correct
   when `v_max = 1.0`, which is what `reference.cfg` said when the ladder was
   written. P7 changed `v_max` to 0.64 for stability, and 0.04/(0.64/25) =
   **1.5625**, so *every rung* was off-node in variance. The solver was
   warning about it on stderr the whole time, but the script called
   `subprocess.run(..., capture_output=True)`, which captured the warning and
   never printed it.
2. **The strike was in the middle of a cell on the coarse rungs.** The payoff
   has a kink at K = 5250. For the kink to be discretised the same way on
   every rung, ds must divide it, but at ds = 200 the strike sits at node
   26.25. A kink that moves relative to the grid adds a sawtooth on top of the
   trend, which is why an early re-run gave an error sequence that went *down
   and then back up*. ds must divide gcd(5200, 5250) = 50, so the ladder now
   starts at ns = 421.
3. **The study was measured against itself.** Defining the error as
   |price − price_finest| means the finest rung gets no error bar, the most
   expensive rung produces a reference instead of a data point, and the
   reference's own error is silently subtracted from every other rung. The
   external Fourier price fixes all three.

`scripts/convergence_plot.py` now works out the alignment requirement from the
config at run time and refuses to start if any rung misses it. It re-raises
the solver's off-node warning, and it reports the order between every pair of
neighbouring rungs. The figure is `results/convergence.png`, with the raw
numbers alongside as CSV. Regenerate it with `make convergence`.

**The reference grid, finally run (P7, 2026-08-09).** `config/reference.cfg`
had never actually produced a number. It specified `nt = 2000`, which is 527×
over the explicit stability bound, so it printed `nan,nan,nan,nan`. The P2
acceptance criterion "price is positive and finite" had been ticked without
anyone running it. With `nt = 980000` and `v_max = 0.64` it gives these
results.

| Source | Price | Delta |
|---|---|---|
| PDE, 2048×512×980000 | **196.1683699** | 0.56717178 |
| semi-analytic (Fourier inversion) | 196.1692 | n/a |
| Monte Carlo, 400 000 paths | 196.1345 ± 0.5524 | 0.566829 ± 0.005377 |

The PDE differs from the semi-analytic price by **−0.0008 dollars, a relative
error of 4.2e-6**. Three routes to the same number agree to four parts per
million. One is a deterministic PDE, one is a deterministic integral
transform, and one uses actual random paths. This is the accuracy claim to
make on camera. It is credible because the Fourier price and the Monte Carlo
price are independent of each other and of the grid. If the PDE had
disagreed, the grid would have been the first suspect.

**Do not quote a runtime for this run on camera.** It was measured on the
laptop, and there is no rangpur `sbatch` job or stamped file behind it, so it
must not appear beside cluster numbers. A single reference solve
exceeds the 15-minute QOS limit by a wide margin and has never been run
there. If a runtime is wanted, derive it from the measured rangpur throughput
and say that you are deriving it. That gives about 2.1 hours for the
Milestone 1 optimised solver and about 5.4 hours for the baseline
(RESULTS.md §8).

There are two footnotes worth knowing if asked. First, the reference grid is
**not** node-aligned (`ds = 21000/2047` puts spot at 506.9 cells), so this
number depends on the P7 bilinear interpolation. With the old nearest-node
readout it would have been off by roughly delta × half a cell, about $2.9.
That code path is now tested. `make check-heston` leg 2 runs a deliberately
off-node grid and checks the interpolated price against the analytic one.
Until 2026-08-10 every automated test used node-aligned grids, where the
blend weights are exactly zero and the interpolation code never ran. Second,
cutting the domain off at `v_max = 0.64` contributes about 2e-4 dollars,
three orders of magnitude below the agreement above. That level is 80%
annualised volatility, the highest VIX reading in modern history, and
P(v_T > 0.64) is about 2e-15 under this calibration.

---

## Part III. The Computer Side (memory, speed and benchmarking)

This part covers the ideas behind how the program runs fast and how we
measure it. The code itself, file by file and snippet by snippet, is in
[CODE_GUIDE.md](CODE_GUIDE.md).

### 10b. Bytes, MiB and memory pages

Memory is counted in bytes. One byte is 8 bits, and one `double` takes 8
bytes. Bigger units go up in steps of 1024, because memory hardware is built
in powers of two, and 1024 is 2 to the power 10.

| Unit | Bytes | In this project |
|---|---|---|
| byte | 1 | one `double` is 8 bytes |
| cache line | 64 | 8 doubles, the unit the caches load and store |
| KiB, kibibyte | 1024 | 4 KiB is one normal memory page |
| MiB, mebibyte | 1024 KiB = 1,048,576 | one buffer is 8 MiB |
| GiB, gibibyte | 1024 MiB | |
| TiB, tebibyte | 1024 GiB | keeping the grid of every timestep would take 7.5 TiB |

The i in KiB and MiB marks the 1024 version. KB and MB strictly mean 1000 and
1,000,000, so 1 MiB is about 5% more than 1 MB. Chip makers often write MB for
cache sizes when they mean MiB, so the 30 MB L3 cache on an r730 socket is 30
MiB.

The sizes to know by heart.

```text
one row     2048 doubles × 8 bytes          = 16,384 bytes    = 16 KiB
one buffer  2048 × 512 doubles × 8 bytes    = 8,388,608 bytes = 8 MiB
two buffers current and next                                  = 16 MiB
L3 cache    one r730 socket                                   = 30 MiB, so both buffers fit
```

A memory page is the fixed-size chunk the operating system hands out memory
in. A program never sees real RAM addresses. Every address it uses is a
virtual address, and the operating system keeps a page table that maps each
page of virtual addresses to a page of real RAM. On Linux on x86 a normal
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

Each TLB entry covers one page whatever its size, so the page size decides how
much memory the TLB covers. 1024 entries of 4 KiB cover only 4 MiB. 1024
entries of 2 MiB huge pages cover 2 GiB.

| Page size | Pages needed for the 16 MiB of buffers | Fits in 1024 TLB entries? |
|---|---|---|
| 4 KiB, normal page | 4096 | no, so the core keeps doing page walks |
| 2 MiB, huge page | 8 | yes, with room to spare |

That is where the 1.66 times in the page-size trap (§13) comes from. The same
binary on the same grid ran at 134.8 Mcell/s on huge pages and 81.1 on 4 KiB
pages. Nothing in the code changed, only the page size the memory happened to
get. Linux gives a block huge pages by itself (transparent huge pages) when
the block comes fresh from the operating system through `mmap`.
`MALLOC_MMAP_THRESHOLD_=1048576` makes every block over 1 MiB come that way,
so every benchmark sets it.

Pages also decide where memory sits on a two-socket node. Linux places each
page in the memory of the socket whose thread first writes it. This is called
the first-touch policy (§16.8).

### 11. Double buffering

Within one timestep, every new cell must be computed from the *old* buffer
only. The worked example in §7.4 assumed all nine neighbours came from
step 0, and two buffers are how the code makes that true.

**What goes wrong with one buffer.** Suppose we updated the grid in place,
walking each row left to right. By the time we reach the cell at S = 5300,
its west neighbour at S = 5250 has already been overwritten with its step 1
value, and so has the whole row below. Its east neighbour is still at step 0.
The stencil now mixes two moments in time.

```text
  one buffer, walking left to right along a row

     S = 5200     S = 5250     S = 5300     S = 5350
   [ step 1 ]   [ step 1 ]   [  ????  ]   [ step 0 ]
                     ^            ^            ^
                     |       being updated     |
               already new                 still old
```

Running this on the test grid for one step gives three different answers for
the same cell (S = 5300, v = 0.04), depending only on the loop order.

```text
two buffers                     50.0007471
one buffer, left to right       50.0007951
one buffer, right to left       50.0007465
```

After a single step the worst cell is already 3 cents out, and that error
compounds over 56000 steps. With threads it gets worse, because which
neighbours are new would depend on which thread got there first. That is a
race condition, and the answer would change from run to run.

**What two buffers do.** One buffer is only read during a step, and the other
is only written. At the end of the step the two buffers swap pointers, so the
buffer just written is the one read in the next step.

```text
  step 1    read  A (payoff)  -->  update every cell  -->  write  B
            swap: A and B change roles
  step 2    read  B           -->  update every cell  -->  write  A   (payoff overwritten, no longer needed)
            swap
  step 3    read  A           -->  update every cell  -->  write  B
  ...
```

We keep only two buffers because each step needs only the grid from the step
just before it. Keeping the grid of every timestep on the reference grid
would take 980000 copies of 8 MiB each, which is about 7.5 TiB.

The swap at the end of each step exchanges two pointers and copies none of
the million doubles. CODE_GUIDE.md §5.4 shows the code, and §7.2 there
explains why the solver's names for the two buffers stay correct after each
swap.

Nothing ever reads a half-written buffer. Because of this rule the later
parallel versions are correct by design. Threads only read data that nobody is
changing, and each thread writes its own separate cells.

### 12. Memory layout, and why `idx = j*Ns + i`

In the code these have descriptive names. The index is
`index(stock_i, var_j) = var_j * num_stock_nodes + stock_i`, the spacings are
`stock_spacing` and `variance_spacing`, and the buffers are `current()` and
`next()`. The maths shorthand `i`, `j` and `Ns` below means the same thing.

The grid is stored as ONE flat array. A Python `list[list]` would scatter the
rows all over the heap, which is bad for the cache.
This is the layout.

```text
index(i, j) = j * Ns + i      i = stock index, j = variance index
```

This block is the formula `Grid::index` uses to turn a (column, row) pair
into a position in the flat array. We put the row number on the outside
because the inner loop runs along a row, so consecutive cells in the loop are
consecutive in memory.

```text
  the grid as a picture (Ns = 4 columns, 3 rows)      the same grid in memory

  row 2 |  8   9  10  11 |                            index  0  1  2  3 | 4  5  6  7 | 8  9 10 11
  row 1 |  4   5   6   7 |                                   ---row 0-- | --row 1--- | --row 2---
  row 0 |  0   1   2   3 |
          i=0 i=1 i=2 i=3                              index(i, j) = j * 4 + i, so index(2, 1) = 6
```

This picture labels each cell with its position in memory. East and west are
one place apart, while north and south are one row (`Ns` places) apart, which
is why `row_above = row_centre + num_stock_nodes` in the serial code
(CODE_GUIDE.md §8.5).

So an entire variance row (`j` fixed, `i` running) is **contiguous in
memory**. The inner loop walks `i`, so it walks memory in order, which is
what the CPU's caches and prefetcher work best with. This one decision is the
foundation of every optimisation that follows, including SIMD in Milestone 2,
because vector instructions need adjacent data.

This is the cache picture to have for the interview. RAM is slow, about
100 ns per access. Caches are fast, about 1 to 10 ns, but small, and data
moves in 64-byte lines, which is 8 doubles at a time. If the loop reads
memory in order, 7 of every 8 loads hit a cache line that is already loaded.
If it reads random addresses, every load can miss and wait for main
memory.

**The measurements forced a correction here.** This section used to end by
claiming that the cell update is only about 20 flops, so memory access sets
this solver's speed. P7 built the roofline to back that
claim up, and it did not survive.

| | Mcell-updates/s | GFLOP/s (×31 flops/cell) | share of measured scalar peak |
|---|---|---|---|
| baseline, 4 KiB pages | 44.8 | 1.39 | 27% |
| baseline, huge pages | 52.6 | 1.63 | 32% |
| level 6, 4 KiB pages | 81.1 | 2.51 | 49% |
| level 6, huge pages | 134.8 | 4.18 | **82%** |
| level 6 @4096×1024 | 156.9 | 4.86 | **95%** |

The peak in the last column is **measured on this node**. It is 5.117 GFLOP/s
from a multiply-add microbenchmark that stays in registers, compiled exactly
like the solver (`slurm/bench_peak.sh`). The memory side of the roofline is
nowhere near the limit. The optimised kernel asks for about 3 GB/s while
moving 24 bytes per cell, on a machine with four DDR4 channels.

So the accurate statement has three parts.

- **The baseline is limited by division throughput.** It runs at between 42
  and 58 Mcell-updates/s whether the working set is 1 MiB or 64 MiB
  (RESULTS.md Table 4), so memory size makes little difference to it. Five
  `divsd` instructions per cell, at about 13 to 20 cycles each, set its speed.
- **After optimisation the kernel runs almost AT the scalar compute peak.** It
  reaches 82% at the reference grid and 95% at the largest one, against a
  ceiling measured on the same machine. There is almost nothing left for any
  memory technique to recover, which is why levels 3, 4 and 5 came back with
  no gain.
- **Memory still matters at the margin, through latency.** The two memory
  costs P7 actually measured were *page size* (1.66×, from TLB pressure) and
  *traversal order* (2.1×). Both are latency effects. At the reference size
  the whole 16 MiB working set fits inside this node's 30 MB L3 cache. The
  optimised kernel asks for about 3 GB/s from a machine with four DDR4
  channels, so bandwidth is nowhere near the limit.

Layout and loop order still matter, and the flat contiguous array is still the
right decision. It makes SIMD possible in Milestone 2, and the `ctl-order`
control shows that getting it wrong has a real cost. But "memory-bound" was an
assumption taken from the shape of the problem, and it was never measured. The
measurement gives a different conclusion. **Optimising the arithmetic moved
this kernel from division-bound to compute-bound, so the next gain has to come
from doing more work per instruction, which SIMD provides.** That is the
bridge into Milestone 2.

### 12b. Two pictures that make the memory story checkable

Claims like "memory-bound" and "loop order dominates" are easy to make and
easy to doubt. P7 added two tools whose job is to make them checkable. Both
are worth being able to talk through.

**`scripts/memory_anim.py --style cache` is the traffic counter.** Two panels
sweep the *same* nine-point stencil across the *same* flat array in the two
possible orders, with a simulated fully-associative LRU cache of 64-byte lines
running underneath. The picture shows the difference before any number does.
The row-major panel lights up a **horizontal strip three rows high**, because
three rows stay in cache as it walks along them. The swapped-loop panel lights
up a **vertical column**, because three lines stay in cache as it walks down a
16 KB stride.

The headline number is bytes moved per cell update, and the hit rate is
deliberately left out of it. The hit rate barely moves, 98.5% against 92.8%,
because consecutive cells in the swapped order still share six of their nine
accesses vertically. The bytes moved per cell update are **8.5 B against
41.2 B**, a 4.7× difference for identical arithmetic. The roofline uses this
number, and it measures what a cache line costs. A miss moves all 64 bytes
of the line from memory, however few of the eight doubles the loop uses
before the line is evicted.

Be clear about the model's limits on camera. It simulates **capacity only**.
The measured slowdown on real hardware is larger than 4.7×, and the gap comes
from hardware prefetching (a contiguous stream is predictable and a 16 KB
stride is not) and from TLB reach, neither of which the model simulates.
Quoting the model's 4.7× as if it explained the whole measured slowdown would
be overclaiming, and the animation's caption says so.

**`scripts/roofline.py` is the test the memory-bound claim failed.** Flops
per cell are counted **term by term from the level-6 kernel**.
There are 15 for the five finite differences (counting `2.0*V` once, since it
is shared) and 16 for assembling the PDE and stepping, which makes **31 flops
per interior cell**. Bytes are counted the same way, assuming perfect
row-major reuse. Each cell reads one new double of `current` (8 B) and writes
one double of `next` (8 B). There are another 8 B because the store misses
and its line is *read* before being overwritten, which is called
write-allocate and is the detail people forget. That makes **24 bytes per
cell**.

The arithmetic intensity is 31/24 ≈ **1.3 flop/byte**. Achieved GFLOP/s is the
measured cell updates per second times 31. The second ceiling on the plot is
an **empirical** streaming ceiling derived from our own `ctl-order` control.
That kernel uses roughly one double per 64-byte line, so its measured
throughput times 256 B per cell is a bandwidth the machine has been shown to
sustain. It is a lower bound, it is our own measurement, and it needs no
datasheet.

**`scripts/memory_anim.py --style buffers`** covers the other half of §11 and
§13. It shows the double buffering (read `current`, write `next`, then
`std::swap` exchanging two pointers, so no doubles are copied). It then shows
the Milestone 1 ladder's data structures appearing rung by rung, which were
the row constants, the S and S² tables, the row offsets, the row pointers and
the four-wide unroll. It is the clearest way to show that the ladder is
*cumulative*.

### 13. The serial optimisation story (Milestone 1's 25%)

**This section is Milestone 1 history.** `OptSolver`, `--opt-level`,
`solver_opt_kernels.h` and `test_opt_matches` were removed at the start of
Milestone 2 and are kept at git tag `milestone1`. The results below are still
what you quote when asked how the serial rung got its speed. This table says
where each M1 technique is in the M2 code.

| M1 level | Technique | Where it is in the M2 code |
|---|---|---|
| L1 | hoisting, lookup tables | `StencilConstants`, built once before the time loop |
| L2 | reciprocals instead of division | the `inv_*` fields of `StencilConstants` |
| L3 | row-major traversal | `var_j` outer loop in `SerialSolver::solve`, `stock_i` inner loop in `update_cells_scalar` |
| L4 | loop splitting | `update_zero_row` and `apply_stock_edges` are separate from the interior |
| L5 | induction variables | `row_centre`, `row_above`, `row_below`, `out_row` pointers |
| L6 | unrolling by four | dropped, it measured inside run-to-run spread |

`test_ladder_matches` replaced `test_opt_matches` and checks every M2 rung
against the baseline.

The baseline solver is written to be simple on purpose. It recomputes all 9
stencil weights for every cell at every step. It is correct and slow, and it
is the before picture. In Milestone 1, `OptSolver` then applied the
techniques below one at a time, measuring each, and `test_opt_matches` proved
every step gave the same answers. These are the course's own named
techniques (L03), so use these names in the video and interview.

1. **Loop-invariant code motion (hoisting) and a lookup table.** The weights
   depend only on the row (`v`) and the column (`S`), and they do not change
   with time. They are computed into tables before the time loop, so thousands
   of repeated per-cell multiplications become one table lookup.
2. **Strength reduction.** There is no division in the hot loop. A division
   costs 5 to 14 cycles against about 1 for a multiplication, so the
   reciprocals (`1/ds²` and so on) are computed once and multiplied.
3. **Cache-friendly loop order.** The outer loop is `j` and the inner loop is
   `i`, a contiguous walk that matches the memory layout. The lecturer
   measured an 8× difference from loop order alone on a large 2-D array.
4. **Hoisting conditionals out of loops, or loop splitting.** The interior
   loop runs without branches (`i = 1..Ns−2`), and the edges are handled in
   small separate loops, so there is no `if` per cell.
5. **Induction-variable simplification.** The per-cell `j*Ns + i` index
   arithmetic is replaced by a running index that just increments. The
   lecture uses exactly this 2D-to-1D address case as its example.
6. **Loop unrolling** of the inner loop, kept only if it measured faster.

We deliberately did *not* use tiling or blocking (outside the course content),
fast-math (it would weaken the correctness story) or float32 (saved as an M2
experiment, since it halves the memory traffic and doubles the SIMD width).

**How we proved which technique did what is called an ablation ladder.** An
ablation is a controlled experiment borrowed from science. You change ONE
thing, measure, and repeat. `--opt-level 0..6` turned the list above into a
cumulative ladder. Level 0 was the baseline algorithm, level 1 added the
lookup table, level 2 added strength reduction on top, and so on. Each rung
was benchmarked separately, so the bar chart shows the *extra* gain of each
named technique. Two rules kept the experiment fair.

- The level was chosen ONCE, before the time loop, by picking a kernel
  function for that level (`step_level3()` and so on). Checking
  `if (level >= 3)` inside the hot loop instead would slow every level down
  and corrupt the measurement, which is the classic mistake of the instrument
  disturbing the experiment.
- A rung that shows **no speedup is reported anyway**. It usually means the
  compiler already did that optimisation at `-O2`, which is likely for
  induction variables and unrolling. The task sheet explicitly gives credit
  for techniques you tried that didn't work, because a null result shows
  method.

The interview line was this. *"Level N enables the first N techniques
cumulatively. Diffing two adjacent kernels shows exactly one technique, and
`test_opt_matches` proves every level still produces the baseline answer. It
is bit-for-bit through level 1, and within a measured 3.3e-15 relative from
level 2 up, where the reciprocals change the last bits."*

The metric is **cell updates per second**, which is `Ns·Nv·nt / seconds`. It
can be compared across grid sizes, which raw seconds cannot.

**How the Milestone 1 ladder was wired.**

- **The kernel was chosen once, through a function pointer.** `StepKernel`
  (in `solver_opt_kernels.h`) was a *function pointer*, a plain variable whose
  value is which function to call. `OptSolver::solve` asked
  `kernel_for_level(cfg.opt_level)` for it ONCE, before the time loop, so the
  hot loop never ran a per-cell check of the level. It is similar to assigning
  a function to a variable in Python, except that the call compiles to an
  indirect call through the stored address, with no dictionary lookup.
- **A kernel owned the interior and the v = 0 row.** The S = 0 and S = max
  columns and the v = vmax row were shared code in `OptSolver::solve`, the
  same for every level and for the baseline, so the ladder timings compared
  kernels and nothing else.
- **Level 0 was the anchor.** `step_level0` was line for line the baseline
  arithmetic. Its output was bit-identical, and it measured 44.6 against
  44.8 Mcell-updates/s for the baseline on rangpur, which showed the harness
  itself cost nothing.
- The opt solver's CSV label carried the level (`opt-L0`, `opt-L3`). The two
  controls were labelled `opt-ctl-order` and `opt-ctl-branch`, so a benchmark
  row could never be misread as an `opt-L7` or `opt-L8` that continued the
  ladder upward.

#### The measured ladder (rangpur r730-2, Xeon E5-2670 v3, g++ 8.5.0 -O2)

These are 2048×512 runs, median of 5 reps, in the huge-page regime (the page
size section below explains why that matters).

| Level | Technique added | Mcell-updates/s | Step | vs baseline |
|---|---|---|---|---|
| baseline | reference solver | 52.6 | n/a | 1.00× |
| L0 | none (anchor) | 52.6 | −0.1% | 1.00× |
| L1 | hoisting, lookup tables, common subexpressions | 52.7 | +0.3% | 1.00× |
| **L2** | **strength reduction** | **131.2** | **+148.7%** | **2.49×** |
| L3 | traversal order | 131.5 | +0.3% | 2.50× |
| L4 | loop splitting | 131.3 | −0.2% | 2.50× |
| L5 | induction variable | 134.2 | +2.2% | 2.55× |
| L6 | unrolling ×4 | 134.8 | +0.5% | **2.56×** |

**Only one rung out of six gave a measurable gain.** The headline is 2.56× end
to end, and essentially all of it is level 2.

**This is the sentence to have ready for why.** GCC 8.5 at `-O2` had already
done levels 1, 3, 4, 5 and 6 for us. Hoisting loop-invariant expressions,
removing common subexpressions, simplifying induction variables and unrolling
counted loops all leave the result unchanged, so the compiler is free to apply
them and does. **Level 2 is the only rung that changes the result.** Replacing
`x / (2*ds)` with `x * (1/(2*ds))` produces a different double, because
`1/(2*ds)` is itself rounded. Without `-ffast-math` the compiler is *not
allowed* to make that change, however much faster the result would be. The
programmer has to decide the accuracy cost is acceptable and make the change
by hand.

So the rule the ladder demonstrates is that **at `-O2`, the optimisations
worth doing by hand are the ones the compiler is not allowed to do for you.**
This rule is only visible because every rung was measured separately.

The cost is real too, and it is measured. Level 2 moves the answer by 3.3e-15
relative (about 15 ulp) over 5200 timesteps, and the matching test reported
it on every run. The video should state the trade out loud, which is **2.5×
the throughput for 15 units in the last place.**

#### The page-size trap (the most surprising thing P7 found)

The ladder above says huge-page regime because there are two regimes, and
the difference between them is 1.66×.

`Grid` allocates two 8 MiB buffers per solve and frees them at the end of
each `--bench` rep. This is the measured behaviour on this node, reproducible
to three digits.

| | Rep 1 | Reps 2 to 5 |
|---|---|---|
| level 6 @2048×512 | 134.8 Mcell-updates/s | 81.1 |

The two obvious suspects were both checked, and neither is the cause.

- **The clock was ruled out.** `slurm/bench_peak.sh` samples `/proc/cpuinfo`
  *during* an actual solve and finds the clock pinned at 3100 MHz the whole
  way through, so there is no turbo ramp-down.
- **The cache was ruled out.** The working set does not change between reps.

The cause is the **memory allocator**. The first allocation is served by `mmap`,
and with `/sys/kernel/mm/transparent_hugepage/enabled = [always]` it is
backed by 2 MiB transparent huge pages. The 16 MiB working set then needs
about **8 TLB entries**. Afterwards glibc raises its *dynamic* `mmap`
threshold past 8 MiB, so later reps are served from the reused main heap on
4 KiB pages. The same 16 MiB then needs **4096 entries**, far more than this
core's L2 TLB holds, and the stencil's loads trigger page walks.

We proved this by setting the threshold by hand (`slurm/diag_alloc.sh`).

| Setting | Rep 1 | Reps 2 to 5 |
|---|---|---|
| default | 134.8 | 81.1 |
| `MALLOC_MMAP_THRESHOLD_=1048576` (force mmap) | 134.3 | **134.2, the drop disappears** |
| `MALLOC_MMAP_THRESHOLD_=134217728` (force heap) | 134.8 | 81.1 |

**This matters for what gets reported.** A real production solve is a fresh
process that allocates once, so it runs in the *fast* regime. A plain
median of 5 understates real performance by 40% and understates the speedup
(1.81× instead of 2.56×). Both regimes are measured and both are in
`results/bench_summary.md`. The huge-page one is the headline because it is
the one a user experiences.

This gives a good interview answer. **The median of 5 reps is the right
statistic only when the reps sample the same thing.** Here rep 1 and reps 2 to 5
were sampling two different memory regimes, and averaging them would have
hidden both. The fix was to find out *why* the reps disagreed.

#### The negative controls, which show what levels 3 and 4 are actually worth

Levels 3 and 4 could not show a gain because the baseline never had the
loss. So the Milestone 1 ladder shipped two kernels that *do* have the loss,
measured in the same job, at the same `nt`, on the same node. A comparison
across jobs would mix the difference in run length into the result.

| Kernel | What it removes | Mcell-updates/s | vs paired reference |
|---|---|---|---|
| level 6 (paired reference) | n/a | 135.3 | 1.00× |
| `ctl-branch` | the v = 0 row fused back in behind a per-cell `if` | 135.3 | 1.0×, free |
| `ctl-order` | loops swapped, giving an `ns`-element stride | 63.6 | **2.1× slower** |

**Both controls came back smaller than a simple reading of the cache model
would predict, and the reason is the most useful thing this project learned
about the machine.** The reference grid is two 8 MiB buffers, which is
**16 MiB against this node's 30 MB L3**. The whole problem fits in the
last-level cache. So the wrong traversal still finds its data in cache, just
one level further away. Swapping the loops moves roughly 4.8× more bytes
(measured by the cache model in §12b) and costs only 2.1× the time, because
at this size the extra bytes are read from L3, which is much faster than
DRAM.

That gap between *traffic* and *time* is worth saying out loud. A 4.8×
traffic penalty that costs only 2.1× is evidence that the kernel is not
limited by bandwidth here, which is the same conclusion the roofline reaches
from the other direction.

`ctl-branch` being free is expected, for a different reason. The branch is
`var_j == 0`, which is false for 510 of 511 rows, so the branch predictor gets
it right essentially every time. "Hoist the conditional out of the loop" was
good advice for hardware without speculative execution. Today it buys about
1%, and it is *still* worth doing because it costs nothing and makes the loop
bodies simpler.

**Be careful with local numbers.** On the Mac, `ctl-order` measured about 30×
slower. That figure was wrong. The reference solve was running at the same
time, and a bandwidth-hungry kernel suffers most from a neighbour that
competes for bandwidth. The rangpur number was taken on an exclusively
allocated node and is the one that counts. This is why `--exclusive` is in
every benchmark job, and it shows how easy it is to measure something real
that answers a different question from the one you asked.

#### Sweep C, which checks whether any of this is the compiler's auto-vectoriser

| Solver | `-O2` | `-O2 -fno-tree-vectorize` | Vectoriser's share |
|---|---|---|---|
| baseline | 44.8 | 44.8 | 1.00× |
| opt level 6 | 81.1 | 81.1 | 1.00× |

Both columns are default-allocator numbers, because sweep C rebuilt and
re-ran with the allocator left alone. Pairing it against the huge-page
ladder would have credited the vectoriser with the 1.66× page effect. That
pairing mistake is easy to make, and the summary script now guards against
it.

The two columns are identical to three digits. The question matters because
Milestone 1 forbids SIMD and Milestone 2 is *about* SIMD, so the M1 number has
to be a true scalar baseline for M2 to improve on. But the control changed
nothing for a reason that makes it vacuous. GCC 8.5 enables
`-ftree-vectorize` at `-O3` and leaves it off at `-O2`, so there was nothing to
switch off, and RESULTS.md §6 records sweep C as measuring nothing. Sweep G is
the valid control. It built at `-O3` with and without `-fno-tree-vectorize` in
one job, and level 6 ran at 135.0 Mcell-updates/s with the vectoriser off
against 134.8 in the headline sweep, so none of the speedup is hidden
automatic SIMD. The same sweep showed the vectoriser lifting the *baseline*
by 1.92×, because division-heavy code gives it real work.

#### Sweep B, which shows how the speedup scales with problem size

These are the pinned re-run (job 556834, RESULTS.md Table 4), with
`MALLOC_MMAP_THRESHOLD_` set so every rep of a row is served the same way.

| Grid | Working set | Page regime (measured) | Baseline | Level 6 | Speedup |
|---|---|---|---|---|---|
| 512×128 | 1 MiB | 4 KiB | 44.4 | 155.4 | 3.50× |
| 1024×256 | 4 MiB | 4 KiB | 41.8 | 104.4 | 2.49× |
| 2048×512 | 16 MiB | 2 MiB (THP) | 52.6 | 135.1 | 2.57× |
| 4096×1024 | 64 MiB | 2 MiB (THP) | 57.8 | 157.5 | 2.72× |

There are three things to read off this table.

- **The baseline is nearly flat.** It stays between 41.8 and 57.8
  Mcell-updates/s across a 64× range of working-set sizes. A kernel whose
  speed barely changes with where its data sits in the memory hierarchy is
  not memory-bound. The baseline's five divisions per cell set its pace.
- **Level 6 is not flat.** It falls from 155.4 at 1 MiB to 104.4 at 4 MiB,
  where RESULTS.md §3 reads the working set as moving past L2 into the shared
  L3, and recovers to 135.1 and 157.5 on the two large grids, whose buffers
  get huge pages. With the division bottleneck removed, the kernel is fast
  enough to show the memory hierarchy, so the optimisation *moved* the
  bottleneck.
- **The rows still span two page regimes.** Transparent huge pages need an
  aligned 2 MiB stretch of memory, and only the 8 MiB and 32 MiB buffers of
  the two large grids get one. The two small grids therefore run on 4 KiB
  pages whatever the allocator is told, and the within-regime comparisons are
  the clean ones. The first run of this sweep (job 556506) did not pin the
  allocator, mixed the two regimes within single rows and put its 2048×512
  row at 44.8 and 81.1, so its curve read backwards and it is retired to
  `results/superseded/`. It shows why an unexplained bump in the data is
  worth chasing.

There are two smaller observations, worth a sentence each.

- **The Mac and the Xeon disagree about unrolling.** On Apple Silicon (through
  Rosetta) level 6 was a *23% slowdown*. Four cells in flight need more live
  registers than were available there, so the compiler had to store some
  values to the stack and load them back (register spilling). On the Xeon it
  was +0.5%, which is no change. The same source gave opposite results,
  because the effect of an optimisation depends on the code, the compiler and
  the machine together. This is exactly why the rubric insists numbers come
  from the target machine.
- **Level 1 gave no gain, and it was still needed.** It is the rung that
  makes levels 2 to 6 *possible to write*. The reciprocals of
  level 2 have to be hoisted somewhere, and the S and S² tables are what let
  the inner loop become a pure streaming loop. The compiler reproduced its
  speed effect, but it did not write the code that follows.

### 14. Benchmarking properly (the other 25%)

- Numbers come **only from rangpur compute nodes** through `sbatch`. The login
  node is shared and busy, so numbers from it are meaningless.
- Run at least 5 repetitions and report the **median** and the spread, never
  a single run.
- Record the hostname, compiler version, flags and grid size with every
  result.
- Time only the solve loop, never file input or output.
- Flags are fixed at `-O2`. A run at `-O3` with and without
  `-fno-tree-vectorize` (sweep G) shows the serial gains aren't hidden
  automatic SIMD. M1 forbids SIMD, and the lectures describe `-O3` as the
  level that enables vectorisation, which is another reason we stay at `-O2`.
- **Profile before optimising.** With `gprof` you compile with `-pg`, run the
  program, then run `gprof prog gmon.out`. Find where the time actually goes
  and optimise that. Don't polish code that takes 1% of the runtime.
- **Use your results.** Print the price. The lecturer showed a benchmark whose
  runtime dropped to zero because the compiler noticed the result was never
  used and deleted the whole loop.
- **The problem should be larger than the cache**, otherwise the benchmark
  measures the cache. Our reference grid is two 8 MiB buffers, so say so.
  §13 shows that those 16 MiB do fit in this node's 30 MiB L3, and why that
  matters.
- **Check correctness while benchmarking.** It is very easy to make code run
  fast that doesn't work, so our benchmark job runs `make test` first.
- **Report null results.** If an optimisation doesn't help, the compiler may
  already have done it, and saying so shows method.

**What the benchmark machine actually is, and why it matters (P7).** Rangpur's
`cosc3500` partition mixes 24-core Intel Xeon E5-2670 v3 boxes (`r730-*`,
Haswell, 2014) with 8-core virtual machines (`a100-*`, `vcpu-*`). Three things
followed from that, and all of them are worth being able to explain.

- **Every benchmark job pins to the same physical node**
  (`--nodelist=r730-2`) **and takes it exclusively** (`--exclusive`). Without
  the pin, the ladder might be measured on a Xeon and the controls on a
  virtual machine, and comparing them would be meaningless. Without
  `--exclusive`, another job on a different core of the same socket competes
  for the *shared* memory bandwidth and L3 cache. It never runs on our core,
  and it still changes our numbers.
- **The QOS caps a job at 15 minutes.** A single job covering everything sits
  in the queue forever with reason `QOSMaxWallDurationPerJobLimit`. Slurm
  never rejects the job and never starts it, which is a much more confusing
  failure. The Milestone 1 suite was five sweeps plus two diagnostic jobs,
  each writing its own CSV (removed in Milestone 2, still at git tag
  `milestone1`). The M2 sweeps follow the same one-job-per-sweep rule, and
  COMMANDS.md lists them.
- **That node is slow, and the reason is informative.** It runs the baseline
  at about 4.7e7 cell updates per second, roughly a tenth of a modern laptop
  core, at a measured **3.1 GHz**. `bench_peak.sh` sampled `/proc/cpuinfo`
  *during* a real level-6 solve and found the clock pinned at 3100 MHz
  throughout, so this is not throttling. An earlier draft of this section
  said 1.6 GHz, which was the idle reading taken before that job existed, and
  the cycles-per-cell figure here is corrected for it. That works out to
  about **66 cycles per cell update**. The baseline does **five divisions per
  cell**, and a Haswell `divsd` takes about 13 to 20 cycles and pipelines
  poorly, so those five alone account for 65 to 100 cycles, which is the
  whole budget. So the *baseline* on this machine is limited by division
  before anything else. That is why strength reduction (level 2) gives such a
  large gain here, and the corrected clock makes that argument stronger.

**A useful debugging habit is to check a surprising number before believing
it.** When rangpur first came back 10× slower than the laptop, the obvious
suspect was denormal arithmetic. `bench.cfg` uses a very short maturity, so
much of the grid holds values near 1e-300, and x86 takes a microcode assist
of tens of cycles on every operation that reads or produces a denormal. That
would have made the slow speed a benchmarking artefact. So it was tested with
`slurm/diag_denormal.sh`. The same run was repeated
with `-ffast-math` (which links `crtfastmath.o` and sets FTZ/DAZ so denormals
flush to zero), with the grid full of NaNs, and with values around 1e33. All
of them gave the *same* throughput to within noise. The hypothesis was wrong
for the Milestone 1 kernels, and the hardware is old, which is now a
measured statement. That test was also the reason the short `nt` was trusted,
since for those kernels a 400-step run measured the same *rate* as a
400 000-step one. Milestone 2 overturned this for the faster rungs.
Unflushed denormals cost the baseline 1.07× and the simd rung 1.42×
(RESULTS.md M2.7), so every M2 benchmark flushes them to zero (§18.5).

**Self-test.** Is the solver memory-bound or compute-bound, and what
measurement settles it? Why is `j*Ns + i` the right index formula for a loop
whose inner index is `i`? What does the buffer swap cost? Why take the median
of 5 runs on a compute node? Why does `--exclusive` matter when we only use
one core? What is the difference between a job that is *rejected* and a job
that never starts?

---

## Part IV. The Parallel Future (Milestone 2 preview, know the shape)

You will be asked where the parallelism is, even about Milestone 1. This is
the answer.

**Within one timestep the work is fully parallel.** Every cell of the new
time level depends only on the old one, so all `Ns·Nv ≈ 1M` cells could be
computed at the same time, in any order.

**Across timesteps the work is strictly sequential.** Step k+1 can't start
until step k is complete. Every one of the `nt` timesteps is a point where
all workers must wait for each other.

That combination, with a lot of parallelism inside a step and a barrier
between steps, makes this a *real* parallel-computing problem. Monte Carlo is
embarrassingly parallel because its paths never exchange data, while this
solver has to synchronise every step.

The M2 ladder uses it rung by rung, with the same `Solver` interface each
time, using only what the course teaches (L04 and L05).

- **SIMD uses AVX intrinsics.** These come from `<immintrin.h>`, a `__m256d`
  holds 4 doubles, and FMA does the multiply-adds. The vector code is written
  by hand with intrinsics, which sit between compiler auto-vectorisation
  (which the lectures call "pretty subpar") and assembly. Adjacent cells in a
  row share the same formula and sit next to each other in memory, so one
  vector instruction computes 4 cells. The plan uses unaligned loads
  (`_mm256_loadu_pd`), because the stencil reads `stock_i - 1`, `stock_i` and
  `stock_i + 1`, so at most one of those three loads could ever be aligned.
  Experiment E2 measured whether alignment helps anyway (§15.4 and §15.5). The
  lectures' rule of thumb is **about 4×, below the full vector width**,
  because such kernels are usually memory-bound. **Our measurements show this
  kernel is compute-bound** (§12). The M1 optimised kernel already runs at 82%
  of this node's measured scalar peak while asking for only about 3 GB/s.
  Arithmetic throughput limits its speed, and SIMD raises arithmetic
  throughput. So the plan expected a gain closer to the vector width than the
  rule of thumb predicts, and the measured gain over serial was 2.87×
  (RESULTS.md M2.10). A float32 version (8 wide) became experiment E3. §15 has
  how the rung was built and checked.
- **OpenMP puts `#pragma omp parallel for` on the outer `var_j` loop.** The
  lecture rule is to parallelise the outermost loop first. It uses
  `default(none)` with explicit sharing lists, as the lecture recommends, and
  the thread count comes from `OMP_NUM_THREADS`. Our design needs **no
  synchronisation inside a step**. Threads read the shared `current` buffer,
  which nobody changes during the step, and write separate rows of `next`.
  The only wait is the implicit barrier at the end of the parallel loop. In
  the lecture's warning example, putting `critical` in a hot loop made the
  code **40× slower** than serial, and our structure, which needs no locks,
  is the answer to how we avoided that. Expect a speedup close to the number
  of physical cores (the r730 has no hyper-threading). Each cell is computed
  by the same arithmetic in the same order whichever thread runs it, and
  there are no sums across threads, so the omp rung matches simd bit for bit
  (RESULTS.md M2.1).
- **MPI is hybrid with OpenMP.** Several processes, possibly on different
  nodes, each own a block of rows. A process cannot read another's memory,
  so before each step each process sends its first and last rows to its
  neighbours and receives theirs into its copied rows.
- **CUDA is optional.** It would use one GPU thread per cell, with both
  buffers kept on the GPU so nothing crosses the PCIe bus inside the time
  loop.

**How each rung splits the grid.** Each phase adds the details. This is the
shape to have in your head.

```text
  SIMD, inside one row            OpenMP, rows to threads       MPI, blocks of rows to processes

  row var_j                       row 7  thread 3               +---------------------+
  [c c c c][c c c c][c c c c]     row 6  thread 3               | rank 1  rows 4..7   |
   4 cells per instruction        row 5  thread 2               | copy of row 3       |
                                  row 4  thread 2               +---------------------+
                                  row 3  thread 1               | copy of row 4       |
                                  row 2  thread 1               | rank 0  rows 0..3   |
                                  row 1  thread 0               +---------------------+
```

This picture shows that every rung divides up the same rectangle from §6.1.
It works without locks because of the two buffers in §11. Every worker reads
`current`, which nobody writes during a step, and writes only its own rows of
`next`. One copied row on each side is enough because the stencil reaches
exactly one row up and one row down (§7.2).

Have two laws ready. **Amdahl's law** says the speedup is capped by the serial
fraction of the work, and our per-step barrier is exactly such a fraction. The
lecture example is that 90% parallel code can reach at most 10×. **Gustafson's
law** says that if you scale the *problem* up with the machine, efficiency
recovers, so finer grids justify more cores. The task sheet names Gustafson's
law, but the lectures do not teach it, so the video cites the paper, J. L.
Gustafson, Reevaluating Amdahl's law, Communications of the ACM 31(5), 532 to
533, 1988. With serial share f of the scaled run on N cores, the scaled
speedup is f + (1 - f)N. Our measured f = 0.092 on 24 threads gives 21.9,
against 16.3 measured (RESULTS.md M2.3).

The video stays with what the course teaches. Where a slide uses an idea from
outside the lectures, it carries a source line, and M2_SCRIPT.md lists them.
The roofline is from Williams, Waterman and Patterson (2009). Fitting
Amdahl's serial fraction to measured speedups is the Karp and Flatt metric
(1990), which `amdahl_fit` in `scripts/m2_figures.py` does by least squares
over every thread count. Page placement on NUMA nodes is from Drepper (2007).
Denormal numbers are from Goldberg (1991).

---

## Part IVb. Milestone 2, rung by rung

Each section below was written when its rung was finished, and says what the
rung does, why it works, how it was checked, and what it measured. The code
itself is walked through in CODE_GUIDE.md.

### 15. SIMD with AVX2 (the `simd` rung)

#### 15.1 What SIMD is

SIMD stands for single instruction, multiple data. A normal (scalar)
instruction adds one pair of numbers. A SIMD instruction adds several pairs in
one instruction, held side by side in one wide register. The CPU has to have
the wide registers and the instructions, and they come in generations.

| Name | Register | Width | Doubles per register | Register names |
|---|---|---|---|---|
| SSE2 | 128 bit | 16 bytes | 2 | `xmm0` to `xmm15` |
| AVX, AVX2 | 256 bit | 32 bytes | 4 | `ymm0` to `ymm15` |
| AVX-512 | 512 bit | 64 bytes | 8 | `zmm0` to `zmm31` |

The rangpur r730 nodes (Xeon E5-2670 v3, Haswell) have AVX2 and FMA but not
AVX-512. So the widest useful register holds four doubles.

In C++ a 256-bit register of four doubles has the type `__m256d`, and every
operation on it is a function from `<immintrin.h>` called an intrinsic. Each
intrinsic becomes one machine instruction.

| Intrinsic | Does, in every lane | Python with numpy, on 4-element arrays |
|---|---|---|
| `_mm256_set1_pd(x)` | fill all four lanes with `x` | `np.full(4, x)` |
| `_mm256_loadu_pd(p)` | load `p[0]` to `p[3]` | `a = p[i:i+4]` |
| `_mm256_storeu_pd(p, a)` | store four lanes to `p[0]` to `p[3]` | `p[i:i+4] = a` |
| `_mm256_add_pd(a, b)` | `a + b` | `a + b` |
| `_mm256_sub_pd(a, b)` | `a - b` | `a - b` |
| `_mm256_mul_pd(a, b)` | `a * b` | `a * b` |
| `_mm256_fmadd_pd(a, b, c)` | `a * b + c`, rounded once | `a * b + c` |
| `_mm256_fnmadd_pd(a, b, c)` | `c - a * b`, rounded once | `c - a * b` |

The names read as `_mm256` (256-bit), then the operation, then `_pd` (packed
doubles, meaning four doubles side by side). The
[Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html)
lists every one.

SIMD is the opposite of what the next two rungs do. OpenMP and MPI are MIMD
(multiple instruction, multiple data). Each core runs its own instruction
stream on its own data. SIMD happens inside one core, and the two multiply,
so 24 cores each doing four cells per instruction is 96 cells in flight.

#### 15.2 Why our loop vectorises cleanly

Three facts about the interior row update make it a good fit.

1. Every cell in a row runs exactly the same formula. There is no `if` that
   depends on the cell, so four lanes can follow one instruction stream.
2. Neighbouring cells in a row are neighbours in memory (`idx = j*Ns + i`,
   §12). Four consecutive cells are one 32-byte load.
3. Within one step no cell reads another cell's new value (§11). Lane 2
   never needs lane 1's answer, so the four lanes are independent.

The per-row coefficients (`stock_diffusion`, `cross_diffusion`,
`variance_diffusion`, `variance_drift`) are the same for every cell in the
row, so each is filled into a register once per row with `set1`. S and S²
change with the column, so they are loaded from the tables four at a time
like the grid values.

#### 15.3 How the nine neighbours are loaded

For four cells starting at column `i`, the AVX row loads nine registers, one
per neighbour direction. Each is the same four-wide load shifted by one
column or one row.

```text
                    columns    i-1   i    i+1  i+2  i+3  i+4
  row_above (j+1)              [ north_west          ]
                                    [ north              ]
                                         [ north_east          ]
  row_centre (j)               [ west                ]
                                    [ centre             ]
                                         [ east                ]
  row_below (j-1)              [ south_west          ]
                                    [ south              ]
                                         [ south_east          ]
```

Lane 0 of `east` is the east neighbour of lane 0 of `centre`, lane 1 of
`east` is the east neighbour of lane 1, and so on. So once the nine loads are
done, the rest of the formula is the scalar formula written with intrinsics,
and every line does four cells. The loads overlap, which means most values
are loaded more than once, but they come from L1 cache and the CPU can do
two loads per cycle.

#### 15.4 Why unaligned loads

An aligned load (`_mm256_load_pd`) needs the address to be a multiple of 32
bytes, and crashes otherwise. Our west, centre and east loads start at
`i - 1`, `i` and `i + 1`, which are 8 bytes apart. At most one of any three
addresses 8 bytes apart can be a multiple of 32. So the three loads cannot
all be aligned, and the shipped row uses `_mm256_loadu_pd` everywhere. On
Haswell and later, an unaligned load that happens to be aligned runs at full
speed, and one that crosses a 64-byte cache line costs a little more.

Experiment E2 measured the alternatives (RESULTS.md M2.6, the E2b rows with
unsplit loads, one thread). The aligned load instruction ran at exactly the
same speed as the unaligned one at the same address (477 and 477
Mcell-updates/s). Where the data sits still matters. A row whose centre load
is 32 byte aligned ran at 477, one aligned on the west or east load at 439,
and one with nothing aligned at 412, because every unaligned 32 byte load
crosses a cache line on every other iteration. Building east and west from
the aligned centre by register shuffles (460) was no faster. §15.5 covers the
compiler setting that E2 found first.

#### 15.5 What FMA changes

FMA is fused multiply-add. `_mm256_fmadd_pd(a, b, c)` computes `a * b + c`
as one instruction. It has two effects.

- **Speed.** A multiply and an add cost one instruction instead of two. In
  the AVX row eight multiply and add pairs became single FMA instructions
  (two in the second derivatives, five in summing the terms, one in the
  final step), out of 24 arithmetic instructions per four cells.
- **Rounding.** A normal multiply rounds `a * b` to the nearest double, then
  the add rounds again. FMA keeps the exact product and rounds only once at
  the end. The FMA answer is usually slightly more accurate, but it is
  different in the last bits, so the SIMD rung is not bit-identical to the
  serial rung.

That second effect is why the AVX2 flags are given to one file only. With
`-mfma` set, g++ is allowed to fuse any `a * b + c` it finds into an FMA, in
any file. If every file were built with `-mfma`, the baseline and serial
rungs would quietly change their answers and their speed too, and the
ladder would no longer measure one change per rung. So the Makefile builds
`stencil_avx.cpp` with `-mavx2 -mfma` and every other file without.

We use `-mavx2 -mfma` and leave out `-march=native`, because `-march=native`
enables every feature of the build machine, so a binary built on a newer node
could use AVX-512 and crash on the r730s. The check for this is `objdump -d
heston | grep zmm`, which must print nothing, because `zmm` registers only
exist in AVX-512. Some AVX-512 intrinsics also have names starting with
`_mm256_`, so the name alone does not prove the code is AVX2.

The AVX2 file gets two more flags, `-mno-avx256-split-unaligned-load` and
`-mno-avx256-split-unaligned-store` (`NOSPLIT` in the Makefile). Given only
`-mavx2`, g++ 8.5 tunes for a generic x86 and compiles every unaligned 256 bit
load as two 128 bit loads joined by `vinsertf128`, and every unaligned store
as two 128 bit stores. Experiment E2 found the split in the disassembly. On
the r730s' Haswell cores the split version is slower, and E2b measured simd
1.19 times and omp 1.12 times faster with the two flags, with the same
answers. `-march=native` would also keep the 256 bit loads unsplit, because it
tunes for Haswell, but it ties the binary to the build machine for the reason
above.

#### 15.6 The remainder

A row has `ns - 2` interior cells. The AVX loop takes them four at a time
while four whole cells are left, and the last zero to three are done by the
scalar function `update_cells_scalar`, the same code the serial rung uses. For
the 2048-wide grid that is 511 groups of four and two scalar cells per row,
so the remainder is 0.1% of the work.

#### 15.7 How it was checked

The Mac is arm64 and has no AVX2, so there are three layers.

1. **On the Mac, the fallback.** `stencil_avx.cpp` is wrapped in
   `#if defined(__AVX2__) && defined(__FMA__)`. Without the flags the
   compiler takes the `#else` branch, where `update_row_avx` calls the scalar
   row. So `make test` runs on the Mac and the `simd` rung gives the serial
   answer, with a warning on stderr.
2. **On the Mac, under Rosetta.** `make avx-check` builds an Intel (x86_64)
   copy of the row test with `clang++ -arch x86_64 -mavx2 -mfma`. macOS runs
   it through Rosetta, which translates AVX2 into ARM instructions. The
   answers are valid, and the timings are meaningless.
3. **On rangpur.** `slurm/debug.sh` builds with the real flags, runs
   `make test` (which includes `simd` in `test_ladder_matches`) and
   `make avx-check`, and counts `zmm` instructions.

| Check | Where | Result |
|---|---|---|
| `test_avx_row`, AVX2 row vs scalar row, random buffer, row lengths 130 to 133 and 2048 | Rosetta and rangpur | worst gap 2.3e-15 of the largest value, tolerance 1e-13 |
| `test_ladder_matches`, simd vs baseline, 5200 steps, price and three Greeks | rangpur job 651064 | 1.06e-15 call, 0 put, tolerance 1e-11 |
| `test_ladder_matches`, serial vs baseline, same run | rangpur job 651064 | 3.3e-15 call, 0 put, unchanged |
| `zmm` instructions in `heston` | rangpur job 651064 | 0 |
| `ymm` instructions in `stencil_avx.o` | rangpur job 651064 | 57 |

The row test was also checked against a planted bug. Changing the
discounting line from `fnmadd` to `fmadd` made it fail, so it does catch a
wrong formula.

Rosetta and rangpur gave exactly the same numbers, which is expected because
both run the same x86_64 instructions on IEEE doubles.

#### 15.8 What we expect it to be worth

PLAN §5 hypothesis 1 predicts 2 to 3.5 times over the serial rung. The M1
serial kernel was compute-bound at 82 to 95% of the node's scalar peak, so a
four-wide register with FMA should carry most of its width through. The
debug job's smoke run on a shared `vcpu-2` node showed about 2.4 times. The
benchmark on r730-2 (job 659238, denormals flushed, §18.5, loads kept whole)
measured serial 159 and simd 456 Mcell-updates/s, 2.87 times, inside the
predicted range. §18.6 explains why it falls short of 4.

#### 15.9 Float instead of double (experiment E3)

A `double` is 8 bytes and keeps about 16 significant digits. A `float` is 4
bytes and keeps about 7. The gap between 1 and the next number up is called
machine epsilon, and it is about 1.1e-16 for double (2^-53, half the gap,
is the worst rounding) and about 6e-8 for float (2^-24).

The same 256-bit register holds four doubles (`__m256d`) or eight floats
(`__m256`), so a float row does twice the cells per instruction, and each
buffer takes half the memory. That gives two ways it could be faster. On one
thread, where the kernel is compute-bound, a float register has eight lanes
where a double register has four. On 24 threads on a large grid, where the
kernel is memory-bound, float moves half the bytes per cell.

The cost is accuracy, and it is not the same on every grid. Each time step
rounds every cell once, so the error builds up step after step. If each
rounding is random, after `nt` steps the error is about `sqrt(nt)` roundings,
and if they all lean the same way it is `nt` roundings. The more steps, the
worse it gets, and the reference grid has 980000.

```text
grid        steps     price error from float    gamma error from float
bench        2000     2e-7                      7e-7
demo        24000     9e-5                      1.6e-4
reference  980000     2.2e-4                    11.9%
```

Each row compares `omp_float` with `omp` on the same grid. The bench and
demo rows are local runs on the Mac, and the reference row is from rangpur
job 652647, where delta moved by 0.46% and vega by 2.3% as well. Gamma
suffers most because it is a second difference, `(east - 2 centre + west)
/ ds^2`, which subtracts nearly equal numbers and then divides by a small
one, so the rounding error in each value is multiplied many times over.

The Fourier price for the reference contract is 196.1692. The double
solver gives 196.1684, $0.0008 off. The float solver gives 196.1256, $0.044
off. That is worse than the double answer on the much coarser 841 by 65
grid with 224000 steps, which is $0.026 off (RESULTS §7). So float loses
the accuracy that the fine grid's extra cost was spent on.

And the speed, from rangpur job 659247 on r730-2, in Mcell-updates/s.

```text
run                          double     float    float / double
1 thread, 2048 x 512            453       871    1.92x
24 threads, 2048 x 512         6375     11043    1.73x
24 threads, 4096 x 1024        2975     14353    4.82x
24 threads, reference grid     6856     11676    1.70x
```

On one thread the gain is close to the full 2x, because float doubles the
lanes per instruction and halves the bytes per cell at the same time. On 24
threads at 2048 by 512 both versions fit in the 60 MB of L3 (16 MB of
buffers for double, 8 MB for float), and the gain is 1.73x. At 4096 by 1024
the double buffers take 64 MB and no longer fit in L3, so they are read from
main memory, while the float buffers take 32 MB and still fit in the two L3
caches. So float halves the bytes and also keeps the whole grid in cache, and
the gain rises to 4.82x. Figure F4 shows the same drop for the double grid
as it grows.

So float is a real speedup, and on a large grid a big one. It was not taken
into the ladder, because the price it gives on the reference grid is
further from the true answer than a far cheaper double grid gets. The
error comes from storing every cell as a float after every step, so doing
only the arithmetic in double would not remove it.

#### 15.10 Self-test

- **What is a `__m256d`?** A 256-bit register holding four doubles side by
  side. One AVX instruction operates on all four lanes.
- **Why can the row be vectorised at all?** Every cell runs the same formula,
  neighbouring cells sit next to each other in memory, and no cell reads
  another cell's new value inside a step.
- **Why does the row use `loadu`?** The west, centre and east loads start 8
  bytes apart, so at most one of the three per row can be 32-byte aligned,
  and the aligned load crashes on the others.
- **What does FMA change?** It does a multiply and an add as one instruction
  with one rounding. That is faster, and it moves the last bits of the
  answer, so the SIMD rung matches the baseline to about 1e-15 and is not
  bit-identical to it.
- **Why are only some files built with `-mfma`?** With `-mfma` the compiler
  may fuse multiply-adds in any file, which would change the baseline and
  serial rungs too.
- **Why not `-march=native`?** It could enable AVX-512 on a newer build
  machine, and the binary would then crash on the r730 nodes. `objdump` shows
  no `zmm` registers, which proves there is no AVX-512.
- **What happens to the last few cells of a row?** The cells that do not
  fill a register are done by the scalar function the serial rung uses.
- **Why are the boundaries not vectorised?** They are a few thousand cells
  per step against about a million in the interior, so they are not where
  the time goes.
- **How did you test AVX2 code on a Mac that does not have it?** An x86_64
  build of the row test runs under Rosetta, then the full suite runs on
  rangpur. Both gave the same numbers.
- **Why not use float and get eight lanes?** It was measured (E3). Float
  rounds at about 6e-8 instead of 1e-16, and an explicit scheme rounds every
  cell every step, so on the 980000 step reference grid the price moved by
  $0.044 and gamma by 12%. That is more error than the grid itself has.

---

### 16. OpenMP threads (the `omp` rung)

#### 16.1 Threads and shared memory

A process is one running program with its own memory. A thread is one stream
of instructions inside a process. All the threads of a process share the same
memory, so every thread can read and write the same `current` and `next`
vectors without copying anything. Each thread gets its own stack, so local
variables inside the loop body (like `centre` and `dV_dS`) are private to the
thread that runs them.

Python has threads too, but the global interpreter lock lets only one of them
run Python code at a time. C++ threads have no such lock, so 24 threads on 24
cores really do run at the same time.

OpenMP is a set of compiler pragmas and a small runtime library that creates
and manages the threads for you. You mark a loop, and the compiler turns it
into code that splits the iterations across threads. The program must be
built with `-fopenmp` (g++) or `-Xpreprocessor -fopenmp` plus libomp (Apple
clang). GCC 8.5 on rangpur implements OpenMP 4.5, whose
[specification](https://www.openmp.org/wp-content/uploads/openmp-4.5.pdf) is
the reference for every clause used here.

#### 16.2 Fork and join, once per timestep

```text
  step k    main ──┬── thread 0: rows   1..128 ──┬── main: boundaries, swap
                   ├── thread 1: rows 129..256 ──┤
                   ├── thread 2: rows 257..383 ──┤
                   └── thread 3: rows 384..510 ──┘
                 fork                          join (barrier)
  step k+1  main ──┬── ...
```

At `#pragma omp parallel for` the main thread wakes a team of threads (the
fork). Each does its rows. At the end of the loop every thread waits at an
implicit barrier until the last one finishes (the join). Then only the main
thread carries on and does the boundaries and the swap.

The runtime keeps the threads alive between loops, so a fork does not create
new operating system threads every step. It still costs a few microseconds
to wake the team and to wait at the barrier, and it happens once per step.
Experiment E4 measures whether keeping one parallel region open for the
whole time loop is faster.

#### 16.3 The pragma, word by word

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
```

| Part | Meaning |
|---|---|
| `parallel` | start a team of threads for the next statement |
| `for` | split the iterations of the following `for` loop across the team |
| `schedule(static)` | cut the iterations into one equal contiguous block per thread, decided before the loop starts |
| `default(none)` | every variable used inside must have its sharing stated, or the compile fails |
| `shared(...)` | these variables are the same object for every thread |

The loop counter `var_j` is made private automatically, so each thread has
its own copy. `default(none)` is there as a safety check. Without it, any
variable you forgot about is shared by default, and a variable that should
have been private silently becomes a data race. With it, forgetting a
variable is a compile error.

`schedule(static)` suits this loop because every row has the same number of
cells and the same formula, so equal blocks mean equal work. `dynamic`
scheduling hands out rows on demand, which helps when iterations cost
different amounts but adds a little bookkeeping per chunk. Experiment E5
measures the difference.

When the rows do not divide evenly, the first threads get one extra row each.
510 interior rows over 4 threads gives 128, 128, 127 and 127. Both libgomp
(g++) and libomp (clang) split it this way.

#### 16.4 Why the threads split rows

Each thread gets a contiguous block of whole rows. That has three good
properties.

1. It parallelises the outermost loop, so there is one fork per step, where
   parallelising the inner loop would fork once per row.
2. Each thread's rows are one contiguous block of memory, so each thread
   streams through its own part of the arrays.
3. The inner loop along a row is unchanged, so it is still the AVX2 row from
   §15 and the two speedups multiply.

Splitting along columns instead would cut every row into short pieces, would
need a fork inside the row loop or a 2-D split, and would break the long
contiguous AVX2 loop.

#### 16.5 Why no locks are needed

A data race is when two threads access the same memory location at the same
time and at least one of them writes. Our loop has none.

- Every thread reads `current`, but nobody writes `current` during the loop.
  Reading the same memory from many threads is always safe.
- Every thread writes `next`, but only its own rows of it, and no two threads
  have the same row.
- `constants` is only read.

The two-buffer design from §11 helps here a second time. With a
single buffer, a thread would read a neighbour row that another thread might
already have overwritten, and the answer would depend on timing.

So there is no `critical`, `atomic` or lock in the hot loop. The lecture
example where `critical` inside a loop made the code about 40 times slower
than serial is the reason to avoid needing them. The only synchronisation is
the barrier at the end of the loop, once per step.

#### 16.6 Why the answer is bit-identical to the SIMD rung

The lectures warn that a parallel sum can give a different answer from a
serial one, because floating-point addition is not associative and threads
add in a different order. That applies to a reduction, where several
threads add into one total. This loop has no reduction. Every cell is
computed by exactly the same `update_row_avx` instructions it would run in
the SIMD rung, only on a different core. So the OpenMP answer is the SIMD
answer bit for bit, at any thread count.

`test_ladder_matches` checks this with a tolerance of exactly zero at 1, 2
and 7 threads. Seven threads do not divide the 30 interior rows of the test
grid, so the blocks are uneven. If a row were skipped or computed twice at a
block edge, the result would change, and a tolerance like 1e-11 could hide
that.

#### 16.7 The serial part, and Amdahl's law

The v = 0 row, the two stock edges and the top row are done by the main
thread after the barrier. That is about 2048 + 2 × 512 + 2048 cells per step
against about a million interior cells, roughly 0.5% of the work. The fork,
the barrier and the swap are also serial.

Amdahl's law says that if a fraction `f` of the run cannot be parallelised,
the speedup on `p` threads is at most

```text
  S(p) = 1 / ( f + (1 - f) / p )
```

With f = 0.005, 24 threads could give at most about 21.5 times. The strong
scaling sweep measures S(p) and fits `f`. The fitted `f` will include
everything that does not shrink with more threads, such as the barrier and
memory bandwidth limits, so it is expected to be larger than the 0.5% count
of boundary cells. Experiment E8 parallelises the boundaries to see if that
term is even visible.

#### 16.8 Thread count and placement

`OMP_NUM_THREADS` sets the team size. The rung records what it actually got
with `omp_get_max_threads()`, and that number goes in the CSV.

Setting more threads than cores is called oversubscription. The threads then
take turns on the cores, and since every step ends in a barrier, the whole
team waits for whichever thread was paused. So oversubscription makes the
run slower.

The r730 nodes have two sockets of 12 cores. Each socket has its own 30 MiB
L3 cache and its own memory. Two settings decide where threads run.

| Variable | Values | Meaning |
|---|---|---|
| `OMP_PLACES` | `cores` | a place is one physical core |
| `OMP_PROC_BIND` | `false` | the operating system may move threads around |
| | `close` | thread n goes on the place after thread n-1 in the place list |
| | `spread` | threads are spread out evenly over the place list |

What that means in sockets depends on how the machine numbers its CPUs. On
the r730 nodes the even CPUs are socket 0 and the odd ones socket 1, and the
`cores` place list is in CPU order, so `close` alternates sockets from the
second thread on, and `spread` with 6 or 12 threads uses only even CPUs,
all on socket 0. That is the opposite of the usual picture, and job 651500
(`slurm/diag_placement.sh`) measured it by printing each thread's CPU.

Linux places each page in the memory of the socket whose thread first writes
it (§10b). `Grid` fills both buffers on the main thread, so all the pages are
placed in socket 0's memory, and threads on socket 1 read them across the
link between sockets. Experiments E6 (binding) and E7 (page placement)
measure how much this matters, which should only be visible once the grid no
longer fits in L3.

#### 16.9 How it was checked

| Check | Where | Result |
|---|---|---|
| `test_ladder_matches`, omp vs baseline | Mac, 1, 2, 7 threads | 2.35e-15 call, 2.24e-15 put, tolerance 1e-11 |
| `test_ladder_matches`, omp vs simd | Mac, 1, 2, 7 threads | exactly 0 |
| GCC 8.5 accepts the `default(none) shared(...)` list | rangpur login node, compile only | no errors or warnings |
| `make test` with real AVX2 at 1, 2, 7 threads | rangpur debug job 651077 | all pass, omp vs simd exactly 0 at 1, 2, 7 and 8 threads, omp vs baseline 1.06e-15 call and 0 put |

#### 16.10 Self-test

- **What is the difference between a thread and a process?** Threads of one
  process share its memory, processes do not. OpenMP uses threads, MPI uses
  processes.
- **What does `default(none)` do and why use it?** It forces every variable
  in the loop to have its sharing stated, so a forgotten variable is caught
  as a compile error.
- **Why `schedule(static)`?** Every row costs the same, so equal blocks
  decided up front give equal work with no runtime bookkeeping.
- **Why split rows across the threads?** It is the outer loop, so there is
  one fork per step, each thread gets contiguous memory, and the AVX2 inner
  loop is kept.
- **Why are no locks needed?** Threads only read `current`, which nobody
  writes during the step, and each writes only its own rows of `next`.
- **Why is the answer exactly the same as the SIMD rung?** There is no
  reduction. Each cell runs the same instructions, only on a different core.
- **What is the barrier for?** The boundaries and the swap must not start
  until every interior row is written, especially `copy_top_row`, which
  reads row `nv - 2`.
- **What happens if `OMP_NUM_THREADS` is more than the number of cores?** The
  threads take turns on the cores, and the barrier makes everyone wait for
  the slowest one, so it gets slower.
- **What does Amdahl's law say about this code?** The serial boundary work is
  about 0.5% of a step, so the cap from that alone is about 21.5 times on 24
  threads. The measured fit will be larger because it also catches the
  barrier and memory limits.
- **Where does Linux place the grid's memory pages?** Linux places each page
  in the memory of the socket whose thread first writes it. Our grid is
  filled by the main thread, so all of it is placed in socket 0's memory.

---

### 17. MPI and the hybrid rung (the `mpi` rung)

#### 17.1 Processes that share nothing

OpenMP threads share one process's memory. MPI (Message Passing Interface)
works with separate processes, usually called ranks. Each rank is a full copy
of the program with its own memory, and it may run on a different machine.
One rank cannot read another rank's variables. A rank moves data to another
rank by sending it as a message, which the other rank receives.

That is more work for the programmer, and it is needed to use more than one
machine, because two nodes have no shared memory at all.

`mpiexec -n 4 ./heston ...` starts four copies of the same program. Each one
calls the same four functions.

| Call | Does |
|---|---|
| `MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level)` | joins this process to the group started by `mpiexec`, and promises that only the main thread makes MPI calls |
| `MPI_Comm_rank(MPI_COMM_WORLD, &rank)` | which copy am I, 0 to n-1 |
| `MPI_Comm_size(MPI_COMM_WORLD, &num_ranks)` | how many copies there are |
| `MPI_Finalize()` | leaves the group, must be called before the program ends |

`MPI_COMM_WORLD` is the communicator holding every rank. Started without
`mpiexec`, the program is a group of one (a singleton) and `rank` is 0, which
is how the other rungs and `make test` run inside the same MPI binary.

#### 17.2 Splitting the grid into blocks of rows

Each rank owns a contiguous block of variance rows. With 512 rows on 2 ranks,
rank 0 owns rows 0 to 255 and rank 1 owns 256 to 511. When the rows do not
divide evenly the first ranks get one extra row. Rank 0 owns the v = 0 row and
the last rank owns the top row, so every boundary has exactly one owner. In
the code a block of rows is a `RowRange` called `owned_rows`, made by
`rows_for_rank`.

Rows are the right unit here for the same reasons as in OpenMP (§16.4), and
for one more. A row is contiguous in memory (`idx = j*Ns + i`), so a whole
row goes as one message straight out of the array. Splitting by columns
would need the edge column copied into a separate buffer every step, because
a column's values are `Ns` doubles apart.

#### 17.3 Copied rows

The stencil for row `j` reads rows `j - 1`, `j` and `j + 1` (§7.2). So a
rank's first row needs the row just below its block of rows, and its last row
needs the row just above. Those rows belong to the neighbours. Each step, each
rank sends its first and last rows to its neighbours and receives theirs into
its copied rows. The copied rows are one row deep because the stencil reaches
exactly one row up and one row down. The code names them by position,
`row_below` and `row_above`, and the function that sends the edge rows and
receives the copied rows is `exchange_edge_rows`.

```text
  rank 1 owns rows 256..511                    rank 0 owns rows 0..255

  row 511  top row (rank 1 owns it)
   ...
  row 256  first row    ─── sent to rank 0 ──►  row 256  row_above (copy)
  row 255  row_below (copy)   ◄── sent to rank 1 ───  row 255  last row
                                                 ...
                                                row 0  v = 0 row (rank 0 owns it)
```

Per step each rank sends at most two rows and receives at most two rows. At
2048 columns a row is 16 KiB, so a rank in the middle moves 64 KiB per step.
The compute per step is the whole block of rows, about half a million cells
for two ranks, so the message is small next to the work. That ratio is PLAN §5
hypothesis 6.

#### 17.4 Non-blocking messages and why the receives go first

`MPI_Send` and `MPI_Recv` are blocking. `MPI_Send` may not return until the
other side has started receiving. If two neighbours both call `MPI_Send`
first, each waits for the other to receive, and neither ever does. That is a
deadlock. Small messages often complete because MPI sends them with the
eager protocol and buffers them at the receiver, so the bug can stay hidden
until the rows grow past the eager limit.

The rung uses the non-blocking versions.

| Call | Does |
|---|---|
| `MPI_Irecv(buffer, count, MPI_DOUBLE, source, tag, comm, &request)` | starts a receive and returns immediately |
| `MPI_Isend(buffer, count, MPI_DOUBLE, dest, tag, comm, &request)` | starts a send and returns immediately |
| `MPI_Waitall(n, requests, MPI_STATUSES_IGNORE)` | blocks until all n sends and receives have completed |

Nothing blocks until `MPI_Waitall`, and by then every send has a matching
receive posted, so a deadlock is impossible. The receives are posted before
the sends so that when a row arrives, MPI already knows where it goes and can
copy it straight into the array instead of into a temporary buffer. The
rows being sent and received must not be read or written between the `I`
call and the wait. Here nothing reads or writes them, because the exchange
happens after the step's update and before the swap.

The `tag` is a number that labels a message. Both directions between a pair
of ranks use tag 0, which is safe because MPI matches a receive by source as
well as tag, and the two messages have different sources.

#### 17.5 Where the exchange sits in the step

```text
  each step on every rank
    1. update own interior rows      OpenMP threads, reads current, writes next
    2. own boundaries                v = 0 row on rank 0, top row on the last rank
    3. exchange edge rows of next    MPI_Irecv, MPI_Isend, MPI_Waitall
    4. swap_buffers                  next becomes current
```

The exchange is done on `next` at the end of the step, so when the following
step starts each rank already has fresh copied rows in `current`. At the
first step the copied rows are already correct without an exchange, because
every rank filled the whole grid with the payoff.

The wait in step 3 is also a synchronisation point. A rank cannot start step
k + 1 until its neighbours have finished step k, which keeps neighbouring
ranks in step with each other without any global barrier.

The time spent in step 3 is added up in `exchange_seconds` so the benchmarks can
show communication against compute. Experiment E9 overlaps it with compute,
by updating and sending the edge rows first and then doing the middle rows
while the messages are in transit.

#### 17.6 Hybrid, MPI between and OpenMP within

Each rank runs the OpenMP loop over its own rows, so the rung is MPI between
ranks and OpenMP threads inside each one. The usual layout is one rank per
node or one per socket, with a thread per core inside it.

| Layout on two r730 nodes | Ranks | Threads per rank |
|---|---|---|
| one rank per node | 2 | 24 |
| one rank per socket | 4 | 12 |
| one rank per core, pure MPI | 48 | 1 |

Hybrid needs fewer ranks than pure MPI, so there are fewer copied rows and fewer
messages. One rank per socket also keeps each rank's memory on its own
socket, which PLAN §5 hypothesis 5 says should beat one rank of 24 threads
whose pages are all in socket 0's memory (§16.8).

When launching, the rank must not be pinned to one core, or all its threads
would share that core. `--bind-to none` stops Open MPI pinning, or
`--map-by ppr:1:socket:PE=12 --bind-to core` gives each rank the 12 cores of
one socket.

#### 17.7 Every rank keeps the whole grid, a trade-off

Every rank allocates both full buffers and only updates its own rows. The
rest of its copy goes stale and is never read, apart from the two copied
rows.

That costs memory. At 2048 by 512 the two buffers are 16 MiB per rank, which
is nothing on a node with hundreds of GiB. In return, every row has the same
index on every rank, so `update_row_avx`, `apply_stock_edges` and the
boundary functions are exactly the code the other rungs run, and rank 0 can
reuse the same `extract_result` after one gather.

For a grid too big for one machine, each rank would instead allocate only its
block of rows plus one copied row on each side, and every row index would be
shifted by `owned_rows.begin - 1`. The maths would not change. That is the
textbook layout, and this rung would switch to it if memory ran out.

#### 17.8 Collecting the answer with MPI_Gatherv

After the last step, `MPI_Gatherv` copies each rank's rows into rank 0's spare
buffer. A plain `MPI_Gather` needs every rank to send the same amount, and
`Gatherv` (v for vector) takes a count and an offset per rank, which we need
because blocks of rows can differ by a row. Then rank 0 swaps the two buffer
pointers so the gathered buffer is `current`, and reads the price and Greeks
as every other rung does. The gather is outside the timed loop.

#### 17.9 Why it is bit-identical too

Each cell is computed by the same `update_row_avx` instructions as in the SIMD
rung. The copied rows hold exactly the neighbour's values, because a message
carries the 64-bit doubles unchanged. Nothing is summed across ranks. So the
MPI answer is the SIMD answer bit for bit, at any number of ranks and threads,
and `test_ladder_matches` checks this with zero tolerance.

With the edge row exchange removed on purpose, the two-rank answer moved by
2e-9. The stale copied rows (still the payoff) spread their error one row
further into the block each step, and the quote point is far enough away
that the error is small. A loose tolerance might have let that through, and
the zero tolerance caught it.

#### 17.10 How it was checked

| Check | Where | Result |
|---|---|---|
| `test_ladder_matches` with mpi as a singleton at 1, 2, 7 threads | Mac | mpi vs baseline 2.35e-15 call, 2.24e-15 put, mpi vs simd 0 |
| `make test-mpi`, 2 and 3 ranks, 2 threads each | Mac | mpi vs simd exactly 0 |
| demo grid price at 1, 2 and 4 ranks | Mac | 196.1054098 every time |
| edge row exchange removed on purpose | Mac | mpi vs simd 2.0e-9, test fails as it should |
| `test_ladder_matches` at 2, 3 and 4 ranks over two nodes, 2 threads each | rangpur job 651092, a100-a and a100-b | mpi vs simd exactly 0, mpi vs baseline 1.06e-15 call and 0 put |
| `make test` and `make test-mpi` on one node | rangpur job 651091, r730-0 | all pass, mpi vs simd exactly 0 as a singleton at 1, 2, 7, 8 threads and at 2 and 3 ranks |

#### 17.11 Self-test

- **What is the difference between MPI and OpenMP?** OpenMP threads share
  one process's memory. MPI ranks are separate processes with separate
  memory, possibly on different machines, and move data only by messages.
- **Why is one copied row on each side enough?** The stencil reads exactly one
  row above and one below, so a rank needs exactly one row from each
  neighbour.
- **Why are the blocks made of rows?** A row is contiguous in memory,
  so it goes as one message straight from the array. A column would have to
  be copied into a buffer first.
- **Why `Irecv` before `Isend`?** When the data arrives MPI already knows
  where to put it, so it is received directly into the array without an
  extra copy.
- **Why non-blocking at all?** Two neighbours doing blocking sends to each
  other at the same time can deadlock. With non-blocking calls nothing blocks
  until every send and receive has been posted.
- **What does `MPI_Waitall` protect?** The rows being sent and received. They
  must not be read or written until the transfers have finished.
- **Why does every rank hold the whole grid?** It keeps every row index the
  same on every rank, so the row functions and `extract_result` are reused
  unchanged. The cost is memory, 16 MiB per rank here. For a huge grid each
  rank would hold only its block of rows plus two copied rows.
- **Why `Gatherv`?** Blocks of rows can differ in size by a row, and
  `Gatherv` takes a count and an offset per rank, where `MPI_Gather` needs
  the same count from every rank.
- **Why is the answer bit-identical?** Every cell runs the same instructions,
  copied rows are exact copies, and nothing is summed across ranks.
- **What does the hybrid layout gain over pure MPI?** Fewer ranks means fewer
  exchange messages, and one rank per socket keeps each rank's memory on its
  own socket.
- **What does `--bind-to none` do?** It stops Open MPI pinning each rank to a
  single core, which would put all of that rank's threads on that one core.

---

### 18. Measuring the parallel rungs

#### 18.1 The protocol, and why each rule is there

Every number comes from an exclusive r730 node through `sbatch`. Each rule
below exists because breaking it would give a wrong or unrepeatable number.

| Rule | Why |
|---|---|
| exclusive node, pinned to r730-2 (and r730-3 for two nodes) | another job on the node would compete for cores, cache and memory bandwidth |
| build from clean, run every test, refuse to time if one fails | a fast wrong answer is worth nothing |
| `MALLOC_MMAP_THRESHOLD_=1048576` | M1 found the same binary runs 1.66 times faster on huge pages, and without this only the first rep got them |
| `OMP_PLACES=cores`, `OMP_PROC_BIND=close` | threads that migrate between cores lose the data in their cache, and the binding decides which socket's memory they use |
| median of 5 reps, min and max kept | one rep can be unlucky, the spread shows how much to trust a difference |
| every rep runs for several seconds | short runs are dominated by start-up and timer noise |
| flush denormals to zero with `LD_PRELOAD` (§18.5) | without it, short runs on big grids are full of very slow tiny numbers and different runs measure different regimes |
| every job builds in its own copy of the tree | two jobs sharing one directory can delete each other's files mid run, which happened once |
| `maturity_for` keeps `dt` at 70% of the stability limit | an unstable run fills the grid with huge numbers and infinities, which can change the speed |
| provenance header in every CSV | anyone can see which node, compiler, flags, binding and page regime produced a number |

The rate is reported as cell updates per second (`ns * nv * nt / seconds`),
so grids of different sizes and runs of different length compare fairly.

#### 18.2 Strong scaling and Amdahl's law

Strong scaling keeps the problem fixed (2048 by 512) and adds threads. The
speedup on `p` threads is `S(p) = T(1) / T(p)`, and the parallel efficiency
is `S(p) / p`, which is 1 when every thread is fully used.

Amdahl's law (§16.7) says a serial fraction `f` caps the speedup.

```text
  S(p) = 1 / ( f + (1 - f) / p )
```

The script fits `f` to the measured speedups. Rearranged, the law is a
straight line in `f`.

```text
  1/S - 1/p = f * (1 - 1/p)
```

So with `x = 1 - 1/p` and `y = 1/S - 1/p` for each thread count, the least
squares slope through the origin is `f = sum(x*y) / sum(x*x)`. No fitting
library is needed. The fitted `f` lumps together everything that does not
shrink with more threads, the serial boundaries, the fork and barrier each
step, and any shared limit such as memory bandwidth. It is a measured
number, and comparing it with the 0.5% boundary share says how much of the
loss comes from something other than the boundaries.

The figure draws the AVX2 omp rung and the scalar omp rung. Each line is
normalised to its own one-thread run, so the two lines show how well each
scales, and the gap between their absolute rates shows SIMD and threads
multiplying.

#### 18.3 Weak scaling and Gustafson's law

Weak scaling keeps the work per thread fixed and grows the problem with the
threads. Each thread gets 21 interior rows, so the grid has `21p + 2` rows on
`p` threads. If scaling were perfect, the time per step would stay flat.

The variance axis is stretched with the rows (`--vmax`), so the spacing
between rows stays the same, and every row is the same kind of work as in the
strong runs. The stretched axis makes the stability limit smaller, so the
script shortens the maturity to match. That changes nothing about the speed.

Gustafson's law looks at the same fraction from the other side. If a run on
`p` threads spends a fraction `f` of its time in serial work, the same job
on one thread would take `f + p(1 - f)` as long, so the scaled speedup is

```text
  S_scaled(p) = p - f (p - 1)
```

which grows with `p` without a cap. The point is that a bigger machine is
usually used for a bigger problem, and then the serial part matters less.

The two laws use different fractions. Amdahl's `f` is the serial share of
the one-thread run, and Gustafson's `f` is the serial share of the
`p`-thread run. So the figure does not reuse the strong-scaling fit. It
takes Gustafson's `f` straight from the 24-thread weak run, as
`boundary_seconds / seconds`, and plots the measured scaled speedup
`p * t(1) / t(p)` against that line.

One effect to know before reading it. On 1 thread the grid has only 23 rows,
about 370 KiB per buffer, which is larger than the 256 KiB L2 of one core.
On 24 threads it is 506 rows, 8 MiB per buffer. Both fit in L3, so both
ends of the sweep work from L3. The one thread has a whole socket's L3 and
memory bandwidth to itself, while the 24 threads share them, so part of the
rise in step time comes from that sharing and 9.2 percent from the serial
boundaries (RESULTS.md M2.3).

#### 18.4 The roofline, with more than one roof

A roofline plots achieved GFLOP/s against arithmetic intensity, the flops
done per byte moved to or from memory. A kernel cannot go above either roof.
The flat roof is the compute peak, and the sloped roof is bandwidth times
intensity.

The kernel does 31 flops per cell update and moves 24 bytes from DRAM in the
best case (read the cell once, write it once, and the write-allocate read),
so its intensity is 31/24 = 1.29 flops per byte, the same as in M1. FMA does
not change the count, because one FMA is still a multiply and an add.

M2 adds roofs that `bench_peak.sh` measures on the node itself.

| Roof | What it measures |
|---|---|
| scalar peak, 1 core | eight independent scalar multiply-add chains, the M1 roof |
| AVX2 FMA peak, 1 core | ten independent `_mm256_fmadd_pd` chains, 4 lanes, 2 flops each |
| AVX2 FMA peak, 24 cores | the same on every core simultaneously |
| triad bandwidth, 1 thread and 24 threads | `a[i] = b[i] + 3 * c[i]` over three 256 MiB arrays, counted as 24 bytes per element |

Each rung's point sits on the same vertical line, because the maths is the
same, and moves up as it gets faster. The question the figure answers is
which roof stops it. SIMD lifts the compute roof by four lanes and FMA,
threads lift it by 24 cores, but the bandwidth roof only grows with the
number of memory channels in use. A point above the DRAM roof is possible
when the whole grid fits in cache, because then most traffic never reaches
DRAM.

#### 18.5 Denormals, the trap the first round fell into

A double stores a number as a fraction times a power of two. The smallest
power it can hold for a normal number is 2^-1022, about 2.2e-308. Below
that, the format switches to denormal (also called subnormal) numbers, which
give up precision to reach down to about 4.9e-324 before reaching zero.

The hardware does not have a fast path for them. On this Xeon, an add or a
multiply that reads or produces a denormal is handed to a microcode routine
that takes on the order of a hundred cycles instead of a few.

Our grid makes them naturally. A call is worth exactly zero below the
strike at expiry, and each step the diffusion spreads a little value one
cell further into that zero region (§7.6). A cell far from the strike
therefore holds values like 1e-200, then 1e-310, then 1e-320, then
exactly 0. How many cells hold denormals depends on how many steps have run
and on the grid, and benchmark runs vary both.

The first size sweep showed it. One grid ran 2.4 times slower than the grids
either side of it on every solver, and the prices printed were as small as
1e-43. A direct test then showed even the standard grid lost 17% to
denormals at 381 steps.

The fix is the same idea as pinning huge pages. Every timed run sets two bits
in the MXCSR register, the control register for SSE and AVX arithmetic.

| Bit | Name | Effect |
|---|---|---|
| 15 | FTZ, flush to zero | a result too small to be a normal double is written as 0 |
| 6 | DAZ, denormals are zero | a denormal input is read as 0 |

They are set by a tiny library loaded with `LD_PRELOAD`, so the program and
the tests are not changed. The values it throws away are below 1e-308 in a
grid whose prices are in the hundreds, so the printed price did not change in
its first ten digits. `-ffast-math` turns on the same bits, which is one
reason it can make code faster, and one reason it is not allowed to slip in
unnoticed. Experiment E12 measures what denormals cost each rung.

#### 18.6 What the measurements showed, and why

The full tables with job ids are in `docs/RESULTS.md`, Milestone 2. These are
the results to be able to explain, each with its reason.

**The ladder.** baseline 57, serial 159, simd 456, omp on 24 threads 6367,
mpi as 2 ranks of 12 on one node 6984 Mcell-updates/s. That is 123 times the
baseline on one node, and every rung gives the same answer.

**SIMD gave 2.9 times, and four lanes set the upper limit at 4.** The loop
issues nine overlapping loads per four cells, and the unaligned ones cross a
cache line on every other iteration. The last cells of each row are computed
one at a time. At 456 Mcell-updates/s one core already moves about 11 GB/s,
close to the 12.4 GB/s it can stream from memory. Before the build kept the
256 bit loads unsplit, the rung gave only 2.4 times, because g++ split every
unaligned 256 bit load into two 128 bit loads (E2 and E2b).

**Threads gave 13.9 times on 24 cores.** Three things account for the rest.
The serial boundaries grow from 0.44% of a one-thread step to 9.1% of a
24-thread step, because the main thread has to fetch rows that other cores
have just written. The fork and barrier cost the same every step however short
the step gets. And at 24 threads the AVX2 kernel moves about 151 GB/s through
the caches, so it runs out of shared bandwidth first. That is why the scalar
row, whose threads each need data at a third of the rate, scales better, 18.5
times.

**Amdahl's fit gives f = 2.4%.** That is more than five times the boundary
share of a one-thread run, so most of the loss is the per-step overhead and
the bandwidth that grow with the thread count. The scalar row fits 1.6%. A
faster row makes the fixed costs a larger share of each step, which is also
why f doubled from 1.2% when the 256 bit loads were unsplit.

**Small grids scale badly, big grids hit memory.** 512 by 128 on 24 threads
is only 8 times serial, because each step is too short to amortise the fork
and barrier. 2048 by 512 fits in the two L3 caches and gets 40 times. 4096 by
1024 does not fit, sits on the measured DRAM roof, and drops to 19 times.

**Page placement matters.** The grid is written first by the main thread,
so Linux places every page in socket 0's memory. On the big grid, spreading
the pages over both sockets with `numactl --interleave=all` gave 1.74 times
in the latest run and 1.28 times in the two before it, because the runs
without it come out near either 2800 or 4000. One MPI rank per socket, each
with its own copy of the grid in its own socket's memory, gave 2.28 times
one 24-thread rank.

**Two nodes gave 1.21 times on the small grid and 1.64 on the big one.** The
exchange is only 64 KiB per rank per step, but a step on the small grid is
only 0.126 ms, and each exchange takes about 50 microseconds, mostly latency
and waiting for the neighbour. The bigger grid has longer steps, so the same
exchange is a smaller share, and it reaches 82% efficiency. The small grid's
two-node layouts also vary between jobs, so its gain lies between 1.21 and
1.33 times.

**Three things outside the code moved the numbers more than some
optimisations.** Huge pages (1.6 to 1.8 times, E11), denormals (1.07 to 1.42
times, E12), and page placement (up to 2.28 times). That is why the
protocol pins all three, and why each is reported as its own experiment.

**What was tried and rejected.** Dynamic and guided scheduling lost 22 to
43%, because every row costs the same and handing rows out on demand only
adds bookkeeping. One parallel region per solve (E4) and parallel boundaries
(E8) were small wins, 2% and 6% on the benchmark grid and up to 13% on the
small grid, and were not folded into the ladder so the ladder keeps one
technique per rung.

#### 18.7 Profiling and optimisation levels

**What a flat profile is.** Building with `-pg` makes the program note
which function it is in every 10 milliseconds and write the tally to
`gmon.out` when it exits. `gprof -b -p ./heston gmon.out` prints the flat
profile, one line per function with its share of the run time, biggest
first. It is called flat because it does not say who called whom, only
where the time was spent. It answers the first question of any
optimisation, which is where the time goes, before anything is changed.

**What it showed** (job 652697, one thread, 2048 by 512). The baseline
spends 99.9% in `BaselineSolver::solve`, because the whole algorithm is
that one loop nest. The serial rung spends 99.4% in `update_cells_scalar`
and the simd rung 98.9% in `update_row_avx`. The boundaries are under half
a percent. So the interior update is the hot spot, and it is the only
thing worth making faster or splitting between cores. A profile only
works at the level of functions, which is one reason the later rungs put
each piece of work in its own function.

**What the -O levels do.** `-O0` does no optimisation, so every variable is
kept in memory, loaded before each use and stored after each change. `-O1`
keeps values in registers and removes repeated work. `-O2` adds instruction
scheduling and more, and on GCC 8.5 it leaves the auto-vectoriser off. `-O3`
turns on the auto-vectoriser, which rewrites a loop to work on several cells
per instruction when it can prove that is safe.

```text
level   baseline   serial      Mcell-updates/s, job 652696
-O0       10.7      44.8
-O1       57.6     156.1
-O2       57.7     158.4
-O3      112.5     238.6
```

Every level gave exactly the `-O2` answer, because none of them is allowed
to reorder floating point arithmetic. `-O3` gains 1.95 times on the
baseline and 1.51 times on the serial rung, using two-lane SSE2, since the
build has no `-mavx2`. The hand-written AVX2 rung reaches 456, 1.91 times
the `-O3` serial rung, because it uses four lanes and FMA, which the
compiler is not given. The project builds at `-O2` so that each rung of
the ladder measures one change of ours, with the compiler's optimisations
held fixed. Milestone 1's sweeps F and G (§13, RESULTS.md §6) found the same
`-O3` effect on the baseline.

#### 18.8 Self-test

- **Strong or weak scaling, which keeps the problem fixed?** Strong. Weak
  keeps the work per thread fixed and grows the problem.
- **How is Amdahl's serial fraction fitted without a library?** `1/S - 1/p`
  is `f` times `1 - 1/p`, a straight line through the origin, so `f` is the
  least squares slope `sum(x*y) / sum(x*x)`.
- **Why are Amdahl's and Gustafson's f not the same number?** Amdahl's is the
  serial share of the one-thread run, Gustafson's the serial share of the
  many-thread run.
- **Why does the one-thread weak scaling point look good?** Its 23 row grid
  fits in one core's L2 cache.
- **What is the arithmetic intensity of the kernel and does FMA change it?**
  31 flops over 24 bytes, 1.29. FMA does not change it, because one FMA is
  still two flops.
- **How can a point sit above the DRAM roof?** When the grid fits in cache,
  most of its traffic never reaches DRAM.
- **What is a denormal and why did it matter?** A denormal is a double smaller
  than about 2.2e-308. Arithmetic on one is far slower, and short runs on big
  grids were full of them, so different runs were measured in different
  regimes.
- **What do FTZ and DAZ do?** FTZ writes tiny results as zero, DAZ reads tiny
  inputs as zero, so the slow path is never taken.
- **Why did SIMD give 2.9 times with four lanes?** The row issues nine
  overlapping loads per four cells, the unaligned ones cross a cache line on
  every other iteration, the scalar remainder runs one cell at a time, and
  one core already moves about 11 GB/s. Before the 256 bit loads were
  unsplit it was 2.4 times.
- **Why 13.9 times on 24 threads?** The serial boundaries grow to 9% of a
  step, the fork and barrier cost the same every step, and the AVX2 kernel
  runs out of shared cache and memory bandwidth.
- **Why does one rank per socket beat one rank of 24 threads?** Each rank's
  main thread writes its own copy of the grid first, so Linux places those
  pages in that rank's own socket's memory. One 24-thread process has every
  page in socket 0's memory.
- **Why do two nodes not give twice the speed?** Each step ends with an edge row
  exchange that costs about 50 microseconds of latency and waiting, against
  a step of only 0.126 ms on the small grid.
- **Which three settings outside the code mattered most?** Huge pages,
  flushing denormals, and where memory pages are placed.
- **What does a flat profile of the baseline show, and why is it not more
  useful?** 99.9% in one function, because the baseline is one loop nest.
  The later rungs split the work into functions, and their profiles show
  99% in the row update and under half a percent in the boundaries.
- **Why build at -O2 when -O3 is faster?** `-O3` auto-vectorises with two
  SSE2 lanes and gives 1.5 times on the serial rung. Building every rung at
  `-O2` means each rung's gain comes from one change of ours. The hand AVX2
  rung at `-O2` is still 1.9 times the `-O3` serial rung.

---

## Part V. Rapid-fire interview drill

Answer these out loud. If any of them trip you up, reread its section.

**The P2 checkpoint drill (2026-08-09).** All five of these tripped you up the
first time. Re-explain them at the start of the next session, then re-quiz
until they stick.

- **Why does `init_payoff` write into `current_`?** The time loop always
  READS `current()` and WRITES `next()`. The payoff is the first grid the
  loop ever reads, so it must start in the read buffer. (§11)
- **With one buffer instead of two, what goes wrong in the middle of a
  timestep?** Updating in place means some of the neighbours you read were
  already updated *this* step, so the stencil mixes time step n with n+1. The
  results then depend on the loop order and are wrong. With two buffers, the
  buffer being read does not change during the step. (§11)
- **Why use a FORWARD difference for V_v on the v = 0 row?** A central
  difference needs a node at v = −dv, which doesn't exist and
  can't, since negative variance is meaningless. Also, the drift kappa·theta
  > 0 points INTO the grid, so taking the derivative from above, which is
  where the information comes from, is the stable and physical choice.
  Imposing a fixed value would over-determine the reduced equation. (§9)
- **When dt crosses the stability bound, which cells blow up first, and
  why?** The high-S, high-v corner blows up first. The diffusion coefficients
  (½vS² and ½ξ²v) are largest there, so that is where the centre weight goes
  negative first, and the error flips sign from cell to cell each step, which
  is the checkerboard. (§7.5 and §8)
- **Feller is violated for our parameters, so why not fix it, and what does it
  actually cost us?** 2κθ ≥ ξ² is a property of the market-calibrated
  parameters, which our code cannot change, and equity calibrations usually
  violate it. The consequence is that the variance process puts real
  probability near v = 0, where our uniform v grid is relatively coarse, so
  that region is slightly under-resolved. The scheme stays valid, and we state
  it as a known limitation, with a non-uniform v grid as future work. (§9)

**Finance and the model**
1. What's the difference between a call and a put? What about European and
   American?
2. Why does Heston exist when Black-Scholes has a formula? (Volatility moves.)
3. What does each Heston parameter do? Why is rho negative?
4. What are delta, gamma and vega, and who cares about them? (Hedgers.)

**Numerics**
5. Why solve backwards from expiry? (That's where the answer is known.)
6. Why does the stencil need 9 points where the heat equation needs 5? (The
   cross term, from the correlation.)
7. What makes the explicit scheme blow up, where does it start, and why
   there?
8. Explain the v = 0 boundary to a child. Then explain the Feller condition.
9. How do you *know* the price is right without market data? (Three tests.)
10. What does the convergence plot prove that the other tests don't?

**Code and HPC**
11. Who owns the grid memory, and when is it freed? (RAII, at the end of the
    scope.)
12. Why use two buffers and a swap instead of updating in place?
13. Why is a flat array with a contiguous inner loop faster than a list of
    lists?
14. What was the first serial optimisation, and why did it win? (Hoisting the
    repeated per-cell weight computation, together with replacing divisions,
    which is where the speedup came from.)
15. Why benchmark on a compute node, and why take the median of several runs?
16. Where will the parallelism come from in M2, and what limits it? (The
    per-step barrier, which brings in Amdahl's law.)

**Milestone 2 build decisions**
17. Why does `update_cells_scalar` copy every constant into a local variable
    before the loop? (The loop writes through a `double*` into the next
    buffer. The compiler cannot prove that pointer never points into the
    `StencilConstants` struct, so it would have to reload each constant from
    memory after every store. A local variable cannot be reached through
    that pointer, so it stays in a register. This is called pointer
    aliasing.)
18. Why is every rung, even `baseline` and `serial`, built with `mpicxx`?
    (There is one binary for the whole ladder, and the `mpi` rung needs the
    MPI library. `mpicxx` is g++ with the MPI include and library paths added,
    so the serial rungs compile exactly as with g++. `main` calls
    `MPI_Init_thread` once, and run without `mpiexec` it becomes a single
    process MPI program.)
19. Why is each experiment its own solver class, with no flag inside
    `OmpSolver`? (Each one is the `omp` loop with one thing changed, so the
    rung the video shows stays one pragma, and each experiment goes through
    `test_ladder_matches` like a rung before it is timed. `omp_runtime`
    covers every schedule with one class by reading `OMP_SCHEDULE`.)

**The grid, the stencil and the buffers (work these on paper)**
- In the first timestep, which time level do the nine values come from? (The
  expiry grid, §6.3.)
- Is the nine-point update an average of the nine cells? (No. It is a
  weighted sum with very uneven weights that add up to `1 − r·dt`, §6.4.)
- Which of the nine cells does each derivative read, and why is the cross
  term the only one that uses the corners? (§7.2)
- Work out the strike cell's value after one step on the test grid. (§7.4,
  the answer is 0.04959.)
- Why must the centre weight stay positive, and how does that give the
  stability limit the solver prints? (§7.5)
- Why does today's spot cell stay at zero for the first step? (§7.6)
- Give a concrete example of what goes wrong with one buffer. (§11)
- Why do the `current` and `next` references in the solver still point at
  the right buffers after `swap_buffers`? (CODE_GUIDE.md §7.2)

CODE_GUIDE.md §14 has the practice questions on the code itself.

**Limitations to have ready**
17. The uniform grid under-resolves the region near v = 0 (Feller is
    violated). With more time we would use a non-uniform grid.
18. The explicit scheme's dt ∝ ds² limit makes fine grids expensive. Implicit
    or ADI schemes remove it, at the cost of solving linear systems.
19. Cutting the domain off at S_max and v_max adds a small, measured boundary
    error. Doubling the domain doesn't move the price.
