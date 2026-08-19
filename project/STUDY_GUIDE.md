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

### 5a. Where is the strike? (the model and the contract are separate)

Look at the two random equations, or at the PDE above, and the strike `K` is
nowhere in them. If you are used to the Black-Scholes formula, where `S` and
`K` sit side by side, this is disorienting. The resolution is that pricing
always has two separate ingredients:

1. **The model** — how the share price behaves. That is what the two random
   equations (and the PDE they become) describe. They apply to *any*
   contract on this share and know nothing about yours.
2. **The contract** — what you are owed at expiry. That is the payoff
   `max(S−K, 0)` for a call, and it is the *terminal condition*: the values
   the grid is filled with at `t = T` before the solver steps backwards
   (section 9's first row).

So the strike is in the program, but it lives in `init_payoff` and the
boundary values, not in the update stencil. The Black-Scholes *formula*
feels different only because it is the special case where both ingredients
could be merged into one closed expression. Heston does not build on that
formula; it relaxes one of its assumptions (constant volatility) and
contains Black-Scholes as the `xi = 0` special case, which is exactly what
the collapse test exploits.

Quiz yourself: if the same model prices a put instead of a call, what
changes in the program? Only the payoff and the `S`-boundary values. The
PDE, the stencil, and every optimisation stay identical.

### 5b. Where did the randomness go? (the question every marker will ask)

Search the C++ for a random number generator and you will find **nothing** —
no `<random>`, no `rand()`, no `mt19937`, not one seed. The model is defined
by two *random* equations, and the code that solves it is completely
deterministic. That is not a bug. It is the whole point of the PDE route,
and you need to be able to say why in one breath.

**The one-sentence answer:** we never simulate the randomness because we
solved for its *average* in advance — the PDE is the equation that the
expectation over all random paths satisfies.

**The theorem that licenses it (Feynman–Kac).** The fair price is by
definition an expectation over random futures:

```text
V(S, v, t) = E^Q[ e^(−r(T−t)) · payoff(S_T)  |  S_t = S, v_t = v ]
```

Feynman–Kac says that this expectation, viewed as a function of today's
state `(S, v, t)`, satisfies exactly the PDE in §5. So there are two honest
ways to get the number:

| route | what it does | cost |
|---|---|---|
| Monte Carlo | draw millions of random paths, average the payoffs | error shrinks like 1/√N — slow |
| **PDE (ours)** | solve the equation the average obeys | error shrinks like the grid spacing — fast, and gives every `(S, v)` at once |

**The intuition behind the squaring (the version the video uses).** Over one
small step the noise shocks `S` up or down by an amount of size `sqrt(v)·S`,
with both directions equally likely. Average the value over the two cases.
The slope contributions cancel, because what the up shock gains the down
shock loses, and the first surviving term of the Taylor expansion is
`½·(shock size)²·V_SS`. Squaring the shock size `sqrt(v)·S` is exactly what
produces `½·v·S²·V_SS`. This is the plain-English shadow of the "quadratic
variation" entries in the table below: `dW` itself never rewrites into the
formula, it gets *averaged out*, and the curvature term is what the
averaging leaves behind.

**Where each piece of noise ended up.** This is the part worth memorising,
because you can point at the code and name the term:

| in the SDE | in the PDE (`solver_baseline.cpp:84–85`) | why |
|---|---|---|
| `sqrt(v)·S·dW1` | `0.5 · v · S² · V_SS` | quadratic variation of S: the diffusion coefficient gets **squared**, so `sqrt(v)² = v` |
| `xi·sqrt(v)·dW2` | `0.5 · xi² · v · V_vv` | same, for variance |
| `corr(dW1, dW2) = rho` | `rho · xi · v · S · V_Sv` | cross-variation `d⟨S,v⟩ = rho·(sqrt(v)S)·(xi·sqrt(v)) dt` |
| drift of S under Q | `(r − q) · S · V_S` | risk-neutral: the real-world drift `mu` is irrelevant |
| drift of v | `kappa · (theta − v) · V_v` | mean reversion |

Notice the tell: **`sqrt(v)` appears nowhere in the numerics — only `v`.**
That squaring *is* the fingerprint of the noise having been integrated out.
If you ever see `sqrt(v)` in a diffusion coefficient in a PDE, something is
wrong.

Two details worth knowing before someone asks:

- **The cross-term has no ½.** Itô's lemma contributes
  `½·(V_Sv + V_vS)·d⟨S,v⟩ = V_Sv·d⟨S,v⟩` — the mixed partial appears twice,
  so the ½ cancels. A spurious `0.5` there is the single most common Heston
  implementation bug, and ours does not have it.
- **`kappa` and `theta` are already risk-neutral.** The variance drift
  carries no market-price-of-volatility-risk term `lambda`. That is standard
  and correct *for calibrated inputs*, but it is a silent assumption: feed in
  `kappa`/`theta` estimated from historical (real-world) data and you would
  be mixing measures.

**How we proved all this rather than asserting it** (`scripts/monte_carlo_check.py`,
run 2026-08-09 on `config/demo.cfg`):

```text
PDE (deterministic grid)          196.1054
Monte Carlo (400k random paths)   196.1345  ± 0.5524   <- actually simulates dW1, dW2
semi-analytic (Fourier)           196.1692             <- no randomness, no mesh
PDE delta 0.566967   vs   MC delta 0.566829
sample corr(z1, z2) = −0.6991     (target rho = −0.7000)
```

Three completely different machines, one number. The Monte Carlo *does*
carry a random number generator, correlates its two Brownian increments with
the Cholesky factor `z2 = rho·z1 + sqrt(1−rho²)·z_indep`, and lands on the
grid's answer. That is the empirical proof that the deterministic grid really
is computing an average over random paths.

**Self-test:** Where is `dW1` in the code? Why is the coefficient `v` and not
`sqrt(v)`? Why does the cross-term have no ½? What would break if you fed in
real-world `kappa` and `theta`?

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

**P6 rendering notes (2026-08-09) — what the four clips show and why they
look the way they do:**

- `weather_map.mp4` — raw V(S, v), inferno heat ramp on a DARK stage. The
  first cut used blue-on-white and was unreadable: the low-value half of a
  light-to-dark ramp vanishes into a white figure background. Dark stage +
  black→red→yellow means cold recedes and heat glows.
- `weather_map_timevalue.mp4` — **the actual heat-transfer shot.** Raw V
  spends the whole colour range on the boring linear deep-ITM ramp (up to
  ~$15,750) while the PDE's real work happens in a ~$400 band near the
  strike. Subtracting the expiry payoff (time value = V − payoff) cancels
  the ramp: what remains is a flame anchored at K = 5250, widening as
  variance grows, pinching to a point at v = 0 — where diffusion dies (all
  diffusion terms carry a factor v) and only transport survives. Time value
  IS what diffuses; this clip is the spreadsheet metaphor made visible.
- `surface.mp4` / `surface_timevalue.mp4` — the same data as a 3-D sheet
  with a slow camera drift: the kinked payoff ramp smoothing out, and the
  time-value sail growing out of the flat expiry sheet. The "matrix being
  computed", frame by frame.
- Mechanics: one fixed colour scale across all frames (else the animation
  "breathes"), 98th-percentile cap so late-time high-v glow doesn't
  compress early structure, NaN cells magenta on dark / black on white so
  broken cells are unmistakable.
- **Marker-friendly framing (author feedback):** headlines count calendar
  time ("46 of 91 days before expiry"), with the raw solver step demoted
  to a smaller second line; every flat frame carries a dashed line at the
  strike, a ring on today's market point (S = 5200, v = 0.04) labelled
  "the price we quote", and a one-line caption saying what the sheet is
  and which direction the solver moves. Smooth run dumps every 100 steps
  (200 frames) so the motion reads as flow, not slideshow. One trap fixed
  on the way: the step total must be read from the FULL dump list before
  `--max-frames` truncation, or every title mislabels the fraction.

**P6b — the clips were measurably motionless, and what fixed it
(2026-08-09).** Author feedback was "it's hard to see anything change".
That turned out to be literally true, not a matter of taste. Measured on
`results/smooth` in the default render (raw V, linear scale):

- The colour bar spans $15,337 (the deep-ITM ramp), of which the entire
  time evolution occupies **5.66%**.
- **Zero** of the 199 frame-to-frame transitions contained a single pixel
  that moved by even one 8-bit colour level. Largest per-frame change:
  0.51 of one level; median frame: 0.05.
- 83.6% of the frame area never changed perceptibly (CIELAB ΔE < 2).

Three scale mismatches stack up, and each has a fix:

1. **Wrong colour range.** The payoff ramp eats it all. → `--time-value`
   (already existed) plus `--gamma 0.5` (a power-law scale, because heat
   grows like `sqrt(tau)` so a linear ramp leaves the first half of the
   clip black) and `--cap 99.9` (the old 98th-percentile cap *clipped* the
   final frame flat).
2. **Wrong viewport.** With `s_max = 4K` and `v_max = 1.0`, the payoff
   smoothing band is ~5% of the S axis and the reachable variance band is
   12.8% of the v axis — `v_max = 1.0` sits **32.7 standard deviations**
   above where the variance process actually goes. → `--zoom`, which crops
   to ±3 CIR standard deviations of `v` and the matching lognormal spread
   in `S`. Not a guess: it uses the exact CIR variance formula.
3. **Wrong time sampling.** Diffusion widths grow like `sqrt(tau)`, so
   uniform frames front-load all the motion — 75% of a uniform clip lives
   in the second half of the `sqrt(tau)` progression. → `--warp sqrt`, and
   fewer frames overall (the largest single win: at 200 frames the median
   frame's visibly-moving area is 0.00%; at 128 it is 10.5%).

Two additions that make the motion unmistakable rather than merely visible:
`--envelope` draws the ±1σ cone `S = K·exp(±sqrt(v·tau))` — the set of
scenarios whose outcome is still in doubt — which visibly flares open every
frame; and `--price-trace` puts a live dollar readout and an accumulating
curve of the quoted price in the corner, so numbers move even where pixels
don't.

Result on `results/demo`: **~10% of pixels change visibly every frame,
up from 0.00%.**

**The deeper fix was the grid, not the renderer — `config/demo.cfg`.**
Cropping in the renderer exposed that the smoke grid only had **9 variance
rows** inside the reachable band, so the "zoomed" frame was blocky. The
honest fix is to stop spending nodes on variance levels the market cannot
reach. Measured (holding `dv` fixed at 0.002, `ns = 421`): `v_max = 0.4`
gives 195.8869, `v_max = 0.2` gives 195.8876 — a **0.0004%** difference.
Two-thirds of the frame, and two-thirds of the flops, were being spent on
fantasy. Combined with node alignment (spot on stock node 130, `v0` on
variance node 20), `demo.cfg` reaches **196.105 in 2.84 s**, against the
P4 Richardson estimate of ~196.16 — better than the P4 ladder's
421×101×56000 rung, which needed 4.92 s to reach 195.925. Putting the nodes
where the physics is beat spending more of them. This is a genuine HPC
point, not a plotting one, and it belongs in the video.

**P6c — the two "Heston becomes Black-Scholes" animations (2026-08-09).**
Both answer the same question visually, and they are complementary:

- `results/bs_collapse.mp4` — `scripts/bs_collapse_anim.py`. Set `xi = 0`,
  `rho = 0`, `v0 = theta`, and the Heston PDE reduces *exactly* to
  Black-Scholes: `xi = 0` kills both `0.5·xi²·v·V_vv` and the cross-term,
  and `v = theta` makes `kappa(theta − v)` vanish, so the `v = v0` row
  decouples completely and runs the pure Black-Scholes stencil. The clip
  animates that row against the closed-form curve at the same `tau`, with a
  dollar-error panel underneath. Measured: final-frame max error **$0.0375**
  (1.85 basis points of the at-the-money price); halving `ds` cuts it to
  $0.0099 (a 3.78× ratio — clean second order), and the refined quoted price
  is **202.808075 vs closed-form 202.8079942, a relative error of 4.0e−07**.
  This is the strongest single validation number the project has.
  (`v0` must land exactly on a variance node or the collapse is only
  approximate — the script hard-fails rather than animate the wrong row.)
- `results/xi_smile.mp4` — `scripts/xi_sweep_smile.py`. 30 solver runs
  sweeping `xi` from 0.35 down to 0; each final sheet is converted into
  **Black-Scholes implied volatility per strike**. At `xi = 0.35` you get a
  skewed smile spanning **12.58 vol points**; at `xi = 0` it collapses to a
  flat line at 20.0% (residual spread **0.076 vol points**, which is the
  numerical noise floor of the whole pipeline, since the exact answer is
  flat 20%). This is the picture of the sentence "Heston contains
  Black-Scholes as a special case". Strikes come from the homogeneity of a
  European call (rescaling `K` maps the grid onto itself because `ds` scales
  with `K`) — verified against re-solving 5 strikes from scratch, worst
  disagreement **3.9e−10**, i.e. the solver's print precision.

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

**P5 measured boundary (2026-08-09, smoke grid 256×64, T = 0.25):** the
solver prints dt_stable ≈ 1.526e-5. Bracketing empirically by sweeping nt:
nt = 15000 (dt = 1.092× the bound) still finishes finite (price 205.3226 —
note the 4th decimal already drifting: a marginal mode decaying slowly);
nt = 14900 (1.100×) ends at 2.8e+69 — caught mid-explosion; nt = 14600
(1.122×) is NaN. So the *real* cliff sits ~9% above the printed estimate.
Why conservative: the printed bound plugs in the single worst cell's
coefficients (the far corner), but an unstable mode occupies a *region*,
so its growth rate averages in neighbouring, smaller coefficients. A safety
bound that errs safe is exactly what you want — say this in the video.
The horror-show config (`unstable.cfg`) runs at 4× the bound: checkerboard
erupts at the high-S/high-v corner and overflows to NaN by ~step 450 of
4096; `weather_map.py --diverging` films it (red/blue = the +/- oscillation,
black = overflowed cells).

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
4. **Fourier reference — DONE (2026-08-09).** Heston has a semi-analytical
   price via its characteristic function. `scripts/monte_carlo_check.py`
   carries one (numpy only, Albrecher "little Heston trap" form) *and* a
   Monte Carlo that simulates the SDEs directly, so any disagreement can be
   arbitrated by a third opinion instead of argued about.

The iron rule: **no version is benchmarked until it matches the reference.**
A fast wrong answer is worth nothing.

**The blind spot tests 1–3 share, and why test 4 mattered.** Tests 1 and 2
compare the solver *against itself* at the same grid cell, and
`extract_result` used to read the price at the cell **nearest** (spot, v0)
with no interpolation. A shared quote-point offset therefore cancels exactly
in put-call parity, and appears identically on both sides of the BS-collapse
comparison — so both tests pass while the absolute price is wrong. That is
exactly what happened on `config/smoke.cfg`: `dv = 1/63 = 0.015873`, so
`v0 = 0.04` snapped to node 3 = **0.047619**, a 19% jump in variance, worth
about +$16 at the reported vega; the stock snap was worth about −$6.6. Net
+$9.2 — and smoke's quoted 205.32 was indeed **+4.67%** above the true
196.1692. It was not instability and not discretisation error: the same run
matched the semi-analytic price *evaluated at its own snapped node* to
within $0.29. Node-aligned configs (the P4 ladder, `config/demo.cfg`) never
had this problem — on `demo.cfg` the three routes agree to 0.03%. The
lesson worth carrying into Milestone 2: **a validation suite made only of
self-consistency checks cannot see a systematic offset. One comparison
against an external absolute reference can.**

**The fix (P7): bilinear interpolation at the quote point.** `extract_result`
no longer snaps. `Grid` gained `stock_cell()`/`stock_weight()` and the
variance pair: the cell just below the quote point, plus how far across it
the point sits. The price is a blend of the four surrounding cells, and each
arm of the Greek stencils (east/west/north/south) is blended the *same* way,
so the reported delta really is the derivative of the reported price rather
than a derivative of some other surface.

Three details worth being able to defend on camera:

- **Node alignment still had to keep working, bit-for-bit.** The whole P3/P4
  record — `demo.cfg`'s 196.1054098, the four convergence rungs, both test
  grids — is quoted from aligned grids. On an aligned grid the weight is
  exactly 0, and the blend is written so that `weight = 0` collapses to
  `1.0*V + 0.0 + 0.0 + 0.0`, which returns the nodal double unchanged. That
  was *checked*, not assumed: every one of those numbers is still identical
  to the last digit after the change. Anything else would have meant the
  interpolation was wrong.
- **The weight is snapped to 0 within 1e-9 of a node.** `0.04 / (0.2/100)`
  happens to come out as exactly 20.0 in IEEE double on these configs, but
  relying on that would be relying on luck; a grid that landed on
  19.999999999999996 would silently interpolate between nodes 19 and 20 and
  break the bit-identity guarantee. The snap makes the guarantee structural.
- **The cell index is clamped to `[1, n-3]`, not `[1, n-2]`.** The blend
  reads `cell` and `cell+1`; the gamma stencil reads one further out on each
  side. So `cell-1` and `cell+2` must both exist, which costs one more node
  of margin at each edge than the old nearest-node readout needed.

The solver now also **says when a config is off-node**: one stderr line
naming the offset in cells, with the same `[baseline]`/`[opt-L6]` prefix as
the Feller line. Node alignment stays the discipline — the warning exists so
that breaking it is a visible choice rather than a silent one.

Measured effect: `smoke.cfg` went from 205.3242902 to **195.247466**. The
remaining gap to 196.1692 is genuine discretisation error on a 256×64 grid
(`dv = 1/63` is very coarse), not a readout offset — and the +4.67% snap
error is gone.

**P3 implementation notes — what the two test binaries actually do:**

- **The normal CDF without approximation.** Black-Scholes needs N(x), the
  standard normal CDF. `<cmath>` has no `norm_cdf`, but it has `erfc`, the
  *complementary error function* — and N(x) = ½·erfc(−x/√2) is an exact
  identity, not an approximation. One line, no polynomial fits to defend.
- **Node-aligned test grid.** `extract_result` reads the cell *nearest*
  (spot, v0) — no interpolation. On the smoke grid the nearest stock node is
  $11.80 away from spot 5200; at delta ≈ 0.56 that's ~$6.60 of price error
  from the snap alone — bigger than everything we're trying to measure. The
  tests instead pick node counts so spot and v0 land *exactly* on nodes:
  `ns = 421` → spacing 50 → 5200 is node 104 (and strike 5250 is node 105);
  `nv = 51` → spacing 0.02 → v0 = 0.04 is node 2.
- **Why parity is near-exact on the grid, not just on paper.** The parity
  portfolio (long call, short put) has payoff S − K — *linear* in S. Central
  differences are exact on linear functions, every boundary formula treats
  it exactly, and it has no v-dependence, so the scheme propagates it with
  only "compounding" error: after nt steps the grid holds (1 − r·dt)^nt
  where the formula says e^(−rT) — a gap of order 1e-5 dollars here. So a
  dollar-tight tolerance is fair for parity, while the BS-collapse test
  carries genuine discretisation error and needs a measured tolerance.
- **Tolerances are measured, not guessed.** Each test was run once to record
  its actual error, then the tolerance was fixed just above it. Measured
  (2026-08-09, node-aligned 421×51 grid, nt = 56000): BS-collapse rel err
  1.16e-3 (call) / 1.12e-3 (put) → tolerance 2.5e-3 (~2× measured); parity
  gap 5.4e-6 dollars → tolerance 1e-4 (20× measured — headroom for compiler
  FP differences on rangpur, still ~5 orders below bug-sized). A real bug
  moves the answer by dollars, so the margins trip on regressions, not noise.
  The plan's "rel err < 1e-3" target is for the *reference-sized* grid
  (spacing ~10 vs 50 here); the P4 convergence study demonstrates it.

**P4 measured convergence (re-done 2026-08-10 — the earlier study was
invalid; see below).** Node-aligned ladder, dt quarters per rung, error
measured against the **semi-analytic Fourier price 196.1692090** rather than
against the study's own finest grid:

| grid | dt | price | error vs analytic |
|---|---|---|---|
| 421×33 | 4.464e-6 | 195.635524 | 0.533685 |
| 841×65 | 1.116e-6 | 196.143695 | 0.025514 |
| 1681×129 | 2.790e-7 | 196.155137 | 0.014072 |

Successive error ratios 20.9 and 1.81, so the observed order in dt is **2.19
then 0.43**. Read that honestly: the error falls fast and then *flattens onto
a floor* rather than converging at a single clean rate. Quoting the mean (1.31)
would describe a straight line that is not what the data does. The floor is
consistent with the `v = 0` row being first-order accurate (a one-sided
difference) while the interior is second-order, so the v-direction error stops
shrinking as fast as everything else. Confirming that needs the fourth rung
(3361×257, ~1.3 h), which is in `LADDER` and has not been run.

**Three things were wrong with the earlier version of this study, and all
three are worth understanding.**

1. **The ladder had gone stale against its own config.** It picked `nv` from
   multiples of 25 so that v0 = 0.04 lands on a node — correct when
   `v_max = 1.0`, which is what `reference.cfg` said when the ladder was
   written. P7 changed `v_max` to 0.64 for stability, and 0.04/(0.64/25) =
   **1.5625**, so *every rung* was off-node in variance. The solver was
   warning about it on stderr the whole time, and the script was calling
   `subprocess.run(..., capture_output=True)`, which swallowed the warning.
2. **The strike was mid-cell on the coarse rungs.** The payoff has a kink at
   K = 5250. For the kink to be discretised the same way on every rung, ds
   must divide it: at ds = 200 the strike sits at node 26.25. A kink that
   moves relative to the grid leaves a sawtooth on top of the trend, which is
   exactly why an early re-run produced an error sequence that went *down and
   then back up*. ds must divide gcd(5200, 5250) = 50, so the ladder now
   starts at ns = 421.
3. **The study was anchored on itself.** Defining the error as
   |price − price_finest| means the finest rung gets no error bar, one rung of
   expensive compute produces a reference instead of a data point, and the
   anchor's own error is silently subtracted from every other rung. The
   external Fourier price fixes all three.

`scripts/convergence_plot.py` now derives the alignment requirement from the
config at run time and refuses to start if any rung misses it, re-raises the
solver's off-node warning instead of discarding it, and reports every pairwise
order rather than only the mean. Figure: `results/convergence.png`, raw numbers
alongside as CSV. Regenerate with `make convergence`.

**P7 — the reference grid, finally run (2026-08-09).** `config/reference.cfg`
had never actually produced a number: it specified `nt = 2000`, which is 527×
over the explicit stability bound, so it printed `nan,nan,nan,nan` and the P2
acceptance criterion "price is positive and finite" had been ticked without
anyone running it. With `nt = 980000` and `v_max = 0.64` it gives:

| source | price | delta |
|---|---|---|
| PDE, 2048×512×980000 | **196.1683699** | 0.56717178 |
| semi-analytic (Fourier inversion) | 196.1692 | — |
| Monte Carlo, 400 000 paths | 196.1345 ± 0.5524 | 0.566829 ± 0.005377 |

PDE against the semi-analytic price: **−0.0008 dollars, a relative error of
4.2e-6**. Three routes to the same number — a deterministic PDE, a
deterministic integral transform, and actual random paths — agreeing to four
parts per million. This is the accuracy claim to make on camera, and note
what makes it *credible*: the Fourier price and the Monte Carlo price are
independent of each other and of the grid, so if the PDE had disagreed, the
mesh would have been the thing to blame, not the arithmetic.

**Do not quote a runtime for this run on camera.** It was measured on the
laptop, not on rangpur, and there is no `sbatch` job or stamped file behind it,
so it must not appear beside cluster numbers. A single reference solve exceeds
the 15-minute QOS limit by a wide margin and has never been run there. If a
runtime is wanted, derive it from the measured rangpur throughput and say you
are deriving it: about 2.1 hours at level 6, about 5.4 hours for the baseline
(RESULTS.md §8).

Two footnotes worth knowing if asked. The reference grid is **not**
node-aligned (`ds = 21000/2047` puts spot at 506.9 cells), so this number
depends on the P7 bilinear interpolation — with the old nearest-node readout
it would have been off by roughly delta × half a cell ≈ $2.9. That path is now
gated: `make check-heston` leg 2 runs a deliberately off-node grid and checks
the interpolated price against the analytic one, because until 2026-08-10 every
automated test used node-aligned grids where the blend weights are exactly zero
and the bilinear code never ran. And the domain
truncation at `v_max = 0.64` (80% annualised volatility, the highest VIX print
in modern history; P(v_T > 0.64) ≈ 2e-15 under this calibration) contributes
about 2e-4 dollars, three orders of magnitude below the agreement above.

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

Python analogies for the C++ constructs (in code, comments explain each new
construct to a beginner in plain English first, ending "similar to X in
Python" only where it helps — never a labelled "Python analogy:" comment):

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

(In the code these carry descriptive names: `index(stock_i, var_j) =
var_j * num_stock_nodes + stock_i`, spacing `stock_spacing`/`variance_spacing`,
sheets `current()`/`next()`. The maths shorthand `i`, `j`, `Ns` below means
the same thing.)

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
randomly and every load pays full price.

**A correction the measurements forced, and the most honest thing in this
guide.** This section used to end "cell-update arithmetic is only ~20 flops,
so this solver's speed is governed by memory access patterns, not arithmetic."
P7 built the roofline to back that claim up and it did not survive:

| | Mcell-updates/s | GFLOP/s (×31 flops/cell) | share of measured scalar peak |
|---|---|---|---|
| baseline, 4 KiB pages | 44.8 | 1.39 | 27% |
| baseline, huge pages | 52.6 | 1.63 | 32% |
| level 6, 4 KiB pages | 81.1 | 2.51 | 49% |
| level 6, huge pages | 134.8 | 4.18 | **82%** |
| level 6 @4096×1024 | 156.9 | 4.86 | **95%** |

The peak in that last column is **measured, not quoted**: 5.117 GFLOP/s from a
register-only multiply-add microbenchmark compiled exactly like the solver
(`slurm/bench_peak.sh`). And the memory side of the roofline is nowhere near
binding — the optimised kernel demands about 3 GB/s while moving 24 bytes per
cell, on a machine with four DDR4 channels.

So the accurate statement is:

- **The baseline is bound by division throughput**, not by memory. It runs at
  the *same* ~45 Mcell-updates/s whether the working set is 1 MiB or 64 MiB —
  a completely memory-insensitive profile. Five `divsd` per cell at ~13–20
  cycles each is what sets its speed.
- **After optimisation the kernel is essentially AT scalar compute peak.** 82%
  at the reference grid and 95% at the largest one, against a ceiling measured
  on the same machine. There is almost nothing left for any memory technique
  to recover, which is exactly why levels 3, 4 and 5 came back null.
- **Memory still matters, but at the margin and not through bandwidth.** The
  two costs P7 actually measured are *page size* (1.66×, from TLB pressure)
  and *traversal order* (2.1×) — both latency effects, not bandwidth ones,
  because at the reference footprint the whole 16 MiB working set fits inside
  this node's 30 MB L3. The optimised kernel demands ~3 GB/s from a machine
  with four DDR4 channels; bandwidth is nowhere near binding.

Layout and loop order still matter, and the flat contiguous array is still the
right decision — it is what makes SIMD possible in Milestone 2, and the
`ctl-order` control shows the cost of getting it wrong is real. But "memory-
bound" was an assumption inherited from the shape of the problem rather than a
measurement, and the correct version is more interesting: **optimising the
arithmetic moved this kernel from division-bound to compute-bound, and the
next win has to come from doing more work per instruction — which is precisely
what SIMD is.** That is the natural bridge into Milestone 2.

### 14b. Two pictures that make the memory story provable

Claims like "memory-bound" and "loop order dominates" are easy to assert and
easy to disbelieve. P7 added two artefacts whose job is to make them
checkable. Both are worth being able to narrate.

**`scripts/memory_anim.py --style cache` — the traffic counter.** Two panels
sweep the *same* nine-point stencil across the *same* flat array in the two
possible orders, with a modelled fully-associative LRU cache of 64-byte lines
running underneath. The visual tells the story before any number does: the
row-major panel lights up a **horizontal band** (three rows resident, walking
along them), the swapped-loop panel lights up a **vertical column** (three
lines resident, walking down a 16 KB stride).

The headline metric is deliberately *not* hit rate. Hit rate barely moves —
98.5% versus 92.8% — because consecutive cells in the swapped order still
share six of their nine accesses vertically. The metric that matters is
**bytes moved per cell update: 8.5 B versus 41.2 B**, a 4.7× difference for
byte-identical arithmetic. That is the number the roofline consumes, and it is
the honest way to say what a cache line costs you: a miss drags all 64 bytes
across the bus however few of the eight doubles you use before eviction.

Be honest about the model's limits on camera: it simulates **capacity only**.
The measured slowdown on real hardware is larger than 4.7×, and the gap is
hardware prefetching (a contiguous stream is predictable; a 16 KB stride is
not) and TLB reach, neither of which the model simulates. Quoting the model's
4.7× as if it explained the whole measured slowdown would be overclaiming; the
animation's caption says exactly this.

**`scripts/roofline.py` — the test the memory-bound claim failed.** Flops per cell are counted
**term by term off the level-6 kernel**, not estimated: 15 for the five finite
differences (with `2.0*V` counted once, since it is a shared subexpression),
16 for assembling the PDE and stepping — **31 flops per interior cell**. Bytes
are counted the same way, with perfect row-major reuse: read one new double of
`cur` (8 B), write one double of `next` (8 B), plus 8 B because the store
misses and its line is *read* before being overwritten (write-allocate — the
detail people forget). **24 bytes per cell.**

Arithmetic intensity = 31/24 ≈ **1.3 flop/byte**. Achieved GFLOP/s is just the
measured cell-updates/sec times 31. The second ceiling on the plot is the nice
part: rather than quoting a vendor bandwidth figure nobody measured, it derives
an **empirical** streaming ceiling from our own `ctl-order` control — that
kernel uses roughly one double per 64-byte line, so its measured throughput
times 256 B/cell is a bandwidth the machine *demonstrably sustained*. It is a
lower bound, it is ours, and it needs no datasheet.

**`scripts/memory_anim.py --style buffers`** covers the other half of §13 and
§15: the ping-pong (read `current`, write `next`, `std::swap` exchanging two
pointers in O(1) with not one double copied), then the ladder's data
structures appearing rung by rung — row constants, the S and S² tables, row
base offsets, row pointers, the four-wide unroll — which is the clearest way
to show that the ladder is *cumulative*.

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

**How we prove which technique did what: the ablation ladder.** An
"ablation" is a controlled experiment borrowed from science: change ONE
thing, measure, repeat. `--opt-level 0..6` turns the list above into a
cumulative ladder — level 0 is the baseline algorithm, level 1 adds the
lookup table, level 2 adds strength reduction on top, and so on. Each rung
is benchmarked separately, so the bar chart shows the *marginal* gain of
each named technique, not just one before/after pair. Two rules make the
experiment honest:

- The level choice happens ONCE, before the time loop, by picking a
  per-level kernel function (`step_level3()` etc.). If instead we checked
  `if (level >= 3)` inside the hot loop, the checking itself would slow
  every level down and corrupt the measurement — the classic
  instrument-disturbs-the-experiment mistake.
- A rung that shows **no speedup is reported anyway**. It usually means the
  compiler already did that optimisation at `-O2` (likely for induction
  variables and unrolling). The task sheet explicitly gives credit for
  "techniques you tried that didn't work" — a null result demonstrates
  methodology.

Interview line: *"Level N enables the first N techniques cumulatively;
diffing two adjacent kernels shows exactly one technique, and
`test_opt_matches` proves every level still produces the baseline answer —
bit-for-bit through level 1, and to a measured 3.3e-15 relative from level 2
up, where reciprocals legitimately move the last bits."*

The metric: **cell-updates per second** = `Ns·Nv·nt / seconds` — comparable
across grid sizes, which raw seconds are not.

**How the ladder is wired:**

- **Dispatch once, via a function pointer.** `StepKernel` (in
  `solver_opt_kernels.h`) is a *function pointer* — a plain variable whose
  value is "which function to call". `OptSolver::solve` asks
  `kernel_for_level(cfg.opt_level)` for it ONCE, before the time loop, so
  the hot loop never pays a per-cell "which level?" branch. Similar to
  assigning a function to a variable in Python — but resolved to a raw
  jump address, no dictionary lookup per call.
- **A kernel owns the interior + the v=0 row.** The S=0/S=max columns and
  the v=vmax Neumann row are shared plumbing in `OptSolver::solve`,
  identical for every level and for baseline — ladder timings compare
  kernels, nothing else.
- **Level 0 is the anchor.** `step_level0` is line-for-line the baseline
  arithmetic: bit-identical output, and measured 44.6 vs 44.8 Mcell-updates/s
  against baseline on rangpur — proof the harness itself costs nothing.
- The opt solver's CSV label carries the level (`opt-L0`, `opt-L3`), and the
  two controls are labelled `opt-ctl-order` / `opt-ctl-branch` rather than
  `opt-L7/L8`, so a benchmark row can never be misread as the ladder
  continuing upward.

#### The measured ladder (rangpur r730-2, Xeon E5-2670 v3, g++ 8.5.0 -O2)

2048×512, median of 5 reps, huge-page regime (see "the page-size trap" below
for why that qualifier matters). **This is the result, and it is not the
result anyone expects:**

| level | technique added | Mcell-updates/s | step | vs baseline |
|---|---|---|---|---|
| baseline | reference solver | 52.6 | — | 1.00× |
| L0 | none (anchor) | 52.6 | −0.1% | 1.00× |
| L1 | hoisting + lookup tables + CSE | 52.7 | +0.3% | 1.00× |
| **L2** | **strength reduction** | **131.2** | **+148.7%** | **2.49×** |
| L3 | traversal order | 131.5 | +0.3% | 2.50× |
| L4 | loop splitting | 131.3 | −0.2% | 2.50× |
| L5 | induction variable | 134.2 | +2.2% | 2.55× |
| L6 | unrolling ×4 | 134.8 | +0.5% | **2.56×** |

**One rung out of six did all the work.** The honest headline is 2.56× end to
end, and essentially all of it is level 2.

**Why — and this is the sentence to have ready.** GCC 8.5 at `-O2` had already
done levels 1, 3, 4, 5 and 6 for us. Hoisting loop-invariant expressions,
eliminating common subexpressions, strength-reducing induction variables and
unrolling counted loops are all *value-preserving* transformations, so the
compiler is free to apply them and does. **Level 2 is the only rung that is
not.** Replacing `x / (2*ds)` with `x * (1/(2*ds))` produces a different
double — `1/(2*ds)` is itself rounded — so without `-ffast-math` the compiler
is *forbidden* from making that change, no matter how obviously profitable it
is. The programmer has to decide the accuracy cost is acceptable and make it
by hand.

So the rule the ladder actually demonstrates is: **at `-O2`, the optimisations
worth doing by hand are the ones the compiler is not allowed to do for you.**
That is a far more useful takeaway than "we got 1.8×", and it is only visible
because every rung was measured separately instead of shipping one
before/after pair.

The cost is real, too, and it is measured rather than hand-waved: level 2
moves the answer by 3.3e-15 relative (about 15 ulp) over 5200 timesteps —
`test_opt_matches` reports it every run. That is the trade the video should
state out loud: **2.5× throughput for 15 units in the last place.**

#### The page-size trap (the most surprising thing P7 found)

The ladder above says "huge-page regime" because there are two of them, and
the difference is 1.66×.

`Grid` allocates two 8 MiB buffers per solve and frees them at the end of each
`--bench` rep. Measured behaviour on this node, reproducible to three digits:

| | rep 1 | reps 2–5 |
|---|---|---|
| level 6 @2048×512 | 134.8 Mcell-updates/s | 81.1 |

The obvious suspects are both wrong, and both were checked rather than
assumed:

- **Not the clock.** `slurm/bench_peak.sh` samples `/proc/cpuinfo` *during* an
  actual solve and finds 3100 MHz pinned the whole way through. No turbo
  ramp-down.
- **Not the cache.** The working set does not change between reps.

It is the **allocator**. The first allocation is served by `mmap`, and with
`/sys/kernel/mm/transparent_hugepage/enabled = [always]` it is backed by 2 MiB
transparent huge pages — 16 MiB then needs about **8 TLB entries**. glibc
afterwards raises its *dynamic* `mmap` threshold past 8 MiB, so later reps are
served from the reused main heap on 4 KiB pages — the same 16 MiB now needs
**4096 entries**, far past this core's L2 TLB, and the stencil pays page walks.

Proved by pinning the threshold by hand (`slurm/diag_alloc.sh`):

| setting | rep 1 | reps 2–5 |
|---|---|---|
| default | 134.8 | 81.1 |
| `MALLOC_MMAP_THRESHOLD_=1048576` (force mmap) | 134.3 | **134.2 — the drop disappears** |
| `MALLOC_MMAP_THRESHOLD_=134217728` (force heap) | 134.8 | 81.1 |

**Why it matters for what gets reported.** A real production solve is a fresh
process that allocates once — so it lives in the *fast* regime, and a plain
median-of-5 understates real performance by 40% and understates the speedup
(1.81× instead of 2.56×). Both regimes are measured and both are in
`results/bench_summary.md`; the huge-page one is the headline because it is
the one a user experiences.

The transferable lesson, and a good interview answer: **"median of 5 reps" is
not automatically the honest statistic.** It is honest when the reps are
samples of the same thing. Here rep 1 and reps 2–5 were sampling two different
memory regimes, and averaging them would have hidden both. The fix was not a
better statistic — it was finding out *why* the reps disagreed.

#### The negative controls: what levels 3 and 4 are actually worth

Levels 3 and 4 could not show a gain because the baseline never had the loss.
So the ladder ships two kernels that *do* have the loss, measured in the same
job, at the same `nt`, on the same node (a cross-job comparison would measure
run length, not technique):

| kernel | what it ablates | Mcell-updates/s | vs paired reference |
|---|---|---|---|
| level 6 (paired reference) | — | 135.3 | 1.00× |
| `ctl-branch` | v=0 row fused back in behind a per-cell `if` | 135.3 | 1.0× — free |
| `ctl-order` | loops swapped → `ns`-element stride | 63.6 | **2.1× slower** |

**Both controls came back smaller than a naive reading of the cache model
would predict, and the reason is the most useful thing this project learned
about the machine.** The reference grid is two 8 MiB buffers — **16 MiB
against this node's 30 MB L3**. The whole problem fits in last-level cache. So
the "wrong" traversal still finds its data in cache; it just finds it one
level further away. Swapping the loops moves roughly 4.8× more bytes (measured
by the cache model in §14b) and costs 2.1× the time, not 4.8×, because at this
footprint the extra bytes are served by L3 rather than DRAM.

That gap between *traffic* and *time* is worth saying out loud rather than
hiding: a 4.8× traffic penalty costing 2.1× is evidence that the kernel is not
bandwidth-limited here, which is the same conclusion the roofline reaches from
the other direction.

`ctl-branch` being free is the expected result for different reasons: the
branch is `var_j == 0`, false for 510 of 511 rows, so the predictor gets it
right essentially always. "Hoist the conditional out of the loop" was
excellent advice for hardware without speculative execution. Reporting that it
buys ~1% today is the honest answer, and it is *still* worth doing because it
costs nothing and makes the loop bodies simpler.

**A caution about local numbers.** On the Mac, `ctl-order` measured ~39×
slower. That figure was wrong: the reference solve was running concurrently
and a bandwidth-hungry kernel is exactly what suffers most from a co-tenant.
The rangpur number was taken on an exclusively-allocated node and is the one
that counts. This is why `--exclusive` is in every benchmark job — and it is a
good story about how easy it is to measure something real that answers a
different question than the one you asked.

#### Sweep C: is any of this just the compiler's auto-vectoriser?

| solver | `-O2` | `-O2 -fno-tree-vectorize` | vectoriser's share |
|---|---|---|---|
| baseline | 44.8 | 44.8 | 1.00× |
| opt level 6 | 81.1 | 81.1 | 1.00× |

(Both columns are default-allocator numbers, because sweep C rebuilt and re-ran
with the allocator left alone — pairing it against the huge-page ladder would
have credited the vectoriser with the 1.66× page effect. Getting that pairing
wrong is an easy mistake and the summary script now guards against it.)

Identical to three digits. None of the speedup is hidden auto-SIMD — which
matters because Milestone 1 forbids SIMD and Milestone 2 is *about* SIMD, so
the M1 number has to be an honest scalar baseline for M2 to improve on. Note
*why* the control is a no-op rather than a discovery: GCC 8.5 enables
`-ftree-vectorize` at `-O3`, not `-O2`, so there was nothing to switch off.
That is still worth running and reporting — "we assumed it was off" and "we
measured that it was off" are different claims.

#### Sweep B: how the speedup scales with problem size

| grid | working set | baseline | level 6 | speedup |
|---|---|---|---|---|
| 512×128 | 1 MiB | 44.5 | 143.6 | 3.23× |
| 1024×256 | 4 MiB | 41.9 | 100.4 | 2.39× |
| 2048×512 | 16 MiB | 44.8 | 81.1 | 1.81× |
| 4096×1024 | 64 MiB | 57.8 | 156.9 | 2.71× |

Two things to read off this, and one trap:

- **The baseline is flat.** ~42–58 Mcell-updates/s across a 64× range of
  working-set sizes. A kernel whose speed does not care whether its data fits
  in cache is not memory-bound; the baseline's five divisions per cell set its
  pace and nothing else gets a look in.
- **Level 6 is not flat.** It falls from 143.6 to 81.1 as the working set
  grows from 1 MiB to 16 MiB. Having removed the division bottleneck, the
  kernel is now fast enough to notice the memory hierarchy — the optimisation
  *moved* the bottleneck rather than just raising the ceiling.
- **The trap:** the 4096×1024 row goes back UP, which looks like the memory
  effect reversing. It is not. Those buffers are 32 MiB each, which is at
  glibc's `DEFAULT_MMAP_THRESHOLD_MAX`, so *every* rep is served by a fresh
  `mmap` and gets huge pages, while the three smaller grids fall back to 4 KiB
  pages for reps 2–5. That row is measured in a different memory regime and is
  annotated as such on the figure. It is the page-size effect again, wearing a
  disguise — and a good example of why an unexplained non-monotonicity is
  worth chasing rather than shrugging at.

Two secondary observations worth a sentence each:

- **The Mac and the Xeon disagree about unrolling.** On Apple Silicon (via
  Rosetta) level 6 was a *23% regression* — four cells in flight need more
  live registers than the architecture exposed there, so the unroll spilled.
  On the Xeon it is +0.8%, i.e. a null. Same source, opposite sign: an
  optimisation is a property of a code/compiler/machine triple, never of the
  code alone. This is exactly why the rubric insists numbers come from the
  target machine.
- **Level 1 being null is not the same as level 1 being pointless.** It is
  the rung that makes levels 2–6 *expressible* — the reciprocals of level 2
  have to be hoisted somewhere, and the S and S² tables are what let the
  inner loop become a pure streaming operation. The compiler reproduced its
  effect; it did not make the code that follows possible.

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

**What the benchmark machine actually is, and why it matters (P7).** rangpur's
`cosc3500` partition is not one machine — it mixes 24-core Intel Xeon
E5-2670 v3 boxes (`r730-*`, Haswell, 2014) with 8-core virtual machines
(`a100-*`, `vcpu-*`). Three things followed from that, all of which are worth
being able to explain:

- **Every benchmark job pins to the same physical node** (`--nodelist=r730-2`)
  **and takes it exclusively** (`--exclusive`). Without the pin, the ladder
  might be measured on a Xeon and the controls on a VM, and the comparison
  between them would be meaningless. Without `--exclusive`, another job on a
  different core of the same socket competes for the *shared* memory
  bandwidth and L3 — it never touches our core, and it still changes our
  numbers.
- **The QOS caps a job at 15 minutes.** A single job covering everything sits
  in the queue forever with reason `QOSMaxWallDurationPerJobLimit` — it is not
  rejected, it just never starts, which is a much more confusing failure. The
  suite is five sweeps plus two diagnostic jobs, each writing its own CSV.
  RUNNING.md lists them and how to submit them.
- **That node is slow, and slow in an informative way.** It turns the baseline
  over at ~4.7e7 cell-updates/s, roughly a tenth of a modern laptop core, at a
  measured **3.1 GHz** — `bench_peak.sh` sampled `/proc/cpuinfo` *during* a
  real level-6 solve and found the clock pinned at 3100 MHz throughout, so
  this is not a throttling story. (An earlier draft of this section said
  1.6 GHz, which was the idle floor read before that job existed; the
  cycles-per-cell figure below is corrected accordingly.) That is
  ~**66 cycles per cell update**. The baseline performs **five divisions per
  cell**, and a Haswell `divsd` is ~13–20 cycles and poorly pipelined — five
  of them alone account for 65–100 cycles, which is the whole budget. So the
  *baseline* on this machine is division-bound before it is anything else,
  which is exactly why strength reduction (level 2) pays so well here, and the
  corrected clock makes that argument stronger rather than weaker.

**A debugging habit worth stealing: check the surprising number before
believing it.** When rangpur first came back 10× slower than the laptop, the
obvious suspect was denormal arithmetic — `bench.cfg` uses a very short
maturity, so much of the grid holds values near 1e-300, and x86 takes a
microcode assist of tens of cycles on every operation that touches a denormal.
That would have been a benchmarking artefact, not a result. So it was tested
(`slurm/diag_denormal.sh`) rather than assumed: the same run with `-ffast-math`
(which links `crtfastmath.o` and sets FTZ/DAZ so denormals flush to zero), with
the grid full of NaNs, and with values around 1e33, all produced the *same*
throughput to within noise. The hypothesis was wrong, the hardware is simply
old — and now that is a measured statement instead of a guess. The same
diagnostic also justifies the short `nt` used throughout the benchmarks: since
the data content provably does not affect throughput, a 400-step run measures
the same *rate* as a 400 000-step one.

**Self-test:** Is the solver memory-bound or compute-bound — and what
measurement settles it?
Why is `j*Ns + i` the right index formula and not `i*Nv + j`? What does the
buffer swap cost? Why median of 5 runs on a compute node? Why does
`--exclusive` matter when we only use one core? What is the difference between
a job that is *rejected* and a job that never starts?

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
  `aligned_alloc`. The lectures' rule of thumb is **~4×, not the full vector
  width**, because such kernels are usually memory-bound. **Careful: we
  measured that and it is not true here** (§14). Level 6 already runs at 82%
  of this node's measured scalar peak while asking for only ~3 GB/s, so this
  kernel is compute-bound and the thing standing between us and more speed is
  arithmetic throughput, which is exactly what SIMD widens. So expect closer
  to the vector width than the rule of thumb predicts — and then measure it,
  because that is how the rule of thumb got caught out in the first place. A
  float32 variant (8-wide) is the natural follow-on experiment.
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

**P2 checkpoint drill (2026-08-09 — all five stumbled; re-explain at next
session start, then re-quiz until they stick)**

- **Q: Why does `init_payoff` write into `current_` and not `next_`?**
  A: The time loop always READS `current()` and WRITES `next()`. The payoff
  is the first sheet the loop ever reads, so it must start in the read
  buffer. (§13)
- **Q: With one buffer instead of two, what goes wrong mid-timestep?**
  A: Updating in place means some neighbours you read were already updated
  *this* step — the stencil mixes time level n with n+1. Results then depend
  on loop order and are simply wrong. Two buffers freeze the read sheet for
  the whole step. (§13)
- **Q: Why a FORWARD difference for V_v on the v=0 row — why not central?**
  A: Central needs a node at v = −dv, which doesn't exist and can't (negative
  variance is meaningless). And the drift kappa*theta > 0 points INTO the
  domain, so taking the derivative from above (upwind — from where the
  information flows) is the stable, physical choice. Imposing a value
  (Dirichlet) would over-determine the degenerate equation. (§9, PLAN §1c)
- **Q: When dt crosses the stability bound, which cells blow up first, why?**
  A: The high-S, high-v corner. The diffusion coefficients (½vS², ½ξ²v) are
  largest there, so the explicit update's amplification factor crosses 1
  there first; the error flips sign cell-to-cell each step — the
  checkerboard. (§8)
- **Q: Feller is violated for our parameters — why not "fix" it, and what
  does it actually cost us?**
  A: 2κθ ≥ ξ² is a property of the market-calibrated parameters, not of our
  code — equity calibrations typically violate it. Consequence: the variance
  process puts real probability mass near v=0, where our uniform v-grid is
  relatively coarse → slight under-resolution there. The scheme stays valid;
  we state it as a known limitation (non-uniform v-grid = future work).
  (§9, PLAN §1c)

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
