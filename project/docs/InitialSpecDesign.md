# Pricing an Option When Volatility Won't Sit Still — COSC3500 Project Brief

**Course:** COSC3500 High-Performance Computing
**Student:** [Name and student number]
**Status:** Approved

## 1. What this project actually does

The goal is simple to state: given one stock option contract — for example, a call option on the S&P 500 with a strike of 5,250 that expires in three months — compute what it is fairly worth today, and how sensitive that value is to changes in the market (the "Greeks": delta, gamma, and vega). The catch is the model. In the classic Black-Scholes world, volatility is a fixed number. In the real world, volatility itself moves around — calm markets become wild ones overnight. This project uses a Heston-style stochastic-volatility model, where volatility is a second random quantity of its own, and prices the option numerically.

This is a good high-performance computing project because the numerical method is a large grid of cell updates repeated thousands of times, where every cell depends on its neighbours from the previous step. It is heavy enough to be worth optimising, and — unlike plain Monte Carlo, which is embarrassingly parallel — the dependencies between cells mean the parallel versions have to be earned. The implementation ladder is: a correct serial C++ solver, then a serial-optimised version, then SIMD, then OpenMP, with CUDA as an optional extension.

## 2. The model: two random processes tied together

The Heston model describes how two things evolve at the same time: the stock price `S` and its instantaneous variance `v` (variance is volatility squared). Written in plain text:

```text
dS = (r - q) * S * dt  +  sqrt(v) * S * dW1
```

In words: over a tiny time interval, the stock drifts upward at the interest rate `r` (minus the dividend yield `q`), and jiggles randomly — and the size of the jiggle depends on the current variance `v`. When `v` is high, the stock moves violently; when `v` is low, it barely moves.

```text
dv = kappa * (theta - v) * dt  +  xi * sqrt(v) * dW2

correlation between dW1 and dW2 = rho
```

In words: variance gets pulled back toward its long-run normal level `theta` at speed `kappa`, while being kicked around randomly with strength `xi`. The correlation `rho` is typically negative, which captures a well-known market fact: stock crashes and volatility spikes tend to arrive together.

The four "personality parameters" of the model, one sentence each:

- `kappa` — how fast volatility gets pulled back to normal after a spike.
- `theta` — the long-run variance level the market keeps returning to.
- `xi` — how violently variance itself jumps around (the "vol of vol").
- `rho` — how tightly stock drops and volatility spikes are linked (negative in practice).

This second random process is exactly why the computation gets bigger. Under Black-Scholes, volatility is a known constant, so the solver only needs to ask "what if the stock were at various prices?" — a single row of what-ifs. Under Heston, volatility is *itself* a what-if, so the solver needs a whole spreadsheet of scenarios: every combination of "what if the stock were here" and "what if the market were this calm or this wild".

![One simulated path of the two coupled Heston processes: the stock price and its volatility](../assets/heston-paths.gif)

*One simulated three-month path of the two coupled processes, using the reference parameters from section 8. The lower panel is the point: volatility is not a constant but a jagged random path of its own, and with `rho = -0.70` its climbs line up with the stock's slides.*

## 3. The pricing equation, gently

From the two equations above, standard financial mathematics derives one bookkeeping rule that the option's value `V` must obey at every point in time and in every scenario. It is a partial differential equation (PDE); here `d` stands for a partial derivative:

```text
dV/dt + (r - q) * S * dV/dS + kappa * (theta - v) * dV/dv
      + (1/2) * v * S^2 * d2V/dS2
      + (1/2) * xi^2 * v * d2V/dv2
      + rho * xi * v * S * d2V/dSdv
      = r * V
```

Each term has a plain meaning:

- `dV/dt` — how the option's value decays as time passes.
- `(r - q) * S * dV/dS` — the drift of the stock pushing value along the price direction.
- `kappa * (theta - v) * dV/dv` — the pull of variance back to normal pushing value along the volatility direction.
- `(1/2) * v * S^2 * d2V/dS2` — curvature in the price direction (this is where gamma lives).
- `(1/2) * xi^2 * v * d2V/dv2` — curvature in the volatility direction.
- `rho * xi * v * S * d2V/dSdv` — the cross-term linking price moves and volatility moves; this is the term that makes diagonal neighbours matter in the numerical method.
- `r * V` — discounting: money later is worth less than money now.

We never solve this equation symbolically. Its only job in this project is to supply the blending weights for the numerical method below. From here on, the document speaks entirely in spreadsheet terms.

## 4. How the answer is computed: filling in the spreadsheet backwards

Picture one big spreadsheet. Each **column** is a candidate stock price: "what if the S&P were at 4,800? 4,900? 5,000? ...". Each **row** is a candidate variance level: "what if the market were very calm? normal? panicking?". Each **cell** holds the answer to one question: "what would this one option be worth if the market were in that scenario?". The whole sheet describes a single fixed contract — one underlying, one strike, one expiry — under many hypothetical states of the world.

There is one such sheet for every point in time between today and expiry, and the trick is that we fill them in **backwards**:

1. **The expiry sheet is free.** At the moment the option expires there is no uncertainty left: a call is worth `max(stock price - strike, 0)` and a put is worth `max(strike - stock price, 0)`. Every cell of the final sheet can be written down instantly, no model required.
2. **Step back one small time interval.** To fill in a cell of the earlier sheet, the solver takes a weighted blend of nine cells from the later sheet: the same cell and its eight immediate neighbours. Left and right neighbours are slightly lower and higher stock prices; the cells above and below are slightly calmer and wilder variance levels; the diagonal neighbours carry the price–volatility correlation from the cross-term in the PDE. The blend weights come directly from discretising the pricing equation, and they are the same for a given row, so they can be precomputed once.
3. **Repeat** — thousands of small steps, each producing a full earlier sheet from a full later sheet — until the solver arrives at today.
4. **Read off the answer.** Find the cell closest to today's actual stock price and today's actual variance: that cell is the option's fair price. The slope between neighbouring cells in the price direction gives delta, the curvature gives gamma, and the slope in the variance direction gives vega. The Greeks come for free from the finished sheet.

```text
   tomorrow's sheet (already filled in)         today's sheet (being filled in)

              $90   $95  $100  $105  $110
   wilder   |     |     |     |     |     |
   market   |     | [*] | [*] | [*] |     |
            |     | [*] | [*] | [*] |     |    ...  [ C ]  <- one cell
   calmer   |     | [*] | [*] | [*] |     |
   market   |     |     |     |     |     |

   To fill in cell C ("stock at $100, medium volatility, today"),
   blend the nine starred cells from tomorrow's sheet.

   expiry sheet  <── ... <── tomorrow's sheet <── today's sheet
   (known for free)         (computation moves right to left)
```

![The option value surface evolving backwards in time from the jagged payoff at expiry to the smooth surface today](../assets/value-surface.gif)

*The whole spreadsheet drawn as a surface: stock price and volatility along the base, option value as height and colour. The animation runs the way the solver runs — it starts at expiry, where the surface is the sharp hockey-stick payoff and volatility is irrelevant, and steps backwards to today, where the kink has been smoothed out and higher volatility visibly lifts the value. (Illustrative prototype from a small Python reference solver; the Milestone 1 "weather map" deliverable regenerates this from the C++ engine.)*

In memory the solver keeps exactly two sheets at once. The later sheet is read-only during an update; the earlier sheet is written, each cell exactly once; and only when every cell is finished are the two buffers swapped. Nothing ever reads a half-finished sheet.

One honest caveat: this "explicit" scheme is only stable if the timesteps are small enough relative to the cell spacing. Take steps that are too big and the sheet blows up into oscillating garbage. That is not a flaw to hide — measuring exactly where the stable/unstable boundary sits, and the accuracy-versus-runtime trade-off around it, is part of the research plan. An implicit or ADI scheme is a possible extension, not required scope.

![The error field exploding into a checkerboard when the timestep exceeds the stability limit](../assets/instability.gif)

*The instability, previewed: the same solver run at 4× the stable timestep, plotted as its difference from a stable reference run. For a while nothing seems wrong — then a checkerboard mode erupts in the high-volatility, high-price corner (where the diffusion coefficients are largest) and swallows the sheet. This is the failure the stability measurements in the research plan map out.*

## 5. Where the parallelism lives

The dependency structure of the method is the whole story. Within one timestep, every cell of today's sheet reads only from tomorrow's finished sheet, so **all cells can be computed simultaneously, in any order**. But across timesteps, today's sheet cannot begin until tomorrow's sheet is complete, so **the timesteps themselves are strictly sequential**. That combination — massive parallelism inside each step, a hard synchronisation between steps — is what makes this a genuine parallel-computing problem rather than an embarrassingly parallel one.

The implementation ladder climbs that structure one rung at a time:

**Serial, then serial-optimised.** The first version is a straightforward correct C++ solver. The optimised version stores each sheet as one flat row-major array so that walking along a row walks contiguously through memory, precomputes the nine blend weights per row instead of recomputing them per cell, chooses the loop order to match the memory layout, and pulls the boundary cells out of the hot inner loop. Evidence: profiling, plus before/after cell-updates-per-second.

**SIMD.** The interior cells of a row sit next to each other in memory and all apply the same blend formula, so one vector instruction can compute four to eight adjacent cells at once. Evidence: bit-level or tolerance-level match against the serial answer, then vector speedup and where memory bandwidth caps it.

**OpenMP.** Each thread takes a band of rows of today's sheet. Because every thread only reads the shared read-only later sheet and writes its own disjoint cells, no locks are needed — just one barrier per timestep before the buffer swap. Evidence: strong scaling and parallel efficiency versus thread count.

**CUDA (optional).** One GPU thread per cell, with both sheets resident on the device so nothing crosses the PCIe bus inside the time loop. Evidence: CPU/GPU crossover point by grid size. MPI slab decomposition with halo exchange is a further optional extension.

At every rung the same discipline applies: prove the output matches the serial reference first, then measure.

## 6. How we know it's right

Three named correctness checks anchor the project:

1. **Collapse to Black-Scholes.** Set the vol-of-vol `xi` to zero and start variance at its long-run level `theta`, and the Heston model degenerates into Black-Scholes — which has an exact textbook formula. The solver's price must converge to it.
2. **Put-call parity.** For European options, `C - P = S * e^(-q*T) - K * e^(-r*T)` must hold regardless of the model. Pricing a call and a put on identical settings and checking this identity catches whole classes of bugs.
3. **Grid refinement.** Halving the cell spacing and the timestep must make the price settle toward a stable value; the convergence plot is a deliverable, not an afterthought.

And one iron rule: **no parallel version is ever benchmarked until it matches the serial answer to a stated tolerance.** A fast wrong answer is worth nothing.

## 7. Making it fun

**The weather map (Milestone 1 video — required scope).** Colour the spreadsheet by option value and it becomes a weather map of option prices. Dumping the grid every few hundred timesteps and stitching the frames with matplotlib/ffmpeg gives an animation of the value surface flowing backwards from the jagged payoff at expiry to the smooth surface today. It is cheap to build, makes the ten-minute video visually alive, and doubles as a debugging tool — instabilities are visible on the map before they show up in the numbers.

**The instability horror show (both videos and the interview).** Deliberately exceed the stable timestep on camera and watch the surface explode into noise, then show the measured stability boundary. It costs one config file and is the most memorable thirty seconds either video will contain — and it demonstrates real understanding of explicit schemes.

**Market-episode replay (Milestone 2 narrative).** No data pipeline needed: hand-pick two parameter sets, one for calm February 2020 (low spot volatility, `v0 ≈ 0.02`) and one for the COVID crash of March 2020 (spot down 30%, `v0 ≈ 0.36`, `rho` more negative), re-run the engine, and narrate how the price and the Greeks lurch. It turns a parameter sweep into a story.

**The stress dial (interview demo — build only after Milestone 1 is safe).** A tiny terminal loop that re-solves the grid live as spot or volatility is nudged, printing solves-per-second — first with the serial solver (sluggish), then with OpenMP (snappy). It converts the speedup number into something the interviewer can feel rather than read off a chart.

## 8. Milestones, scope, and guardrails

**Milestone 1 (serial only, 24%, ten-minute video).** Deliverables: the correct serial solver with boundaries handled, the Black-Scholes and put-call-parity validation, the grid-refinement convergence plot, the serial before/after optimisation numbers, and the weather-map animation. No SIMD — the course counts it as parallelism.

**Milestone 2 (parallel, ten-minute video plus interview).** Deliverables: the SIMD and OpenMP versions with their correctness matches, the scaling and efficiency story, the stress-dial demo, and CUDA only if the CPU evidence is already complete.

All benchmark inputs are synthetic and version-controlled, so every run is reproducible and grid size can be scaled freely. A reference configuration:

```yaml
underlying: SPX          # spot 5200, rate 0.045, dividend yield 0.013
option:    { strike: 5250.0, maturity_years: 0.25, type: call }
heston:    { v0: 0.04, kappa: 1.5, theta: 0.04, xi: 0.35, rho: -0.70 }
grid:      { stock_nodes: 2048, variance_nodes: 512, time_steps: 2000 }
```

**What this project deliberately does not do:** European calls and puts only — no American, barrier, or multi-asset options; no calibration to a real volatility surface; no backtester or trading strategy; no live data feed. A small historical end-of-day comparison against real option quotes is optional garnish, attempted only if lawfully accessible data is already in hand, and it is the first thing cut. If work must be dropped, the cut order is fixed: MPI first, then CUDA, then the historical study — never the solver, the validation, or the benchmarks. If the 2D boundary treatment proves too risky, the documented fallback is a 1D local-volatility PDE solver, which preserves the entire optimisation and parallelisation story on a simpler grid.
