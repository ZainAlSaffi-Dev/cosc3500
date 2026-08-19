# Milestone 1 Build Plan — Heston PDE Option Pricer (Serial)

Spec: `docs/InitialSpecDesign.md`. Rubric: Intro/Background 20%, Optimisation 25%,
Benchmarking 25%, Presentation 20%, Reflection 10%, Code submission pass/fail.
Milestone 1 is **serial only** — no SIMD intrinsics, no OpenMP, no CUDA.

---

## 0. The one-paragraph mental model

We price one European option under Heston by filling a 2-D grid of "what-if"
scenarios (stock price × variance) backwards in time from expiry to today.
Each backward step blends, for every interior cell, nine cells of the later
sheet (3×3 stencil — the diagonal entries carry the price/vol correlation).
Two buffers, swap per step. Price = cell nearest (spot, v0); Greeks = finite
differences around it. Explicit scheme ⇒ conditionally stable ⇒ the stability
boundary is itself a deliverable, not a bug.

## 1. Numerical scheme (decisions locked)

| Decision | Choice | Rationale |
|---|---|---|
| Scheme | Explicit Euler in time, central differences in S and v, central cross-term | Simplest correct; stability study is a deliverable |
| Grid layout | One flat `std::vector<double>`, **row-major with S contiguous**: `idx = j*Ns + i` (`i` = stock index, `j` = variance index) | Inner loop walks contiguous memory; sets up SIMD/OpenMP for M2 |
| Buffers | Exactly two (`cur` read-only, `next` write-only), pointer swap per step | Spec §4; no half-read sheets |
| S domain | `[0, s_max_mult * strike]`, uniform. Default `s_max_mult = 4` | Standard practice; far boundary ≈ linear payoff |
| v domain | `[0, v_max]`, uniform. Default `v_max = 1.0` (vol 100%) | Covers crash scenarios (v0=0.36 in M2 replay) |
| Boundary S=0 | Call: `V=0`. Put: `V = K*exp(-r*(T-t))` | Exact — stock absorbed at 0 |
| Boundary S=Smax | Call: `V = S*exp(-q*τ) - K*exp(-r*τ)` (deep ITM). Put: `V=0` | Asymptotic payoff |
| Boundary v=0 | **Degenerate-PDE row — see §1c.** No condition imposed; solve `V_t + (r−q)S·V_S + kappa·theta·V_v − rV = 0` with central S-diff + forward (upwind) v-diff. Never Dirichlet, never central-in-v | in't Hout & Foulon treatment; drift kappa·theta > 0 flows inward |
| Boundary v=vmax | `dV/dv = 0` (Neumann): copy row `j=Nv-1` update from `j=Nv-2` curvature-free form | Vega saturates for huge v |
| Greeks | delta/gamma: central differences in S at nearest cell to (spot, v0); vega: central difference in v, chain rule ×(2*sqrt(v)) noted in write-up | Free from finished sheet |
| Timestep | User-set `nt`; solver also computes and reports the explicit stability estimate `dt_stable ≈ min over grid of 1/( v*S²/ΔS² + xi²*v/Δv² + r )`-style bound | Enables the instability demo and the stability-boundary study |
| Precision | `double` everywhere | Correctness first; float experiments are M2 territory |
| Config format | Plain `key = value` text file (`config/*.cfg`), hand-rolled parser | Zero dependencies on rangpur |
| Compiler flags | Baseline benchmark: `-O2`. Report flags in every benchmark. Add `-fno-tree-vectorize` variant to prove serial-opt gains are not hidden auto-SIMD | M1 forbids SIMD; keep the story clean |
| Timing | `std::chrono::steady_clock` around the time loop only (exclude I/O); metric = **cell-updates per second** = Ns*Nv*nt / seconds | Comparable across grid sizes |

Weights: for a given variance row `j`, the nine stencil coefficients depend on
(j, i) only through simple per-row and per-column factors. Baseline recomputes
per cell; optimised solver precomputes per row (see §4).

## 1b. Code style — "OOP-shaped modern C++, memory-safe" (locked)

Author knows Python OOP (some Java, some C), is learning C++. Rules for ALL
project C++, with the Python mapping stated so the code doubles as a tutorial:

| C++ construct | Python equivalent | Where used |
|---|---|---|
| `class` with `private:` fields (`name_`), public getters/methods | class with `_attrs` — but privacy compiler-enforced | `Grid` |
| abstract base: `virtual m() = 0;` + `virtual ~T() = default;` | `ABC` + `@abstractmethod` | `Solver` |
| `class B : public A` + `override` | subclassing; `override` catches method-name typos at compile | `BaselineSolver`, `OptSolver` |
| `const` method suffix | no analogue — "doesn't mutate self", compiler-checked | all getters |
| plain `struct` (all public data, no methods) | `@dataclass` | `Config`, `SolveResult` |
| constructor + RAII (destructor auto-frees members) | GC/`with`-block, but deterministic at scope exit | `Grid` owns its vectors |
| `std::unique_ptr<Solver>` via `make_solver` factory | factory function; ownership explicit (Python GC hides it) | polymorphic dispatch only |
| `const Config&` parameter | object passing, but never-null + read-only enforced | everywhere |

Two concepts Python has no instinct for — flag them in comments wherever
they appear:
1. **Value vs reference semantics.** Python variables are all references to
   heap objects. C++ `Grid g` IS the object, on the stack, destroyed at
   scope exit. References are opt-in (`Grid&`). This is why RAII works and
   why the memory is contiguous (fast).
2. **Explicit ownership.** Every byte has exactly one owner: a vector or a
   unique_ptr. Scope exit frees it. No GC needed, no leaks representable.

Memory-safety rules (the point):
- **No raw `new`/`delete` anywhere.** `std::make_unique` in the factory is
  the only allocation of a polymorphic object; `std::vector<double>` owns
  all grid memory.
- **No raw owning pointers, ever.** `.data()` raw pointers appear only
  inside OptSolver's hot kernel, as locals, never stored.
- References over pointers for parameters; nothing is nullable.
- Exceptions only at the edges: `load_config` / `make_solver` throw
  `std::runtime_error`; `main` has the single try/catch.
- No templates of our own, no lambdas, no `std::function`, no inheritance
  deeper than `Solver` → implementations. `auto` only where the type is
  obvious on the same line.
- `std::string` for paths/names; `printf`-family for output (CSV lines).
- One class per header; files stay under ~200 lines.

Naming (locked at P2 start): descriptive identifiers over abbreviations —
`num_stock_nodes` not `ns`, `stock_spacing` not `ds`, `current()`/`next()`
sheets, loop indices `stock_i`/`var_j`. Docs and cfg/CLI keys keep the short
maths shorthand (`ns`, `nv`, `nt`, `idx = j*ns + i`); `load_config` maps the
short keys onto the descriptive fields. Exception: Heston parameters keep
their maths symbols (`kappa`, `theta`, `xi`, `rho`, `v0`) — they must stay
recognisable against the PDE and the literature.

Narration comments (so the author can follow along line-by-line):
- **Every serial-optimisation technique gets a one-line plain-English label
  on the line (or block) where it is applied** — e.g.
  `// loop-invariant hoisting: coefficients depend only on j, computed once per row`,
  `// strength reduction: replace divide with precomputed reciprocal`,
  `// row-major traversal: j outer, i inner, so memory access is contiguous`.
- **Every boundary-condition handling site gets a one-line comment** naming
  which boundary it is and what condition is imposed — e.g.
  `// v=0 boundary: one-sided derivative, diffusion term vanishes here`.
- **Every non-obvious decision gets a one-line comment stating the choice
  and the reason** — e.g. `// two buffers, swap pointers: avoids copying the
  whole grid each timestep`.
- **New-construct comments explain to a beginner first, in plain English**
  — say what the construct does on its own terms; at most end with
  "similar to X in Python". Never a labelled "Python analogy:" comment.
- Keep these to one line each, plain English, on/above the relevant line.
  They are teaching aids for the 10-min video and interview prep, same
  spirit as the new-construct comments above.

Why this is safe: every byte is owned by exactly one vector or unique_ptr,
and scope exit frees it — leaks and double-frees are unrepresentable.

## 1c. The v=0 boundary, spelled out (the trap boundary)

At v=0 every diffusion term carries a factor v and vanishes; the PDE
degenerates to pure transport on that row:

```text
V_t + (r−q)·S·V_S + kappa·theta·V_v − r·V = 0
```

The v-drift at the boundary is kappa·theta > 0 — pointing INTO the domain —
so the correct treatment is to solve this degenerate equation on the j=0
row, not to impose a value:

- 4-point stencil: self, S-left, S-right (central), one v-neighbour above
  (forward one-sided difference for V_v; upwind matches the positive drift).
- FORBIDDEN: Dirichlet at v=0 (over-determines), central v-difference
  (needs a v<0 ghost point — meaningless).

**Feller check (report at startup):** 2·kappa·theta ≥ xi² decides whether
variance can touch zero. Reference params: 0.1200 < 0.1225 — VIOLATED
(marginally), as is typical for equity calibrations. Scheme stays valid;
consequence is solution mass near v=0 that a uniform v-grid slightly
under-resolves. State as a known limitation (video beat); non-uniform
v-grid is a stretch goal, not M1. Solver should print the Feller status
alongside dt_stable_estimate.

Stability note: the v=0 row is the least stiff (its coefficients vanished);
the explicit dt bound binds at the high-v/high-S corner — which is exactly
where the instability checkerboard erupts first.

## 1d. Working agreement — learning by building (locked)

The author wants to finish the project AND understand every line (interview
prep for M2, video narration for M1). Rules for every coding session:

- **Explain before code.** Before any file is written/changed, state in 2–3
  plain-English sentences what it will do and which STUDY_GUIDE section it
  maps to. New concept → STUDY_GUIDE gets a section before the code lands.
- **Author-writes list.** These core pieces the author types themself, with
  guidance (skeleton + hints first, review after): payoff initialisation,
  the interior stencil update, ONE boundary (v=0 — the trap one), and ONE
  opt-level kernel (level 1, hoisting — the biggest idea). Everything else
  (plumbing, CLI, tests, scripts, remaining boundaries/levels) can be
  scaffolded by the assistant, but always walked through afterwards.
  - **WAIVER, P7 session (2026-08-09):** the author explicitly waived the
    `step_level1` reservation for that session, so all of levels 1–6 plus
    the two negative controls were assistant-written. The condition attached
    to the waiver was that the kernels be commented well enough to learn
    them from afterwards — see the block comment at the top of
    `src/solver_opt_kernels.cpp` and the bit-identity rule spelled out above
    `step_level1`. Everything else on this list still stands.
- **Narration comments** (§1b) are mandatory everywhere — they are the
  follow-along thread through the code.
- **The author's comments are theirs.** In author-written pieces, comments
  are personal study notes — review the code, not the comments. Conceptual
  corrections go in conversation, never as requested comment edits.
- **Checkpoint quiz.** At the end of each phase (P1–P8), 3–5 rapid-fire
  questions on what was just built (same style as STUDY_GUIDE Part V).
  Wrong answer → that topic gets re-explained and added to the drill list.
- **STUDY_GUIDE stays in sync.** Any decision change lands in the guide in
  the same session (already a repo rule; restated here because it's the
  learning backbone).

## 2. Repo layout (this scaffold)

```
project/
├── PLAN.md                  ← this file
├── Makefile                 ← build + test + sync + remote-bench targets
├── .gitignore               ← results/, binaries, frames
├── config/
│   ├── reference.cfg        ← spec §8 reference run (SPX 5250 call)
│   ├── demo.cfg             ← node-aligned animation grid (P6); the quoted 196.105
│   ├── bench.cfg            ← reference FOOTPRINT, short maturity: throughput only
│   ├── smoke.cfg            ← tiny grid, seconds to run, for tests/CI
│   └── unstable.cfg         ← 4× stable dt — instability horror show
├── include/
│   ├── params.h             ← Option/Market/Heston/Grid/Config structs + parser decl
│   ├── grid.h               ← flat 2-buffer grid, idx(), swap
│   ├── solver.h             ← SolveResult, solver interface (baseline & opt)
│   ├── solver_opt_kernels.h ← StepKernel function-pointer type + kernel_for_level()
│   ├── black_scholes.h      ← closed form + norm_cdf (validation only)
│   └── io.h                 ← grid snapshot dumps, result CSV line
├── src/
│   ├── main.cpp             ← CLI: --config --solver --dump-every --bench
│   ├── params.cpp
│   ├── grid.cpp
│   ├── solver.cpp           ← shared extract_result() readout + make_solver() factory
│   ├── solver_baseline.cpp  ← correct-first serial solver
│   ├── solver_opt.cpp       ← optimised serial (same interface, same answers); dispatches on --opt-level
│   ├── solver_opt_kernels.cpp ← step_level1()..step_level6() ladder kernels (§4b)
│   ├── black_scholes.cpp
│   └── io.cpp
├── tests/
│   ├── test_bs_collapse.cpp ← xi=0, v0=theta ⇒ match Black-Scholes formula
│   ├── test_parity.cpp      ← C - P = S*exp(-qT) - K*exp(-rT)
│   └── test_opt_matches.cpp ← every --opt-level 0..6 + both controls ≡ baseline
│                              (bit-identical to L1; measured 3.3e-15 rel above)
├── scripts/
│   ├── weather_map.py       ← snapshots → coloured frames → gif/mp4 (ffmpeg)
│   ├── convergence_plot.py  ← refinement study → convergence figure
│   ├── bench_plot.py        ← bench CSVs → opt-level ladder bars + scaling figure
│   ├── memory_anim.py       ← cache-line sweep + double-buffer ping-pong clips
│   ├── roofline.py          ← flops/bytes per cell → roofline; TESTS the bound
│   ├── monte_carlo_check.py ← independent MC + Fourier price (external reference)
│   ├── bs_collapse_anim.py  ← xi→0 collapse to Black-Scholes, animated
│   └── xi_sweep_smile.py    ← implied-vol smile flattening as xi→0
├── slurm/                   ← how to run these: RUNNING.md §2–§3
│   ├── debug_serial.sh      ← 5-min sanity job (build + smoke run)
│   ├── bench_common.sh      ← shared header/run helpers, sourced by the five sweeps
│   ├── bench_ladder.sh      ← sweep A (ladder, default allocator)
│   ├── bench_scaling.sh     ← sweep B (baseline vs level 6 across four grids)
│   ├── bench_novec.sh       ← sweep C (-fno-tree-vectorize rebuild)
│   ├── bench_controls.sh    ← sweep D (ctl-order, ctl-branch + paired reference)
│   ├── bench_pages.sh       ← sweep E (ladder again, huge-page regime — headline 2.56×)
│   ├── bench_peak.sh        ← measured scalar FLOP ceiling for the roofline
│   ├── diag_alloc.sh        ← diagnostic: proved the page-size effect
│   └── diag_denormal.sh     ← diagnostic: ruled out denormal stalls
└── results/                 ← gitignored; snapshots, CSVs, frames land here
                               (RESULTS.md carries the tracked copy of the numbers)
```

One binary `heston`; tests are three tiny extra binaries reusing the same
objects. No external C++ dependencies. Python (matplotlib, ffmpeg) local only —
never needed on rangpur.

## 3. CLI contract (fixed now so scripts can be written against it)

```
./heston --config config/reference.cfg [options]
  --solver baseline|opt        (default baseline)
  --opt-level 0..6             (opt solver only; cumulative technique ladder, default 6 = all)
              7 | 8            (NEGATIVE CONTROLS: ctl-order, ctl-branch — §4b)
  --type call|put              (overrides config)
  --dump-every N               (snapshot every N steps → --dump-dir)
  --dump-dir results/run1      (default results/)
  --bench R                    (repeat solve R times, print per-rep timings CSV)
  --nt/--ns/--nv X             (override grid, for convergence & stability sweeps)
  --maturity Y                 (override option.maturity_years; benchmark knob —
                                dt = maturity/nt, so this is how a bigger grid is
                                held inside its stability bound without changing
                                the footprint being timed. Never for pricing.)
Output (stdout, one CSV line):
  price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec
```

The `solver` column self-describes: `baseline`, `opt-L0`..`opt-L6`, and
`opt-ctl-order` / `opt-ctl-branch` for the controls — deliberately NOT
`opt-L7/L8`, so a CSV row can never be misread as the ladder continuing up.
The benchmark scripts prefix two more columns, `sweep,build`, ahead of this
line; see `slurm/bench_common.sh`.

## 4. Serial optimisation menu (the 25% Optimisation story)

Aligned with the techniques the course actually teaches (L03 — use the
lecturer's names for them in the video). Baseline is written deliberately
straightforward (weights recomputed per cell, generic indexing). Optimised
version applies, in order, each measured separately:

1. **Loop-invariant code motion / hoisting + lookup table**: precompute the
   nine stencil weights per variance-row into a table before the time loop —
   hoists thousands of per-cell multiplications into a per-row lookup.
   (L03 techniques #7 hoisting + #9 lookup tables + #13 common
   subexpression elimination, all in one move.)
2. **Strength reduction**: no divisions in the hot loop — precompute
   reciprocals (1/ds², 1/dv², …) once; division throughput is 5–14 cycles
   vs 0–1 for multiply (L03 #15).
3. **Cache-friendly loop order**: outer j, inner i, contiguous walk matching
   the row-major layout — spatial/temporal locality, 64-byte cache lines
   (L03 #22; lecturer measured 8x from loop order alone).
4. **Hoist conditionals out of the loop**: interior loop runs `i=1..Ns-2`,
   `j=1..Nv-2` branch-free; boundaries in separate small loops (L03 #18/#20
   loop splitting).
5. **Induction-variable simplification**: replace per-cell `j*ns+i` index
   math with an incremented running index/pointer (L03 #14 — explicitly
   framed in lecture as the 2D→1D address-calculation case).
6. (Optional, measure) **loop unrolling** of the inner loop (L03 #17).

NOT doing: tiling/blocking (not course content), restrict tricks beyond
local raw pointers, fast-math (correctness story comes first). Precision
stays double for M1; float32 is noted as an M2 experiment (L03 #1 pairs it
with SIMD width).

### 4b. Ablation ladder — `--opt-level 0..6` (locked)

The "measured separately" promise is delivered as a **cumulative ladder**,
re-runnable by anyone (including the markers) with one sbatch sweep:

| Level | Adds technique | L03 name | Expectation |
|---|---|---|---|
| 0 | none — identical algorithm to BaselineSolver (sanity anchor) | — | null by design |
| 1 | per-row constants, S table, named spacing products, CSE | hoisting + lookup + CSE | real gain |
| 2 | no divisions in hot loop, S² table | strength reduction | real gain; first rung whose last bits move |
| 3 | contiguous traversal (outer j, inner i), row bases hoisted | loop order / locality | **null by construction** |
| 4 | branch-free interior, v=0 row peeled, bounds hoisted | loop splitting | **null by construction** |
| 5 | row pointers instead of `j*ns+i` per access | induction variable | measure |
| 6 | inner-loop unrolling ×4 | unrolling (optional) | measure |

**Negative controls (`--opt-level 7` and `8`, labelled `opt-ctl-order` and
`opt-ctl-branch` — never `opt-L7/L8`).** Levels 3 and 4 are null *because the
baseline was written correctly from the start*: it already loops outer-variance
/ inner-stock, and already peels the v=0 row into its own loop. A rung cannot
recover a loss that never happened. The honest way to show what those
techniques are worth is not to sabotage an earlier rung but to build the wrong
version and time it:

| Control | Built from | Changes exactly one thing | Gives meaning to |
|---|---|---|---|
| `ctl-order` | level 5 | loops swapped → `ns`-element stride | level 3 |
| `ctl-branch` | level 5 | v=0 row fused back in behind a per-cell `if` | level 4 |

Both are arithmetically identical to levels 2–6, so `test_opt_matches` holds
them to the same tolerance — a control that got a *different* answer would be
measuring a bug, not a technique. `ctl-order` is also the ground truth behind
the cache animation (`scripts/memory_anim.py --style cache`).

**Floating-point contract across the ladder.** Two claims, two checks:
- **Levels 0–1 are BIT-IDENTICAL to the baseline** — asserted as
  `max_abs_diff == 0.0` exactly, which is stronger than any epsilon. Level 1
  may give a subexpression a *name* but never change the order it was
  multiplied in (`0.5*v*S*S` groups as `((0.5*v)*S)*S`, so the hoisted form is
  `(half_v*S)*S`), and every true division stays a division.
- **Levels 2–8 agree to a MEASURED relative tolerance.** Level 2 replaces
  divisions with reciprocal multiplies, which genuinely changes the last bits.
  One absolute tolerance cannot cover price ≈ 196, delta ≈ 0.57, gamma ≈ 7.9e-4
  and vega ≈ 2190, so each is checked against its own magnitude. Measured worst
  case over 5200 steps: **3.3e-15** (call leg; exactly 0.0 on the put leg),
  tolerance set at 1e-11. Levels 2–8 also agree with *each other* bit-for-bit,
  which is what proves they change addressing and structure but not arithmetic.

Design rules:
- **Dispatch once, outside the time loop**: OptSolver picks a per-level
  `step_levelN()` kernel function up front. NEVER per-cell `if (level >= k)`
  checks — that would poison the very timings the ladder exists to measure.
- Kernels live in their own file (`solver_opt_kernels.cpp`) since each is a
  short (~40-line) variant; duplication is accepted — the ladder IS the
  ablation instrument, and diffing adjacent kernels shows exactly one
  technique.
- `test_opt_matches` loops over ALL levels 0–6 AND both controls, each vs
  baseline, under the two-tolerance contract above.
- Benchmark output: 7-bar ladder figure (cell-updates/sec per level)
  replaces the plain before/after bar. Baseline-vs-level-6 remains the
  headline number.
- **Honest null results**: if a rung shows no gain (compiler already did
  it — likely for levels 5/6 at -O2), report it as-is; the task sheet
  explicitly credits "techniques you tried that didn't work" under
  Optimisation (25%), and L03 says null results show methodology.

Profile with **gprof** (the course-taught tool: `-pg`, flat + call-graph
profile) plus manual `<chrono>` instrumentation; report where time goes
before and after. Evidence = cell-updates/sec ladder table on rangpur,
identical answers proven by `test_opt_matches` at every level.

## 5. Phases, order, acceptance criteria

Dependency chain: P1 → P2 → P3 → (P4, P5, P6 in any order) → P7 → P8.

**P1 — Plumbing.** `params` parser, `grid`, CLI skeleton, `io` result line.
✓ `./heston --config config/smoke.cfg` runs end-to-end printing zeros.

**P2 — Baseline solver.** Payoff init, interior stencil, all four boundaries,
time loop, price + Greeks readout.
✓ Price for `config/reference.cfg` is positive, finite and between intrinsic
  (0 for this call, since spot 5200 < strike 5250) and spot. **This criterion
  was recorded as passed at P2 but was actually FAILING**: reference.cfg then
  specified `nt = 2000`, which is 527× over the explicit stability bound, and
  the run printed `nan,nan,nan,nan`. Fixed in P7 by `nt = 980000` and
  `v_max = 0.64`; the config header carries the arithmetic and the measured
  justification for both. The lesson is in the criterion itself — "positive,
  finite" is only a real gate if somebody actually runs the config.
  **Measured 2026-08-09: 196.1683699**, delta 0.56717178, in 2295 s serial —
  against the semi-analytic Fourier price of 196.1692 that is a relative error
  of **4.2e-6**, and Monte Carlo over 400 000 real random paths brackets both
  at 196.1345 ± 0.5524. That is the strongest absolute-accuracy number the
  project has, and it exists only because the config was fixed.
✓ Weather-map dump of smoke run looks smooth (first visual debug).

**P3 — Validation gate.** `test_bs_collapse` (tolerance: rel err < 1e-3 on
reference-sized grid; record actual), `test_parity` (abs err < 1e-3·K·1e-3? —
measure, then fix tolerance and justify in video).
✓ Both tests green. **Nothing downstream starts until this passes.**

**P4 — Convergence study.** Script sweeps (Ns,Nv,nt) doublings via CLI
overrides; `convergence_plot.py` produces the deliverable figure.
✓ Price differences shrink at expected first-order-in-dt rate; plot saved.

**P5 — Instability demo.** `unstable.cfg` + snapshot dumps + weather map of
the checkerboard blow-up; log measured stable/unstable dt boundary vs the
solver's printed estimate.
✓ 30-second horror-show clip material exists.

**P6 — Weather map (deliverable quality).** `--dump-every` + `weather_map.py`
→ H.264 mp4 of value surface flowing backwards.
✓ Animation renders from C++ engine output (replaces Python prototype).

**P7 — Serial optimisation + rangpur benchmarking.** Implement §4 steps as
the §4b opt-level ladder, `test_opt_matches` green after each level lands,
then the `slurm/bench_*.sh` sweeps run solver=baseline plus opt levels 0–6
on a compute node across grid sizes. (An earlier single `bench_serial.sh`
was deleted: it could never have run inside the 15-minute QOS limit. The
five sweeps that replaced it are listed in §6 and RUNNING.md §3.)
✓ Ladder cell-updates/sec table + 7-bar figure from rangpur numbers
  (baseline vs level 6 = headline before/after).
✓ Flags, hostname, grid sizes, opt level all recorded in the CSV.

**P8 — Video assembly.** Script: problem (weather map) → method (spreadsheet
metaphor) → validation (tests) → convergence → optimisation story →
instability finale → reflection. 10 minutes, H.264.

## 6. Benchmark protocol (rubric: Benchmarking 25%)

- Only compute-node numbers count (`sbatch`, partition `cosc3500`).
- ≥5 reps per configuration; report median + min/max; `--bench R` handles it.
- **The node is pinned.** `--exclusive` (no co-tenant stealing memory
  bandwidth) and `--nodelist=r730-2`, because the partition mixes 24-core
  Xeon E5-2670 v3 physical machines with 8-core VMs, and sweeps measured on
  different hardware cannot be compared with each other.
- **The QOS caps wall time at 00:15:00** (`sacctmgr show qos` → cosc3500
  MaxWall). One job covering everything sits in the queue forever with reason
  `QOSMaxWallDurationPerJobLimit`, so the suite is **five sweeps** — A ladder,
  B scaling, C no-vectorise, D controls, E ladder-on-huge-pages — plus two
  one-off diagnostic jobs (`bench_peak.sh`, `diag_alloc.sh`, `diag_denormal.sh`),
  each writing its own CSV. `make remote-bench` submits the five; `make
  remote-diag` submits the diagnostics. **RUNNING.md is the single home for
  the operational detail**; this section states the protocol and why.
- Sweep grid sizes: 512×128, 1024×256, 2048×512 (reference), 4096×1024.
  **nt is scaled per grid so that both (a) dt ≤ 0.7·dt_stable, and (b) every
  grid does ~4e8 cell updates**, making the four points equal-work. The
  maturity moves with nt (`--maturity`), because dt = maturity/nt and the
  stability bound depends on the GRID, not on the contract:

  | grid | dt_stable (v_max 0.64) | nt | maturity | dt/dt_stable | cell updates |
  |---|---|---|---|---|---|
  | 512×128 | 5.875e-6 | 6000 | 0.0246 | 0.698 | 3.93e8 |
  | 1024×256 | 1.466e-6 | 1500 | 0.00153 | 0.696 | 3.93e8 |
  | 2048×512 | 3.661e-7 | 400 | 0.0001 | 0.683 | 4.19e8 |
  | 4096×1024 | 9.147e-8 | 100 | 0.0000064 | 0.700 | 4.19e8 |

  This is legitimate because cell-updates/sec is a RATE: a 400-step run
  measures the same throughput as a 400 000-step one. Verified rather than
  assumed — `slurm/diag_denormal.sh` showed throughput is unchanged when the
  grid is full of denormals, full of NaNs, or full of values around 1e33.
- Opt-level ladder (§4b): full 0–6 sweep plus `ctl-branch` at the reference
  grid size; other grid sizes need only baseline + level 6 (scaling story);
  `ctl-order` runs in the controls job at a shorter time loop because it is
  ~30× slower than the ladder top.
- Record: hostname, CPU model, compiler version, flags, date, job id and
  `lscpu` cache sizes → CSV header comments, reproduced on every figure.
- Course benchmarking rules (L03, follow verbatim):
  - **Consume the results** (print the price) — otherwise the compiler can
    delete the benchmarked loop entirely (lecturer demoed runtime → 0).
  - **Problem must not fit in cache** — reference grid is 2×8 MB buffers,
    comfortably exceeding L2/L3 per-core; state this in the video.
  - **Verify correctness while benchmarking** — `make test` runs before
    every bench job, wired into `build_and_validate()` in
    `slurm/bench_common.sh`, which exits rather than benchmark a build that
    has not just re-proved itself against the reference solver.
  - Timer resolution: solve takes seconds; `<chrono>` steady_clock is ample.
  - Don't be surprised by null results — compiler may have already done the
    optimisation; report honest null results, they show methodology.
- Amdahl/Gustafson framing saved for M2; M1 story = redundant-work elimination
  (strength reduction) plus two measured memory effects that are *latency*,
  not bandwidth: page size and traversal order.

**Two measurement traps found on this node, both now encoded in the scripts:**

1. **Page size, worth 1.66×.** `Grid` allocates two 8 MiB buffers per solve
   and frees them per rep. Rep 1 is served by `mmap` and, with
   `transparent_hugepage=[always]`, backed by 2 MiB pages — 16 MiB then needs
   ~8 TLB entries. glibc afterwards raises its dynamic `mmap` threshold past
   8 MiB, so reps 2–5 come from the reused heap on 4 KiB pages: 4096 entries,
   far past the core's L2 TLB, and every stencil access risks a page walk.
   Result: rep 1 at 1.35e8 cell-updates/s, reps 2–5 at 8.1e7, reproducible to
   three digits. Proved by pinning `MALLOC_MMAP_THRESHOLD_` in
   `slurm/diag_alloc.sh` — forcing `mmap` makes all five reps fast.
   **Consequence for reporting:** a single production solve is a fresh process
   and gets the fast regime, so sweep A's median understates real performance.
   Sweep E (`slurm/bench_pages.sh`) re-measures the whole ladder with the
   threshold pinned, and both are reported.
2. **The clock is NOT the explanation** — `slurm/bench_peak.sh` samples
   `/proc/cpuinfo` during an actual solve and finds 3100 MHz pinned
   throughout, so the rep-1 advantage is not turbo. That job also measures the
   machine's scalar double peak (5.117 GFLOP/s) so the roofline uses a
   measured ceiling rather than a datasheet figure.

## 6b. Milestone 2 technique alignment (locked now so M1 sets it up)

From the lectures, the exact M2 toolkit (nothing fancier):

- **SIMD = AVX intrinsics, `<immintrin.h>`** — NOT compiler auto-vectorisation
  (L04: "typically gives pretty subpar results") and NOT assembly (L04:
  discouraged). `__m256d` (4 doubles) on our data; optional float32 variant
  doubles the width to 8 (ties to L03 precision technique). Use
  `aligned_alloc` to 64-byte (cache-line) boundaries — exactly why the grid
  is a flat contiguous array with S contiguous. FMA (`_mm256_fmadd_pd`) for
  the stencil blend. Remainder loop for ns not divisible by width.
  Realistic expectation: L04 quotes ~4x rather than the full vector width,
  on the grounds that such kernels are usually memory-bound. **That reasoning
  does not apply here and the measurement says so** — see §6 and RESULTS.md,
  where level 6 reaches 82% of this node's measured scalar peak while asking
  for ~3 GB/s on a four-channel machine. This kernel is compute-bound, so
  widening the arithmetic is exactly the right next move and the ceiling is
  the vector width rather than the bus. Expect closer to the width than L04's
  rule of thumb, and measure it rather than assume either way.
- **OpenMP** — `#pragma omp parallel for` on the OUTER j loop (L05:
  "parallelise outermost loops first"), `default(none)` with explicit
  shared/private lists (L05 recommendation), thread counts via
  OMP_NUM_THREADS. Our design needs ZERO synchronisation inside a step
  (threads read shared immutable cur, write disjoint rows of next) — one
  implicit barrier per step at the end of the parallel for. Never
  `critical` in the hot loop (L05 demo: 40x SLOWDOWN). Expect speedup ≈
  physical cores, not hyper-threads; accept tiny FP differences from
  reordering (L05: non-associativity), bounded by test tolerance.
- **MPI/CUDA** — optional extensions only, first to cut (unchanged).

## 7. Risks & cut order

| Risk | Mitigation |
|---|---|
| v=0 boundary subtly wrong | BS-collapse test catches it (xi=0 exercises v-drift); compare against textbook Heston prices if doubt remains |
| ~~Explicit dt too small ⇒ reference nt=2000 unstable~~ **RESOLVED P7** | It did fire — reference.cfg printed `nan,nan,nan,nan`. Fixed by `v_max` 1.0→0.64 and `nt` 2000→980000; the solver's printed stability estimate is what diagnosed it. The instability became P5 content rather than a failure. |
| Rangpur queue congestion near deadline | Benchmarks are P7 but slurm scripts tested (debug_serial.sh) in P1 week |
| Scope creep | Cut order fixed by spec: historical study first, then extras; never solver/validation/benchmarks |

## 8. Definition of done (Milestone 1)

- [x] `make test` green locally and on rangpur (9 of 9 ladder levels + controls)
- [x] Convergence figure committed
- [x] Instability clip + measured boundary
- [x] Weather-map animation from C++ output
- [x] Reference config actually runs: 196.1684 vs semi-analytic 196.1692 (4.2e-6)
- [x] Rangpur opt-level ladder benchmark table + figure (headline baseline→L6
      = **2.56× huge pages / 1.81× default allocator**; both reported, and the
      1.66× page-size effect that separates them is its own result)
- [x] Negative controls measured, so the two null rungs have a number attached
- [x] Roofline + cache-traffic figures TESTING the "memory-bound" claim
      (it failed: level 6 reaches ~81% of measured scalar peak — compute-bound)
- [ ] 10-min H.264 video, submission format per rubric
- [ ] Repo clean: `make clean && make` works from fresh clone on rangpur
