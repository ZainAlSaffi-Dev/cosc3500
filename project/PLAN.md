# Milestone 1 Build Plan — Heston PDE Option Pricer (Serial)

Spec: `InitialSpecDesign.md`. Rubric: Intro/Background 20%, Optimisation 25%,
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
| Boundary v=0 | PDE degenerates: diffusion-in-v and cross terms vanish; solve reduced 1-D equation with **one-sided (upwind) dV/dv** since kappa*theta > 0 pushes inward | Standard Heston treatment |
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

Why this is safe: every byte is owned by exactly one vector or unique_ptr,
and scope exit frees it — leaks and double-frees are unrepresentable.

## 2. Repo layout (this scaffold)

```
project/
├── PLAN.md                  ← this file
├── Makefile                 ← build + test + sync + remote-bench targets
├── .gitignore               ← results/, binaries, frames
├── config/
│   ├── reference.cfg        ← spec §8 reference run (SPX 5250 call)
│   ├── smoke.cfg            ← tiny grid, seconds to run, for tests/CI
│   └── unstable.cfg         ← 4× stable dt — instability horror show
├── include/
│   ├── params.h             ← Option/Market/Heston/Grid/Config structs + parser decl
│   ├── grid.h               ← flat 2-buffer grid, idx(), swap
│   ├── solver.h             ← SolveResult, solver interface (baseline & opt)
│   ├── black_scholes.h      ← closed form + norm_cdf (validation only)
│   └── io.h                 ← grid snapshot dumps, result CSV line
├── src/
│   ├── main.cpp             ← CLI: --config --solver --dump-every --bench
│   ├── params.cpp
│   ├── grid.cpp
│   ├── solver_baseline.cpp  ← correct-first serial solver
│   ├── solver_opt.cpp       ← optimised serial (same interface, same answers)
│   ├── black_scholes.cpp
│   └── io.cpp
├── tests/
│   ├── test_bs_collapse.cpp ← xi=0, v0=theta ⇒ match Black-Scholes formula
│   ├── test_parity.cpp      ← C - P = S*exp(-qT) - K*exp(-rT)
│   └── test_opt_matches.cpp ← solver_opt ≡ solver_baseline to 1e-12
├── scripts/
│   ├── weather_map.py       ← snapshots → coloured frames → gif/mp4 (ffmpeg)
│   ├── convergence_plot.py  ← refinement study → convergence figure
│   └── bench_plot.py        ← bench CSVs → before/after bar + scaling figure
├── slurm/
│   ├── debug_serial.sh      ← 5-min sanity job (build + smoke run)
│   └── bench_serial.sh      ← the real benchmark job (baseline vs opt sweep)
└── results/                 ← gitignored; snapshots, CSVs, frames land here
```

One binary `heston`; tests are three tiny extra binaries reusing the same
objects. No external C++ dependencies. Python (matplotlib, ffmpeg) local only —
never needed on rangpur.

## 3. CLI contract (fixed now so scripts can be written against it)

```
./heston --config config/reference.cfg [options]
  --solver baseline|opt        (default baseline)
  --type call|put              (overrides config)
  --dump-every N               (snapshot every N steps → --dump-dir)
  --dump-dir results/run1      (default results/)
  --bench R                    (repeat solve R times, print per-rep timings CSV)
  --nt/--ns/--nv X             (override grid, for convergence & stability sweeps)
Output (stdout, one CSV line):
  price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec
```

## 4. Serial optimisation menu (the 25% Optimisation story)

Baseline is written deliberately straightforward (weights recomputed per cell,
generic indexing). Optimised version applies, in order, each measured
separately for the before/after narrative:

1. **Precompute stencil weights per row** (hoist all j-dependent and
   i-dependent coefficient work out of the inner loop).
2. **Loop order matches layout**: outer j, inner i, contiguous walk.
3. **Boundary cells out of the hot loop**: interior loop runs `i=1..Ns-2`,
   `j=1..Nv-2` branch-free; boundaries handled in separate small loops.
4. **Restrict-qualified raw pointers** for cur/next inside the kernel.
5. (Stretch) blocked/tiled traversal — measure, keep only if it wins.

Profile with `gprof` or `perf stat` locally; report where time goes before and
after. Evidence = cell-updates/sec table on rangpur, identical answers proven
by `test_opt_matches`.

## 5. Phases, order, acceptance criteria

Dependency chain: P1 → P2 → P3 → (P4, P5, P6 in any order) → P7 → P8.

**P1 — Plumbing.** `params` parser, `grid`, CLI skeleton, `io` result line.
✓ `./heston --config config/smoke.cfg` runs end-to-end printing zeros.

**P2 — Baseline solver.** Payoff init, interior stencil, all four boundaries,
time loop, price + Greeks readout.
✓ Price for reference.cfg is positive, finite, between intrinsic and spot.
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

**P7 — Serial optimisation + rangpur benchmarking.** Implement §4 steps,
`test_opt_matches` green after each, then `slurm/bench_serial.sh` runs the
baseline-vs-opt sweep on a compute node across grid sizes.
✓ Before/after cell-updates/sec table + figure from rangpur numbers.
✓ Flags, hostname, grid sizes all recorded in the CSV.

**P8 — Video assembly.** Script: problem (weather map) → method (spreadsheet
metaphor) → validation (tests) → convergence → optimisation story →
instability finale → reflection. 10 minutes, H.264.

## 6. Benchmark protocol (rubric: Benchmarking 25%)

- Only compute-node numbers count (`sbatch`, partition `cosc3500`).
- ≥5 reps per configuration; report median + min/max; `--bench R` handles it.
- Sweep grid sizes: 512×128, 1024×256, 2048×512 (reference), 4096×1024,
  nt scaled to keep the scheme stable.
- Record: hostname, compiler version, flags, date → CSV header comment.
- Amdahl/Gustafson framing saved for M2; M1 story = memory-access + redundant
  work elimination.

## 7. Risks & cut order

| Risk | Mitigation |
|---|---|
| v=0 boundary subtly wrong | BS-collapse test catches it (xi=0 exercises v-drift); compare against textbook Heston prices if doubt remains |
| Explicit dt too small ⇒ reference nt=2000 unstable | Solver prints stability estimate; adjust nt or shrink v_max; instability is content, not failure |
| Rangpur queue congestion near deadline | Benchmarks are P7 but slurm scripts tested (debug_serial.sh) in P1 week |
| Scope creep | Cut order fixed by spec: historical study first, then extras; never solver/validation/benchmarks |

## 8. Definition of done (Milestone 1)

- [ ] `make test` green locally and on rangpur
- [ ] Convergence figure committed
- [ ] Instability clip + measured boundary
- [ ] Weather-map animation from C++ output
- [ ] Rangpur before/after benchmark table + figure
- [ ] 10-min H.264 video, submission format per rubric
- [ ] Repo clean: `make clean && make` works from fresh clone on rangpur
