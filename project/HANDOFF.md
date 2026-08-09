# HANDOFF — next session pickup (updated 2026-08-09, end of P2)

## Kickoff prompt (paste this to start the next session)

> Read `project/HANDOFF.md`, then `project/PLAN.md` (§1b, §1c, §1d, §4b).
> P2 (baseline solver) is COMPLETE and validated; everything is
> uncommitted — commit first (author-only, NO Claude trailers). Then,
> BEFORE any new code: re-explain the five P2 drill topics in
> STUDY_GUIDE Part V ("P2 checkpoint drill") one at a time and re-quiz
> me until I can answer them — I couldn't at the end of last session.
> Only then start P3 (validation tests): you scaffold
> `black_scholes.cpp` + both tests per §1d and walk me through them.

## Where the project is

- **P1 + P2 done, working end-to-end, NOT committed.** Last commit is
  `47c3970` (plan/docs); ALL code since is uncommitted working tree.
  Commit it as the P2 milestone first thing.
- `./heston --config config/smoke.cfg` runs: price **205.3243** (BS
  ballpark ≈203 ✓, between intrinsic 0 and spot 5200 ✓), delta 0.559,
  gamma 0.00073, ~4.9e8 cell-updates/s local, 0.67 s. Feller prints
  VIOLATED (0.1200 vs 0.1225) as PLAN §1c predicts. Snapshot dumps are
  smooth (zero checkerboard flips; call value monotone in S except
  ~1e-7-dollar deep-OTM dips — harmless explicit-FD artifact, one video
  sentence, don't "fix").
- **Author wrote (per §1d)**: payoff init, interior 9-point stencil,
  v=0 degenerate row (had a V−V ≡ 0 bug both derivatives; author fixed
  after review). Remaining author piece: level-1 hoisting kernel (P7).
- **Assistant scaffolded**: grid plumbing, solve() skeleton
  (dt/stability/Feller/timing/boundaries S=0, Smax, vmax),
  extract_result, `load_config` parser, `main` two-pass CLI, `io`
  (snapshots + CSV line). `smoke.cfg` nt raised 500→20000 (500 was ~30×
  over the stability bound).

## Quiz status (§1d) — THE open learning debt

P2 checkpoint quiz asked; **author could answer none of the five**.
Questions + full answers are in STUDY_GUIDE Part V under "P2 checkpoint
drill". Next session: re-explain each topic conversationally (don't just
point at the guide), then re-quiz. Topics: payoff→current_ reasoning,
two-buffer rationale, v=0 forward/upwind difference, where instability
erupts first and why, Feller meaning/consequence.

## Decisions locked this session (all codified in PLAN §1b/§1d)

1. **Comment style** — new C++ constructs explained to a beginner in
   plain English; may end "similar to X in Python"; NEVER a labelled
   "Python analogy:" comment.
2. **Descriptive naming** — `num_stock_nodes` not `ns`, `stock_spacing`
   not `ds` ("spacing" not "step" — avoids timestep collision),
   `current()` not `cur()`, indices `stock_i`/`var_j`. Cfg/CLI keys stay
   short (parser maps). Heston maths symbols (kappa, xi, rho, …) exempt.
3. **Author's comments are theirs** — in author-written pieces, review
   code only; corrections go in conversation, never comment edits.

## Open decision (needed by P7, not before): reference.cfg stability

Solver prints dt_stable ≈ 2.368e-7 for the 2048×512 reference grid;
reference.cfg nt=2000 → dt=1.25e-4, ~500× over. A stable full solve
needs nt ≈ 1.3M ≈ 45–60 min/rep — kills the ≥5-rep benchmark protocol.
Options: **(A)** benchmark reference grid at nt=2000 with maturity
shortened so dt ≈ 0.8×dt_stable — finite numbers, honest timings,
correctness proven separately by P3 tests (recommended); **(B)** full
stable solves, 1 rep; **(C)** shrink the reference grid (deviates from
spec §8). Also: unstable.cfg is ~2000× over (comment claims ~4×) — tune
its nt in P5 using the printed bound.

## Immediate next steps (in order)

1. Commit everything (author-only, no trailers) as the P2 milestone.
2. Re-explain + re-quiz the five drill topics (see Quiz status).
3. P3 validation gate: implement `black_scholes.cpp` (norm_cdf +
   closed form), `test_bs_collapse` (xi=0, v0=theta ⇒ match BS; note
   xi=0 makes Feller "OK" and dt_stable larger), `test_parity`
   (C − P = S·e^(−qT) − K·e^(−rT)). Use smoke-sized STABLE grids
   (respect dt_stable). Record actual errors, then fix tolerances.
   **Nothing downstream starts until both tests are green.**
4. Then P4 (convergence), P5 (instability demo), P6 (weather map),
   P7 (opt ladder — author writes level-1 kernel; reference.cfg
   decision due), P8 (video).

## Standing repo rules (also CLAUDE.md — do not relearn)

- Develop/validate locally; rangpur is benchmark-only via sbatch
  (login-node numbers invalid). Sync with rsync, not git.
- No SIMD in M1. No parallel version benchmarked until it matches
  serial reference.
- OOP-shaped modern C++ per PLAN §1b: RAII, unique_ptr factory, no raw
  new/delete; descriptive names; beginner-first comments (rules 1–3
  above); narration comments mandatory.
- The author writes the §1d core pieces themself (one remains: level-1
  kernel); assistant scaffolds the rest, explains before code, quizzes
  at phase boundaries, keeps STUDY_GUIDE in sync.
- `.env` holds rangpur/PDF passwords — NEVER print/echo/commit.
