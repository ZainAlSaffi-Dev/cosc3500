# HANDOFF — next session pickup (updated 2026-08-09, end of P4–P7-scaffold)

## Kickoff prompt (paste this to start the next session)

> Read `project/HANDOFF.md`, then PLAN §4/§4b and STUDY_GUIDE §15.
> P1–P6 are DONE and committed; P7 is scaffolded and waiting on MY
> level-1 kernel (`step_level1` in `src/solver_opt_kernels.cpp` — the
> hint skeleton is in place). First: the P2 drill topics 2–5 plus the
> P3/P4/P5 checkpoint quizzes are still owed — run them. Then guide me
> (hints, review after — don't write it) through step_level1; after it's
> green in test_opt_matches, you scaffold levels 2–6 one commit each.

## Where the project is

- **P1–P6 complete and committed; P7 scaffolded.** `make test` fully
  green locally (BS collapse, parity, opt_matches with L0).
- **P3 gate:** collapse rel err 1.16e-3 (tol 2.5e-3), parity gap 5.4e-6
  (tol 1e-4), node-aligned 421×51 grid — STUDY_GUIDE §10 P3 notes.
- **P4 convergence:** prices 190.996 → 195.605 → 195.925 → 196.102 on
  the aligned ladder; observed order ≈ 1.2 in dt; Richardson true price
  ≈ 196.16. `results/convergence.png` + `.csv`; regenerate with
  `python3 scripts/convergence_plot.py`.
- **P5 instability:** smoke-grid printed bound dt_stable = 1.526e-5;
  measured cliff: finite at 1.092× bound, 2.8e69 at 1.100×, NaN at
  1.122× → estimate ~9% conservative (why: corner-cell bound vs
  region-occupying mode — STUDY_GUIDE §8 P5 notes). unstable.cfg now
  smoke-sized at 4× bound (old comment claimed 4× but was ~2000× over on
  the reference grid; also reference frames are 13 MB each, unfilmable).
- **P6 weather map:** `weather_map.py` renders mp4 (+gif) from snapshot
  dumps; author reviewed the first cut (Blues on white) as unreadable →
  now inferno on a dark stage, plus `--time-value` (V − payoff: the real
  heat-transfer shot, flame at the strike pinching to a point at v=0)
  and `--style surface` (3-D sheet, slow camera drift). `--diverging`
  RdBu on white for the blow-up clip. Fixed colour scale, 98th-pct cap,
  NaN magenta-on-dark/black-on-white, `--max-frames` caps the NaN tail.
  Rendered: `results/weather_map.{mp4,gif}`, `weather_map_timevalue.mp4`
  (best single clip), `surface.mp4`, `surface_timevalue.mp4`,
  `blowup.mp4`. Snapshot dirs: results/smooth, results/blowup.
  STUDY_GUIDE §6 P6 rendering notes explain all four choices.
- **P7 scaffold:** `--opt-level 0..6` CLI → `Config::opt_level`;
  OptSolver mirrors baseline, dispatches a `StepKernel` function pointer
  ONCE before the time loop; kernels own interior + v=0 row, solve()
  owns shared boundaries; CSV label `opt-L<k>`. `step_level0` = exact
  baseline copy, measured 0.0e+00 diff, same throughput. Levels 1–6
  return nullptr → solver throws at the edge → test SKIPs honestly
  (all-SKIP = FAIL, so the test can't go green by doing nothing).

## THE author-writes piece next: step_level1 (hoisting ladder rung)

Working agreement §1d: the author types `step_level1` themself. Hint
skeleton is in `src/solver_opt_kernels.cpp` (hoist v-dependent scalars
per row; S/S² lookup tables O(ns) per step; named reciprocal spacings;
arithmetic must stay 1e-12-identical). After it lands and matches,
assistant scaffolds levels 2–6, one technique per commit, diffable.

## Open learning debt (quiz backlog)

- P2 drill: topic 1 (payoff→current_, buffer ping-pong) re-taught this
  session — author followed the two-whiteboard walkthrough but hasn't
  been re-quizzed to completion; topics 2–5 untouched (two-buffer
  rationale, v=0 upwind, where instability erupts + why, Feller).
  Questions + answers: STUDY_GUIDE Part V "P2 checkpoint drill".
- New quiz material since: P3 (erfc, node alignment, discrete parity,
  measured tolerances), P4 (order ~1.2, why it oscillates around 1),
  P5 (why the printed bound is ~9% conservative), P7 scaffold
  (function-pointer dispatch, why dispatch-once, why level 0 exists).

## Open decision (needed at P7 benchmarking): reference.cfg stability

Unchanged: reference.cfg nt=2000 is ~500× over its stability bound
(dt_stable ≈ 2.37e-7 → stable full solve ≈ 1.3M steps ≈ 45–60 min/rep).
Options: (A) shorten maturity so dt ≈ 0.8× bound at benchmarkable nt —
honest timings, correctness proven by P3 tests (recommended); (B) full
stable solve, 1 rep; (C) shrink grid (deviates from spec §8). Decide
before writing `slurm/bench_serial.sh` sweeps.

## Remaining to Milestone 1 done

1. Quiz backlog (above).
2. Author writes step_level1 → assistant levels 2–6 → ladder green.
3. reference.cfg decision → bench_serial.sh 0–6 sweep on rangpur
   (≥5 reps, compute node only) → bench_plot.py 7-bar ladder figure.
4. P8: 10-min video assembly (weather map → method → validation →
   convergence → opt ladder → instability finale → reflection).

## Local viz environment note

Plot/animation scripts need `pip install -r scripts/requirements.txt`.
This Mac's user site-packages currently mixes arm64 and x86_64 wheels
(pre-existing arm64 numpy + an x86_64 matplotlib installed 2026-08-09
under a Rosetta shell) — this session sidestepped it with a scratchpad
`pip --target` install; a native-terminal `pip3 install --user
--force-reinstall matplotlib` will fix it properly if plots fail.

## Standing repo rules (also CLAUDE.md — do not relearn)

- Develop/validate locally; rangpur is benchmark-only via sbatch
  (login-node numbers invalid). Sync with rsync, not git.
- No SIMD in M1. No parallel version benchmarked until it matches
  serial reference (test_opt_matches enforces this).
- OOP-shaped modern C++ per PLAN §1b; descriptive names; beginner-first
  narration comments; author's own comments are theirs — never edit.
- Author writes §1d core pieces (step_level1 is the last one);
  assistant scaffolds the rest, explains before code, quizzes at phase
  boundaries, keeps STUDY_GUIDE in sync.
- Commits authored by user alone — no Claude trailers/footers.
- `.env` holds rangpur/PDF passwords — NEVER print/echo/commit.
