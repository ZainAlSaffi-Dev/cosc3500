# HANDOFF — next session pickup (updated 2026-08-09, end of P3)

## Kickoff prompt (paste this to start the next session)

> Read `project/HANDOFF.md`, then `project/PLAN.md` (§1b, §1d, §5).
> P3 (validation gate) is GREEN and committed. FIRST: re-explain the
> five P2 drill topics in STUDY_GUIDE Part V ("P2 checkpoint drill")
> one at a time and re-quiz me until I can answer them — still owed
> from two sessions ago. Then P4 (convergence study) and P6 (weather
> map), both assistant-scaffolded per §1d.

## Where the project is

- **P1 + P2 + P3 done and committed.** P2 milestone commit `d70876e`;
  P3 (this session) on top.
- **P3 validation gate GREEN** (both tests, run locally):
  - `tests/test_bs_collapse`: xi=0, v0=theta vs closed form.
    Measured rel err 1.161e-3 (call), 1.116e-3 (put); tolerance 2.5e-3.
  - `tests/test_parity`: full Heston, C−P vs S·e^(−qT)−K·e^(−rT).
    Measured gap 5.4e-6 dollars; tolerance 1e-4.
  - Both use a **node-aligned** stable grid: ns=421 (spacing 50 → spot
    5200 = node 104), nv=51 (spacing 0.02 → v0 = node 2), nt=56000
    (dt ≈ 0.79× the printed stability bound). Alignment matters because
    `extract_result` reads the NEAREST cell, no interpolation — the
    smoke grid's snap alone is ~$7 of price (and v snaps 0.04→0.0476).
  - Rationale + measured numbers recorded in STUDY_GUIDE §10 "P3
    implementation notes".
- `make test` still fails at `test_opt_matches` — intentional, that's
  the P7 stub. Run the two P3 binaries directly until P7 lands.
- Plan's "rel err < 1e-3" was written for the reference-sized grid;
  1.16e-3 at spacing 50 is consistent with it — P4 convergence study
  will demonstrate the trend.

## Quiz status (§1d) — THE open learning debt (unchanged)

P2 checkpoint quiz: author could answer none of the five. Questions +
answers in STUDY_GUIDE Part V "P2 checkpoint drill". Re-explain each
conversationally, then re-quiz. Topics: payoff→current_ reasoning,
two-buffer rationale, v=0 forward/upwind difference, where instability
erupts first and why, Feller meaning/consequence. A P3 checkpoint quiz
(norm_cdf/erfc, node alignment, why parity is near-exact discretely,
measured-tolerance method) is also owed.

## Open decision (needed by P7, not before): reference.cfg stability

Unchanged from last session — see PLAN §5 P7 and previous discussion:
reference.cfg nt=2000 is ~500× over the stability bound; options are
(A) shortened-maturity honest benchmark (recommended), (B) 1-rep full
solve, (C) shrink grid. Also unstable.cfg nt is ~2000× over, not ~4×
as its comment claims — tune in P5 using the printed bound.

## Immediate next steps (in order)

1. Re-explain + re-quiz the five P2 drill topics; quiz P3 topics too.
2. P4 convergence study: sweep (ns, nv, nt) doublings via CLI
   overrides, `convergence_plot.py` → deliverable figure.
3. P5 instability demo (tune unstable.cfg nt from printed bound).
4. P6 weather-map animation (dump-every + weather_map.py → mp4).
5. P7 opt ladder — author writes level-1 kernel; reference.cfg
   benchmark decision due then.

## Standing repo rules (also CLAUDE.md — do not relearn)

- Develop/validate locally; rangpur is benchmark-only via sbatch
  (login-node numbers invalid). Sync with rsync, not git.
- No SIMD in M1. No parallel version benchmarked until it matches
  serial reference.
- OOP-shaped modern C++ per PLAN §1b: RAII, unique_ptr factory, no raw
  new/delete; descriptive names; beginner-first comments; narration
  comments mandatory; author's own comments are theirs — never edit.
- The author writes the §1d core pieces themself (one remains: level-1
  kernel); assistant scaffolds the rest, explains before code, quizzes
  at phase boundaries, keeps STUDY_GUIDE in sync.
- Commits authored by user alone — no Claude trailers/footers.
- `.env` holds rangpur/PDF passwords — NEVER print/echo/commit.
