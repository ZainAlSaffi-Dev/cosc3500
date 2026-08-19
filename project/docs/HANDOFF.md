# HANDOFF — next session pickup (updated 2026-08-09)

**P1–P7 are done.** What remains for Milestone 1 is **P8 (the video)** and the
**learning debt** at the bottom of this file. Everything P7 produced is listed
below with the numbers, so the video script can be written straight from here.

## Kickoff prompt (paste this to start the next session)

> Read `project/docs/HANDOFF.md` end to end, then `STUDY_GUIDE.md` §14b, §15 and
> §16 and `PLAN.md` §4b/§6. P1–P7 are done; what is left is P8 (the 10-minute
> video) and the quiz backlog.
>
> Start by confirming `make test` is green and that
> `./heston --config config/demo.cfg` still prints `196.1054098` — that number
> is the interpolation regression gate and every figure downstream assumes it.
>
> Standing rules still apply: benchmark numbers only from a rangpur compute
> node via `sbatch` (credentials in `project/.env` — never print, echo or
> commit them); `make test` green after every change; run python as
> `arch -arm64 .venv/bin/python`; commits authored by the user alone.

---

## What P7 produced

The full results narrative now lives in **[RESULTS.md](RESULTS.md)**, which is
tracked and carries the provenance beside each number. The rangpur operating
detail moved to **[RUNNING.md](RUNNING.md)** and the video plan to
**[PRESENTATION.md](../presentation/PRESENTATION.md)**. This file keeps only what is not
recorded anywhere else: the pickup prompt above, the outstanding work below and
the learning debt.

Headline, for orientation: reference price 196.1683699, agreeing with Fourier
inversion to 4.2e-6; serial speedup baseline to level 6 of 2.56x on huge pages
and 1.81x on the default allocator, essentially all of it from the single rung
the compiler is forbidden to perform.

## What is left

### P8 — the video (the only remaining deliverable)

The beat sheet, time budget and submission checklist are in
**[PRESENTATION.md](../presentation/PRESENTATION.md)**. Do not rewrite the spine here.

### Open items found in the 2026-08-10 audit (status updated 2026-08-11)

Ordered by what they cost. The first is the only one that is pass/fail.

1. **Commit the work — still open, still pass/fail.** The P7 sources, all
   slurm scripts, the five deliverable docs, the presentation layer
   (`presentation/`, `scripts/code_shots.py`, `scripts/payoff_diagram.py`)
   and the 2026-08-11 fixes are untracked or modified. A fresh clone of HEAD
   builds a pre-P7 program with one ladder level and a vacuously-passing
   test suite. The Code Submission gate is pass/fail and a zip built from
   git would fail it.
2. ~~Sweep B mixes memory regimes~~ **Closed 2026-08-11.** Re-run pinned as
   job 556834 with a deep-in-the-money readout for the largest grid
   (`config/bench_itm.cfg`); old CSV retired to `results/superseded/`.
   RESULTS.md §3 carries the new table and the remaining physical caveat
   (buffers smaller than one huge page cannot be THP-backed).
3. **The stability boundary has no artifact.** The nt=15000 / 14900 / 14600
   bracket exists only as prose in STUDY_GUIDE §P5 and in a commit message.
   Every other quantitative claim has a stamped file behind it.
4. ~~`bench_flags.sh` has not been run~~ **Closed 2026-08-11.** Ran as job
   556833. The `-O2` leg replicated sweep E to three digits; `-O3` lifted
   the baseline ×1.92 (all of it the vectoriser, per the new sweep G, job
   556839, `slurm/bench_novec_o3.sh`) while L2/L6 did not move; the
   `-ffast-math` leg was refused by the iron rule because level 1 lost
   bit-identity — reported as such in RESULTS.md §2/§6/§10, and the old
   sweep C is recorded there as a vacuous control (GCC 8.5 has no
   vectoriser at `-O2` to disable).
5. **The refinement ladder co-varies ns, nv and nt**, so no single factor is
   isolated. Defensible for a convergence study, worth saying out loud
   (now said in RESULTS.md §7 and §10).
6. **`dt_stable` omits the advective and cross terms.** It is conservative in
   practice, but the printed number is not the full CFL condition (recorded
   in RESULTS.md §10).

### Learning debt (quiz backlog)

Owed to the author, questions and answers in STUDY_GUIDE Part V:

- **P2 drill:** topic 1 (payoff → `current_`, buffer ping-pong) was re-taught
  but never re-quizzed to completion; topics 2–5 untouched (two-buffer
  rationale, v=0 upwind, where instability erupts and why, Feller).
- **P3:** erfc, node alignment, discrete parity, measured tolerances.
- **P4 (rewritten 2026-08-10):** why the ladder has to put the STRIKE on a
  node and not just the spot; why an external anchor beats using your own
  finest grid; what it means when the observed order falls from ~2.2 to ~0.4
  across the ladder (the error is landing on a floor, not converging faster
  then slower).
- **P5:** why the printed bound is ~9% conservative.
- **P7 (new):** why levels 0–1 must be bit-identical but level 2 need not be;
  why `(half_v*S)*S` and not `half_v*(S*S)`; why levels 3 and 4 are null by
  construction and what the controls measure; why the baseline is
  division-bound and level 6 is not; what changed when the bottleneck moved.
- **§5b — the highest-value set, most likely interview question:** where the
  randomness went. Feynman–Kac, why there is no RNG, which PDE term each piece
  of noise became, why the coefficient is `v` and not `sqrt(v)`, why the
  cross-term has no ½.

---

## Environment

The Mac's user site-packages mix arm64 and x86_64 wheels and this shell runs
under Rosetta, so bare `python3 -c "import numpy"` fails. A pinned arm64 venv
lives at `project/.venv` (git-ignored); every script runs as
`arch -arm64 .venv/bin/python ...`, which is what the Makefile's `PY` does.
Rebuild with `arch -arm64 /usr/bin/python3 -m venv .venv && arch -arm64
.venv/bin/python -m pip install -r scripts/requirements.txt`.

Note that the local shell's `g++` targets `x86_64` under Rosetta, so local
timings are sanity checks only and are never quotable. They also disagree with
rangpur about unrolling — level 6 is a 23% regression locally and a null on
the Xeon. Same source, opposite sign; the target machine is the only one that
counts.

## Standing rules (also CLAUDE.md — do not relearn)

- Develop and validate locally; rangpur is benchmark-only via `sbatch`
  (login-node numbers are invalid). Sync with rsync, not git.
- No SIMD in M1. No parallel version is benchmarked until it matches serial
  (`test_opt_matches` enforces this).
- OOP-shaped modern C++ per PLAN §1b; descriptive names; beginner-first
  narration comments; **the author's own comments are theirs — never edit them.**
- PLAN §1d's "author writes level 1" was explicitly waived for the P7 session
  and the waiver is recorded in PLAN §1d. Everything else on that list stands.
- Commits authored by the user alone — no Claude trailers or footers.
- `.env` holds rangpur/PDF passwords — never print, echo or commit them.
