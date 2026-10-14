# HANDOFF, next session pickup (updated 2026-10-14, Milestone 2)

## Milestone 2 state

P0 to P8 of `PLAN.md` §10 are done, and so are the P9 plan and script. All five
rungs (baseline, serial, simd, omp, mpi) exist, match the baseline, and the
parallel ones match simd bit for bit. Every sweep has run on r730-2 (and
r730-3 for two nodes) with denormals flushed and with unsplit AVX2 loads, in
the Phase 13 rerun (jobs 659238 to 659249 and 659436). `make m2-figures`
draws F1 to F18, and `make code-shots` draws the code images into
`results/m2/code_shots/`. The numbers are in `docs/RESULTS.md` under
Milestone 2, the decisions in `CHANGELOG.md` Phases 1 to 13, the concepts in
`STUDY_GUIDE.md` §15 to §18 and the code in `CODE_GUIDE.md` §8.10, §9a to
§9e. `COMMANDS.md` §5 is the interview demo runbook.

The current state of the work is listed below.

- The shipped build keeps unaligned 256 bit loads and stores unsplit
  (`NOSPLIT` in the Makefile, experiment E2b), and every M2 sweep was rerun
  with it in Phase 13. The ladder is 57, 159, 456, 6367 and 6984 Mcell/s,
  123x on one node, and two nodes reach 8310, 147x.
- E1 (jobs 659248 and 659249) and E9 (job 659243) were rerun with the flags.
  Still open and untested are why `-march=native` slows the auto-vectorised
  serial rung, why the fast-math serial rung is slower with unsplit loads, and
  whether raising `OMPI_MCA_btl_vader_eager_limit` lets E9 overlap on the
  2 x 12 layouts.
- The small grid's two-node layouts vary between jobs, 8310 in job 659242
  and 9364 in the E9 job 659243. The headline quotes 659242, and appendix
  A22 of the deck shows both.
- omp on 24 threads at 2896 x 724 fell to about half its earlier rate in
  both Phase 13 size jobs (3288, and 3267 in the control job 659437 with the
  split back on). It is the
  2^24 address effect of RESULTS M2.8, and the flags are not the cause.
- `presentation/INTERVIEW.md` is the interview preparation, with the format,
  a friend's account of it, the plain-English answer, the numbers to know,
  the unusual C++ with an answer for each, a code tour and the questions to
  expect.
- P9 is the recording. The deck (`presentation/deck_m2.pdf`, 35 recorded
  frames, the divider and 32 backup slides, 68 pages) and the script
  (`presentation/M2_SCRIPT.md`, 9:43 at 140 words a minute) carry the Phase
  13 numbers. Every rung that changed the code has a before and after, and
  frame 34 is the reflection. The bottom-right corner of every recorded
  frame is kept empty for the webcam. Only the recording is left.
- `HPC_JUSTIFICATIONS.md` holds the reason for every HPC choice. Update it
  whenever a decision changes.
- P10 is the interview rehearsal and the final pass of both guides,
  including the older M1 sections that still have em dashes. CHANGELOG and
  PLAN still use the old word for a block of rows in places.
- P11 is the CUDA bonus, to be done only after everything else.
- Milestone 2 is committed up to the E1 and E9 setup (fa00544). Everything
  since, including Phase 13, the deck rebuild and INTERVIEW.md, is not
  committed yet. Commits are authored by the user alone, with no co-author
  trailer.

These notes cover building on the WSL machine.

- `make m2-figures`, `make code-shots` and `make slides-m2` work here. `PY`
  selects `.venv/bin/python` on any machine other than the Mac, and
  `slides-m2` falls back to pdflatex twice because tectonic is not
  installed.
- `results/` is git-ignored, so the figures exist only on disk and are
  included in the zip that `make package` writes.
- WSL has no MPI compiler, so `make` cannot build `heston` here. A conda Open
  MPI in a scratch folder built it once for the snapshots that F14 and F18
  read (`results/collapse_work` and `results/demo`). Installing
  `libopenmpi-dev` would let `make collapse-anim` and `make demo-snapshots`
  run directly.
- `project/CLAUDE.md` is local only (git-ignored) and holds the writing,
  code and workflow rules.

These rules were learned this session.

- Never run `make sync` while a job is running unless the job builds in its
  own work copy. Every slurm script now does (`~/cosc3500/work/job_<id>`),
  and `sync` excludes `slurm-*.out` and the test binaries.
- Every benchmark loads `slurm/flush_denormals.c` with `LD_PRELOAD`. Without
  it rates move by up to 42% with `nt` and grid size.
- On r730 nodes `OMP_PROC_BIND=close` alternates sockets (the
  even-numbered CPUs are on socket 0).
- `perf stat` does not count properly inside a slurm job on rangpur.
- Validation takes about 2.5 of a job's 14 minutes. If it grows, check the
  CSV header's `# date` against the job start before blaming the benchmark.
- The rangpur login node is an AMD EPYC 7542, so `-march=native` there does
  not match the r730 Xeons. Use `-march=haswell` for login node disassembly.
- The docs use the terms edge row exchange, copied rows and exchange share
  for the MPI exchange (CHANGELOG Phase 12), block of rows for a rank's
  share of the grid, and buffer for each of the two arrays.
  `project/CLAUDE.md` lists the full set of terms.
- `docs/E1_E9_FINDINGS.md` collects the first E1 and E9 results. Its numbers
  predate Phase 13, and RESULTS M2.6 has the current ones.
- The project venv needed rebuilding (`.venv/bin/python` was missing).

---

# Milestone 1 handoff (kept for history)


**P1 to P7 are done.** What remains for Milestone 1 is **P8 (the video)** and
the **learning debt** at the bottom of this file. Everything P7 produced is
listed below with the numbers, so the video script can be written directly
from here.

## Kickoff prompt (paste this to start the next session)

> Read `project/docs/HANDOFF.md` end to end, then `STUDY_GUIDE.md` §12b, §13 and
> §14 and `PLAN.md` §4b/§6. P1 to P7 are done, and what is left is P8 (the
> 10-minute video) and the quiz backlog.
>
> Start by confirming `make test` passes and that
> `./heston --config config/demo.cfg` still prints `196.1054098`. That number
> is the regression check for the interpolation, and every figure downstream
> assumes it.
>
> The standing rules still apply. Benchmark numbers come only from a rangpur
> compute node via `sbatch`. The credentials are in `project/.env`, and they
> must never be printed, echoed or committed. `make test` must pass after
> every change. Python runs as `arch -arm64 .venv/bin/python`. Commits are
> authored by the user alone.

---

## What P7 produced

The full results narrative is now in **[RESULTS.md](RESULTS.md)**, which is
tracked and carries the provenance beside each number. The rangpur operating
detail moved to **[RUNNING.md](RUNNING.md)** and the video plan to
**[PRESENTATION.md](../presentation/PRESENTATION.md)**. This file keeps only
what is recorded nowhere else, which is the pickup prompt above, the
outstanding work below and the learning debt.

For orientation, the reference price is 196.1683699, which agrees with
Fourier inversion to 4.2e-6. The serial speedup from the baseline to level 6
is 2.56x on huge pages and 1.81x on the default allocator, and almost all of
it comes from the single rung the compiler is forbidden to perform.

## What is left

### P8, the video (the only remaining deliverable)

The slide-by-slide plan, time budget and submission checklist are in
**[PRESENTATION.md](../presentation/PRESENTATION.md)**. Do not repeat the
video outline here.

### Open items found in the 2026-08-10 audit (status updated 2026-08-11)

They are ordered by what they cost. The first is the only one that is
pass/fail.

1. **Commit the work. This is still open and still pass/fail.** The P7
   sources, all slurm scripts, the five deliverable docs, the presentation layer
   (`presentation/`, `scripts/code_shots.py`, `scripts/payoff_diagram.py`)
   and the 2026-08-11 fixes are untracked or modified. A fresh clone of HEAD
   builds a pre-P7 program with one ladder level and a test suite that
   passes without testing anything. The Code Submission gate is pass/fail,
   and a zip built from git would fail it.
2. ~~Sweep B mixes memory regimes~~ **Closed 2026-08-11.** It was re-run
   pinned as job 556834 with a deep-in-the-money readout for the largest
   grid (`config/bench_itm.cfg`), and the old CSV was retired to
   `results/superseded/`. RESULTS.md §3 has the new table and the remaining
   physical caveat (buffers smaller than one huge page cannot be backed by
   transparent huge pages).
3. **The stability boundary has no saved output.** The nt=15000 / 14900 /
   14600 bracket exists only as prose in STUDY_GUIDE §P5 and in a commit
   message. Every other quantitative claim has a file with a provenance
   header behind it.
4. ~~`bench_flags.sh` has not been run~~ **Closed 2026-08-11.** It ran as job
   556833. The `-O2` leg replicated sweep E to three digits. `-O3` raised
   the baseline ×1.92, all of it from the vectoriser according to the new
   sweep G (job 556839, `slurm/bench_novec_o3.sh`), while L2/L6 did not
   change. The `-ffast-math` leg was stopped by the validation rule because
   level 1 lost bit-identity. RESULTS.md §2/§6/§10 reports this and records
   the old sweep C as a vacuous control (GCC 8.5 has no vectoriser at `-O2`
   to disable).
5. **The refinement ladder co-varies ns, nv and nt**, so no single factor is
   isolated. This is defensible for a convergence study and worth saying out
   loud, and RESULTS.md §7 and §10 now say it.
6. **`dt_stable` omits the advective and cross terms.** It is conservative in
   practice, but the printed number leaves out part of the CFL condition
   (recorded in RESULTS.md §10).

### Learning debt (quiz backlog)

These quizzes are owed to the author, with the questions and answers in
STUDY_GUIDE Part V.

- **P2 drill.** Topic 1 (payoff → `current_`, and the two buffers swapping
  pointers each step) was taught again but never quizzed again to
  completion. Topics 2 to 5 have not been covered (why there are two
  buffers, the upwind difference at v=0, where instability starts and why,
  and the Feller condition).
- **P3.** The topics are erfc, node alignment, discrete parity and measured
  tolerances.
- **P4 (rewritten 2026-08-10).** Why does the refinement ladder have to put
  the strike on a node as well as the spot? Why is an external anchor better
  than the solver's own finest grid? What does it mean when the observed
  order falls from ~2.2 to ~0.4 across the ladder? (The error is reaching a
  floor that refinement cannot reduce.)
- **P5.** Why is the printed bound ~9% conservative?
- **P7 (new).** Why must levels 0 and 1 be bit-identical when level 2 need
  not be? Why is it written `(half_v*S)*S` instead of `half_v*(S*S)`? Why
  are levels 3 and 4 null by construction, and what do the controls
  measure? Why is the baseline division-bound when level 6 is not? What
  changed when the bottleneck moved?
- **§5b, the most valuable set and the most likely interview question.**
  Where did the randomness go? This covers Feynman-Kac, why there is no
  random number generator, which PDE term each source of noise became, why
  the coefficient is `v` instead of `sqrt(v)`, and why the cross term has
  no ½.

---

## Environment

The Mac's user site-packages mix arm64 and x86_64 wheels and this shell runs
under Rosetta, so bare `python3 -c "import numpy"` fails. A pinned arm64 venv
is at `project/.venv` (git-ignored), and every script runs as
`arch -arm64 .venv/bin/python ...`, the same command the Makefile's `PY`
uses.
Rebuild with `arch -arm64 /usr/bin/python3 -m venv .venv && arch -arm64
.venv/bin/python -m pip install -r scripts/requirements.txt`.

The local shell's `g++` targets `x86_64` under Rosetta, so local timings are
sanity checks only and are never quoted. They also disagree with rangpur
about unrolling. Level 6 is a 23% regression locally and makes no difference
on the Xeon, from the same source code, so only timings on the target
machine are used.

## Standing rules (also in CLAUDE.md)

- Develop and validate locally. Rangpur is for benchmarks only, through
  `sbatch`, and numbers from the login node are invalid. Sync the code to
  rangpur with rsync.
- Milestone 1 uses no SIMD. No parallel version is benchmarked until it
  matches the serial version (`test_opt_matches` enforces this).
- The C++ is modern and object-oriented as PLAN §1b describes, with
  descriptive names and narration comments written for a beginner. **The
  author's own comments are theirs, and they are never edited.**
- PLAN §1d's "author writes level 1" was explicitly waived for the P7 session
  and the waiver is recorded in PLAN §1d. Everything else on that list stands.
- Commits are authored by the user alone, with no Claude trailers or footers.
- `.env` holds the rangpur and PDF passwords. Never print, echo or commit
  them.
