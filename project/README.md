# Heston PDE Option Pricer — COSC3500 Milestone 1 (serial)

A serial C++17 solver that prices a European option under the Heston
stochastic-volatility model by marching a finite-difference grid backwards in
time from expiry, plus the validation, benchmarking and optimisation study
around it.

```bash
make          # builds ./heston and the three test binaries, no dependencies
make test     # the validation gates; prints ALL TESTS PASS
./heston --config config/reference.cfg
```

---

## Where things are

| document | what it is for |
|---|---|
| **[WALKTHROUGH.md](docs/WALKTHROUGH.md)** | **Start here.** The whole program explained in plain English, from the problem down to the individual loops. No finance or C++ background assumed. |
| [RUNNING.md](docs/RUNNING.md) | How to build, test, benchmark on rangpur, and reproduce every number. The single home for the slurm and cluster mechanics. |
| [RESULTS.md](docs/RESULTS.md) | The measured findings, with the provenance behind each number. Tracked, unlike `results/`. |
| [PRESENTATION.md](presentation/PRESENTATION.md) | The Milestone 1 video plan: rubric, time budget, beat sheet, submission checklist. |
| [presentation/](presentation/) | The recording side of that plan. A timestamped shot table, notes on how each asset is shot, and the slide deck the video is recorded over. What to say stays in PRESENTATION.md. |
| [PLAN.md](PLAN.md) | The design contract: numerical scheme, code style rules, optimisation ladder design, benchmark protocol. |
| [STUDY_GUIDE.md](STUDY_GUIDE.md) | The long-form explainer, Parts I to V, from the finance through to the parallel future. Interview preparation. |
| [InitialSpecDesign.md](docs/InitialSpecDesign.md) | The original submitted brief. Frozen deliverable, deliberately not updated to match later decisions. |
| [HANDOFF.md](docs/HANDOFF.md) | Session pickup notes and the outstanding learning backlog. |

To avoid the same fact drifting apart in several places, each topic has one
owner: the scheme and style live in PLAN.md, cluster operations in RUNNING.md,
measured numbers in RESULTS.md, explanation in WALKTHROUGH.md and
STUDY_GUIDE.md. Everything else links rather than restates.

---

## Layout

```
include/   headers: params, grid, solver interface, kernels, io, black_scholes
src/       the solver: baseline, optimised ladder, grid, io, main
docs/      the written record: walkthrough, running, results, spec, handoff
presentation/  the video: plan, shot list, timeline, slide deck
tests/     three validation binaries, run by `make test`
config/    reference.cfg (the real contract), bench.cfg (timing only),
           demo.cfg (animations), smoke.cfg, unstable.cfg
slurm/     the sbatch jobs; five benchmark sweeps and two diagnostics
scripts/   figures and animations (Python, local only, never on rangpur)
results/   gitignored: CSVs, figures, snapshots, animations
```

One binary, `heston`. Three test binaries reusing the same objects. No external
C++ dependencies; the Python scripts are local-only and listed in
`scripts/requirements.txt`.

---

## What the project demonstrates

- A correct explicit finite-difference solver for the Heston PDE, including the
  `v = 0` transport boundary that most implementations get wrong by imposing a
  value instead of solving the reduced equation.
- Four independent correctness checks: collapse to Black-Scholes, put-call
  parity, grid refinement, and two external prices from Fourier inversion and
  Monte Carlo.
- A seven-rung optimisation ablation in which each rung adds exactly one
  technique, plus two negative controls that measure what the rungs which came
  out null would have been worth.
- Benchmarking on a rangpur compute node with the node reserved exclusively and
  pinned, provenance stamped into every result file, and a build that refuses to
  be benchmarked until it has re-proved itself against the reference solver.

The headline is not the two-and-a-half-times speedup. It is that only one of
the six techniques did anything, and it was the only one the compiler is
forbidden to perform on your behalf.
