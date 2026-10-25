# Heston PDE option pricer, COSC3500 Milestone 2

This program prices a European call option under the Heston model, where the
share price and its volatility both move randomly. It solves the Heston
partial differential equation on a grid of share price by variance, stepping
back in time from expiry, where every value is known, to today. Each step
rewrites every cell from its nine neighbours in the previous step, and no cell
in a step depends on another cell in the same step. The five versions of the
solver use this independence, each adding one technique to the version
before it, and every version gives the same price as the first.

## The result

Million cell updates per second on the 2048 x 512 grid, one exclusive rangpur
r730 node (2 x 12 cores), each the median of 5 runs. The job id of every
number is in [docs/RESULTS.md](docs/RESULTS.md).

| Rung | What it adds | Rate | Over the rung below | Over the baseline |
|---|---|---|---|---|
| `baseline` | The original loop, with divisions inside it | 57 | | 1x |
| `serial` | Divisions turned into multiplications, worked out once per solve | 159 | 2.8x | 2.8x |
| `simd` | Four cells per AVX2 instruction, with fused multiply-add | 456 | 2.9x | 8.0x |
| `omp` | The rows of each step shared between 24 threads | 6367 | 14.0x | 112x |
| `mpi` | One process of 12 threads per socket, edge rows exchanged | 6984 | 1.10x | 123x |

Across two nodes, two processes of 12 threads on each node reach 8310,
which is 147 times the baseline and 1.21 times one node in the same job.
On the 4096 x 1024 grid, which does not fit in the L3 caches, the second
node gives 1.64 times. Every rung matches the baseline to within 3.3e-15
on price, delta, gamma and vega, and the OpenMP and MPI rungs match the
SIMD rung bit for bit.

## Build, test and run

The build needs an MPI C++ compiler (`mpicxx`) with OpenMP. On rangpur that
is Open MPI 4.1 over g++ 8.5, loaded as in `slurm/bench_common.sh`.

```bash
make            # builds ./heston and the test programs
make test       # every test, prints ALL TESTS PASS
make test-mpi   # the ladder test on 2 and 3 MPI processes
make avx-check  # the AVX2 row against the scalar row

./heston --config config/bench.cfg --solver serial
OMP_NUM_THREADS=24 ./heston --config config/bench.cfg --solver omp
OMP_NUM_THREADS=12 mpiexec -n 2 ./heston --config config/bench.cfg --solver mpi
```

Each run prints one CSV line with the price, the Greeks, the time and the
rate. [COMMANDS.md](COMMANDS.md) has every command, the rangpur workflow and
the benchmark sweeps.

## Where everything is

| Path | What it holds |
|---|---|
| `include/grid.h` | The grid, `index()` and the two buffers |
| `src/solver_baseline.cpp` | Rung 0, the original loop nest. Every rung is checked against it |
| `src/stencil.cpp` | The constants worked out once, the scalar row and the boundaries |
| `src/stencil_avx.cpp` | The AVX2 row, the only file built with the AVX2 flags |
| `src/solver_serial.cpp`, `solver_simd.cpp`, `solver_omp.cpp`, `solver_mpi.cpp` | Rungs 1 to 4, one time loop each |
| `src/solver.cpp` | `make_solver`, the stability report and reading the price off the grid |
| `src/solver_experiments.cpp`, `solver_align.cpp`, `solver_float.cpp` | The experiments, each the OpenMP rung with one change |
| `tests/` | Four test programs, run before every benchmark |
| `slurm/` | One job script per benchmark sweep, sharing `bench_common.sh` |
| `scripts/` | Every figure, drawn from the CSVs the jobs wrote |
| `presentation/` | The video deck `deck_m2.html` (from `deck_m2.qmd`) with backup slides after the AI declaration, and the script `M2_SCRIPT.md` |
| `docs/RESULTS.md` | Every measured number, with its job id |

## Documents

| Document | What it is for |
|---|---|
| [docs/RESULTS.md](docs/RESULTS.md) | The measured findings, with the provenance behind each number |
| [CHANGELOG.md](CHANGELOG.md) | Every Milestone 2 change in order, with the reason for each decision |
| [COMMANDS.md](COMMANDS.md) | Every command to build, test, run and benchmark |
| [HPC_JUSTIFICATIONS.md](HPC_JUSTIFICATIONS.md) | Why each HPC technique was used, with the number behind it |
| [CODE_GUIDE.md](CODE_GUIDE.md) | The code, file by file, with what each block does and why |
| [STUDY_GUIDE.md](STUDY_GUIDE.md) | The maths and ideas behind the code |
| [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md) | The whole program in plain English |

## What the project shows

- A correct explicit finite-difference solver for the Heston PDE, checked
  against the Black-Scholes closed form, put-call parity and Fourier
  inversion.
- A ladder of five rungs, each adding one technique, each checked against
  the baseline's answer before it is timed.
- Every benchmark run on an exclusive compute node through `sbatch`, after
  the tests pass in the same job, with the node, compiler, flags and binding
  written into the header of every result file.
- The measured ceilings behind the ladder. The serial rung runs at the
  scalar peak, the SIMD rung below the AVX2 ceiling, and the threaded rung
  at the memory roof for grids larger than the L3 caches.
- The experiments that were tried and measured, including the ones that were
  rejected, in [docs/RESULTS.md](docs/RESULTS.md) M2.6.
