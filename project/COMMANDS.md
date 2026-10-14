# Commands

This file lists every command used to build, test, run and benchmark the
pricer, locally and on rangpur. Commands are run from `project/`. Section 5 is the five minute
live demo for the interview.

## 1. Local build and tests

The build needs Homebrew `open-mpi` and `libomp` on the Mac. On rangpur the
same Makefile uses Open MPI over g++ 8.5.

```bash
make            # builds ./heston and the three test binaries
make test       # runs the tests, prints ALL TESTS PASS
make test-mpi   # the ladder test as 2 and 3 MPI processes
make avx-check  # AVX2 row against the scalar row, x86_64 build
make smoke      # one tiny solve to check the binary runs
make clean      # removes the binary, tests and object files
```

`make test` runs these in order.

| Binary | What it checks | Tolerance |
|---|---|---|
| `tests/test_bs_collapse` | with `xi = 0` and `v0 = theta` the baseline matches closed-form Black-Scholes | 2.5e-3 relative |
| `tests/test_parity` | put-call parity holds on the full parameter set | 1e-4 dollars |
| `tests/test_ladder_matches` | every rung matches the baseline on price, delta, gamma and vega, call and put, and `omp` matches `simd` exactly | 1e-11 relative, 0 for omp vs simd |

`make test` then runs `test_ladder_matches` again with `OMP_NUM_THREADS` set
to 1, 2 and 7.

`make avx-check` builds `tests/test_avx_row` with `-mavx2 -mfma` and runs it.
On the arm64 Mac it uses `clang++ -arch x86_64`, so the test is an Intel
program that runs under Rosetta. On rangpur it builds natively. It prints the
worst difference for five row lengths and fails if the AVX2 row was not built.
`make test` on the Mac still runs the `simd` rung, but there it uses the
scalar fallback and prints a warning, so only `make avx-check` and rangpur
test the real AVX2 code.

`make test` cannot check the cross term or the vol of vol, because the
collapse test switches them off and parity holds in every model. This target
prices the full parameter set against Fourier inversion and fails if they
disagree. It needs numpy, so it runs locally only, after any change to the
stencil or the readout.

```bash
make check-heston
```

If the project venv is broken, any python with numpy and scipy works.

```bash
make check-heston PY=/usr/bin/python3
```

## 2. Running the pricer

```bash
./heston --config config/bench.cfg --solver serial
```

This takes about 2 seconds on the Mac. `config/reference.cfg` is the
accurate 980000 step grid, about 1e12 cell updates, which takes about 3
minutes with `omp` on the Mac's 8 cores and about 16 minutes with `serial`.

It prints one CSV line, with the columns described below the flag table.
The command line overrides the config file.

| Flag | Meaning |
|---|---|
| `--config PATH` | required, the `key = value` config file |
| `--solver baseline\|serial\|simd\|omp\|mpi` | which rung to run, default `baseline` |
| `--solver omp_scalar\|omp_persistent\|omp_runtime\|omp_edges\|omp_float` | the experiment solvers, each the `omp` rung with one thing changed (CODE_GUIDE §9d and §9e) |
| `--type call\|put` | contract type |
| `--ns / --nv / --nt` | grid nodes in stock, variance and time |
| `--maturity Y` | years to expiry, shortened by benchmarks so a big grid stays stable |
| `--vmax V` | top of the variance axis, the weak scaling sweep grows it with the rows |
| `--bench R` | repeat the solve R times, one CSV line each |
| `--dump-every N` | write `snap_NNNNNN.csv` every N steps for the animations |
| `--dump-dir DIR` | where the snapshots go |

The `threads` and `ranks` columns say how many OpenMP threads and MPI
processes did the work. The full header is

```
price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec,threads,ranks,exchange_seconds,boundary_seconds
```

`boundary_seconds` is the part of `seconds` the omp and mpi rungs spent on
the serial boundary work.

`exchange_seconds` is the part of `seconds` that rank 0 of the `mpi` rung
spent on the edge row exchange. It is 0 for every other rung.

The stability and Feller line and any off-node warning go to stderr, so
stdout stays one clean CSV line.

### 2a. OpenMP settings

The `omp` rung reads these environment variables. They go in front of the
command, or are exported in a slurm script.

```bash
OMP_NUM_THREADS=8 ./heston --config config/bench.cfg --solver omp
OMP_NUM_THREADS=24 OMP_PLACES=cores OMP_PROC_BIND=close ./heston --config config/bench.cfg --solver omp
```

| Variable | Meaning |
|---|---|
| `OMP_NUM_THREADS=N` | team size. When it is unset, the runtime uses every CPU the process may run on, which under slurm is `--cpus-per-task` |
| `OMP_PLACES=cores` | a place is one physical core, so two threads never share one core's hyper-threads |
| `OMP_PROC_BIND=close` | threads go on consecutive places in the place list. On the r730 nodes the cores are numbered even on socket 0 and odd on socket 1, so consecutive threads alternate sockets (RESULTS M2.2) |
| `OMP_PROC_BIND=spread` | threads are spread evenly across both sockets |
| `OMP_PROC_BIND=false` | the operating system may move threads between cores |
| `OMP_DISPLAY_ENV=true` | the runtime prints the settings it actually used to stderr at start up |

The CSV `threads` column comes from `omp_get_max_threads()`, so it shows the
team size the program actually got.

### 2b. MPI runs

The `mpi` rung is started with `mpiexec`. Every other rung can be too, but
then every rank runs the whole solve on its own and only rank 0 prints.

```bash
# Mac, 2 processes with 4 threads each
OMP_NUM_THREADS=4 mpiexec --oversubscribe --bind-to none -n 2 ./heston --config config/bench.cfg --solver mpi

# rangpur, one process per node on two nodes, 24 threads each
mpiexec -n 2 --map-by node --bind-to none -x OMP_NUM_THREADS=24 ./heston --config config/bench.cfg --solver mpi

# rangpur, one process per socket, 12 cores each
mpiexec -n 4 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12 ./heston --config config/bench.cfg --solver mpi
```

| Flag | Meaning |
|---|---|
| `-n N` | start N processes (ranks) |
| `--map-by node` | place ranks round robin over the nodes, so 2 ranks on 2 nodes is one each |
| `--map-by ppr:1:socket:PE=12` | one rank per socket, each given 12 cores (processing elements) |
| `--bind-to none` | do not pin a rank to one core, so its threads can use every core |
| `--bind-to core` | pin each rank to the cores it was given by `--map-by` |
| `--report-bindings` | print which cores each rank was bound to, on stderr |
| `--oversubscribe` | allow more ranks than slots, needed on the Mac and inside a small slurm job |
| `-x VAR=value` | set an environment variable on every rank, including ranks on other nodes |

On rangpur the slurm scripts also export these two before `mpiexec`, as the
matrixass MPI runs needed. They stop Open MPI trying network transports that
fail to start on these nodes.

```bash
export OMPI_MCA_mtl=^ofi,psm,psm2,psm3
export OMPI_MCA_btl=^uct,ofi
```

The `mpi` rung refuses `--dump-every`, and refuses to run when there are fewer
than two variance rows per rank.

| Config | Purpose |
|---|---|
| `config/reference.cfg` | the real contract, gives a correct price |
| `config/bench.cfg` | timing only, maturity shortened so a rep takes seconds, same grid footprint as reference |
| `config/bench_itm.cfg` | timing at the largest grid, read out deep in the money |
| `config/demo.cfg` | the animation run |
| `config/smoke.cfg` | tiny grid, under a second |
| `config/unstable.cfg` | timestep over the stability limit, for the blow-up clip |

## 3. Rangpur

Off campus, connect the UQ VPN first. The Makefile uses the ssh alias
`rangpur` with connection sharing (see `CLAUDE.md`).

```bash
make sync           # rsync the source up, excludes .git, results and .env
make remote-debug   # sbatch slurm/debug.sh, build and test on a compute node
make remote-peak    # sbatch slurm/bench_peak.sh, the node's compute ceiling
make queue          # squeue for this user
make fetch          # rsync results and slurm-*.out back into results/
```

`slurm/debug.sh` asks for 8 cores, builds from clean, runs one smoke solve
per rung (`omp` on 8 threads, `mpi` as 2 ranks of 4 threads), runs
`make test`, `make test-mpi` and `make avx-check`, and counts AVX-512
(`zmm`) instructions in the binary, which must be 0, and AVX2 (`ymm`) instructions in
`stencil_avx.o`, which must not be. It does not request a specific node, so
its timings are not benchmarks.

These two commands check the AVX2 instructions by hand after a build on
rangpur.

```bash
objdump -d heston | grep -c zmm             # must print 0
objdump -d src/stencil_avx.o | grep vfmadd  # the FMA instructions
```

`slurm/debug_mpi.sh` asks for two nodes with two tasks each, builds once on
the shared home directory, and runs `test_ladder_matches` as 2 ranks (one
per node), 3 ranks and 4 ranks, then prints the host of each rank.

```bash
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/debug_mpi.sh'
```

`sbatch` and `squeue` are in `/opt/slurm/bin`, which a non-interactive ssh
does not put on `PATH`. The Makefile targets add it. This command submits
one script by hand.

```bash
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/debug.sh'
```

Every slurm script loads the compiler and Open MPI modules and puts the MPI
wrappers on `PATH`, because the module does not do that by itself.

```bash
module load compiler-rt/latest
module add mpi/openmpi-x86_64
export PATH=/usr/lib64/openmpi/bin:$PATH
```

The `cosc3500` QOS caps a job at 15 minutes, so each benchmark sweep is its
own job. Benchmark numbers only ever come from a compute node through
`sbatch`, never from the login node.

## 3a. Benchmark sweeps

Every sweep is its own job on an exclusive r730 node, because the QOS caps a
job at 15 minutes. Each one sources `slurm/bench_common.sh`, builds from
clean, runs `make test`, `make test-mpi` and `make avx-check`, and refuses to
time anything if a test fails. Results are written to
`results/bench_<sweep>_<jobid>.csv` with a provenance header, and stderr in
the matching `.err`.

```bash
make sync
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/bench_m2_ladder.sh'
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/bench_m2_strong.sh omp'
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/bench_m2_strong.sh omp_scalar'
make queue
make fetch
```

| Script | Node | What it runs |
|---|---|---|
| `bench_m2_probe.sh` | r730-2 | core layout, MPI binding check, one rep of each rung, used to size the others |
| `bench_m2_ladder.sh` | r730-2 | baseline, serial, simd, omp on 24 threads, mpi as 2 ranks of 12, at 2048 by 512 |
| `bench_m2_strong.sh SOLVER` | r730-2 | SOLVER (`omp` or `omp_scalar`) at 1 to 24 threads, 2048 by 512 |
| `bench_m2_weak.sh` | r730-2 | omp at 1 to 24 threads, 21 interior rows per thread |
| `bench_m2_size.sh` | r730-2 | serial, simd, omp on 24 threads, 512 by 128 up to 8192 by 1024 |
| `bench_m2_exp_a.sh` | r730-2 | E4 persistent region, E5 schedule, E8 parallel boundaries |
| `bench_m2_exp_c.sh` | r730-2 | E11 page regime, pinned against glibc's default, on 24 threads and serial |
| `bench_m2_exp_float.sh` | r730-2 | E3 float buffers against double, on 1 and 24 threads, and the float error on the reference grid |
| `bench_m2_denormal.sh` | r730-2 | E12, every rung with and without flushing denormals to zero |
| `bench_m2_optlevels.sh` | r730-2 | baseline and serial built at -O0 to -O3, one thread, each level checked for the -O2 answer first |
| `profile_m2.sh` | r730-2 | gprof flat profiles of baseline, serial and simd, written to `results/profile_<solver>_<jobid>.txt` |
| `bench_m2_exp_b.sh` | r730-2 | E6 thread binding, E7 NUMA page placement with `numactl --interleave=all` |
| `bench_m2_mpi.sh` | r730-2 and r730-3 | mpi as 1x24 and 2x12 on one node, 1x24 and 2x12 per node on two nodes |
| `bench_peak.sh` | r730-2 | scalar peak, AVX2 FMA peak on 1, 12, 24 cores, triad bandwidth on 1, 12, 24 cores |

`bench_common.sh` pins the page regime with `MALLOC_MMAP_THRESHOLD_=1048576`,
the Milestone 1 headline regime, and sets `OMP_PLACES=cores` and
`OMP_PROC_BIND=close` unless a script overrides them. Its `maturity_for ns
nv vmax nt` prints the maturity that puts `dt` at 70% of the stability
limit, so a sweep can change `nt` or the grid and stay stable. `run` takes an
optional `RUN_PREFIX`, used for `numactl --interleave=all` and for
`env -u MALLOC_MMAP_THRESHOLD_`.

`--vmax V` on the command line sets the top of the variance axis. The weak
sweep uses it to grow the grid in variance at a fixed spacing.

## 4. Figures and animations

These need the Python venv.

```bash
arch -arm64 /usr/bin/python3 -m venv .venv
arch -arm64 .venv/bin/python -m pip install -r scripts/requirements.txt
```

| Target | Produces |
|---|---|
| `make demo-anim` | `results/demo_timevalue.mp4`, `demo_surface.mp4`, the value map |
| `make collapse-anim` | `results/bs_collapse.mp4`, `xi_smile.mp4`, Heston collapsing onto Black-Scholes |
| `make mc-check` | a Monte Carlo price, a third check alongside the PDE and Fourier prices |
| `make convergence` | `results/convergence.png`, the grid refinement study |
| `make payoff` | `results/payoff.png` |
| `make m2-figures` | F1 to F14 and `summary.txt` in `results/m2/`, drawn from the final sweep CSVs the Makefile names by job id, and the grid and Black-Scholes figures from `results/collapse_work/dump` |
| `make code-shots` | the Milestone 2 code screenshots in `results/m2/code_shots/`, the data layout, one per rung, the test gate, the CSV header and the test output |

## 5. Interview demo runbook

The demo takes about five minutes. Connect to the UQ VPN and open the shared
`rangpur` ssh connection before the interview starts, so no command has to
wait for a password.

Start the compute node job first, because it queues and then takes a few
minutes. It builds from clean on a compute node, runs every test and the
AVX2 check, and counts AVX-512 instructions, which must be 0.

```bash
make remote-debug
make queue
```

While it runs, build and test on the Mac. `make test` takes about 20
seconds and ends with `ALL TESTS PASS`. It runs the ladder test at 1, 2
and 7 threads, and 7 does not divide the rows evenly.

```bash
make clean && make && make test
```

Run each rung on the same grid. The price, delta, gamma and vega columns
are the same in every line, and the `seconds` column decreases. On the Mac
`simd` prints a warning and runs the scalar row, because the Mac has no AVX2.

```bash
./heston --config config/bench.cfg --solver baseline
./heston --config config/bench.cfg --solver serial
./heston --config config/bench.cfg --solver simd
OMP_NUM_THREADS=8 ./heston --config config/bench.cfg --solver omp
OMP_NUM_THREADS=4 mpiexec --oversubscribe --bind-to none -n 2 ./heston --config config/bench.cfg --solver mpi
```

Show where a benchmark number comes from. The header names the node, the
job, the compiler, the flags and the binding, and the rows below it are
the five repetitions of each rung.

```bash
head -30 results/bench_m2_ladder_651295.csv
```

Fetch the debug job's output and show the test lines and the `zmm` count.

```bash
make fetch
ls -t results/logs | head -1
grep -E 'PASS|FAIL|zmm' results/logs/slurm-<jobid>.out
```

If the queue is slow, `results/logs/slurm-651295.out` is the full log of
the ladder job, with the build, every test and the AVX2 check on r730-2.

Have these files open in the editor, in this order. They are
`src/stencil.cpp` (`update_cells_scalar`), `src/stencil_avx.cpp` (the FMA
lines), `src/solver_omp.cpp` (the pragma), `src/solver_mpi.cpp`
(`exchange_edge_rows`) and `results/m2/f1_ladder.png`.
