# RUNNING.md — how to build, validate and reproduce every number

Everything in this project can be rebuilt from source with `make` and a C++17
compiler. Every performance number in `results/` came from a rangpur **compute
node** via `sbatch`, never from the login node and never from a laptop. This
file is the recipe for reproducing each one.

Student: s4800836. Cluster: `rangpur.eait.uq.edu.au`, partition `cosc3500`,
benchmark node `r730-2` (Intel Xeon E5-2670 v3, Haswell, 24 cores, 30 MB L3).

---

## 1. Local (macOS or Linux) — build and validate

No external dependencies. C++17 only.

```bash
cd project
make            # builds ./heston and the three test binaries
make test       # the three validation gates; prints ALL TESTS PASS
make smoke      # one tiny solve, sanity check that the binary runs
make clean
```

`make test` runs, in order:

| binary | what it proves | tolerance |
|---|---|---|
| `tests/test_bs_collapse` | with `xi = 0` and `v0 = theta` the Heston solver reproduces the closed-form Black-Scholes price | absolute, in dollars |
| `tests/test_parity` | put-call parity `C - P = S e^{-qT} - K e^{-rT}` holds for the full parameter set, including the cross term | fractions of a cent |
| `tests/test_opt_matches` | every optimisation level 0-6 and both negative controls agree with the baseline solver | L0/L1 bit-identical, L2+ relative 1e-11 |

The last one is the gate that makes the benchmarks meaningful: no optimised
kernel is ever timed unless it has just proved it still gives the reference
answer.

### The check that `make test` cannot do

```bash
make check-heston    # ~40 s, needs the Python venv
```

`make test` contains no gate on the two terms that make this a Heston solver
rather than a Black-Scholes solver. `test_bs_collapse` sets `xi = 0`, which
zeroes both the cross term and the vol-of-vol diffusion. `test_parity` checks
an identity that holds in every model and is measurably blind to both, as
RESULTS.md §7 shows. `test_opt_matches` only compares the solver against
itself.

`make check-heston` prices the full parameter set against Fourier inversion and
fails if they disagree. Its second leg uses a deliberately off-node grid, so the
price it checks is the bilinearly interpolated one, which is the code path the
headline number actually goes through.

It is not part of `make test` because it needs numpy, and `make test` has to run
inside every benchmark job on rangpur where there is no virtual environment.
Run it locally after any change to the stencil or the readout.

```bash
make convergence     # ~6 min, the grid-refinement study and its figure
```

### Running the pricer by hand

```bash
./heston --config config/reference.cfg
```

Prints one CSV line:

```
price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec
```

Flags (`src/main.cpp`; the command line overrides the config file):

| flag | meaning |
|---|---|
| `--config PATH` | required, the `key = value` config file |
| `--solver baseline\|opt` | reference solver, or the optimised one |
| `--opt-level 0..6` | which rung of the optimisation ladder to run |
| `--opt-level 7` | negative control: loops swapped (`opt-ctl-order`) |
| `--opt-level 8` | negative control: `v = 0` row fused back in (`opt-ctl-branch`) |
| `--type call\|put` | contract type |
| `--ns / --nv / --nt` | grid nodes in stock, variance, time |
| `--maturity Y` | years to expiry; used by the benchmarks to shorten the run without changing the memory footprint |
| `--bench R` | repeat the whole solve R times, one CSV line per rep |
| `--dump-every N` | write `snap_NNNNNN.csv` every N steps, for the animations |
| `--dump-dir DIR` | where those snapshots go |

### Config files

| file | purpose |
|---|---|
| `config/reference.cfg` | the real contract. This is the one that produces a correct price. |
| `config/bench.cfg` | **timing only.** Maturity is shortened to about four and a half hours so a rep takes seconds instead of forty minutes. The grid footprint and access pattern are byte-for-byte identical to `reference.cfg`, which is all that cell-updates-per-second depends on. The price column in a benchmark CSV is not a meaningful option price and is not read as one. |
| `config/demo.cfg` | the animation run, sized so the surface moves visibly. |
| `config/smoke.cfg` | tiny grid, runs in under a second. |
| `config/unstable.cfg` | deliberately over-long timestep, used to film the instability blow-up. |

---

## 2. Rangpur — the benchmark loop

Off campus, connect the UQ VPN first. The Makefile assumes an ssh alias named
`rangpur` with connection multiplexing (see the repo `CLAUDE.md`).

```bash
make sync          # rsync the source tree up (excludes .git, results, .env)
make remote-bench  # sbatch the four benchmark jobs
make queue         # squeue -u $USER
make fetch         # rsync results/ and slurm-*.out back down
```

`sbatch` and `squeue` live in `/opt/slurm/bin`, which rangpur only puts on
`PATH` for an interactive login shell. A plain `ssh host cmd` runs
non-interactively and reports `sbatch: command not found`, so the Makefile
targets prepend that directory explicitly.

To submit one job by hand:

```bash
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/bench_ladder.sh'
```

### Why the work is split across several jobs

The `cosc3500` QOS caps wall time at **15 minutes** (`sacctmgr show qos` ->
MaxWall). One job covering the ladder, the scaling sweep, the controls and the
no-vectorise rebuild needs about forty minutes, so it would sit in the queue
forever with reason `QOSMaxWallDurationPerJobLimit`. Each script below is
sized to finish comfortably inside the limit and writes its own CSV.
`scripts/bench_plot.py` takes any number of CSVs and merges them on the sweep
column.

### Benchmark hygiene, shared by every job

Set in each `#SBATCH` header and in `slurm/bench_common.sh`:

- `--exclusive` — no other job shares the node. The kernel is memory-bound and
  divide-bound, so a co-tenant on another core changes our numbers even though
  it never touches our core.
- `--nodelist=r730-2` — the `cosc3500` partition mixes 24-core r730 physical
  machines with 8-core a100/vcpu VMs. Every benchmark pins to the *same*
  physical node, otherwise the ladder, the scaling sweep and the controls
  cannot be compared with one another at all.
- `build_and_validate()` runs `make clean && make && make test` and **exits
  before benchmarking if any test fails.** No unvalidated build is ever timed.
- `write_header()` stamps `host`, `date`, `jobid`, compiler version, exact
  compiler flags, rep count, CPU model and all four cache sizes into the CSV as
  `#` comment lines. That header is the provenance proof: `host` reads
  `r730-2.compute.eait.uq.edu.au`, which is a compute node, not the login node.
- `REPS=5` by default, and the summary reports the **median** with the observed
  min-max range rather than a best-of.

---

## 3. What each job produces

| script | sweep | what it measures | output |
|---|---|---|---|
| `slurm/bench_ladder.sh` | A | the ablation ladder, baseline through level 6, at 2048x512, default allocator | `results/bench_ladder_<jobid>.csv` |
| `slurm/bench_scaling.sh` | B | baseline vs level 6 across 512x128, 1024x256, 2048x512, 4096x1024 (1 MiB to 64 MiB working set); `MALLOC_MMAP_THRESHOLD_` pinned since job 556834, largest grid reads out deep in the money via `config/bench_itm.cfg` (the unpinned first run is retired to `results/superseded/`) | `results/bench_scaling_<jobid>.csv` |
| `slurm/bench_novec.sh` | C | the headline pair rebuilt with `-O2 -fno-tree-vectorize` — **kept as a record, but vacuous as a control**: GCC 8.5 has no vectoriser on at `-O2` to disable (RESULTS.md §6); the valid control is sweep G | `results/bench_novec_<jobid>.csv` |
| `slurm/bench_controls.sh` | D | the two negative controls against their paired reference, in the same job at the same `nt` | `results/bench_controls_<jobid>.csv` |
| `slurm/bench_pages.sh` | E | the full ladder again with `MALLOC_MMAP_THRESHOLD_` pinned so every rep gets huge-page-backed memory | `results/bench_pages_<jobid>.csv` |
| `slurm/bench_flags.sh` | F | the compiler-flag control: rebuilds at `-O3` and at `-O2 -ffast-math` to test the claim that level 2 is the only rung the compiler is forbidden to perform. Ran as job 556833; the `-ffast-math` leg was refused by validation (level 1 loses bit-identity), which RESULTS.md §2 reports as a finding | `results/bench_flags_<jobid>.csv` |
| `slurm/bench_novec_o3.sh` | G | the vectorisation control sweep C should have been: `-O3` vs `-O3 -fno-tree-vectorize` in one job, attributing the `-O3` baseline lift (ran as job 556839) | `results/bench_novec_o3_<jobid>.csv` |
| `slurm/bench_peak.sh` | — | the node's achievable scalar double FLOP rate, measured rather than quoted from a datasheet, for the roofline's ceiling | `results/peak_<jobid>.txt` |
| `slurm/diag_alloc.sh` | — | diagnostic: proves the rep-1-is-fast effect is glibc's dynamic mmap threshold and not the kernel | `results/alloc_<jobid>.txt` |
| `slurm/diag_denormal.sh` | — | diagnostic: proves the low absolute throughput is the hardware and the five divisions per cell, not denormal numbers | job stdout |
| `slurm/debug_serial.sh` | — | five-minute smoke job, submitted early to shake out module and toolchain problems | job stdout |

Job stdout lands in `slurm-<jobid>.out` in the remote project directory, and
`make fetch` copies it into `results/logs/`.

---

## 4. Figures and animations (local, after `make fetch`)

These need the Python virtual environment:

```bash
arch -arm64 /usr/bin/python3 -m venv .venv
arch -arm64 .venv/bin/python -m pip install -r scripts/requirements.txt
```

| target | produces | what it is for |
|---|---|---|
| `make bench-plot` | `results/bench_ladder.png`, `bench_ladder_hugepage.png`, `bench_scaling.png` | the optimisation and benchmarking figures |
| `make roofline` | `results/roofline.png` | tests the memory-bound claim; the ceiling is the measured 5.117 GFLOP/s from `bench_peak.sh`, not a datasheet figure |
| `make memory-anim` | `results/memory_cache.mp4`, `memory_buffers.mp4` | why loop order and the double buffer matter; the caption quotes the measured slowdown straight from the CSV so it cannot drift from the table |
| `make demo-anim` | `results/demo_timevalue.mp4`, `demo_surface.mp4` | the weather map, the value surface marching backwards from expiry |
| `make collapse-anim` | `results/bs_collapse.mp4`, `xi_smile.mp4` | validation as a picture: Heston collapsing onto Black-Scholes, and the volatility smile flattening as `xi` goes to zero |
| `make mc-check` | console output | an independent third opinion on the price from Monte Carlo paths |
| `make code-shots` | `results/code_shots/*.png` | slide-ready screenshots of the load-bearing code (L1-vs-L2 diff, the swapped-loop control, the iron rule, a CSV provenance header, the matching-test output); needs pygments + pillow |
| `make payoff` | `results/payoff.png` | the payoff diagram for the video's first Introduction beat |

`scripts/bench_plot.py` and `scripts/roofline.py` **exit non-zero rather than
plot laptop numbers**, so a figure in `results/` can only have come from a
cluster CSV.

---

## 5. Reproducing the published tables end to end

```bash
cd project
make sync
make remote-bench                              # sweeps A, B, C, D
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/bench_pages.sh'
ssh rangpur 'export PATH=/opt/slurm/bin:$PATH; cd cosc3500/project && sbatch slurm/bench_peak.sh'
make queue                                     # wait until empty
make fetch
make bench-plot roofline
```

`results/bench_summary.md` is the write-up of those CSVs, and its provenance
block lists the job id and timestamp behind every table in it.
