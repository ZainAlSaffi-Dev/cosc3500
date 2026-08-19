#!/bin/bash -l
#
# Sweep A — the §4b ablation ladder at the reference footprint. This is the
# job the Optimisation mark hangs on: nine bars from baseline to level 6, all
# at 2048x512 (two 8 MB buffers, far larger than this node's 30 MB L3 once
# both are counted, and far larger than its 256 KB L2).
#
# ctl-order is NOT here: it runs ~30x slower than the ladder top and gets its
# own job at a shorter time loop (bench_controls.sh).
#
# SIZING NOTE. r730-2 is an Intel Xeon E5-2670 v3 (Haswell, 2014) clocked at
# a measured 3.1 GHz (bench_peak.sh sampled the clock during a real solve and
# found it pinned at 3100 MHz, so this is not throttling), and it turns this
# kernel over at ~4.7e7 cell-updates/s for the
# baseline — about a tenth of a modern laptop core. Measured, not assumed: a
# diagnostic job (slurm/diag_denormal.sh) confirmed the throughput is
# identical with denormals flushed, with the grid full of NaNs, and with
# values around 1e33, so it is the hardware and the five divisions per cell,
# not the data.
#
# The consequence is that nt has to be small to fit the 15-minute QOS limit.
# That is harmless: cell-updates/sec is a RATE, so a 400-step run measures the
# same throughput as a 400000-step one. Maturity is scaled with nt to keep
# dt/dt_stable at ~0.7, so every solve is stable and prints a finite price
# (PLAN §6: consume the results or the compiler may delete the loop).
#
# Expected runtime ~9 min, of which build + `make test` is ~5 (the
# validation tests solve a 421x101x56000 grid twice on each leg, and on this
# CPU that is not cheap — but the iron rule is the iron rule).
#
#SBATCH --job-name=heston_ladder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
# Benchmark hygiene, not politeness:
#   --exclusive  no other job shares this node. The kernel is memory- and
#                divide-bound, so a co-tenant on another core changes our
#                numbers even though it never touches our core.
#   --nodelist   the cosc3500 partition mixes 24-core r730 physical machines
#                with 8-core a100/vcpu VMs. Every benchmark job pins to the
#                SAME physical node (Xeon E5-2670 v3) so the ladder, the
#                scaling sweep, the controls and the no-vectorise control can
#                be compared with each other at all.
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00

module load compiler-rt/latest
set -u

SWEEP_NAME=ladder
source slurm/bench_common.sh

build_and_validate "$BASE_FLAGS"
write_header

# nt=400 at 2048x512 -> 4.2e8 cell updates per rep: ~9 s baseline, ~4 s opt.
LADDER_ARGS="--config config/bench.cfg --nt 400 --maturity 0.0001"

echo "=== Sweep A: ablation ladder at 2048x512 ==="
run A O2 $LADDER_ARGS --solver baseline
for LEVEL in 0 1 2 3 4 5 6; do
  run A O2 $LADDER_ARGS --solver opt --opt-level "$LEVEL"
done
# ctl-branch is cheap enough to keep beside the rung it ablates.
run A O2 $LADDER_ARGS --solver opt --opt-level 8

echo "done: $OUT"
wc -l "$OUT"
