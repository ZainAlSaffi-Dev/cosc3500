#!/bin/bash -l
#
# Sweep D — the negative controls, in their own job because ctl-order is slow
# enough to blow the QOS limit on its own.
#
# Levels 3 (traversal order) and 4 (loop splitting) are null BY CONSTRUCTION:
# BaselineSolver was already written with the right loop order and with the
# v=0 row already peeled out, so those rungs have nothing to recover. Rather
# than manufacture a gain by sabotaging an earlier rung, this job measures
# what those techniques are WORTH by building the wrong version and timing it:
#   ctl-order   — level 5 with the loops swapped (ns-element stride)
#   ctl-branch  — level 5 with the v=0 row fused back in behind a per-cell if
#
# Everything runs at nt=100 on the same node in the same job, so each control
# has a PAIRED reference measured under identical conditions.
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
#SBATCH --job-name=heston_controls
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

SWEEP_NAME=controls
source slurm/bench_common.sh

build_and_validate "$BASE_FLAGS"
write_header

# nt=100 -> 1.05e8 cell updates: ~2 s for a ladder rung, ~35 s for ctl-order.
SHORT="--config config/bench.cfg --nt 100 --maturity 0.000025"

echo "=== Sweep D: negative controls vs their paired reference ==="
run D O2 $SHORT --solver baseline
run D O2 $SHORT --solver opt --opt-level 5
run D O2 $SHORT --solver opt --opt-level 6
run D O2 $SHORT --solver opt --opt-level 8   # ctl-branch
# Slowest item in the whole suite, and that cost IS the result.
run D O2 $SHORT --solver opt --opt-level 7   # ctl-order

echo "done: $OUT"
wc -l "$OUT"
