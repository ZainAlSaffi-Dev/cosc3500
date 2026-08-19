#!/bin/bash -l
#
# Sweep B — scaling. Baseline vs level 6 across the four grid sizes PLAN §6
# fixes, so speedup can be plotted against problem size and the memory-bound
# claim can be checked as the working set grows past L2 (256 KB) and L3 (30 MB).
#
# Per-grid nt is chosen so every grid does ~4e8 cell updates, making the four
# points equal-work and directly comparable:
#   grid        dt_stable    nt     maturity    dt/dt_stable   cell updates
#   512x128     5.875e-6     6000   0.0246      0.698          3.93e8
#   1024x256    1.466e-6     1500   0.00153     0.696          3.93e8
#   2048x512    3.661e-7      400   0.0001      0.683          4.19e8
#   4096x1024   9.147e-8      100   0.0000064   0.700          4.19e8
#
# RE-RUN NOTE (supersedes job 556506). The first run of this sweep let the
# allocator do what it liked, and what it liked differed by row: the 64 MiB
# grid was above glibc's mmap ceiling so every rep got fresh huge-page-backed
# mappings, while the smaller grids fell back to the reused 4 KiB heap after
# rep 1. That made the curve confound cache scaling with page size — the
# whole subject of sweep E — and the 4096x1024 row read FASTER than 2048x512,
# which is backwards for a cache story. Two fixes here:
#   1. MALLOC_MMAP_THRESHOLD_ is pinned low (same trick as bench_pages.sh),
#      so every row that CAN take huge pages takes them on every rep. One
#      honest exception: 512x128's buffers are 0.5 MiB each, smaller than a
#      single 2 MiB huge page, so that row physically cannot join the
#      huge-page regime — but at 256 total 4 KiB pages it fits the STLB
#      (1024 entries on Haswell) outright, so no page-walk pressure exists
#      at that size and the regime label does not matter.
#   2. The 4096x1024 row reads out deep in the money (config/bench_itm.cfg,
#      spot moved, grid identical) because at maturity 6.4e-6 the standard
#      $50-out-of-the-money readout underflows to ~2e-16 — a printed zero
#      still defeats dead-code elimination, but a ~$5,150 price that baseline
#      and opt-L6 must agree on is a witness with digits in it.
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
#SBATCH --job-name=heston_scaling
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

SWEEP_NAME=scaling
source slurm/bench_common.sh

build_and_validate "$BASE_FLAGS"
write_header

# Pin the allocator so every rep of every row is served by fresh mmap (and
# huge pages where the buffer is big enough to hold one). 256 KiB, not the
# 1 MiB bench_pages.sh uses, because the smallest grid's buffers are 0.5 MiB
# and must still clear the threshold.
export MALLOC_MMAP_THRESHOLD_=262144
echo "# MALLOC_MMAP_THRESHOLD_=${MALLOC_MMAP_THRESHOLD_} (single page regime across rows and reps)" >> "$OUT"

echo "=== Sweep B: scaling, baseline vs level 6, pinned page regime ==="
for SPEC in "512 128 6000 0.0246 config/bench.cfg" \
            "1024 256 1500 0.00153 config/bench.cfg" \
            "2048 512 400 0.0001 config/bench.cfg" \
            "4096 1024 100 0.0000064 config/bench_itm.cfg"; do
  set -- $SPEC
  NS=$1; NV=$2; NT=$3; MAT=$4; CFG=$5
  run B O2-hugepage --config "$CFG" --ns "$NS" --nv "$NV" --nt "$NT" \
      --maturity "$MAT" --solver baseline
  run B O2-hugepage --config "$CFG" --ns "$NS" --nv "$NV" --nt "$NT" \
      --maturity "$MAT" --solver opt --opt-level 6
done

echo "done: $OUT"
wc -l "$OUT"
