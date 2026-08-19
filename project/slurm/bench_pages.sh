#!/bin/bash -l
#
# Sweep E — the ablation ladder re-measured with huge pages, which is the
# regime a REAL run actually experiences.
#
# Why this job exists. slurm/diag_alloc.sh established that on this node:
#   rep 1 of every configuration runs at ~1.35e8 cell-updates/s;
#   reps 2-5 run at ~8.1e7 — a 66% drop, reproducible to three digits.
# The cause is not the CPU (the clock is pinned at 3100 MHz throughout, see
# slurm/bench_peak.sh) and not the cache. It is the allocator. Grid allocates
# two 8 MiB buffers per solve and frees them at the end of each rep:
#   - the FIRST allocation is served by mmap, and with
#     /sys/kernel/mm/transparent_hugepage/enabled = [always] it is backed by
#     2 MiB transparent huge pages -> 16 MiB needs ~8 TLB entries;
#   - glibc then RAISES its dynamic mmap threshold past 8 MiB, so later reps
#     are served from the reused main heap on 4 KiB pages -> 4096 entries for
#     the same 16 MiB, far past this core's L2 TLB, and every stencil access
#     pays a page walk.
# Setting MALLOC_MMAP_THRESHOLD_ below the buffer size pins the good regime:
# diag_alloc.sh showed all five reps then land at 1.34e8.
#
# A single production solve is a fresh process that allocates once, so it gets
# the huge-page regime. Sweep A's median-of-5 therefore UNDERSTATES real
# performance, and understates the speedup: 1.81x there versus 2.56x here.
# Both are reported; this one is the one a user would experience.
#
#SBATCH --job-name=heston_pages
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00

module load compiler-rt/latest
set -u

SWEEP_NAME=pages
source slurm/bench_common.sh

build_and_validate "$BASE_FLAGS"
write_header

# Pin the allocator into the huge-page regime for every rep, so all five reps
# measure the same thing instead of one good rep and four degraded ones.
export MALLOC_MMAP_THRESHOLD_=1048576
echo "# MALLOC_MMAP_THRESHOLD_=${MALLOC_MMAP_THRESHOLD_} (forces mmap + THP for the grid buffers)" >> "$OUT"

LADDER_ARGS="--config config/bench.cfg --nt 400 --maturity 0.0001"

echo "=== Sweep E: ablation ladder, huge-page regime ==="
run E O2-hugepage $LADDER_ARGS --solver baseline
for LEVEL in 0 1 2 3 4 5 6; do
  run E O2-hugepage $LADDER_ARGS --solver opt --opt-level "$LEVEL"
done
run E O2-hugepage $LADDER_ARGS --solver opt --opt-level 8   # ctl-branch
run E O2-hugepage $LADDER_ARGS --solver opt --opt-level 7   # ctl-order

echo "done: $OUT"
wc -l "$OUT"
