#!/bin/bash -l
#
# Sweep G — the vectorisation control that sweep C should have been.
#
# Sweep C paired -O2 against -O2 -fno-tree-vectorize and measured 1.00x, which
# turned out to be a tautology: GCC 8.5 does not enable the tree vectoriser at
# -O2 at all (that default arrived in GCC 12), so the flag disabled a pass that
# was already off. RESULTS.md §6 records the correction.
#
# Sweep F then showed the fact that makes the valid control worth running: at
# -O3 the BASELINE nearly doubles (5.25e7 -> 1.01e8 cell-updates/s) while the
# hand-optimised kernels do not move. This job asks WHICH -O3 pass did that.
# Pairing -O3 against -O3 -fno-tree-vectorize inside one job isolates the
# vectoriser's own share:
#   - if the -fno-tree-vectorize leg falls back to ~5.3e7, the whole -O3 lift
#     was the vectoriser packing the five divisions two-per-divpd, and the
#     ladder's story gains a sharp footnote: the compiler CAN halve the
#     division cost, but only by vectorising — which Milestone 1 excludes as
#     parallelism — and it still lands 31% short of scalar level 2;
#   - if it stays near 1.01e8, the lift was other -O3 passes and the
#     vectoriser is as irrelevant at -O3 as it was absent at -O2.
# Either answer completes §6. Both legs run baseline, L2 and L6 so the
# optimised kernels' indifference to -O3 is measured in the same job too.
#
# Huge-page regime pinned, matching sweeps E and F, so all three sweeps'
# numbers sit on one footing.
#
# Expected runtime ~13 min: two builds with validation (~5 min each) plus six
# short benchmark runs. Tight against the 14:00 limit but inside it; sweep F
# fitted three builds in the same window.
#
#SBATCH --job-name=heston_novec_o3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
# Same hygiene as every other benchmark job: exclusive so no co-tenant steals
# memory bandwidth, and pinned to r730-2 so this is comparable with sweeps A-F.
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00

module load compiler-rt/latest
set -u

SWEEP_NAME=novec_o3
source slurm/bench_common.sh

O3_FLAGS="-std=c++17 -O3 -Wall -Wextra -Iinclude"
O3_NOVEC_FLAGS="-std=c++17 -O3 -fno-tree-vectorize -Wall -Wextra -Iinclude"

export MALLOC_MMAP_THRESHOLD_=1048576
LADDER_ARGS="--config config/bench.cfg --nt 400 --maturity 0.0001"

build_and_validate "$O3_FLAGS"
write_header
echo "# flags_o3='${O3_FLAGS}'" >> "$OUT"
echo "# flags_o3_novec='${O3_NOVEC_FLAGS}'" >> "$OUT"
echo "# MALLOC_MMAP_THRESHOLD_=${MALLOC_MMAP_THRESHOLD_} (huge-page regime, matches sweeps E and F)" >> "$OUT"

echo "=== Sweep G leg 1: -O3, vectoriser on ==="
run G O3 $LADDER_ARGS --solver baseline
run G O3 $LADDER_ARGS --solver opt --opt-level 2
run G O3 $LADDER_ARGS --solver opt --opt-level 6

echo "=== Sweep G leg 2: -O3 -fno-tree-vectorize ==="
build_and_validate "$O3_NOVEC_FLAGS"
run G O3-novec $LADDER_ARGS --solver baseline
run G O3-novec $LADDER_ARGS --solver opt --opt-level 2
run G O3-novec $LADDER_ARGS --solver opt --opt-level 6

# Leave the tree in the shipped configuration.
make clean
make all

echo "done: $OUT"
wc -l "$OUT"
