#!/bin/bash -l
#
# Sweep F — the compiler-flag control. This job exists to TEST the headline
# claim of the whole optimisation section rather than assert it.
#
# The claim (RESULTS.md §2): "level 2 is the only rung on the ladder that the
# compiler is forbidden to perform, which is why it is the only one that pays."
# Every other sweep builds at -O2 and so cannot distinguish that explanation
# from "the other rungs happen not to help on this machine". Two rebuilds
# separate them decisively.
#
#   -O3            turns on the transformations -O2 leaves off, notably
#                  vectorisation and more aggressive unrolling. If the null
#                  rungs really were already done by the compiler at -O2, then
#                  -O3 should not move the ladder much either. If -O3 lifts the
#                  BASELINE towards level 6, that is the compiler doing by hand
#                  what we did by hand, which is a result worth reporting.
#
#   -O2 -ffast-math  is the decisive one. It licenses exactly the reciprocal
#                  substitution that level 2 performs, x/(2*ds) -> x*(1/(2*ds)),
#                  by dropping strict IEEE semantics. The prediction is sharp
#                  and falsifiable: the BASELINE at -ffast-math should jump to
#                  roughly level 2's throughput, and level 2's advantage over
#                  it should collapse to nothing. If that happens, the claim is
#                  demonstrated. If the baseline does NOT move, the claim as
#                  worded is too strong, and the honest finding becomes the
#                  more interesting one: the compiler was GIVEN permission and
#                  still did not take it, so hand optimisation beat it on
#                  merit rather than on legal grounds.
#
# LOCAL PILOT, and the reason this job is worth queueing. Run on the Mac
# (clang, x86-64 under Rosetta, so indicative only and never quotable), the
# baseline moved 4.57e8 -> 4.77e8 cell-updates/s under -ffast-math, about 4%,
# while level 2 moved 8.41e8 -> 9.14e8. The baseline did NOT close the gap.
# That points at the second outcome above rather than the first. It has to be
# confirmed on r730-2 with GCC 8.5 before it is claimed, because the local
# toolchain already disagrees with the Xeon about unrolling (level 6 is a 23%
# regression locally and a null there), and the whole point of this job is that
# the target machine is the only one that counts.
#
# NOTE ON CORRECTNESS. -ffast-math changes floating-point semantics, so the
# usual iron rule needs a word. build_and_validate still runs the full test
# suite on each build and still refuses to benchmark a failing one, but
# test_opt_matches demands bit-identical answers from levels 0 and 1, and
# -ffast-math can break that legitimately rather than because of a bug. If this
# job dies in validation on the fast-math leg, that is informative and should
# be reported as such, not worked around by loosening the test — the -O3 leg
# and the -O2 reference leg still stand on their own.
#
# These builds are CONTROLS, not the shipped configuration. Nothing in
# RESULTS.md is measured at -O3 or -ffast-math, and the project's own numbers
# stay at -O2 throughout (PLAN §1), because -ffast-math would invalidate the
# put-call parity and Black-Scholes-collapse tolerances the validation rests on.
#
# Expected runtime ~11 min: three builds at ~1.5 min each including validation,
# plus six short benchmark runs.
#
#SBATCH --job-name=heston_flags
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
# Same hygiene as every other benchmark job: exclusive so no co-tenant steals
# memory bandwidth, and pinned to r730-2 so this is comparable with sweeps A-E.
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00

module load compiler-rt/latest
set -u

SWEEP_NAME=flags
source slurm/bench_common.sh

O3_FLAGS="-std=c++17 -O3 -Wall -Wextra -Iinclude"
FAST_FLAGS="-std=c++17 -O2 -ffast-math -Wall -Wextra -Iinclude"

# The huge-page regime, matching sweep E, so these rows are comparable with the
# headline table rather than with the 4 KiB one.
export MALLOC_MMAP_THRESHOLD_=1048576
LADDER_ARGS="--config config/bench.cfg --nt 400 --maturity 0.0001"

build_and_validate "$BASE_FLAGS"
write_header
echo "# flags_o3='${O3_FLAGS}'" >> "$OUT"
echo "# flags_fastmath='${FAST_FLAGS}'" >> "$OUT"
echo "# MALLOC_MMAP_THRESHOLD_=${MALLOC_MMAP_THRESHOLD_} (huge-page regime, matches sweep E)" >> "$OUT"

echo "=== Sweep F leg 1: -O2 reference, same job so the comparison is fair ==="
run F O2 $LADDER_ARGS --solver baseline
run F O2 $LADDER_ARGS --solver opt --opt-level 2
run F O2 $LADDER_ARGS --solver opt --opt-level 6

echo "=== Sweep F leg 2: -O3 ==="
build_and_validate "$O3_FLAGS"
run F O3 $LADDER_ARGS --solver baseline
run F O3 $LADDER_ARGS --solver opt --opt-level 2
run F O3 $LADDER_ARGS --solver opt --opt-level 6

echo "=== Sweep F leg 3: -O2 -ffast-math (the decisive one) ==="
echo "If validation fails here, REPORT THAT — -ffast-math may legitimately"
echo "break the bit-identical requirement on levels 0 and 1."
build_and_validate "$FAST_FLAGS"
run F O2-fastmath $LADDER_ARGS --solver baseline
run F O2-fastmath $LADDER_ARGS --solver opt --opt-level 2
run F O2-fastmath $LADDER_ARGS --solver opt --opt-level 6

# Leave the tree in the shipped configuration.
make clean
make all

echo "done: $OUT"
wc -l "$OUT"
