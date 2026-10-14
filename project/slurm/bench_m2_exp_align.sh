#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Experiment E2, where the rows sit in their cache lines, and aligned loads against unaligned.
# omp_align_loadu runs the shipped row with column 1 moved 0 to 7 doubles past a 64 byte boundary.
# Shifts 0 and 4 align the centre, 1 and 5 the west, 3 and 7 the east window, and 2 and 6 none.
# Shifts 4 to 7 repeat 0 to 3, so they show the noise.
# omp_align_load swaps in aligned loads and stores, omp_align_shuffle builds east and west in registers.
#
# The cosc3500 QOS allows 15 minutes, so the sweep runs as two jobs.
#   sbatch slurm/bench_m2_exp_align.sh one    one thread
#   sbatch slurm/bench_m2_exp_align.sh many   24 threads
#
#SBATCH --job-name=heston_m2_exp_align
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
PART=${1:-}
if [ "$PART" != "one" ] && [ "$PART" != "many" ]; then
  echo "usage: sbatch slurm/bench_m2_exp_align.sh one|many"
  exit 1
fi
SWEEP_NAME=m2_exp_align_${PART}
# E2 measured the build that split unaligned loads, so these turn the split back on. They come after the shipped no-split flags, so they win
export AVXEXTRA="-mavx256-split-unaligned-load -mavx256-split-unaligned-store"
source slurm/bench_common.sh
build_and_validate
write_header

# The machine code of the three rows, to show which loads g++ 8.5 made aligned.
# g++ folds some loads into the arithmetic instruction, and those never check alignment
objdump -d -C --no-show-raw-insn src/stencil_avx.o > "$SUBMIT_DIR/results/e2_disasm_${SLURM_JOB_ID}.txt"

if [ "$PART" = "one" ]; then

# One thread on the 2048 by 512 grid. The shipped simd rung first and last, to catch drift
NT=1000
ONE="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"
export OMP_NUM_THREADS=1
run E2 simd_1_2048 $ONE --solver simd
for SHIFT in 0 1 2 3 4 5 6 7; do
  run E2 loadu_s${SHIFT}_1_2048 $ONE --solver omp_align_loadu --align-shift $SHIFT
done
for SHIFT in 0 4; do
  run E2 load_s${SHIFT}_1_2048 $ONE --solver omp_align_load --align-shift $SHIFT
  run E2 shuffle_s${SHIFT}_1_2048 $ONE --solver omp_align_shuffle --align-shift $SHIFT
done
run E2 simd_1_2048_end $ONE --solver simd

else

# 24 threads, the sheets fit the two L3 caches
export OMP_NUM_THREADS=24
NT=20000
BIG="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"
run E2 omp_24_2048 $BIG --solver omp
for SHIFT in 0 1 2 3; do
  run E2 loadu_s${SHIFT}_24_2048 $BIG --solver omp_align_loadu --align-shift $SHIFT
done
run E2 load_s0_24_2048 $BIG --solver omp_align_load --align-shift 0
run E2 shuffle_s0_24_2048 $BIG --solver omp_align_shuffle --align-shift 0

# 24 threads, the sheets spill to memory
NT=4000
HUGE="--config config/bench.cfg --ns 4096 --nv 1024 --nt $NT --maturity $(maturity_for 4096 1024 0.64 $NT)"
run E2 omp_24_4096 $HUGE --solver omp
for SHIFT in 0 1 2 3; do
  run E2 loadu_s${SHIFT}_24_4096 $HUGE --solver omp_align_loadu --align-shift $SHIFT
done
run E2 load_s0_24_4096 $HUGE --solver omp_align_load --align-shift 0
run E2 shuffle_s0_24_4096 $HUGE --solver omp_align_shuffle --align-shift 0

fi

echo "done: $OUT"
wc -l "$OUT"
