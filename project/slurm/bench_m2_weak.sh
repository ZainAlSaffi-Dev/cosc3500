#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Weak scaling, 21 interior rows per thread on 1 to 24 threads.
# vmax grows with the rows so the variance spacing, and the work per cell, stays the same.
#
#SBATCH --job-name=heston_m2_weak
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_weak
source slurm/bench_common.sh
build_and_validate
write_header

ROWS_PER_THREAD=21
NT=20000
for THREADS in 1 2 4 6 8 10 12 16 20 24; do
  # 21 interior rows per thread plus the v = 0 row and the top row
  NV=$((ROWS_PER_THREAD * THREADS + 2))
  # Same spacing as the 512 row grid, 0.64 / 511 per row
  VMAX=$(awk -v nv="$NV" 'BEGIN { printf "%.8f", 0.64 * (nv - 1) / 511 }')
  OMP_NUM_THREADS=$THREADS run W "omp_${THREADS}" --config config/bench.cfg --solver omp --nv "$NV" --vmax "$VMAX" --nt "$NT" --maturity $(maturity_for 2048 "$NV" "$VMAX" "$NT")
done

echo "done: $OUT"
wc -l "$OUT"
