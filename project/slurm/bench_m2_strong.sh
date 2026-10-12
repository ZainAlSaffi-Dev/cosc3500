#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Strong scaling, the 2048 by 512 grid on 1 to 24 threads.
# Pass the solver as the first argument, omp for the AVX2 row or omp_scalar for the scalar row.
#
#SBATCH --job-name=heston_m2_strong
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SOLVER=${1:-omp}
SWEEP_NAME=m2_strong_${SOLVER}
source slurm/bench_common.sh
build_and_validate
write_header

# Timesteps per thread, so every rep takes about the same few seconds
if [ "$SOLVER" = "omp_scalar" ]; then
  STEPS_PER_THREAD=400
else
  STEPS_PER_THREAD=1000
fi

for THREADS in 1 2 4 6 8 10 12 16 20 24; do
  NT=$((STEPS_PER_THREAD * THREADS))
  OMP_NUM_THREADS=$THREADS run S "${SOLVER}_${THREADS}" --config config/bench.cfg --solver "$SOLVER" --nt "$NT" --maturity $(maturity_for 2048 512 0.64 "$NT")
done

echo "done: $OUT"
wc -l "$OUT"
