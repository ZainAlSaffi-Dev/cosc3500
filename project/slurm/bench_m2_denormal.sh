#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Experiment E12, what denormal numbers cost each rung.
# Every ladder rung at its ladder timestep count, with and without flush to zero.
#
#SBATCH --job-name=heston_m2_denormal
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_denormal
source slurm/bench_common.sh
build_and_validate
write_header
export OMP_NUM_THREADS=24

# E12, every ladder rung with and without flushing denormals, at its ladder timestep count
ladder_args() {
  echo "--config config/bench.cfg --nt $1 --maturity $(maturity_for 2048 512 0.64 $1)"
}
for FLUSH in 1 0; do
  FLUSH_DENORMALS=$FLUSH OMP_NUM_THREADS=1 run E12 "baseline_ftz$FLUSH" $(ladder_args 400) --solver baseline
  FLUSH_DENORMALS=$FLUSH OMP_NUM_THREADS=1 run E12 "serial_ftz$FLUSH" $(ladder_args 1000) --solver serial
  FLUSH_DENORMALS=$FLUSH OMP_NUM_THREADS=1 run E12 "serial_short_ftz$FLUSH" $(ladder_args 381) --solver serial
  FLUSH_DENORMALS=$FLUSH OMP_NUM_THREADS=1 run E12 "simd_ftz$FLUSH" $(ladder_args 2000) --solver simd
  FLUSH_DENORMALS=$FLUSH run E12 "omp24_ftz$FLUSH" $(ladder_args 20000) --solver omp
done

echo "done: $OUT"
wc -l "$OUT"
