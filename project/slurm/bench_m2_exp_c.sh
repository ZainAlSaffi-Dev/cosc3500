#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Experiment E11, the pinned page regime against glibc's default, threaded and serial.
# Split from bench_m2_exp_a.sh, which ran out of time before reaching it.
#
#SBATCH --job-name=heston_m2_exp_c
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_exp_c
source slurm/bench_common.sh
build_and_validate
write_header

export OMP_NUM_THREADS=24
NT=20000
BIG="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"

# E11, the pinned page regime against glibc's default
run E11 pinned_2048 $BIG --solver omp
RUN_PREFIX="env -u MALLOC_MMAP_THRESHOLD_" run E11 default_2048 $BIG --solver omp
run E11 pinned_serial $BIG --solver serial --nt 1000 --maturity $(maturity_for 2048 512 0.64 1000)
RUN_PREFIX="env -u MALLOC_MMAP_THRESHOLD_" run E11 default_serial $BIG --solver serial --nt 1000 --maturity $(maturity_for 2048 512 0.64 1000)

echo "done: $OUT"
wc -l "$OUT"
