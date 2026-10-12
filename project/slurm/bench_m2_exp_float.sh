#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Experiment E3, float sheets against double sheets.
# Speed on one thread and on 24, on a grid that fits the two L3 caches and one that does not.
# Accuracy on the reference grid, whose Fourier price is 196.1692.
#
#SBATCH --job-name=heston_m2_exp_float
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_exp_float
source slurm/bench_common.sh
build_and_validate
write_header

# One thread, simd against the float row on the 2048 by 512 grid
NT=1000
ONE="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"
OMP_NUM_THREADS=1 run E3 simd_1_2048 $ONE --solver simd
OMP_NUM_THREADS=1 run E3 float_1_2048 $ONE --solver omp_float

# 24 threads, the double sheets fill 16 MB of the 60 MB of L3, the float sheets 8 MB
NT=20000
BIG="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"
OMP_NUM_THREADS=24 run E3 omp_24_2048 $BIG --solver omp
OMP_NUM_THREADS=24 run E3 float_24_2048 $BIG --solver omp_float

# 24 threads, the double sheets take 64 MB and spill to memory, the float sheets 32 MB
NT=4000
HUGE="--config config/bench.cfg --ns 4096 --nv 1024 --nt $NT --maturity $(maturity_for 4096 1024 0.64 $NT)"
OMP_NUM_THREADS=24 run E3 omp_24_4096 $HUGE --solver omp
OMP_NUM_THREADS=24 run E3 float_24_4096 $HUGE --solver omp_float

# Accuracy, one run each on the reference grid with its own 980000 steps
OMP_NUM_THREADS=24 REPS=1 run E3 omp_reference --config config/reference.cfg --solver omp
OMP_NUM_THREADS=24 REPS=1 run E3 float_reference --config config/reference.cfg --solver omp_float

echo "done: $OUT"
wc -l "$OUT"
