#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# The Milestone 2 ladder, every rung at the 2048 by 512 benchmark grid.
# Each rung gets enough timesteps that one rep takes several seconds.
#
#SBATCH --job-name=heston_m2_ladder
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_ladder
source slurm/bench_common.sh
build_and_validate
write_header

ARGS="--config config/bench.cfg"
# Runs one rung with nt timesteps and the maturity that keeps dt stable
ladder_run() {
  local label=$1
  local nt=$2
  shift 2
  run L "$label" $ARGS --nt "$nt" --maturity $(maturity_for 2048 512 0.64 "$nt") "$@"
}

ladder_run baseline 400 --solver baseline
ladder_run serial 1000 --solver serial
ladder_run simd 2000 --solver simd
OMP_NUM_THREADS=24 ladder_run omp24 20000 --solver omp
# Two processes on one node, one per socket, twelve threads each
run_mpi L mpi2x12 "-n 2 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12" $ARGS --solver mpi --nt 20000 --maturity $(maturity_for 2048 512 0.64 20000)

echo "done: $OUT"
wc -l "$OUT"
