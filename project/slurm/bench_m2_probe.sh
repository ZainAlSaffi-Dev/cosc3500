#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Short look at the benchmark node before the real sweeps.
# Prints the core layout, checks the MPI binding, and times one rep of each rung.
#
#SBATCH --job-name=heston_probe
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:10:00

set -u
SWEEP_NAME=probe
REPS=1
source slurm/bench_common.sh

# The core layout, sockets and NUMA nodes
lscpu
numactl --hardware

build_and_validate
write_header

ARGS="--config config/bench.cfg"
run P baseline $ARGS --solver baseline --nt 100 --maturity $(maturity_for 2048 512 0.64 100)
run P serial $ARGS --solver serial --nt 200 --maturity $(maturity_for 2048 512 0.64 200)
run P simd $ARGS --solver simd --nt 400 --maturity $(maturity_for 2048 512 0.64 400)
OMP_NUM_THREADS=12 run P omp12 $ARGS --solver omp --nt 2000 --maturity $(maturity_for 2048 512 0.64 2000)
OMP_NUM_THREADS=24 run P omp24 $ARGS --solver omp --nt 2000 --maturity $(maturity_for 2048 512 0.64 2000)
run_mpi P mpi2x12 "-n 2 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12" $ARGS --solver mpi --nt 2000 --maturity $(maturity_for 2048 512 0.64 2000)

cat "$OUT"
echo "--- stderr, bindings ---"
grep -i -E 'bound|binding|MCW' "$ERR" | head -20
