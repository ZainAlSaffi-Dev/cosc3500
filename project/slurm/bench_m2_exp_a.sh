#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# OpenMP experiments E4, E5 and E8 at 24 threads.
# Persistent region, schedule, and parallel boundaries.
#
#SBATCH --job-name=heston_m2_exp_a
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_exp_a
source slurm/bench_common.sh
build_and_validate
write_header

export OMP_NUM_THREADS=24
NT=20000
BIG="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"
# The small grid is where a per step cost shows most
SMALL_NT=200000
SMALL="--config config/bench.cfg --ns 512 --nv 128 --nt $SMALL_NT --maturity $(maturity_for 512 128 0.64 $SMALL_NT)"

# E4, one fork per step against one parallel region for the whole loop
run E4 omp_2048 $BIG --solver omp
run E4 persistent_2048 $BIG --solver omp_persistent
run E4 omp_512 $SMALL --solver omp
run E4 persistent_512 $SMALL --solver omp_persistent

# E5, the schedule, read from OMP_SCHEDULE by the omp_runtime solver
for SCHEDULE in static dynamic,1 dynamic,8 guided; do
  OMP_SCHEDULE=$SCHEDULE run E5 "${SCHEDULE/,/_}" $BIG --solver omp_runtime
done

# E8, boundaries serial against boundaries inside the parallel loop
run E8 omp_2048 $BIG --solver omp
run E8 edges_2048 $BIG --solver omp_edges
run E8 omp_512 $SMALL --solver omp
run E8 edges_512 $SMALL --solver omp_edges

echo "done: $OUT"
wc -l "$OUT"
