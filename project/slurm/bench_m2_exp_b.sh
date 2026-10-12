#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# OpenMP experiments E6 and E7, thread binding and NUMA page placement.
# Both grids, one inside the two L3 caches and one larger than both.
#
#SBATCH --job-name=heston_m2_exp_b
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_exp_b
source slurm/bench_common.sh
build_and_validate
write_header

# Cell updates per rep at each thread count, so reps take a few seconds
grid_args() {
  local ns=$1
  local nv=$2
  local threads=$3
  local nt=$((threads * 600000000 / (ns * nv)))
  echo "--config config/bench.cfg --ns $ns --nv $nv --nt $nt --maturity $(maturity_for $ns $nv 0.64 $nt)"
}

# E6, where the threads run
for GRID in "2048 512" "4096 1024"; do
  set -- $GRID
  for THREADS in 6 12 24; do
    for BIND in false close spread; do
      OMP_NUM_THREADS=$THREADS OMP_PROC_BIND=$BIND run E6 "${BIND}_${THREADS}_$1x$2" $(grid_args $1 $2 $THREADS) --solver omp
    done
  done
done

# E7, where the pages live, all on the main thread's socket or spread over both
for GRID in "2048 512" "4096 1024"; do
  set -- $GRID
  OMP_NUM_THREADS=24 run E7 "default_24_$1x$2" $(grid_args $1 $2 24) --solver omp
  OMP_NUM_THREADS=24 RUN_PREFIX="numactl --interleave=all" run E7 "interleave_24_$1x$2" $(grid_args $1 $2 24) --solver omp
done

echo "done: $OUT"
wc -l "$OUT"
