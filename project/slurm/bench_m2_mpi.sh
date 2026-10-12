#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Hybrid MPI and OpenMP layouts on one and two r730 nodes, experiment E10.
# One rank per node or one per socket, threads filling the rest of the cores.
#
#SBATCH --job-name=heston_m2_mpi
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2,r730-3
#SBATCH --time=0-00:14:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_mpi
source slurm/bench_common.sh
build_and_validate
write_header

# Cell updates per rep scale with the cores used, so reps take a few seconds
grid_args() {
  local ns=$1
  local nv=$2
  local cores=$3
  local nt=$((cores * 400000000 / (ns * nv)))
  echo "--config config/bench.cfg --ns $ns --nv $nv --nt $nt --maturity $(maturity_for $ns $nv 0.64 $nt) --solver mpi"
}

for GRID in "2048 512" "4096 1024"; do
  set -- $GRID
  SIZE="$1x$2"
  # One node, one rank of 24 threads
  run_mpi M "1node_1x24_$SIZE" "-n 1 --map-by ppr:1:node:PE=24 --bind-to core --report-bindings -x OMP_NUM_THREADS=24" $(grid_args $1 $2 24)
  # One node, one rank per socket
  run_mpi M "1node_2x12_$SIZE" "-n 2 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12" $(grid_args $1 $2 24)
  # Two nodes, one rank per node
  run_mpi M "2node_1x24_$SIZE" "-n 2 --map-by ppr:1:node:PE=24 --bind-to core --report-bindings -x OMP_NUM_THREADS=24" $(grid_args $1 $2 48)
  # Two nodes, one rank per socket
  run_mpi M "2node_2x12_$SIZE" "-n 4 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12" $(grid_args $1 $2 48)
done

echo "done: $OUT"
wc -l "$OUT"
