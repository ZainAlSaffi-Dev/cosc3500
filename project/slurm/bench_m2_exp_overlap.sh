#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Experiment E9, overlapping the edge row exchange with the update (the L06 overlap pattern).
# The mpi rung posts its non-blocking sends and receives after the update and waits for them at once,
# so the whole exchange sits on the critical path. On two nodes that was 38% of a step (M2.9).
# The mpi_overlap rung updates the two rows the neighbours need first, starts sending them,
# updates the rest of its rows while the messages travel, and waits only at the end.
# Every layout runs both rungs back to back. The exchange column holds the time the exchange
# still costs, the time to post plus the time spent waiting.
# The one rank layout has no neighbours, so it measures only the cost of the reordered loop.
#
#SBATCH --job-name=heston_m2_exp_overlap
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2,r730-3
#SBATCH --time=0-00:14:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_exp_overlap
source slurm/bench_common.sh
build_and_validate
write_header

# Cell updates per rep scale with the cores used, so reps take a few seconds.
# The same sizes as the hybrid sweep, so the mpi rows can be checked against M2.9
grid_args() {
  local ns=$1
  local nv=$2
  local cores=$3
  local nt=$((cores * 400000000 / (ns * nv)))
  echo "--config config/bench.cfg --ns $ns --nv $nv --nt $nt --maturity $(maturity_for $ns $nv 0.64 $nt)"
}

ONE_NODE_1x24="-n 1 --map-by ppr:1:node:PE=24 --bind-to core --report-bindings -x OMP_NUM_THREADS=24"
ONE_NODE_2x12="-n 2 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12"
TWO_NODE_1x24="-n 2 --map-by ppr:1:node:PE=24 --bind-to core --report-bindings -x OMP_NUM_THREADS=24"
TWO_NODE_2x12="-n 4 --map-by ppr:1:socket:PE=12 --bind-to core --report-bindings -x OMP_NUM_THREADS=12"

for GRID in "2048 512" "4096 1024"; do
  set -- $GRID
  SIZE="$1x$2"
  for SOLVER in mpi mpi_overlap; do
    run_mpi E9 "1node_1x24_${SIZE}_$SOLVER" "$ONE_NODE_1x24" $(grid_args $1 $2 24) --solver $SOLVER
  done
  for SOLVER in mpi mpi_overlap; do
    run_mpi E9 "1node_2x12_${SIZE}_$SOLVER" "$ONE_NODE_2x12" $(grid_args $1 $2 24) --solver $SOLVER
  done
  for SOLVER in mpi mpi_overlap; do
    run_mpi E9 "2node_1x24_${SIZE}_$SOLVER" "$TWO_NODE_1x24" $(grid_args $1 $2 48) --solver $SOLVER
  done
  for SOLVER in mpi mpi_overlap; do
    run_mpi E9 "2node_2x12_${SIZE}_$SOLVER" "$TWO_NODE_2x12" $(grid_args $1 $2 48) --solver $SOLVER
  done
done

echo "done: $OUT"
wc -l "$OUT"
