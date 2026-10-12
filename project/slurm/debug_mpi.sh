#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Checks the mpi rung across two separate machines.
# Two nodes, two processes per node, two threads per process.
#
#SBATCH --job-name=heston_debug_mpi
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --time=0-00:05:00

# Compiler runtime and Open MPI, which provides mpicxx and mpiexec
module load compiler-rt/latest
module add mpi/openmpi-x86_64
export PATH=/usr/lib64/openmpi/bin:$PATH
# Skip the network transports rangpur's Open MPI cannot start
export OMPI_MCA_mtl=^ofi,psm,psm2,psm3
export OMPI_MCA_btl=^uct,ofi

# Build in a private copy of the tree, so a running benchmark is never disturbed
WORK_DIR=$HOME/cosc3500/work/job_${SLURM_JOB_ID}
mkdir -p "$WORK_DIR"
rsync -a --exclude results --exclude '*.o' --exclude heston --exclude 'slurm-*.out' ./ "$WORK_DIR/"
cd "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT

# The home directory is shared, so one build serves both nodes
make clean
make all
mpicxx --version | head -1
echo "nodes: $SLURM_JOB_NODELIST"

export OMP_NUM_THREADS=2

# One process per node, so the halo rows cross the network
echo "== 2 processes, one per node"
mpiexec -n 2 --map-by node --bind-to none --report-bindings ./tests/test_ladder_matches 2>/dev/null

# Three processes over two nodes, an uneven split of nodes and of rows
echo "== 3 processes over two nodes"
mpiexec -n 3 --map-by node --bind-to none ./tests/test_ladder_matches 2>/dev/null

# Four processes, two per node
echo "== 4 processes, two per node"
mpiexec -n 4 --map-by node --bind-to none ./tests/test_ladder_matches 2>/dev/null

# Which host each process ran on, to prove two machines were used
mpiexec -n 4 --map-by node --bind-to none hostname
