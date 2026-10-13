#!/bin/bash -l
#
#SBATCH --job-name=MM_BenchMPI
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --time=0-00:14:00
module load compiler-rt/latest
module add mkl/latest
module add mpi/openmpi-x86_64
module load cuda/11.1

export PATH=/opt/local/stow/cuda-11.1/bin:$PATH
export PATH=/usr/lib64/openmpi/bin:$PATH

make clean
make all
hostname
# Two real nodes at graded size, no nodelist pin, so the pair is whatever is
# free and the interconnect may differ from a100-a/a100-b. The final MPI
# number comes from the Phase 7 JudgementDay run.
# -map-by node puts one rank per node, -bind-to none stops OpenMPI pinning
# each rank to a single core (which would kill OpenMP threading).
mpiexec -n 2 -map-by node -bind-to none -x OMPI_MCA_mtl=^ofi,psm,psm2,psm3 -x OMPI_MCA_btl=^uct,ofi ./Assignment1_GradeBot 2048 4 0 0 1
