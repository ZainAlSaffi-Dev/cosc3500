#!/bin/bash -l
#
#SBATCH --job-name=MM_BenchEPYC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=vcpu-4
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
lscpu | head -30
# vcpu-4 has the same EPYC 7542 host CPU as the graded a100 pair, so this
# calibrates ratios on the graded silicon without touching a100-a/b.
./Assignment1_GradeBot 2048 4 1 0 0
./Assignment1_GradeBot 2048 4 1 0 0 384 96
./Assignment1_GradeBot 2048 4 1 0 0 256 192
./Assignment1_GradeBot 2048 4 1 0 0 256 96 1
./Assignment1_GradeBot 2048 4 1 0 0 384 96 1
