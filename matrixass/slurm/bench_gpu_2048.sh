#!/bin/bash -l
#
#SBATCH --job-name=MM_BenchGPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres gpu:1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --time=0-00:10:00
module load compiler-rt/latest
module add mkl/latest
module add mpi/openmpi-x86_64
module load cuda/11.1

export PATH=/opt/local/stow/cuda-11.1/bin:$PATH
export PATH=/usr/lib64/openmpi/bin:$PATH

make clean
make all
hostname
nvidia-smi
# Full GPU (gres gpu:1), unlike the DebugGPU script's shard:1 shared slice.
# Only timings from this script count as real GPU numbers.
./Assignment1_GradeBot 2048 4 0 1 0
# Kernel A/B at the graded size, arg picks the kernel, 1 naive, 2 shared
# tiled. The deltas bracket how the time splits between kernel and PCIe
# transfers.
./Assignment1_GradeBot 2048 4 0 1 0 1
./Assignment1_GradeBot 2048 4 0 1 0 2
# Default kernel at other sizes. Transfers scale as N^2 and compute as N^3,
# so the scaling separates the two.
./Assignment1_GradeBot 1024 4 0 1 0
./Assignment1_GradeBot 4096 4 0 1 0
