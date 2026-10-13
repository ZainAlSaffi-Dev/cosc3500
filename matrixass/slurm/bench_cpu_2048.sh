#!/bin/bash -l
#
#SBATCH --job-name=MM_BenchCPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
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
lscpu | head -30
# Graded size, CPU benchmark only, runs on whatever node is free.
# Real CPU numbers come from here, not from the 5-minute debug script.
./Assignment1_GradeBot 2048 4 1 0 0
