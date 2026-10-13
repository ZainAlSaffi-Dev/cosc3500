#!/bin/bash -l
#
#SBATCH --job-name=MM_Sweep
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
# CPU tuning sweep. args[0]=Kc panel depth, args[1]=Mc row band,
# args[2]=1 selects the 3M scheme. Each invocation is a full 8-rep
# benchmark. The graded run passes NO args, defaults must match the winner.
for MC in 48 96 192; do
    for KC in 192 256 384; do
        echo "=== Kc=$KC Mc=$MC 4M ==="
        ./Assignment1_GradeBot 2048 4 1 0 0 $KC $MC
    done
done
echo "=== Kc=256 Mc=96 3M ==="
./Assignment1_GradeBot 2048 4 1 0 0 256 96 1
echo "=== Kc=384 Mc=96 3M ==="
./Assignment1_GradeBot 2048 4 1 0 0 384 96 1
# GPU variant, uncomment and move to a gres gpu:1 job before using.
# for FLAG in 0 1; do ./Assignment1_GradeBot 2048 4 0 1 0 $FLAG; done
