#!/bin/bash -l
#
# 5-minute sanity job: does it build and run on a compute node at all?
# Submit early and often (PLAN §7 risk table) — do not wait for P7.
#
#SBATCH --job-name=heston_debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --time=0-00:05:00

module load compiler-rt/latest

make clean
make all
hostname
g++ --version | head -1

./heston --config config/smoke.cfg
make test
