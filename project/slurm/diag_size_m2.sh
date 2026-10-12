#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Looks at the slow simd point at 2896 by 724 in the size sweep.
# Changes the row width and the row count separately to see which one matters.
#
#SBATCH --job-name=heston_diag_size
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-3
#SBATCH --time=0-00:12:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

set -u
SWEEP_NAME=diag_size
REPS=3
source slurm/bench_common.sh
make clean >/dev/null
make heston >/dev/null
write_header

# One thread, the same work per rep as the size sweep
simd_run() {
  local ns=$1
  local nv=$2
  local nt=$((700000000 / (ns * nv)))
  OMP_NUM_THREADS=1 run X "$3_${ns}x${nv}" --config config/bench.cfg --solver $3 --ns $ns --nv $nv --nt $nt --maturity $(maturity_for $ns $nv 0.64 $nt)
}

# Grids whose sheet plus malloc's 16 byte header fills exactly 16 MiB of pages,
# so the two sheets sit exactly 2^24 bytes apart, against near neighbours that do not
simd_run 2896 724 simd
simd_run 1448 1448 simd
simd_run 2044 1026 simd
simd_run 2044 1027 simd
simd_run 2880 724 simd
# The same grid on 2 and 4 threads
OMP_NUM_THREADS=2 run X "omp2_2896x724" --config config/bench.cfg --solver omp --ns 2896 --nv 724 --nt 666 --maturity $(maturity_for 2896 724 0.64 666)
OMP_NUM_THREADS=2 run X "omp2_2880x724" --config config/bench.cfg --solver omp --ns 2880 --nv 724 --nt 666 --maturity $(maturity_for 2880 724 0.64 666)
grep -v '^#' "$OUT" | cut -d, -f2,12
