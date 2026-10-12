#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Throughput against grid size for serial, simd and omp on 24 threads.
# The sizes run from one that fits in L2 to one larger than both L3 caches.
#
#SBATCH --job-name=heston_m2_size
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_size
source slurm/bench_common.sh
build_and_validate
write_header

# Cell updates per rep for each solver, so each rep takes a few seconds
size_run() {
  local solver=$1
  local label=$2
  local updates=$3
  local ns=$4
  local nv=$5
  local nt=$((updates / (ns * nv)))
  run Z "${label}_${ns}x${nv}" --config config/bench.cfg --solver "$solver" --ns "$ns" --nv "$nv" --nt "$nt" --maturity $(maturity_for "$ns" "$nv" 0.64 "$nt")
}

for GRID in "512 128" "1024 256" "2048 512" "2896 724" "4096 1024" "8192 1024"; do
  set -- $GRID
  size_run serial serial 300000000 $1 $2
  size_run simd simd 700000000 $1 $2
  OMP_NUM_THREADS=24 size_run omp omp24 8000000000 $1 $2
done

echo "done: $OUT"
wc -l "$OUT"
