#!/bin/bash -l
#
# The real Milestone 1 benchmark (PLAN §6): baseline vs opt across grid
# sizes, >=5 reps each, on a compute node. Output: results/bench_<jobid>.csv
# (one PLAN §3 CSV line per rep, prefixed with grid label) — pull with
# `make fetch`, plot with scripts/bench_plot.py.
#
#SBATCH --job-name=heston_bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --time=0-01:00:00

module load compiler-rt/latest

make clean
make all
make test   # iron rule: no benchmark on an unvalidated build

OUT=results/bench_${SLURM_JOB_ID}.csv
mkdir -p results
{
  echo "# host=$(hostname) date=$(date -Is) flags='-std=c++17 -O2' job=${SLURM_JOB_ID}"
  echo "# grid,solver_rep_lines_follow"
} > "$OUT"

REPS=5
# ns nv nt triples — nt scaled with grid so the scheme stays stable (PLAN §6)
for SPEC in "512 128 2000" "1024 256 2000" "2048 512 2000" "4096 1024 4000"; do
  set -- $SPEC
  NS=$1; NV=$2; NT=$3
  for SOLVER in baseline opt; do
    echo "grid=${NS}x${NV}x${NT} solver=${SOLVER}"
    ./heston --config config/reference.cfg \
             --solver "$SOLVER" --ns "$NS" --nv "$NV" --nt "$NT" \
             --bench "$REPS" >> "$OUT"
  done
done

echo "done: $OUT"
