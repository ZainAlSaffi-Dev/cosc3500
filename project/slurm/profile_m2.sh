#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# gprof flat profiles of the baseline, serial and simd rungs on one thread.
# Shows which functions the time goes into before anything is parallelised.
#
#SBATCH --job-name=heston_m2_profile
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_profile
source slurm/bench_common.sh
build_and_validate

export OMP_NUM_THREADS=1

# -pg makes the program record where it is every 10 ms into gmon.out
make clean > /dev/null
make heston CXXFLAGS="-std=c++17 -O2 -Wall -Wextra -Iinclude -pg" > /dev/null || { echo "FATAL: -pg build failed"; exit 1; }

for SOLVER in baseline serial simd; do
  if [ "$SOLVER" = "baseline" ]; then
    NT=400
  else
    NT=1000
  fi
  PROFILE=$SUBMIT_DIR/results/profile_${SOLVER}_${SLURM_JOB_ID}.txt
  rm -f gmon.out
  {
    echo "# host=$(hostname)"
    echo "# jobid=${SLURM_JOB_ID}"
    echo "# compiler=$(mpicxx --version | head -1)"
    echo "# flags=-std=c++17 -O2 -pg, plus -mavx2 -mfma on stencil_avx.cpp"
    echo "# run=./heston --config config/bench.cfg --solver $SOLVER --nt $NT, 2048 x 512, one thread, denormals flushed"
    env "$(preload_setting)" ./heston --config config/bench.cfg --solver "$SOLVER" --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT) 2>/dev/null | sed 's/^/# csv=/'
    # -b drops the long explanation, -p prints only the flat profile
    gprof -b -p ./heston gmon.out
  } > "$PROFILE"
  echo "wrote $PROFILE"
  head -20 "$PROFILE"
done
