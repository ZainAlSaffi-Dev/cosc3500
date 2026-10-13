#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Experiment E1, what the compiler can do on its own.
# The whole program is rebuilt under several sets of flags, and each build is timed on one thread.
# The question is whether the compiler could have written the hand-made AVX2 rung for us.
# The earlier -O3 bar (M2.11) used only 128-bit SSE2, because no -mavx2 was given, so it did not answer that.
#
# Each build records three things.
#   Speed of the baseline, serial and simd rungs on the 2048 by 512 grid.
#   How far the answer moves from the -O2 answer, printed to 17 digits on config/smoke.cfg.
#   Whether the compiler vectorised the scalar row, from its own report and from the machine code.
#
# The cosc3500 QOS allows 15 minutes, so the sweep runs as two jobs.
#   sbatch slurm/bench_m2_exp_compiler.sh vector     -O2, -O3, -O3 with AVX2 and FMA, -O3 -march=native
#   sbatch slurm/bench_m2_exp_compiler.sh fastmath   -O2, then -ffast-math on top of -O2, -O3, and -O3 with AVX2 and FMA
#
#SBATCH --job-name=heston_m2_exp_compiler
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
PART=${1:-}
if [ "$PART" != "vector" ] && [ "$PART" != "fastmath" ]; then
  echo "usage: sbatch slurm/bench_m2_exp_compiler.sh vector|fastmath"
  exit 1
fi
SWEEP_NAME=m2_exp_compiler_${PART}
source slurm/bench_common.sh
build_and_validate
write_header
echo "# note=each row was built with the flags named in its label, see the BUILD lines in the job log" >> "$OUT"

BASE="-std=c++17 -Wall -Wextra -Iinclude"
# Each entry is a short label, a colon, then the flags that replace -O2
if [ "$PART" = "vector" ]; then
  BUILDS=("O2:-O2" "O3:-O3" "O3avx2:-O3 -mavx2 -mfma" "O3native:-O3 -march=native")
else
  BUILDS=("O2:-O2" "O2fast:-O2 -ffast-math" "O3fast:-O3 -ffast-math" "O3avx2fast:-O3 -mavx2 -mfma -ffast-math")
fi

export OMP_NUM_THREADS=1
BASELINE_NT=200
SERIAL_NT=400
SIMD_NT=1000
BASELINE_ARGS="--config config/bench.cfg --nt $BASELINE_NT --maturity $(maturity_for 2048 512 0.64 $BASELINE_NT)"
SERIAL_ARGS="--config config/bench.cfg --nt $SERIAL_NT --maturity $(maturity_for 2048 512 0.64 $SERIAL_NT)"
SIMD_ARGS="--config config/bench.cfg --nt $SIMD_NT --maturity $(maturity_for 2048 512 0.64 $SIMD_NT)"

# The answers of each build on the real contract, one line per build and rung
DRIFT=$SUBMIT_DIR/results/e1_answers_${PART}_${SLURM_JOB_ID}.csv
echo "build,solver,price,delta,gamma,vega" > "$DRIFT"
# What the compiler vectorised in each build
VECTOR_LOG=$SUBMIT_DIR/results/e1_vectorised_${PART}_${SLURM_JOB_ID}.txt
: > "$VECTOR_LOG"

# Counts the instructions of one function in the disassembly that use 256-bit (ymm) registers,
# that work on packed doubles in 128-bit (xmm) registers, and that work on one double at a time
count_instructions() {
  local name=$1
  objdump -d -C --no-show-raw-insn heston | awk -v name="$name" '
    /^[0-9a-f]+ <.*>:$/ { inside = (index($0, "<" name) > 0) }
    inside && /\t/ {
      if ($0 ~ /ymm/) ymm++
      else if ($0 ~ /pd[ \t]/ && $0 ~ /xmm/) xmm_packed++
      else if ($0 ~ /sd[ \t]/) scalar++
    }
    END { printf "  %-28s ymm %d, packed xmm %d, scalar double %d\n", name, ymm, xmm_packed, scalar }'
}

for ENTRY in "${BUILDS[@]}"; do
  BUILD=${ENTRY%%:*}
  FLAGS=${ENTRY#*:}
  echo "BUILD $BUILD = $FLAGS"
  make clean > /dev/null
  make heston CXXFLAGS="$BASE $FLAGS" > /dev/null || { echo "FATAL: build $BUILD failed"; exit 1; }

  # The compiler's own report of the loops it vectorised in the two hot files
  {
    echo "== $BUILD ($FLAGS)"
    for FILE in src/stencil.cpp src/solver_baseline.cpp; do
      mpicxx $BASE $FLAGS -DOMPI_SKIP_MPICXX -fopenmp -fopt-info-vec-optimized -c $FILE -o /dev/null 2>&1 | sed "s|^|  |"
    done
    count_instructions "update_cells_scalar"
    count_instructions "BaselineSolver::solve"
    count_instructions "update_row_avx("
  } >> "$VECTOR_LOG"

  # The answer to 17 digits, the same contract for every build
  for SOLVER in baseline serial simd; do
    ANSWER=$(env "$(preload_setting)" ./heston --config config/smoke.cfg --solver $SOLVER --digits 17 2>/dev/null | cut -d, -f1-4)
    echo "$BUILD,$SOLVER,$ANSWER" >> "$DRIFT"
  done

  run E1 "baseline_$BUILD" $BASELINE_ARGS --solver baseline
  run E1 "serial_$BUILD" $SERIAL_ARGS --solver serial
  run E1 "simd_$BUILD" $SIMD_ARGS --solver simd
done

# The largest relative change of the price and the three Greeks against the -O2 build of the same rung
echo "relative drift against -O2:"
awk -F, 'NR > 1 {
  key = $2
  if ($1 == "O2") { for (k = 3; k <= 6; k++) reference[key, k] = $k }
  worst = 0
  for (k = 3; k <= 6; k++) {
    a = $k; b = reference[key, k]
    scale = (a < 0 ? -a : a); if ((b < 0 ? -b : b) > scale) scale = (b < 0 ? -b : b)
    diff = a - b; if (diff < 0) diff = -diff
    if (scale > 0 && diff / scale > worst) worst = diff / scale
  }
  printf "  %-11s %-9s %.3e\n", $1, $2, worst
}' "$DRIFT"

echo "done: $OUT"
wc -l "$OUT"
