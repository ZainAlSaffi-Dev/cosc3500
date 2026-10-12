#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# The baseline and serial rungs built at -O0, -O1, -O2 and -O3.
# Shows how much of the speed comes from the compiler before any change of ours.
#
#SBATCH --job-name=heston_m2_optlevels
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:14:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12

set -u
SWEEP_NAME=m2_optlevels
source slurm/bench_common.sh
build_and_validate
write_header
# The flags line above is the default build. The label of each row says which -O level it ran
echo "# note=each row was built with CXXFLAGS=-std=c++17 -O<level> -Wall -Wextra -Iinclude, the level is in the label" >> "$OUT"

export OMP_NUM_THREADS=1

# The -O2 answer on a small grid, which every other level must give exactly
SMALL="--config config/smoke.cfg"
REFERENCE_ANSWER=$(LD_PRELOAD=$FTZ_LIB ./heston $SMALL --solver baseline 2>/dev/null | cut -d, -f1-4)
echo "reference answer at -O2: $REFERENCE_ANSWER"

BASELINE_NT=200
SERIAL_NT=400
BASELINE_ARGS="--config config/bench.cfg --nt $BASELINE_NT --maturity $(maturity_for 2048 512 0.64 $BASELINE_NT)"
SERIAL_ARGS="--config config/bench.cfg --nt $SERIAL_NT --maturity $(maturity_for 2048 512 0.64 $SERIAL_NT)"

for LEVEL in O0 O1 O2 O3; do
  # Rebuild only the program, with the one flag changed
  make clean > /dev/null
  make heston CXXFLAGS="-std=c++17 -$LEVEL -Wall -Wextra -Iinclude" > /dev/null || { echo "FATAL: -$LEVEL build failed"; exit 1; }

  # Refuse to time a build whose answer moved
  ANSWER=$(LD_PRELOAD=$FTZ_LIB ./heston $SMALL --solver baseline 2>/dev/null | cut -d, -f1-4)
  echo "answer at -$LEVEL: $ANSWER"
  if [ "$ANSWER" != "$REFERENCE_ANSWER" ]; then
    echo "FATAL: -$LEVEL changed the answer, refusing to benchmark it"
    exit 1
  fi

  run OPT "baseline_$LEVEL" $BASELINE_ARGS --solver baseline
  run OPT "serial_$LEVEL" $SERIAL_ARGS --solver serial
done

echo "done: $OUT"
wc -l "$OUT"
