#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Checks whether denormal numbers slow down the short benchmark runs on big grids.
# Runs each case twice, once as normal and once with flush to zero switched on.
#
#SBATCH --job-name=heston_diag_denormal
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-3
#SBATCH --time=0-00:12:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

set -u
SWEEP_NAME=diag_denormal
REPS=1
source slurm/bench_common.sh
make clean >/dev/null
make heston >/dev/null
write_header

# A tiny library that sets FTZ and DAZ in the MXCSR register when heston starts.
# FTZ writes denormal results as zero, DAZ reads denormal inputs as zero
cat > /tmp/ftz_$$.c <<'C'
#include <xmmintrin.h>
__attribute__((constructor)) static void flush_denormals(void) { _mm_setcsr(_mm_getcsr() | 0x8040); }
C
gcc -O2 -shared -fPIC /tmp/ftz_$$.c -o /tmp/ftz_$$.so

for CASE in "2048 512 381" "2048 512 3810" "2896 724 190" "2896 724 1900" "4096 1024 95" "4096 1024 950"; do
  set -- $CASE
  ARGS="--config config/bench.cfg --solver serial --ns $1 --nv $2 --nt $3 --maturity $(maturity_for $1 $2 0.64 $3)"
  run D "normal_$1x$2_nt$3" $ARGS
  RUN_PREFIX="env LD_PRELOAD=/tmp/ftz_$$.so" run D "ftz_$1x$2_nt$3" $ARGS
done
rm -f /tmp/ftz_$$.c /tmp/ftz_$$.so
grep -v '^#' "$OUT"
