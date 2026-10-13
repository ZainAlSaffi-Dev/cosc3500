#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Experiment E2b, the E2 sweep again with g++ told not to split unaligned 256 bit loads and stores.
# Under -mavx2 alone g++ 8.5 tunes for a generic x86 and turns every _mm256_loadu_pd into a
# 128 bit load plus a vinsertf128, and every _mm256_storeu_pd into two 128 bit stores (E2, jobs 657421 and 657422).
# Aligned loads are never split, which is why load and shuffle beat loadu in E2.
#
# The cosc3500 QOS allows 15 minutes, so the sweep runs as two jobs.
#   sbatch slurm/bench_m2_exp_nosplit.sh one    one thread
#   sbatch slurm/bench_m2_exp_nosplit.sh many   24 threads
#
#SBATCH --job-name=heston_m2_exp_nosplit
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
if [ "$PART" != "one" ] && [ "$PART" != "many" ]; then
  echo "usage: sbatch slurm/bench_m2_exp_nosplit.sh one|many"
  exit 1
fi
SWEEP_NAME=m2_exp_nosplit_${PART}

# Only the AVX2 file gets these. The Makefile adds AVXEXTRA to AVXFLAGS
export AVXEXTRA="-mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store"

source slurm/bench_common.sh
build_and_validate
write_header
echo "# avx_extra=${AVXEXTRA}" >> "$OUT"

# Keep the machine code, to check the loads really are whole 256 bit loads now
objdump -d -C --no-show-raw-insn src/stencil_avx.o > "$SUBMIT_DIR/results/e2b_disasm_${SLURM_JOB_ID}.txt"

if [ "$PART" = "one" ]; then

# One thread on the 2048 by 512 grid, the same runs as E2
NT=1000
ONE="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"
export OMP_NUM_THREADS=1
run E2b simd_1_2048 $ONE --solver simd
for SHIFT in 0 1 2 3 4 5 6 7; do
  run E2b loadu_s${SHIFT}_1_2048 $ONE --solver omp_align_loadu --align-shift $SHIFT
done
for SHIFT in 0 4; do
  run E2b load_s${SHIFT}_1_2048 $ONE --solver omp_align_load --align-shift $SHIFT
  run E2b shuffle_s${SHIFT}_1_2048 $ONE --solver omp_align_shuffle --align-shift $SHIFT
done
run E2b simd_1_2048_end $ONE --solver simd

else

# 24 threads, the sheets fit the two L3 caches
export OMP_NUM_THREADS=24
NT=20000
BIG="--config config/bench.cfg --nt $NT --maturity $(maturity_for 2048 512 0.64 $NT)"
run E2b omp_24_2048 $BIG --solver omp
for SHIFT in 0 1 2 3; do
  run E2b loadu_s${SHIFT}_24_2048 $BIG --solver omp_align_loadu --align-shift $SHIFT
done
run E2b load_s0_24_2048 $BIG --solver omp_align_load --align-shift 0
run E2b shuffle_s0_24_2048 $BIG --solver omp_align_shuffle --align-shift 0

# 24 threads, the sheets spill to memory. Fewer runs, E2 took 13:57 of the 14 minutes
NT=4000
HUGE="--config config/bench.cfg --ns 4096 --nv 1024 --nt $NT --maturity $(maturity_for 4096 1024 0.64 $NT)"
run E2b omp_24_4096 $HUGE --solver omp
run E2b loadu_s0_24_4096 $HUGE --solver omp_align_loadu --align-shift 0
run E2b shuffle_s0_24_4096 $HUGE --solver omp_align_shuffle --align-shift 0

fi

echo "done: $OUT"
wc -l "$OUT"
