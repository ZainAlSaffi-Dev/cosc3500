#!/bin/bash -l
#
# Diagnostic, not a result. In the scaling sweep, 2048x512 showed rep 1 at
# 1.35e8 cell-updates/s and reps 2-5 at 8.1e7 — a 66% drop — while 4096x1024
# showed all five reps identical to three digits. That is backwards for a
# cache explanation (the bigger problem should be the unstable one), and it
# is far too large to be clock behaviour.
#
# Hypothesis: glibc's DYNAMIC mmap threshold. Grid allocates two buffers per
# solve and frees them at the end of each rep.
#   - 4096x1024 -> 32 MiB per buffer, at or above DEFAULT_MMAP_THRESHOLD_MAX,
#     so every rep gets fresh mmap'd (and transparent-huge-page eligible)
#     memory. All reps identical.
#   - 2048x512  -> 8 MiB per buffer. Rep 1 is served by mmap; when it is
#     freed, glibc RAISES its threshold past 8 MiB so later reps are served
#     from the reused main heap instead — same size, worse page backing,
#     more TLB pressure.
#
# The test: pin the threshold by hand. Forcing mmap should make every rep
# look like rep 1; forcing the heap should make rep 1 look like the others.
# If both fire, the effect is the allocator and nothing to do with the kernel.
#
#SBATCH --job-name=heston_alloc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:12:00

module load compiler-rt/latest
set -u
mkdir -p results
OUT=results/alloc_${SLURM_JOB_ID}.txt

{
echo "# host=$(hostname) date=$(date -Is) jobid=${SLURM_JOB_ID}"
echo -n "# transparent_hugepage: "
cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null || echo "unreadable"
echo

make clean >/dev/null && make heston >/dev/null
ARGS="--config config/bench.cfg --nt 400 --maturity 0.0001 --solver opt --opt-level 6 --bench 5"

echo "### default allocator behaviour (what the benchmarks measured)"
./heston $ARGS 2>/dev/null | awk -F, '{printf "  rep sec=%-8.3f cups=%s\n", $9, $10}'

echo
echo "### MALLOC_MMAP_THRESHOLD_=1048576  (force mmap for the 8 MiB buffers)"
MALLOC_MMAP_THRESHOLD_=1048576 ./heston $ARGS 2>/dev/null \
  | awk -F, '{printf "  rep sec=%-8.3f cups=%s\n", $9, $10}'

echo
echo "### MALLOC_MMAP_THRESHOLD_=134217728  (force the heap for everything)"
MALLOC_MMAP_THRESHOLD_=134217728 ./heston $ARGS 2>/dev/null \
  | awk -F, '{printf "  rep sec=%-8.3f cups=%s\n", $9, $10}'
} 2>&1 | tee "$OUT"

echo "done: $OUT"
