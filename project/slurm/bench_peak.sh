#!/bin/bash -l
#
# Measures the two things a roofline needs and nobody should quote from a
# datasheet: this node's SUSTAINED clock under our own load, and its
# achievable SCALAR double-precision FLOP rate.
#
# Why this job exists. The roofline built from the ladder data says level 6
# demands ~1.9 GB/s of memory bandwidth. A Xeon E5-2670 v3 has four DDR4
# channels and tens of GB/s, so if that number is right the kernel is NOT
# memory-bound — which contradicts the assumption the whole optimisation
# story was written around. Before rewriting that story, the compute roof has
# to be measured rather than assumed.
#
# The microbenchmark below is deliberately dumb: eight INDEPENDENT
# multiply-add chains on registers, no memory traffic, no divisions, compiled
# the same way the solver is (-O2, no -march, so plain SSE2 scalar — the same
# instruction set the kernel gets). Eight chains is enough to hide the ~5
# cycle latency of an addsd so the loop measures throughput, not latency.
# The result is consumed and printed so it cannot be optimised away.
#
#SBATCH --job-name=heston_peak
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:10:00

module load compiler-rt/latest
set -u
mkdir -p results
OUT=results/peak_${SLURM_JOB_ID}.txt

{
echo "# host=$(hostname) date=$(date -Is) jobid=${SLURM_JOB_ID}"
echo "# compiler=$(g++ --version | head -1)"
lscpu | grep -E 'Model name|CPU max MHz|CPU min MHz'

cat > /tmp/peak_$$.cpp <<'CPP'
#include <cstdio>
#include <chrono>
int main() {
    // Eight independent accumulators: enough in-flight work that the loop
    // measures issue THROUGHPUT rather than the latency of one chain.
    double a0=1.0,a1=1.1,a2=1.2,a3=1.3,a4=1.4,a5=1.5,a6=1.6,a7=1.7;
    const double m = 1.0000000001, c = 1e-9;
    const long iters = 200000000L;
    auto t0 = std::chrono::steady_clock::now();
    for (long i = 0; i < iters; ++i) {
        a0 = a0*m + c; a1 = a1*m + c; a2 = a2*m + c; a3 = a3*m + c;
        a4 = a4*m + c; a5 = a5*m + c; a6 = a6*m + c; a7 = a7*m + c;
    }
    double secs = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    // 8 chains x 2 flops (one multiply, one add) per iteration.
    double flops = (double)iters * 8.0 * 2.0;
    // Consume the results so the compiler cannot delete the loop.
    std::printf("checksum %.6f\n", a0+a1+a2+a3+a4+a5+a6+a7);
    std::printf("scalar_peak_gflops %.3f  (%.3f s)\n", flops/secs/1e9, secs);
    return 0;
}
CPP
g++ -std=c++17 -O2 -o /tmp/peak_$$ /tmp/peak_$$.cpp
echo
echo "### achievable scalar double FLOP rate (-O2, no -march: SSE2 scalar)"
/tmp/peak_$$
/tmp/peak_$$

echo
echo "### sustained core clock DURING an actual level-6 solve"
make clean >/dev/null && make heston >/dev/null
./heston --config config/bench.cfg --nt 2000 --maturity 0.00051 \
        --solver opt --opt-level 6 >/tmp/solve_$$.out 2>/dev/null &
SOLVE=$!
sleep 3
for i in 1 2 3 4 5 6 7 8; do
  grep -m1 'cpu MHz' /proc/cpuinfo | sed 's/^/  /'
  sleep 1
done
wait $SOLVE
echo "  solve result: $(cat /tmp/solve_$$.out)"
rm -f /tmp/peak_$$ /tmp/peak_$$.cpp /tmp/solve_$$.out
} 2>&1 | tee "$OUT"

echo "done: $OUT"
