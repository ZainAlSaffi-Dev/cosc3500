#!/bin/bash -l
#
# Measures the two things a roofline needs and nobody should quote from a
# datasheet: this node's SUSTAINED clock under our own load, and its
# achievable SCALAR double-precision FLOP rate.
#
# Why this job exists. The roofline built from the ladder data says level 6
# demands ~1.9 GB/s of memory bandwidth. A Xeon E5-2670 v3 has four DDR4
# channels and tens of GB/s, so if that number is right the kernel is NOT
# memory-bound, which contradicts the assumption the whole optimisation
# story was written around. Before rewriting that story, the compute roof has
# to be measured rather than assumed.
#
# The microbenchmark below is deliberately dumb: eight INDEPENDENT
# multiply-add chains on registers, no memory traffic, no divisions, compiled
# the same way the solver is (-O2, no -march, so plain SSE2 scalar, the same
# instruction set the kernel gets). Eight chains is enough to hide the ~5
# cycle latency of an addsd so the loop measures throughput, not latency.
# The result is consumed and printed so it cannot be optimised away.
#
#SBATCH --job-name=heston_peak
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:10:00

module load compiler-rt/latest
module add mpi/openmpi-x86_64
export PATH=/usr/lib64/openmpi/bin:$PATH
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
echo "### sustained core clock DURING an actual serial solve"
make clean >/dev/null && make heston >/dev/null
./heston --config config/bench.cfg --nt 2000 --maturity 0.00051 \
        --solver serial >/tmp/solve_$$.out 2>/dev/null &
SOLVE=$!
sleep 3
for i in 1 2 3 4 5 6 7 8; do
  grep -m1 'cpu MHz' /proc/cpuinfo | sed 's/^/  /'
  sleep 1
done
wait $SOLVE
echo "  solve result: $(cat /tmp/solve_$$.out)"
rm -f /tmp/peak_$$ /tmp/peak_$$.cpp /tmp/solve_$$.out

# Milestone 2 ceilings. The AVX2 FMA peak on one core and on every core,
# and the memory bandwidth one core and every core can pull.
cat > /tmp/peak_m2_$$.cpp <<'CPP'
#include <immintrin.h>
#include <omp.h>
#include <cstdio>
#include <memory>

// Ten independent FMA chains, enough to hide the five cycle FMA latency on both FMA ports
static double fma_chains(long iters) {
    __m256d a0 = _mm256_set1_pd(1.0), a1 = _mm256_set1_pd(1.1), a2 = _mm256_set1_pd(1.2), a3 = _mm256_set1_pd(1.3), a4 = _mm256_set1_pd(1.4);
    __m256d a5 = _mm256_set1_pd(1.5), a6 = _mm256_set1_pd(1.6), a7 = _mm256_set1_pd(1.7), a8 = _mm256_set1_pd(1.8), a9 = _mm256_set1_pd(1.9);
    __m256d m = _mm256_set1_pd(1.0000000001);
    __m256d c = _mm256_set1_pd(1e-9);
    for (long i = 0; i < iters; ++i) {
        a0 = _mm256_fmadd_pd(a0, m, c); a1 = _mm256_fmadd_pd(a1, m, c); a2 = _mm256_fmadd_pd(a2, m, c); a3 = _mm256_fmadd_pd(a3, m, c); a4 = _mm256_fmadd_pd(a4, m, c);
        a5 = _mm256_fmadd_pd(a5, m, c); a6 = _mm256_fmadd_pd(a6, m, c); a7 = _mm256_fmadd_pd(a7, m, c); a8 = _mm256_fmadd_pd(a8, m, c); a9 = _mm256_fmadd_pd(a9, m, c);
    }
    __m256d sum = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(a0, a1), _mm256_add_pd(a2, a3)), _mm256_add_pd(_mm256_add_pd(a4, a5), _mm256_add_pd(_mm256_add_pd(a6, a7), _mm256_add_pd(a8, a9))));
    double lanes[4];
    _mm256_storeu_pd(lanes, sum);
    return lanes[0] + lanes[1] + lanes[2] + lanes[3];
}

int main() {
    const long iters = 100000000L;
    // Ten chains, four lanes, two flops per FMA
    const double flops_per_thread = (double)iters * 10.0 * 4.0 * 2.0;
    double checksum = 0.0;
    int threads = omp_get_max_threads();
    double t0 = omp_get_wtime();
    #pragma omp parallel reduction(+:checksum)
    checksum += fma_chains(iters);
    double secs = omp_get_wtime() - t0;
    std::printf("checksum %.6f\n", checksum);
    std::printf("avx2_fma_peak_gflops %.3f  threads %d  (%.3f s)\n", flops_per_thread * threads / secs / 1e9, threads, secs);

    // Triad a = b + s * c on three arrays of 256 MiB, the STREAM convention of 24 bytes per element
    const long n = 1L << 25;
    // new double[n] leaves the memory untouched, so the parallel loop below places each page.
    // A std::vector would zero every page on the main thread first, putting it all on socket 0
    std::unique_ptr<double[]> a(new double[n]), b(new double[n]), c(new double[n]);
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < n; ++i) { a[i] = 0.0; b[i] = 1.0; c[i] = 2.0; }
    double best = 1e30;
    for (int rep = 0; rep < 10; ++rep) {
        double start = omp_get_wtime();
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < n; ++i) { a[i] = b[i] + 3.0 * c[i]; }
        double took = omp_get_wtime() - start;
        if (took < best) { best = took; }
    }
    std::printf("check %.1f\n", a[n / 2]);
    std::printf("triad_bandwidth_gbs %.2f  threads %d\n", 24.0 * n / best / 1e9, threads);
    return 0;
}
CPP
g++ -std=c++17 -O3 -mavx2 -mfma -fopenmp -o /tmp/peak_m2_$$ /tmp/peak_m2_$$.cpp
echo
echo "### AVX2 FMA peak and triad bandwidth (-O3 -mavx2 -mfma -fopenmp, OMP_PLACES=cores OMP_PROC_BIND=close)"
for THREADS in 1 12 24; do
  OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close /tmp/peak_m2_$$
done
echo "### the same at 24 threads with spread binding"
OMP_NUM_THREADS=24 OMP_PLACES=cores OMP_PROC_BIND=spread /tmp/peak_m2_$$
rm -f /tmp/peak_m2_$$ /tmp/peak_m2_$$.cpp
} 2>&1 | tee "$OUT"

echo "done: $OUT"
