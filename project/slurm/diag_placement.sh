#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Prints which CPU, and so which socket, each OpenMP thread runs on for each binding.
# On r730 nodes the even CPUs are socket 0 and the odd CPUs are socket 1.
#
#SBATCH --job-name=heston_diag_placement
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --exclusive
#SBATCH --nodelist=r730-2
#SBATCH --time=0-00:03:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

module load compiler-rt/latest
cat > /tmp/where_$$.cpp <<'CPP'
#include <omp.h>
#include <sched.h>
#include <cstdio>
int main() {
    int cpu_of[64];
    #pragma omp parallel
    cpu_of[omp_get_thread_num()] = sched_getcpu();
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        std::printf("%d%s", cpu_of[t], t + 1 < omp_get_max_threads() ? " " : "\n");
    }
}
CPP
g++ -O2 -fopenmp /tmp/where_$$.cpp -o /tmp/where_$$
for BIND in close spread; do
  for THREADS in 6 12 24; do
    echo -n "OMP_PLACES=cores OMP_PROC_BIND=$BIND threads $THREADS -> cpus: "
    OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=$BIND /tmp/where_$$
  done
done
rm -f /tmp/where_$$ /tmp/where_$$.cpp
