#!/bin/bash -l
#
# Diagnostic, not a result. rangpur measured 4.1e7 cell-updates/s on
# config/bench.cfg — ten times slower than a laptop, which is not a plausible
# hardware gap for this kernel. The suspect is DENORMAL arithmetic: bench.cfg
# uses a very short maturity, so most of the sub-strike half of the grid holds
# values around 1e-300, and x86 takes a microcode assist (tens of cycles) on
# every operation that consumes or produces a denormal.
#
# The test: -ffast-math makes GCC link crtfastmath.o, which sets FTZ/DAZ so
# denormals are flushed to zero in hardware. If throughput jumps by ~10x, the
# hypothesis is confirmed and bench.cfg needs a maturity that keeps the value
# field in normal range. Nothing here is quoted as a benchmark number.
#
#SBATCH --job-name=heston_diag
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
BASE="-std=c++17 -O2 -Wall -Wextra -Iinclude"

echo "### CPU"
lscpu | grep -E 'Model name|^CPU\(s\)|MHz|L1d|L2|L3'

echo
echo "### A: default -O2, bench.cfg short maturity (0.00051)"
make clean >/dev/null && make CXXFLAGS="$BASE" heston >/dev/null
./heston --config config/bench.cfg --nt 200 --maturity 0.0000512 --bench 2

echo
echo "### B: -O2 -ffast-math (sets FTZ/DAZ), same config"
make clean >/dev/null && make CXXFLAGS="$BASE -ffast-math" heston >/dev/null
./heston --config config/bench.cfg --nt 200 --maturity 0.0000512 --bench 2

echo
echo "### C: default -O2, LONG maturity so values stay normal (nt scaled)"
make clean >/dev/null && make CXXFLAGS="$BASE" heston >/dev/null
./heston --config config/bench.cfg --nt 200 --maturity 0.25 --bench 2

echo
echo "### D: default -O2, moderate maturity 0.01"
./heston --config config/bench.cfg --nt 200 --maturity 0.01 --bench 2
