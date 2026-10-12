#!/bin/bash -l
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Five minute check that the code builds and passes the tests on a compute node.
# Eight cores so the OpenMP tests at 1, 2 and 7 threads run on real cores.
#
#SBATCH --job-name=heston_debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --time=0-00:05:00

# Compiler runtime and Open MPI, which provides mpicxx
module load compiler-rt/latest
module add mpi/openmpi-x86_64
# The module does not put the MPI wrappers on PATH by itself
export PATH=/usr/lib64/openmpi/bin:$PATH
# Skip the network transports rangpur's Open MPI cannot start
export OMPI_MCA_mtl=^ofi,psm,psm2,psm3
export OMPI_MCA_btl=^uct,ofi

# Build in a private copy of the tree, so a running benchmark is never disturbed
WORK_DIR=$HOME/cosc3500/work/job_${SLURM_JOB_ID}
mkdir -p "$WORK_DIR"
rsync -a --exclude results --exclude '*.o' --exclude heston --exclude 'slurm-*.out' ./ "$WORK_DIR/"
cd "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT

# Fresh build so no stale object files are reused
make clean
make all
# Record where this ran and with what compiler
hostname
mpicxx --version | head -1

# One tiny solve per rung, then the full test suite
./heston --config config/smoke.cfg --solver baseline
./heston --config config/smoke.cfg --solver serial
./heston --config config/smoke.cfg --solver simd
OMP_NUM_THREADS=8 ./heston --config config/smoke.cfg --solver omp
OMP_NUM_THREADS=4 mpiexec --oversubscribe --bind-to none -n 2 ./heston --config config/smoke.cfg --solver mpi
make test
make test-mpi

# The AVX2 row against the scalar row on random data
make avx-check

# zmm registers are AVX-512, which this Xeon does not have.
# Some AVX-512 intrinsics also start with _mm256_, so check none slipped in
zmm_count=$(objdump -d heston | grep -c zmm)
echo "zmm instructions in heston: $zmm_count"
# ymm registers are the 256 bit AVX2 ones, so this should not be zero
ymm_count=$(objdump -d src/stencil_avx.o | grep -c ymm)
echo "ymm instructions in stencil_avx.o: $ymm_count"
