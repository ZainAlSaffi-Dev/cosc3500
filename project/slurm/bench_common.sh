# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Shared setup for every benchmark job. Sourced by each bench script, never run on its own.
# Each script sets SWEEP_NAME first, then calls build_and_validate, write_header and run.
#
# CSV columns, the first two are added here so rows from different sweeps never mix up
#   sweep,label,price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec,threads,ranks,exchange_seconds,boundary_seconds

# Compiler runtime and Open MPI, which provides mpicxx and mpiexec
module load compiler-rt/latest
module add mpi/openmpi-x86_64
# The module does not put the MPI wrappers on PATH by itself
export PATH=/usr/lib64/openmpi/bin:$PATH
# Skip the network transports rangpur's Open MPI cannot start
export OMPI_MCA_mtl=^ofi,psm,psm2,psm3
export OMPI_MCA_btl=^uct,ofi

# Pin the page regime to the Milestone 1 headline.
# Buffers over 1 MiB come straight from mmap, so every rep gets 2 MiB huge pages
export MALLOC_MMAP_THRESHOLD_=1048576

# One thread per physical core, filling neighbouring cores first
export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}

# Repetitions per measurement, the plots take the median
REPS=${REPS:-5}

# Build and run in a private copy of the tree, so two jobs never share object files.
# It sits in the home directory, which every node can see, and results go back to the submitted tree
SUBMIT_DIR=$(pwd)
WORK_DIR=$HOME/cosc3500/work/job_${SLURM_JOB_ID}
mkdir -p "$WORK_DIR" "$SUBMIT_DIR/results"
rsync -a --exclude results --exclude '*.o' --exclude heston --exclude 'slurm-*.out' ./ "$WORK_DIR/"
cd "$WORK_DIR"
# Remove the copy when the job ends, however it ends
trap 'rm -rf "$WORK_DIR"' EXIT

# Pin the floating point regime too. Tiny values below about 1e-308 (denormals)
# take a slow path in the CPU, and how many there are depends on the grid and nt.
# The preload library sets flush to zero, set FLUSH_DENORMALS=0 to run without it
FLUSH_DENORMALS=${FLUSH_DENORMALS:-1}
gcc -O2 -shared -fPIC slurm/flush_denormals.c -o "$WORK_DIR/flush_denormals.so"
FTZ_LIB=$WORK_DIR/flush_denormals.so

# Prints the LD_PRELOAD setting for one run, empty when flushing is off
preload_setting() {
  if [ "$FLUSH_DENORMALS" = "1" ]; then
    echo "LD_PRELOAD=$FTZ_LIB"
  else
    echo "LD_PRELOAD="
  fi
}

OUT=$SUBMIT_DIR/results/bench_${SWEEP_NAME}_${SLURM_JOB_ID}.csv
ERR=$SUBMIT_DIR/results/bench_${SWEEP_NAME}_${SLURM_JOB_ID}.err

# Builds from clean and runs every test.
# Refuses to benchmark anything if a single test fails
build_and_validate() {
  make clean
  make all || { echo "FATAL: build failed, refusing to benchmark"; exit 1; }
  make test || { echo "FATAL: make test failed, refusing to benchmark"; exit 1; }
  make test-mpi || { echo "FATAL: make test-mpi failed, refusing to benchmark"; exit 1; }
  make avx-check || { echo "FATAL: make avx-check failed, refusing to benchmark"; exit 1; }
}

# Writes everything a reader needs to reproduce or challenge the numbers.
# Lines start with # so the plotting script skips them
write_header() {
  {
    echo "# sweep_name=${SWEEP_NAME}"
    echo "# host=$(hostname)"
    echo "# nodes=${SLURM_JOB_NODELIST}"
    echo "# date=$(date -Is)"
    echo "# jobid=${SLURM_JOB_ID}"
    echo "# compiler=$(mpicxx --version | head -1)"
    echo "# mpi=$(mpiexec --version | head -1)"
    echo "# flags=$(make -s -p -n 2>/dev/null | grep -m1 '^CXXFLAGS' | sed 's/^CXXFLAGS *[?:]*= *//')"
    echo "# avx_flags=$(make -s -p -n 2>/dev/null | grep -m1 '^AVXFLAGS' | sed 's/^AVXFLAGS *[?:]*= *//')"
    echo "# reps=${REPS}"
    echo "# page_regime=MALLOC_MMAP_THRESHOLD_=${MALLOC_MMAP_THRESHOLD_}, THP $(cat /sys/kernel/mm/transparent_hugepage/enabled)"
    echo "# omp_places=${OMP_PLACES} omp_proc_bind=${OMP_PROC_BIND}"
    echo "# flush_denormals=${FLUSH_DENORMALS} (FTZ and DAZ set by slurm/flush_denormals.c through LD_PRELOAD, unless a row label says otherwise)"
    echo "# cpu=$(lscpu | grep -m1 'Model name' | sed 's/.*: *//')"
    lscpu | grep -E '^(Socket|Core|Thread|NUMA node)' | sed 's/^/# /'
    lscpu | grep -E '^(L1d|L1i|L2|L3) cache' | sed 's/^/# cache_/'
    echo "# sweep,label,price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec,threads,ranks,exchange_seconds,boundary_seconds"
  } > "$OUT"
}

# Prints the maturity that puts dt at 70 percent of the stability limit.
# Arguments are ns nv vmax nt, the formula is the one report_stability uses.
# xi and rate below must match config/bench.cfg
maturity_for() {
  awk -v ns="$1" -v nv="$2" -v vmax="$3" -v nt="$4" 'BEGIN {
    xi = 0.35
    rate = 0.045
    dt_stable = 1.0 / (vmax * (ns - 1) * (ns - 1) + xi * xi * (nv - 1) * (nv - 1) / vmax + rate)
    printf "%.6e", 0.7 * dt_stable * nt
  }'
}

# Runs heston REPS times and tags each CSV line with the sweep and a label.
# Usage is run <sweep> <label> <heston arguments>.
# RUN_PREFIX, if set, goes in front of heston, for example numactl --interleave=all
run() {
  local sweep=$1
  local label=$2
  shift 2
  echo "[$sweep/$label] $(preload_setting) ${RUN_PREFIX:-} ./heston $*"
  env "$(preload_setting)" ${RUN_PREFIX:-} ./heston "$@" --bench "$REPS" 2>>"$ERR" | sed "s/^/${sweep},${label},/" >> "$OUT"
}

# The same as run, but starts heston under mpiexec.
# Usage is run_mpi <sweep> <label> "<mpiexec arguments>" <heston arguments>
run_mpi() {
  local sweep=$1
  local label=$2
  local launcher=$3
  shift 3
  echo "[$sweep/$label] mpiexec $launcher -x $(preload_setting) ./heston $*"
  mpiexec $launcher -x MALLOC_MMAP_THRESHOLD_ -x OMP_PLACES -x OMP_PROC_BIND -x "$(preload_setting)" ./heston "$@" --bench "$REPS" 2>>"$ERR" | sed "s/^/${sweep},${label},/" >> "$OUT"
}
