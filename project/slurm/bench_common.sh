# Shared setup for the three benchmark jobs. Sourced, never executed directly.
#
# Why three jobs and not one: the cosc3500 QOS caps wall time at 00:15:00
# (`sacctmgr show qos` -> MaxWall). A single job covering the ladder, the
# scaling sweep, the negative controls and the no-vectorise rebuild needs
# roughly forty minutes, so it sits in the queue forever with reason
# QOSMaxWallDurationPerJobLimit. Each script below is sized to finish well
# inside the limit and writes its own CSV; scripts/bench_plot.py takes any
# number of CSVs and merges them on the sweep column.
#
# CSV columns (the two-field prefix is added here so rows from different
# sweeps and builds can never be confused):
#   sweep,build,price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec

REPS=${REPS:-5}
BASE_FLAGS="-std=c++17 -O2 -Wall -Wextra -Iinclude"
NOVEC_FLAGS="$BASE_FLAGS -fno-tree-vectorize"

mkdir -p results
OUT=results/bench_${SWEEP_NAME}_${SLURM_JOB_ID}.csv
ERR=results/bench_${SWEEP_NAME}_${SLURM_JOB_ID}.err

build_and_validate() {
  # $1 = CXXFLAGS to build with. Iron rule (PLAN §6): never benchmark a build
  # that has not just proved it still matches the reference solver.
  make clean
  make CXXFLAGS="$1" all
  make CXXFLAGS="$1" test || {
    echo "FATAL: make test failed for flags '$1' — refusing to benchmark"
    exit 1
  }
}

write_header() {
  # Everything a reader needs to reproduce or challenge these numbers.
  # '#' comments so bench_plot.py can skip them.
  {
    echo "# sweep_name=${SWEEP_NAME}"
    echo "# host=$(hostname)"
    echo "# date=$(date -Is)"
    echo "# jobid=${SLURM_JOB_ID}"
    echo "# compiler=$(g++ --version | head -1)"
    echo "# flags_default='${BASE_FLAGS}'"
    echo "# flags_novec='${NOVEC_FLAGS}'"
    echo "# reps=${REPS}"
    echo "# cpu=$(lscpu | grep -m1 'Model name' | sed 's/.*: *//')"
    lscpu | grep -E '^(L1d|L1i|L2|L3) cache' | sed 's/^/# cache_/'
    echo "# sweep,build,price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec"
  } > "$OUT"
}

run() {
  # run <sweep> <build> <heston args...>  — tags each CSV line with its origin.
  local sweep=$1 build=$2
  shift 2
  echo "[$sweep/$build] $*"
  ./heston "$@" --bench "$REPS" 2>>"$ERR" | sed "s/^/${sweep},${build},/" >> "$OUT"
}
