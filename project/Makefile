# Heston PDE option pricer, serial reference through hybrid MPI and OpenMP.
# Builds on macOS (Homebrew Open MPI over Apple clang) and on rangpur
# (Open MPI 4.1 over g++ 8.5). Benchmark flags stay at -O2.

# One binary for every rung, so it is always built with the MPI wrapper
CXX      = mpicxx
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude

# Use only MPI's C functions. Open MPI 4.1's old C++ wrapper header
# makes g++ 8.5 print warnings under -Wextra, and nothing here uses it
MPIFLAGS = -DOMPI_SKIP_MPICXX

# Apple clang needs Homebrew libomp passed by hand, GCC only needs -fopenmp
ifneq ($(shell $(CXX) --version 2>/dev/null | grep -c clang),0)
  LIBOMP   := $(shell brew --prefix libomp 2>/dev/null)
  OMPFLAGS  = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
  OMPLIBS   = -L$(LIBOMP)/lib -lomp
else
  OMPFLAGS  = -fopenmp
  OMPLIBS   = -fopenmp
endif

# AVX2 and FMA only on x86_64 machines such as rangpur.
# Only the SIMD row file gets them, so the other rungs compile exactly as before
ifeq ($(shell uname -m),x86_64)
  AVXFLAGS = -mavx2 -mfma $(AVXEXTRA)
endif
# Extra flags for the AVX2 file only, empty unless a job sets it.
# Experiment E2b sets -mno-avx256-split-unaligned-load and -store
AVXEXTRA ?=

# GCC 8 keeps std::filesystem in a separate library, newer compilers reject it
FSLIB := $(shell echo 'int main(){}' | $(CXX) -x c++ -std=c++17 - -lstdc++fs \
           -o /dev/null 2>/dev/null && echo -lstdc++fs)

# Objects shared by the main binary and the tests
CORE_OBJS = src/params.o src/grid.o src/stencil.o src/stencil_avx.o \
            src/solver.o src/solver_baseline.o src/solver_serial.o \
            src/solver_simd.o src/solver_omp.o src/solver_mpi.o \
            src/solver_experiments.o src/solver_float.o src/solver_align.o \
            src/black_scholes.o src/io.o

TESTS = tests/test_bs_collapse tests/test_parity tests/test_ladder_matches

default: all
all: heston $(TESTS)

heston: src/main.o $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@ $(OMPLIBS) $(FSLIB)

tests/%: tests/%.cpp $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) $^ -o $@ $(OMPLIBS) $(FSLIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -c $< -o $@

# The one file built with AVX2 and FMA
src/stencil_avx.o: src/stencil_avx.cpp
	$(CXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) $(AVXFLAGS) -c $< -o $@

# Header dependencies, listed by hand
src/main.o: include/params.h include/solver.h include/io.h
src/params.o: include/params.h
src/grid.o: include/grid.h include/params.h
src/stencil.o src/stencil_avx.o: include/stencil.h include/grid.h include/params.h
SOLVER_OBJS = src/solver.o src/solver_baseline.o src/solver_serial.o src/solver_simd.o src/solver_omp.o src/solver_mpi.o src/solver_experiments.o src/solver_float.o src/solver_align.o
$(SOLVER_OBJS): include/solver.h include/grid.h include/params.h
src/solver_baseline.o src/solver_serial.o src/solver_simd.o src/solver_omp.o src/solver_mpi.o: include/io.h
src/solver_mpi.o: include/solver_experiments.h
src/solver_serial.o src/solver_simd.o src/solver_omp.o src/solver_mpi.o src/solver_experiments.o src/solver_float.o src/solver_align.o: include/stencil.h include/io.h
src/solver.o src/solver_experiments.o src/solver_float.o src/solver_align.o: include/solver_experiments.h
src/black_scholes.o: include/black_scholes.h
src/io.o: include/io.h include/grid.h include/params.h include/solver.h

# The ladder test runs again at 1, 2 and 7 threads.
# 7 does not divide the row count, so the bands come out uneven
test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done
	@for n in 1 2 7; do echo "== tests/test_ladder_matches with OMP_NUM_THREADS=$$n"; OMP_NUM_THREADS=$$n ./tests/test_ladder_matches 2>/dev/null || exit 1; done
	@echo "ALL TESTS PASS"

# The ladder test again as 2 and 3 processes.
# 3 does not divide the test grid's 32 rows, so the row blocks come out uneven.
# Two threads per process, and no core pinning so the threads can spread out
MPIEXEC = mpiexec --oversubscribe --bind-to none
test-mpi: tests/test_ladder_matches
	@for n in 2 3; do echo "== tests/test_ladder_matches on $$n processes"; OMP_NUM_THREADS=2 $(MPIEXEC) -n $$n ./tests/test_ladder_matches 2>/dev/null || exit 1; done
	@echo "MPI TESTS PASS"

# Checks the AVX2 row against the scalar row on random data.
# On an Apple silicon Mac it is built for x86_64 and runs under Rosetta
ifeq ($(shell uname -s)-$(shell uname -m),Darwin-arm64)
  AVXCHECK_CXX = clang++ -arch x86_64
else
  AVXCHECK_CXX = $(CXX)
endif
AVXCHECK_FLAGS = -std=c++17 -O2 -Wall -Wextra -Iinclude

avx-check:
	$(AVXCHECK_CXX) $(AVXCHECK_FLAGS) -mavx2 -mfma -c src/stencil_avx.cpp -o tests/stencil_avx_check.o
	$(AVXCHECK_CXX) $(AVXCHECK_FLAGS) tests/test_avx_row.cpp src/stencil.cpp src/grid.cpp src/params.cpp tests/stencil_avx_check.o -o tests/test_avx_row
	./tests/test_avx_row

smoke: heston
	./heston --config config/smoke.cfg

# ---- local-only figures and animations (need scripts/requirements.txt) ----
# This Mac's system python mixes arm64 and x86_64 wheels, so the venv is
# pinned to arm64. Rebuild it with:
#   arch -arm64 /usr/bin/python3 -m venv .venv && \
#   arch -arm64 .venv/bin/python -m pip install -r scripts/requirements.txt
PY = arch -arm64 .venv/bin/python

# The weather map, rendered so the motion is actually visible (STUDY_GUIDE
# P6b): time value, cropped to the reachable region, sqrt-warped frames.
demo-anim: heston
	./heston --config config/demo.cfg --dump-every 60 --dump-dir results/demo
	$(PY) scripts/weather_map.py results/demo --config config/demo.cfg \
	  --time-value --zoom --gamma 0.5 --cap 99.9 --warp sqrt --frames 140 \
	  --fps 12 --envelope --price-trace --gif --out results/demo_timevalue.mp4
	$(PY) scripts/weather_map.py results/demo --config config/demo.cfg \
	  --style surface --time-value --zoom --gamma 0.6 --cap 99.9 \
	  --warp sqrt --frames 110 --fps 12 --out results/demo_surface.mp4

# "Heston becomes Black-Scholes", twice: the xi=0 collapse against the
# closed form, and the xi sweep flattening the implied-vol smile.
collapse-anim: heston
	$(PY) scripts/bs_collapse_anim.py --work results/collapse_work \
	  --out results/bs_collapse.mp4 --frames 100 --refine-check --gif
	$(PY) scripts/xi_sweep_smile.py --work results/xi_work --steps 30 \
	  --verify 5 --out results/xi_smile.mp4 --gif

# Third opinion on the price: random paths (MC) and Fourier inversion, both
# independent of the grid. Answers "is the noise modelled correctly?".
mc-check: heston
	./heston --config config/demo.cfg
	$(PY) scripts/monte_carlo_check.py --config config/demo.cfg \
	  --paths 400000 --steps 800 --pde-price 196.1054098

# The grid-refinement study (P4). Every other figure had a target and this one
# did not, which is how its ladder drifted out of alignment with reference.cfg
# without anyone re-running it. ~6 minutes for three rungs.
convergence: heston
	$(PY) scripts/convergence_plot.py --binary ./heston \
	  --config config/reference.cfg --out results/convergence.png \
	  --max-rungs 3

# THE GATE THAT COVERS xi AND rho. `make test` does not: test_bs_collapse sets
# xi=0 (killing both the cross term and the vol-of-vol diffusion), test_parity
# checks an identity that holds in every model and is measurably blind to both,
# and test_ladder_matches only compares solvers with each other. This target
# prices the full parameter set against Fourier inversion and FAILS if they
# disagree. It is Python + numpy, so it cannot live inside `make test` (which
# has to run on rangpur with no venv) — run it locally after any change to the
# stencil. ~40 s.
check-heston: heston
	@echo "=== leg 1: node-aligned grid (526x101) ==="
	@price=$$(./heston --config config/demo.cfg 2>/dev/null | cut -d, -f1); \
	echo "PDE price: $$price"; \
	$(PY) scripts/monte_carlo_check.py --config config/demo.cfg \
	  --paths 40000 --steps 400 --no-delta --pde-price $$price --gate 2e-3
	@echo
	@echo "=== leg 2: deliberately OFF-node grid (527x100) ==="
	@echo "Every other gate uses a node-aligned grid, where extract_result's"
	@echo "blend weights are exactly 0 and the bilinear path is never taken."
	@echo "This leg puts spot 130.25 cells and v0 19.8 cells along, so the"
	@echo "reported price IS the interpolated one, and checks it against the"
	@echo "same external analytic answer. That path priced the headline number"
	@echo "and had no absolute validation before."
	@price=$$(./heston --config config/demo.cfg --ns 527 --nv 100 2>/dev/null | cut -d, -f1); \
	echo "PDE price (interpolated): $$price"; \
	$(PY) scripts/monte_carlo_check.py --config config/demo.cfg \
	  --paths 40000 --steps 400 --no-delta --pde-price $$price --gate 2e-3

# ---- P7 optimisation figures ----------------------------------------------
# Benchmark figures from the rangpur CSVs. Run `make fetch` first — these
# scripts exit non-zero rather than plot laptop numbers as results.
bench-plot:
	$(PY) scripts/bench_plot.py results/bench_*.csv --out results/

# The Milestone 2 figures F1 to F10, from the final rangpur sweeps (denormals flushed).
# Every CSV is named by its job id, so the figures can always be traced to a run
M2_CSVS = results/bench_m2_ladder_651295.csv results/bench_m2_strong_omp_651296.csv \
          results/bench_m2_strong_omp_scalar_651297.csv results/bench_m2_weak_651298.csv \
          results/bench_m2_size_656959.csv results/bench_m2_exp_a_651300.csv \
          results/bench_m2_exp_b_651301.csv results/bench_m2_exp_c_651474.csv \
          results/bench_m2_exp_float_652647.csv results/bench_m2_optlevels_652696.csv \
          results/bench_m2_mpi_651303.csv results/bench_m2_denormal_651294.csv
m2-figures:
	$(PY) scripts/m2_figures.py $(M2_CSVS) --peak results/peak_651302.txt \
	  --test-log results/logs/slurm-652647.out --out results/m2 \
	  --profile results/profile_baseline_652697.txt results/profile_serial_652697.txt results/profile_simd_652697.txt
	$(PY) scripts/decomposition_diagram.py --out results/m2/f10_decomposition.png
	$(PY) scripts/method_figures.py --out results/m2

# Why layout and loop order dominate: the cache-line sweep (row-major vs the
# ctl-order control, with a modelled LRU cache) and the double-buffer
# ping-pong. --bench-csv makes the caption quote the MEASURED slowdown, so it
# can never drift from the number in the table.
memory-anim:
	$(PY) scripts/memory_anim.py --style cache --gif \
	  --bench-csv results/bench_*.csv --out results/memory_cache.mp4
	$(PY) scripts/memory_anim.py --style buffers \
	  --out results/memory_buffers.mp4

# The figure that TESTS the "memory-bound" claim: flops and bytes counted term
# by term from the level-6 kernel, achieved point from measured throughput.
# --peak-flops is not a datasheet number — it is what slurm/bench_peak.sh
# measured on r730-2 (5.117 GFLOP/s, scalar SSE2 double, eight independent
# multiply-add chains, the same -O2 no -march build the solver gets).
roofline:
	$(PY) scripts/roofline.py results/bench_*.csv --peak-flops 5.117 \
	  --out results/roofline.png

# ---- presentation assets ---------------------------------------------------
# Slide-ready screenshots of the load-bearing code (needs pygments + pillow,
# see scripts/requirements.txt) and the Introduction beat's payoff diagram.
# Shot specs and where each lands in the video: presentation/TIMELINE.md.
code-shots:
	$(PY) scripts/code_shots.py --out results/m2/code_shots

payoff:
	$(PY) scripts/payoff_diagram.py --out results/payoff.png

# The Beamer deck. A PDF cannot play the MP4 animations, so those slides
# carry a still frame (stills) and the clip itself is spliced in during the
# edit, per presentation/SHOTLIST.md. tectonic fetches beamer and the fonts
# on first run (brew install tectonic if missing).
stills:
	$(PY) scripts/slide_stills.py --out presentation/stills

slides: stills
	tectonic presentation/deck.tex

# The Milestone 2 deck, one frame per on-screen item in presentation/M2_SCRIPT.md.
slides-m2:
	tectonic presentation/deck_m2.tex

clean:
	rm -f heston $(TESTS) src/*.o tests/*.o tests/test_avx_row

# ---- submission packaging --------------------------------------------------
# The rubric wants ONE zip named Milestone2_<student>.zip holding the code and
# the video, and the video must be H.264 under 100 MB.
#
# This target exists because a plain `zip -r` of this directory would ship
# project/.env, which holds the rangpur password, and would also carry ~570 MB
# of venv and rendered frames. The exclusions below are the point of the
# target, not incidental tidiness.
STUDENT = 48008361
PKG     = Milestone2_$(STUDENT)
VIDEO   = $(PKG).mp4

# The video check comes BEFORE clean, so a missing video does not leave the
# tree unbuilt as a parting gift.
package:
	@test -f "$(VIDEO)" || { \
	  echo "MISSING: $(VIDEO) — record and export it first (H.264, <100 MB)."; \
	  echo "See presentation/PRESENTATION.md for the beat sheet and the format gate."; \
	  exit 1; }
	$(MAKE) clean
	@rm -rf /tmp/$(PKG) && mkdir -p /tmp/$(PKG)
	rsync -a --exclude '.env' --exclude '.venv' --exclude '.git' \
	         --exclude 'graphify-out' --exclude '__pycache__' \
	         --exclude '.DS_Store' --exclude '*.o' --exclude 'heston' \
	         --exclude 'results/demo' --exclude 'results/blowup' \
	         --exclude 'results/smooth' --exclude 'results/smoke_dump' \
	         --exclude 'results/collapse_work' --exclude 'results/xi_work' \
	         ./ /tmp/$(PKG)/code/
	cp "$(VIDEO)" /tmp/$(PKG)/
	@echo "--- verifying no secrets made it in ---"
	@! grep -rqI --exclude-dir=.git '[P]ASSWORD\|[p]assword *=' /tmp/$(PKG)/code/ \
	  || { echo "REFUSING: something password-shaped is in the payload"; exit 1; }
	@test ! -e /tmp/$(PKG)/code/.env || { echo "REFUSING: .env present"; exit 1; }
	cd /tmp && rm -f $(PKG).zip && zip -qr $(PKG).zip $(PKG)
	@mv /tmp/$(PKG).zip .
	@rm -rf /tmp/$(PKG)
	@echo "wrote $(PKG).zip  ($$(du -h $(PKG).zip | cut -f1))"
	@echo "Check it is under Blackboard's 100 MB limit before uploading."

# ---- rangpur targets (ssh alias 'rangpur', see ~/CLAUDE.md workflow) ----
REMOTE_DIR = cosc3500/project

# --delete never touches excluded paths, so job logs and built tests on rangpur survive.
# The test sources are included first, so only the test binaries are excluded
sync:
	rsync -az --delete \
	  --exclude '.git' --exclude 'graphify-out' --exclude 'results' \
	  --exclude '*.o' --exclude 'heston' --exclude '.env' \
	  --exclude 'slurm-*.out' --include 'tests/test_*.cpp' --exclude 'tests/test_*' \
	  ./ rangpur:$(REMOTE_DIR)/

# sbatch/squeue live in /opt/slurm/bin, which rangpur only puts on PATH for an
# interactive LOGIN shell — a plain `ssh host cmd` runs non-interactively and
# reports "sbatch: command not found". Prepend the directory explicitly.
SLURM_BIN = /opt/slurm/bin

remote-debug: sync
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/debug.sh'

# The node's measured compute ceiling for the roofline
remote-peak: sync
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/bench_peak.sh'

queue:
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; squeue -u $$USER'

fetch:
	mkdir -p results/logs
	rsync -az rangpur:$(REMOTE_DIR)/results/ results/
	rsync -az 'rangpur:$(REMOTE_DIR)/slurm-*.out' results/logs/ || true

.PHONY: default all test test-mpi avx-check smoke clean sync remote-debug remote-peak \
        queue fetch demo-anim collapse-anim mc-check check-heston \
        package \
        convergence bench-plot memory-anim roofline m2-figures slides-m2
