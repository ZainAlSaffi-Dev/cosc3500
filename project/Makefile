# COSC3500 Milestone 1 — serial Heston PDE pricer.
# Local dev on macOS (clang++ via g++ alias) and rangpur (g++ 8.x) — plain
# C++17, no external deps. Benchmark flags fixed at -O2 (PLAN §1); report
# any change in every benchmark.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude

# GCC 8 (rangpur ships 8.5.0) keeps std::filesystem in a SEPARATE library and
# fails to link io.cpp without it; clang/libc++ on the Mac and GCC 9+ have it
# in the main runtime and reject the flag. Probe once per make run rather than
# hardcoding either, so the same Makefile builds in both places.
FSLIB := $(shell echo 'int main(){}' | $(CXX) -x c++ -std=c++17 - -lstdc++fs \
           -o /dev/null 2>/dev/null && echo -lstdc++fs)

# Objects shared by the main binary and the tests (everything except main.o)
CORE_OBJS = src/params.o src/grid.o src/solver.o src/solver_baseline.o \
            src/solver_opt.o src/solver_opt_kernels.o src/black_scholes.o \
            src/io.o

TESTS = tests/test_bs_collapse tests/test_parity tests/test_opt_matches

default: all
all: heston $(TESTS)

heston: src/main.o $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(FSLIB)

tests/%: tests/%.cpp $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(FSLIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Header deps (small project — list by hand, keep honest)
src/main.o: include/params.h include/solver.h include/io.h
src/params.o: include/params.h
src/grid.o: include/grid.h include/params.h
src/solver.o src/solver_baseline.o src/solver_opt.o: include/solver.h include/grid.h include/params.h
src/solver.o: include/io.h
src/solver_opt.o src/solver_opt_kernels.o: include/solver_opt_kernels.h
src/black_scholes.o: include/black_scholes.h
src/io.o: include/io.h include/grid.h include/params.h include/solver.h

test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done
	@echo "ALL TESTS PASS"

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
# and test_opt_matches only compares the solver against itself. This target
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
	$(PY) scripts/code_shots.py --out results/code_shots

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

clean:
	rm -f heston $(TESTS) src/*.o

# ---- submission packaging --------------------------------------------------
# The rubric wants ONE zip named Milestone1_<student>.zip holding the code and
# the video, and the video must be H.264 under 100 MB.
#
# This target exists because a plain `zip -r` of this directory would ship
# project/.env, which holds the rangpur password, and would also carry ~570 MB
# of venv and rendered frames. The exclusions below are the point of the
# target, not incidental tidiness.
STUDENT = 48008361
PKG     = Milestone1_$(STUDENT)
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
	@! grep -rqI --exclude-dir=.git 'PASSWORD\|password *=' /tmp/$(PKG)/code/ \
	  || { echo "REFUSING: something password-shaped is in the payload"; exit 1; }
	@test ! -e /tmp/$(PKG)/code/.env || { echo "REFUSING: .env present"; exit 1; }
	cd /tmp && rm -f $(PKG).zip && zip -qr $(PKG).zip $(PKG)
	@mv /tmp/$(PKG).zip .
	@rm -rf /tmp/$(PKG)
	@echo "wrote $(PKG).zip  ($$(du -h $(PKG).zip | cut -f1))"
	@echo "Check it is under Blackboard's 100 MB limit before uploading."

# ---- rangpur targets (ssh alias 'rangpur', see ~/CLAUDE.md workflow) ----
REMOTE_DIR = cosc3500/project

sync:
	rsync -az --delete \
	  --exclude '.git' --exclude 'graphify-out' --exclude 'results' \
	  --exclude '*.o' --exclude 'heston' --exclude '.env' \
	  ./ rangpur:$(REMOTE_DIR)/

# sbatch/squeue live in /opt/slurm/bin, which rangpur only puts on PATH for an
# interactive LOGIN shell — a plain `ssh host cmd` runs non-interactively and
# reports "sbatch: command not found". Prepend the directory explicitly.
SLURM_BIN = /opt/slurm/bin

remote-debug: sync
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/debug_serial.sh'

# Several jobs, not one: the cosc3500 QOS caps wall time at 15 minutes and each
# job spends ~5 of those on the mandatory build-and-validate. Each writes its
# own CSV; bench_plot.py and roofline.py merge them on the sweep column.
#
# All five sweeps, in one target. bench_pages.sh (sweep E) used to be a manual
# step, which meant anyone reproducing the repo got sweep A's 1.81x and never
# found the 2.56x headline. RUNNING.md §3 says what each job measures.
remote-bench: sync
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/bench_ladder.sh'
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/bench_scaling.sh'
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/bench_controls.sh'
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/bench_novec.sh'
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/bench_pages.sh'

# The roofline's ceiling and the two diagnostics that explained the surprising
# numbers. Separate from remote-bench because these are measured once and then
# cited, rather than re-run every time the ladder changes.
remote-diag: sync
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/bench_peak.sh'
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/diag_alloc.sh'
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; cd $(REMOTE_DIR) && sbatch slurm/diag_denormal.sh'

queue:
	ssh rangpur 'export PATH=$(SLURM_BIN):$$PATH; squeue -u $$USER'

fetch:
	mkdir -p results/logs
	rsync -az rangpur:$(REMOTE_DIR)/results/ results/
	rsync -az 'rangpur:$(REMOTE_DIR)/slurm-*.out' results/logs/ || true

.PHONY: default all test smoke clean sync remote-debug remote-bench \
        remote-diag queue fetch demo-anim collapse-anim mc-check check-heston \
        package \
        convergence bench-plot memory-anim roofline
