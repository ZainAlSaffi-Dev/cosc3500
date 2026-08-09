# COSC3500 Milestone 1 — serial Heston PDE pricer.
# Local dev on macOS (clang++ via g++ alias) and rangpur (g++ 8.x) — plain
# C++17, no external deps. Benchmark flags fixed at -O2 (PLAN §1); report
# any change in every benchmark.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude

# Objects shared by the main binary and the tests (everything except main.o)
CORE_OBJS = src/params.o src/grid.o src/solver_baseline.o src/solver_opt.o \
            src/solver_opt_kernels.o src/black_scholes.o src/io.o

TESTS = tests/test_bs_collapse tests/test_parity tests/test_opt_matches

default: all
all: heston $(TESTS)

heston: src/main.o $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

tests/%: tests/%.cpp $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Header deps (small project — list by hand, keep honest)
src/main.o: include/params.h include/solver.h include/io.h
src/params.o: include/params.h
src/grid.o: include/grid.h include/params.h
src/solver_baseline.o src/solver_opt.o: include/solver.h include/grid.h include/params.h
src/solver_opt.o src/solver_opt_kernels.o: include/solver_opt_kernels.h
src/black_scholes.o: include/black_scholes.h
src/io.o: include/io.h include/grid.h include/params.h include/solver.h

test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done
	@echo "ALL TESTS PASS"

smoke: heston
	./heston --config config/smoke.cfg

clean:
	rm -f heston $(TESTS) src/*.o

# ---- rangpur targets (ssh alias 'rangpur', see ~/CLAUDE.md workflow) ----
REMOTE_DIR = cosc3500/project

sync:
	rsync -az --delete \
	  --exclude '.git' --exclude 'graphify-out' --exclude 'results' \
	  --exclude '*.o' --exclude 'heston' --exclude '.env' \
	  ./ rangpur:$(REMOTE_DIR)/

remote-debug: sync
	ssh rangpur 'cd $(REMOTE_DIR) && sbatch slurm/debug_serial.sh'

remote-bench: sync
	ssh rangpur 'cd $(REMOTE_DIR) && sbatch slurm/bench_serial.sh'

queue:
	ssh rangpur 'squeue -u $$USER'

fetch:
	rsync -az rangpur:$(REMOTE_DIR)/results/ results/
	rsync -az 'rangpur:$(REMOTE_DIR)/slurm-*.out' results/logs/ || true

.PHONY: default all test smoke clean sync remote-debug remote-bench queue fetch
