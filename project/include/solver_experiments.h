// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_SOLVER_EXPERIMENTS_H
#define HESTON_SOLVER_EXPERIMENTS_H

#include <string>

#include "solver.h"

// The OpenMP solver with the scalar row instead of the AVX2 row.
// Used to show that threads and SIMD speedups multiply.
class OmpScalarSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_scalar"; }
};

// The OpenMP solver with one parallel region around the whole time loop.
// Tests whether starting the threads once instead of every step is faster.
class OmpPersistentSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_persistent"; }
};

// The OpenMP solver with the schedule read from OMP_SCHEDULE at run time.
// Lets one build try static, dynamic and guided scheduling.
class OmpRuntimeScheduleSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_runtime"; }
};

// The OpenMP solver with the v = 0 row and the stock edges done by the threads.
// Tests whether the serial boundary work is big enough to matter.
class OmpParallelEdgesSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_edges"; }
};

// The OpenMP solver with both sheets in float instead of double.
// Trades accuracy for eight cells per instruction and half the memory traffic.
class OmpFloatSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_float"; }
};

#endif  // HESTON_SOLVER_EXPERIMENTS_H
