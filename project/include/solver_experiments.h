// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_SOLVER_EXPERIMENTS_H
#define HESTON_SOLVER_EXPERIMENTS_H

#include <string>

#include "solver.h"

// The OpenMP solver with the scalar row in place of the AVX2 row.
// Shows that the thread speedup and the SIMD speedup multiply together.
class OmpScalarSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_scalar"; }
};

// The OpenMP solver with one parallel region around the whole time loop.
// Tests whether starting the threads once for the whole solve is faster.
class OmpPersistentSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_persistent"; }
};

// The OpenMP solver with the schedule read from OMP_SCHEDULE when it runs.
// One build can try static, dynamic and guided scheduling.
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

// The OpenMP solver with both buffers stored as float.
// Eight cells per instruction and half the memory traffic, with less accuracy.
class OmpFloatSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp_float"; }
};

// Experiment E2. The omp solver on buffers placed at a chosen alignment.
// --align-shift K puts column 1 of every row K doubles past a 64 byte boundary.
// The kernel says which row function runs.
//   loadu    the shipped row, unaligned loads, any K. Only the address moves.
//   load     aligned loads and stores where the address allows, K = 0 or 4.
//   shuffle  one aligned load per row, east and west built in registers, K = 0 or 4.
enum class AlignKernel { loadu, load, shuffle };

class OmpAlignSolver : public Solver {
public:
    // explicit stops a plain AlignKernel turning into a solver by accident.
    explicit OmpAlignSolver(AlignKernel kernel) : kernel_(kernel) {}
    SolveResult solve(const Config& cfg) override;
    std::string name() const override;

private:
    AlignKernel kernel_;
};

// Experiment E9. The mpi solver with the edge row messages overlapped with the interior update.
// The edge rows are updated first and sent, the other rows are updated while the messages travel.
class MpiOverlapSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "mpi_overlap"; }
};

#endif  // HESTON_SOLVER_EXPERIMENTS_H
