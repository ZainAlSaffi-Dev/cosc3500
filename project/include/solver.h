// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_SOLVER_H
#define HESTON_SOLVER_H

#include <memory>
#include <string>

#include "grid.h"
#include "params.h"

// What a solver hands back after pricing the option.
// Holds the price, the three Greeks and the timing.
struct SolveResult {
    // Option price today.
    double price = 0.0;
    // How much the price moves when the stock moves.
    double delta = 0.0;
    // How much delta moves when the stock moves.
    double gamma = 0.0;
    // How much the price moves when the variance moves.
    double vega = 0.0;
    // Seconds spent inside the time loop.
    double seconds = 0.0;
    // Cells updated per second.
    double cell_updates_per_sec = 0.0;
    // Largest stable timestep.
    double dt_stable_estimate = 0.0;
    // OpenMP threads per process.
    int threads = 1;
    // MPI processes.
    int ranks = 1;
    // Seconds the first process spent exchanging edge rows with its neighbours. Included in seconds.
    double exchange_seconds = 0.0;
    // Seconds spent on the serial boundary work. Included in seconds.
    double boundary_seconds = 0.0;
};

// The common shape every solver follows, like an abstract base class in Python.
// Each solver writes its own solve() and name() and shares the helpers below.
class Solver {
public:
    // Lets a subclass clean up when deleted through a Solver pointer.
    virtual ~Solver() = default;
    // virtual means the subclass's version runs.
    // The = 0 means every subclass must write its own version.
    virtual SolveResult solve(const Config& cfg) = 0;
    virtual std::string name() const = 0;

protected:
    // Helpers shared by every subclass, written in solver.cpp.
    void extract_result(const Grid& grid, const Config& cfg, SolveResult& out) const;
    void report_stability(const Grid& grid, const Config& cfg, double dt, SolveResult& out) const;
};

// The first solver, with the whole PDE written out inside the loop.
// Every other solver is checked against its answers.
// The : public Solver part makes it a subclass of Solver.
class BaselineSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    // override marks this as replacing the Solver version.
    std::string name() const override { return "baseline"; }
};

// The single core solver with the constants worked out before the time loop.
// Gives the same answers as the baseline with less work per cell.
class SerialSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "serial"; }
};

// The single core solver that updates four cells per instruction with AVX2.
// Same time loop as SerialSolver with a different row update.
class SimdSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "simd"; }
};

// The solver that splits the rows between the cores of one machine with OpenMP.
// Same time loop as SimdSolver with one pragma above the row loop.
class OmpSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp"; }
};

// The solver that splits the rows between MPI processes, with threads inside each.
// Each process owns a block of rows and exchanges edge rows with its neighbours every step.
class MpiSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "mpi"; }
};

// Creates the solver whose name matches, for example "serial".
std::unique_ptr<Solver> make_solver(const std::string& name);

#endif  // HESTON_SOLVER_H
