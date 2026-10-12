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
// The price, the three Greeks, and how long the solve took.
struct SolveResult {
    // Option price today
    double price = 0.0;
    // How much the price moves when the stock moves
    double delta = 0.0;
    // How much delta moves when the stock moves
    double gamma = 0.0;
    // How much the price moves when the variance moves
    double vega = 0.0;
    // Time spent inside the time loop
    double seconds = 0.0;
    // Cells updated per second, so different grid sizes compare fairly
    double cell_updates_per_sec = 0.0;
    // Largest timestep the scheme can take without blowing up
    double dt_stable_estimate = 0.0;
    // OpenMP threads per process that did the work
    int threads = 1;
    // MPI processes that did the work
    int ranks = 1;
    // Time the first process spent swapping halo rows, part of seconds
    double halo_seconds = 0.0;
    // Time spent on the serial boundary work, part of seconds
    double boundary_seconds = 0.0;
};

// The common shape every solver follows, like a Python abstract base class.
// Each solver must provide solve() and name(), and shares the helpers below.
class Solver {
public:
    // Virtual so the right cleanup runs when deleting through a Solver pointer
    virtual ~Solver() = default;
    // The = 0 means every subclass has to write its own version
    virtual SolveResult solve(const Config& cfg) = 0;
    virtual std::string name() const = 0;

protected:
    // Helpers every solver can call, written once in solver.cpp
    void extract_result(const Grid& grid, const Config& cfg, SolveResult& out) const;
    void report_stability(const Grid& grid, const Config& cfg, double dt, SolveResult& out) const;
};

// The simple reference solver that every other solver must match.
// Recomputes everything per cell, written to be correct rather than fast.
class BaselineSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "baseline"; }
};

// The single core solver with constants worked out once up front.
// Same answers as the baseline, about twice as fast.
class SerialSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "serial"; }
};

// The single core solver that updates four cells per instruction with AVX2.
// Same time loop as SerialSolver, only the row update is different.
class SimdSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "simd"; }
};

// The solver that splits the rows between all the cores of one machine with OpenMP.
// Same time loop as SimdSolver, with one pragma above the row loop.
class OmpSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "omp"; }
};

// The solver that splits the rows between several processes with MPI, threads inside each.
// Each process owns one slab of rows and swaps one edge row with each neighbour per step.
class MpiSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "mpi"; }
};

// Creates the solver whose name matches, for example "serial".
std::unique_ptr<Solver> make_solver(const std::string& name);

#endif  // HESTON_SOLVER_H
