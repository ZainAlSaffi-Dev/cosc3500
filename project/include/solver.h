#ifndef HESTON_SOLVER_H
#define HESTON_SOLVER_H

#include <memory>
#include <string>

#include "grid.h"
#include "params.h"

// Result of one full backward solve.
// similar to python dataclass
struct SolveResult {  
    double price = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    //vega is the sensitivity of the option price to changes in the volatility of the underlying asset
    double vega = 0.0;
    //seconds is the time taken to solve the option
    double seconds = 0.0;
    //cell_updates_per_sec is the number of cells updated per second
    double cell_updates_per_sec = 0.0;
    //dt_stable_estimate is the estimated stability of the time step
    double dt_stable_estimate = 0.0;
};

// Solver is an abstract class with an abstract method solve()
class Solver {
public:
// virtual destructor is needed for polymorphism
    virtual ~Solver() = default;
    // solve is an abstract method that must be implemented by all solvers
    virtual SolveResult solve(const Config& cfg) = 0;
    virtual std::string name() const = 0;

protected:
    // Shared by all implementations
    // read price + Greeks off the finished sheet by finite differences.
    void extract_result(const Grid& g, const Config& cfg, SolveResult& out) const;
};

// Correct-first reference implementation. Weights recomputed per cell.
class BaselineSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "baseline"; }
};

// Optimised serial. Same answers as BaselineSolver.
// tests/test_opt_matches.cpp before any benchmarking.
class OptSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "opt"; }
};

// Factory
// unique_ptr = explicit sole owner; C++
// makes the owner a named variable, and scope exit deletes the solver.
// throws std::runtime_error on unknown name.
std::unique_ptr<Solver> make_solver(const std::string& name);

#endif  // HESTON_SOLVER_H
