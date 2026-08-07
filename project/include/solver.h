#ifndef HESTON_SOLVER_H
#define HESTON_SOLVER_H

#include <memory>
#include <string>

#include "grid.h"
#include "params.h"

// Result of one full backward solve.
struct SolveResult {  // plain data bag — Python @dataclass. No behaviour.
    double price = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double vega = 0.0;   // dV/dv at (spot, v0); chain rule to per-vol in write-up
    double seconds = 0.0;              // time loop only, no I/O
    double cell_updates_per_sec = 0.0; // ns*nv*nt / seconds
    double dt_stable_estimate = 0.0;   // explicit-scheme stability bound (PLAN §1)
};

// Python analogy: class Solver(ABC) with @abstractmethod solve().
// `virtual ... = 0` = abstract method; `override` in subclasses catches
// method-name typos at compile time (Python would silently add a new method).
// Milestone 2 adds SimdSolver / OpenMpSolver behind this same interface.
class Solver {
public:
    virtual ~Solver() = default;  // rule: interfaces need a virtual destructor
    virtual SolveResult solve(const Config& cfg) = 0;
    virtual std::string name() const = 0;

protected:
    // Shared by all implementations (Python: a normal method on the ABC):
    // read price + Greeks off the finished sheet by finite differences.
    void extract_result(const Grid& g, const Config& cfg, SolveResult& out) const;
};

// Correct-first reference implementation. Weights recomputed per cell.
class BaselineSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "baseline"; }
};

// Optimised serial (PLAN §4). Same answers as BaselineSolver — enforced by
// tests/test_opt_matches.cpp before any benchmarking.
class OptSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "opt"; }
};

// Factory (Python: a function returning whichever subclass fits the name).
// unique_ptr = explicit sole owner — Python's GC does this invisibly; C++
// makes the owner a named variable, and scope exit deletes the solver.
// This is the only place polymorphism forces a pointer.
// Throws std::runtime_error on unknown name.
std::unique_ptr<Solver> make_solver(const std::string& name);

#endif  // HESTON_SOLVER_H
