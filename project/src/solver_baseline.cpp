// BaselineSolver — correct-first, deliberately straightforward:
// stencil weights recomputed per cell, no hoisting, no cleverness.
// This is the "before" of the optimisation story and the reference
// every other solver must match. Resist optimising this file.

#include <stdexcept>

#include "io.h"
#include "solver.h"

SolveResult BaselineSolver::solve(const Config& cfg) {
    // TODO(P2), in this order:
    //  1. Grid g(cfg.grid, cfg.option); g.init_payoff(cfg.option);
    //  2. dt = T/nt; compute dt_stable_estimate (PLAN §1), store in result.
    //  3. Time loop n = 1..nt (stepping backwards from expiry):
    //     a. interior cells (j=1..nv-2 outer, i=1..ns-2 inner): 9-point
    //        blend of g.cur() into g.next() — weights from discretised PDE.
    //     b. boundaries: S=0, S=Smax (Dirichlet), v=0 (one-sided reduced
    //        equation), v=vmax (Neumann copy-form). All per PLAN §1 table.
    //     c. g.swap_buffers();
    //     d. if (cfg.dump_every > 0 && n % cfg.dump_every == 0)
    //          dump_snapshot(cfg.dump_dir, n, g);   // outside timed region
    //  4. Timing: steady_clock around the loop, minus snapshot time
    //     (simplest: time each step, sum; or dump in a separate pass).
    //  5. extract_result(g, cfg, result);
    (void)cfg;
    return SolveResult{};
}

void Solver::extract_result(const Grid& g, const Config& cfg,
                            SolveResult& out) const {
    // TODO(P2): i0 = g.nearest_i(spot), j0 = g.nearest_j(v0);
    // price = g.cur()[g.idx(i0, j0)];
    // delta = (V[i0+1] - V[i0-1]) / (2*ds); gamma second difference;
    // vega  = (V[j0+1] - V[j0-1]) / (2*dv)  (per-variance, note chain rule).
    (void)g;
    (void)cfg;
    (void)out;
}

std::unique_ptr<Solver> make_solver(const std::string& name) {
    // std::make_unique wraps `new` so ownership is never loose — the
    // returned unique_ptr deletes the solver when it goes out of scope.
    if (name == "baseline") return std::make_unique<BaselineSolver>();
    if (name == "opt") return std::make_unique<OptSolver>();
    throw std::runtime_error("unknown solver: " + name);
}
