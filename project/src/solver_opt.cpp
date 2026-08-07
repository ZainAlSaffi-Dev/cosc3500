// OptSolver — PLAN §4 menu, applied one measured step at a time:
//   1. per-row precomputed stencil weights (hoisted out of inner loop)
//   2. loop order matches layout (outer j, inner contiguous i)
//   3. boundary handling outside the branch-free interior hot loop
//   4. const double* in = g.cur().data(); double* out = g.next().data();
//      (raw pointers ONLY inside the kernel, never stored in members)
//   5. (stretch) tiling — keep only if it measures faster
// MUST match BaselineSolver (tests/test_opt_matches.cpp) before benchmarking.

#include "solver.h"

SolveResult OptSolver::solve(const Config& cfg) {
    // TODO(P7): start as a copy of the finished BaselineSolver::solve, then
    // apply §4 steps one commit each, re-running test_opt_matches after each.
    (void)cfg;
    return SolveResult{};
}
