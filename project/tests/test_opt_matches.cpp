// The iron rule (spec §6): no optimised version is benchmarked until it
// matches the reference. solve_opt vs solve_baseline on the smoke grid:
// max abs difference over the final sheet-readouts (price, delta, gamma,
// vega) must be < 1e-12 — same arithmetic, same order where possible.

#include <cstdio>

#include "params.h"
#include "solver.h"

int main() {
    // TODO(P7): run both solvers on smoke.cfg (call and put);
    // compare all SolveResult fields; printf diffs; return 0 iff all < 1e-12.
    std::printf("test_opt_matches: NOT IMPLEMENTED (P7)\n");
    return 1;
}
