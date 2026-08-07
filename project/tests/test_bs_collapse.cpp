// Validation gate 1 (PLAN P3): with xi = 0 and v0 = theta, Heston degenerates
// to Black-Scholes with sigma = sqrt(theta). Solver price must converge to the
// closed form. Tolerance: rel err < 1e-3 on the reference-sized grid —
// measure the actual number and record it in the video script.

#include <cstdio>

#include "black_scholes.h"
#include "params.h"
#include "solver.h"

int main() {
    // TODO(P3): cfg = reference values; heston.xi = 0; heston.v0 = heston.theta;
    // run solve_baseline for call AND put; compare each against bs_price;
    // printf both prices + rel errs; return 0 iff both within tolerance.
    std::printf("test_bs_collapse: NOT IMPLEMENTED (P3)\n");
    return 1;
}
