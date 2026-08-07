// Validation gate 2 (PLAN P3): put-call parity, model-independent:
//   C - P = S*exp(-q*T) - K*exp(-r*T)
// Full Heston parameters (xi != 0) — this catches whole bug classes the
// BS-collapse test cannot (cross-term, v-boundaries).

#include <cstdio>

#include "params.h"
#include "solver.h"

int main() {
    // TODO(P3): reference cfg; solve call and put on identical grid;
    // check |(C - P) - (S*e^{-qT} - K*e^{-rT})| < tolerance (measure, then
    // fix and justify); printf all four numbers; return 0 iff within.
    std::printf("test_parity: NOT IMPLEMENTED (P3)\n");
    return 1;
}
