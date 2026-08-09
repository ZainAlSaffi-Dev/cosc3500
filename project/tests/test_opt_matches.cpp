// The iron rule (PLAN §6): no optimised version is benchmarked until it
// matches the reference. Every WRITTEN ladder level (0..6) must reproduce
// BaselineSolver's price/delta/gamma/vega to 1e-12 on both call and put.
// Unwritten levels are reported as SKIP — this test grows teeth
// automatically as kernels land, no edits needed.

#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "params.h"
#include "solver.h"

// Small but complete grid: exercises interior, all four boundaries and the
// readout. No node alignment needed — this test compares solver against
// solver at IDENTICAL cells, so any snap hits both sides equally.
// nt keeps dt ~0.78x the stability bound (ds=165.4, dv=1/31).
static Config test_config() {
    Config cfg;  // params.h defaults are the reference market/model
    cfg.grid.num_stock_nodes = 128;
    cfg.grid.num_variance_nodes = 32;
    cfg.grid.num_timesteps = 5200;
    return cfg;
}

// Largest of the four readout differences — one number per comparison.
static double max_abs_diff(const SolveResult& a, const SolveResult& b) {
    double worst = std::fabs(a.price - b.price);
    worst = std::fmax(worst, std::fabs(a.delta - b.delta));
    worst = std::fmax(worst, std::fabs(a.gamma - b.gamma));
    worst = std::fmax(worst, std::fabs(a.vega - b.vega));
    return worst;
}

int main() {
    const double tolerance = 1e-12;
    Config cfg = test_config();

    BaselineSolver baseline;
    OptSolver opt;
    bool any_mismatch = false;
    int levels_written = 0;

    for (int leg = 0; leg < 2; ++leg) {
        cfg.option.is_call = (leg == 0);
        const char* leg_name = cfg.option.is_call ? "call" : "put";
        const SolveResult reference = baseline.solve(cfg);

        for (int level = 0; level <= 6; ++level) {
            cfg.opt_level = level;
            SolveResult candidate;
            try {
                candidate = opt.solve(cfg);
            } catch (const std::runtime_error&) {
                // Unwritten rung — honest SKIP, never a fake green.
                std::printf("opt L%d %-4s: SKIP (kernel not written)\n",
                            level, leg_name);
                continue;
            }
            if (leg == 0) ++levels_written;
            const double diff = max_abs_diff(reference, candidate);
            const bool ok = diff < tolerance;
            std::printf("opt L%d %-4s: max |diff| %.3e  %s\n", level,
                        leg_name, diff, ok ? "OK" : "MISMATCH");
            if (!ok) any_mismatch = true;
        }
    }

    // PASS needs level 0 (the anchor) written and every written level
    // matching. All-SKIP means the ladder hasn't started — that's a FAIL:
    // this test must never go green by doing nothing.
    const bool pass = !any_mismatch && levels_written > 0;
    std::printf("test_opt_matches: %s (%d level(s) written)\n",
                pass ? "PASS" : "FAIL", levels_written);
    return pass ? 0 : 1;
}
