// PLAN §6 says that nothing optimised gets benchmarked until it matches the
// reference, and this test is what enforces that. Every ladder level that has
// been written, along with both negative controls, has to reproduce
// BaselineSolver's price, delta, gamma and vega on both the call and the put.
// A level whose kernel does not exist yet reports SKIP, so this test covers
// more of the ladder as kernels land without needing to be edited.
//
// Two different claims are being checked, so there are two tolerances.
// Levels 0 and 1 have to be bit-identical, because those rungs only give
// subexpressions a name, so the check is a plain comparison against zero.
// Levels 2 upwards introduce reciprocals and so differ in the last bits, and
// they are checked relative to the size of each quantity, since one absolute
// tolerance cannot cover a price near 196 and a gamma near 0.00079 at once.

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "params.h"
#include "solver.h"

// A small but complete grid that still exercises the interior, all four
// boundaries and the readout. It does not need to be node-aligned, because
// this test compares one solver against another at the same cells, so any
// interpolation weight applies equally to both sides. The timestep count
// keeps dt at roughly 0.78 of the stability bound.
static Config test_config() {
    Config cfg;  // params.h defaults are the reference market/model
    cfg.grid.num_stock_nodes = 128;
    cfg.grid.num_variance_nodes = 32;
    cfg.grid.num_timesteps = 5200;
    return cfg;
}

// Highest ladder level plus the two controls.
static const int kMaxLevel = 8;
// Levels at or below this must be bit-identical to the baseline.
static const int kBitIdenticalThrough = 1;
// The worst relative difference allowed from level 2 upwards, measured before
// it was chosen. On this grid the worst case over 5200 steps is 3.3e-15 on
// the call and exactly zero on the put, so 1e-11 leaves room for a different
// compiler while staying far tighter than anything that could hide a bug.
//
// Every level from 2 upwards reports the same difference, which this number
// also checks, because those kernels change addressing and loop structure but
// never the arithmetic.
static const double kRelativeTolerance = 1e-11;

// Largest of the four readout differences, one number per comparison.
static double max_abs_diff(const SolveResult& a, const SolveResult& b) {
    double worst = std::fabs(a.price - b.price);
    worst = std::fmax(worst, std::fabs(a.delta - b.delta));
    worst = std::fmax(worst, std::fabs(a.gamma - b.gamma));
    worst = std::fmax(worst, std::fabs(a.vega - b.vega));
    return worst;
}

// |a - b| divided by the size of the quantity itself, so gamma at 7.9e-4 and
// vega at 2190 get held to the same standard.
static double relative_diff(double a, double b) {
    const double scale = std::fmax(std::fabs(a), std::fabs(b));
    if (scale == 0.0) return 0.0;  // both exactly zero: no disagreement
    return std::fabs(a - b) / scale;
}

static double max_relative_diff(const SolveResult& a, const SolveResult& b) {
    double worst = relative_diff(a.price, b.price);
    worst = std::fmax(worst, relative_diff(a.delta, b.delta));
    worst = std::fmax(worst, relative_diff(a.gamma, b.gamma));
    worst = std::fmax(worst, relative_diff(a.vega, b.vega));
    return worst;
}

// "L0".."L6" for the ladder, plus names for the two controls. These match
// the CSV labels in io.cpp so a failure here points at the right bench row.
static std::string level_name(int level) {
    if (level == 7) return "ctl-order";
    if (level == 8) return "ctl-branch";
    return "L" + std::to_string(level);
}

int main() {
    Config cfg = test_config();

    BaselineSolver baseline;
    OptSolver opt;
    bool any_mismatch = false;
    int levels_written = 0;

    for (int leg = 0; leg < 2; ++leg) {
        cfg.option.is_call = (leg == 0);
        const char* leg_name = cfg.option.is_call ? "call" : "put";
        const SolveResult reference = baseline.solve(cfg);

        for (int level = 0; level <= kMaxLevel; ++level) {
            cfg.opt_level = level;
            SolveResult candidate;
            try {
                candidate = opt.solve(cfg);
            } catch (const std::runtime_error&) {
                // No kernel for this rung yet. Skip it rather than let it
                // count as a pass.
                std::printf("opt %-10s %-4s: SKIP (kernel not written)\n",
                            level_name(level).c_str(), leg_name);
                continue;
            }
            if (leg == 0) ++levels_written;

            bool ok = false;
            if (level <= kBitIdenticalThrough) {
                const double diff = max_abs_diff(reference, candidate);
                ok = (diff == 0.0);
                std::printf("opt %-10s %-4s: max |diff| %.3e   "
                            "(bit-identical required)  %s\n",
                            level_name(level).c_str(), leg_name, diff,
                            ok ? "OK" : "MISMATCH");
            } else {
                const double rel = max_relative_diff(reference, candidate);
                ok = (rel < kRelativeTolerance);
                std::printf("opt %-10s %-4s: max rel  %.3e   (tol %.0e)"
                            "            %s\n",
                            level_name(level).c_str(), leg_name, rel,
                            kRelativeTolerance, ok ? "OK" : "MISMATCH");
            }
            if (!ok) any_mismatch = true;
        }
    }

    // Passing requires that the anchor level exists and that every level
    // which does exist matches. If everything skipped then the ladder has not
    // been started, and that counts as a failure, because this test should
    // never report success for having done nothing.
    const bool pass = !any_mismatch && levels_written > 0;
    std::printf("test_opt_matches: %s (%d of %d level(s) written)\n",
                pass ? "PASS" : "FAIL", levels_written, kMaxLevel + 1);
    return pass ? 0 : 1;
}
