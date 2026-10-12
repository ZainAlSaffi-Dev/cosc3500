// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "params.h"
#include "solver.h"

// Solvers that must give the same answers as the baseline
static const std::vector<std::string> kRungs = {"serial", "simd", "omp", "mpi", "omp_scalar", "omp_persistent", "omp_runtime", "omp_edges"};

// Solvers that must match simd bit for bit.
// They split the same row updates between workers and never add up across workers
static const std::vector<std::string> kSameBitsAsSimd = {"omp", "mpi", "omp_persistent", "omp_runtime", "omp_edges"};

// Worst relative difference allowed.
// Measured at most 3.3e-15 on serial with clang and g++ 8.5
static const double kRelativeTolerance = 1e-11;

// The float experiment keeps only about 7 digits, so it gets its own looser limit.
// Measured at most 8.3e-6 on this grid, on the put
static const double kFloatTolerance = 1e-4;

// A small grid that still reaches every boundary.
// The timestep is about 0.78 of the stability limit, so the run is stable.
static Config test_config() {
    Config cfg;
    cfg.grid.num_stock_nodes = 128;
    cfg.grid.num_variance_nodes = 32;
    cfg.grid.num_timesteps = 5200;
    return cfg;
}

// The difference between two numbers divided by the size of the larger one.
// This way gamma near 0.0008 and vega near 2000 are judged by the same rule.
static double relative_diff(double a, double b) {
    double scale = std::fmax(std::fabs(a), std::fabs(b));
    // Both exactly zero means they agree
    if (scale == 0.0) {
        return 0.0;
    }
    return std::fabs(a - b) / scale;
}

// The worst relative difference across price, delta, gamma and vega.
// One number per comparison makes the pass or fail check simple.
static double max_relative_diff(const SolveResult& a, const SolveResult& b) {
    double worst = relative_diff(a.price, b.price);
    worst = std::fmax(worst, relative_diff(a.delta, b.delta));
    worst = std::fmax(worst, relative_diff(a.gamma, b.gamma));
    worst = std::fmax(worst, relative_diff(a.vega, b.vega));
    return worst;
}

// Runs the baseline and every other solver on the call and the put.
// Fails if any solver's answer drifts from the baseline by more than the tolerance.
int main(int argc, char** argv) {
    // The mpi rung needs MPI running, even when this is one process
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Config cfg = test_config();
    BaselineSolver baseline;
    bool all_ok = true;

    for (int leg = 0; leg < 2; ++leg) {
        // First leg is the call, second leg is the put
        cfg.option.is_call = (leg == 0);
        const char* leg_name = "put";
        if (cfg.option.is_call) {
            leg_name = "call";
        }
        SolveResult reference = baseline.solve(cfg);
        SolveResult simd_reference = make_solver("simd")->solve(cfg);

        for (const std::string& rung : kRungs) {
            // The factory hands back the only owner of the solver
            std::unique_ptr<Solver> solver = make_solver(rung);
            SolveResult candidate = solver->solve(cfg);
            // Every process runs every solver, but only the first one holds the mpi answer
            if (rank != 0) {
                continue;
            }
            double worst_diff = max_relative_diff(reference, candidate);
            bool ok = worst_diff < kRelativeTolerance;
            if (!ok) {
                all_ok = false;
            }
            std::printf("%-8s %-4s: worst relative diff %.3e  (tol %.0e)  %s\n",
                        rung.c_str(), leg_name, worst_diff, kRelativeTolerance, ok ? "OK" : "MISMATCH");
        }

        // The float experiment, held to its own looser limit
        SolveResult float_candidate = make_solver("omp_float")->solve(cfg);
        if (rank == 0) {
            double float_diff = max_relative_diff(reference, float_candidate);
            bool float_ok = float_diff < kFloatTolerance;
            if (!float_ok) {
                all_ok = false;
            }
            std::printf("%-8s %-4s: worst relative diff %.3e  (tol %.0e, float)  %s\n",
                        "omp_float", leg_name, float_diff, kFloatTolerance, float_ok ? "OK" : "MISMATCH");
        }

        for (const std::string& rung : kSameBitsAsSimd) {
            SolveResult candidate = make_solver(rung)->solve(cfg);
            if (rank != 0) {
                continue;
            }
            double simd_diff = max_relative_diff(simd_reference, candidate);
            // Any difference at all means a bug, not rounding
            bool ok = simd_diff == 0.0;
            if (!ok) {
                all_ok = false;
            }
            std::printf("%-8s %-4s: diff from simd %.3e  (must be exactly 0, %d threads, %d processes)  %s\n",
                        rung.c_str(), leg_name, simd_diff, candidate.threads, candidate.ranks, ok ? "OK" : "MISMATCH");
        }
    }

    MPI_Finalize();

    // Only the first process judged the results
    if (rank != 0) {
        return 0;
    }
    if (all_ok) {
        std::printf("test_ladder_matches: PASS\n");
        return 0;
    } else {
        std::printf("test_ladder_matches: FAIL\n");
        return 1;
    }
}
