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

// Solvers that must give the same answers as the baseline.
static const std::vector<std::string> solver_names = {"serial", "simd", "omp", "mpi", "omp_scalar", "omp_persistent", "omp_runtime", "omp_edges"};

// Solvers that must match simd bit for bit.
// They split the same row updates between workers and never sum across workers.
static const std::vector<std::string> same_bits_as_simd = {"omp", "mpi", "omp_persistent", "omp_runtime", "omp_edges"};

// The worst relative difference allowed.
// Measured at most 3.3e-15 on serial with clang and g++ 8.5.
static const double relative_tolerance = 1e-11;

// The float experiment keeps about 7 digits and gets a looser limit.
// Measured at most 8.3e-6 on this grid, on the put.
static const double float_tolerance = 1e-4;

// A small grid that still reaches every boundary.
// The timestep is about 0.78 of the stability limit.
static Config test_config() {
    Config cfg;
    cfg.grid.num_stock_nodes = 128;
    cfg.grid.num_variance_nodes = 32;
    cfg.grid.num_timesteps = 5200;
    return cfg;
}

// The difference between two numbers divided by the size of the larger one.
// Gamma near 0.0008 and vega near 2000 are judged by the same rule.
static double relative_diff(double a, double b) {
    double scale = std::fmax(std::fabs(a), std::fabs(b));
    // Both being exactly zero means they agree.
    if (scale == 0.0) {
        return 0.0;
    }
    return std::fabs(a - b) / scale;
}

// The worst relative difference across price, delta, gamma and vega.
// One number per comparison.
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
    // The mpi solver needs MPI running, even on one process.
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Config cfg = test_config();
    BaselineSolver baseline;
    bool all_ok = true;

    for (int leg = 0; leg < 2; ++leg) {
        // First leg is the call, second leg is the put.
        cfg.option.is_call = (leg == 0);
        const char* leg_name = "put";
        if (cfg.option.is_call) {
            leg_name = "call";
        }
        SolveResult reference = baseline.solve(cfg);
        SolveResult simd_reference = make_solver("simd")->solve(cfg);

        // Works like for solver_name in solver_names in Python.
        for (const std::string& solver_name : solver_names) {
            // The factory hands back the only owner of the solver.
            std::unique_ptr<Solver> solver = make_solver(solver_name);
            SolveResult candidate = solver->solve(cfg);
            // Every process runs every solver. Only the first one holds the mpi answer.
            if (rank != 0) {
                continue;
            }
            double worst_diff = max_relative_diff(reference, candidate);
            // Written as less than so a NaN fails.
            bool ok = worst_diff < relative_tolerance;
            const char* verdict = "OK";
            if (!ok) {
                verdict = "MISMATCH";
                all_ok = false;
            }
            std::printf("%-8s %-4s: worst relative diff %.3e  (tol %.0e)  %s\n",
                        solver_name.c_str(), leg_name, worst_diff, relative_tolerance, verdict);
        }

        // The float experiment is held to its own looser limit.
        SolveResult float_candidate = make_solver("omp_float")->solve(cfg);
        if (rank == 0) {
            double float_diff = max_relative_diff(reference, float_candidate);
            // Written as less than so a NaN fails.
            bool ok = float_diff < float_tolerance;
            const char* verdict = "OK";
            if (!ok) {
                verdict = "MISMATCH";
                all_ok = false;
            }
            std::printf("%-8s %-4s: worst relative diff %.3e  (tol %.0e, float)  %s\n",
                        "omp_float", leg_name, float_diff, float_tolerance, verdict);
        }

        for (const std::string& solver_name : same_bits_as_simd) {
            SolveResult candidate = make_solver(solver_name)->solve(cfg);
            if (rank != 0) {
                continue;
            }
            double simd_diff = max_relative_diff(simd_reference, candidate);
            // Any difference at all means a bug.
            const char* verdict = "OK";
            if (simd_diff != 0.0) {
                verdict = "MISMATCH";
                all_ok = false;
            }
            std::printf("%-8s %-4s: diff from simd %.3e  (must be exactly 0, %d threads, %d processes)  %s\n",
                        solver_name.c_str(), leg_name, simd_diff, candidate.threads, candidate.ranks, verdict);
        }
    }

    MPI_Finalize();

    // Only the first process judged the results.
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
