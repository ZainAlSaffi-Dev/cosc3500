// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <mpi.h>

#include <cstdio>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>

#include "io.h"
#include "params.h"
#include "solver.h"

// Prints how to use the program.
// Shown when the config file is missing or a flag is not recognised.
static void usage() {
    std::fprintf(stderr, "usage: heston --config <path> [--solver baseline|serial|simd|omp|mpi]\n");
    std::fprintf(stderr, "              [--type call|put] [--dump-every N] [--dump-dir DIR]\n");
    std::fprintf(stderr, "              [--bench R] [--ns X] [--nv X] [--nt X] [--maturity Y] [--vmax V]\n");
    std::fprintf(stderr, "              [--align-shift K]  (experiment E2, the omp_align solvers only)\n");
    std::fprintf(stderr, "              [--digits D]  (significant digits of the price and Greeks)\n");
}

// Reads the settings, runs the chosen solver and prints one CSV line for each solve.
// Only the first process prints. Errors are caught and printed here.
static int run(int argc, char** argv, int rank) {
    // Works like try and except in Python. A throw anywhere inside jumps to the catch below.
    try {
        // The config file comes first, as in heston --config file.cfg.
        if (argc < 3 || std::string(argv[1]) != "--config") {
            usage();
            return 1;
        }
        Config cfg = load_config(argv[2]);

        // Every other flag is a pair, the flag then its value.
        // Each flag overwrites the value from the config file.
        for (int i = 3; i < argc; i = i + 2) {
            std::string flag = argv[i];
            // A flag typed last with no value after it.
            if (i + 1 >= argc) {
                throw std::runtime_error(flag + " needs a value");
            }
            std::string value = argv[i + 1];

            if (flag == "--solver") {
                cfg.solver = value;
            } else if (flag == "--type") {
                if (value != "call" && value != "put") {
                    throw std::runtime_error("--type must be call or put");
                }
                cfg.option.is_call = (value == "call");
            } else if (flag == "--dump-every") {
                cfg.dump_every = std::stoi(value);
            } else if (flag == "--dump-dir") {
                cfg.dump_dir = value;
            } else if (flag == "--bench") {
                cfg.bench_reps = std::stoi(value);
            } else if (flag == "--ns") {
                cfg.grid.num_stock_nodes = std::stoi(value);
            } else if (flag == "--nv") {
                cfg.grid.num_variance_nodes = std::stoi(value);
            } else if (flag == "--nt") {
                cfg.grid.num_timesteps = std::stoi(value);
            } else if (flag == "--vmax") {
                // Used by the weak scaling sweep to grow the variance axis.
                cfg.grid.variance_max = std::stod(value);
            } else if (flag == "--maturity") {
                // A shorter maturity keeps a big grid stable in the benchmarks.
                cfg.option.maturity_years = std::stod(value);
            } else if (flag == "--align-shift") {
                // Where column 1 sits in its cache line, for the alignment experiment.
                cfg.align_shift = std::stoi(value);
            } else if (flag == "--digits") {
                // Prints the price and the Greeks to this many significant digits.
                cfg.print_digits = std::stoi(value);
            } else {
                // Anything else is a typo.
                usage();
                throw std::runtime_error("unknown flag: " + flag);
            }
        }

        // The unique_ptr owns the solver and deletes it at the end of run.
        std::unique_ptr<Solver> solver = make_solver(cfg.solver);

        // --bench R repeats the solve R times.
        int reps = 1;
        if (cfg.bench_reps > 0) {
            reps = cfg.bench_reps;
        }
        for (int rep = 0; rep < reps; ++rep) {
            SolveResult result = solver->solve(cfg);
            // With several processes only the first one prints.
            if (rank == 0) {
                print_result_csv(result, cfg);
            }
        }
        return 0;
    } catch (const std::exception& e) {
        // Errors from the config, the solver choice or the process count land here.
        std::fprintf(stderr, "heston: %s\n", e.what());
        return 1;
    }
}

// Starts MPI, runs the program, and shuts MPI down again.
// MPI_Finalize must run even when the program fails, so the work is in run().
int main(int argc, char** argv) {
    // Every solver runs inside MPI. The non MPI solvers run as one process.
    // FUNNELED promises MPI that only the main thread calls it, even while OpenMP threads run.
    // The mpi_overlap solver needs that, because its master thread posts messages inside a parallel region.
    int thread_level = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int exit_code = run(argc, argv, rank);
    MPI_Finalize();
    return exit_code;
}
