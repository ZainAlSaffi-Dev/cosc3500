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
}

// Returns the value that follows a flag, such as 1024 in --ns 1024.
// Moves i forward past the value so the caller skips it.
static std::string flag_value(int argc, char** argv, int& i) {
    if (i + 1 >= argc) {
        throw std::runtime_error(std::string(argv[i]) + " needs a value");
    }
    i = i + 1;
    return argv[i];
}

// Reads the settings, runs the chosen solver and prints one CSV line per solve.
// Only the first process prints, and any error is caught and printed instead of crashing.
static int run(int argc, char** argv, int rank) {
    try {
        // First pass finds only the config file.
        // That way the other flags can override it in any order
        std::string config_path;
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--config") {
                config_path = flag_value(argc, argv, i);
            }
        }
        if (config_path.empty()) {
            usage();
            return 1;
        }
        Config cfg = load_config(config_path);

        // Second pass applies every other flag on top of the config file
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--config") {
                // Already read in the first pass, just skip its value
                flag_value(argc, argv, i);
            } else if (arg == "--solver") {
                cfg.solver = flag_value(argc, argv, i);
            } else if (arg == "--type") {
                std::string type = flag_value(argc, argv, i);
                if (type != "call" && type != "put") {
                    throw std::runtime_error("--type must be call or put");
                }
                cfg.option.is_call = (type == "call");
            } else if (arg == "--dump-every") {
                cfg.dump_every = std::stoi(flag_value(argc, argv, i));
            } else if (arg == "--dump-dir") {
                cfg.dump_dir = flag_value(argc, argv, i);
            } else if (arg == "--bench") {
                cfg.bench_reps = std::stoi(flag_value(argc, argv, i));
            } else if (arg == "--ns") {
                cfg.grid.num_stock_nodes = std::stoi(flag_value(argc, argv, i));
            } else if (arg == "--nv") {
                cfg.grid.num_variance_nodes = std::stoi(flag_value(argc, argv, i));
            } else if (arg == "--nt") {
                cfg.grid.num_timesteps = std::stoi(flag_value(argc, argv, i));
            } else if (arg == "--vmax") {
                // The weak scaling sweep grows the grid in variance at a fixed spacing
                cfg.grid.variance_max = std::stod(flag_value(argc, argv, i));
            } else if (arg == "--maturity") {
                // A shorter maturity keeps a big grid stable for benchmarks
                cfg.option.maturity_years = std::stod(flag_value(argc, argv, i));
            } else {
                // Anything else is a typo, so stop
                usage();
                throw std::runtime_error("unknown flag: " + arg);
            }
        }

        // The unique_ptr owns the solver and deletes it at the end of run
        std::unique_ptr<Solver> solver = make_solver(cfg.solver);

        // --bench R repeats the solve so the scripts can take a median
        int reps = 1;
        if (cfg.bench_reps > 0) {
            reps = cfg.bench_reps;
        }
        for (int rep = 0; rep < reps; ++rep) {
            SolveResult result = solver->solve(cfg);
            // With several processes only the first one prints the result
            if (rank == 0) {
                print_result_csv(result, cfg);
            }
        }
        return 0;
    } catch (const std::exception& e) {
        // Reading the config, picking the solver or a bad process count can throw
        std::fprintf(stderr, "heston: %s\n", e.what());
        return 1;
    }
}

// Starts MPI, runs the program, and shuts MPI down again.
// MPI_Finalize must run even when the program fails, so the work is in run().
int main(int argc, char** argv) {
    // Every rung runs inside MPI, the non MPI ones simply as one process
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int exit_code = run(argc, argv, rank);
    MPI_Finalize();
    return exit_code;
}
