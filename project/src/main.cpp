// The serial Heston PDE option pricer for COSC3500 Milestone 1. PLAN.md §3
// defines the command line. This file stays thin and only parses the flags,
// picks a solver and prints the result.

#include <cstdio>
#include <exception>
#include <memory>
#include <string>

#include "io.h"
#include "params.h"
#include "solver.h"

static void usage() {
    std::fprintf(stderr,
        "usage: heston --config <path> [--solver baseline|opt] [--type call|put]\n"
        "              [--opt-level 0..6 | 7=ctl-order | 8=ctl-branch]\n"
        "              [--dump-every N] [--dump-dir DIR]\n"
        "              [--bench R] [--ns X] [--nv X] [--nt X] [--maturity Y]\n");
}

// Fetches the value that has to follow a flag such as "--ns 1024".
static std::string flag_value(int argc, char** argv, int& i) {
    if (i + 1 >= argc)
        throw std::runtime_error(std::string(argv[i]) + " needs a value");
    ++i;
    return argv[i];
}

int main(int argc, char** argv) {
    try {
        // The first pass looks only for the config file, so that the file's
        // values load first and every other flag can then override them no
        // matter what order they were typed in.
        std::string config_path;
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--config")
                config_path = flag_value(argc, argv, i);
        }
        if (config_path.empty()) {
            usage();
            return 1;
        }
        Config cfg = load_config(config_path);

        // The second pass applies the command-line overrides on top.
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--config") {
                flag_value(argc, argv, i);  // already handled in pass 1
            } else if (arg == "--solver") {
                cfg.solver = flag_value(argc, argv, i);
            } else if (arg == "--opt-level") {
                cfg.opt_level = std::stoi(flag_value(argc, argv, i));
                // Levels 0 to 6 are the ladder and 7 and 8 are the two
                // negative controls. The accepted range has only ever grown,
                // so older scripts still work.
                if (cfg.opt_level < 0 || cfg.opt_level > 8)
                    throw std::runtime_error(
                        "--opt-level must be 0..6 (ladder) or 7/8 (controls)");
            } else if (arg == "--type") {
                const std::string type = flag_value(argc, argv, i);
                if (type != "call" && type != "put")
                    throw std::runtime_error("--type must be call or put");
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
                cfg.grid.num_variance_nodes =
                    std::stoi(flag_value(argc, argv, i));
            } else if (arg == "--nt") {
                cfg.grid.num_timesteps = std::stoi(flag_value(argc, argv, i));
            } else if (arg == "--maturity") {
                // This exists for benchmarking rather than pricing. Since dt
                // is the maturity divided by nt, shrinking the maturity keeps
                // a larger grid inside its stability bound without changing
                // the memory footprint or the access pattern being measured.
                cfg.option.maturity_years =
                    std::stod(flag_value(argc, argv, i));
            } else {
                usage();
                throw std::runtime_error("unknown flag: " + arg);
            }
        }

        // The unique_ptr is the sole owner of the solver and deletes it when
        // it goes out of scope. This is like Python's garbage collector,
        // except that you know exactly when it happens.
        std::unique_ptr<Solver> solver = make_solver(cfg.solver);
        // Passing --bench R repeats the whole solve R times and prints one
        // CSV line each, so the scripts can take a median across reps.
        const int reps = cfg.bench_reps > 0 ? cfg.bench_reps : 1;
        for (int rep = 0; rep < reps; ++rep) {
            const SolveResult result = solver->solve(cfg);
            print_result_csv(result, cfg);
        }
        return 0;
    } catch (const std::exception& e) {
        // This is the only catch in the program. Parsing and the factory can
        // throw, and the numerics never do.
        std::fprintf(stderr, "heston: %s\n", e.what());
        return 1;
    }
}
