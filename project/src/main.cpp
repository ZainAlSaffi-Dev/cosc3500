// heston — serial Heston PDE option pricer (COSC3500 Milestone 1).
// CLI contract in PLAN.md §3. Keep this file thin: parse flags, dispatch, print.

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
        "              [--opt-level 0..6] [--dump-every N] [--dump-dir DIR]\n"
        "              [--bench R] [--ns X] [--nv X] [--nt X]\n");
}

// Fetch the value that must follow a flag like "--ns 1024".
static std::string flag_value(int argc, char** argv, int& i) {
    if (i + 1 >= argc)
        throw std::runtime_error(std::string(argv[i]) + " needs a value");
    ++i;
    return argv[i];
}

int main(int argc, char** argv) {
    try {
        // Pass 1: find the config file only, so file values load first and
        // every other flag can override them afterwards, regardless of order.
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

        // Pass 2: CLI overrides on top of the file values.
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--config") {
                flag_value(argc, argv, i);  // already handled in pass 1
            } else if (arg == "--solver") {
                cfg.solver = flag_value(argc, argv, i);
            } else if (arg == "--opt-level") {
                cfg.opt_level = std::stoi(flag_value(argc, argv, i));
                if (cfg.opt_level < 0 || cfg.opt_level > 6)
                    throw std::runtime_error("--opt-level must be 0..6");
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
            } else {
                usage();
                throw std::runtime_error("unknown flag: " + arg);
            }
        }

        // unique_ptr: sole owner of the polymorphic solver — deleted
        // automatically at scope exit, like Python GC but deterministic.
        std::unique_ptr<Solver> solver = make_solver(cfg.solver);
        // --bench R repeats the whole solve R times, one CSV line each,
        // so scripts get per-rep timings (median/min/max downstream).
        const int reps = cfg.bench_reps > 0 ? cfg.bench_reps : 1;
        for (int rep = 0; rep < reps; ++rep) {
            const SolveResult result = solver->solve(cfg);
            print_result_csv(result, cfg);
        }
        return 0;
    } catch (const std::exception& e) {
        // The single catch (PLAN §1b): parsing and the factory throw; nothing
        // in the numerics does.
        std::fprintf(stderr, "heston: %s\n", e.what());
        return 1;
    }
}
