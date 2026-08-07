// heston — serial Heston PDE option pricer (COSC3500 Milestone 1).
// CLI contract in PLAN.md §3. Keep this file thin: parse flags, dispatch, print.

#include <cstdio>
#include <exception>
#include <memory>

#include "io.h"
#include "params.h"
#include "solver.h"

static void usage() {
    std::fprintf(stderr,
        "usage: heston --config <path> [--solver baseline|opt] [--type call|put]\n"
        "              [--dump-every N] [--dump-dir DIR] [--bench R]\n"
        "              [--ns X] [--nv X] [--nt X]\n");
}

int main(int argc, char** argv) {
    try {
        // TODO(P1): walk argv: load_config on --config, then apply CLI
        // overrides onto the Config fields.
        // TODO(P1): std::unique_ptr<Solver> solver = make_solver(cfg.solver);
        // --bench R loops solver->solve(cfg) R times; print_result_csv each.
        (void)argc;
        (void)argv;
        usage();
        return 1;  // P1 exit: real dispatch replaces this
    } catch (const std::exception& e) {
        // Only load_config and make_solver throw (PLAN §1b). One catch, here.
        std::fprintf(stderr, "heston: %s\n", e.what());
        return 1;
    }
}
