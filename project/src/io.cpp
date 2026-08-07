#include "io.h"

void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid) {
    // TODO(P2): write dump_dir/snap_<step, zero-padded>.csv — nv rows, ns
    // columns, plain doubles. Create dump_dir if missing. Keep it simple;
    // this runs outside the timed loop.
    (void)dump_dir;
    (void)step;
    (void)grid;
}

void print_result_csv(const SolveResult& r, const Config& cfg) {
    // TODO(P1): exactly the PLAN §3 line, no extra whitespace —
    // scripts/bench_plot.py and slurm logs parse this.
    (void)r;
    (void)cfg;
}
