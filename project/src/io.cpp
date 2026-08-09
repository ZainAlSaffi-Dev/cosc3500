#include "io.h"

#include <cstdio>      // printf-family for the machine-readable result line
#include <filesystem>  // create_directories — similar to os.makedirs(..., exist_ok=True)
#include <fstream>     // std::ofstream: RAII file handle, closes at scope exit
#include <string>

void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid) {
    // Feeds weather_map.py; only ever called outside the timed region.
    std::filesystem::create_directories(dump_dir);
    // Zero-padded step number so filenames sort in time order.
    char name[32];
    std::snprintf(name, sizeof name, "snap_%06d.csv", step);
    const std::string path = dump_dir + "/" + name;
    std::ofstream out(path);
    if (!out) {
        std::fprintf(stderr, "dump_snapshot: cannot write %s\n", path.c_str());
        return;
    }
    out.precision(10);
    // Row = variance level, column = stock node (matches weather_map.py).
    for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
        for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
            if (stock_i > 0) out << ',';
            out << grid.current()[grid.index(stock_i, var_j)];
        }
        out << '\n';
    }
}

void print_result_csv(const SolveResult& r, const Config& cfg) {
    // Exactly the PLAN §3 line — bench_plot.py and slurm logs parse this,
    // so no extra whitespace or columns.
    std::printf("%.10g,%.8g,%.8g,%.8g,%d,%d,%d,%s,%.6f,%.6g\n",
                r.price, r.delta, r.gamma, r.vega, cfg.grid.num_stock_nodes,
                cfg.grid.num_variance_nodes, cfg.grid.num_timesteps,
                cfg.solver.c_str(), r.seconds, r.cell_updates_per_sec);
}
