#include "io.h"

#include <cstdio>      // printf-family for the machine-readable result line
#include <filesystem>  // create_directories, similar to os.makedirs(..., exist_ok=True)
#include <fstream>     // std::ofstream: RAII file handle, closes at scope exit
#include <string>

void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid) {
    // Feeds weather_map.py. Only ever called outside the timed region.
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

std::string solver_label(const Config& cfg) {
    if (cfg.solver != "opt") return cfg.solver;
    // Levels 7 and 8 are the negative controls rather than ladder rungs, so
    // they get their own names. Calling them opt-L7 and opt-L8 would let a
    // CSV row be misread as the ladder carrying on upwards.
    if (cfg.opt_level == 7) return "opt-ctl-order";
    if (cfg.opt_level == 8) return "opt-ctl-branch";
    return "opt-L" + std::to_string(cfg.opt_level);
}

void print_result_csv(const SolveResult& r, const Config& cfg) {
    // This is exactly the line PLAN §3 specifies. bench_plot.py and the slurm
    // logs parse it, so it must not gain extra whitespace or columns. The opt
    // solver's label carries its ladder level so that a benchmark CSV records
    // which kernel produced each row.
    const std::string label = solver_label(cfg);
    std::printf("%.10g,%.8g,%.8g,%.8g,%d,%d,%d,%s,%.6f,%.6g\n",
                r.price, r.delta, r.gamma, r.vega, cfg.grid.num_stock_nodes,
                cfg.grid.num_variance_nodes, cfg.grid.num_timesteps,
                label.c_str(), r.seconds, r.cell_updates_per_sec);
}
