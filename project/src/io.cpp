// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "io.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

// Saves the current sheet to dump_dir/snap_<step>.csv.
// One line per variance row and one column per stock node, for weather_map.py.
void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid) {
    // Make the folder if it is missing, like os.makedirs with exist_ok=True
    std::filesystem::create_directories(dump_dir);

    // Zero padded step number so the files sort in time order
    char file_name[32];
    std::snprintf(file_name, sizeof file_name, "snap_%06d.csv", step);
    std::string path = dump_dir + "/" + file_name;

    // The file closes itself when out goes out of scope
    std::ofstream out(path);
    if (!out) {
        std::fprintf(stderr, "dump_snapshot: cannot write %s\n", path.c_str());
        return;
    }
    out.precision(10);

    for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
        for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
            // Comma before every value except the first in a row
            if (stock_i > 0) {
                out << ',';
            }
            out << grid.current()[grid.index(stock_i, var_j)];
        }
        out << '\n';
    }
}

// Prints price, Greeks, grid size, solver name, timing, threads, processes, halo and boundary time as one CSV line.
// The benchmark scripts parse this line, so new columns only ever go on the end.
void print_result_csv(const SolveResult& result, const Config& cfg) {
    std::printf("%.10g,%.8g,%.8g,%.8g,%d,%d,%d,%s,%.6f,%.6g,%d,%d,%.6f,%.6f\n",
                result.price, result.delta, result.gamma, result.vega,
                cfg.grid.num_stock_nodes, cfg.grid.num_variance_nodes, cfg.grid.num_timesteps,
                cfg.solver.c_str(), result.seconds, result.cell_updates_per_sec,
                result.threads, result.ranks, result.halo_seconds, result.boundary_seconds);
}
