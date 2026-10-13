// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "io.h"

#include <cstdio>
#include <filesystem>
#include <string>

// Saves the current buffer to dump_dir/snap_<step>.csv.
// One line per variance row and one column per stock node.
void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid) {
    // Make the folder if it is missing, like os.makedirs in Python.
    std::filesystem::create_directories(dump_dir);

    // Zero padded step number, so the files sort in time order.
    char file_name[32];
    std::snprintf(file_name, sizeof file_name, "snap_%06d.csv", step);
    std::string path = dump_dir + "/" + file_name;

    // c_str() gives the C string that fopen expects.
    FILE* out = std::fopen(path.c_str(), "w");
    if (out == NULL) {
        std::fprintf(stderr, "dump_snapshot: cannot write %s\n", path.c_str());
        return;
    }

    for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
        for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
            // Comma before every value except the first in a row.
            if (stock_i > 0) {
                std::fprintf(out, ",");
            }
            // %.10g writes 10 significant digits.
            std::fprintf(out, "%.10g", grid.current()[grid.index(stock_i, var_j)]);
        }
        std::fprintf(out, "\n");
    }
    std::fclose(out);
}

// Prints the price, Greeks, grid size, solver, timing, threads, processes, row exchange and boundary time as one CSV line.
// The scripts read the columns by position, so new columns go on the end.
void print_result_csv(const SolveResult& result, const Config& cfg) {
    std::printf("%.10g,%.8g,%.8g,%.8g,%d,%d,%d,%s,%.6f,%.6g,%d,%d,%.6f,%.6f\n",
                result.price, result.delta, result.gamma, result.vega,
                cfg.grid.num_stock_nodes, cfg.grid.num_variance_nodes, cfg.grid.num_timesteps,
                cfg.solver.c_str(), result.seconds, result.cell_updates_per_sec,
                result.threads, result.ranks, result.exchange_seconds, result.boundary_seconds);
}
