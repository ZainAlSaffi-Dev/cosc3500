// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "grid.h"
#include "params.h"
#include "stencil.h"

// Worst difference allowed, as a fraction of the largest value.
// FMA rounds once where the scalar row rounds twice.
static const double tolerance = 1e-13;

// The same check for the float row, scaled to float precision.
// One float rounding step is about 6e-8 of the value.
static const double float_tolerance = 1e-5;

// Runs the scalar row and the AVX row on the same random grid values and compares them.
// Returns the worst difference divided by the largest grid value.
static double compare_rows(int num_stock_nodes, std::mt19937& random_engine) {
    Config cfg;
    cfg.grid.num_stock_nodes = num_stock_nodes;
    cfg.grid.num_variance_nodes = 32;
    Grid grid(cfg.grid, cfg.option);
    StencilConstants constants(grid, cfg, 1e-5);

    // A call payoff with random noise on top, so every derivative is non zero.
    int num_cells = num_stock_nodes * cfg.grid.num_variance_nodes;
    std::vector<double> current(num_cells);
    // noise(random_engine) gives a random number between -1 and 1.
    std::uniform_real_distribution<double> noise(-1.0, 1.0);
    for (int cell = 0; cell < num_cells; ++cell) {
        int stock_i = cell % num_stock_nodes;
        double payoff = std::fmax(constants.stock_price[stock_i] - cfg.option.strike, 0.0);
        current[cell] = payoff + noise(random_engine);
    }

    // Both outputs start with the same marker to catch stray writes.
    std::vector<double> next_scalar(num_cells, -7.0);
    std::vector<double> next_avx(num_cells, -7.0);
    for (int var_j = 1; var_j < cfg.grid.num_variance_nodes - 1; ++var_j) {
        update_row_scalar(constants, current, next_scalar, var_j);
        update_row_avx(constants, current, next_avx, var_j);
    }

    // Compare every cell, edges included.
    double largest = 0.0;
    double worst = 0.0;
    for (int cell = 0; cell < num_cells; ++cell) {
        largest = std::fmax(largest, std::fabs(next_scalar[cell]));
        worst = std::fmax(worst, std::fabs(next_scalar[cell] - next_avx[cell]));
    }
    return worst / largest;
}

// Runs the float scalar row and the float AVX row on the same random grid values and compares them.
// Returns the worst difference divided by the largest grid value.
static double compare_float_rows(int num_stock_nodes, std::mt19937& random_engine) {
    Config cfg;
    cfg.grid.num_stock_nodes = num_stock_nodes;
    cfg.grid.num_variance_nodes = 32;
    Grid grid(cfg.grid, cfg.option);
    StencilConstants constants(grid, cfg, 1e-5);
    FloatStencilConstants float_constants(constants);

    // A call payoff with random noise on top, as in compare_rows.
    int num_cells = num_stock_nodes * cfg.grid.num_variance_nodes;
    std::vector<float> current(num_cells);
    std::uniform_real_distribution<float> noise(-1.0f, 1.0f);
    for (int cell = 0; cell < num_cells; ++cell) {
        int stock_i = cell % num_stock_nodes;
        float payoff = std::fmax(float_constants.stock_price[stock_i] - (float)cfg.option.strike, 0.0f);
        current[cell] = payoff + noise(random_engine);
    }

    // Both outputs start with the same marker to catch stray writes.
    std::vector<float> next_scalar(num_cells, -7.0f);
    std::vector<float> next_avx(num_cells, -7.0f);
    for (int var_j = 1; var_j < cfg.grid.num_variance_nodes - 1; ++var_j) {
        update_cells_scalar_float(float_constants, current, next_scalar, var_j, 1, num_stock_nodes - 1);
        update_row_avx_float(float_constants, current, next_avx, var_j);
    }

    // Compare every cell, edges included.
    double largest = 0.0;
    double worst = 0.0;
    for (int cell = 0; cell < num_cells; ++cell) {
        largest = std::fmax(largest, std::fabs(next_scalar[cell]));
        worst = std::fmax(worst, std::fabs(next_scalar[cell] - next_avx[cell]));
    }
    return worst / largest;
}

// Checks the AVX2 rows against the scalar rows for several row lengths.
// The lengths cover every remainder left after the groups of four doubles or eight floats.
int main() {
    // Fail if the AVX2 row was not built.
    if (!avx_row_is_vectorised()) {
        std::printf("test_avx_row: FAIL, built without AVX2 and FMA\n");
        return 1;
    }

    // Random number generator with a fixed seed, like random.seed(3500) in Python.
    // Every run uses the same random grid values.
    std::mt19937 random_engine(3500);
    std::vector<int> row_lengths = {130, 131, 132, 133, 2048};
    bool all_ok = true;

    // Works like for num_stock_nodes in row_lengths in Python.
    for (int num_stock_nodes : row_lengths) {
        double worst = compare_rows(num_stock_nodes, random_engine);
        // Written as less than so a NaN fails.
        bool ok = worst < tolerance;
        const char* verdict = "OK";
        if (!ok) {
            verdict = "MISMATCH";
            all_ok = false;
        }
        // Interior cells left over after the groups of four.
        int remainder = (num_stock_nodes - 2) % 4;
        std::printf("ns %5d (remainder %d): worst diff %.3e of the largest value  (tol %.0e)  %s\n",
                    num_stock_nodes, remainder, worst, tolerance, verdict);
    }

    // The float row does eight cells per pass, so test remainders up to seven.
    std::vector<int> float_row_lengths = {130, 131, 132, 133, 134, 135, 136, 137, 2048};
    for (int num_stock_nodes : float_row_lengths) {
        double worst = compare_float_rows(num_stock_nodes, random_engine);
        // Written as less than so a NaN fails.
        bool ok = worst < float_tolerance;
        const char* verdict = "OK";
        if (!ok) {
            verdict = "MISMATCH";
            all_ok = false;
        }
        int remainder = (num_stock_nodes - 2) % 8;
        std::printf("float ns %5d (remainder %d): worst diff %.3e of the largest value  (tol %.0e)  %s\n",
                    num_stock_nodes, remainder, worst, float_tolerance, verdict);
    }

    if (all_ok) {
        std::printf("test_avx_row: PASS\n");
        return 0;
    } else {
        std::printf("test_avx_row: FAIL\n");
        return 1;
    }
}
