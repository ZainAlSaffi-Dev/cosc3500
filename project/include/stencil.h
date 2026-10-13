// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_STENCIL_H
#define HESTON_STENCIL_H

#include <vector>

#include "grid.h"
#include "params.h"

// Every number the update needs that stays the same for the whole solve.
// Filled once before the time loop.
struct StencilConstants {
    // The grid size.
    int num_stock_nodes = 0;
    int num_variance_nodes = 0;

    // The size of one timestep in years.
    double dt = 0.0;
    // The risk free rate.
    double rate = 0.0;
    // The dividend yield.
    double div_yield = 0.0;
    // Rate minus dividend yield. Multiplies S times dV/dS.
    double stock_drift = 0.0;

    // One over each spacing, so the loop multiplies.
    // ds is the stock spacing and dv is the variance spacing.
    // 1 / (2 ds).
    double inv_two_ds = 0.0;
    // 1 / ds squared.
    double inv_ds_sq = 0.0;
    // 1 / (2 dv).
    double inv_two_dv = 0.0;
    // 1 / dv squared.
    double inv_dv_sq = 0.0;
    // 1 / (4 ds dv).
    double inv_four_ds_dv = 0.0;
    // 1 / dv.
    double inv_dv = 0.0;

    // Contract details for the stock edges.
    bool is_call = true;
    double strike = 0.0;
    double max_stock_price = 0.0;

    // Stock price and its square, one entry per stock column.
    std::vector<double> stock_price;
    std::vector<double> stock_price_sq;

    // PDE coefficients, one entry per variance row.
    // 0.5 v. Multiplies S squared times d2V/dS2.
    std::vector<double> stock_diffusion;
    // rho xi v. Multiplies S times d2V/dSdv.
    std::vector<double> cross_diffusion;
    // 0.5 xi squared v. Multiplies d2V/dv2.
    std::vector<double> variance_diffusion;
    // kappa (theta - v). Multiplies dV/dv.
    std::vector<double> variance_drift;

    // Fills every field above, in stencil.cpp.
    StencilConstants(const Grid& grid, const Config& cfg, double timestep);
};

// Updates the cells first_stock to end_stock - 1 of one variance row, one cell at a time.
// The & passes each vector without a copy. const means read only.
void update_cells_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j, int first_stock, int end_stock);

// The same with each buffer passed as a pointer to its cell 0.
// The vector versions call these, so the maths lives in one place.
// The alignment experiment needs them, because its buffers start part way into a vector.
void update_cells_scalar(const StencilConstants& constants, const double* current, double* next, int var_j, int first_stock, int end_stock);

// Updates one interior variance row of the next buffer, one cell at a time.
void update_row_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j);

// Updates one interior variance row four cells at a time with AVX2.
// Defined in stencil_avx.cpp, the only file built with the AVX2 and FMA flags.
void update_row_avx(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j);

// The same row with every buffer passed as a pointer.
// stock_price and stock_price_sq point at column 0 of the S and S squared tables.
// Every load and store is unaligned (loadu and storeu).
void update_row_avx(const StencilConstants& constants, const double* current, double* next, int var_j, const double* stock_price, const double* stock_price_sq);

// Experiment E2. The same row with aligned loads and stores where the address allows.
// Column 1 of every row and of both tables must sit on a 32 byte boundary, or the program crashes.
// Centre, north, south, S, S squared and the store are aligned. The six east and west windows cannot be.
void update_row_avx_load(const StencilConstants& constants, const double* current, double* next, int var_j, const double* stock_price, const double* stock_price_sq);

// Experiment E2. One aligned load per row per four cells.
// The east and west windows are built in registers from neighbouring blocks with permutes and blends.
// Needs the same 32 byte alignment as update_row_avx_load.
void update_row_avx_shuffle(const StencilConstants& constants, const double* current, double* next, int var_j, const double* stock_price, const double* stock_price_sq);

// True when update_row_avx was built with AVX2 and FMA.
bool avx_row_is_vectorised();

// Updates the bottom row, where the variance is zero.
void update_zero_row(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next);
void update_zero_row(const StencilConstants& constants, const double* current, double* next);

// Writes the known values at S = 0 and at the largest S.
void apply_stock_edges(const StencilConstants& constants, std::vector<double>& next, int row_begin, int row_end, double time_to_expiry);
void apply_stock_edges(const StencilConstants& constants, double* next, int row_begin, int row_end, double time_to_expiry);

// Copies the second top row into the top row.
void copy_top_row(const StencilConstants& constants, std::vector<double>& next);
void copy_top_row(const StencilConstants& constants, double* next);

// The constants the float experiment needs, rounded down from the double ones.
// A float holds about 7 significant digits and a double about 16.
struct FloatStencilConstants {
    // The grid size.
    int num_stock_nodes = 0;
    int num_variance_nodes = 0;

    // The same scalars as StencilConstants, as floats.
    float dt = 0.0f;
    float rate = 0.0f;
    float stock_drift = 0.0f;
    float inv_two_ds = 0.0f;
    float inv_ds_sq = 0.0f;
    float inv_two_dv = 0.0f;
    float inv_dv_sq = 0.0f;
    float inv_four_ds_dv = 0.0f;
    float inv_dv = 0.0f;

    // The same tables as StencilConstants, as floats.
    std::vector<float> stock_price;
    std::vector<float> stock_price_sq;
    std::vector<float> stock_diffusion;
    std::vector<float> cross_diffusion;
    std::vector<float> variance_diffusion;
    std::vector<float> variance_drift;

    // Rounds every field of the double version, in stencil.cpp.
    FloatStencilConstants(const StencilConstants& constants);
};

// The float version of update_cells_scalar.
void update_cells_scalar_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j, int first_stock, int end_stock);

// Updates one interior variance row eight cells at a time with AVX2 floats.
// Defined in stencil_avx.cpp next to the double version.
void update_row_avx_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j);

#endif  // HESTON_STENCIL_H
