// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "stencil.h"

#include <vector>

// The AVX2 code is only built with -mavx2 -mfma.
// On the arm64 Mac the scalar fallback at the bottom of this file is used.
#if defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>

// Says the real AVX2 row was built into this program.
// The SIMD solver uses this to warn when it is running the scalar fallback.
bool avx_row_is_vectorised() {
    return true;
}

// Updates every interior cell in one variance row, four cells per instruction.
// A __m256d register holds four doubles, so each line works on four cells.
// current and next point at cell 0 of each buffer, the two tables at column 0.
void update_row_avx(const StencilConstants& constants, const double* current, double* next, int var_j, const double* stock_price_table, const double* stock_price_sq_table) {
    int num_stock_nodes = constants.num_stock_nodes;

    // Copy every constant into all four lanes of a register.
    // set1 fills all four slots with the same number.
    __m256d dt = _mm256_set1_pd(constants.dt);
    __m256d rate = _mm256_set1_pd(constants.rate);
    __m256d stock_drift = _mm256_set1_pd(constants.stock_drift);
    __m256d inv_two_ds = _mm256_set1_pd(constants.inv_two_ds);
    __m256d inv_ds_sq = _mm256_set1_pd(constants.inv_ds_sq);
    __m256d inv_two_dv = _mm256_set1_pd(constants.inv_two_dv);
    __m256d inv_dv_sq = _mm256_set1_pd(constants.inv_dv_sq);
    __m256d inv_four_ds_dv = _mm256_set1_pd(constants.inv_four_ds_dv);
    __m256d two = _mm256_set1_pd(2.0);

    // This row's four PDE coefficients, the same in every lane.
    __m256d stock_diffusion = _mm256_set1_pd(constants.stock_diffusion[var_j]);
    __m256d cross_diffusion = _mm256_set1_pd(constants.cross_diffusion[var_j]);
    __m256d variance_diffusion = _mm256_set1_pd(constants.variance_diffusion[var_j]);
    __m256d variance_drift = _mm256_set1_pd(constants.variance_drift[var_j]);

    // Pointers to the start of the rows to read and write.
    const double* row_centre = current + var_j * num_stock_nodes;
    const double* row_above = row_centre + num_stock_nodes;
    const double* row_below = row_centre - num_stock_nodes;
    double* out_row = next + var_j * num_stock_nodes;

    // The interior runs from column 1 to column num_stock_nodes - 2.
    int last_interior = num_stock_nodes - 2;
    int stock_i = 1;

    // Four cells per pass, for as long as four whole cells are left.
    while (stock_i + 3 <= last_interior) {
        // Load four neighbouring cells and their eight neighbours.
        // loadu allows an address that is not a multiple of 32 bytes.
        __m256d centre = _mm256_loadu_pd(row_centre + stock_i);
        __m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
        __m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);
        __m256d north = _mm256_loadu_pd(row_above + stock_i);
        __m256d south = _mm256_loadu_pd(row_below + stock_i);
        __m256d north_east = _mm256_loadu_pd(row_above + stock_i + 1);
        __m256d north_west = _mm256_loadu_pd(row_above + stock_i - 1);
        __m256d south_east = _mm256_loadu_pd(row_below + stock_i + 1);
        __m256d south_west = _mm256_loadu_pd(row_below + stock_i - 1);

        // First derivatives.
        __m256d east_minus_west = _mm256_sub_pd(east, west);
        __m256d dV_dS = _mm256_mul_pd(east_minus_west, inv_two_ds);
        __m256d north_minus_south = _mm256_sub_pd(north, south);
        __m256d dV_dv = _mm256_mul_pd(north_minus_south, inv_two_dv);

        // Second derivatives, the two neighbours minus twice the centre.
        // fnmadd(a, b, c) works out c minus a times b in one instruction.
        __m256d east_plus_west = _mm256_add_pd(east, west);
        __m256d stock_curvature = _mm256_fnmadd_pd(two, centre, east_plus_west);
        __m256d d2V_dS2 = _mm256_mul_pd(stock_curvature, inv_ds_sq);
        __m256d north_plus_south = _mm256_add_pd(north, south);
        __m256d variance_curvature = _mm256_fnmadd_pd(two, centre, north_plus_south);
        __m256d d2V_dv2 = _mm256_mul_pd(variance_curvature, inv_dv_sq);

        // Cross derivative from the four corner cells.
        __m256d corners = _mm256_sub_pd(north_east, north_west);
        corners = _mm256_sub_pd(corners, south_east);
        corners = _mm256_add_pd(corners, south_west);
        __m256d d2V_dS_dv = _mm256_mul_pd(corners, inv_four_ds_dv);

        // Load S and S squared for these four columns.
        __m256d stock_price = _mm256_loadu_pd(stock_price_table + stock_i);
        __m256d stock_price_sq = _mm256_loadu_pd(stock_price_sq_table + stock_i);

        // Multiply each PDE coefficient by the S factor its term needs.
        __m256d stock_diffusion_factor = _mm256_mul_pd(stock_diffusion, stock_price_sq);
        __m256d cross_factor = _mm256_mul_pd(cross_diffusion, stock_price);
        __m256d stock_drift_factor = _mm256_mul_pd(stock_drift, stock_price);

        // Add up the five PDE terms and the discounting.
        // fmadd(a, b, c) works out a times b plus c with one rounding.
        __m256d change = _mm256_mul_pd(stock_diffusion_factor, d2V_dS2);
        change = _mm256_fmadd_pd(cross_factor, d2V_dS_dv, change);
        change = _mm256_fmadd_pd(variance_diffusion, d2V_dv2, change);
        change = _mm256_fmadd_pd(stock_drift_factor, dV_dS, change);
        change = _mm256_fmadd_pd(variance_drift, dV_dv, change);
        change = _mm256_fnmadd_pd(rate, centre, change);

        // Step back one dt and store the four new values.
        __m256d new_value = _mm256_fmadd_pd(dt, change, centre);
        _mm256_storeu_pd(out_row + stock_i, new_value);

        stock_i = stock_i + 4;
    }

    // The scalar row does the last zero to three cells.
    update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
}

// The same row with the buffers passed as vectors, as every rung calls it.
// .data() is the address of cell 0.
void update_row_avx(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    update_row_avx(constants, current.data(), next.data(), var_j, constants.stock_price.data(), constants.stock_price_sq.data());
}

// Experiment E2. The row above with the aligned load and store wherever the address allows.
// _mm256_load_pd and _mm256_store_pd crash on an address that is not a multiple of 32 bytes.
// Column 1 is on a 32 byte boundary, so the windows starting at stock_i are aligned and
// the windows starting at stock_i - 1 and stock_i + 1 are 8 bytes off and stay loadu.
void update_row_avx_load(const StencilConstants& constants, const double* current, double* next, int var_j, const double* stock_price_table, const double* stock_price_sq_table) {
    int num_stock_nodes = constants.num_stock_nodes;

    // The constants, the same as update_row_avx.
    __m256d dt = _mm256_set1_pd(constants.dt);
    __m256d rate = _mm256_set1_pd(constants.rate);
    __m256d stock_drift = _mm256_set1_pd(constants.stock_drift);
    __m256d inv_two_ds = _mm256_set1_pd(constants.inv_two_ds);
    __m256d inv_ds_sq = _mm256_set1_pd(constants.inv_ds_sq);
    __m256d inv_two_dv = _mm256_set1_pd(constants.inv_two_dv);
    __m256d inv_dv_sq = _mm256_set1_pd(constants.inv_dv_sq);
    __m256d inv_four_ds_dv = _mm256_set1_pd(constants.inv_four_ds_dv);
    __m256d two = _mm256_set1_pd(2.0);
    __m256d stock_diffusion = _mm256_set1_pd(constants.stock_diffusion[var_j]);
    __m256d cross_diffusion = _mm256_set1_pd(constants.cross_diffusion[var_j]);
    __m256d variance_diffusion = _mm256_set1_pd(constants.variance_diffusion[var_j]);
    __m256d variance_drift = _mm256_set1_pd(constants.variance_drift[var_j]);

    const double* row_centre = current + var_j * num_stock_nodes;
    const double* row_above = row_centre + num_stock_nodes;
    const double* row_below = row_centre - num_stock_nodes;
    double* out_row = next + var_j * num_stock_nodes;

    int last_interior = num_stock_nodes - 2;
    int stock_i = 1;

    while (stock_i + 3 <= last_interior) {
        // Aligned, the three windows that start at stock_i.
        __m256d centre = _mm256_load_pd(row_centre + stock_i);
        __m256d north = _mm256_load_pd(row_above + stock_i);
        __m256d south = _mm256_load_pd(row_below + stock_i);
        // Unaligned, the six windows one cell either side.
        __m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
        __m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);
        __m256d north_east = _mm256_loadu_pd(row_above + stock_i + 1);
        __m256d north_west = _mm256_loadu_pd(row_above + stock_i - 1);
        __m256d south_east = _mm256_loadu_pd(row_below + stock_i + 1);
        __m256d south_west = _mm256_loadu_pd(row_below + stock_i - 1);

        // From here to the store, the same sums in the same order as update_row_avx.
        __m256d east_minus_west = _mm256_sub_pd(east, west);
        __m256d dV_dS = _mm256_mul_pd(east_minus_west, inv_two_ds);
        __m256d north_minus_south = _mm256_sub_pd(north, south);
        __m256d dV_dv = _mm256_mul_pd(north_minus_south, inv_two_dv);

        __m256d east_plus_west = _mm256_add_pd(east, west);
        __m256d stock_curvature = _mm256_fnmadd_pd(two, centre, east_plus_west);
        __m256d d2V_dS2 = _mm256_mul_pd(stock_curvature, inv_ds_sq);
        __m256d north_plus_south = _mm256_add_pd(north, south);
        __m256d variance_curvature = _mm256_fnmadd_pd(two, centre, north_plus_south);
        __m256d d2V_dv2 = _mm256_mul_pd(variance_curvature, inv_dv_sq);

        __m256d corners = _mm256_sub_pd(north_east, north_west);
        corners = _mm256_sub_pd(corners, south_east);
        corners = _mm256_add_pd(corners, south_west);
        __m256d d2V_dS_dv = _mm256_mul_pd(corners, inv_four_ds_dv);

        // The tables start at the same alignment as the rows, so these are aligned too.
        __m256d stock_price = _mm256_load_pd(stock_price_table + stock_i);
        __m256d stock_price_sq = _mm256_load_pd(stock_price_sq_table + stock_i);

        __m256d stock_diffusion_factor = _mm256_mul_pd(stock_diffusion, stock_price_sq);
        __m256d cross_factor = _mm256_mul_pd(cross_diffusion, stock_price);
        __m256d stock_drift_factor = _mm256_mul_pd(stock_drift, stock_price);

        __m256d change = _mm256_mul_pd(stock_diffusion_factor, d2V_dS2);
        change = _mm256_fmadd_pd(cross_factor, d2V_dS_dv, change);
        change = _mm256_fmadd_pd(variance_diffusion, d2V_dv2, change);
        change = _mm256_fmadd_pd(stock_drift_factor, dV_dS, change);
        change = _mm256_fmadd_pd(variance_drift, dV_dv, change);
        change = _mm256_fnmadd_pd(rate, centre, change);

        // The store starts at stock_i, so it is aligned.
        __m256d new_value = _mm256_fmadd_pd(dt, change, centre);
        _mm256_store_pd(out_row + stock_i, new_value);

        stock_i = stock_i + 4;
    }

    update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
}

// Experiment E2. Each row is read once, in aligned blocks of four, and never at stock_i - 1 or stock_i + 1.
// The west window [c-1, c, c+1, c+2] is the last cell of the block before plus the first three of this block.
// The east window [c+1, c+2, c+3, c+4] is the last three of this block plus the first of the block after.
// permute4x64 rotates the four lanes of one register. blend picks each lane from one of two registers.
void update_row_avx_shuffle(const StencilConstants& constants, const double* current, double* next, int var_j, const double* stock_price_table, const double* stock_price_sq_table) {
    int num_stock_nodes = constants.num_stock_nodes;

    // The constants, the same as update_row_avx.
    __m256d dt = _mm256_set1_pd(constants.dt);
    __m256d rate = _mm256_set1_pd(constants.rate);
    __m256d stock_drift = _mm256_set1_pd(constants.stock_drift);
    __m256d inv_two_ds = _mm256_set1_pd(constants.inv_two_ds);
    __m256d inv_ds_sq = _mm256_set1_pd(constants.inv_ds_sq);
    __m256d inv_two_dv = _mm256_set1_pd(constants.inv_two_dv);
    __m256d inv_dv_sq = _mm256_set1_pd(constants.inv_dv_sq);
    __m256d inv_four_ds_dv = _mm256_set1_pd(constants.inv_four_ds_dv);
    __m256d two = _mm256_set1_pd(2.0);
    __m256d stock_diffusion = _mm256_set1_pd(constants.stock_diffusion[var_j]);
    __m256d cross_diffusion = _mm256_set1_pd(constants.cross_diffusion[var_j]);
    __m256d variance_diffusion = _mm256_set1_pd(constants.variance_diffusion[var_j]);
    __m256d variance_drift = _mm256_set1_pd(constants.variance_drift[var_j]);

    const double* row_centre = current + var_j * num_stock_nodes;
    const double* row_above = row_centre + num_stock_nodes;
    const double* row_below = row_centre - num_stock_nodes;
    double* out_row = next + var_j * num_stock_nodes;

    int last_interior = num_stock_nodes - 2;
    int stock_i = 1;

    // The lane orders for permute4x64. _MM_SHUFFLE lists the source lane for lanes 3, 2, 1, 0.
    // Rotate right gives [x3, x0, x1, x2]. Rotate left gives [x1, x2, x3, x0].
    const int rotate_right = _MM_SHUFFLE(2, 1, 0, 3);
    const int rotate_left = _MM_SHUFFLE(0, 3, 2, 1);

    // Before the first pass, load this block of each row and rotate it both ways.
    // The block before is rotated right, and only its lane 0 (the cell at stock_i - 1) is used.
    // At var_j = 1 the row below is row 0, so that block reads 3 cells before cell 0.
    // The solver leaves padding there.
    __m256d centre_block = _mm256_load_pd(row_centre + stock_i);
    __m256d above_block = _mm256_load_pd(row_above + stock_i);
    __m256d below_block = _mm256_load_pd(row_below + stock_i);
    __m256d centre_left = _mm256_permute4x64_pd(centre_block, rotate_left);
    __m256d above_left = _mm256_permute4x64_pd(above_block, rotate_left);
    __m256d below_left = _mm256_permute4x64_pd(below_block, rotate_left);
    __m256d centre_prev_right = _mm256_permute4x64_pd(_mm256_load_pd(row_centre + stock_i - 4), rotate_right);
    __m256d above_prev_right = _mm256_permute4x64_pd(_mm256_load_pd(row_above + stock_i - 4), rotate_right);
    __m256d below_prev_right = _mm256_permute4x64_pd(_mm256_load_pd(row_below + stock_i - 4), rotate_right);

    while (stock_i + 3 <= last_interior) {
        // The block after, the only new loads this pass.
        // On the last pass of the top row this reads up to 3 cells past the grid. The solver leaves padding there.
        __m256d centre_next = _mm256_load_pd(row_centre + stock_i + 4);
        __m256d above_next = _mm256_load_pd(row_above + stock_i + 4);
        __m256d below_next = _mm256_load_pd(row_below + stock_i + 4);

        // Rotate this block right and the next block left.
        __m256d centre_right = _mm256_permute4x64_pd(centre_block, rotate_right);
        __m256d above_right = _mm256_permute4x64_pd(above_block, rotate_right);
        __m256d below_right = _mm256_permute4x64_pd(below_block, rotate_right);
        __m256d centre_next_left = _mm256_permute4x64_pd(centre_next, rotate_left);
        __m256d above_next_left = _mm256_permute4x64_pd(above_next, rotate_left);
        __m256d below_next_left = _mm256_permute4x64_pd(below_next, rotate_left);

        // West: lane 0 from the previous block, lanes 1 to 3 from this one. Bit i of the mask picks the second register.
        // East: lanes 0 to 2 from this block, lane 3 from the next one.
        __m256d centre = centre_block;
        __m256d north = above_block;
        __m256d south = below_block;
        __m256d west = _mm256_blend_pd(centre_right, centre_prev_right, 0b0001);
        __m256d east = _mm256_blend_pd(centre_left, centre_next_left, 0b1000);
        __m256d north_west = _mm256_blend_pd(above_right, above_prev_right, 0b0001);
        __m256d north_east = _mm256_blend_pd(above_left, above_next_left, 0b1000);
        __m256d south_west = _mm256_blend_pd(below_right, below_prev_right, 0b0001);
        __m256d south_east = _mm256_blend_pd(below_left, below_next_left, 0b1000);

        // From here to the store, the same sums in the same order as update_row_avx.
        __m256d east_minus_west = _mm256_sub_pd(east, west);
        __m256d dV_dS = _mm256_mul_pd(east_minus_west, inv_two_ds);
        __m256d north_minus_south = _mm256_sub_pd(north, south);
        __m256d dV_dv = _mm256_mul_pd(north_minus_south, inv_two_dv);

        __m256d east_plus_west = _mm256_add_pd(east, west);
        __m256d stock_curvature = _mm256_fnmadd_pd(two, centre, east_plus_west);
        __m256d d2V_dS2 = _mm256_mul_pd(stock_curvature, inv_ds_sq);
        __m256d north_plus_south = _mm256_add_pd(north, south);
        __m256d variance_curvature = _mm256_fnmadd_pd(two, centre, north_plus_south);
        __m256d d2V_dv2 = _mm256_mul_pd(variance_curvature, inv_dv_sq);

        __m256d corners = _mm256_sub_pd(north_east, north_west);
        corners = _mm256_sub_pd(corners, south_east);
        corners = _mm256_add_pd(corners, south_west);
        __m256d d2V_dS_dv = _mm256_mul_pd(corners, inv_four_ds_dv);

        __m256d stock_price = _mm256_load_pd(stock_price_table + stock_i);
        __m256d stock_price_sq = _mm256_load_pd(stock_price_sq_table + stock_i);

        __m256d stock_diffusion_factor = _mm256_mul_pd(stock_diffusion, stock_price_sq);
        __m256d cross_factor = _mm256_mul_pd(cross_diffusion, stock_price);
        __m256d stock_drift_factor = _mm256_mul_pd(stock_drift, stock_price);

        __m256d change = _mm256_mul_pd(stock_diffusion_factor, d2V_dS2);
        change = _mm256_fmadd_pd(cross_factor, d2V_dS_dv, change);
        change = _mm256_fmadd_pd(variance_diffusion, d2V_dv2, change);
        change = _mm256_fmadd_pd(stock_drift_factor, dV_dS, change);
        change = _mm256_fmadd_pd(variance_drift, dV_dv, change);
        change = _mm256_fnmadd_pd(rate, centre, change);

        __m256d new_value = _mm256_fmadd_pd(dt, change, centre);
        _mm256_store_pd(out_row + stock_i, new_value);

        // Move one block along. This block's right rotation is the next pass's previous block,
        // and the next block's left rotation is the next pass's own.
        centre_prev_right = centre_right;
        above_prev_right = above_right;
        below_prev_right = below_right;
        centre_block = centre_next;
        above_block = above_next;
        below_block = below_next;
        centre_left = centre_next_left;
        above_left = above_next_left;
        below_left = below_next_left;

        stock_i = stock_i + 4;
    }

    update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
}

// Updates every interior cell in one variance row, eight float cells per instruction.
// A __m256 register holds eight floats in 256 bits.
void update_row_avx_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j) {
    int num_stock_nodes = constants.num_stock_nodes;

    // The _ps functions are the float versions of the _pd ones.
    __m256 dt = _mm256_set1_ps(constants.dt);
    __m256 rate = _mm256_set1_ps(constants.rate);
    __m256 stock_drift = _mm256_set1_ps(constants.stock_drift);
    __m256 inv_two_ds = _mm256_set1_ps(constants.inv_two_ds);
    __m256 inv_ds_sq = _mm256_set1_ps(constants.inv_ds_sq);
    __m256 inv_two_dv = _mm256_set1_ps(constants.inv_two_dv);
    __m256 inv_dv_sq = _mm256_set1_ps(constants.inv_dv_sq);
    __m256 inv_four_ds_dv = _mm256_set1_ps(constants.inv_four_ds_dv);
    __m256 two = _mm256_set1_ps(2.0f);

    // This row's four PDE coefficients, the same in every lane.
    __m256 stock_diffusion = _mm256_set1_ps(constants.stock_diffusion[var_j]);
    __m256 cross_diffusion = _mm256_set1_ps(constants.cross_diffusion[var_j]);
    __m256 variance_diffusion = _mm256_set1_ps(constants.variance_diffusion[var_j]);
    __m256 variance_drift = _mm256_set1_ps(constants.variance_drift[var_j]);

    // Pointers to the start of the rows to read and write.
    const float* stock_price_table = constants.stock_price.data();
    const float* stock_price_sq_table = constants.stock_price_sq.data();
    const float* row_centre = current.data() + var_j * num_stock_nodes;
    const float* row_above = row_centre + num_stock_nodes;
    const float* row_below = row_centre - num_stock_nodes;
    float* out_row = next.data() + var_j * num_stock_nodes;

    int last_interior = num_stock_nodes - 2;
    int stock_i = 1;

    // Eight cells per pass, for as long as eight whole cells are left.
    while (stock_i + 7 <= last_interior) {
        // Load eight neighbouring cells and their neighbours.
        __m256 centre = _mm256_loadu_ps(row_centre + stock_i);
        __m256 east = _mm256_loadu_ps(row_centre + stock_i + 1);
        __m256 west = _mm256_loadu_ps(row_centre + stock_i - 1);
        __m256 north = _mm256_loadu_ps(row_above + stock_i);
        __m256 south = _mm256_loadu_ps(row_below + stock_i);
        __m256 north_east = _mm256_loadu_ps(row_above + stock_i + 1);
        __m256 north_west = _mm256_loadu_ps(row_above + stock_i - 1);
        __m256 south_east = _mm256_loadu_ps(row_below + stock_i + 1);
        __m256 south_west = _mm256_loadu_ps(row_below + stock_i - 1);

        // First derivatives.
        __m256 east_minus_west = _mm256_sub_ps(east, west);
        __m256 dV_dS = _mm256_mul_ps(east_minus_west, inv_two_ds);
        __m256 north_minus_south = _mm256_sub_ps(north, south);
        __m256 dV_dv = _mm256_mul_ps(north_minus_south, inv_two_dv);

        // Second derivatives, the two neighbours minus twice the centre.
        __m256 east_plus_west = _mm256_add_ps(east, west);
        __m256 stock_curvature = _mm256_fnmadd_ps(two, centre, east_plus_west);
        __m256 d2V_dS2 = _mm256_mul_ps(stock_curvature, inv_ds_sq);
        __m256 north_plus_south = _mm256_add_ps(north, south);
        __m256 variance_curvature = _mm256_fnmadd_ps(two, centre, north_plus_south);
        __m256 d2V_dv2 = _mm256_mul_ps(variance_curvature, inv_dv_sq);

        // Cross derivative from the four corner cells.
        __m256 corners = _mm256_sub_ps(north_east, north_west);
        corners = _mm256_sub_ps(corners, south_east);
        corners = _mm256_add_ps(corners, south_west);
        __m256 d2V_dS_dv = _mm256_mul_ps(corners, inv_four_ds_dv);

        // Load S and S squared for these eight columns.
        __m256 stock_price = _mm256_loadu_ps(stock_price_table + stock_i);
        __m256 stock_price_sq = _mm256_loadu_ps(stock_price_sq_table + stock_i);

        // Multiply each PDE coefficient by the S factor its term needs.
        __m256 stock_diffusion_factor = _mm256_mul_ps(stock_diffusion, stock_price_sq);
        __m256 cross_factor = _mm256_mul_ps(cross_diffusion, stock_price);
        __m256 stock_drift_factor = _mm256_mul_ps(stock_drift, stock_price);

        // Add up the five PDE terms and the discounting.
        __m256 change = _mm256_mul_ps(stock_diffusion_factor, d2V_dS2);
        change = _mm256_fmadd_ps(cross_factor, d2V_dS_dv, change);
        change = _mm256_fmadd_ps(variance_diffusion, d2V_dv2, change);
        change = _mm256_fmadd_ps(stock_drift_factor, dV_dS, change);
        change = _mm256_fmadd_ps(variance_drift, dV_dv, change);
        change = _mm256_fnmadd_ps(rate, centre, change);

        // Step back one dt and store the eight new values.
        __m256 new_value = _mm256_fmadd_ps(dt, change, centre);
        _mm256_storeu_ps(out_row + stock_i, new_value);

        stock_i = stock_i + 8;
    }

    // The scalar float row does the last zero to seven cells.
    update_cells_scalar_float(constants, current, next, var_j, stock_i, last_interior + 1);
}

#else

// Says the real AVX2 row was not built into this program.
// The SIMD solver uses this to warn when it is running the scalar fallback.
bool avx_row_is_vectorised() {
    return false;
}

// Stands in for the AVX2 row on machines without AVX2.
// Gives the serial answer so the tests run on the Mac.
void update_row_avx(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    update_row_scalar(constants, current, next, var_j);
}

// The three pointer rows of experiment E2 also fall back to the scalar row.
void update_row_avx(const StencilConstants& constants, const double* current, double* next, int var_j, const double*, const double*) {
    update_cells_scalar(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}

void update_row_avx_load(const StencilConstants& constants, const double* current, double* next, int var_j, const double*, const double*) {
    update_cells_scalar(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}

void update_row_avx_shuffle(const StencilConstants& constants, const double* current, double* next, int var_j, const double*, const double*) {
    update_cells_scalar(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}

// Stands in for the AVX2 float row on machines without AVX2.
// Runs the float scalar row over the whole row.
void update_row_avx_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j) {
    update_cells_scalar_float(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}

#endif
