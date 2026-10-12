// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "stencil.h"

#include <vector>

// The AVX2 code is only built when the compiler was given -mavx2 -mfma.
// On the arm64 Mac neither is set, so the scalar fallback at the bottom is used
#if defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>

// Says the real AVX2 row was built into this program.
// The SIMD solver uses this to warn when it is running the scalar fallback.
bool avx_row_is_vectorised() {
    return true;
}

// Updates every interior cell in one variance row, four cells per instruction.
// A __m256d register holds four doubles, so each line below works on four neighbouring cells.
void update_row_avx(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    int num_stock_nodes = constants.num_stock_nodes;

    // Copy every constant into all four lanes of a register.
    // set1 is like filling a four slot array with the same number
    __m256d dt = _mm256_set1_pd(constants.dt);
    __m256d rate = _mm256_set1_pd(constants.rate);
    __m256d stock_drift = _mm256_set1_pd(constants.stock_drift);
    __m256d inv_two_ds = _mm256_set1_pd(constants.inv_two_ds);
    __m256d inv_ds_sq = _mm256_set1_pd(constants.inv_ds_sq);
    __m256d inv_two_dv = _mm256_set1_pd(constants.inv_two_dv);
    __m256d inv_dv_sq = _mm256_set1_pd(constants.inv_dv_sq);
    __m256d inv_four_ds_dv = _mm256_set1_pd(constants.inv_four_ds_dv);
    __m256d two = _mm256_set1_pd(2.0);

    // This row's four PDE coefficients, the same in every lane
    __m256d stock_diffusion = _mm256_set1_pd(constants.stock_diffusion[var_j]);
    __m256d cross_diffusion = _mm256_set1_pd(constants.cross_diffusion[var_j]);
    __m256d variance_diffusion = _mm256_set1_pd(constants.variance_diffusion[var_j]);
    __m256d variance_drift = _mm256_set1_pd(constants.variance_drift[var_j]);

    // Pointers to the start of the rows we read and write, as in the scalar row
    const double* stock_price_table = constants.stock_price.data();
    const double* stock_price_sq_table = constants.stock_price_sq.data();
    const double* row_centre = current.data() + var_j * num_stock_nodes;
    const double* row_above = row_centre + num_stock_nodes;
    const double* row_below = row_centre - num_stock_nodes;
    double* out_row = next.data() + var_j * num_stock_nodes;

    // The interior runs from column 1 to column num_stock_nodes - 2
    int last_interior = num_stock_nodes - 2;
    int stock_i = 1;

    // Four cells per pass, for as long as four whole cells are left
    while (stock_i + 3 <= last_interior) {
        // Load four neighbouring cells and their eight neighbours.
        // loadu means the address does not have to be a multiple of 32 bytes
        __m256d centre = _mm256_loadu_pd(row_centre + stock_i);
        __m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
        __m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);
        __m256d north = _mm256_loadu_pd(row_above + stock_i);
        __m256d south = _mm256_loadu_pd(row_below + stock_i);
        __m256d north_east = _mm256_loadu_pd(row_above + stock_i + 1);
        __m256d north_west = _mm256_loadu_pd(row_above + stock_i - 1);
        __m256d south_east = _mm256_loadu_pd(row_below + stock_i + 1);
        __m256d south_west = _mm256_loadu_pd(row_below + stock_i - 1);

        // First derivatives, the gap across the cell times one over twice the spacing
        __m256d east_minus_west = _mm256_sub_pd(east, west);
        __m256d dV_dS = _mm256_mul_pd(east_minus_west, inv_two_ds);
        __m256d north_minus_south = _mm256_sub_pd(north, south);
        __m256d dV_dv = _mm256_mul_pd(north_minus_south, inv_two_dv);

        // Second derivatives, the two neighbours minus twice the centre.
        // fnmadd(a, b, c) works out c minus a times b in one instruction
        __m256d east_plus_west = _mm256_add_pd(east, west);
        __m256d stock_curvature = _mm256_fnmadd_pd(two, centre, east_plus_west);
        __m256d d2V_dS2 = _mm256_mul_pd(stock_curvature, inv_ds_sq);
        __m256d north_plus_south = _mm256_add_pd(north, south);
        __m256d variance_curvature = _mm256_fnmadd_pd(two, centre, north_plus_south);
        __m256d d2V_dv2 = _mm256_mul_pd(variance_curvature, inv_dv_sq);

        // Cross derivative from the four corner cells
        __m256d corners = _mm256_sub_pd(north_east, north_west);
        corners = _mm256_sub_pd(corners, south_east);
        corners = _mm256_add_pd(corners, south_west);
        __m256d d2V_dS_dv = _mm256_mul_pd(corners, inv_four_ds_dv);

        // S and S squared for these four columns, from the tables
        __m256d stock_price = _mm256_loadu_pd(stock_price_table + stock_i);
        __m256d stock_price_sq = _mm256_loadu_pd(stock_price_sq_table + stock_i);

        // Each PDE coefficient times the S factor its term needs
        __m256d stock_diffusion_factor = _mm256_mul_pd(stock_diffusion, stock_price_sq);
        __m256d cross_factor = _mm256_mul_pd(cross_diffusion, stock_price);
        __m256d stock_drift_factor = _mm256_mul_pd(stock_drift, stock_price);

        // Add up the five PDE terms and the discounting.
        // fmadd(a, b, c) works out a times b plus c with only one rounding
        __m256d change = _mm256_mul_pd(stock_diffusion_factor, d2V_dS2);
        change = _mm256_fmadd_pd(cross_factor, d2V_dS_dv, change);
        change = _mm256_fmadd_pd(variance_diffusion, d2V_dv2, change);
        change = _mm256_fmadd_pd(stock_drift_factor, dV_dS, change);
        change = _mm256_fmadd_pd(variance_drift, dV_dv, change);
        change = _mm256_fnmadd_pd(rate, centre, change);

        // Step back one dt and store the four new values
        __m256d new_value = _mm256_fmadd_pd(dt, change, centre);
        _mm256_storeu_pd(out_row + stock_i, new_value);

        stock_i = stock_i + 4;
    }

    // The last zero to three cells do not fill a register, so do them one at a time
    update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
}

// Updates every interior cell in one variance row, eight float cells per instruction.
// A __m256 register holds eight floats in the same 256 bits that hold four doubles.
void update_row_avx_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j) {
    int num_stock_nodes = constants.num_stock_nodes;

    // The _ps functions are the float versions of the _pd ones above
    __m256 dt = _mm256_set1_ps(constants.dt);
    __m256 rate = _mm256_set1_ps(constants.rate);
    __m256 stock_drift = _mm256_set1_ps(constants.stock_drift);
    __m256 inv_two_ds = _mm256_set1_ps(constants.inv_two_ds);
    __m256 inv_ds_sq = _mm256_set1_ps(constants.inv_ds_sq);
    __m256 inv_two_dv = _mm256_set1_ps(constants.inv_two_dv);
    __m256 inv_dv_sq = _mm256_set1_ps(constants.inv_dv_sq);
    __m256 inv_four_ds_dv = _mm256_set1_ps(constants.inv_four_ds_dv);
    __m256 two = _mm256_set1_ps(2.0f);

    // This row's four PDE coefficients, the same in every lane
    __m256 stock_diffusion = _mm256_set1_ps(constants.stock_diffusion[var_j]);
    __m256 cross_diffusion = _mm256_set1_ps(constants.cross_diffusion[var_j]);
    __m256 variance_diffusion = _mm256_set1_ps(constants.variance_diffusion[var_j]);
    __m256 variance_drift = _mm256_set1_ps(constants.variance_drift[var_j]);

    // Pointers to the start of the rows we read and write
    const float* stock_price_table = constants.stock_price.data();
    const float* stock_price_sq_table = constants.stock_price_sq.data();
    const float* row_centre = current.data() + var_j * num_stock_nodes;
    const float* row_above = row_centre + num_stock_nodes;
    const float* row_below = row_centre - num_stock_nodes;
    float* out_row = next.data() + var_j * num_stock_nodes;

    int last_interior = num_stock_nodes - 2;
    int stock_i = 1;

    // Eight cells per pass, for as long as eight whole cells are left
    while (stock_i + 7 <= last_interior) {
        // Load eight neighbouring cells and their neighbours
        __m256 centre = _mm256_loadu_ps(row_centre + stock_i);
        __m256 east = _mm256_loadu_ps(row_centre + stock_i + 1);
        __m256 west = _mm256_loadu_ps(row_centre + stock_i - 1);
        __m256 north = _mm256_loadu_ps(row_above + stock_i);
        __m256 south = _mm256_loadu_ps(row_below + stock_i);
        __m256 north_east = _mm256_loadu_ps(row_above + stock_i + 1);
        __m256 north_west = _mm256_loadu_ps(row_above + stock_i - 1);
        __m256 south_east = _mm256_loadu_ps(row_below + stock_i + 1);
        __m256 south_west = _mm256_loadu_ps(row_below + stock_i - 1);

        // First derivatives
        __m256 east_minus_west = _mm256_sub_ps(east, west);
        __m256 dV_dS = _mm256_mul_ps(east_minus_west, inv_two_ds);
        __m256 north_minus_south = _mm256_sub_ps(north, south);
        __m256 dV_dv = _mm256_mul_ps(north_minus_south, inv_two_dv);

        // Second derivatives, the two neighbours minus twice the centre
        __m256 east_plus_west = _mm256_add_ps(east, west);
        __m256 stock_curvature = _mm256_fnmadd_ps(two, centre, east_plus_west);
        __m256 d2V_dS2 = _mm256_mul_ps(stock_curvature, inv_ds_sq);
        __m256 north_plus_south = _mm256_add_ps(north, south);
        __m256 variance_curvature = _mm256_fnmadd_ps(two, centre, north_plus_south);
        __m256 d2V_dv2 = _mm256_mul_ps(variance_curvature, inv_dv_sq);

        // Cross derivative from the four corner cells
        __m256 corners = _mm256_sub_ps(north_east, north_west);
        corners = _mm256_sub_ps(corners, south_east);
        corners = _mm256_add_ps(corners, south_west);
        __m256 d2V_dS_dv = _mm256_mul_ps(corners, inv_four_ds_dv);

        // S and S squared for these eight columns
        __m256 stock_price = _mm256_loadu_ps(stock_price_table + stock_i);
        __m256 stock_price_sq = _mm256_loadu_ps(stock_price_sq_table + stock_i);

        // Each PDE coefficient times the S factor its term needs
        __m256 stock_diffusion_factor = _mm256_mul_ps(stock_diffusion, stock_price_sq);
        __m256 cross_factor = _mm256_mul_ps(cross_diffusion, stock_price);
        __m256 stock_drift_factor = _mm256_mul_ps(stock_drift, stock_price);

        // Add up the five PDE terms and the discounting
        __m256 change = _mm256_mul_ps(stock_diffusion_factor, d2V_dS2);
        change = _mm256_fmadd_ps(cross_factor, d2V_dS_dv, change);
        change = _mm256_fmadd_ps(variance_diffusion, d2V_dv2, change);
        change = _mm256_fmadd_ps(stock_drift_factor, dV_dS, change);
        change = _mm256_fmadd_ps(variance_drift, dV_dv, change);
        change = _mm256_fnmadd_ps(rate, centre, change);

        // Step back one dt and store the eight new values
        __m256 new_value = _mm256_fmadd_ps(dt, change, centre);
        _mm256_storeu_ps(out_row + stock_i, new_value);

        stock_i = stock_i + 8;
    }

    // The last zero to seven cells do not fill a register, so do them one at a time
    update_cells_scalar_float(constants, current, next, var_j, stock_i, last_interior + 1);
}

#else

// Says the real AVX2 row was not built into this program.
// The SIMD solver uses this to warn when it is running the scalar fallback.
bool avx_row_is_vectorised() {
    return false;
}

// Stands in for the AVX2 row on machines without AVX2, such as the arm64 Mac.
// It gives exactly the serial answer, so the tests still run here.
void update_row_avx(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    update_row_scalar(constants, current, next, var_j);
}

// Stands in for the AVX2 float row on machines without AVX2.
// It runs the float scalar row over the whole row.
void update_row_avx_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j) {
    update_cells_scalar_float(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}

#endif
