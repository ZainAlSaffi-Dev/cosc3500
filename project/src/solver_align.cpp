// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <omp.h>

#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "io.h"
#include "solver_experiments.h"
#include "stencil.h"

// Spare doubles kept before cell 0 and after the last cell.
// The shuffle row reads up to 3 cells past each end of the grid.
// The front also gives up to 7 more doubles of room to move cell 0 to the chosen alignment.
static const int front_padding = 16;
static const int back_padding = 8;

// One buffer whose cell 0 starts part way into a vector, so column 1 lands where we choose.
// The vector owns the memory, start is cell 0.
struct ShiftedBuffer {
    std::vector<double> storage;
    double* start = nullptr;
};

// Makes a zeroed buffer of num_cells doubles with column 1 of row 0 align_shift doubles past a 64 byte boundary.
// Every row starts at the same place in its cache line when the row length is a multiple of 8 doubles.
static void make_shifted_buffer(ShiftedBuffer& buffer, int num_cells, int align_shift) {
    // The vector is big, so like the Grid's buffers it comes from mmap with the same page regime.
    buffer.storage.assign(num_cells + front_padding + back_padding, 0.0);

    // Try the 8 places from double 8 on until column 1 is the right number of bytes into its cache line.
    // uintptr_t is an integer as wide as an address, so % works on it.
    std::uintptr_t wanted_byte = 8 * (std::uintptr_t)align_shift;
    for (int offset = 8; offset < 16; ++offset) {
        double* cell_zero = buffer.storage.data() + offset;
        std::uintptr_t column_one = (std::uintptr_t)(cell_zero + 1);
        if (column_one % 64 == wanted_byte) {
            buffer.start = cell_zero;
            return;
        }
    }
    throw std::runtime_error("could not place the buffer, the vector is not 8 byte aligned");
}

// How many bytes past a 64 byte boundary an address sits.
static int byte_in_cache_line(const double* address) {
    return (int)((std::uintptr_t)address % 64);
}

std::string OmpAlignSolver::name() const {
    if (kernel_ == AlignKernel::load) {
        return "omp_align_load";
    }
    if (kernel_ == AlignKernel::shuffle) {
        return "omp_align_shuffle";
    }
    return "omp_align_loadu";
}

// Prices the option like the omp solver, on buffers whose alignment is chosen by --align-shift.
// Same row updates, same boundaries, same order, so the answer matches simd bit for bit.
SolveResult OmpAlignSolver::solve(const Config& cfg) {
    SolveResult result;
    result.threads = omp_get_max_threads();

    if (cfg.dump_every > 0) {
        throw std::runtime_error("--dump-every is not supported by the " + name() + " solver");
    }

    // Every row must start at the same place in its cache line, so the row length must be a multiple of 8 doubles (64 bytes).
    int align_shift = cfg.align_shift;
    if (cfg.grid.num_stock_nodes % 8 != 0) {
        throw std::runtime_error(name() + " needs --ns to be a multiple of 8");
    }
    if (align_shift < 0 || align_shift > 7) {
        throw std::runtime_error("--align-shift must be 0 to 7");
    }
    // The aligned loads crash unless column 1 sits on a 32 byte boundary, 0 or 4 doubles past a 64 byte one.
    if (kernel_ != AlignKernel::loadu && align_shift % 4 != 0) {
        throw std::runtime_error(name() + " needs --align-shift 0 or 4");
    }

    // The payoff and the constants, the same as the omp solver.
    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;
    int total_cells = num_stock_nodes * num_variance_nodes;

    // The two sheets, placed at the chosen alignment, with the payoff copied in.
    ShiftedBuffer current_buffer;
    ShiftedBuffer next_buffer;
    make_shifted_buffer(current_buffer, total_cells, align_shift);
    make_shifted_buffer(next_buffer, total_cells, align_shift);
    for (int cell = 0; cell < total_cells; ++cell) {
        current_buffer.start[cell] = grid.current()[cell];
    }

    // The S and S squared tables, moved to the same alignment as the rows.
    ShiftedBuffer price_buffer;
    ShiftedBuffer price_sq_buffer;
    make_shifted_buffer(price_buffer, num_stock_nodes, align_shift);
    make_shifted_buffer(price_sq_buffer, num_stock_nodes, align_shift);
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        price_buffer.start[stock_i] = constants.stock_price[stock_i];
        price_sq_buffer.start[stock_i] = constants.stock_price_sq[stock_i];
    }
    const double* stock_price_table = price_buffer.start;
    const double* stock_price_sq_table = price_sq_buffer.start;

    // Say where column 1 sits for this run and where the shipped omp rung's Grid put it.
    std::fprintf(stderr, "[%s] column 1 at byte %d of its cache line, the shipped Grid has it at byte %d\n",
                 name().c_str(), byte_in_cache_line(current_buffer.start + 1), byte_in_cache_line(grid.current().data() + 1));

    // Pick the row function once, before the time loop.
    // A function pointer holds the address of a function, so row_update(...) calls whichever was picked.
    void (*row_update)(const StencilConstants&, const double*, double*, int, const double*, const double*) = update_row_avx;
    if (kernel_ == AlignKernel::load) {
        row_update = update_row_avx_load;
    }
    if (kernel_ == AlignKernel::shuffle) {
        row_update = update_row_avx_shuffle;
    }

    double* current = current_buffer.start;
    double* next = next_buffer.start;

    for (int step = 1; step <= num_timesteps; ++step) {
        double step_start = omp_get_wtime();

        // The same pragma as the omp solver.
        #pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes, row_update, stock_price_table, stock_price_sq_table)
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            row_update(constants, current, next, var_j, stock_price_table, stock_price_sq_table);
        }

        // Boundaries on the main thread, top row last.
        double boundary_start = omp_get_wtime();
        update_zero_row(constants, current, next);
        apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
        copy_top_row(constants, next);

        // Swap the two pointers. No numbers are copied.
        std::swap(current, next);
        double step_end = omp_get_wtime();
        result.boundary_seconds += step_end - boundary_start;
        result.seconds += step_end - step_start;
    }

    // Copy the answer into the grid and make it current, so extract_result can read it.
    std::vector<double>& answer = grid.next();
    for (int cell = 0; cell < total_cells; ++cell) {
        answer[cell] = current[cell];
    }
    grid.swap_buffers();
    extract_result(grid, cfg, result);

    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
