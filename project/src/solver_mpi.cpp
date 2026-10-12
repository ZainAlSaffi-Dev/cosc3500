// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <vector>

#include "io.h"
#include "solver.h"
#include "stencil.h"

// The block of variance rows one process owns, from begin up to end - 1.
// A plain record, like a Python tuple with names.
struct RowSlab {
    int begin = 0;
    int end = 0;
};

// Works out which rows a process owns when the rows are shared out evenly.
// When they do not divide evenly the first processes get one extra row each.
static RowSlab slab_for_rank(int num_rows, int rank, int num_ranks) {
    int rows_per_rank = num_rows / num_ranks;
    int extra_rows = num_rows % num_ranks;

    RowSlab slab;
    if (rank < extra_rows) {
        // This process is one of the first few, so it takes one extra row
        slab.begin = rank * (rows_per_rank + 1);
        slab.end = slab.begin + rows_per_rank + 1;
    } else {
        // Skip past the processes that took an extra row
        slab.begin = extra_rows * (rows_per_rank + 1) + (rank - extra_rows) * rows_per_rank;
        slab.end = slab.begin + rows_per_rank;
    }
    return slab;
}

// Swaps the edge rows of one sheet with the processes above and below.
// Afterwards this process holds a fresh copy of the row just outside each end of its slab.
static void exchange_halo_rows(std::vector<double>& sheet, int num_stock_nodes, RowSlab slab, int rank, int num_ranks) {
    // Each row is num_stock_nodes doubles in a row, so it goes as one message
    double* first_row = sheet.data() + slab.begin * num_stock_nodes;
    double* last_row = sheet.data() + (slab.end - 1) * num_stock_nodes;
    double* halo_below = sheet.data() + (slab.begin - 1) * num_stock_nodes;
    double* halo_above = sheet.data() + slab.end * num_stock_nodes;
    int halo_tag = 0;

    // One request per message, so we can wait for all of them at the end
    std::vector<MPI_Request> requests;

    // The bottom process has nobody below it
    if (rank > 0) {
        // Post the receive first so the message can land straight in place
        MPI_Request receive_below;
        MPI_Irecv(halo_below, num_stock_nodes, MPI_DOUBLE, rank - 1, halo_tag, MPI_COMM_WORLD, &receive_below);
        requests.push_back(receive_below);

        MPI_Request send_below;
        MPI_Isend(first_row, num_stock_nodes, MPI_DOUBLE, rank - 1, halo_tag, MPI_COMM_WORLD, &send_below);
        requests.push_back(send_below);
    }

    // The top process has nobody above it
    if (rank < num_ranks - 1) {
        MPI_Request receive_above;
        MPI_Irecv(halo_above, num_stock_nodes, MPI_DOUBLE, rank + 1, halo_tag, MPI_COMM_WORLD, &receive_above);
        requests.push_back(receive_above);

        MPI_Request send_above;
        MPI_Isend(last_row, num_stock_nodes, MPI_DOUBLE, rank + 1, halo_tag, MPI_COMM_WORLD, &send_above);
        requests.push_back(send_above);
    }

    // Nothing may touch these rows until every send and receive has finished
    MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

// Prices the option on several processes, each owning one slab of variance rows.
// Each process runs the OpenMP loop over its own rows and swaps one edge row with each neighbour per step.
SolveResult MpiSolver::solve(const Config& cfg) {
    SolveResult result;

    // Which process this is, and how many there are
    int rank = 0;
    int num_ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    result.ranks = num_ranks;
    result.threads = omp_get_max_threads();

    // Snapshots would only hold this process's rows, so they are not supported
    if (cfg.dump_every > 0) {
        throw std::runtime_error("--dump-every is not supported by the mpi solver");
    }

    // Every process keeps the whole grid but only updates its own rows
    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    // Only the first process prints the stability line
    if (rank == 0) {
        report_stability(grid, cfg, dt, result);
    }

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;

    // The top row copies the row below it, so the last slab needs at least two rows
    if (num_variance_nodes / num_ranks < 2) {
        throw std::runtime_error("too many processes, each one needs at least two variance rows");
    }
    // This process's rows, out of all rows from 0 to num_variance_nodes - 1
    RowSlab slab = slab_for_rank(num_variance_nodes, rank, num_ranks);

    // The interior rows inside this slab, skipping row 0 and the top row
    int first_interior = std::max(slab.begin, 1);
    int end_interior = std::min(slab.end, num_variance_nodes - 1);
    bool owns_zero_row = (rank == 0);
    bool owns_top_row = (rank == num_ranks - 1);

    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        auto step_start = std::chrono::steady_clock::now();

        // This process's interior rows, split between its threads as in the omp solver
        #pragma omp parallel for schedule(static) default(none) shared(constants, current, next, first_interior, end_interior)
        for (int var_j = first_interior; var_j < end_interior; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }

        // Each boundary is done only by the process that owns it
        auto boundary_start = std::chrono::steady_clock::now();
        if (owns_zero_row) {
            update_zero_row(constants, current, next);
        }
        apply_stock_edges(constants, next, slab.begin, slab.end, step * dt);
        if (owns_top_row) {
            copy_top_row(constants, next);
        }
        auto boundary_end = std::chrono::steady_clock::now();
        result.boundary_seconds += std::chrono::duration<double>(boundary_end - boundary_start).count();

        // Send our edge rows to the neighbours and receive theirs.
        // The next step reads one row past each end of the slab, so it needs these
        auto halo_start = std::chrono::steady_clock::now();
        exchange_halo_rows(next, num_stock_nodes, slab, rank, num_ranks);
        auto halo_end = std::chrono::steady_clock::now();
        result.halo_seconds += std::chrono::duration<double>(halo_end - halo_start).count();

        grid.swap_buffers();

        auto step_end = std::chrono::steady_clock::now();
        result.seconds += std::chrono::duration<double>(step_end - step_start).count();
    }

    // Collect every slab into the first process's spare sheet.
    // Gatherv is the version of gather where each process can send a different amount
    std::vector<int> counts(num_ranks);
    std::vector<int> offsets(num_ranks);
    for (int other_rank = 0; other_rank < num_ranks; ++other_rank) {
        RowSlab other_slab = slab_for_rank(num_variance_nodes, other_rank, num_ranks);
        counts[other_rank] = (other_slab.end - other_slab.begin) * num_stock_nodes;
        offsets[other_rank] = other_slab.begin * num_stock_nodes;
    }
    const double* my_rows = current.data() + slab.begin * num_stock_nodes;
    int my_count = (slab.end - slab.begin) * num_stock_nodes;
    MPI_Gatherv(my_rows, my_count, MPI_DOUBLE, next.data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // The gathered sheet becomes current, so extract_result reads it
    grid.swap_buffers();
    if (rank == 0) {
        extract_result(grid, cfg, result);
    }

    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
