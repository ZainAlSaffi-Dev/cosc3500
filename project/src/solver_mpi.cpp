// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "io.h"
#include "solver.h"
#include "stencil.h"

// The variance rows one process owns, from begin up to end - 1.
// A plain record, like a Python tuple with names.
struct RowRange {
    int begin = 0;
    int end = 0;
};

// Works out which rows a process owns when the rows are shared out evenly.
// The first few processes get one extra row each when they do not divide evenly.
static RowRange rows_for_rank(int num_rows, int rank, int num_ranks) {
    int rows_per_rank = num_rows / num_ranks;
    int extra_rows = num_rows % num_ranks;

    RowRange rows;
    if (rank < extra_rows) {
        // One of the first few processes, with one extra row.
        rows.begin = rank * (rows_per_rank + 1);
        rows.end = rows.begin + rows_per_rank + 1;
    } else {
        // Skip past the processes that took an extra row.
        rows.begin = extra_rows * (rows_per_rank + 1) + (rank - extra_rows) * rows_per_rank;
        rows.end = rows.begin + rows_per_rank;
    }
    return rows;
}

// Sends this process's first and last rows to the processes below and above, and receives theirs.
// The received rows land just outside this process's rows. HPC texts call them halo or ghost rows.
static void exchange_edge_rows(std::vector<double>& buffer, int num_stock_nodes, RowRange owned_rows, int rank, int num_ranks) {
    // Each row is num_stock_nodes doubles next to each other, sent as one message.
    double* first_row = buffer.data() + owned_rows.begin * num_stock_nodes;
    double* last_row = buffer.data() + (owned_rows.end - 1) * num_stock_nodes;
    double* row_below = buffer.data() + (owned_rows.begin - 1) * num_stock_nodes;
    double* row_above = buffer.data() + owned_rows.end * num_stock_nodes;
    int message_tag = 0;

    // One request per message.
    // push_back adds to the end, like append in Python.
    std::vector<MPI_Request> requests;

    // The bottom process has nobody below it.
    if (rank > 0) {
        // Post the receive first.
        MPI_Request receive_below;
        MPI_Irecv(row_below, num_stock_nodes, MPI_DOUBLE, rank - 1, message_tag, MPI_COMM_WORLD, &receive_below);
        requests.push_back(receive_below);

        MPI_Request send_below;
        MPI_Isend(first_row, num_stock_nodes, MPI_DOUBLE, rank - 1, message_tag, MPI_COMM_WORLD, &send_below);
        requests.push_back(send_below);
    }

    // The top process has nobody above it.
    if (rank < num_ranks - 1) {
        MPI_Request receive_above;
        MPI_Irecv(row_above, num_stock_nodes, MPI_DOUBLE, rank + 1, message_tag, MPI_COMM_WORLD, &receive_above);
        requests.push_back(receive_above);

        MPI_Request send_above;
        MPI_Isend(last_row, num_stock_nodes, MPI_DOUBLE, rank + 1, message_tag, MPI_COMM_WORLD, &send_above);
        requests.push_back(send_above);
    }

    // Wait for every send and receive to finish.
    MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

// Prices the option on several processes, each owning a block of variance rows.
// Each process runs the OpenMP loop over its rows and swaps edge rows with its neighbours every step.
SolveResult MpiSolver::solve(const Config& cfg) {
    SolveResult result;

    // Which process this is, and how many there are.
    int rank = 0;
    int num_ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    result.ranks = num_ranks;
    result.threads = omp_get_max_threads();

    // Snapshots are not supported with MPI.
    if (cfg.dump_every > 0) {
        throw std::runtime_error("--dump-every is not supported by the mpi solver");
    }

    // Every process keeps the whole grid and updates only its own rows.
    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    // Only the first process prints the stability line.
    if (rank == 0) {
        report_stability(grid, cfg, dt, result);
    }

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;

    // Each process needs at least two rows.
    if (num_variance_nodes / num_ranks < 2) {
        throw std::runtime_error("too many processes, each one needs at least two variance rows");
    }
    // The rows this process owns.
    RowRange owned_rows = rows_for_rank(num_variance_nodes, rank, num_ranks);

    // Interior rows this process owns, without row 0 and the top row.
    int first_interior = std::max(owned_rows.begin, 1);
    int end_interior = std::min(owned_rows.end, num_variance_nodes - 1);
    bool owns_zero_row = (rank == 0);
    bool owns_top_row = (rank == num_ranks - 1);

    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        double step_start = omp_get_wtime();

        // Split this process's interior rows between its threads.
        #pragma omp parallel for schedule(static) default(none) shared(constants, current, next, first_interior, end_interior)
        for (int var_j = first_interior; var_j < end_interior; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }

        // Only the owning process updates each boundary.
        double boundary_start = omp_get_wtime();
        if (owns_zero_row) {
            update_zero_row(constants, current, next);
        }
        apply_stock_edges(constants, next, owned_rows.begin, owned_rows.end, step * dt);
        if (owns_top_row) {
            copy_top_row(constants, next);
        }
        double boundary_end = omp_get_wtime();
        result.boundary_seconds += boundary_end - boundary_start;

        // Send the edge rows to the neighbours and receive theirs.
        // The next step reads one row past each end of this process's rows.
        double exchange_start = omp_get_wtime();
        exchange_edge_rows(next, num_stock_nodes, owned_rows, rank, num_ranks);
        double exchange_end = omp_get_wtime();
        result.exchange_seconds += exchange_end - exchange_start;

        grid.swap_buffers();

        double step_end = omp_get_wtime();
        result.seconds += step_end - step_start;
    }

    // Collect the rows of every process into the next buffer of the first process.
    // Gatherv lets each process send a different amount.
    std::vector<int> counts(num_ranks);
    std::vector<int> offsets(num_ranks);
    for (int other_rank = 0; other_rank < num_ranks; ++other_rank) {
        RowRange other_rows = rows_for_rank(num_variance_nodes, other_rank, num_ranks);
        counts[other_rank] = (other_rows.end - other_rows.begin) * num_stock_nodes;
        offsets[other_rank] = other_rows.begin * num_stock_nodes;
    }
    const double* my_values = current.data() + owned_rows.begin * num_stock_nodes;
    int my_count = (owned_rows.end - owned_rows.begin) * num_stock_nodes;
    MPI_Gatherv(my_values, my_count, MPI_DOUBLE, next.data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Make the gathered buffer current for extract_result.
    grid.swap_buffers();
    if (rank == 0) {
        extract_result(grid, cfg, result);
    }

    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
