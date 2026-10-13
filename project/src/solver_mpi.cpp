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
#include "solver_experiments.h"
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

// Starts sending this process's first and last rows to the processes below and above, and starts receiving theirs.
// The received rows land just outside this process's rows. HPC texts call them halo or ghost rows.
// The calls return at once. The rows are only safe to touch after MPI_Waitall on the requests.
static void post_edge_rows(std::vector<double>& buffer, int num_stock_nodes, RowRange owned_rows, int rank, int num_ranks,
                           std::vector<MPI_Request>& requests) {
    // Each row is num_stock_nodes doubles next to each other, sent as one message.
    double* first_row = buffer.data() + owned_rows.begin * num_stock_nodes;
    double* last_row = buffer.data() + (owned_rows.end - 1) * num_stock_nodes;
    double* row_below = buffer.data() + (owned_rows.begin - 1) * num_stock_nodes;
    double* row_above = buffer.data() + owned_rows.end * num_stock_nodes;
    int message_tag = 0;

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
}

// Sends the edge rows and receives the neighbours' rows, and returns only when every message has finished.
static void exchange_edge_rows(std::vector<double>& buffer, int num_stock_nodes, RowRange owned_rows, int rank, int num_ranks) {
    // One request per message.
    // push_back adds to the end, like append in Python.
    std::vector<MPI_Request> requests;
    post_edge_rows(buffer, num_stock_nodes, owned_rows, rank, num_ranks, requests);

    // Wait for every send and receive to finish.
    MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

// Collects the rows of every process into the next buffer of the first process.
// The answer of the last step is in the current buffer of every process.
static void gather_rows(Grid& grid, int num_stock_nodes, int num_variance_nodes, int rank, int num_ranks) {
    // Gatherv lets each process send a different amount.
    std::vector<int> counts(num_ranks);
    std::vector<int> offsets(num_ranks);
    for (int other_rank = 0; other_rank < num_ranks; ++other_rank) {
        RowRange other_rows = rows_for_rank(num_variance_nodes, other_rank, num_ranks);
        counts[other_rank] = (other_rows.end - other_rows.begin) * num_stock_nodes;
        offsets[other_rank] = other_rows.begin * num_stock_nodes;
    }
    RowRange owned_rows = rows_for_rank(num_variance_nodes, rank, num_ranks);
    const double* my_values = grid.current().data() + owned_rows.begin * num_stock_nodes;
    int my_count = (owned_rows.end - owned_rows.begin) * num_stock_nodes;
    MPI_Gatherv(my_values, my_count, MPI_DOUBLE, grid.next().data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

    gather_rows(grid, num_stock_nodes, num_variance_nodes, rank, num_ranks);

    // Make the gathered buffer current for extract_result.
    grid.swap_buffers();
    if (rank == 0) {
        extract_result(grid, cfg, result);
    }

    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}

// Experiment E9. The mpi solver with the edge row messages overlapped with the interior update.
// Each step updates the two rows the neighbours need first, starts sending them,
// updates the remaining rows while the messages travel, and only then waits for them.
// The row updates and boundaries are the same as the mpi solver's, so the answer matches simd bit for bit.
SolveResult MpiOverlapSolver::solve(const Config& cfg) {
    SolveResult result;

    int rank = 0;
    int num_ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    result.ranks = num_ranks;
    result.threads = omp_get_max_threads();

    // One thread calls MPI from inside the parallel region.
    // MPI allows that only when it was started at thread level FUNNELED or higher.
    int thread_level = MPI_THREAD_SINGLE;
    MPI_Query_thread(&thread_level);
    if (thread_level < MPI_THREAD_FUNNELED) {
        throw std::runtime_error("the mpi_overlap solver needs MPI started with MPI_THREAD_FUNNELED");
    }

    if (cfg.dump_every > 0) {
        throw std::runtime_error("--dump-every is not supported by the mpi_overlap solver");
    }

    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    if (rank == 0) {
        report_stability(grid, cfg, dt, result);
    }

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;

    if (num_variance_nodes / num_ranks < 2) {
        throw std::runtime_error("too many processes, each one needs at least two variance rows");
    }
    RowRange owned_rows = rows_for_rank(num_variance_nodes, rank, num_ranks);

    int first_interior = std::max(owned_rows.begin, 1);
    int end_interior = std::min(owned_rows.end, num_variance_nodes - 1);
    int num_interior_rows = end_interior - first_interior;
    bool owns_zero_row = (rank == 0);
    bool owns_top_row = (rank == num_ranks - 1);

    // The rows the neighbours need, or -1 where there is no neighbour.
    // Every process owns at least two rows, so both are always interior rows.
    int send_below_row = -1;
    if (rank > 0) {
        send_below_row = owned_rows.begin;
    }
    int send_above_row = -1;
    if (rank < num_ranks - 1) {
        send_above_row = owned_rows.end - 1;
    }

    // The owned rows that are not sent. They get their stock edges after the parallel region.
    // A row may not be written while MPI is still sending it, so the sent rows are left out.
    int unsent_begin = owned_rows.begin;
    if (send_below_row >= 0) {
        unsent_begin = unsent_begin + 1;
    }
    int unsent_end = owned_rows.end;
    if (send_above_row >= 0) {
        unsent_end = unsent_end - 1;
    }

    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();
    std::vector<MPI_Request> requests;

    for (int step = 1; step <= num_timesteps; ++step) {
        double step_start = omp_get_wtime();
        double time_to_expiry = step * dt;
        requests.clear();

        #pragma omp parallel default(none) shared(constants, current, next, first_interior, num_interior_rows, send_below_row, send_above_row, time_to_expiry, num_stock_nodes, owned_rows, rank, num_ranks, requests, result)
        {
            // The same split of rows between threads that schedule(static) makes.
            // It is worked out by hand so each thread knows its rows before it starts.
            RowRange my_rows = rows_for_rank(num_interior_rows, omp_get_thread_num(), omp_get_num_threads());
            my_rows.begin = my_rows.begin + first_interior;
            my_rows.end = my_rows.end + first_interior;

            // Phase 1. A thread holding a row the neighbours need updates that row first.
            // Every other thread updates its own first row, so no thread waits at the barrier with nothing done.
            int early_rows[2] = {-1, -1};
            int num_early_rows = 0;
            if (my_rows.begin <= send_below_row && send_below_row < my_rows.end) {
                early_rows[num_early_rows] = send_below_row;
                num_early_rows = num_early_rows + 1;
            }
            if (my_rows.begin <= send_above_row && send_above_row < my_rows.end) {
                early_rows[num_early_rows] = send_above_row;
                num_early_rows = num_early_rows + 1;
            }
            if (num_early_rows == 0 && my_rows.begin < my_rows.end) {
                early_rows[0] = my_rows.begin;
                num_early_rows = 1;
            }
            for (int early = 0; early < num_early_rows; ++early) {
                update_row_avx(constants, current, next, early_rows[early]);
            }

            // Both rows to send are finished once every thread reaches this line.
            #pragma omp barrier

            // The master thread starts the messages. The other threads do not wait for it,
            // because master, unlike single, has no barrier at its end.
            #pragma omp master
            {
                double post_start = omp_get_wtime();
                // A sent row needs its stock edges before it leaves.
                if (send_below_row >= 0) {
                    apply_stock_edges(constants, next, send_below_row, send_below_row + 1, time_to_expiry);
                }
                if (send_above_row >= 0) {
                    apply_stock_edges(constants, next, send_above_row, send_above_row + 1, time_to_expiry);
                }
                post_edge_rows(next, num_stock_nodes, owned_rows, rank, num_ranks, requests);
                result.exchange_seconds += omp_get_wtime() - post_start;
            }

            // Phase 2. The rest of this thread's rows, while the messages travel.
            for (int var_j = my_rows.begin; var_j < my_rows.end; ++var_j) {
                if (var_j == early_rows[0] || var_j == early_rows[1]) {
                    continue;
                }
                update_row_avx(constants, current, next, var_j);
            }
        }

        // The boundaries of the rows that were not sent, while the messages may still be travelling.
        double boundary_start = omp_get_wtime();
        if (owns_zero_row) {
            update_zero_row(constants, current, next);
        }
        apply_stock_edges(constants, next, unsent_begin, unsent_end, time_to_expiry);
        if (owns_top_row) {
            copy_top_row(constants, next);
        }
        double boundary_end = omp_get_wtime();
        result.boundary_seconds += boundary_end - boundary_start;

        // Only now wait for the messages. The time spent here is the part of the exchange the update did not hide.
        MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        double wait_end = omp_get_wtime();
        result.exchange_seconds += wait_end - boundary_end;

        grid.swap_buffers();

        double step_end = omp_get_wtime();
        result.seconds += step_end - step_start;
    }

    gather_rows(grid, num_stock_nodes, num_variance_nodes, rank, num_ranks);

    grid.swap_buffers();
    if (rank == 0) {
        extract_result(grid, cfg, result);
    }

    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
