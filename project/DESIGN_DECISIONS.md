# Design Decisions and Their Alternatives

This file lists every code and design decision in the project, what the
alternative would have looked like, why we went the way we did, and whether
the course taught it. It is for the interview question "why did you do it
this way and not that way". `HPC_JUSTIFICATIONS.md` has the measured numbers
behind the parallel choices, and this file points there instead of repeating
them.

Course coverage is checked against the transcripts in `docs/lecs/` (L04
SIMD, L05 OpenMP, L06 MPI, tutorials T04 OpenMP, T05 MPI, T06 CUDA) and the
course's starter code for the matrix assignment (`hpcfiles/matmul/`). The
starter code is the course's own style. It sets sizes with `#define`,
allocates with `malloc`, and reads `argv` with `atoi`. Wherever our code does
something different from that, it needs a reason, and the reason is here.

Every entry ends with a plain answer to whether the simpler version would
have done the job. Where it would, say so in the interview. Saying "the
simpler version works too, and here is the one thing ours adds" is a better
answer than defending something you do not need.

---

## Summary

| # | Decision | Main alternative | Taught in the course? | Simpler alternative good enough? |
|---|---|---|---|---|
| 1 | Settings read from a config file at run time | `#define` or `const` constants | No, starter code uses `#define` | No, the sweeps need it |
| 2 | Command-line flags override the file | Config file only, or flags only | No | Flags only would work |
| 3 | `--config` first, then flags read in pairs | Two passes so flags can come in any order | No | Ours is now the simple one |
| 4 | Unknown config keys and flags stop the program | Ignore them | No | No, typos would go unnoticed |
| 5 | Settings grouped in plain structs | Loose global variables | No | Globals would work, less clear |
| 6 | Grid as a class with private fields | Plain struct or raw arrays in solve | No | Yes, a struct would work |
| 7 | One flat array per sheet | 2-D array or vector of vectors | Partly, flat arrays in starter code | No, flat is already the simple one |
| 8 | `std::vector` for memory | `malloc`/`free` or `new[]`/`delete[]` | No, starter code uses `malloc` | No, vector is shorter and safer |
| 9 | Two sheets swapped each step | Copy the new sheet back each step | No | No, copying is slower and no simpler |
| 10 | Solver base class with one subclass per rung | One plain function per rung | No | Yes, plain functions would work |
| 11 | `make_solver` factory returning `unique_ptr` | `if` chain calling plain functions | No | Yes, goes with 10 |
| 12 | Errors thrown and caught once in `run` | Print and `exit(1)` on the spot | No | Mostly, but MPI must still finish |
| 13 | Each rung has its own copy of the time loop | One shared loop with a row-function parameter | No | Shared loop is shorter but harder to read |
| 14 | Constants worked out once in `StencilConstants` | Compute everything per cell like the baseline | No | No, it is the main serial speedup |
| 15 | Row pointers inside the kernel | `grid.index()` for every neighbour | Partly, pointers in L04 | No, this is part of the speed |
| 16 | Code split into headers and `.cpp` files | One big file | No | One file would work, harder to navigate |
| 17 | Include guards (`#ifndef`) | `#pragma once` | No | Either is fine |
| 18 | `static` helper functions | Plain functions | No | Yes, `static` only hides them |
| 19 | One CSV line on stdout, messages on stderr | Human-readable printout | No | No, the scripts parse the line |
| 20 | `printf` for output | `std::cout` | No | Either is fine |
| 21 | `omp_get_wtime()` for timing | `std::chrono::steady_clock` or `MPI_Wtime` | No | Ours is the simple one |
| 22 | Time accumulated inside each step | One timer around the whole loop | No | One timer works if snapshots are off |
| 23 | One binary for every rung, always inside MPI | A separate binary per rung | No | Separate binaries are more build work |
| 24 | Hand-written Makefile | CMake | No | Makefile is the simpler one |
| 25 | Explicit finite-difference scheme | Implicit scheme or Monte Carlo | No | Explicit is the simplest |
| 26 | `double` throughout | `float` | No | No, float was measured and rejected |
| 27 | AVX2 intrinsics written by hand | `#pragma omp simd` or `-O3` auto-vectorisation | Yes, intrinsics in L04 | `-O3` alone got half the gain |
| 28 | AVX2 flags on one file, scalar fallback with `#if` | `-march=native` on everything | Flags yes, the split no | Needed so the Mac build works |
| 29 | Unaligned loads (`loadu`) | Aligned memory plus `load` | Yes, L04 says no penalty now | No, aligned needs more code |
| 30 | Leftover cells done by the scalar function | Pad rows to a multiple of four | No | No, padding changes the grid |
| 31 | `parallel for` over variance rows | Over stock columns, or `collapse(2)` | Yes, L05 and T04 | No |
| 32 | `schedule(static)` | `dynamic` or `guided` | Yes, T04 | Static is the simple one |
| 33 | `default(none)` with an explicit `shared` list | Leave the default (everything shared) | Yes, L05 recommends it | Default works, less safe |
| 34 | Boundaries done serially after the loop | Boundaries inside the parallel loop | No | Yes, the serial version is simpler |
| 35 | New thread team each step | One persistent parallel region | Partly, `single` in L05 | Yes, gain was 3.5% |
| 36 | MPI split by slabs of rows | Blocks in both directions, or columns | No | Rows are the simplest |
| 37 | Every rank holds the whole grid | Each rank holds only its slab and two halo rows | No | No, whole grid is the simpler code |
| 38 | `MPI_Isend`/`MPI_Irecv` with `MPI_Waitall` | `MPI_Sendrecv`, or blocking send and receive | Non-blocking yes, `Waitall` and `Sendrecv` no | `MPI_Sendrecv` would be shorter |
| 39 | `MPI_Gatherv` at the end | `MPI_Gather` | Gather yes, Gatherv no | No, slabs can be uneven |
| 40 | Only rank 0 prints | Every rank prints | No | No |
| 41 | Price and Greeks read off the final grid with bilinear blending | Require the quote point on a node | No | Node-only would work for our configs |
| 42 | Greeks from neighbouring grid values | Re-solve with a bumped spot or variance | No | No, bumping costs extra solves |
| 43 | Tests as small programs with hand-written checks | A test framework such as GoogleTest | No | Ours is the simple one |

---

## Part 1. Settings and the command line

### 1. Why a config file read at run time, not `#define` constants?

The program reads every setting from a text file such as `config/bench.cfg`
when it starts.

```cpp
// config/bench.cfg
grid.ns         = 2048
grid.nv         = 512
grid.nt         = 2000

// src/params.cpp, load_config, one key at a time
} else if (key == "grid.ns") {
    cfg.grid.num_stock_nodes = std::stoi(value);
}
```

The course's starter code does it with `#define`, which would look like this.

```cpp
#define NUM_STOCK_NODES 2048
#define NUM_VARIANCE_NODES 512
#define NUM_TIMESTEPS 2000
```

A modern C++ version of the same idea is a constant, `const int
kNumStockNodes = 2048;`. Both are fixed when the program is compiled.

We read the settings at run time because the benchmarks change them on every
run. The slurm scripts call the same binary with `--ns`, `--nv`, `--nt`,
`--maturity` and `--solver` for the size sweep, the weak scaling sweep and
every rung of the ladder. With `#define`, each grid size would need its own
compile, so one sweep would build dozens of binaries, and a number in the
results could not be traced to one build. The config files also keep named
parameter sets (`reference.cfg` for the real price, `bench.cfg` for timing,
`smoke.cfg` for a quick check), each with comments explaining its numbers.

There is one real cost. With `#define` the compiler knows the grid size and
can sometimes make the loop a little faster. We did not measure that.

The course does not teach config files. The starter code uses `#define`.

Would the simpler version have done the job? No. Fixed constants would have
needed a rebuild for every point of every sweep.

### 2. Why do command-line flags override the config file?

```cpp
// src/main.cpp, the flags after --config <file>
} else if (flag == "--ns") {
    cfg.grid.num_stock_nodes = std::stoi(value);
}
```

The alternatives are a config file only (one file per grid size) or flags
only (every setting typed on the command line, with defaults in
`params.h`). Flags only would look like the starter code.

```cpp
int num_stock_nodes = 2048;
if (argc > 1) {
    num_stock_nodes = std::atoi(argv[1]);
}
```

We kept both because the file holds the 17 settings that describe the
problem and the flags change the two or three things a sweep varies. A file
for every sweep point would be hundreds of near-identical files, and flags
only would make every command 17 settings long.

The course does not cover this.

Would the simpler version have done the job? Flags only with defaults in
`params.h` would work. You would lose the commented, named parameter sets.

### 3. Why must `--config` come first, with the flags read in pairs?

```cpp
// The config file comes first, as in heston --config file.cfg.
if (argc < 3 || std::string(argv[1]) != "--config") {
    usage();
    return 1;
}
Config cfg = load_config(argv[2]);

// Every other flag is a pair, the flag then its value.
for (int i = 3; i < argc; i = i + 2) {
    std::string flag = argv[i];
    std::string value = argv[i + 1];
    ...
}
```

The earlier version went through the command line twice. The first pass
searched for `--config` anywhere, and the second applied the other flags
through a helper, `flag_value`, that took the loop counter by reference
(`int& i`) so it could move it on. That let flags come in any order.

We changed it to the simple version because every script already put
`--config` first, and the pair loop has nothing in it that needs explaining
beyond `i = i + 2`. A library such as `getopt` is the other alternative. It
handles flags in any order and flags without values, which this program
never has.

Would the simpler version have done the job? Yes, and it is now the one in
the code.

### 4. Why stop on an unknown config key or flag?

```cpp
// src/params.cpp, the end of apply
} else {
    // Nothing matched, so the key is a typo.
    throw std::runtime_error(where + ": unknown key '" + key + "'");
}
```

The alternative is to ignore lines it does not understand.

We stop because a typo such as `grid.nss = 4096` would otherwise run the
default grid size without telling you, and a benchmark would be labelled
with the wrong size. A value that is not a number at all, such as `abc`,
also stops the program with the file and line number.

```cpp
// src/params.cpp, load_config
try {
    apply(cfg, key, value, where);
} catch (const std::invalid_argument&) {
    throw std::runtime_error(where + ": not a number '" + value + "'");
}
```

The earlier version also refused text such as `1.5abc`, through a helper
that checked `std::stod` had used every character. We dropped that helper
to keep the code short. `std::stod` now reads `1.5abc` as 1.5, the same as
C's `atof` in the starter code.

Would the simpler version have done the job? It would run, but a mistake in
a config could produce wrong results that look right.

### 5. Why group the settings in plain structs?

```cpp
// include/params.h
struct OptionSpec {
    double strike = 5250.0;
    double maturity_years = 0.25;
    bool is_call = true;
};
struct Config {
    OptionSpec option;
    MarketSpec market;
    HestonParams heston;
    GridSpec grid;
    ...
};
```

The alternative is loose global variables, as in the starter code.

```cpp
double strike = 5250.0;
double maturity_years = 0.25;
int num_stock_nodes = 2048;
```

A struct is a named bundle of variables, like a Python dataclass. One
`Config` is passed into every solver, so a function's arguments say exactly
what it depends on, and the tests can build a `Config`, change one field and
run two solvers on it side by side. With globals, the tests would have to
change shared state between runs, and every function could read anything.

A `struct` and a `class` are the same thing in C++, except that a struct's
fields are public by default. We use `struct` for plain data with no rules
to protect.

Would the simpler version have done the job? Globals would work for a single
run. The tests that run several configurations in one program would be
messier.

---

## Part 2. Memory and the grid

### 6. Why is `Grid` a class with private fields and getters?

```cpp
class Grid {
public:
    int num_stock_nodes() const { return num_stock_nodes_; }
    const std::vector<double>& current() const { return current_; }
    std::vector<double>& next() { return next_; }
    void swap_buffers();
private:
    int num_stock_nodes_ = 0;
    std::vector<double> current_;
    std::vector<double> next_;
};
```

The alternative is a struct with public fields, or just two arrays and a
few `int`s declared at the top of each `solve`.

The class does one useful thing beyond organising the code. `current()`
returns a read-only reference, so writing into the sheet being read is a
compile error. That turns the double-buffer rule into something the compiler
checks. It also keeps the grid spacing and the index formula in one place.
The class shape also matches the OOP style you know from Python.

Would the simpler version have done the job? Yes. A struct with public
fields would compute the same answers. You would lose the compiler check on
the read-only sheet.

### 7. Why one flat array per sheet?

```cpp
// Position of node (stock_i, var_j) in the flat buffer.
int index(int stock_i, int var_j) const { return var_j * num_stock_nodes_ + stock_i; }
```

The alternatives are a vector of vectors, `std::vector<std::vector<double>>
grid(nv, std::vector<double>(ns));`, or a fixed 2-D array, `double
grid[512][2048];`.

A vector of vectors puts every row in its own separate block of memory, so
moving from one row to the next is a jump to an unrelated address. A fixed
2-D array is laid out flat like ours, but its size must be known when
compiling and an 8 MiB array is too big for the stack. One flat block is
also what MPI needs, since a whole row or slab can be sent as one pointer and
a count. The starter code's matrices are flat arrays as well.

`HPC_JUSTIFICATIONS.md` section 2 covers the layout and loop order.

Would the simpler version have done the job? The flat array is already the
simple choice.

### 8. Why `std::vector` and not `malloc` or `new[]`?

```cpp
current_.assign(total_cells, 0.0);
next_.assign(total_cells, 0.0);
```

The starter code allocates like this.

```cpp
double* current = (double*)std::malloc(total_cells * sizeof(double));
...
std::free(current);
```

and plain C++ without vector would use `new double[total_cells]` with a
matching `delete[]`.

A vector frees its own memory when it goes out of scope, so there is no
`free` to forget, no double free, and nothing leaks if an error is thrown
halfway through. It also knows its own size. `std::swap` on two vectors
swaps their internal pointers, so the sheet swap costs nothing. The project
rule is no raw `new` or `delete`.

There is one measured downside. `assign` fills the whole sheet with zeros on
the main thread, so on a two-socket node every page lands on socket 0's
memory (first touch). The peak bandwidth benchmark had to switch to `new
double[n]` touched in parallel to get full bandwidth. The solver fixes the
same problem with one MPI rank per socket instead. See
`HPC_JUSTIFICATIONS.md`, why page placement matters.

The course does not teach vector. The starter code uses `malloc`.

Would the simpler version have done the job? No. `malloc` needs the same
amount of code plus a matching `free` on every path out.

### 9. Why two sheets that swap, instead of copying back?

```cpp
void Grid::swap_buffers() {
    std::swap(current_, next_);
}
```

The alternative writes into a second sheet and then copies it back over the
first every step.

```cpp
for (int c = 0; c < total_cells; ++c) {
    current[c] = next[c];
}
```

Updating one sheet in place is not an option, because cells would read
neighbours that were already updated this step and the answer would be
wrong. Copying back is correct but moves 8 MiB per step for nothing. The
swap exchanges two pointers.

Would the simpler version have done the job? Copying back gives the same
answer, about as much code, and a slower program. The swap is the better
choice on both counts.

---

## Part 3. How the program is structured

### 10. Why a `Solver` base class with one subclass per rung?

```cpp
class Solver {
public:
    virtual ~Solver() = default;
    virtual SolveResult solve(const Config& cfg) = 0;
    virtual std::string name() const = 0;
protected:
    void extract_result(const Grid& grid, const Config& cfg, SolveResult& out) const;
    void report_stability(const Grid& grid, const Config& cfg, double dt, SolveResult& out) const;
};

class SerialSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "serial"; }
};
```

The simpler alternative is one plain function per rung, with the shared
helpers as plain functions too.

```cpp
SolveResult solve_baseline(const Config& cfg);
SolveResult solve_serial(const Config& cfg);
SolveResult solve_simd(const Config& cfg);
SolveResult solve_omp(const Config& cfg);
SolveResult solve_mpi(const Config& cfg);

void extract_result(const Grid& grid, const Config& cfg, SolveResult& out);
```

Both give the same speed and the same answers. The class version was chosen
for three reasons. Every rung is forced to provide the same `solve` and
`name`, so the tests and `main` treat them all alike. The shared helpers
can use `name()` to label their warnings. And the shape matches the Python
abstract base class pattern you already know, which the project style asked
for. `virtual ... = 0` is the C++ way of writing `@abstractmethod`.

The course does not cover classes.

Would the simpler version have done the job? Yes. Plain functions would work
equally well. If asked, say the class hierarchy is an organisational choice,
not a performance one, and that it maps onto the Python pattern of an
abstract base class with one subclass per variant.

### 11. Why a `make_solver` factory that returns a `unique_ptr`?

```cpp
std::unique_ptr<Solver> make_solver(const std::string& name) {
    if (name == "serial") {
        return std::make_unique<SerialSolver>();
    }
    ...
    throw std::runtime_error("unknown solver: " + name);
}

// main.cpp
std::unique_ptr<Solver> solver = make_solver(cfg.solver);
SolveResult result = solver->solve(cfg);
```

The solver is picked by name at run time from `--solver`, so the program
needs something that turns a string into a solver. `unique_ptr` is a pointer
that deletes what it points to when it goes out of scope, so there is no
`delete` to forget. Returning a plain `Solver` by value is not possible,
because the base class is abstract and copying a subclass into a base would
lose the subclass.

The alternative, if the rungs were plain functions, is a plain `if` chain.

```cpp
SolveResult run_solver(const std::string& name, const Config& cfg) {
    if (name == "serial") {
        return solve_serial(cfg);
    }
    ...
    throw std::runtime_error("unknown solver: " + name);
}
```

The course does not cover `unique_ptr`.

Would the simpler version have done the job? Yes, together with decision 10.
The factory exists because the rungs are classes.

### 12. Why throw errors and catch them in one place?

```cpp
// src/main.cpp
static int run(int argc, char** argv, int rank) {
    try {
        ...
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "heston: %s\n", e.what());
        return 1;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    ...
    int exit_code = run(argc, argv, rank);
    MPI_Finalize();
    return exit_code;
}
```

The alternative prints a message and exits on the spot.

```cpp
if (!file) {
    std::fprintf(stderr, "cannot open config: %s\n", path.c_str());
    std::exit(1);
}
```

With MPI running, `exit` skips `MPI_Finalize`, and `mpiexec` then reports
the job as crashed, sometimes leaving other ranks hanging. Throwing means
every error, from a missing file to a typo to too many processes, travels up
to the single `catch` in `run`, which prints it and returns, so
`MPI_Finalize` always runs. An exception is like a Python `raise`, and the
`catch` is like `except`.

The course does not cover exceptions.

Would the simpler version have done the job? For a program without MPI,
yes. With MPI, every error path would need `MPI_Abort` or careful cleanup,
which is more code than one `catch`.

### 13. Why does every rung have its own copy of the time loop?

`solver_serial.cpp`, `solver_simd.cpp` and `solver_omp.cpp` are nearly the
same 50 lines. The only line that changes between serial and simd is
`update_row_scalar` versus `update_row_avx`, and omp adds one pragma.

The alternative is one shared loop that takes the row function as a
parameter.

```cpp
using RowFunction = void (*)(const StencilConstants&, const std::vector<double>&, std::vector<double>&, int);

SolveResult time_loop(const Config& cfg, RowFunction update_row) {
    ...
    for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
        update_row(constants, current, next, var_j);
    }
    ...
}
```

The copies are deliberate. Each rung reads from top to bottom in one file,
so in the video and the interview you can open one file and show the whole
rung, and the difference between two rungs is a visible one-line change. A
shared loop would also not cover the OpenMP pragma or the MPI halo exchange,
which change the loop itself, not just the row function.

Would the simpler version have done the job? A shared loop would be fewer
lines, but each rung would be harder to show on its own. For explaining the
code, the copies are simpler.

### 14. Why work out the constants once in `StencilConstants`?

```cpp
// src/stencil.cpp, once before the time loop
inv_two_ds = 1.0 / (2.0 * stock_spacing);
stock_price_sq[stock_i] = price * price;
stock_diffusion[var_j] = 0.5 * variance;
```

The alternative is what the baseline does, which works everything out again
in every cell of every step.

```cpp
double dV_dS = (east - west) / (2.0 * grid.stock_spacing());
double updated = centre + dt * (0.5 * variance * stock_price * stock_price * d2V_dS2 + ...);
```

Division is several times slower than multiplication, and the baseline does
five per cell. Precomputing the reciprocals and tables is the main reason
the serial rung is 2.81 times the baseline. The baseline is kept as written
because it is the reference every other rung must match.

Would the simpler version have done the job? No, the baseline is the simpler
version and it is the slow one. See `HPC_JUSTIFICATIONS.md` section 3.

### 15. Why raw row pointers inside the kernel?

```cpp
const double* row_centre = current.data() + var_j * num_stock_nodes;
const double* row_above = row_centre + num_stock_nodes;
double east = row_centre[stock_i + 1];
```

The alternative is calling `grid.index(stock_i + 1, var_j)` for each of the
nine neighbours, as the baseline does.

The pointers are set once per row, so each neighbour is a fixed offset and
the compiler does no multiplications to find it. They do not own memory, the
vectors do, so they cannot leak. The SIMD lecture uses pointers into arrays
in the same way for loads.

Would the simpler version have done the job? `index()` gives the same
answer. It is part of what makes the baseline slow.

### 16 and 17. Why headers plus `.cpp` files, and include guards?

Each part of the program has a header in `include/` that says what exists
and a `.cpp` file in `src/` that does the work. Every header starts with an
include guard.

```cpp
#ifndef HESTON_GRID_H
#define HESTON_GRID_H
...
#endif
```

The alternative is one large `.cpp` file, which needs no headers.

Splitting lets the tests reuse the solver files without `main.cpp`, lets
only `stencil_avx.cpp` be compiled with the AVX2 flags (decision 28), and
keeps each file about one thing. The guard stops a header being pasted in
twice when two files include it. `#pragma once` at the top of a header does
the same in one line and every compiler we use supports it. The guard is the
older standard form.

Would the simpler version have done the job? One file would compile, but the
AVX2 flags would then apply to every rung, and the tests would need their
own copy of the solvers.

### 18. Why are some helper functions marked `static`?

```cpp
// src/params.cpp
static void apply(Config& cfg, const std::string& key, const std::string& value, const std::string& where) { ... }
```

`static` on a function outside a class means only that file can call it,
like a leading underscore in Python but enforced. Without it the program
still works. It stops two files clashing if both have a helper called
`apply`.

Would the simpler version have done the job? Yes.

---

## Part 4. Output and timing

### 19 and 20. Why one CSV line on stdout, messages on stderr, and `printf`?

```cpp
std::printf("%.10g,%.8g,%.8g,%.8g,%d,%d,%d,%s,%.6f,%.6g,%d,%d,%.6f,%.6f\n", ...);
std::fprintf(stderr, "[%s] feller %s ...\n", ...);
```

The alternative is a human-readable printout such as `Price: 196.17`.

Every benchmark is run by a script that appends the output to a CSV file and
plots it. One line of comma-separated numbers on stdout can be appended
directly, and anything meant for a person goes to stderr so it never mixes
into the data. New columns only go on the end so old CSV files still parse.

`printf` is used over `std::cout` because one format string controls the
number of digits for every column in one line. `std::cout` would need a
`std::setprecision` for each. The snapshot files are written the same way,
with C's `fopen`, `fprintf("%.10g")` and `fclose`, rather than a C++
`std::ofstream`, so all output in the project uses the one `printf` style.

Would the simpler version have done the job? A human printout would mean
parsing text in every script. `std::cout` would work as well as `printf`.

### 21 and 22. Why `omp_get_wtime` timing, summed step by step?

```cpp
double step_start = omp_get_wtime();
...
double step_end = omp_get_wtime();
result.seconds += step_end - step_start;

// Write snapshots after the timer stops.
```

The alternatives are C++'s `std::chrono::steady_clock` or `MPI_Wtime()`, and
timing the whole loop with one start and one stop. The chrono version looks
like this.

```cpp
auto step_start = std::chrono::steady_clock::now();
...
auto step_end = std::chrono::steady_clock::now();
result.seconds += std::chrono::duration<double>(step_end - step_start).count();
```

`omp_get_wtime()` returns the wall clock time in seconds as a plain
`double`, so the time for a step is one subtraction. The chrono version
needs `auto`, a time point type and a `.count()` conversion to get the same
number. Every file is built with OpenMP, so the baseline and serial solvers
call it too and every solver is timed the same way. On Linux both
`omp_get_wtime` (GCC's libgomp) and `steady_clock` read `CLOCK_MONOTONIC`, a
clock that never jumps if the system time is adjusted. The benchmark numbers
in `docs/RESULTS.md` were taken with the chrono version, and the switch reads
the same clock. Timing each step and adding up keeps the optional snapshot
file writing out of the measurement. The setup before the loop and the
result extraction after it are not timed either, because every rung does
them the same way.

The lectures time code but never name a timer. `omp_get_wtime` is defined in
the OpenMP API specification, runtime library routines, timing routines.

Would the simpler version have done the job? The timer is already the
simple one. One timer around the whole loop would be correct whenever snapshots are off, which they are for every benchmark. The
persistent-region experiment does use one timer, because its threads stay
inside the loop the whole time.

### 23. Why one binary for every rung, always started inside MPI?

```cpp
int main(int argc, char** argv) {
    // Every solver runs inside MPI. The non MPI solvers run as one process.
    MPI_Init(&argc, &argv);
```

The alternative is a separate program per rung, with only the MPI one
linked against MPI.

One binary means every rung is compiled with the same compiler, flags and
code, so a speed difference can only come from the rung itself. The tests
compare all the rungs in one program. Running a serial rung as one MPI
process costs nothing measurable.

Would the simpler version have done the job? Separate programs would mean
more Makefile targets and more chances for the builds to differ.

### 24. Why a hand-written Makefile?

The alternative is CMake, which generates the Makefile. The course's starter
code uses a plain Makefile, and ours follows it. The extras in ours are the
AVX2 flags for one file, Apple clang versus g++ for OpenMP, and the test and
rangpur targets.

Would the simpler version have done the job? The Makefile is already the
simpler one.

---

## Part 5. The numerical method

### 25. Why an explicit finite-difference scheme?

```cpp
out_row[stock_i] = centre + dt * change;
```

The alternatives are an implicit scheme such as Crank-Nicolson or ADI,
which solves a linear system every step, and Monte Carlo, which simulates
random paths.

The explicit update is one line per cell and every cell is independent,
which is the easiest possible thing to vectorise, thread and split with MPI.
The cost is that the timestep must stay small to be stable, so the reference
solve needs 980000 steps. An implicit scheme allows big steps but couples
every cell in a row through a tridiagonal solve, which is far harder to
parallelise. Monte Carlo parallelises easily but gives a noisy answer and
noisy Greeks. We use Monte Carlo and a Fourier price only as independent
checks (`make check-heston`).

See `HPC_JUSTIFICATIONS.md`, why an explicit scheme.

### 26. Why `double` and not `float`?

The `omp_float` experiment kept both sheets in `float`. It was 1.75 times
faster at 24 threads, and 4.46 times on the big grid, but the price moved by
$0.044 and gamma by 11.9 percent on the reference solve. That error is
bigger than the error from the grid itself, so `double` stayed.

---

## Part 6. SIMD

### 27. Why AVX2 intrinsics written by hand?

```cpp
__m256d centre = _mm256_loadu_pd(row_centre + stock_i);
__m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
...
change = _mm256_fmadd_pd(cross_factor, d2V_dS_dv, change);
```

The alternatives are adding `#pragma omp simd` above the scalar loop and
letting the compiler vectorise it, or compiling with `-O3`, which turns the
auto-vectoriser on.

```cpp
#pragma omp simd
for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {
    ...
}
```

The SIMD lecture teaches intrinsics, and the rung is meant to show the
technique by hand. Writing it out means we know exactly which instructions
run, and the gain can be credited to SIMD and not to the compiler. Measured,
`-O3` alone roughly doubles the baseline (57.7 to 112.5) and takes the
serial rung from 158.4 to 238.6, while the hand-written AVX2 rung reaches
379. Benchmarks stay at `-O2` so the ladder measures our changes.

The lecture covers intrinsics, `set1`, FMA and aligned versus unaligned
loads.

Would the simpler version have done the job? `-O3` would get part of the
way with no extra code, which is a fair point to make if asked. The hand
version is about 60 percent faster than that.

### 28. Why AVX2 flags on one file, with a scalar fallback?

```cpp
// src/stencil_avx.cpp
#if defined(__AVX2__) && defined(__FMA__)
... the AVX2 row ...
#else
void update_row_avx(...) {
    update_row_scalar(constants, current, next, var_j);
}
#endif
```

The alternative is `-march=native` or `-mavx2 -mfma` on every file.

Two reasons. The Mac has an ARM chip with no AVX2, so without a fallback the
project would not build there, and all development and testing happens on
the Mac. And if every file got the flags, the compiler could quietly use
AVX2 inside the serial rung, so the serial and SIMD rungs would no longer be
a clean comparison.

Would the simpler version have done the job? On rangpur alone, flags on
everything would compile. The Mac build and the clean serial baseline need
the split.

### 29. Why unaligned loads?

```cpp
__m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
```

The alternative is allocating with 32-byte alignment (`_mm_malloc` or
`std::aligned_alloc`) and using `_mm256_load_pd`.

The neighbour loads start one cell to the left and right of each other, so
at most one of the nine loads per group could ever be aligned, whatever the
allocation. The SIMD lecture says the speed penalty for unaligned loads has
been greatly reduced or removed on modern CPUs. Aligned loads were listed as
experiment E2 and not run.

Would the simpler version have done the job? `loadu` is the simpler one.

### 30. Why let the scalar function do the leftover cells?

```cpp
while (stock_i + 3 <= last_interior) {
    ... four cells ...
    stock_i = stock_i + 4;
}
// The scalar row does the last zero to three cells.
update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
```

The alternatives are padding every row to a multiple of four or using masked
loads (`_mm256_maskload_pd`).

Padding changes the grid's shape and the index formula everywhere. Masked
loads are more intrinsics to explain. Reusing the scalar function costs two
cells per row on the standard grid and is already tested.

---

## Part 7. OpenMP

### 31. Why `parallel for` over the variance rows?

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
```

The alternatives are splitting the stock columns inside each row, or
`collapse(2)` over both loops, which the matrix assignment used.

Rows give each thread long runs of memory in order, reuse the AVX2 row
function without changes, and make one thread's work per step large compared
with the cost of starting the threads. `collapse(2)` would cut rows into
pieces and break the AVX2 row function.

The lecture and tutorial cover `parallel for`.

### 32. Why `schedule(static)`?

The alternatives are `dynamic` and `guided`, which hand out rows while the
loop runs. Every row costs the same, so an even split decided up front is
ideal. Measured, `dynamic,1` ran at 0.60 of static. The tutorial covers
schedules.

### 33. Why `default(none)`?

The alternative is to leave it out, in which case every outside variable is
shared. `default(none)` makes the compiler refuse any variable that is not
listed, so nothing is shared by accident. The OpenMP lecture recommends it
in almost those words.

Would the simpler version have done the job? The code would give the same
answer without it. It is a safety check that costs one line.

### 34. Why are the boundaries done by one thread after the loop?

```cpp
// Only the main thread runs the boundaries.
update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
```

The alternative, tested as `omp_edges`, moves row 0 and the edges into the
parallel loop. It was 5 percent faster at 24 threads, because the serial
share had grown to 8.2 percent. We kept the serial version on the main rung
because it is simpler and the top row must wait for the finished rows
anyway.

Would the simpler version have done the job? The serial version is the
simpler one, at a 5 percent cost.

### 35. Why start a new thread team every step?

The alternative, tested as `omp_persistent`, opens one parallel region
around the whole time loop and uses `omp for` and `omp single` inside it. It
was 3.5 percent faster on the standard grid. The main rung keeps the one
pragma because it is the easiest to explain, and its correctness does not
depend on reasoning about barriers inside a region. The lecture covers
`single` and barriers.

---

## Part 8. MPI

### 36. Why split the grid into slabs of rows?

The alternatives are columns, or blocks in both directions.

Rows are contiguous in memory, so a halo is one row sent as one message with
no packing. Each rank's slab is processed by the same OpenMP row loop as the
omp rung. Blocks in both directions would need halos on four sides, with
column halos that are not contiguous.

### 37. Why does every rank keep the whole grid?

```cpp
// Every process keeps the whole grid and updates only its own rows.
Grid grid(cfg.grid, cfg.option);
```

The alternative is giving each rank only its own rows plus one halo row on
each side, with its own local row numbers.

Keeping the whole grid means `index()`, the row functions, the boundaries
and the final gather all use the same global row numbers as every other
rung, so no code had to change to handle local numbering. It costs memory,
16 MiB per rank on the standard grid, which is small next to a node's memory.

Would the simpler version have done the job? The whole grid is the simpler
code. The slab-only version would only matter for grids too big for one
node's memory.

### 38. Why `MPI_Isend`, `MPI_Irecv` and `MPI_Waitall`?

```cpp
MPI_Irecv(row_below, num_stock_nodes, MPI_DOUBLE, rank - 1, message_tag, MPI_COMM_WORLD, &receive_below);
MPI_Isend(first_row, num_stock_nodes, MPI_DOUBLE, rank - 1, message_tag, MPI_COMM_WORLD, &send_below);
...
MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
```

There are two alternatives. Blocking `MPI_Send` and `MPI_Recv`, as in the
tutorial, can deadlock if two neighbours both send first and wait for the
other to receive, so they need an ordering rule such as even ranks send
first. `MPI_Sendrecv` does a send and a receive together without deadlock.

```cpp
MPI_Sendrecv(first_row, num_stock_nodes, MPI_DOUBLE, below, 0,
             row_below, num_stock_nodes, MPI_DOUBLE, below, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Sendrecv(last_row, num_stock_nodes, MPI_DOUBLE, above, 0,
             row_above, num_stock_nodes, MPI_DOUBLE, above, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
```

With `MPI_PROC_NULL` as the neighbour of the end ranks, those two calls
would replace the whole of `exchange_edge_rows`.

We used non-blocking calls because the MPI lecture teaches them as the way
to let every process send and receive at once without deadlock. Posting the
receive first lets the message land straight in place. `MPI_Waitall` waits
for all four messages at the end.

The lecture covers blocking and non-blocking messages. `MPI_Waitall` and
`MPI_Sendrecv` are not named in the transcripts.

Would the simpler version have done the job? `MPI_Sendrecv` would be shorter
and just as fast here, since we wait straight away anyway. It is a fair
thing to say in the interview. The non-blocking version is what the lecture
showed, and it is the starting point for overlapping messages with
computation, which was planned as E9 and not done.

### 39. Why `MPI_Gatherv` and not `MPI_Gather`?

```cpp
MPI_Gatherv(my_values, my_count, MPI_DOUBLE, next.data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
```

`MPI_Gather`, covered in the tutorial, needs every rank to send the same
amount. When the rows do not divide evenly, for example 512 rows on 3 ranks,
the slabs are 171, 171 and 170 rows. `MPI_Gatherv` takes a count and a
starting position for each rank.

The alternative that skips the gather entirely is to have only the rank that
owns today's stock price and variance send the few values `extract_result`
needs. That moves less data but is more code, and the gather runs once per
solve, so it does not matter for speed.

Would the simpler version have done the job? Plain `Gather` would work only
when the rows divide evenly. The tests use 3 ranks on 32 rows on purpose to
make sure uneven slabs work.

### 40. Why does only rank 0 print?

Every rank runs the same program. Without the check, a run on 4 ranks would
print 4 CSV lines and the scripts would record the result 4 times. Only rank
0 has the gathered grid, so only it can extract the price anyway.

---

## Part 9. Reading the answer and testing it

### 41. Why bilinear blending to read the price?

```cpp
out.price = bilinear_interpolate(grid, values, base_i, base_j, stock_weight, variance_weight);
```

The alternative is to require that today's stock price and variance land
exactly on grid nodes and read the value directly.

The test grids and `demo.cfg` are set up so they do land on nodes, and then
the blend weights are zero and the blend returns the node value. Blending lets any
grid size be run, which the size sweep needs. A warning goes to stderr
whenever the point is between nodes.

Would the simpler version have done the job? For the reference and test
configs, yes. The sweeps over grid sizes would then be limited to sizes
where the spot lands on a node.

### 42. Why take the Greeks from neighbouring grid values?

```cpp
out.delta = (east - west) / (2.0 * grid.stock_spacing());
```

The alternative is bump and revalue, solving again with the spot moved
slightly and taking the difference in price. The grid already holds the
value at every stock price and variance, so the Greeks come free from
neighbouring cells. Bumping would need four extra full solves.

### 43. Why tests as small programs and not a test framework?

```cpp
// tests/test_parity.cpp
if (abs_gap < abs_tolerance) {
    std::printf("test_parity: PASS\n");
    return 0;
}
```

The alternative is a framework such as GoogleTest or Catch2. Each test here
is a short `main` that runs solvers, compares numbers and returns 0 or 1,
which `make test` checks. That needs nothing installed on rangpur.

Would the simpler version have done the job? Ours is the simple version.

---

## If you want the code simpler

None of the choices above changes the answers or the timings by more than a
few percent, apart from the ones that are the actual optimisations (14, 15,
27, 31).

The settings code (decisions 3 and 4) has already been simplified. It now
reads the config file as three words per line with `std::istringstream`
and plain `std::stod`, and reads the flags in pairs. The helpers `trim`,
`parse_double`, `parse_int` and `flag_value` are gone, and the price,
Greeks and tests are unchanged.

Two more swaps would remove unfamiliar ideas without touching performance.

- Decisions 10 and 11. Replace the `Solver` classes and the factory with one
  plain function per rung and an `if` chain. This removes `virtual`,
  `override` and `protected`. It is the largest piece of C++ that is about
  organisation, not speed.
- Decision 38. Replace `exchange_edge_rows` with two `MPI_Sendrecv` calls.

Both live in code that runs during the timed loop or shapes every solver,
so they would need the tests, the guides and a rangpur re-run.

---

## Documentation for the techniques the course did not teach

The HPC techniques below are used in the code but are not in the lectures
or tutorials. Each one is defined in official documentation, which is what
to point to if asked where it comes from. The C++ standard library calls
are all described on cppreference, https://en.cppreference.com.

| Technique | Where in the code | Documentation |
|---|---|---|
| `MPI_Isend`, `MPI_Irecv`, `MPI_Waitall` | `solver_mpi.cpp`, `exchange_edge_rows` | MPI Standard 4.0, chapter 3, nonblocking communication, https://www.mpi-forum.org/docs/ |
| `MPI_Gatherv` | `solver_mpi.cpp`, end of `solve` | MPI Standard 4.0, chapter 6, gather, https://www.mpi-forum.org/docs/ |
| `MPI_Sendrecv`, `MPI_PROC_NULL` (alternative only) | not used | MPI Standard 4.0, chapter 3, https://www.mpi-forum.org/docs/ |
| Halo (ghost) rows between neighbouring slabs | `solver_mpi.cpp` | W. Gropp, E. Lusk and A. Skjellum, Using MPI, 3rd edition, MIT Press, 2014, the Poisson example with ghost points |
| `mpiexec --map-by ppr:1:socket`, `--bind-to core` | `slurm/bench_m2_mpi.sh` | Open MPI 4.1 `mpirun` manual page, https://www.open-mpi.org/doc/v4.1/man1/mpirun.1.php |
| `schedule(runtime)` and `OMP_SCHEDULE` | `solver_experiments.cpp` | OpenMP API specification, worksharing-loop construct and environment variables, https://www.openmp.org/specifications/ |
| `OMP_PLACES=cores`, `OMP_PROC_BIND=close` | `slurm/bench_common.sh` | OpenMP API specification, environment variables, https://www.openmp.org/specifications/ |
| Implicit barriers at the end of `omp for` and `omp single` | `solver_experiments.cpp`, persistent solver | OpenMP API specification, worksharing constructs, https://www.openmp.org/specifications/ |
| `omp_get_wtime()` | every `solver_*.cpp`, the step timer | OpenMP API specification, runtime library routines, timing routines, https://www.openmp.org/specifications/ |
| First-touch page placement on NUMA nodes | explains decision 8 and the MPI layout | Linux kernel documentation, NUMA memory policy, https://www.kernel.org/doc/html/latest/admin-guide/mm/numa_memory_policy.html, and U. Drepper, What every programmer should know about memory, 2007, section 5 |
| `numactl --interleave=all` | experiment E7 | numactl manual page, https://man7.org/linux/man-pages/man8/numactl.8.html |
| Transparent huge pages | page-size results | Linux kernel documentation, transparent hugepage support, https://www.kernel.org/doc/html/latest/admin-guide/mm/transhuge.html |
| `MALLOC_MMAP_THRESHOLD_` | `slurm/bench_common.sh` | glibc `mallopt` manual page, `M_MMAP_THRESHOLD`, https://man7.org/linux/man-pages/man3/mallopt.3.html |
| Flush to zero and denormals are zero (`_mm_setcsr`) | `slurm/flush_denormals.c` | Intel 64 and IA-32 Architectures Software Developer's Manual, volume 1, the MXCSR register, and D. Goldberg, What every computer scientist should know about floating-point arithmetic, 1991 |
| `-mavx2 -mfma` and the `__AVX2__`, `__FMA__` macros | `Makefile`, `stencil_avx.cpp` | GCC manual, x86 options, https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html |
| `_mm256_fnmadd_pd`, `_mm256_loadu_pd` and the other intrinsics | `stencil_avx.cpp` | Intel Intrinsics Guide, https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html |
| Roofline model | results and slides | S. Williams, A. Waterman and D. Patterson, Roofline, Communications of the ACM 52(4), 2009 |
| Weak scaling | results and slides | J. L. Gustafson, Reevaluating Amdahl's law, Communications of the ACM 31(5), 1988 |
| Serial fraction fitted from measured speedups | results | A. H. Karp and H. P. Flatt, Measuring parallel processor performance, Communications of the ACM 33(5), 1990 |
| Triad bandwidth benchmark | `slurm/bench_peak.sh` | J. D. McCalpin, STREAM benchmark, https://www.cs.virginia.edu/stream/ |
| gprof flat profile | `slurm/profile_m2.sh` | GNU gprof manual, https://sourceware.org/binutils/docs/gprof/ |
| `sbatch`, exclusive nodes | `slurm/` | Slurm `sbatch` documentation, https://slurm.schedmd.com/sbatch.html |
