# Code Guide. The Heston Pricer, Snippet by Snippet

This guide goes through the C++ one block at a time, in the order the program
runs. It explains what each piece of code does and why it is written that
way. The maths behind it, such as what the PDE means, how finite differences
work and why there are two buffers, is in [STUDY_GUIDE.md](STUDY_GUIDE.md).
This guide points there whenever a block depends on an idea from it.

**Interactive version.** The
[Heston Code Walkthrough](https://claude.ai/artifact/5ggg9HfwJiXwvRRPwHHSWh)
shows this guide one snippet at a time, with progress ticks, a mode that
hides each explanation until you have tried it yourself, and the practice
questions as flashcards. It is built from this file, so rebuild it when this
file changes.

**How the snippets work.** Every snippet is copied from the real file, mostly
with its comments removed so it stays short. The file and function are named
above it. Under each one, the first sentence says what the block does and the
second says why we do it that way.

**Reading order.** If you are short on time, read §1 for the map, then §7
(the baseline, where the whole algorithm is in one function), then §8 and §9
(the serial rung). Everything else is the setup and plumbing around them.

---

## 1. The map

### 1.1 Which file does what

Each `.h` file declares names, and its `.cpp` file holds the bodies.
`#include` only brings in the declarations. The Makefile compiles every
`.cpp` on its own, and the linker joins them into one program. A Python
`import` gives you both halves at once, which is why this feels unfamiliar.

```text
main.cpp
 ├─ params.h   Config, load_config()                      body in params.cpp
 ├─ io.h       print_result_csv(), dump_snapshot()        body in io.cpp
 └─ solver.h   Solver, BaselineSolver, SerialSolver,
               SimdSolver, make_solver()
                 └─ grid.h   Grid (two sheets, index, swap)  body in grid.cpp

solver.cpp           Solver::extract_result, Solver::report_stability, make_solver
solver_baseline.cpp  BaselineSolver::solve, the whole algorithm in one function
solver_serial.cpp    SerialSolver::solve, only the time loop
solver_simd.cpp      SimdSolver::solve, the same loop calling the AVX2 row
solver_omp.cpp       OmpSolver::solve, the SIMD loop with one OpenMP pragma
solver_mpi.cpp       MpiSolver::solve, the OpenMP loop over one slab per process
solver_experiments.cpp  the omp loop with one thing changed, for experiments E4, E5, E8
solver_float.cpp     OmpFloatSolver::solve, the omp loop on float sheets (E3)
                       └─ stencil.h  StencilConstants and the row functions
                                     scalar bodies in stencil.cpp
                                     AVX2 row in stencil_avx.cpp
black_scholes.cpp    closed-form formula, used only by the tests
```

This block is the include map, with each header on the left and the file
holding its function bodies on the right. It is split this way so the
Makefile can compile each `.cpp` on its own and only rebuild what changed.

### 1.2 What happens when you run it

This is what happens when you run
`./heston --config config/smoke.cfg --solver serial`.

1. `main` calls `load_config` to build a `Config`, then calls
   `make_solver("serial")`, which hands back a `unique_ptr<SerialSolver>`.
2. `solver->solve(cfg)` runs `SerialSolver::solve`. It builds the `Grid` and
   fills in the payoff, prints the stability line with `report_stability`,
   and builds `StencilConstants` once.
3. Each timestep calls `update_row_scalar` for every interior row, then
   `update_zero_row`, `apply_stock_edges` and `copy_top_row`, then
   `swap_buffers`.
4. `extract_result` reads the price and Greeks off the last sheet, and `main`
   prints them with `print_result_csv`.

### 1.3 One timestep at a glance

Both solvers do these five things, in this order, every step.

```text
  for each step:
    1. rows 1 .. nv-2      full 9-point update      reads current, writes next
    2. row 0               drift-only update        reads current, writes next
    3. columns 0 and ns-1  known edge values        writes next
    4. row nv-1            copy of row nv-2         reads next, so it must come after 1
    5. swap_buffers        next becomes current
```

This block is the order of work inside one timestep. The order matters only
for step 4, which reads values that step 1 has just written. STUDY_GUIDE §6.1
has a picture of which cells each step writes.

### 1.4 Baseline and stencil side by side

The baseline does every step inline in one function. The serial rung moves
each step into its own function in `stencil.cpp`, so the SIMD, OpenMP and MPI
rungs can reuse them. This table maps one to the other.

| Block in `solver_baseline.cpp` | Same job in `stencil.cpp` |
|---|---|
| spacings, S, S² and the v terms recomputed per cell | `StencilConstants` constructor, built once |
| the inner `stock_i` loop over interior cells | `update_row_scalar`, one call per row |
| the v = 0 loop | `update_zero_row` |
| the S = 0 and S = max edge loop | `apply_stock_edges` |
| the top row copy | `copy_top_row` |

The names follow the PDE. `centre` is the cell being updated, and `east`,
`west`, `north`, `south`, `north_east` and so on are its neighbours. East and
west are one stock node apart, and north and south are one variance row
apart. `dV_dS`, `d2V_dS2`, `dV_dv`, `d2V_dv2` and `d2V_dS_dv` are the five
derivatives. The four per-row coefficients are named after the PDE term they
multiply, so `stock_diffusion` is 0.5v, `cross_diffusion` is ρξv,
`variance_diffusion` is 0.5ξ²v and `variance_drift` is κ(θ − v).
`stock_drift` is r − q. STUDY_GUIDE §5 has the full table from PDE term to
variable.

Splitting the stencil into its own file costs nothing measurable. The call
into `stencil.cpp` happens once per row, which is about 510 calls per step
against about a million cell updates, and the inner loop compiles the same
either way.

---

## 2. The C++ you need

### 2.1 Two ideas Python never taught you

**Value versus reference semantics.** In Python, every variable is a
reference to an object somewhere on the heap. In C++, `Grid g;` **is** the
object. It lives right there on the stack, and when the scope ends, it is
destroyed. If you want reference behaviour, you ask for it with `Grid&`. The
result is that you always know where your data physically is, which is what
high-performance code needs.

**Ownership.** Python's garbage collector quietly works out when memory can
be freed. C++ makes ownership explicit. Every byte belongs to exactly one
owner, and when the owner is destroyed, the memory is freed at that exact
point. This is **RAII** (Resource Acquisition Is Initialization). The `Grid`
constructor allocates its two vectors, and the destructor, which runs
automatically at the end of the scope, frees them. Our safety rules are
these. There is no raw `new` or `delete` anywhere, `std::vector` owns the
sheets, `std::unique_ptr` owns the solver object, and raw pointers appear
only as short-lived local variables inside the hot loop.

Why can't we hold a `Solver` by value? A variable of type `Solver` must
really be *some concrete subclass*, and its size isn't known from the
interface. That is the same reason interface variables in Java and Python are
references under the hood. C++ just makes you say it, with
`std::unique_ptr<Solver>`.

### 2.2 C++ constructs and their Python equivalents

This table gives a Python equivalent for each C++ construct in the code. In
the code itself, comments explain each new construct in plain English first
and only add "similar to X in Python" where it helps.

| C++ | Python equivalent | Watch out |
|---|---|---|
| `class Grid { private: ... }` | class with `_attrs` | privacy is *enforced* by the compiler |
| `class Solver` + `virtual solve() = 0` | `ABC` + `@abstractmethod` | you can't create the base class itself |
| `BaselineSolver : public Solver` + `override` | subclassing | `override` catches typos at compile time |
| `struct Config` | `@dataclass` | just data, all public |
| `std::vector<double>` | `list`, but typed and contiguous | the main container, and it owns its memory |
| `const Config& cfg` parameter | passing an object | can't be None and can't be changed |
| `std::unique_ptr<Solver>` | garbage-collected ownership, made explicit | the only owner, deleted automatically at the end of its scope |
| `const` after a method | (nothing) | a promise that the method doesn't change the object |
| `static` on a function in a `.cpp` | a leading underscore, `_helper` | only this file can call it, and the compiler enforces it |

---

## 3. The entry point (`main.cpp`)

**3.1 Reading a flag's value** (`flag_value`)

```cpp
static std::string flag_value(int argc, char** argv, int& i) {
    if (i + 1 >= argc) {
        throw std::runtime_error(std::string(argv[i]) + " needs a value");
    }
    i = i + 1;
    return argv[i];
}
```

This block returns the word after a flag, such as `1024` in `--ns 1024`, and
moves `i` past it. `i` is passed by reference (`int&`) so the caller's loop
counter moves too, which stops the value from being read again as if it were
a flag.

**3.2 First pass, find the config file**

```cpp
std::string config_path;
for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--config") {
        config_path = flag_value(argc, argv, i);
    }
}
if (config_path.empty()) {
    usage();
    return 1;
}
Config cfg = load_config(config_path);
```

This block looks only for `--config` and loads that file. We read the config
first so that every other flag can override it, whatever order the flags are
typed in.

**3.3 Second pass, apply the other flags**

```cpp
for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--config") {
        flag_value(argc, argv, i);
    } else if (arg == "--solver") {
        cfg.solver = flag_value(argc, argv, i);
    } else if (arg == "--ns") {
        cfg.grid.num_stock_nodes = std::stoi(flag_value(argc, argv, i));
    ...
    } else {
        usage();
        throw std::runtime_error("unknown flag: " + arg);
    }
}
```

This block writes each flag's value over the matching `Config` field. An
unknown flag stops the program, because a typo such as `--nss` would
otherwise be silently ignored and you would benchmark the wrong grid.

**3.4 Make the solver and run it**

```cpp
std::unique_ptr<Solver> solver = make_solver(cfg.solver);

int reps = 1;
if (cfg.bench_reps > 0) {
    reps = cfg.bench_reps;
}
for (int rep = 0; rep < reps; ++rep) {
    SolveResult result = solver->solve(cfg);
    if (rank == 0) {
        print_result_csv(result, cfg);
    }
}
```

This block creates the chosen solver and runs it once, or `--bench R` times.
Repeating inside one process lets the benchmark scripts take a median, and
`run` never needs to know which solver it has, because every solver has the
same `solve` method (§6.2). Under `mpiexec` there are several copies of the
program, and only rank 0 prints, so the output is still one CSV line per
solve.

**3.5 Catch errors at the top**

```cpp
} catch (const std::exception& e) {
    std::fprintf(stderr, "heston: %s\n", e.what());
    return 1;
}
```

This block catches any error thrown while reading the config, picking the
solver or setting up the MPI slabs, prints it and returns 1. Catching it once
here keeps the rest of the code free of error checks, and the non-zero exit
code tells a Slurm script the run failed.

**3.6 Start and stop MPI** (`main`)

```cpp
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int exit_code = run(argc, argv, rank);
    MPI_Finalize();
    return exit_code;
}
```

This block is the whole of `main`. Everything from §3.1 to §3.5 lives in
`run`, and `main` only starts MPI, asks which process this is, runs the
program and shuts MPI down. We do it this way because `MPI_Finalize` has to
run on every exit path, including errors, and keeping the work in `run`
means there is one place it is called. There is one binary for every rung,
so the non-MPI rungs also start MPI. Run without `mpiexec`, MPI simply
starts as one process (a singleton) and `rank` is 0.

---

## 4. The settings (`params.h`, `params.cpp`)

**4.1 The settings structs** (`params.h`)

```cpp
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
    int dump_every = 0;
    std::string dump_dir = "results";
    int bench_reps = 0;
    std::string solver = "baseline";
};
```

This block defines the settings as plain structs with default values, one
struct per topic and `Config` holding them all. Structs work like Python
dataclasses here. Grouping the settings means every function takes one
`const Config&` instead of twenty separate arguments.

**4.2 Parsing a number safely** (`parse_double`)

```cpp
static double parse_double(const std::string& value, const std::string& where) {
    std::size_t used = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error(where + ": not a number: '" + value + "'");
    }
    if (used != value.size()) {
        throw std::runtime_error(where + ": trailing text in '" + value + "'");
    }
    return parsed;
}
```

This block turns text like `0.045` into a double and rejects anything else.
`std::stod` alone would read `1.5abc` as 1.5 and ignore the rest, so we check
that it used every character. The error message includes the file name and
line number so a bad config is easy to find.

**4.3 Matching a key to a field** (`apply`)

```cpp
if (key == "option.strike") {
    cfg.option.strike = parse_double(value, where);
} else if (key == "grid.ns") {
    cfg.grid.num_stock_nodes = parse_int(value, where);
...
} else {
    throw std::runtime_error(where + ": unknown key '" + key + "'");
}
```

This block stores one `key = value` line in the matching field. It is a plain
`if` / `else if` chain so it reads top to bottom, and an unknown key is an
error for the same reason as an unknown flag in §3.3.

**4.4 Reading the file** (`load_config`)

```cpp
std::ifstream file(path);
...
while (std::getline(file, line)) {
    line_number = line_number + 1;
    std::size_t hash_position = line.find('#');
    if (hash_position != std::string::npos) {
        line = line.substr(0, hash_position);
    }
    line = trim(line);
    if (line.empty()) {
        continue;
    }
    std::size_t equals_position = line.find('=');
    ...
    std::string key = trim(line.substr(0, equals_position));
    std::string value = trim(line.substr(equals_position + 1));
    apply(cfg, key, value, where);
}
```

This block reads the config one line at a time, drops comments and blank
lines, and splits each line at the `=`. The file closes itself when `file`
goes out of scope at the end of the function, which is RAII doing the job of
Python's `with open(...)`.

---

## 5. The grid (`grid.h`, `grid.cpp`)

**5.1 The class layout** (`grid.h`)

```cpp
class Grid {
public:
    Grid(const GridSpec& spec, const OptionSpec& opt);
    int num_stock_nodes() const { return num_stock_nodes_; }
    double stock_spacing() const { return stock_spacing_; }
    double stock_price(int stock_i) const { return stock_i * stock_spacing_; }
    double variance(int var_j) const { return var_j * variance_spacing_; }
    int index(int stock_i, int var_j) const { return var_j * num_stock_nodes_ + stock_i; }
    const std::vector<double>& current() const { return current_; }
    std::vector<double>& next() { return next_; }
    void swap_buffers();
    void init_payoff(const OptionSpec& opt);
    ...
private:
    int num_stock_nodes_ = 0;
    int num_variance_nodes_ = 0;
    double stock_spacing_ = 0.0;
    double variance_spacing_ = 0.0;
    std::vector<double> current_;
    std::vector<double> next_;
};
```

This block declares the grid, with private data and public getters. The data
is private so nothing outside the class can resize a sheet or change the
spacing by mistake. The small getters are written inside the class so the
compiler can inline them, which makes `grid.index(...)` as fast as writing
the multiplication by hand.

**5.2 The constructor** (`Grid::Grid`)

```cpp
num_stock_nodes_ = spec.num_stock_nodes;
num_variance_nodes_ = spec.num_variance_nodes;
stock_spacing_ = spec.stock_max_multiplier * opt.strike / (num_stock_nodes_ - 1);
variance_spacing_ = spec.variance_max / (num_variance_nodes_ - 1);

std::size_t total_cells = (std::size_t)num_stock_nodes_ * num_variance_nodes_;
current_.assign(total_cells, 0.0);
next_.assign(total_cells, 0.0);
```

This block works out the spacing and makes both sheets, filled with zeros.
The spacing divides by `N - 1` because N nodes have N − 1 gaps between them.
Each sheet is one flat block of doubles rather than a list of rows, which
STUDY_GUIDE §12 explains.

**5.3 The index formula**

```cpp
int index(int stock_i, int var_j) const { return var_j * num_stock_nodes_ + stock_i; }
```

This block turns a (column, row) pair into a position in the flat sheet. The
row number is multiplied because a whole row of `num_stock_nodes_` cells sits
before it, so cells next to each other in a row are next to each other in
memory. STUDY_GUIDE §12 has a picture.

**5.4 The swap** (`Grid::swap_buffers`)

```cpp
void Grid::swap_buffers() {
    std::swap(current_, next_);
}
```

This block exchanges which sheet is read and which is written. We use
`std::swap` because a `std::vector` is a small object holding a pointer to
its block of memory, so swapping two vectors swaps two pointers and copies
none of the million doubles. STUDY_GUIDE §11 explains why there are two
sheets at all.

```text
  before swap                            after swap
  current_ : [pointer] ---> block A      current_ : [pointer] ---> block B
  next_    : [pointer] ---> block B      next_    : [pointer] ---> block A
```

**5.5 The payoff** (`Grid::init_payoff`)

```cpp
for (int var_j = 0; var_j < num_variance_nodes_; ++var_j) {
    for (int stock_i = 0; stock_i < num_stock_nodes_; ++stock_i) {
        double price = stock_price(stock_i);
        double payoff = 0.0;
        if (opt.is_call) {
            payoff = std::max(price - opt.strike, 0.0);
        } else {
            payoff = std::max(opt.strike - price, 0.0);
        }
        current_[index(stock_i, var_j)] = payoff;
    }
}
```

This block fills the first sheet with what the option pays at expiry. It
writes into `current_`, not `next_`, because the time loop always reads
`current` first, and the payoff is the first sheet it reads. Every row is the
same because the payoff does not depend on variance.

**5.6 Finding today's market on the grid** (`locate`)

```cpp
double position = std::clamp(coord / spacing, 0.0, (double)(num_nodes - 1));
double nearest_node = std::round(position);
bool on_node = std::fabs(position - nearest_node) < kNodeSnapCells;
...
int highest_cell = std::max(1, num_nodes - 3);
out.cell = std::clamp((int)lower_node, 1, highest_cell);
if (on_node && out.cell == (int)lower_node) {
    out.weight = 0.0;
} else {
    out.weight = std::clamp(position - out.cell, 0.0, 1.0);
}
```

This block finds the cell just below today's spot (or v0), and how far across
that cell the point sits. The weight is snapped to exactly 0 when the point
is within 1e-9 of a node, so an aligned grid always returns the node's value
unchanged. The cell is kept at least one node in from each edge, because the
Greeks read one cell further out on each side (STUDY_GUIDE §10).

---

## 6. The solver family (`solver.h`, `solver.cpp`)

**6.1 What a solver returns** (`SolveResult`)

```cpp
struct SolveResult {
    double price = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double vega = 0.0;
    double seconds = 0.0;
    double cell_updates_per_sec = 0.0;
    double dt_stable_estimate = 0.0;
    int threads = 1;
    int ranks = 1;
    double halo_seconds = 0.0;
    double boundary_seconds = 0.0;
};
```

This block is the record every solver hands back. Returning one struct keeps
`main`, the CSV printer and the tests the same for every solver. `threads`
and `ranks` say how many OpenMP threads and MPI processes did the work. They
default to 1, so only the parallel rungs have to set them. `halo_seconds` is
the part of `seconds` that the MPI rung spent swapping halo rows, and
`boundary_seconds` the part the OpenMP and MPI rungs spent on the serial
boundary work. Together they let the benchmarks split a step into its
pieces.

**6.2 The base class** (`Solver`)

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
```

This block is the interface every solver follows, like a Python abstract base
class. `= 0` means each subclass must write its own `solve` and `name`, and
`virtual` means a call through a `Solver` pointer runs the subclass's
version. The two helpers are `protected`, so subclasses can call them but
`main` cannot, and they are written once so every rung reads results the
same way.

**6.3 The solvers**

```cpp
class BaselineSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "baseline"; }
};

class SerialSolver : public Solver {
public:
    SolveResult solve(const Config& cfg) override;
    std::string name() const override { return "serial"; }
};
```

`SimdSolver`, `OmpSolver` and `MpiSolver` follow the same pattern with the
names `simd`, `omp` and `mpi`.

This block declares the rungs. `override` makes the
compiler check that the method really matches one in the base class, so a
typo in the signature is a compile error instead of a silent new method.
Each new rung adds one class like these.

**6.4 Picking a solver by name** (`make_solver`)

```cpp
std::unique_ptr<Solver> make_solver(const std::string& name) {
    if (name == "baseline") {
        return std::make_unique<BaselineSolver>();
    }
    if (name == "serial") {
        return std::make_unique<SerialSolver>();
    }
    if (name == "simd") {
        return std::make_unique<SimdSolver>();
    }
    if (name == "omp") {
        return std::make_unique<OmpSolver>();
    }
    if (name == "mpi") {
        return std::make_unique<MpiSolver>();
    }
    throw std::runtime_error("unknown solver: " + name);
}
```

This block turns the `--solver` text into a solver object. It returns a
`unique_ptr` to the `Solver` base class, so `main` can call `solve` without
knowing which rung it has, and the object is deleted when `main` ends. Each
new rung adds one `if` here.

**6.5 The stability report** (`Solver::report_stability`)

```cpp
out.dt_stable_estimate = 1.0 / (max_variance * max_stock_price * max_stock_price / (stock_spacing * stock_spacing) + xi * xi * max_variance / (variance_spacing * variance_spacing) + rate);

bool feller_ok = 2.0 * kappa * theta >= xi * xi;
...
std::fprintf(stderr, "[%s] feller %s (2*kappa*theta=%.4f, xi^2=%.4f)  dt=%.3e  dt_stable~%.3e\n",
             name().c_str(), feller_text, 2.0 * kappa * theta, xi * xi, dt, out.dt_stable_estimate);
```

This block works out the largest stable timestep and prints it next to the
one we are using, along with whether the Feller condition holds. The formula
is the point where the centre weight of the update reaches zero at the far
corner of the grid (STUDY_GUIDE §7.5). It prints to stderr so stdout stays
one clean CSV line for the scripts.

**6.6 Reading the answer** (`Solver::extract_result`)

```cpp
AxisPosition spot_position = grid.stock_position(cfg.market.spot);
AxisPosition v0_position = grid.variance_position(cfg.heston.v0);
...
out.price = blend_cell(grid, values, base_i, base_j, stock_weight, variance_weight);
double east = blend_cell(grid, values, base_i + 1, base_j, stock_weight, variance_weight);
double west = blend_cell(grid, values, base_i - 1, base_j, stock_weight, variance_weight);
double north = blend_cell(grid, values, base_i, base_j + 1, stock_weight, variance_weight);
double south = blend_cell(grid, values, base_i, base_j - 1, stock_weight, variance_weight);

out.delta = (east - west) / (2.0 * grid.stock_spacing());
out.gamma = (east - 2.0 * out.price + west) / (grid.stock_spacing() * grid.stock_spacing());
out.vega = (north - south) / (2.0 * grid.variance_spacing());
```

This block reads the price at today's spot and variance, and gets the Greeks
from the same finite differences as STUDY_GUIDE §7. Every value goes through
the same blend, so the reported delta is the slope of the reported price and
not of some other surface.

**6.7 The blend** (`blend_cell`)

```cpp
double stock_lower = 1.0 - stock_weight;
double variance_lower = 1.0 - variance_weight;
double bottom_left = values[grid.index(base_i, base_j)];
double bottom_right = values[grid.index(base_i + 1, base_j)];
double top_left = values[grid.index(base_i, base_j + 1)];
double top_right = values[grid.index(base_i + 1, base_j + 1)];
return variance_lower * (stock_lower * bottom_left + stock_weight * bottom_right) + variance_weight * (stock_lower * top_left + stock_weight * top_right);
```

This block is bilinear interpolation. It blends along the stock axis on the
bottom and top rows, then blends those two results along the variance axis.
When both weights are 0 it reduces to `1.0 * bottom_left` plus zeros, which
returns the node value exactly, so aligned grids give bit-identical answers.

---

## 7. The baseline (`solver_baseline.cpp`, `BaselineSolver::solve`)

The baseline has the whole algorithm in one function, with nothing
precomputed. It is slow on purpose, and it is the fixed reference every other
rung is tested against.

**7.1 Set up the grid**

```cpp
Grid grid(cfg.grid, cfg.option);
grid.init_payoff(cfg.option);

double dt = cfg.option.maturity_years / num_timesteps;
report_stability(grid, cfg, dt, result);
```

This block makes the two sheets, fills the first one with the payoff and
works out the timestep. We print the stability report before any maths runs
because a `dt` above the limit makes the whole run blow up (STUDY_GUIDE §8),
and it is easier to spot in the first line of output than as a NaN at the
end.

**7.2 Name the two sheets**

```cpp
const std::vector<double>& current = grid.current();
std::vector<double>& next = grid.next();
```

This block gives short names to the sheet we read and the sheet we write.
`current` is `const` because nothing in a step is allowed to change it, and
the compiler enforces that. These names stay correct after every swap
because `current` always refers to the `current_` member and `next` to the
`next_` member. The swap changes what is inside those members, not which
member the names refer to.

**7.3 The time loop and the interior loops**

```cpp
for (int step = 1; step <= num_timesteps; ++step) {
    auto step_start = std::chrono::steady_clock::now();
    for (int var_j = 1; var_j < grid.num_variance_nodes() - 1; ++var_j) {
        for (int stock_i = 1; stock_i < grid.num_stock_nodes() - 1; ++stock_i) {
```

This block visits every interior cell once per step, with variance rows on
the outside and stock columns on the inside. The inner loop runs along
`stock_i` because cells in the same row sit next to each other in memory, so
the processor reads them in order. Both loops start at 1 and stop one short
of the end because the edge cells have no neighbour on one side.

**7.4 Read the nine cells**

```cpp
double stock_price = grid.stock_price(stock_i);
double variance = grid.variance(var_j);

double centre = current[grid.index(stock_i, var_j)];
double east = current[grid.index(stock_i + 1, var_j)];
double west = current[grid.index(stock_i - 1, var_j)];
double north = current[grid.index(stock_i, var_j + 1)];
double south = current[grid.index(stock_i, var_j - 1)];
double north_east = current[grid.index(stock_i + 1, var_j + 1)];
double north_west = current[grid.index(stock_i - 1, var_j + 1)];
double south_east = current[grid.index(stock_i + 1, var_j - 1)];
double south_west = current[grid.index(stock_i - 1, var_j - 1)];
```

This block loads the cell and its eight neighbours from the current sheet,
laid out as in the diagram in STUDY_GUIDE §7.2. Everything is read from
`current` and nothing from `next`, which is the rule that makes the update
correct (STUDY_GUIDE §11).

**7.5 The five derivatives**

```cpp
double dV_dS = (east - west) / (2.0 * grid.stock_spacing());
double d2V_dS2 = (east - 2.0 * centre + west) / (grid.stock_spacing() * grid.stock_spacing());
double dV_dv = (north - south) / (2.0 * grid.variance_spacing());
double d2V_dv2 = (north - 2.0 * centre + south) / (grid.variance_spacing() * grid.variance_spacing());
double d2V_dS_dv = (north_east - north_west - south_east + south_west) / (4.0 * grid.stock_spacing() * grid.variance_spacing());
```

This block turns the nine values into the five derivatives the PDE needs,
using the formulas in STUDY_GUIDE §7.2. The baseline divides by the spacing
every time, five divisions per cell per step. That is slow, and the serial
rung removes it (§8.2), but here it keeps each line identical to the maths.

**7.6 One step of the PDE**

```cpp
double updated = centre + dt * (0.5 * variance * stock_price * stock_price * d2V_dS2
                                + rho * xi * variance * stock_price * d2V_dS_dv
                                + 0.5 * xi * xi * variance * d2V_dv2
                                + (rate - div_yield) * stock_price * dV_dS
                                + kappa * (theta - variance) * dV_dv
                                - rate * centre);
next[grid.index(stock_i, var_j)] = updated;
```

This block is the update from STUDY_GUIDE §7.3, with one PDE term per line,
written into the next sheet. The baseline keeps it as one expression, rather
than the named terms the serial rung uses, so that its output stays
bit-identical to Milestone 1.

**7.7 The v = 0 row**

```cpp
for (int stock_i = 1; stock_i < grid.num_stock_nodes() - 1; ++stock_i) {
    double stock_price = grid.stock_price(stock_i);
    double centre = current[grid.index(stock_i, 0)];
    double dV_dS = (current[grid.index(stock_i + 1, 0)] - current[grid.index(stock_i - 1, 0)]) / (2.0 * grid.stock_spacing());
    double dV_dv = (current[grid.index(stock_i, 1)] - centre) / grid.variance_spacing();
    double updated = centre + dt * ((rate - div_yield) * stock_price * dV_dS + kappa * theta * dV_dv - rate * centre);
    next[grid.index(stock_i, 0)] = updated;
}
```

This block updates the bottom row with only the drift and discount terms. The
diffusion terms all have a factor of v, so at v = 0 they are zero, and
`dV_dv` looks one row up because there is no row below zero (STUDY_GUIDE §9).

**7.8 The stock edges**

```cpp
double time_to_expiry = step * dt;
double zero_stock_value = 0.0;
double max_stock_value = 0.0;
if (cfg.option.is_call) {
    zero_stock_value = 0.0;
    max_stock_value = max_stock_price * std::exp(-div_yield * time_to_expiry) - strike * std::exp(-rate * time_to_expiry);
} else {
    zero_stock_value = strike * std::exp(-rate * time_to_expiry);
    max_stock_value = 0.0;
}
for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
    next[grid.index(0, var_j)] = zero_stock_value;
    next[grid.index(grid.num_stock_nodes() - 1, var_j)] = max_stock_value;
}
```

This block writes the known option value into the first and last stock
column of every row. We can write them straight in because the option's value
at S = 0 and at a very large S is known without solving anything, and it only
depends on how long is left until expiry.

**7.9 The top row**

```cpp
int top_row = grid.num_variance_nodes() - 1;
for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
    next[grid.index(stock_i, top_row)] = next[grid.index(stock_i, top_row - 1)];
}
```

This block copies the second-top row into the top row. It reads from `next`,
not `current`, so it has to run after the interior update has finished. The
copy says the value stops changing with variance at the top of the grid,
which is the Neumann condition from STUDY_GUIDE §9.

**7.10 Swap, time and snapshot**

```cpp
grid.swap_buffers();

auto step_end = std::chrono::steady_clock::now();
result.seconds += std::chrono::duration<double>(step_end - step_start).count();

if (cfg.dump_every > 0 && step % cfg.dump_every == 0) {
    dump_snapshot(cfg.dump_dir, step, grid);
}
```

This block makes the new sheet the current one, adds this step's time to the
total and sometimes writes the sheet to a file. The timer stops before the
snapshot so that writing files never counts as compute time in a benchmark.

**7.11 Read off the answer**

```cpp
extract_result(grid, cfg, result);

double total_cell_updates = (double)grid.num_stock_nodes() * grid.num_variance_nodes() * num_timesteps;
result.cell_updates_per_sec = total_cell_updates / result.seconds;
```

This block reads the price and Greeks off the final sheet and works out the
throughput. We report cell updates per second rather than seconds because it
lets grids of different sizes be compared on the same chart. The `(double)`
stops the multiplication overflowing an `int` on big grids.

---

## 8. The shared stencil (`stencil.h`, `stencil.cpp`)

The serial rung does exactly what §7 does. The difference is where the work
happens. Anything that does not change during a solve is worked out once in
`StencilConstants`, and the per-row work lives in functions that the SIMD,
OpenMP and MPI rungs will reuse.

**8.1 The constants struct** (`stencil.h`)

```cpp
struct StencilConstants {
    int num_stock_nodes = 0;
    int num_variance_nodes = 0;
    double dt = 0.0;
    double rate = 0.0;
    double stock_drift = 0.0;
    double inv_two_ds = 0.0;
    double inv_ds_sq = 0.0;
    ...
    std::vector<double> stock_price;
    std::vector<double> stock_price_sq;
    std::vector<double> stock_diffusion;
    std::vector<double> cross_diffusion;
    std::vector<double> variance_diffusion;
    std::vector<double> variance_drift;
    StencilConstants(const Grid& grid, const Config& cfg, double timestep);
};
```

This block holds every number the update needs that stays the same for the
whole solve. It is a plain struct because it is just a record of values, and
building it once before the time loop means the loop only ever reads it.

**8.2 The reciprocals** (`StencilConstants` constructor)

```cpp
double stock_spacing = grid.stock_spacing();
double variance_spacing = grid.variance_spacing();
inv_two_ds = 1.0 / (2.0 * stock_spacing);
inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
inv_two_dv = 1.0 / (2.0 * variance_spacing);
inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
inv_four_ds_dv = 1.0 / (4.0 * stock_spacing * variance_spacing);
inv_dv = 1.0 / variance_spacing;
```

This block does each division once, before the time loop, and stores the
answer. We do it because a division costs about ten times as much as a
multiplication, and the derivatives can then multiply by these instead of
dividing. In Milestone 1 this one change was a 2.5× speedup
(STUDY_GUIDE §13).

**8.3 The tables**

```cpp
for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
    double price = stock_i * stock_spacing;
    stock_price[stock_i] = price;
    stock_price_sq[stock_i] = price * price;
}

for (int var_j = 0; var_j < num_variance_nodes; ++var_j) {
    double variance = var_j * variance_spacing;
    stock_diffusion[var_j] = 0.5 * variance;
    cross_diffusion[var_j] = (rho * xi) * variance;
    variance_diffusion[var_j] = ((0.5 * xi) * xi) * variance;
    variance_drift[var_j] = kappa * (theta - variance);
}
```

This block stores S and S² for every column, and the four PDE coefficients
for every row. We do it because these numbers depend only on where a cell is,
never on time, so the baseline was recomputing the same values in every one
of the thousands of timesteps.

**8.4 Copy the constants into locals** (`update_cells_scalar`)

```cpp
double dt = constants.dt;
double rate = constants.rate;
double stock_drift = constants.stock_drift;
double inv_two_ds = constants.inv_two_ds;
...
double stock_diffusion = constants.stock_diffusion[var_j];
double cross_diffusion = constants.cross_diffusion[var_j];
double variance_diffusion = constants.variance_diffusion[var_j];
double variance_drift = constants.variance_drift[var_j];
```

This block copies everything the row needs out of the struct and into local
variables. We do it because the loop writes through a pointer, and the
compiler cannot prove that pointer does not point into the struct. Without
the copies it would reload every constant from memory after every write.
Locals live in registers, so the loop never goes back to memory for them.

**8.5 Row pointers**

```cpp
const double* row_centre = current.data() + var_j * num_stock_nodes;
const double* row_above = row_centre + num_stock_nodes;
const double* row_below = row_centre - num_stock_nodes;
double* out_row = next.data() + var_j * num_stock_nodes;
```

This block points at the start of the three rows we read and the one row we
write. With these, `row_above[stock_i]` replaces
`current[grid.index(stock_i, var_j + 1)]`, so no index multiplication happens
inside the loop. The vectors still own the memory, and the pointers only look
into it.

```text
  current (one flat array)
  ... | row var_j-1 ............ | row var_j ............. | row var_j+1 ........... | ...
        ^ row_below                ^ row_centre              ^ row_above
          = row_centre - ns                                    = row_centre + ns
```

**8.6 Neighbours, derivatives and named terms**

```cpp
double centre = row_centre[stock_i];
double east = row_centre[stock_i + 1];
...
double dV_dS = (east - west) * inv_two_ds;
double d2V_dS2 = (east - 2.0 * centre + west) * inv_ds_sq;
...
double stock_diffusion_term = stock_diffusion * stock_price_sq[stock_i] * d2V_dS2;
double cross_term = cross_diffusion * stock_price[stock_i] * d2V_dS_dv;
double variance_diffusion_term = variance_diffusion * d2V_dv2;
double stock_drift_term = stock_drift * stock_price[stock_i] * dV_dS;
double variance_drift_term = variance_drift * dV_dv;

double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
out_row[stock_i] = centre + dt * change;
```

This block is §7.4 to §7.6 again, with multiplications by the reciprocals and
table lookups instead of recomputing. Each PDE term gets its own name so it
can be matched to the table in STUDY_GUIDE §5. The answer differs from the
baseline only in the last digits (about 3e-15 relative), because `x * (1/h)`
rounds slightly differently from `x / h`, and `test_ladder_matches` checks
this on every build (§12.3).

**8.7 The v = 0 row** (`update_zero_row`)

```cpp
const double* row_zero = current.data();
const double* row_one = current.data() + constants.num_stock_nodes;
double* out_row = next.data();
double variance_drift = constants.variance_drift[0];

for (int stock_i = 1; stock_i < constants.num_stock_nodes - 1; ++stock_i) {
    double centre = row_zero[stock_i];
    double dV_dS = (row_zero[stock_i + 1] - row_zero[stock_i - 1]) * constants.inv_two_ds;
    double dV_dv = (row_one[stock_i] - centre) * constants.inv_dv;
    ...
    out_row[stock_i] = centre + constants.dt * change;
}
```

This block is §7.7 with the same row-pointer and reciprocal changes as the
interior. `variance_drift[0]` is κ(θ − 0), which is κθ, so the table built in
§8.3 already holds the value this row needs.

**8.8 The stock edges** (`apply_stock_edges`)

```cpp
void apply_stock_edges(const StencilConstants& constants, std::vector<double>& next, int row_begin, int row_end, double time_to_expiry) {
    double strike_today = constants.strike * std::exp(-constants.rate * time_to_expiry);
    ...
    for (int var_j = row_begin; var_j < row_end; ++var_j) {
        next[var_j * constants.num_stock_nodes] = zero_stock_value;
        next[var_j * constants.num_stock_nodes + last_stock] = max_stock_value;
    }
}
```

This block is §7.8, but it only writes rows `row_begin` to `row_end - 1`. The
serial rung passes the whole grid. The row range is there for the MPI rung,
where each process only owns some of the rows.

**8.9 The top row** (`copy_top_row`)

```cpp
int top_row = constants.num_variance_nodes - 1;
for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
    next[top_row * num_stock_nodes + stock_i] = next[(top_row - 1) * num_stock_nodes + stock_i];
}
```

This block is §7.9 written with the index formula by hand. It must still run
last, because it copies values the interior update has just written.

**8.10 The column range** (`update_cells_scalar`, `update_row_scalar`)

```cpp
void update_cells_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j, int first_stock, int end_stock) {
    ...
    for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {
    ...
}

void update_row_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    update_cells_scalar(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}
```

This block splits the scalar row in two. `update_cells_scalar` holds the
maths from §8.4 to §8.6 and updates only the columns it is given, and
`update_row_scalar` asks it for every interior column. We do it because the
AVX2 row (§9a) handles cells in groups of four and needs the scalar maths for
the one to three cells left at the end of the row. Sharing one function means
those leftover cells are computed by exactly the code the serial rung uses,
and the serial rung's answers did not change.

---

## 9. The serial rung (`solver_serial.cpp`, `SerialSolver::solve`)

```cpp
Grid grid(cfg.grid, cfg.option);
grid.init_payoff(cfg.option);
double dt = cfg.option.maturity_years / num_timesteps;
report_stability(grid, cfg, dt, result);

StencilConstants constants(grid, cfg, dt);

for (int step = 1; step <= num_timesteps; ++step) {
    for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
        update_row_scalar(constants, current, next, var_j);
    }
    update_zero_row(constants, current, next);
    apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
    copy_top_row(constants, next);
    grid.swap_buffers();
}
```

This block is the five-step timestep from §1.3, one function call per step.
The setup before the loop and the timing, snapshot and `extract_result`
after it are the same as in the baseline (§7.10 and §7.11). We split it this
way because the next rungs only change this loop. The SIMD rung swaps in a
vector version of the row function, OpenMP adds one pragma above the `var_j`
loop, and MPI runs the same loop over only its own rows.

---

## 9a. The SIMD rung (`stencil_avx.cpp`, `solver_simd.cpp`)

The SIMD rung does the same maths as §8.6, but on four neighbouring cells at
once. STUDY_GUIDE §15 explains what SIMD and FMA are. This section walks
through the code.

**9a.1 Build only with AVX2** (`stencil_avx.cpp`)

```cpp
#if defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>
...
#else
bool avx_row_is_vectorised() {
    return false;
}
void update_row_avx(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    update_row_scalar(constants, current, next, var_j);
}
#endif
```

This block picks one of two versions of the file when it is compiled. The
compiler defines `__AVX2__` and `__FMA__` only when it was given
`-mavx2 -mfma`, which the Makefile does on x86_64 machines such as rangpur
(§13.5). On the arm64 Mac neither is defined, so `update_row_avx` simply
calls the scalar row and the program still builds and passes its tests. We
do it this way because the Mac cannot run AVX2 natively, and one `#if`
around the whole file is easier to read than one around every line.

**9a.2 Fill a register with one constant**

```cpp
__m256d dt = _mm256_set1_pd(constants.dt);
__m256d rate = _mm256_set1_pd(constants.rate);
...
__m256d two = _mm256_set1_pd(2.0);
__m256d stock_diffusion = _mm256_set1_pd(constants.stock_diffusion[var_j]);
...
```

This block copies each constant into all four slots of a 256-bit register.
`__m256d` is the type of one such register, holding four doubles, and
`set1` puts the same number in every slot. We do it because every AVX
instruction works on whole registers, so to multiply four cells by `dt` we
need `dt` four times. The per-row coefficients are the same for every cell in
a row, so they are filled once per row and not once per cell. These are local
variables for the same reason as §8.4.

**9a.3 Four cells per pass**

```cpp
int last_interior = num_stock_nodes - 2;
int stock_i = 1;
while (stock_i + 3 <= last_interior) {
    ...
    stock_i = stock_i + 4;
}
update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
```

This block walks along the row four cells at a time, `stock_i` to
`stock_i + 3`, for as long as four whole interior cells are left. The cells
that do not fill a register are finished by the scalar maths from §8.10. We
do it because the number of interior cells, `ns - 2`, is usually not a
multiple of four. For the 2048-wide benchmark grid there are 2046 interior
cells, which is 511 groups of four and then two left over.

```text
  row var_j   col: 0 | 1  2  3  4 | 5  6  7  8 | ... | 2041 .. 2044 | 2045 2046 | 2047
                edge | AVX pass 1 | AVX pass 2 |     | AVX pass 511 |  scalar   | edge
```

**9a.4 Loading the neighbours**

```cpp
__m256d centre = _mm256_loadu_pd(row_centre + stock_i);
__m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
__m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);
__m256d north = _mm256_loadu_pd(row_above + stock_i);
__m256d south = _mm256_loadu_pd(row_below + stock_i);
...
```

This block loads four consecutive doubles into each register. `east` is the
same four-wide window as `centre` shifted one column right, so lane 0 of
`east` is the east neighbour of lane 0 of `centre`, and so on for every lane.
We use `loadu`, the unaligned load, because the windows start at
`stock_i - 1`, `stock_i` and `stock_i + 1`. Those three addresses are 8 bytes
apart, so at most one of them can ever sit on a 32-byte boundary, and the
aligned load `_mm256_load_pd` would crash on the others.

```text
  row_centre:  ... | w0 | c0 | c1 | c2 | c3 | e3 | ...
  west   = [w0 c0 c1 c2]       starts at stock_i - 1
  centre =    [c0 c1 c2 c3]    starts at stock_i
  east   =       [c1 c2 c3 e3] starts at stock_i + 1
```

**9a.5 The derivatives**

```cpp
__m256d east_minus_west = _mm256_sub_pd(east, west);
__m256d dV_dS = _mm256_mul_pd(east_minus_west, inv_two_ds);
...
__m256d east_plus_west = _mm256_add_pd(east, west);
__m256d stock_curvature = _mm256_fnmadd_pd(two, centre, east_plus_west);
__m256d d2V_dS2 = _mm256_mul_pd(stock_curvature, inv_ds_sq);
...
__m256d corners = _mm256_sub_pd(north_east, north_west);
corners = _mm256_sub_pd(corners, south_east);
corners = _mm256_add_pd(corners, south_west);
__m256d d2V_dS_dv = _mm256_mul_pd(corners, inv_four_ds_dv);
```

This block is the five central differences from §8.6, one intrinsic per
arithmetic operation. `_mm256_sub_pd(a, b)` is `a - b` in all four lanes,
`_mm256_mul_pd` is `a * b`, and `_mm256_fnmadd_pd(a, b, c)` is `c - a * b`
done as one instruction with one rounding. We name every intermediate so
each line does one thing, which keeps the AVX version readable next to the
scalar version.

**9a.6 The PDE terms with FMA**

```cpp
__m256d stock_price = _mm256_loadu_pd(stock_price_table + stock_i);
__m256d stock_price_sq = _mm256_loadu_pd(stock_price_sq_table + stock_i);
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
_mm256_storeu_pd(out_row + stock_i, new_value);
```

This block adds up the five PDE terms and the discounting, then steps back
one `dt`, for four cells. S and S² differ per column, so they are loaded from
the tables, while the coefficients are the registers from §9a.2.
`_mm256_fmadd_pd(a, b, c)` is `a * b + c` in one instruction, so each line
adds one term to the running total `change`. We use FMA because it does a
multiply and an add for the cost of one instruction and rounds once instead
of twice. That single rounding is also why the SIMD answers differ from the
serial ones in the last bits, which the tests measure (§12.3, §12.4).

**9a.7 The SIMD time loop** (`solver_simd.cpp`)

```cpp
if (!avx_row_is_vectorised()) {
    std::fprintf(stderr, "[simd] built without AVX2 and FMA, using the scalar row\n");
}
...
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
```

This block is the serial time loop from §9 with one line changed,
`update_row_scalar` became `update_row_avx`. The warning goes to stderr so
nobody mistakes a Mac run of the fallback for an AVX2 result. The boundaries
stay scalar because they are a few thousand cells per step against about a
million in the interior.

---

## 9b. The OpenMP rung (`solver_omp.cpp`, `OmpSolver::solve`)

The OpenMP rung is the SIMD time loop with one pragma added. STUDY_GUIDE §16
explains threads, OpenMP and why no locks are needed.

**9b.1 Record the thread count**

```cpp
#include <omp.h>
...
result.threads = omp_get_max_threads();
```

This block asks the OpenMP runtime how many threads the next parallel loop
will use and stores it in the result, so it ends up in the CSV line. We
record what the program actually got, not what we asked for, because a typo
in `OMP_NUM_THREADS` or a slurm limit can silently give a different number.

**9b.2 The pragma**

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
```

This block is the only real change from §9a.7. `parallel for` starts a team
of threads and gives each one a share of the `var_j` loop. `schedule(static)`
cuts the rows into one equal band per thread before the loop starts, which
suits this loop because every row costs the same. `default(none)` makes the
compiler refuse any variable whose sharing we did not state, and `shared(...)`
lists the four that every thread uses. The loop counter `var_j` is private to
each thread automatically. We do it this way because the only thing that
changes between threads is which rows they own, and each thread writes only
its own rows of `next` while everyone only reads `current`, so no locks are
needed.

```text
  OMP_NUM_THREADS=4, nv = 512, interior rows 1 to 510 (510 rows)

  row 510 ┐
   ...    │ thread 3   rows 384 .. 510   (127 rows)
  row 384 ┘
  row 383 ┐
   ...    │ thread 2   rows 257 .. 383   (127 rows)
  row 257 ┘
  row 256 ┐
   ...    │ thread 1   rows 129 .. 256   (128 rows)
  row 129 ┘
  row 128 ┐
   ...    │ thread 0   rows 1 .. 128     (128 rows)
  row 1   ┘
```

510 rows do not split evenly four ways, so the first two threads get one
extra row each. This split was printed by a small test program under both
Apple clang's libomp and g++ 8.5's libgomp, and they agree.

**9b.3 The barrier and the serial boundaries**

```cpp
auto boundary_start = std::chrono::steady_clock::now();
update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
grid.swap_buffers();
auto boundary_end = std::chrono::steady_clock::now();
result.boundary_seconds += std::chrono::duration<double>(boundary_end - boundary_start).count();
```

This block runs after the parallel loop, on the main thread only. A
`parallel for` ends with an implicit barrier, so no thread gets past the loop
until every row is done, and that is what makes it safe for `copy_top_row`
to read row `nv - 2` and for `swap_buffers` to flip the sheets. The
boundaries are left serial because they are a few thousand cells against a
million. That serial part is the Amdahl term that the strong-scaling sweep
measures. The two clock reads around it add up that serial time in
`boundary_seconds`, so the benchmarks can show it next to the parallel part.
Two clock reads per step cost well under a microsecond against a step of
about a millisecond.

---

## 9c. The MPI rung (`solver_mpi.cpp`, `MpiSolver::solve`)

The MPI rung runs several copies of the program, called processes or ranks,
possibly on different machines. Each rank owns one slab of variance rows and
runs the OpenMP loop over its own rows only. STUDY_GUIDE §17 explains MPI,
halos and the hybrid layout.

**9c.1 Who am I**

```cpp
int rank = 0;
int num_ranks = 1;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
result.ranks = num_ranks;
result.threads = omp_get_max_threads();
```

This block asks MPI which rank this copy is (0, 1, 2 ...) and how many
there are in total. `MPI_COMM_WORLD` means every process started by the same
`mpiexec`. We do it because every later decision, which rows to update and
who to swap with, depends on these two numbers.

**9c.2 Sharing out the rows** (`slab_for_rank`)

```cpp
struct RowSlab {
    int begin = 0;
    int end = 0;
};

static RowSlab slab_for_rank(int num_rows, int rank, int num_ranks) {
    int rows_per_rank = num_rows / num_ranks;
    int extra_rows = num_rows % num_ranks;
    RowSlab slab;
    if (rank < extra_rows) {
        slab.begin = rank * (rows_per_rank + 1);
        slab.end = slab.begin + rows_per_rank + 1;
    } else {
        slab.begin = extra_rows * (rows_per_rank + 1) + (rank - extra_rows) * rows_per_rank;
        slab.end = slab.begin + rows_per_rank;
    }
    return slab;
}
```

This block works out which rows a rank owns, counting every row from 0 to
`nv - 1`. When the rows do not divide evenly the first ranks get one extra
row each, the same rule OpenMP's static schedule uses. We share out every
row, not only the interior ones, so rank 0 also owns the v = 0 row and the
last rank owns the top row, and each boundary has exactly one owner.

```text
  nv = 32 rows, 3 ranks (32 = 11 + 11 + 10)

  rows 22 .. 31   rank 2   owns the top row 31
  rows 11 .. 21   rank 1
  rows  0 .. 10   rank 0   owns the v = 0 row
```

**9c.3 Checks before the loop**

```cpp
if (cfg.dump_every > 0) {
    throw std::runtime_error("--dump-every is not supported by the mpi solver");
}
...
if (num_variance_nodes / num_ranks < 2) {
    throw std::runtime_error("too many processes, each one needs at least two variance rows");
}
```

This block refuses two set-ups the rung cannot handle. A snapshot on rank 0
would hold only its own rows, so it would be a wrong picture. The last rank's
`copy_top_row` reads row `nv - 2`, so the last slab needs at least two rows.
Every rank runs the same check, so they all stop together instead of some
waiting forever for a message.

**9c.4 This rank's part of each step**

```cpp
int first_interior = std::max(slab.begin, 1);
int end_interior = std::min(slab.end, num_variance_nodes - 1);
bool owns_zero_row = (rank == 0);
bool owns_top_row = (rank == num_ranks - 1);
...
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, first_interior, end_interior)
for (int var_j = first_interior; var_j < end_interior; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
if (owns_zero_row) {
    update_zero_row(constants, current, next);
}
apply_stock_edges(constants, next, slab.begin, slab.end, step * dt);
if (owns_top_row) {
    copy_top_row(constants, next);
}
```

This block is the OpenMP step from §9b with the row range cut down to this
rank's slab. `first_interior` and `end_interior` skip row 0 and the top row,
which are not interior rows. The boundaries are done only by the rank that
owns them, and `apply_stock_edges` takes this rank's row range, which is why
it was given one in §8.8. We do it this way so the maths is exactly the code
the other rungs run, and only the rows each copy works on change.

**9c.5 The halo exchange** (`exchange_halo_rows`)

```cpp
double* first_row = sheet.data() + slab.begin * num_stock_nodes;
double* last_row = sheet.data() + (slab.end - 1) * num_stock_nodes;
double* halo_below = sheet.data() + (slab.begin - 1) * num_stock_nodes;
double* halo_above = sheet.data() + slab.end * num_stock_nodes;
std::vector<MPI_Request> requests;

if (rank > 0) {
    MPI_Request receive_below;
    MPI_Irecv(halo_below, num_stock_nodes, MPI_DOUBLE, rank - 1, halo_tag, MPI_COMM_WORLD, &receive_below);
    requests.push_back(receive_below);
    MPI_Request send_below;
    MPI_Isend(first_row, num_stock_nodes, MPI_DOUBLE, rank - 1, halo_tag, MPI_COMM_WORLD, &send_below);
    requests.push_back(send_below);
}
if (rank < num_ranks - 1) {
    ... the same with rank + 1, halo_above and last_row ...
}
MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
```

This block swaps edge rows with the neighbours. Each rank sends its first
row down and its last row up, and receives the row just outside each end of
its slab. The stencil reads one row up and one row down (STUDY_GUIDE §7.2),
so one row from each neighbour is all a rank needs. `MPI_Irecv` and
`MPI_Isend` only start a transfer and return straight away, and
`MPI_Waitall` waits until all of them have finished. The receives are posted
first so an arriving row can be copied straight into place. The four
transfers run at the same time, and no rank can get stuck waiting for
another that is itself waiting, which can happen with blocking sends.

```text
  rank 1's copy of the sheet         rank 0's copy of the sheet

  row 11  first row ── send ──────►  row 11  halo_above
  row 10  halo_below ◄────── send ── row 10  last row
```

**9c.6 Where the exchange happens, and the timer**

```cpp
auto halo_start = std::chrono::steady_clock::now();
exchange_halo_rows(next, num_stock_nodes, slab, rank, num_ranks);
auto halo_end = std::chrono::steady_clock::now();
result.halo_seconds += std::chrono::duration<double>(halo_end - halo_start).count();

grid.swap_buffers();
```

This block swaps the halo rows of `next` at the end of each step, just
before `next` becomes `current`. So when the next step starts, every rank
already holds up-to-date copies of the two rows it reads from its
neighbours. Doing it on `next` also means writing into a sheet the solver is
allowed to change, since `current` is only ever read. The halo time is added
up separately so the benchmarks can show how much of each step is
communication. It is still part of `seconds`, because the step is not done
until the exchange is.

**9c.7 Collect the answer** (`MPI_Gatherv`)

```cpp
for (int other_rank = 0; other_rank < num_ranks; ++other_rank) {
    RowSlab other_slab = slab_for_rank(num_variance_nodes, other_rank, num_ranks);
    counts[other_rank] = (other_slab.end - other_slab.begin) * num_stock_nodes;
    offsets[other_rank] = other_slab.begin * num_stock_nodes;
}
const double* my_rows = current.data() + slab.begin * num_stock_nodes;
int my_count = (slab.end - slab.begin) * num_stock_nodes;
MPI_Gatherv(my_rows, my_count, MPI_DOUBLE, next.data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

grid.swap_buffers();
if (rank == 0) {
    extract_result(grid, cfg, result);
}
```

This block copies every rank's finished rows into rank 0's spare sheet,
then swaps so that sheet becomes `current` for `extract_result`. `Gatherv`
is the gather where each rank can send a different amount, which we need
because the slabs can differ by a row. `counts` says how many doubles come
from each rank and `offsets` where they go. We do it because the price and
Greeks are read off cells that may belong to any rank, and keeping the full
grid on every rank means rank 0 can reuse the same `extract_result` as every
other rung. It happens once, after the timer, so it is not timed.

---

## 9d. The experiment solvers (`solver_experiments.h`, `solver_experiments.cpp`)

These five solvers are not rungs of the ladder. Each is a copy of the
`omp` time loop with exactly one thing changed, so a benchmark can measure
that one thing. They live in their own files so the ladder files stay as
short as the video shows them. Every one is in `test_ladder_matches`.
`omp_scalar` uses the scalar row and `omp_float` uses floats, so neither
can match `simd` bit for bit, and every other one must.

| `--solver` | Changed from `omp` | Experiment |
|---|---|---|
| `omp_scalar` | `update_row_scalar` instead of `update_row_avx` | the scalar line on the strong scaling figure |
| `omp_persistent` | one parallel region around the whole time loop | E4 |
| `omp_runtime` | `schedule(runtime)`, so `OMP_SCHEDULE` picks the schedule | E5 |
| `omp_edges` | the v = 0 row and the stock edges inside the parallel loop | E8 |
| `omp_float` | both sheets in `float`, eight cells per AVX register (`solver_float.cpp`, §9e) | E3 |

**9d.1 One parallel region for the whole loop** (`OmpPersistentSolver`)

```cpp
#pragma omp parallel default(none) shared(constants, current, next, grid, num_variance_nodes, num_timesteps, dt)
{
    for (int step = 1; step <= num_timesteps; ++step) {
        #pragma omp for schedule(static)
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }
        #pragma omp single
        {
            update_zero_row(constants, current, next);
            apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
            copy_top_row(constants, next);
            grid.swap_buffers();
        }
    }
}
```

This block starts the team of threads once, and every thread runs the whole
step loop. `omp for` shares out the rows of each step and ends with a
barrier. `omp single` lets exactly one thread do the boundaries and the swap
while the others wait at its end, which is a second barrier. We try it
because the `omp` rung forks and joins once per step, and this version
swaps that for two barriers per step inside one region. Which is cheaper is
a measurement, not a rule.

**9d.2 Boundaries inside the loop** (`OmpParallelEdgesSolver`)

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes, time_to_expiry)
for (int var_j = 0; var_j < num_variance_nodes - 1; ++var_j) {
    if (var_j == 0) {
        update_zero_row(constants, current, next);
    } else {
        update_row_avx(constants, current, next, var_j);
    }
    apply_stock_edges(constants, next, var_j, var_j + 1, time_to_expiry);
}
copy_top_row(constants, next);
```

This block moves row 0 and the stock edges into the parallel loop. Row 0
reads only `current`, and the edges are columns 0 and `ns - 1`, which the
interior update never writes, so neither needs to wait for any other row.
Only `copy_top_row` stays after the barrier, because it reads the finished
row `nv - 2`. We try it because it shrinks the serial part of each step,
which is the Amdahl term.

## 9e. The float experiment (`solver_float.cpp`, float parts of `stencil.cpp` and `stencil_avx.cpp`)

`--solver omp_float` is the `omp` rung with both sheets stored as `float`
(4 bytes, about 7 significant digits) instead of `double` (8 bytes, about
16). It is experiment E3. Run on one thread it compares with `simd`, and on
24 threads with `omp`. It is not a rung, because it does not match the
baseline to 1e-11 and never can.

**9e.1 The constants, rounded once** (`FloatStencilConstants`, `stencil.cpp`)

```cpp
FloatStencilConstants::FloatStencilConstants(const StencilConstants& constants) {
    num_stock_nodes = constants.num_stock_nodes;
    num_variance_nodes = constants.num_variance_nodes;
    dt = (float)constants.dt;
    rate = (float)constants.rate;
    ...
    stock_price.assign(constants.stock_price.begin(), constants.stock_price.end());
    ...
}
```

This block builds the float constants from the double ones. `(float)x`
rounds a double to the nearest float, and `assign` copies a whole table,
rounding every entry. We start from the double constants so both solvers
use the same numbers up to the last rounding, and any gap in the answer
comes from doing the arithmetic in float, not from different set-up.

**9e.2 The float rows** (`update_cells_scalar_float` in `stencil.cpp`, `update_row_avx_float` in `stencil_avx.cpp`)

```cpp
while (stock_i + 7 <= last_interior) {
    __m256 centre = _mm256_loadu_ps(row_centre + stock_i);
    ...
    change = _mm256_fmadd_ps(cross_factor, d2V_dS_dv, change);
    ...
    _mm256_storeu_ps(out_row + stock_i, new_value);
    stock_i = stock_i + 8;
}
update_cells_scalar_float(constants, current, next, var_j, stock_i, last_interior + 1);
```

This block is the AVX2 row of §9a line for line, with `__m256` in place of
`__m256d` and `_ps` (packed single) in place of `_pd` (packed double). A
256-bit register holds eight floats, so the loop takes eight cells per pass
and leaves zero to seven for the scalar float row. We keep the two versions
line for line the same so the only thing E3 measures is the number type.
The scalar float row lives in `stencil.cpp` for the same reason the double
one does, so it is never built with `-mfma`.

**9e.3 The float sheets** (`OmpFloatSolver::solve`)

```cpp
std::vector<float> current(grid.current().begin(), grid.current().end());
std::vector<float> next(current.size(), 0.0f);
...
std::swap(current, next);
...
std::vector<double>& answer = grid.next();
for (int cell = 0; cell < (int)current.size(); ++cell) {
    answer[cell] = current[cell];
}
grid.swap_buffers();
extract_result(grid, cfg, result);
```

This block keeps the solver's own two float sheets next to the grid's
double ones. The payoff is built in double and rounded into `current`,
the time loop works only on the float sheets, and at the end the answer
is copied back into the grid's spare double sheet and made current, so
`extract_result` reads the price and Greeks exactly as for every other
solver. We do it this way so the grid class, the result code and the
other solvers did not have to change at all.

**9e.4 The float boundaries** (`update_zero_row_float`, `apply_edges_and_top_float`)

```cpp
double strike_today = constants.strike * std::exp(-constants.rate * time_to_expiry);
...
next[var_j * num_stock_nodes] = (float)zero_stock_value;
next[var_j * num_stock_nodes + last_stock] = (float)max_stock_value;
```

This block writes the stock edges and copies the top row on the float
sheet. The edge values come from the double formula of `apply_stock_edges`
and are then rounded, because `std::exp` of a float would round twice. We
copy the boundary code into this file instead of changing the shared
double functions, so the rungs that must match `simd` bit for bit are not
touched by the experiment.

---

## 10. Output (`io.cpp`)

**10.1 The result line** (`print_result_csv`)

```cpp
std::printf("%.10g,%.8g,%.8g,%.8g,%d,%d,%d,%s,%.6f,%.6g,%d,%d,%.6f,%.6f\n",
            result.price, result.delta, result.gamma, result.vega,
            cfg.grid.num_stock_nodes, cfg.grid.num_variance_nodes, cfg.grid.num_timesteps,
            cfg.solver.c_str(), result.seconds, result.cell_updates_per_sec,
            result.threads, result.ranks, result.halo_seconds, result.boundary_seconds);
```

This block prints one CSV line per solve with the price, the Greeks, the grid
size, the solver name, the timing, and the number of threads and processes.
The benchmark and plotting scripts parse these columns by position, so new
columns only ever go on the end, which is where `threads`, `ranks`,
`halo_seconds` and `boundary_seconds` went.

**10.2 Snapshots** (`dump_snapshot`)

```cpp
std::filesystem::create_directories(dump_dir);
char file_name[32];
std::snprintf(file_name, sizeof file_name, "snap_%06d.csv", step);
std::ofstream out(dump_dir + "/" + file_name);
...
for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
    for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
        if (stock_i > 0) {
            out << ',';
        }
        out << grid.current()[grid.index(stock_i, var_j)];
    }
    out << '\n';
}
```

This block writes the whole current sheet to a CSV file, one variance row per
line. The animation scripts read these files, and the step number is
zero-padded so the files sort in time order. The worked example in
STUDY_GUIDE §7.4 was checked against one of these snapshots.

---

## 11. The Black-Scholes formula (`black_scholes.cpp`)

```cpp
double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double bs_price(bool is_call, double spot, double strike, double rate, double div_yield, double sigma, double maturity_years) {
    double sqrt_maturity = std::sqrt(maturity_years);
    double d1 = (std::log(spot / strike) + (rate - div_yield + 0.5 * sigma * sigma) * maturity_years) / (sigma * sqrt_maturity);
    double d2 = d1 - sigma * sqrt_maturity;
    double stock_leg = spot * std::exp(-div_yield * maturity_years);
    double strike_leg = strike * std::exp(-rate * maturity_years);
    if (is_call) {
        return stock_leg * norm_cdf(d1) - strike_leg * norm_cdf(d2);
    } else {
        return strike_leg * norm_cdf(-d2) - stock_leg * norm_cdf(-d1);
    }
}
```

This block is the textbook Black-Scholes price with a dividend yield. Only
the tests use it, as the exact answer the PDE solver must match when the
variance is held fixed. `norm_cdf` uses `erfc` because N(x) = ½·erfc(−x/√2)
is an exact identity, so there is no approximation to defend.

---

## 12. The tests (`tests/`)

All three tests return 0 on a pass and 1 on a fail, so `make test` stops at
the first failure. STUDY_GUIDE §10 explains what each one proves.

**12.1 Black-Scholes collapse** (`test_bs_collapse.cpp`)

```cpp
cfg.heston.xi = 0.0;
cfg.heston.v0 = cfg.heston.theta;
...
SolveResult result = solver.solve(cfg);
double closed_form = bs_price(cfg.option.is_call, cfg.market.spot, cfg.option.strike, cfg.market.rate, cfg.market.div_yield, sigma, cfg.option.maturity_years);
double rel_err = std::fabs(result.price - closed_form) / closed_form;
```

This block switches off the random variance and compares the solver with the
exact formula, for both a call and a put. With `xi = 0` and `v0 = theta` the
variance never moves, so Heston becomes Black-Scholes and the two must agree
to within the measured tolerance of 2.5e-3.

**12.2 Put-call parity** (`test_parity.cpp`)

```cpp
cfg.option.is_call = true;
SolveResult call_result = solver.solve(cfg);
cfg.option.is_call = false;
SolveResult put_result = solver.solve(cfg);

double parity_rhs = cfg.market.spot * std::exp(-cfg.market.div_yield * maturity) - cfg.option.strike * std::exp(-cfg.market.rate * maturity);
double parity_lhs = call_result.price - put_result.price;
double abs_gap = std::fabs(parity_lhs - parity_rhs);
```

This block prices a call and a put on the full model and checks
C − P = S·e^(−qT) − K·e^(−rT). Parity holds in any model, so it tests the
cross term and the variance boundaries that the collapse test switches off.

Both of these tests first check that `dt` is below `dt_stable_estimate`. An
unstable run would make the comparison meaningless, so it fails with an
UNSTABLE message instead.

**12.3 Every rung matches the baseline** (`test_ladder_matches.cpp`)

```cpp
SolveResult reference = baseline.solve(cfg);
for (const std::string& rung : kRungs) {
    std::unique_ptr<Solver> solver = make_solver(rung);
    SolveResult candidate = solver->solve(cfg);
    double worst_diff = max_relative_diff(reference, candidate);
    bool ok = worst_diff < kRelativeTolerance;
```

This block runs the baseline and every other rung on the same small grid and
compares price, delta, gamma and vega. A rung is not benchmarked until this
passes, because a fast wrong answer is worth nothing. Each new rung is added
to `kRungs`.

```cpp
static const std::vector<std::string> kSameBitsAsSimd = {"omp", "mpi", "omp_persistent", "omp_runtime", "omp_edges"};
...
for (const std::string& rung : kSameBitsAsSimd) {
    SolveResult candidate = make_solver(rung)->solve(cfg);
    double simd_diff = max_relative_diff(simd_reference, candidate);
    bool ok = simd_diff == 0.0;
```

This block checks a stronger rule for the OpenMP and MPI rungs and the experiments built on them. They must
give exactly the same bits as the SIMD rung, not just agree to 1e-11. Each cell is still
computed by the same `update_row_avx` code, only by a different thread or
process, and nothing is ever added up across them, so any difference at all
would mean two workers touched the same cell, a row was skipped, or a halo
row was stale. With the halo exchange removed on purpose, the MPI answer on
two ranks moved by only 2e-9, so the exact check gives a far clearer signal
than the tolerance alone.

```cpp
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ...
            SolveResult candidate = solver->solve(cfg);
            if (rank != 0) {
                continue;
            }
    ...
    MPI_Finalize();
```

```cpp
static const double kFloatTolerance = 1e-4;
...
SolveResult float_candidate = make_solver("omp_float")->solve(cfg);
double float_diff = max_relative_diff(reference, float_candidate);
bool float_ok = float_diff < kFloatTolerance;
```

This block holds the float experiment (§9e) to its own limit of 1e-4. A
float keeps about 7 digits, so it can never meet 1e-11, and it is not put
in `kRungs`. On this grid it measures 2.1e-6 on the call and 8.2e-6 on the
put on rangpur (8.3e-6 on the Mac), so 1e-4 leaves room for rounding but still fails a real bug, which
moves the answer by far more.

This block lets the same test run as one process or under `mpiexec`. Every
rank runs every solver, because the MPI rung needs all ranks to take part,
but only rank 0 holds the gathered answer, so only rank 0 compares and
prints.

```cpp
static double relative_diff(double a, double b) {
    double scale = std::fmax(std::fabs(a), std::fabs(b));
    if (scale == 0.0) {
        return 0.0;
    }
    return std::fabs(a - b) / scale;
}
```

This block measures a difference relative to the size of the numbers. We use
a relative difference because gamma is about 0.0008 and vega about 2000, and
one tolerance of 1e-11 has to be fair to both.

**12.4 The AVX2 row against the scalar row** (`test_avx_row.cpp`)

```cpp
if (!avx_row_is_vectorised()) {
    std::printf("test_avx_row: FAIL, built without AVX2 and FMA\n");
    return 1;
}
...
std::vector<int> row_lengths = {130, 131, 132, 133, 2048};
...
std::vector<double> next_scalar(num_cells, -7.0);
std::vector<double> next_avx(num_cells, -7.0);
for (int var_j = 1; var_j < cfg.grid.num_variance_nodes - 1; ++var_j) {
    update_row_scalar(constants, current, next_scalar, var_j);
    update_row_avx(constants, current, next_avx, var_j);
}
```

This block fills a sheet with a call payoff plus random noise, runs both row
functions on it, and compares every cell of the two outputs. The row lengths
130 to 133 leave zero, one, two and three cells for the scalar remainder, so
every path through §9a.3 is tested. Both outputs start filled with -7, so if
the AVX row wrote outside the interior the edges would differ and the test
would fail. It refuses to pass when the AVX2 row was not built, because then
it would only be comparing the scalar row with itself. It runs with
`make avx-check` (§13.5), not `make test`, because on the Mac it needs a
separate x86_64 build.

```cpp
static const double kFloatTolerance = 1e-5;
...
std::vector<int> float_row_lengths = {130, 131, 132, 133, 134, 135, 136, 137, 2048};
...
update_cells_scalar_float(float_constants, current, next_scalar, var_j, 1, num_stock_nodes - 1);
update_row_avx_float(float_constants, current, next_avx, var_j);
```

This block runs the same check on the float row (§9e). The float row takes
eight cells per pass, so the lengths 130 to 137 leave every remainder from
zero to seven. The limit is 1e-5 of the largest value, because one float
rounding is about 6e-8 of a value where one double rounding is about 1e-16.
It measured 6.3e-8 on the short rows and 1.2e-6 on the 2048 row.

---

## 13. The build (`Makefile`)

**13.1 Compiler and flags**

```make
CXX      = mpicxx
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude

ifneq ($(shell $(CXX) --version 2>/dev/null | grep -c clang),0)
  LIBOMP   := $(shell brew --prefix libomp 2>/dev/null)
  OMPFLAGS  = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
  OMPLIBS   = -L$(LIBOMP)/lib -lomp
else
  OMPFLAGS  = -fopenmp
  OMPLIBS   = -fopenmp
endif
```

This block picks the compiler and flags. Everything is built with `mpicxx`
because the final program is one binary for every rung, including MPI.
OpenMP needs different flags on Apple clang (with Homebrew's libomp) and on
GCC, and they sit in their own variable so a benchmark script that overrides
`CXXFLAGS` does not lose them. `-O2` is fixed for every benchmark.

**13.2 Objects and targets**

```make
CORE_OBJS = src/params.o src/grid.o src/stencil.o src/stencil_avx.o \
            src/solver.o src/solver_baseline.o src/solver_serial.o \
            src/solver_simd.o src/black_scholes.o src/io.o

TESTS = tests/test_bs_collapse tests/test_parity tests/test_ladder_matches

heston: src/main.o $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@ $(OMPLIBS) $(FSLIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@
```

This block says which object files make up the program. The program and all
three tests link the same `CORE_OBJS`, so the tests always check exactly the
code that runs. The `%.o: %.cpp` rule compiles any `.cpp` into its `.o`, and
make only reruns it for files that changed.

**13.3 Running the tests**

```make
test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done
	@echo "ALL TESTS PASS"
```

```make
@for n in 1 2 7; do echo "== tests/test_ladder_matches with OMP_NUM_THREADS=$$n"; OMP_NUM_THREADS=$$n ./tests/test_ladder_matches 2>/dev/null || exit 1; done
```

This block builds and runs every test in turn. `|| exit 1` stops at the
first failure, which is why each test returns 1 when it fails. The second
loop runs the ladder test again with 1, 2 and 7 threads. One thread checks
the pragma does nothing harmful, two checks a real split, and seven does not
divide the 30 interior rows of the test grid, so the bands are uneven and a
mistake at a band edge would show.

**13.4 MPI tests**

```make
MPIEXEC = mpiexec --oversubscribe --bind-to none
test-mpi: tests/test_ladder_matches
	@for n in 2 3; do echo "== tests/test_ladder_matches on $$n processes"; OMP_NUM_THREADS=2 $(MPIEXEC) -n $$n ./tests/test_ladder_matches 2>/dev/null || exit 1; done
```

This block runs the ladder test as 2 and then 3 processes, two threads each.
Three ranks do not divide the test grid's 32 rows, so the slabs are 11, 11
and 10 and the uneven case is tested. `--oversubscribe` lets mpiexec start
more processes than slurm or the laptop has slots for, and `--bind-to none`
stops it pinning each process to one core, which would squeeze its OpenMP
threads onto that core.

**13.5 AVX2 flags and the AVX check**

```make
ifeq ($(shell uname -m),x86_64)
  AVXFLAGS = -mavx2 -mfma
endif

src/stencil_avx.o: src/stencil_avx.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(AVXFLAGS) -c $< -o $@

ifeq ($(shell uname -s)-$(shell uname -m),Darwin-arm64)
  AVXCHECK_CXX = clang++ -arch x86_64
else
  AVXCHECK_CXX = $(CXX)
endif

avx-check:
	$(AVXCHECK_CXX) $(AVXCHECK_FLAGS) -mavx2 -mfma -c src/stencil_avx.cpp -o tests/stencil_avx_check.o
	$(AVXCHECK_CXX) $(AVXCHECK_FLAGS) tests/test_avx_row.cpp src/stencil.cpp src/grid.cpp src/params.cpp tests/stencil_avx_check.o -o tests/test_avx_row
	./tests/test_avx_row
```

This block gives `-mavx2 -mfma` to one file only, `stencil_avx.cpp`, and only
on x86_64. We do it because with `-mfma` g++ is allowed to turn any
`a * b + c` into an FMA, so building every file with it would quietly change
the baseline and serial answers and timings. Keeping the flags on one file
means the rungs below SIMD compile exactly as they did before.
`-mavx2 -mfma` is used rather than `-march=native` so the program needs only
the AVX2 features every rangpur r730 node has. The `avx-check` target builds
the row test for x86_64. On an Apple silicon Mac `clang++ -arch x86_64` makes
an Intel program that runs under Rosetta, which translates AVX2 so the
answers can be checked, though Rosetta timings mean nothing. On rangpur the
same target builds natively.

---

## 14. Practice questions on the code

- Why does `main` read the config file in a separate first pass? (§3.2)
- Why does `parse_double` check how many characters `stod` used? (§4.2)
- Why are `Grid`'s fields private, and why are the getters written inside the
  class? (§5.1)
- Why does `init_payoff` write into `current_` and not `next_`? (§5.5)
- What does `virtual ... = 0` mean, and why is `override` worth writing?
  (§6.2, §6.3)
- Why does `make_solver` return a `unique_ptr<Solver>` and not a
  `SerialSolver`? (§6.4)
- Why do the `current` and `next` references still point at the right sheets
  after `swap_buffers`? (§7.2)
- Why does `copy_top_row` have to run after the interior update? (§7.9)
- Why does the timer stop before `dump_snapshot`? (§7.10)
- Why does `update_row_scalar` copy the constants into local variables?
  (§8.4)
- What do the row pointers replace, and who owns the memory they point at?
  (§8.5)
- Why does the serial rung differ from the baseline in the last digits, and
  how do we know that is fine? (§8.6, §12.3)
- Why does `apply_stock_edges` take a row range? (§8.8)
- Why is the test tolerance relative and not absolute? (§12.3)
- Why is everything compiled with `mpicxx`, even the serial rung? (§13.1)
- What does the `#if defined(__AVX2__)` at the top of `stencil_avx.cpp`
  choose between, and why does the Mac build still pass the tests? (§9a.1)
- Why does every constant get copied into four lanes with `set1`? (§9a.2)
- What happens to the last cells of a row when `ns - 2` is not a multiple of
  four? (§9a.3)
- Why `_mm256_loadu_pd` and not `_mm256_load_pd`? (§9a.4)
- What does `_mm256_fmadd_pd(a, b, c)` compute, and why does it change the
  last bits of the answer? (§9a.6)
- Why do only the interior rows use AVX2, and not the boundaries? (§9a.7)
- Why does `test_avx_row` fill both outputs with -7 first? (§12.4)
- Why does only `stencil_avx.cpp` get `-mfma`? (§13.5)
- What does each word of the OpenMP pragma do? (§9b.2)
- Why does the OpenMP rung need no locks or `critical` sections? (§9b.2)
- Why is it safe for `copy_top_row` to run straight after the parallel
  loop? (§9b.3)
- Why must the OpenMP rung match the SIMD rung bit for bit, and why is that
  a stronger test than a tolerance? (§12.3)
- Why are the tests run at 7 threads? (§13.3)
- Why does `main` hand everything to `run`? (§3.6)
- Which rows does each rank own when 32 rows go to 3 ranks? (§9c.2)
- Why does each rank need exactly one row from each neighbour? (§9c.5)
- Why post `MPI_Irecv` before `MPI_Isend`, and why non-blocking at all?
  (§9c.5)
- Why is the halo exchanged on `next` and not on `current`? (§9c.6)
- Why does every rank keep the whole grid, and what would you change for a
  grid too big for one machine? (§9c.7, STUDY_GUIDE §17)
- Why does only rank 0 compare results in `test_ladder_matches`? (§12.3)
- How many barriers does `omp_persistent` have per step, and where are
  they? (§9d.1)
- Why can row 0 and the stock edges join the parallel loop, but the top row
  cannot? (§9d.2)
