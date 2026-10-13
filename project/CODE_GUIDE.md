# Code Guide. The Heston Pricer, Snippet by Snippet

This guide goes through the C++ one block at a time, in the order the program
runs. It explains what each piece of code does and why it is written that
way. The maths behind it, such as what the PDE means, how finite differences
work and why there are two buffers, is in [STUDY_GUIDE.md](STUDY_GUIDE.md).
This guide points there whenever a block depends on an idea from it.

**Interactive version.** The
[Heston Pricer Study Guide](https://claude.ai/artifact/5ggg9HfwJiXwvRRPwHHSWh)
shows this guide one snippet at a time, after the short concept cards in
[docs/STUDY_CARDS.md](docs/STUDY_CARDS.md). It has progress ticks, search, a
mode that hides each answer until you have tried it yourself, and a drill. It
is built from those two files, so rebuild it when either changes.

**How the snippets work.** Every snippet is copied from the real file, mostly
with its comments removed so it stays short. The file and function are named
above it. Under each one, the first sentence says what the block does and the
second says why we do it that way.

**Following the code.** §1.5 traces one run from `main` to the printed
line, with a link to the exact line of every step, and a table of the lines
where each solver differs from the serial one. Start there when you want to
follow the logic in the editor.

**Reading order.** If you are short on time, read §1 for the map, then §7
(the baseline, where the whole algorithm is in one function), then §8 and §9
(the serial rung). Everything else is the setup and plumbing around them.
§15 puts the code before and after each change in the changelog side by
side, so it is the quickest way to see what each rung added.

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
               SimdSolver, OmpSolver, MpiSolver, make_solver()
                 ├─ grid.h   Grid (two buffers, index, swap)  body in grid.cpp
                 └─ solver_experiments.h   the five experiment solvers

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

§1.5 follows the same run line by line, with a link to every function on
the way.

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

### 1.5 Following the logic in the editor

This section is the map for reading the code itself. Every link opens the
file at that line, in VS Code and on GitHub. Keep the code open next to
this guide and follow the links in order. The section number at the end of
each step points to the block that explains those lines.

**The whole run as a call tree**

```text
main                                        src/main.cpp:104
└─ run                                      src/main.cpp:27
   ├─ load_config                           src/params.cpp:59
   │   └─ apply, once per line of the file  src/params.cpp:14
   ├─ make_solver                           src/solver.cpp:95
   ├─ solver->solve, for example SerialSolver::solve   src/solver_serial.cpp:16
   │   ├─ Grid::Grid                        src/grid.cpp:14
   │   ├─ Grid::init_payoff                 src/grid.cpp:37
   │   ├─ report_stability                  src/solver.cpp:70
   │   ├─ StencilConstants::StencilConstants   src/stencil.cpp:12
   │   ├─ the time loop, once per step
   │   │   ├─ update_row_scalar, once per interior row   src/stencil.cpp:127
   │   │   │   └─ update_cells_scalar       src/stencil.cpp:64
   │   │   ├─ update_zero_row               src/stencil.cpp:133
   │   │   ├─ apply_stock_edges             src/stencil.cpp:156
   │   │   ├─ copy_top_row                  src/stencil.cpp:181
   │   │   ├─ Grid::swap_buffers            src/grid.cpp:31
   │   │   └─ dump_snapshot, only with --dump-every   src/io.cpp:13
   │   └─ extract_result                    src/solver.cpp:34
   │       ├─ stock_position and variance_position   src/grid.cpp:110 and 116
   │       │   └─ locate                    src/grid.cpp:57
   │       └─ bilinear_interpolate, five calls   src/solver.cpp:18
   └─ print_result_csv                      src/io.cpp:45
```

The tree reads top to bottom in the order things happen. A line indented
under another is a call made from inside it. Only `solve` changes from
solver to solver. Everything above it and everything it calls for setup
and for the answer is shared.

**Step 1. From the command line to a solver**

1. [main.cpp:104](src/main.cpp#L104) `main` starts MPI, asks which process
   this is, and hands everything to `run`. `MPI_Finalize` has to run on
   every exit path, so all the work sits in `run` and `main` stays four
   lines long. (§3.6)
2. [main.cpp:31](src/main.cpp#L31) checks the command starts with
   `--config <file>`. The file is loaded first so every flag after it can
   overwrite a value from the file. (§3.1)
3. [params.cpp:59](src/params.cpp#L59) `load_config` starts from the
   defaults in [params.h:60](include/params.h#L60), reads the file one line
   at a time and splits each line into three words at
   [params.cpp:79](src/params.cpp#L79). (§4.3)
4. [params.cpp:94](src/params.cpp#L94) passes each `key = value` to
   [`apply`](src/params.cpp#L14), which stores the value in the matching
   field. An unknown key stops the program, so a typo cannot quietly run
   the wrong grid. (§4.2)
5. [main.cpp:39](src/main.cpp#L39) reads the other flags two words at a
   time and overwrites the matching fields. (§3.2)
6. [main.cpp:80](src/main.cpp#L80) calls
   [`make_solver`](src/solver.cpp#L95), which turns the `--solver` text into
   a solver object and returns it in a `unique_ptr<Solver>`. (§6.4)
7. [main.cpp:88](src/main.cpp#L88) calls `solver->solve(cfg)`. `solve` is
   `virtual` ([solver.h:49](include/solver.h#L49)), so this line jumps into
   the `solve` of whichever subclass `make_solver` built. Step 2 follows it
   there. (§6.2)
8. [main.cpp:91](src/main.cpp#L91) prints the result as one CSV line with
   [`print_result_csv`](src/io.cpp#L45). Only rank 0 prints, so an MPI run
   still gives one line per solve. (§10.1)
9. [main.cpp:95](src/main.cpp#L95) catches any error thrown on the way,
   prints it and returns 1, so a Slurm script sees the failure. (§3.5)

**Step 2. Inside one solve** (the serial rung,
[solver_serial.cpp:16](src/solver_serial.cpp#L16))

Every solver's `solve` has the same three parts, set up, time loop, read
off the answer.

1. [solver_serial.cpp:20](src/solver_serial.cpp#L20) builds the `Grid`.
   The constructor at [grid.cpp:14](src/grid.cpp#L14) works out the two
   spacings and makes both buffers full of zeros. (§5.2)
2. [solver_serial.cpp:21](src/solver_serial.cpp#L21) calls
   [`init_payoff`](src/grid.cpp#L37), which writes the payoff at expiry
   into the current buffer. The solve starts at expiry because that is the
   one time the option's value is known exactly. (§5.5)
3. [solver_serial.cpp:24](src/solver_serial.cpp#L24) works out `dt`, and
   [solver_serial.cpp:25](src/solver_serial.cpp#L25) calls
   [`report_stability`](src/solver.cpp#L70), which prints the largest
   stable `dt` next to it before any maths runs. (§6.5)
4. [solver_serial.cpp:28](src/solver_serial.cpp#L28) builds
   `StencilConstants` once. The constructor at
   [stencil.cpp:12](src/stencil.cpp#L12) does every division and every
   per-row and per-column product before the loop starts. (§8.1 to §8.3)
5. [solver_serial.cpp:34](src/solver_serial.cpp#L34) names the two
   buffers `current` and `next`. (§7.2)
6. [solver_serial.cpp:37](src/solver_serial.cpp#L37) is the time loop.
   Each pass moves the whole grid one `dt` closer to today. Step 3 below
   follows one pass.
7. [solver_serial.cpp:63](src/solver_serial.cpp#L63) calls
   [`extract_result`](src/solver.cpp#L34), which finds today's stock price
   and variance on the grid with [`locate`](src/grid.cpp#L57) and reads the
   price and the three Greeks there. (§5.6, §6.6, §6.7)
8. [solver_serial.cpp:66](src/solver_serial.cpp#L66) works out cell
   updates per second, so grids of different sizes compare on one chart.
   (§7.11)

**Step 3. One timestep** (inside the loop at
[solver_serial.cpp:37](src/solver_serial.cpp#L37))

1. [solver_serial.cpp:42](src/solver_serial.cpp#L42) updates every
   interior row with [`update_row_scalar`](src/stencil.cpp#L127). The
   order of the rows does not matter, because every row reads only
   `current` and writes only its own row of `next`. This loop is the one
   every later rung changes. (§8.10)
2. [solver_serial.cpp:47](src/solver_serial.cpp#L47)
   [`update_zero_row`](src/stencil.cpp#L133) updates the v = 0 row with
   the drift terms only. (§8.7)
3. [solver_serial.cpp:48](src/solver_serial.cpp#L48)
   [`apply_stock_edges`](src/stencil.cpp#L156) writes the known values into
   the first and last stock columns. (§8.8)
4. [solver_serial.cpp:49](src/solver_serial.cpp#L49)
   [`copy_top_row`](src/stencil.cpp#L181) copies row `nv - 2` into the top
   row. It goes last because it reads values that step 1 has just written.
   (§8.9)
5. [solver_serial.cpp:52](src/solver_serial.cpp#L52)
   [`swap_buffers`](src/grid.cpp#L31) makes `next` the new `current`, so
   the next pass reads what this pass wrote. (§5.4)
6. [solver_serial.cpp:54](src/solver_serial.cpp#L54) stops the step timer,
   and [solver_serial.cpp:58](src/solver_serial.cpp#L58) writes a snapshot
   only after the timer has stopped. (§7.10)

**Step 4. Inside one row** (where the maths is)

1. [stencil.cpp:127](src/stencil.cpp#L127) `update_row_scalar` asks
   [`update_cells_scalar`](src/stencil.cpp#L64) for columns 1 to
   `ns - 2`. The first and last columns are the stock edges from step 3.3.
2. [stencil.cpp:67](src/stencil.cpp#L67) copies the constants into local
   variables, so they stay in registers. (§8.4)
3. [stencil.cpp:85](src/stencil.cpp#L85) points at the start of the three
   rows it reads and the one it writes. (§8.5)
4. [stencil.cpp:92](src/stencil.cpp#L92) is the loop along the row. For
   each cell it reads the nine values
   ([stencil.cpp:95](src/stencil.cpp#L95)), works out the five derivatives
   ([stencil.cpp:106](src/stencil.cpp#L106)), the five PDE terms
   ([stencil.cpp:113](src/stencil.cpp#L113)) and the new value
   ([stencil.cpp:120](src/stencil.cpp#L120)). (§8.6)

The SIMD row [`update_row_avx`](src/stencil_avx.cpp#L22) does the same four
things for four cells at a time in the loop at
[stencil_avx.cpp:56](src/stencil_avx.cpp#L56), then hands the last zero to
three cells to `update_cells_scalar` at
[stencil_avx.cpp:116](src/stencil_avx.cpp#L116). (§9a)

**Where each solver differs**

Each `solve` is a copy of the serial one with the lines below changed.
Everything else in it, the set-up, the boundaries, the swap, the timer and
`extract_result`, is the same code.

| `--solver` | `solve` starts at | Row function | The lines that differ | Guide |
|---|---|---|---|---|
| `baseline` | [solver_baseline.cpp:16](src/solver_baseline.cpp#L16) | none, the maths is inline | the whole step is written out at [49 to 122](src/solver_baseline.cpp#L49) | §7 |
| `serial` | [solver_serial.cpp:16](src/solver_serial.cpp#L16) | `update_row_scalar` | the reference for the rows below | §9 |
| `simd` | [solver_simd.cpp:17](src/solver_simd.cpp#L17) | `update_row_avx` | the warning at [21](src/solver_simd.cpp#L21) and the row call at [49](src/solver_simd.cpp#L49) | §9a.7 |
| `omp` | [solver_omp.cpp:15](src/solver_omp.cpp#L15) | `update_row_avx` | the thread count at [19](src/solver_omp.cpp#L19), the pragma at [45](src/solver_omp.cpp#L45), the boundary timer at [53](src/solver_omp.cpp#L53) | §9b |
| `mpi` | [solver_mpi.cpp:85](src/solver_mpi.cpp#L85) | `update_row_avx` | the owned rows at [121](src/solver_mpi.cpp#L121), the loop over them at [136](src/solver_mpi.cpp#L136), the owner checks at [143](src/solver_mpi.cpp#L143), the exchange at [156](src/solver_mpi.cpp#L156), the gather at [177](src/solver_mpi.cpp#L177) | §9c |
| `omp_scalar` | [solver_experiments.cpp:17](src/solver_experiments.cpp#L17) | `update_row_scalar` | the row call at [39](src/solver_experiments.cpp#L39) | §9d.1 |
| `omp_persistent` | [solver_experiments.cpp:60](src/solver_experiments.cpp#L60) | `update_row_avx` | one parallel region at [80 to 98](src/solver_experiments.cpp#L80) | §9d.2 |
| `omp_runtime` | [solver_experiments.cpp:111](src/solver_experiments.cpp#L111) | `update_row_avx` | `schedule(runtime)` at [131](src/solver_experiments.cpp#L131) | §9d.3 |
| `omp_edges` | [solver_experiments.cpp:154](src/solver_experiments.cpp#L154) | `update_row_avx` | the boundaries inside the loop at [176 to 184](src/solver_experiments.cpp#L176) | §9d.4 |
| `omp_float` | [solver_float.cpp:67](src/solver_float.cpp#L67) | `update_row_avx_float` | the float buffers at [90 to 95](src/solver_float.cpp#L90), the float row at [103](src/solver_float.cpp#L103), the copy back at [120 to 124](src/solver_float.cpp#L120) | §9e |

**Where the tests start**

Each test is its own small program with its own `main`.

- [test_bs_collapse.cpp:24](tests/test_bs_collapse.cpp#L24) switches off
  the random variance and compares the baseline with the Black-Scholes
  formula at [black_scholes.cpp:17](src/black_scholes.cpp#L17). (§12.1)
- [test_parity.cpp:23](tests/test_parity.cpp#L23) prices a call and a put
  and checks put-call parity. (§12.2)
- [test_ladder_matches.cpp:64](tests/test_ladder_matches.cpp#L64) runs
  every solver and compares it with the baseline, then with `simd`.
  (§12.3)
- [test_avx_row.cpp:100](tests/test_avx_row.cpp#L100) compares the AVX2
  rows with the scalar rows cell by cell. (§12.4)

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
| `std::size_t` | the `int` that `len()` gives you | it can never be negative, so `0 - 1` wraps round to a huge number |
| `omp_get_wtime()` | `time.perf_counter()` | wall clock time in seconds as a `double`, from OpenMP |

**What `std::size_t` is.** `std::size_t` is a whole-number type for sizes
and positions, such as how many cells a sheet holds or where a `#` sits in
a line. It comes from the C++ standard library. It is defined in the header
`<cstddef>`, and `<string>` and `<vector>` include it for you, which is why
`grid.cpp` can use it without including `<cstddef>`
themselves. The `std::` in front means it lives in the standard library's
namespace, the same as `std::vector`.

It differs from `int` in two ways. It is unsigned, so it holds only zero
and positive numbers. It is also as wide as a memory address, which is 64
bits on the Mac and on rangpur, so it can count up to about 1.8e19. A plain
`int` is 32 bits and stops at about 2.1 billion. Every standard container
uses it for sizes, so `value.size()` returns a `std::size_t`. Python never
needs this type because a Python `int` grows as big as it has to.

Our own code never declares a `std::size_t`. It only meets one when it asks
a vector for its size, as in `(int)requests.size()` in the MPI halo
exchange (§9c), where `(int)` turns it back into a plain `int`.

---

## 3. The entry point (`main.cpp`)

**3.1 The config file comes first**

```cpp
if (argc < 3 || std::string(argv[1]) != "--config") {
    usage();
    return 1;
}
Config cfg = load_config(argv[2]);
```

This block checks the command starts with `--config <file>` and loads that
file. `argv` is the list of words typed on the command line, like
`sys.argv` in Python, and `argc` is how many there are. `argv[0]` is the
program name, so the config path is `argv[2]`. Loading the file before any
other flag means every flag after it can overwrite what the file set. Every
script and Makefile target already puts `--config` first.

**3.2 The other flags, two words at a time**

```cpp
for (int i = 3; i < argc; i = i + 2) {
    std::string flag = argv[i];
    if (i + 1 >= argc) {
        throw std::runtime_error(flag + " needs a value");
    }
    std::string value = argv[i + 1];

    if (flag == "--solver") {
        cfg.solver = value;
    } else if (flag == "--type") {
        if (value != "call" && value != "put") {
            throw std::runtime_error("--type must be call or put");
        }
        cfg.option.is_call = (value == "call");
    } else if (flag == "--ns") {
        cfg.grid.num_stock_nodes = std::stoi(value);
    ...
    } else {
        usage();
        throw std::runtime_error("unknown flag: " + flag);
    }
}
```

This block reads the rest of the command line in pairs, such as `--ns`
then `1024`. Every flag takes exactly one value, so the loop moves on two
words at a time with `i = i + 2`. `std::stoi` turns the text `"1024"` into
the int 1024, like `int("1024")` in Python, and `std::stod` does the same
for decimals. A flag typed last with no value after it, or a flag the chain
does not know, stops the program. A typo such as `--nss` would otherwise be
ignored and you would benchmark the wrong grid.

`--type` is the one flag whose value is checked by hand, because it is text
and `std::stoi` cannot catch a wrong word. `(value == "call")` is a `bool`,
true for a call and false for a put, so the one line stores the option
type. `--vmax` and `--maturity` exist for the benchmarks. The weak scaling
sweep grows the variance axis with `--vmax`, and a shorter `--maturity`
keeps a big grid stable without a huge `--nt`.

`throw std::runtime_error(...)` stops the function at that line and jumps
to the `catch` at the bottom of `run` (§3.5), like `raise` in Python.

**3.2b The usage message** (`usage`)

```cpp
static void usage() {
    std::fprintf(stderr, "usage: heston --config <path> [--solver baseline|serial|simd|omp|mpi]\n");
    std::fprintf(stderr, "              [--type call|put] [--dump-every N] [--dump-dir DIR]\n");
    std::fprintf(stderr, "              [--bench R] [--ns X] [--nv X] [--nt X] [--maturity Y] [--vmax V]\n");
}
```

This block prints how to run the program. It writes to `stderr`, the error
stream, so `stdout` only ever holds result lines for the scripts to read.
`static` means only `main.cpp` can call it.

**3.3 Why not a fancier parser**

Libraries such as `getopt` let flags come in any order and handle flags
without values. Every flag here has a value and `--config` always comes
first, so a plain loop and an `if` chain do the job and read top to bottom.
`DESIGN_DECISIONS.md`, decisions 2 and 3.

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

**4.1b The market, the model and the grid** (`params.h`)

```cpp
struct MarketSpec {
    double spot = 5200.0;
    double rate = 0.045;
    double div_yield = 0.013;
};

struct HestonParams {
    double v0 = 0.04;
    double kappa = 1.5;
    double theta = 0.04;
    double xi = 0.35;
    double rho = -0.70;
};

struct GridSpec {
    int num_stock_nodes = 2048;
    int num_variance_nodes = 512;
    int num_timesteps = 2000;
    double stock_max_multiplier = 4.0;
    double variance_max = 1.0;
};
```

This block holds the other three groups of settings. `MarketSpec` is the
market today, the stock price `spot`, the interest rate and the dividend
yield. `HestonParams` is the five numbers of the Heston model, the variance
today `v0`, how fast it pulls back (`kappa`) towards its long run level
(`theta`), how much it jumps around (`xi`) and how it moves with the stock
(`rho`). STUDY_GUIDE §4 explains each one. `GridSpec` is the grid size and
how far each axis reaches. The stock axis runs from 0 to 4 times the
strike, 21000 for a 5250 strike, and the variance axis from 0 to 1.0, which
is 100 percent volatility.

Every field has a default written next to it. `Config cfg;` therefore
gives a complete, runnable set of settings, which is how the tests build
their settings in code with no config file (§12). The grid defaults have
the size of the 2048 by 512 by 2000 benchmark grid, and
`config/bench.cfg` sets `v_max` to 0.64 on top of them.

**4.2 Matching a key to a field** (`apply`)

```cpp
if (key == "option.strike") {
    cfg.option.strike = std::stod(value);
} else if (key == "grid.ns") {
    cfg.grid.num_stock_nodes = std::stoi(value);
...
} else {
    throw std::runtime_error(where + ": unknown key '" + key + "'");
}
```

This block stores one `key = value` line in the matching field. It is a plain
`if` / `else if` chain so it reads top to bottom, and an unknown key is an
error for the same reason as an unknown flag in §3.2. `static` in front of
`apply` means only `params.cpp` can call it, the way a leading underscore
marks a helper as private in Python.

**4.3 Reading the file** (`load_config`)

```cpp
std::ifstream file(path);
if (!file) {
    throw std::runtime_error("cannot open config: " + path);
}

Config cfg;
std::string line;
int line_number = 0;
while (std::getline(file, line)) {
    line_number = line_number + 1;

    std::istringstream words(line);
    std::string key;
    std::string equals;
    std::string value;
    words >> key >> equals >> value;

    if (key.empty() || key[0] == '#') {
        continue;
    }

    std::string where = path + ":" + std::to_string(line_number);
    if (equals != "=" || value.empty()) {
        throw std::runtime_error(where + ": expected key = value with spaces around the =");
    }

    try {
        apply(cfg, key, value, where);
    } catch (const std::invalid_argument&) {
        throw std::runtime_error(where + ": not a number '" + value + "'");
    }
}
return cfg;
```

This block reads the config one line at a time and splits each line into
words. `std::getline` reads one line, like `for line in file` in Python.
`std::istringstream words(line)` treats the line as a small stream of words,
and `words >> key >> equals >> value` takes the first three, like
`key, equals, value = line.split()[:3]` in Python. A blank line leaves `key`
empty, and a comment line has `#` as the first character of its first word,
so both are skipped. Each config line is written `grid.ns = 2048`, so the
middle word must be `=`.

`std::stod` and `std::stoi` throw `std::invalid_argument` when the text is
not a number at all, such as `abc`. The `try` and `catch` turn that into a
message with the file name and line number. Text that starts with a number,
such as `12x`, is read as 12 and the rest is ignored, which is the one check
the simpler code gave up.

`Config cfg;` starts from the defaults written in `params.h`, so a key left
out of the file keeps its default. The file closes itself when `file` goes
out of scope at the end of the function, which is RAII doing the job of
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

int total_cells = num_stock_nodes_ * num_variance_nodes_;
current_.assign(total_cells, 0.0);
next_.assign(total_cells, 0.0);
```

This block works out the spacing and makes both sheets, filled with zeros.
The spacing divides by `N - 1` because N nodes have N − 1 gaps between them.
`assign(n, 0.0)` makes the vector hold `n` zeros, like `np.zeros(n)`. An
`int` stops at about 2.1 billion, and the largest grid we run, 8192 by 1024,
has 8.4 million cells, so a plain `int` is enough. `index()` returns an
`int` too, so the two agree. Each sheet is one flat block of doubles rather than a list of rows, which
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
double position = coord / spacing;
if (position < 0.0) {
    position = 0.0;
}
if (position > num_nodes - 1) {
    position = num_nodes - 1;
}

int lower_node = (int)std::floor(position);
double nearest_node = std::round(position);
bool on_node = std::fabs(position - nearest_node) < node_snap_tolerance;
if (on_node) {
    lower_node = (int)nearest_node;
}

int highest_cell = num_nodes - 3;
if (highest_cell < 1) {
    highest_cell = 1;
}
AxisPosition out;
out.cell = lower_node;
if (out.cell < 1) {
    out.cell = 1;
}
if (out.cell > highest_cell) {
    out.cell = highest_cell;
}

if (on_node && out.cell == lower_node) {
    out.weight = 0.0;
} else {
    out.weight = position - out.cell;
    if (out.weight < 0.0) {
        out.weight = 0.0;
    }
    if (out.weight > 1.0) {
        out.weight = 1.0;
    }
}
return out;
```

This block finds the cell just below today's spot (or v0), and how far across
that cell the point sits. Each limit is a plain pair of `if` statements, the
same as `max(lo, min(x, hi))` in Python. `std::floor` and `std::round` are
the C `floor` and `round` from `<cmath>`. The weight is snapped to exactly 0 when the point
is within 1e-9 of a node, so an aligned grid always returns the node's value
unchanged. The cell is kept at least one node in from each edge, because the
Greeks read one cell further out on each side (STUDY_GUIDE §10).

The block runs in four parts, in this order.

1. `coord / spacing` turns a price or a variance into a position counted in
   cells. With a spacing of 50, a spot of 5200 is position 104.0, and 5225
   would be 104.5. The two `if` statements keep it between 0 and the last
   node.
2. `lower_node` is the node at or below the position. `std::floor` of
   103.9999999999 is 103, so a value that should be exactly on node 104 but
   picked up rounding error would land one cell too low. The snap with
   `nearest_node` fixes that.
3. `out.cell` is moved in to between 1 and `num_nodes - 3`. The Greeks read
   the cell to each side and the one above (§6.6), so cell 0 or a cell at
   the top would read off the grid.
4. `out.weight` is how far past `out.cell` the point sits, from 0 to 1. When
   step 3 moved the cell, the point is no longer inside it, so the weight
   is held between 0 and 1.

```text
  stock axis, spacing 50
  node:      103        104        105
             |----------|----------|
  5200 sits exactly on node 104     cell 104, weight 0.0
  5225 sits halfway to node 105     cell 104, weight 0.5
```

**5.7 The two public lookups** (`Grid::stock_position`, `Grid::variance_position`)

```cpp
AxisPosition Grid::stock_position(double spot) const {
    return locate(spot, stock_spacing_, num_stock_nodes_);
}

AxisPosition Grid::variance_position(double v0) const {
    return locate(v0, variance_spacing_, num_variance_nodes_);
}
```

This block gives each axis its own call into `locate`, passing that axis's
spacing and node count. `locate` is `static` in `grid.cpp`, so code outside
the file cannot call it, and these two methods are the way in. The answer
comes back as an `AxisPosition` ([grid.h:14](include/grid.h#L14)), a struct
with two fields, `cell` and `weight`, like a Python named tuple.

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
    double exchange_seconds = 0.0;
    double boundary_seconds = 0.0;
};
```

This block is the record every solver hands back. Returning one struct keeps
`main`, the CSV printer and the tests the same for every solver. `threads`
and `ranks` say how many OpenMP threads and MPI processes did the work. They
default to 1, so only the parallel rungs have to set them. `exchange_seconds` is
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
double max_stock_price = grid.stock_price(grid.num_stock_nodes() - 1);
double max_variance = grid.variance(grid.num_variance_nodes() - 1);
double stock_spacing = grid.stock_spacing();
double variance_spacing = grid.variance_spacing();
out.dt_stable_estimate = 1.0 / (max_variance * max_stock_price * max_stock_price / (stock_spacing * stock_spacing) + xi * xi * max_variance / (variance_spacing * variance_spacing) + rate);

bool feller_ok = 2.0 * kappa * theta >= xi * xi;
const char* feller_text = "VIOLATED";
if (feller_ok) {
    feller_text = "OK";
}
std::fprintf(stderr, "[%s] feller %s (2*kappa*theta=%.4f, xi^2=%.4f)  dt=%.3e  dt_stable~%.3e\n",
             name().c_str(), feller_text, 2.0 * kappa * theta, xi * xi, dt, out.dt_stable_estimate);
```

This block works out the largest stable timestep and prints it next to the
one we are using, along with whether the Feller condition holds. The formula
is the point where the centre weight of the update reaches zero at the far
corner of the grid (STUDY_GUIDE §7.5). It prints to stderr so stdout stays
one clean CSV line for the scripts.

The far corner is the largest stock price and the largest variance, where
the diffusion terms are biggest, so it is the first cell to go unstable.
The estimate is stored in `out` as well as printed, because the tests read
it to refuse an unstable run (§12.2). `const char*` is a C string, plain
text that is never changed. It starts as `"VIOLATED"` and becomes `"OK"`
when Feller holds. Feller is printed for the write-up only. The solver runs
the same either way, because the v = 0 row (§8.7) handles a variance that
reaches zero.

**6.6 Reading the answer** (`Solver::extract_result`)

```cpp
const std::vector<double>& values = grid.current();

AxisPosition spot_position = grid.stock_position(cfg.market.spot);
AxisPosition v0_position = grid.variance_position(cfg.heston.v0);
int base_i = spot_position.cell;
int base_j = v0_position.cell;
double stock_weight = spot_position.weight;
double variance_weight = v0_position.weight;

if (stock_weight > off_node_tolerance || variance_weight > off_node_tolerance) {
    std::fprintf(stderr, ...);
}

out.price = bilinear_interpolate(grid, values, base_i, base_j, stock_weight, variance_weight);
double east = bilinear_interpolate(grid, values, base_i + 1, base_j, stock_weight, variance_weight);
double west = bilinear_interpolate(grid, values, base_i - 1, base_j, stock_weight, variance_weight);
double north = bilinear_interpolate(grid, values, base_i, base_j + 1, stock_weight, variance_weight);
double south = bilinear_interpolate(grid, values, base_i, base_j - 1, stock_weight, variance_weight);

out.delta = (east - west) / (2.0 * grid.stock_spacing());
out.gamma = (east - 2.0 * out.price + west) / (grid.stock_spacing() * grid.stock_spacing());
out.vega = (north - south) / (2.0 * grid.variance_spacing());
```

This block reads the price at today's spot and variance, and gets the Greeks
from the same finite differences as STUDY_GUIDE §7. Every value goes through
the same blend, so the reported delta is the slope of the reported price and
not of some other surface.

The block runs in four parts.

1. `values` is the current buffer. After the last swap it holds the option
   value today, at every stock price and variance on the grid.
2. `stock_position` and `variance_position` (§5.7) turn today's spot and
   `v0` into a cell and a weight on each axis. `base_i` and `base_j` are the
   cell's column and row.
3. When either weight is above 1e-9, today's market sits between nodes. A
   warning goes to stderr, because the price is then interpolated and
   carries a little extra error. The collapse and parity tests use a grid
   where both land on nodes, so they compare node values. The 2048 by 512
   benchmark grid is off node, which does not matter for timing.
4. The price is read at the point, and the same interpolation is done one
   cell east, west, north and south. Those four values feed the central
   differences, delta and gamma along the stock axis and vega along the
   variance axis.

```text
              north (base_j + 1)
                    |
  west (base_i - 1) - price - east (base_i + 1)
                    |
              south (base_j - 1)
```

**6.7 The blend** (`bilinear_interpolate`)

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
SolveResult result;

Grid grid(cfg.grid, cfg.option);
grid.init_payoff(cfg.option);

double strike = cfg.option.strike;
double rate = cfg.market.rate;
double div_yield = cfg.market.div_yield;
double kappa = cfg.heston.kappa;
double theta = cfg.heston.theta;
double xi = cfg.heston.xi;
double rho = cfg.heston.rho;
int num_timesteps = cfg.grid.num_timesteps;

double dt = cfg.option.maturity_years / num_timesteps;
report_stability(grid, cfg, dt, result);

double max_stock_price = grid.stock_price(grid.num_stock_nodes() - 1);
```

This block makes the two sheets, fills the first one with the payoff and
works out the timestep. `result` starts with every field at its default
from `solver.h`, and the solve fills it in as it goes. The short local
names such as `rate` and `kappa` make the update in §7.6 read like the
maths. `dt` is the maturity split into `num_timesteps` equal steps, for
example 0.25 years over 2000 steps. `max_stock_price` is the last stock
node, which the call's top edge needs in §7.8.

The stability report prints before any maths runs because a `dt` above the limit makes the whole run blow up (STUDY_GUIDE §8),
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
    double step_start = omp_get_wtime();
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

double step_end = omp_get_wtime();
result.seconds += step_end - step_start;

if (cfg.dump_every > 0 && step % cfg.dump_every == 0) {
    dump_snapshot(cfg.dump_dir, step, grid);
}
```

This block makes the new sheet the current one, adds this step's time to the
total and sometimes writes the sheet to a file. `omp_get_wtime()` returns
the wall clock time in seconds as a `double`, so the time for a step is one
subtraction. It comes from OpenMP, and every solver uses it, the serial
ones included, so all of them are timed the same way. The timer stops before the
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

**8.1b The scalars** (`StencilConstants` constructor, [stencil.cpp:12](src/stencil.cpp#L12))

```cpp
StencilConstants::StencilConstants(const Grid& grid, const Config& cfg, double timestep) {
    num_stock_nodes = grid.num_stock_nodes();
    num_variance_nodes = grid.num_variance_nodes();

    dt = timestep;
    rate = cfg.market.rate;
    div_yield = cfg.market.div_yield;
    stock_drift = cfg.market.rate - cfg.market.div_yield;

    is_call = cfg.option.is_call;
    strike = cfg.option.strike;
    max_stock_price = grid.stock_price(grid.num_stock_nodes() - 1);
```

This block is the start of the constructor, which runs once when a solver
writes `StencilConstants constants(grid, cfg, dt);`. It copies the grid
size and the market numbers into the struct, so the row functions need
only `constants` and the two buffers. `stock_drift` is r − q, worked out
once because the update multiplies by it in every cell. `is_call`,
`strike` and `max_stock_price` are there for `apply_stock_edges` (§8.8).
Each field is set with a plain assignment in the body, like `self.x = x` in
a Python `__init__`.

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
stock_price.resize(num_stock_nodes);
stock_price_sq.resize(num_stock_nodes);
for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
    double price = stock_i * stock_spacing;
    stock_price[stock_i] = price;
    stock_price_sq[stock_i] = price * price;
}

double kappa = cfg.heston.kappa;
double theta = cfg.heston.theta;
double xi = cfg.heston.xi;
double rho = cfg.heston.rho;
stock_diffusion.resize(num_variance_nodes);
...
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

`resize(n)` makes a vector hold `n` zeros so the loop can fill it by index,
like `[0.0] * n` in Python. The stock table has one entry per column
(2048) and the coefficient tables one per row (512). The brackets in
`(rho * xi) * variance` show the order C++ multiplies in anyway, left to
right, which is the same order as `rho * xi * variance` in the baseline.

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

    double zero_stock_value = 0.0;
    double max_stock_value = 0.0;
    if (constants.is_call) {
        zero_stock_value = 0.0;
        max_stock_value = constants.max_stock_price * std::exp(-constants.div_yield * time_to_expiry) - strike_today;
    } else {
        zero_stock_value = strike_today;
        max_stock_value = 0.0;
    }

    int last_stock = constants.num_stock_nodes - 1;
    for (int var_j = row_begin; var_j < row_end; ++var_j) {
        next[var_j * constants.num_stock_nodes] = zero_stock_value;
        next[var_j * constants.num_stock_nodes + last_stock] = max_stock_value;
    }
}
```

This block is §7.8, but it only writes rows `row_begin` to `row_end - 1`. The
serial rung passes the whole grid. The row range is there for the MPI rung,
where each process only owns some of the rows.

The two values come from what the option is worth at each end of the stock
axis, with `time_to_expiry` years left (STUDY_GUIDE §9).

- `strike_today` is the strike K paid at expiry, discounted back to now,
  K·e^(−rτ). `std::exp` is e to the power, like `math.exp`.
- A call at S = 0 is worth 0, because the stock can never climb back above
  the strike.
- A call at the top of the axis, 4 times the strike, is almost sure to be
  exercised, so it is worth what a forward contract is worth, the stock
  less its dividends minus the discounted strike, S·e^(−qτ) − K·e^(−rτ).
- A put at S = 0 is sure to pay the full strike, so it is worth
  K·e^(−rτ).
- A put at the top is worth 0.

The values do not depend on the variance, so one pair is worked out and
written into every row. `var_j * num_stock_nodes` is the first cell of row
`var_j`, and adding `last_stock` gives the last cell of the same row.

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

The pointers `row_centre`, `row_above`, `row_below` and `out_row` are made
just before the loop at
[stencil_avx.cpp:44](src/stencil_avx.cpp#L44), the same way as in the
scalar row (§8.5).

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
```

This block is the five central differences from §8.6, one intrinsic per
arithmetic operation. `_mm256_sub_pd(a, b)` is `a - b` in all four lanes,
`_mm256_mul_pd` is `a * b`, and `_mm256_fnmadd_pd(a, b, c)` is `c - a * b`
done as one instruction with one rounding. We name every intermediate so
each line does one thing, which keeps the AVX version readable next to the
scalar version.

Each group matches one line of §8.6. `east_minus_west` times `inv_two_ds`
is (east − west)/(2·ds). `stock_curvature` is east + west − 2·centre, the
top of the second difference, and times `inv_ds_sq` it is d²V/dS². The
variance lines do the same with north and south. `corners` builds
NE − NW − SE + SW one operation at a time, because an intrinsic takes two
inputs.

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
double boundary_start = omp_get_wtime();
update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
grid.swap_buffers();
double boundary_end = omp_get_wtime();
result.boundary_seconds += boundary_end - boundary_start;
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
possibly on different machines. Each rank owns one block of variance rows
(`owned_rows`, called a slab in STUDY_GUIDE) and runs the OpenMP loop over its
own rows only. The rows received from the neighbours are `row_below` and
`row_above`, called halo or ghost rows in HPC texts. STUDY_GUIDE §17 explains MPI,
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

**9c.2 Sharing out the rows** (`rows_for_rank`)

```cpp
struct RowRange {
    int begin = 0;
    int end = 0;
};

static RowRange rows_for_rank(int num_rows, int rank, int num_ranks) {
    int rows_per_rank = num_rows / num_ranks;
    int extra_rows = num_rows % num_ranks;
    RowRange rows;
    if (rank < extra_rows) {
        rows.begin = rank * (rows_per_rank + 1);
        rows.end = rows.begin + rows_per_rank + 1;
    } else {
        rows.begin = extra_rows * (rows_per_rank + 1) + (rank - extra_rows) * rows_per_rank;
        rows.end = rows.begin + rows_per_rank;
    }
    return rows;
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

Grid grid(cfg.grid, cfg.option);
grid.init_payoff(cfg.option);
...
if (rank == 0) {
    report_stability(grid, cfg, dt, result);
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

Every rank builds the whole grid and the whole payoff, the same as the
other rungs. Only rank 0 prints the stability line, so a run on 4 ranks
prints it once.

**9c.4 This rank's part of each step**

```cpp
int first_interior = std::max(owned_rows.begin, 1);
int end_interior = std::min(owned_rows.end, num_variance_nodes - 1);
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
apply_stock_edges(constants, next, owned_rows.begin, owned_rows.end, step * dt);
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

**9c.5 The halo exchange** (`exchange_edge_rows`)

```cpp
double* first_row = buffer.data() + owned_rows.begin * num_stock_nodes;
double* last_row = buffer.data() + (owned_rows.end - 1) * num_stock_nodes;
double* row_below = buffer.data() + (owned_rows.begin - 1) * num_stock_nodes;
double* row_above = buffer.data() + owned_rows.end * num_stock_nodes;
std::vector<MPI_Request> requests;

if (rank > 0) {
    MPI_Request receive_below;
    MPI_Irecv(row_below, num_stock_nodes, MPI_DOUBLE, rank - 1, message_tag, MPI_COMM_WORLD, &receive_below);
    requests.push_back(receive_below);
    MPI_Request send_below;
    MPI_Isend(first_row, num_stock_nodes, MPI_DOUBLE, rank - 1, message_tag, MPI_COMM_WORLD, &send_below);
    requests.push_back(send_below);
}
if (rank < num_ranks - 1) {
    ... the same with rank + 1, row_above and last_row ...
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

  row 11  first row ── send ──────►  row 11  row_above
  row 10  row_below ◄────── send ── row 10  last row
```

**9c.6 Where the exchange happens, and the timer**

```cpp
double exchange_start = omp_get_wtime();
exchange_edge_rows(next, num_stock_nodes, owned_rows, rank, num_ranks);
double exchange_end = omp_get_wtime();
result.exchange_seconds += exchange_end - exchange_start;

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
    RowRange other_rows = rows_for_rank(num_variance_nodes, other_rank, num_ranks);
    counts[other_rank] = (other_rows.end - other_rows.begin) * num_stock_nodes;
    offsets[other_rank] = other_rows.begin * num_stock_nodes;
}
const double* my_values = current.data() + owned_rows.begin * num_stock_nodes;
int my_count = (owned_rows.end - owned_rows.begin) * num_stock_nodes;
MPI_Gatherv(my_values, my_count, MPI_DOUBLE, next.data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

**9d.1 The scalar row on every core** (`OmpScalarSolver`)

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_scalar(constants, current, next, var_j);
}
```

This block is the `omp` loop with `update_row_avx` changed back to
`update_row_scalar`. It runs the threads without SIMD, so the strong
scaling figure can show the thread speedup on its own next to the thread
speedup with SIMD. Comparing the two lines shows whether the two speedups
multiply or get in each other's way, for example when the memory bus is
full.

**9d.2 One parallel region for the whole loop** (`OmpPersistentSolver`)

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

```cpp
double loop_start = omp_get_wtime();
...
double loop_end = omp_get_wtime();
result.seconds = loop_end - loop_start;
```

This solver times the whole loop once, from outside the parallel region.
Inside the region every thread runs the step loop, so a per-step timer
would run once per thread. It has no `boundary_seconds` for the same
reason.

**9d.3 The schedule picked when the program runs** (`OmpRuntimeScheduleSolver`)

```cpp
#pragma omp parallel for schedule(runtime) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
```

This block is the `omp` loop with `schedule(static)` changed to
`schedule(runtime)`. The schedule then comes from the environment variable
`OMP_SCHEDULE`, so one build can try `dynamic,1`, `dynamic,8` or `guided`
without recompiling. `dynamic` hands out rows a few at a time as threads
finish, which helps when rows cost different amounts. Every row here costs
the same, and experiment E5 measured `dynamic,1` at 0.60 times the speed of
`static`, `dynamic,8` at 0.72 and `guided` at 0.84. The rerun gave 0.62,
0.73 and 0.77.

**9d.4 Boundaries inside the loop** (`OmpParallelEdgesSolver`)

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
    stock_price.resize(num_stock_nodes);
    stock_price_sq.resize(num_stock_nodes);
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        stock_price[stock_i] = (float)constants.stock_price[stock_i];
        stock_price_sq[stock_i] = (float)constants.stock_price_sq[stock_i];
    }
    ...
}
```

This block builds the float constants from the double ones. `(float)x`
rounds a double to the nearest float, and the loops copy each table one
entry at a time, rounding every entry. We start from the double constants so both solvers
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
int total_cells = num_stock_nodes * num_variance_nodes;
std::vector<float> current(total_cells, 0.0f);
std::vector<float> next(total_cells, 0.0f);
for (int cell = 0; cell < total_cells; ++cell) {
    current[cell] = (float)grid.current()[cell];
}
...
std::swap(current, next);
...
std::vector<double>& answer = grid.next();
for (int cell = 0; cell < total_cells; ++cell) {
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

`update_zero_row_float` at [solver_float.cpp:18](src/solver_float.cpp#L18)
is `update_zero_row` (§8.7) line for line in float.

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
            result.threads, result.ranks, result.exchange_seconds, result.boundary_seconds);
```

This block prints one CSV line per solve with the price, the Greeks, the grid
size, the solver name, the timing, and the number of threads and processes.
The benchmark and plotting scripts parse these columns by position, so new
columns only ever go on the end, which is where `threads`, `ranks`,
`exchange_seconds` and `boundary_seconds` went.

**10.2 Snapshots** (`dump_snapshot`)

```cpp
std::filesystem::create_directories(dump_dir);
char file_name[32];
std::snprintf(file_name, sizeof file_name, "snap_%06d.csv", step);
std::string path = dump_dir + "/" + file_name;
FILE* out = std::fopen(path.c_str(), "w");
if (out == NULL) {
    std::fprintf(stderr, "dump_snapshot: cannot write %s\n", path.c_str());
    return;
}

for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
    for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
        if (stock_i > 0) {
            std::fprintf(out, ",");
        }
        std::fprintf(out, "%.10g", grid.current()[grid.index(stock_i, var_j)]);
    }
    std::fprintf(out, "\n");
}
std::fclose(out);
```

This block writes the whole current sheet to a CSV file, one variance row per
line. It uses C's `fopen`, `fprintf` and `fclose`, the same calls as in C.
`path.c_str()` turns the C++ string into the C string that `fopen` wants,
and `%.10g` writes each value with 10 significant digits. `fopen` returns
`NULL` when it cannot open the file, for example in a folder with no write
permission. The function then prints a message and returns, so a missing
snapshot never stops a solve. `snprintf` writes the file name into the
32-byte buffer `file_name`, and `sizeof file_name` stops it writing past
the end.
`std::filesystem::create_directories` is the one C++17 call here. It is
like `os.makedirs(dump_dir, exist_ok=True)` and makes every missing folder
on the path, which C's `mkdir` cannot do in one call. The animation scripts read these files, and the step number is
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

**12.0 How each test is built**

```cpp
static Config aligned_test_config() {
    Config cfg;
    cfg.grid.num_stock_nodes = 421;
    cfg.grid.num_variance_nodes = 51;
    cfg.grid.num_timesteps = 56000;
    return cfg;
}
...
BaselineSolver solver;
SolveResult result = solver.solve(cfg);
```

This block is the pattern every test follows. It starts from the defaults
in `params.h` (§4.1b), changes only the grid, runs a solver, prints one
line per check and returns 0 or 1. No config file is read, so a test
cannot pass or fail because a file changed.

The numbers are picked so the answer lands on nodes. 421 stock nodes over
0 to 21000 gives a spacing of 21000 / 420 = 50, so the spot 5200 is node
104 and the strike 5250 is node 105. 51 variance nodes over 0 to 1 gives a
spacing of 0.02, so `v0` = 0.04 is node 2. The weights are then 0 and the
tests compare grid values directly (§6.7). 56000 steps over 0.25 years
gives `dt` = 4.5e-6, about 0.79 of the stability limit for this grid.

`BaselineSolver solver;` makes the solver as a normal local variable on the
stack. The test knows which solver it wants, so it needs no
`make_solver` and no `unique_ptr`, and the solver is freed at the end of
`main`.

`test_ladder_matches` uses a smaller 128 by 32 grid with 5200 steps
([test_ladder_matches.cpp:33](tests/test_ladder_matches.cpp#L33)), because
it runs 16 solves for the call and 16 for the put, and still finishes in
about 7 seconds on the Mac. That grid is off node, so
it also tests the interpolation in §6.7.

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
for (const std::string& solver_name : solver_names) {
    std::unique_ptr<Solver> solver = make_solver(solver_name);
    SolveResult candidate = solver->solve(cfg);
    ...
    double worst_diff = max_relative_diff(reference, candidate);
    bool ok = worst_diff < relative_tolerance;
```

The `for (... : solver_names)` loop is C++'s version of
`for solver_name in solver_names` in Python. The check is written as
`worst_diff < relative_tolerance` rather than `>=` for the fail case, because
a run that blew up gives NaN, and any comparison with NaN is false, so a NaN
fails.

This block runs the baseline and every other rung on the same small grid and
compares price, delta, gamma and vega. A rung is not benchmarked until this
passes, because a fast wrong answer is worth nothing. Each new rung is added
to `solver_names`.

```cpp
static const std::vector<std::string> same_bits_as_simd = {"omp", "mpi", "omp_persistent", "omp_runtime", "omp_edges"};
...
for (const std::string& solver_name : same_bits_as_simd) {
    SolveResult candidate = make_solver(solver_name)->solve(cfg);
    ...
    double simd_diff = max_relative_diff(simd_reference, candidate);
    const char* verdict = "OK";
    if (simd_diff != 0.0) {
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
static const double float_tolerance = 1e-4;
...
SolveResult float_candidate = make_solver("omp_float")->solve(cfg);
...
double float_diff = max_relative_diff(reference, float_candidate);
bool ok = float_diff < float_tolerance;
```

This block holds the float experiment (§9e) to its own limit of 1e-4. A
float keeps about 7 digits, so it can never meet 1e-11, and it is not put
in `solver_names`. On this grid it measures 2.1e-6 on the call and 8.2e-6 on the
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

```cpp
static double max_relative_diff(const SolveResult& a, const SolveResult& b) {
    double worst = relative_diff(a.price, b.price);
    worst = std::fmax(worst, relative_diff(a.delta, b.delta));
    worst = std::fmax(worst, relative_diff(a.gamma, b.gamma));
    worst = std::fmax(worst, relative_diff(a.vega, b.vega));
    return worst;
}
```

This block turns two results into one number, the worst relative
difference across the price and the three Greeks. `std::fmax` is `max` for
doubles. One number per solver keeps each check to one `<` and one printed
line.

**12.4 The AVX2 row against the scalar row** (`test_avx_row.cpp`)

```cpp
if (!avx_row_is_vectorised()) {
    std::printf("test_avx_row: FAIL, built without AVX2 and FMA\n");
    return 1;
}
...
std::mt19937 random_engine(3500);
std::vector<int> row_lengths = {130, 131, 132, 133, 2048};
...
std::uniform_real_distribution<double> noise(-1.0, 1.0);
for (int cell = 0; cell < num_cells; ++cell) {
    int stock_i = cell % num_stock_nodes;
    double payoff = std::fmax(constants.stock_price[stock_i] - cfg.option.strike, 0.0);
    current[cell] = payoff + noise(random_engine);
}

std::vector<double> next_scalar(num_cells, -7.0);
std::vector<double> next_avx(num_cells, -7.0);
for (int var_j = 1; var_j < cfg.grid.num_variance_nodes - 1; ++var_j) {
    update_row_scalar(constants, current, next_scalar, var_j);
    update_row_avx(constants, current, next_avx, var_j);
}
```

This block fills a sheet with a call payoff plus random noise, runs both row
functions on it, and compares every cell of the two outputs. A plain payoff
is a straight line or flat in most places, so most derivatives would be
zero and a wrong term could hide. The noise makes every derivative non
zero. `std::mt19937` is a random number generator, and the fixed seed 3500
makes every run use the same numbers, like `random.seed(3500)` in Python.
`cell % num_stock_nodes` is the column of a cell in the flat buffer. The row lengths
130 to 133 leave zero, one, two and three cells for the scalar remainder, so
every path through §9a.3 is tested. Both outputs start filled with -7, so if
the AVX row wrote outside the interior the edges would differ and the test
would fail. It refuses to pass when the AVX2 row was not built, because then
it would only be comparing the scalar row with itself. It runs with
`make avx-check` (§13.5), not `make test`, because on the Mac it needs a
separate x86_64 build.

```cpp
static const double float_tolerance = 1e-5;
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

- Starting from `main`, name every function a serial run calls, in order.
  (§1.5)
- Which lines of `solver_omp.cpp` differ from `solver_serial.cpp`? (§1.5)
- When `main` calls `solver->solve(cfg)`, how does the program know which
  `solve` to run? (§1.5, §6.2)
- Why must `--config` come first on the command line? (§3.1)
- Why do the tests build a `Config` in code and not read a file? (§4.1b,
  §12.0)
- With a stock spacing of 50, what cell and weight does `locate` give for
  5225? (§5.6)
- What is the call worth at the top of the stock axis, and why? (§8.8)
- Why are the test grids 421 by 51? (§12.0)
- What does `words >> key >> equals >> value` do with a blank line? (§4.3)
- What is `std::size_t`, and where does the code meet one? (§2.2)
- How does a solver time one step, and why does every solver use
  `omp_get_wtime`? (§7.10)
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
  they? (§9d.2)
- Why does `omp_persistent` time the whole loop once and not each step?
  (§9d.2)
- What does `schedule(runtime)` read, and what did E5 measure? (§9d.3)
- Why can row 0 and the stock edges join the parallel loop, but the top row
  cannot? (§9d.4)

---

## 15. Before and after, one change at a time

This section takes every code change in [CHANGELOG.md](CHANGELOG.md), in the
same order, and puts the code before the change next to the code after it.
Each pair is followed by what changed and why. Every After snippet is copied
from the real file with most comments removed. Each Before snippet comes from
the `milestone1` tag or from the rung below. Two of them were never
committed on their own (§15.5 and §15.12), so they are rebuilt from the
changelog's description and say so.

Phases 6, 8 and 9 only changed figures, documents and slides, so they have
no entry here.

### 15.1 Names that read like the PDE (Phase 1, Change 1)

**Before** (`solver_baseline.cpp` at the `milestone1` tag)

```cpp
double S = g.stock_price(stock_i);
double V = cur[g.index(stock_i, var_j)];
double v = g.variance(var_j);
double ne = cur[g.index(stock_i + 1, var_j + 1)];
double V_S = (east - west) / (2.0 * g.stock_spacing());
double V_SS = (east - 2.0 * V + west) / (g.stock_spacing() * g.stock_spacing());
double pde = V + dt * (0.5 * v * S * S * V_SS + rho * xi * v * S * V_Sv + ...);
```

**After** (`solver_baseline.cpp`)

```cpp
double stock_price = grid.stock_price(stock_i);
double centre = current[grid.index(stock_i, var_j)];
double variance = grid.variance(var_j);
double north_east = current[grid.index(stock_i + 1, var_j + 1)];
double dV_dS = (east - west) / (2.0 * grid.stock_spacing());
double d2V_dS2 = (east - 2.0 * centre + west) / (grid.stock_spacing() * grid.stock_spacing());
double updated = centre + dt * (0.5 * variance * stock_price * stock_price * d2V_dS2 + ...);
```

Only the names changed. `V` and `v` differed by one capital letter, so the
cell value is now `centre` and the variance is `variance`. `ne` is spelled
out as `north_east`, and `V_S` became `dV_dS` so it reads like the
derivative it is. No arithmetic moved, so the baseline still gives exactly
the Milestone 1 answer.

### 15.2 The constants move out of the cell loop (Phase 1, Decisions 1 to 3)

**Before** (`solver_baseline.cpp`, inside the loop over every cell)

```cpp
double stock_price = grid.stock_price(stock_i);
double variance = grid.variance(var_j);
double centre = current[grid.index(stock_i, var_j)];
double east = current[grid.index(stock_i + 1, var_j)];
double dV_dS = (east - west) / (2.0 * grid.stock_spacing());
double d2V_dS2 = (east - 2.0 * centre + west) / (grid.stock_spacing() * grid.stock_spacing());
```

**After** (`stencil.cpp`, `StencilConstants` constructor, runs once)

```cpp
inv_two_ds = 1.0 / (2.0 * stock_spacing);
inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
    double price = stock_i * stock_spacing;
    stock_price[stock_i] = price;
    stock_price_sq[stock_i] = price * price;
}
for (int var_j = 0; var_j < num_variance_nodes; ++var_j) {
    double variance = var_j * variance_spacing;
    stock_diffusion[var_j] = 0.5 * variance;
    ...
}
```

**After** (`stencil.cpp`, `update_cells_scalar`, inside the loop over every cell)

```cpp
double dt = constants.dt;
double inv_two_ds = constants.inv_two_ds;
...
const double* row_centre = current.data() + var_j * num_stock_nodes;
const double* row_above = row_centre + num_stock_nodes;
const double* row_below = row_centre - num_stock_nodes;

for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {
    double centre = row_centre[stock_i];
    double east = row_centre[stock_i + 1];
    double dV_dS = (east - west) * inv_two_ds;
    double d2V_dS2 = (east - 2.0 * centre + west) * inv_ds_sq;
```

Three things moved out of the cell loop.

1. The baseline divides by the spacing five times per cell. Division is
   slow, so `StencilConstants` works out each `1 / spacing` once, and the
   loop multiplies by it.
2. The baseline works out S, S² and the four `v` coefficients again for
   every cell. They only depend on the column or the row, so they are now
   tables built once, and the loop looks them up.
3. The baseline calls `grid.index(i, j)` for each of the nine neighbours.
   The new code keeps a pointer to the start of the row, the row above and
   the row below, so a neighbour is just `row_above[stock_i + 1]`.

The constants are also copied into local variables at the top of the row
function. Without that the compiler would have to reload them after every
write through `out_row`, because it cannot prove the write did not change
them. Only hoisting and reciprocals moved the M1 number (53 to 131
Mcell-updates/s), so those are the only M1 techniques the serial rung keeps.

### 15.3 One long expression becomes five named terms (Phase 1, Change 2)

**Before** (`solver_baseline.cpp`, kept this way on purpose)

```cpp
double updated = centre + dt * (0.5 * variance * stock_price * stock_price * d2V_dS2
                                + rho * xi * variance * stock_price * d2V_dS_dv
                                + 0.5 * xi * xi * variance * d2V_dv2
                                + (rate - div_yield) * stock_price * dV_dS
                                + kappa * (theta - variance) * dV_dv
                                - rate * centre);
```

**After** (`stencil.cpp`, `update_cells_scalar`)

```cpp
double stock_diffusion_term = stock_diffusion * stock_price_sq[stock_i] * d2V_dS2;
double cross_term = cross_diffusion * stock_price[stock_i] * d2V_dS_dv;
double variance_diffusion_term = variance_diffusion * d2V_dv2;
double stock_drift_term = stock_drift * stock_price[stock_i] * dV_dS;
double variance_drift_term = variance_drift * dV_dv;

double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
out_row[stock_i] = centre + dt * change;
```

Each PDE term now has its own line and its own name, so you can point at one
term and say which part of the equation it is. The baseline was left as one
expression so its answer stays bit for bit the M1 answer. On the Mac this
split cost about 13%, because Apple clang only fuses a multiply and an add
into one FMA instruction inside a single statement. Rangpur's build has no
`-mfma` for this file, so it is not affected there.

### 15.4 The time loop becomes five function calls (Phase 1, Decision 2)

**Before** (`solver_baseline.cpp`, the body of one timestep, shortened)

```cpp
for (int var_j = 1; var_j < grid.num_variance_nodes() - 1; ++var_j) {
    for (int stock_i = 1; stock_i < grid.num_stock_nodes() - 1; ++stock_i) {
        ... about 30 lines of maths per cell ...
    }
}
for (int stock_i = 1; stock_i < grid.num_stock_nodes() - 1; ++stock_i) {
    ... the v = 0 row ...
}
for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
    next[grid.index(0, var_j)] = zero_stock_value;
    next[grid.index(grid.num_stock_nodes() - 1, var_j)] = max_stock_value;
}
for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
    next[grid.index(stock_i, top_row)] = next[grid.index(stock_i, top_row - 1)];
}
grid.swap_buffers();
```

**After** (`solver_serial.cpp`)

```cpp
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_scalar(constants, current, next, var_j);
}
update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
grid.swap_buffers();
```

The four loops are the same four jobs, each moved into a function in
`stencil.cpp`. The point is that every later rung changes only one of these
lines. Keep this snippet in mind, because §15.6 to §15.9 are all small edits
to it.

### 15.5 The scalar row takes a column range (Phase 2, Decision 3)

**Before** (rebuilt from the changelog, not committed on its own)

```cpp
void update_row_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    ...
    for (int stock_i = 1; stock_i < num_stock_nodes - 1; ++stock_i) {
        ... the cell update ...
    }
}
```

**After** (`stencil.cpp`)

```cpp
void update_cells_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j, int first_stock, int end_stock) {
    ...
    for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {
        ... the same cell update ...
    }
}

void update_row_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    update_cells_scalar(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}
```

The loop body did not change. It only gained a start and an end column. The
AVX2 row does four cells at a time, and when a row's interior is not a
multiple of four the last one to three cells are left over. The AVX2 row
hands those to `update_cells_scalar`, so they are done by exactly the serial
rung's code and not by a second copy of it.

### 15.6 One cell at a time becomes four cells at a time (Phase 2)

**Before** (`stencil.cpp`, `update_cells_scalar`)

```cpp
for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {
    double centre = row_centre[stock_i];
    double east = row_centre[stock_i + 1];
    double west = row_centre[stock_i - 1];

    double dV_dS = (east - west) * inv_two_ds;
    double d2V_dS2 = (east - 2.0 * centre + west) * inv_ds_sq;
    ...
    double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
    out_row[stock_i] = centre + dt * change;
}
```

**After** (`stencil_avx.cpp`, `update_row_avx`)

```cpp
__m256d inv_two_ds = _mm256_set1_pd(constants.inv_two_ds);
__m256d two = _mm256_set1_pd(2.0);

while (stock_i + 3 <= last_interior) {
    __m256d centre = _mm256_loadu_pd(row_centre + stock_i);
    __m256d east = _mm256_loadu_pd(row_centre + stock_i + 1);
    __m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);

    __m256d east_minus_west = _mm256_sub_pd(east, west);
    __m256d dV_dS = _mm256_mul_pd(east_minus_west, inv_two_ds);
    __m256d east_plus_west = _mm256_add_pd(east, west);
    __m256d stock_curvature = _mm256_fnmadd_pd(two, centre, east_plus_west);
    __m256d d2V_dS2 = _mm256_mul_pd(stock_curvature, inv_ds_sq);

    __m256d change = _mm256_mul_pd(stock_diffusion_factor, d2V_dS2);
    change = _mm256_fmadd_pd(cross_factor, d2V_dS_dv, change);
    ...
    change = _mm256_fnmadd_pd(rate, centre, change);
    __m256d new_value = _mm256_fmadd_pd(dt, change, centre);
    _mm256_storeu_pd(out_row + stock_i, new_value);

    stock_i = stock_i + 4;
}
update_cells_scalar(constants, current, next, var_j, stock_i, last_interior + 1);
```

The maths is the same. What changed is the size of each
value. A `double` holds one cell and a `__m256d` holds four cells side by
side, so every line now works on four cells. This table maps one to the
other.

| Scalar (one cell) | AVX2 (four cells) | What it does |
|---|---|---|
| `double x` | `__m256d x` | a value, or four values in one register |
| `double c = constants.c` | `_mm256_set1_pd(constants.c)` | the same constant in all four lanes |
| `row_centre[stock_i + 1]` | `_mm256_loadu_pd(row_centre + stock_i + 1)` | read four cells starting there |
| `a - b`, `a + b`, `a * b` | `_mm256_sub_pd`, `_mm256_add_pd`, `_mm256_mul_pd` | the same operation on all four lanes |
| `c + a * b` | `_mm256_fmadd_pd(a, b, c)` | multiply and add in one instruction, one rounding |
| `c - a * b` | `_mm256_fnmadd_pd(a, b, c)` | the same, with the product taken away |
| `out_row[stock_i] = x` | `_mm256_storeu_pd(out_row + stock_i, x)` | write four cells |
| `++stock_i` | `stock_i = stock_i + 4` | move on four cells |

The scalar `east - 2.0 * centre + west` is written as `east + west` then
`fnmadd(two, centre, ...)`, which takes `2 * centre` off in one instruction.
The leftover cells at the end of the row go to `update_cells_scalar` from
§15.5. Because FMA rounds once where a multiply then an add rounds twice,
the SIMD rung differs from serial in the last bits (about 1e-15), well
inside the 1e-11 test tolerance.

### 15.7 The SIMD time loop (Phase 2)

**Before** (`solver_serial.cpp`)

```cpp
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_scalar(constants, current, next, var_j);
}
```

**After** (`solver_simd.cpp`)

```cpp
if (!avx_row_is_vectorised()) {
    std::fprintf(stderr, "[simd] built without AVX2 and FMA, using the scalar row\n");
}
...
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}
```

One function name changed in the loop. The warning at the top is new. On the
Mac, `update_row_avx` is the scalar fallback from the `#else` branch, and the
warning makes sure nobody mistakes a Mac run for an AVX2 result.

### 15.8 One pragma spreads the rows over the cores (Phase 3, Phase 5 Decision 5)

**Before** (`solver_simd.cpp`)

```cpp
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}

update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
grid.swap_buffers();
```

**After** (`solver_omp.cpp`)

```cpp
result.threads = omp_get_max_threads();
...
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}

double boundary_start = omp_get_wtime();
update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);
grid.swap_buffers();
double boundary_end = omp_get_wtime();
result.boundary_seconds += boundary_end - boundary_start;
```

The loop itself is untouched. The pragma above it tells OpenMP to split the
rows into one equal band per thread (`schedule(static)`) and to make every
variable we did not list a compile error (`default(none)`). Each thread
reads only `current` and writes only its own rows of `next`, so no locks are
needed. The boundaries stay serial after the loop. The two timer lines
around them were added in Phase 5 so the CSV shows how much of each step is
serial work. `result.threads` records the thread count the program really
got.

### 15.9 Each process takes a slab of rows (Phase 4)

**Before** (`solver_omp.cpp`)

```cpp
#pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
    update_row_avx(constants, current, next, var_j);
}

update_zero_row(constants, current, next);
apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
copy_top_row(constants, next);

grid.swap_buffers();
```

**After** (`solver_mpi.cpp`)

```cpp
RowRange owned_rows = rows_for_rank(num_variance_nodes, rank, num_ranks);
int first_interior = std::max(owned_rows.begin, 1);
int end_interior = std::min(owned_rows.end, num_variance_nodes - 1);
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
apply_stock_edges(constants, next, owned_rows.begin, owned_rows.end, step * dt);
if (owns_top_row) {
    copy_top_row(constants, next);
}

exchange_edge_rows(next, num_stock_nodes, owned_rows, rank, num_ranks);

grid.swap_buffers();
```

Four things changed, and the row function is still the same AVX2 one.

1. The loop runs over this process's rows only, `first_interior` to
   `end_interior`, instead of every row. Inside a process the rows are still
   split between threads with the same pragma.
2. Each boundary is done only by the process that owns it. The bottom
   process does the v = 0 row and the top process copies the top row. Every
   process writes the stock edges of its own rows, which is why
   `apply_stock_edges` was given a row range back in Phase 1.
3. `exchange_edge_rows` is new. The next step needs one row just outside
   each end of the process's rows, and those rows belong to the neighbours. So each
   process sends its first and last row and receives the row below and
   above it.
4. After the time loop, `MPI_Gatherv` collects every slab into rank 0 so it
   can read the price, and only rank 0 prints.

### 15.10 The experiment solvers (Phase 5, Decision 1)

Each experiment is a copy of the omp loop in §15.8 with one thing changed.
Two of them only change one word. `omp_scalar` calls `update_row_scalar`
instead of `update_row_avx`, and `omp_runtime` says `schedule(runtime)`
instead of `schedule(static)` so `OMP_SCHEDULE` can pick the schedule. The
other two change the shape of the loop.

**Before** (`solver_omp.cpp`, a new thread team every step)

```cpp
for (int step = 1; step <= num_timesteps; ++step) {
    #pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
    for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
        update_row_avx(constants, current, next, var_j);
    }
    update_zero_row(constants, current, next);
    apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
    copy_top_row(constants, next);
    grid.swap_buffers();
}
```

**After** (`solver_experiments.cpp`, `OmpPersistentSolver`, one team for the whole solve)

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

`parallel for` is really two pragmas in one, `parallel` (start the threads)
and `for` (share out the loop). Here they are split. `parallel` moves
outside the time loop, so the threads are started once, and every thread
runs every step. `omp for` shares out the rows as before, and `omp single`
lets one thread do the boundaries while the others wait at its end. This
tests whether starting the team every step costs anything.

**After** (`solver_experiments.cpp`, `OmpParallelEdgesSolver`, boundaries inside the loop)

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
grid.swap_buffers();
```

The loop now starts at row 0, and row 0 takes the v = 0 update. Each thread
also writes the stock edges of its own rows. Only `copy_top_row` stays after
the loop, because it reads the row below it, which some other thread might
still be writing. This tests how much the serial boundary work holds back
the omp rung.

### 15.11 Denormals flushed to zero (Phase 5, Decision 7)

There is no before and after inside `heston`, which did not change. The fix
is a separate five-line library, `slurm/flush_denormals.c`, that the
benchmark scripts load into every timed run.

```c
__attribute__((constructor)) static void flush_denormals(void) {
    _mm_setcsr(_mm_getcsr() | 0x8040);
}
```

`constructor` makes the function run before `main`. It sets two bits of the
CPU's floating point control register, so a number too small to be a normal
double becomes 0 instead of taking a slow path. Before this, short runs on
big grids left many cells below 1e-300, and the 2896 by 724 grid ran 2.4
times slower than its neighbours for that reason.

### 15.12 The triad arrays stop touching every page on one thread (Phase 5, Decision 8)

**Before** (rebuilt from the changelog, not committed on its own)

```cpp
std::vector<double> a(n), b(n), c(n);
#pragma omp parallel for schedule(static)
for (long i = 0; i < n; ++i) { a[i] = 0.0; b[i] = 1.0; c[i] = 2.0; }
```

**After** (`slurm/bench_peak.sh`, the C++ inside it)

```cpp
std::unique_ptr<double[]> a(new double[n]), b(new double[n]), c(new double[n]);
#pragma omp parallel for schedule(static)
for (long i = 0; i < n; ++i) { a[i] = 0.0; b[i] = 1.0; c[i] = 2.0; }
```

A `std::vector<double>(n)` fills itself with zeros the moment it is made, on
the main thread. Linux places a page of memory on the socket of the thread
that first writes it, so every page ended up on socket 0 and the parallel
loop never got the chance to place them. `new double[n]` leaves the memory
untouched, so the parallel loop is the first write and each page lands next
to the thread that uses it. `unique_ptr` still frees the memory on its own,
so there is no `delete`. Before this change 24 threads measured about the
same bandwidth as 12, which is what one socket's memory gives.

### 15.13 double becomes float (Phase 7)

**Before** (`stencil_avx.cpp`, `update_row_avx`)

```cpp
__m256d dt = _mm256_set1_pd(constants.dt);
while (stock_i + 3 <= last_interior) {
    __m256d centre = _mm256_loadu_pd(row_centre + stock_i);
    ...
    __m256d new_value = _mm256_fmadd_pd(dt, change, centre);
    _mm256_storeu_pd(out_row + stock_i, new_value);
    stock_i = stock_i + 4;
}
```

**After** (`stencil_avx.cpp`, `update_row_avx_float`)

```cpp
__m256 dt = _mm256_set1_ps(constants.dt);
while (stock_i + 7 <= last_interior) {
    __m256 centre = _mm256_loadu_ps(row_centre + stock_i);
    ...
    __m256 new_value = _mm256_fmadd_ps(dt, change, centre);
    _mm256_storeu_ps(out_row + stock_i, new_value);
    stock_i = stock_i + 8;
}
```

A float is half the size of a double, so the same 256-bit register holds
eight of them instead of four. That is the whole change. `__m256d` becomes
`__m256`, every `_pd` (packed double) becomes `_ps` (packed single), and the
loop moves eight cells at a time. In `solver_float.cpp` the sheets become
`std::vector<float>`, and the constants are worked out in double first and
then rounded, so both solvers start from the same numbers. Float was 1.75
times faster at 2048 by 512 on 24 threads, but it was rejected because the
price and gamma moved too far from the double answer.
