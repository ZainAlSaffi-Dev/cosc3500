# Walkthrough. What this program does, in plain English

This file explains the whole program in order, from the problem it solves down
to the individual loops. It is written to be read out loud and understood
without any finance or C++ background. Every claim here can be traced to a
file and a line number in this repository.

---

## Part 1. The problem being solved

An option is a contract that gives someone the right, but not the obligation,
to buy a share at a fixed price on a fixed future date. The fixed price is
called the strike and the fixed date is called expiry. If the share ends up
worth more than the strike, the holder buys at the strike, sells at the market
price and keeps the difference. If the share ends up worth less, the holder
does not exercise it and the contract is worth nothing.

The question this program answers is what that contract is worth today, before
anyone knows where the share price will end up.

The classic answer to that question is the Black-Scholes formula, and it makes
one assumption that is known to be false. It assumes that the size of the
share price's random daily moves, called the volatility, stays the same
forever. Real markets do not behave that way. Calm periods and volatile
periods alternate, and the size of the daily moves is itself random.

The Heston model replaces exactly that assumption. The share price moves
randomly, and the size of those moves also changes randomly over time. There
are now two random quantities instead of one. The share price is the first,
and the variance, which is the square of the volatility, is the second. The
two are allowed to be correlated, and in practice that correlation is
negative, meaning that when prices fall the volatility tends to rise. In this
project that correlation is set to negative zero point seven, which is a
realistic figure for an equity index.

The cost of that extra realism is that the Heston model has no simple formula
you can type into a calculator. To get an answer you have to solve a partial
differential equation numerically, and this program does that.

---

## Part 2. Why the answer becomes a grid

There is a result in mathematics called the Feynman-Kac theorem, and it is the
reason this program contains no random numbers at all. The theorem says that
the average outcome of a random process is the solution of a certain
differential equation. In other words, you can either simulate a
million random share price paths and average the results, or you can solve one
equation and get the same answer exactly. Solving the equation is far cheaper
and far more accurate, so that is the route taken here.

The equation has three inputs, which are the share price, the variance and the
amount of time left until expiry. The program stores the solution as a large
table of numbers called a grid. Every column of the grid is a possible share
price, every row is a possible variance level, and the whole grid holds the
option's value at one instant in time.

The program fills the grid in backwards in time. At expiry the answer is
completely known, because at that moment the option is worth exactly its
payoff. If the share is above the strike the option is worth the difference
and otherwise it is worth zero. There is no uncertainty left, so no modelling
is required. That known grid is where the program starts, and it is created
in `Grid::init_payoff` in `src/grid.cpp:37`.

From there the program steps backwards in time, one small step at a time. At
each step it works out what every cell of the grid must have been worth a
moment earlier, given what the cells around it are worth now. Repeat that
hundreds of thousands of times (980000 steps on the reference grid) and you
arrive back at today, and the cell that matches today's actual share price and
today's actual variance is the answer.

---

## Part 3. How the grid is stored in memory

A grid of numbers is the obvious thing to store as a list of lists, but that
is a poor choice for performance. The program instead stores the entire grid
as one contiguous array of numbers, which `std::vector<double>` provides. This
is the same idea as a NumPy array in Python. The class that owns this is
`Grid` in `include/grid.h`.

To find a particular cell you use the formula in `Grid::index` at
`include/grid.h:38`, which is the variance row number multiplied by the number
of share price columns, plus the share price column number. The important
consequence is that moving one step along the share price direction moves you
one step along in memory, whereas moving one step along the variance direction
is a stride of one whole row. Modern processors read memory far faster in
address order than at scattered addresses, so every loop in this program is
written to run along the share price direction on the inside. Part 8 comes
back to this, because it turned out to be worth a factor of two.

There are two of these arrays, called buffers, and the reason is subtle but
important. When the program computes the new value of a cell it needs the
current values of the cells around it. If it wrote the new value straight back
over the old one, then the next cell along would read a value that had already
been updated, and it would be mixing two different instants in time. So one
buffer is read and the other is written, and at the end of the step the two
buffers swap pointers. That swap happens in `Grid::swap_buffers` at
`src/grid.cpp:31`, and it takes the same tiny time however big the grid is,
because `std::swap` on two vectors exchanges their internal pointers and copies
none of the numbers.

---

## Part 4. What happens in one timestep

The update for one timestep is written out in full in `BaselineSolver::solve`
in `src/solver_baseline.cpp:16`.

For every cell in the middle of the grid the program reads that cell and its
eight neighbours, which is the cell to the left and right, the cell above and
below, and the four diagonals. From those nine numbers it estimates five
slopes and curvatures, which are how fast the value changes as the share price
moves, how that rate of change is itself changing, the same two quantities in
the variance direction, and finally a mixed term that captures how the two
directions interact. Those five estimates are built at
`src/solver_baseline.cpp:68` to `72`, and each one is just a difference between
neighbouring cells divided by the spacing between them, which is the standard
way of approximating a derivative on a grid.

Those five quantities are then used in the Heston equation itself at
`src/solver_baseline.cpp:75`. Each term of that line has a plain meaning. The
first term is the diffusion of the share price, which comes from the random
movement of the share price. The second is the cross term, which is the effect
of the share price and the variance being correlated. The third is the
diffusion of the variance, which comes from the random movement of the
variance itself. The fourth is the drift of the share price under the risk
free rate less the dividend yield. The fifth is the pull of the variance back
towards its long run average, which is the mean reversion. The last term
discounts the whole value by one step of interest.

The reason this is called an explicit scheme is that the new value of each cell
is computed directly from values that are already known. There is no system of
equations to solve, which makes the code simple and makes every cell
independent of every other cell. That independence makes the problem easy to
parallelise later, and it is why this shape of solver was chosen for the
project.

The price of that simplicity is that the timestep cannot be made too large. If
it is, the scheme becomes unstable and the numbers grow without bound within
a few dozen steps. The program computes the largest safe timestep at
`src/solver_baseline.cpp:34` and prints it, and the animation in
`results/blowup.mp4` shows what happens when that bound is deliberately
ignored.

---

## Part 5. The four edges of the grid

The loop above only covers cells that have all nine neighbours available. The
cells on the outside edge of the grid do not, so each edge has to be handled
separately. There are four of them and each has its own reason.

The left edge is where the share price is zero. A share worth nothing stays
worth nothing forever, so at that edge the outcome is certain. A call option is
worthless there and a put option is certain to pay the strike, which just needs
discounting back to today. This is at `src/solver_baseline.cpp:109`.

The right edge is where the share price is four times the strike. The option
is so deeply in the money there that exercising it is effectively guaranteed,
so the value is the share price, discounted by the dividend yield, less the
discounted strike. This is at `src/solver_baseline.cpp:106`.

The top edge is where the variance is at its maximum. By that point the option
value has stopped responding to further increases in volatility, so the program
imposes a flat slope by copying the row underneath. This is at
`src/solver_baseline.cpp:121`, and it deliberately runs last, because the row
it copies has to be finished first.

The bottom edge, where the variance is zero, is the interesting one and it is
the trap that catches most implementations. You might expect to impose a value
by hand there, but that would be wrong. Every diffusion term in the equation
carries a factor of the variance, so when the variance is zero all of them
vanish and what is left behind is a simpler transport equation. The program
solves that reduced equation on the bottom row, at
`src/solver_baseline.cpp:93`.

There is one further subtlety on that bottom row. The slope in the variance
direction has to be measured using the row above only, because there is no row
below to average with. Fortunately that is also the correct choice for a
deeper reason. The mean reversion term pushes the variance upwards away from
zero, so information flows into the grid from that edge, and the one sided
difference points in the direction the information is actually coming from.

---

## Part 6. Reading the answer out

Once the backward time stepping reaches today, the finished grid holds the
option's value for every share price and every variance level. The program
only needs one of those, which is the cell matching today's actual market.
That readout is `Solver::extract_result` in `src/solver.cpp:34`.

One detail of the readout cost real accuracy before it was found. Today's
share price almost never sits exactly on a grid column. Reading the nearest
column instead sounds harmless, but on the small test grid it produced a price
of 205.32 when the true answer was 195.25, which is an error of nearly five
percent. The fix is to blend the four cells surrounding the true market point
in proportion to how close it is to each, which is done by
`bilinear_interpolate` at `src/solver.cpp:18`.

The same function is then reused to compute the sensitivities, which traders
call the Greeks. Delta is how much the option price moves when the share price
moves, gamma is how fast delta itself changes, and vega is how much the price
moves when volatility moves. Because all three are measured off the same
blended surface as the price, the reported delta really is the derivative of
the reported price.

One more detail is worth knowing for questions. When the grid happens to be
aligned so that the market point sits exactly on a node, the blending weight
comes out as exactly zero and the blend returns the original cell value bit
for bit. That is why adding interpolation changed nothing about any previously
recorded number, and that was verified.

---

## Part 7. How we know the answer is right

There are four independent checks, and each one covers something the others
cannot.

The first check is a collapse test, in `tests/test_bs_collapse.cpp`. If you
set the volatility of the variance to zero and start the variance exactly at
its long run average, the variance can never move, and the Heston model
becomes the ordinary Black-Scholes model. Black-Scholes has an exact formula,
so the answer is known in advance. The solver has no idea the formula exists,
so agreement between the two tests the payoff, the stencil, the boundaries and
the readout together.

The second check is put-call parity, in `tests/test_parity.cpp`. There is a
relationship between the price of a call and the price of a put on the same
share that has to hold in absolutely any model, because violating it would let
someone make money with no risk. Checking it exercises the full parameter set,
including the cross term and both variance boundaries that the collapse test
cannot reach because it has to switch the volatility of the variance off.

The third check is convergence, in `results/convergence.csv` and
`results/convergence.png`. If the code is a correct approximation of the
equation, then making the grid finer must make the answer approach a fixed
limit. The solver passes this check. Between the two coarser grids the error
against the Fourier price falls from 0.534 dollars to 0.026, an observed order
of about two point two in the timestep along this refinement path, and on the
finest grid it flattens at about 0.014 dollars, a floor that this refinement
does not remove (`docs/RESULTS.md` section 7).

The fourth check is a completely independent price, from
`scripts/monte_carlo_check.py`. This one throws away the grid entirely and
prices the option two other ways, once by simulating four hundred thousand
random share price paths and averaging the outcomes, and once by Fourier
inversion of the model's characteristic function. The grid says 196.1684, the
Fourier method says 196.1692, which agree to four parts in a million, and the
simulation says 196.1345 give or take 0.5524, which brackets both. Three
methods with almost nothing in common agreeing on the same number is strong
evidence.

There is a fifth check that only matters once optimisation starts, which is
`tests/test_ladder_matches.cpp`. It runs the reference solver and every faster
rung and compares the price, delta, gamma and vega. Every rung must agree with
the reference to one part in ten to the eleventh, and the OpenMP and MPI rungs
must agree with the SIMD rung exactly, bit for bit. No rung is ever benchmarked
until it has passed this, and the benchmark scripts enforce that by refusing
to run if the tests fail.

---

## Part 8. The Milestone 1 optimisation ladder

The obvious way to present an optimisation is to show one slow version and one
fast version and quote the ratio. Milestone 1 deliberately avoided that,
because a single ratio cannot show which change made the difference.

Instead Milestone 1 had seven complete versions of the solver, numbered zero
to six, in `src/solver_opt_kernels.cpp`. That file was removed at the start of
Milestone 2 and is kept at git tag `milestone1`. The techniques that helped are
now in `StencilConstants` in `src/stencil.cpp` and in the serial rung
(STUDY_GUIDE.md section 13). Each version was a full copy of the one below it
with exactly one technique added and nothing else changed. That means the
difference in timing between any two neighbours could only have come from that
one technique. This is why the file was long and repetitive, and the
repetition was deliberate. Sharing setup code between the rungs would have let
one rung change what its neighbour was measuring, which is precisely what this
design ruled out.

Level zero applied no technique and repeated the reference solver line for
line. Its job was to prove that the machinery around the kernels cost nothing
and changed no answers.

Level one named repeated subexpressions and moved calculations out to the
level where they actually vary. The reference version recomputes the share
price for a column inside every one of a million cells, when it really only
changes once per column. Level one was required to stay bit identical, so it
could give a subexpression a name but could never change the order things were
multiplied in.

Level two replaced divisions with multiplications by a precomputed
reciprocal. A division of two doubles costs somewhere between five and
fourteen processor cycles and cannot be pipelined, while a multiplication
costs about one. This is the rung where the answers were allowed to change in
the last few bits, because dividing by two times the spacing and multiplying
by the reciprocal of two times the spacing are genuinely different numbers
once that reciprocal has been rounded.

Level three made the memory friendly traversal order explicit. Level four
split the loops so the special bottom row was handled separately with no
conditional test in the main loop. Both of those were predicted to measure as
nothing, and they did, because the reference solver was already written
correctly in both respects. A confirmed prediction counts as a result.

Level five held a direct pointer to the start of each row instead of
recomputing the address for each of the ten memory accesses a cell makes.
Level six unrolled the innermost loop four times, so the loop bookkeeping was
spread over four cells and the processor had four independent chains of
arithmetic to work on in parallel.

Because levels three and four came out as nothing, there was a separate pair
of deliberately bad versions to measure what they would have been worth. The
first swapped the two loops so that the inner loop moved through memory with a
stride of one whole row every iteration. The second fused the bottom row back
into the main loop so that every single cell had to test whether it was on
that row. Both did arithmetic identical to the good versions, so any
difference in timing was purely structural, and both were held to the same
correctness tolerance, so a wrong answer would have meant a bug in the
control.

---

## Part 9. What the Milestone 1 numbers said

Every benchmark ran on a rangpur compute node through the job scheduler, on a
node reserved exclusively so nothing else could interfere, and pinned to the
same physical machine every time so the sweeps could be compared with each
other. That machine is an Intel Xeon E5-2670 v3.

The Milestone 1 headline was that the optimised solver ran 2.56 times faster
than the reference, and that essentially all of that came from a single rung.

Every technique on the ladder except level two came back inside the run to run
noise, because the compiler had already done all of them at the ordinary
optimisation setting. Level two is the only transformation on the whole list
that the compiler is forbidden to perform, because turning a division into a
multiplication by a reciprocal changes the answer in the last bits and the
code is not compiled with fast maths enabled.

So the conclusion is that at this optimisation level, the hand optimisations
worth doing are the ones the compiler is not allowed to do for you. That
finding is only visible because each rung was measured separately, and reported
null results are explicitly credited by the course.

Two other things came out of the measurements that are worth saying out loud.

The first is that an assumption turned out to be wrong. The solver was assumed
to be limited by memory bandwidth, and the measurements showed that arithmetic
limits it. Counting the arithmetic and the memory traffic straight off the
level six kernel gives thirty one arithmetic operations for every twenty four
bytes moved, and a separate job measured what the node can actually sustain,
so no datasheet figure was needed. Level six runs at eighty two percent of
that measured ceiling while asking for only about three gigabytes per second
of memory bandwidth on a four channel machine. The reference solver sits at
twenty seven percent, and it is limited by the five divisions in every cell.
So the optimisation moved the program from being division limited to being
genuinely compute limited, and from there it can only go faster by doing more
arithmetic per instruction, which vector instructions provide. That led
directly into Milestone 2.

The second is a measurement artefact that was caught and explained. The same
binary on the same node ran at either 134.8 or 81.1 million cell updates per
second depending on nothing but whether the two eight megabyte grid buffers
were backed by large memory pages or small ones. Sixteen megabytes needs eight
address translation entries one way and four thousand and ninety six the other.
The memory allocator switches between the two regimes on its own, so the first
repetition of a benchmark run was getting one regime and the rest were getting
the other, which meant a naive median was silently averaging two different
memory regimes. That was established by elimination, by checking that the
processor clock stayed pinned throughout and then by forcing each regime by
hand and watching the effect appear and disappear on demand.

Milestone 2 built the parallel ladder on top of that serial work. On the same
kind of node the baseline runs at 57 million cell updates per second, the
serial rung at 159, the SIMD rung at 456 and OpenMP on 24 threads at 6367, and
two MPI processes of 12 threads each on one node reach 6984, which is 123
times the baseline (`docs/RESULTS.md` section M2.1).

---

## Part 10. How the code is organised

The program is deliberately small and the pieces have clean responsibilities.

`include/params.h` holds plain records describing the contract, the market, the
model parameters and the grid resolution. These are simple data holders with no
behaviour, which is the same role a dataclass plays in Python.

`src/params.cpp` reads a configuration file of key and value lines and throws
an error if a key is unrecognised or a value does not parse.

`include/grid.h` and `src/grid.cpp` own the two buffers of numbers, the spacing
between nodes, the index arithmetic, the initial payoff and the buffer swap.

`include/solver.h` declares an abstract base class with one method that must be
implemented, which is the same idea as an abstract base class in Python. Each
rung of the ladder is one implementation, `BaselineSolver`, `SerialSolver`,
`SimdSolver`, `OmpSolver` and `MpiSolver`, and `include/solver_experiments.h`
declares the experiment solvers. A factory function, `make_solver`, hands back
whichever one was asked for by name. The pointer that holds it is a unique
pointer, which means there is exactly one named owner and
the object is destroyed automatically when that owner goes out of scope. There
is no manual memory management anywhere in the project.

`src/solver.cpp` holds the parts every solver shares, which are the readout,
the stability report and the factory.

`src/stencil.cpp` holds the constants that are worked out once before the time
loop and the scalar row update, and `src/stencil_avx.cpp` holds the AVX2
version of the same row update, the only file built with the AVX2 flags. Each
rung's time loop has its own file, `src/solver_baseline.cpp`,
`src/solver_serial.cpp`, `src/solver_simd.cpp`, `src/solver_omp.cpp` and
`src/solver_mpi.cpp`, and the experiments are in `src/solver_experiments.cpp`,
`src/solver_align.cpp`, `src/solver_float.cpp` and, for `mpi_overlap`,
`src/solver_mpi.cpp`.

`src/io.cpp` writes grid snapshots for the animations and prints the one line
of comma separated results that every script reads. It is careful to write
warnings to the error stream so that the results stream stays clean.

`src/main.cpp` starts MPI, parses the command line, chooses a solver and
prints the result, from the first process only when several run. It requires
`--config <file>` first, reads that file, and then applies the remaining flags
in pairs on top of it, so a flag always overrides the file.

`src/black_scholes.cpp` is the closed form formula used only by the tests. The
solvers never call it, so it is a genuinely independent check.
