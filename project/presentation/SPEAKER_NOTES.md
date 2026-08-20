# SPEAKER_NOTES.md — scripts for the harder sections

Spoken scripts for the slides that need them, written so that someone new to
the project could deliver them and answer the obvious follow-up question.
Each slide gets a paragraph of what to say in your own words, and where a
follow-up is likely, a second paragraph with the answer. Two sections are
covered, the model-to-equation run in the introduction and the whole
optimisation section. In the recording, levels 3 and 4 share one slide
and one take. Their per-level notes below back the individual slides in
the backup section, for the interview rather than the video.

If you remember nothing else, the section in one sentence. Seven versions of
the same solver, each adding exactly one textbook optimisation. Six did
nothing because the compiler had already done them. The only one that paid,
replacing division with multiplication, is the one the compiler is not
allowed to do, because it changes the answer in the last decimal places.

---

# The model-to-equation run

## From stochastic model to deterministic equation

The fair price of the option is an average, the expected payoff over every
path the market could take between now and expiry, discounted back to today.
The obvious way to compute an average over random paths is to simulate
millions of them and take the mean, and this program never does that. The
Feynman-Kac theorem says that this average, viewed as a function of where
you start, satisfies a deterministic equation. That function is the option
price, written capital V. So instead of simulating the randomness, the
program solves the equation the randomness leaves behind.

## The derivative notation

Capital V is the option price, the thing being solved for, as a function of
three inputs, the share price, the variance, and the time. Say once, out
loud, that capital V and lowercase v are different quantities. V is the
option price and v is the variance, and both appear in the same equation,
so this is the likeliest point of confusion on the slide. A subscript names the
direction of a partial derivative. V with subscript S is how fast the value changes as
the share price moves with everything else held fixed, and doubling the
subscript means the second derivative. On the grid every one
of these becomes a difference between neighbouring cells, which is why a
cell update reads nine neighbours.

## Why the noise terms become second derivatives

This is the one step in the construction that needs an argument, and it fits
in three sentences. Over one small time step, the model shocks the share
price up or down by an amount of size root v times S, with both directions
equally likely. Average the option value over the up case and the down case.
The first-derivative parts cancel, because whatever the up move gains the
down move loses. The remaining term is the curvature term, half the move
size squared times V sub SS. Squaring the move size root v times S produces
the half v S squared coefficient on the slide. The same argument applied to the variance noise, of
size xi root v, gives the V sub vv term, and the correlation between the two
noises produces the mixed V sub Sv term.

A likely question is why the drift terms get no such treatment. The drifts
are not random, so there is nothing to average over. They push value along
an axis and appear with a single first derivative, exactly as they are written in the
model.

## The equation, term by term

Each part of the model maps to one term. The two drifts become first-derivative terms.
The two noises become second-derivative terms through the averaging argument. The
correlation becomes the mixed term, and discounting becomes r times V.
Nothing else goes in, and nothing is left over.

## The equation the program solves

This is the end point of the derivation. Every symbol in the equation except V
is a known number before the solve starts. The interest rate and dividend
yield come from the market. Kappa, theta, xi, rho and the starting variance
are the model's five parameters, chosen up front. The strike is not in the
equation at all, it sits in the expiry payoff. The one unknown is the whole
function V, one value for every combination of share price, variance and
time, and the grid is that function written out cell by cell. The program
fills the grid in backwards from expiry, and the single number reported at
the end is V read out at today's share price and today's variance.

---

## How the optimisation was structured

A single before and after gives one speedup number and no explanation of
where it came from. So there are seven complete solvers. Each is a copy of
the one below it with exactly one technique added and nothing else changed,
which means the difference between neighbouring rungs measures one
technique. No version is benchmarked until it has proved it still produces
the reference answer.

A likely question is why not just profile. A profiler says where time is
spent, not which change moved it. The ladder attributes every improvement to
a named technique, including the ones worth zero, and the zeros turn out to
be the finding.

## Level 0, harness control

Level 0 adds no optimisation at all. It is the baseline arithmetic, line for
line, running through the ladder's kernel selection. It exists to prove that
the selection machinery is free. It measured the same as the baseline, so
every difference higher up the ladder is real.

The machinery in question is a function pointer. The solver picks its
inner-loop function once, before the time loop starts, so the choice of
level never puts a branch inside the measured loop.

## Level 1, hoisting

Level 1 moves each computation to the level where it actually varies. Grid
spacings, their reciprocals, and coefficients that only change per row are
computed once, outside the loops, instead of once per cell. This is the
first optimisation in every textbook. It measured no change. GCC at -O2 had
already hoisted all of it, and the level exists to demonstrate that.

The compiler passes involved are loop-invariant code motion and common
subexpression elimination. Both preserve the answer exactly, so the compiler
is free to do them, and does.

## Level 2, strength reduction

Level 2 replaces each division by a grid constant with a multiplication by
its reciprocal, computed once. The inner loop had five divisions per cell. A
double-precision divide costs five to fourteen cycles and blocks the divide
unit, while a multiply costs about one cycle. This took the solver from 52.6
to 131.2 million cell updates per second, 2.5 times faster, and it is the
only rung that moved.

Here is why it was the only one. Dividing by a number and multiplying by its
reciprocal round differently in the last bits of the answer. The compiler at
-O2 must preserve the answer bit for bit, so it has to leave the divisions
alone. Everything the compiler is allowed to do, it had already done. The
hand optimisations worth making are the ones it is not allowed to make.

On correctness, the change moved the price by 15 units in the last decimal
place of a $196 answer, which is 3.3 parts in 10^15, and an automated test
holds every level to that tolerance.

## Level 3, traversal order

Level 3 makes the loops walk the grid in the same order the cells sit in
memory, so the cache streams data instead of jumping around. It measured
zero, and by construction rather than by accident. The baseline was already
written to walk memory in the right order, so there was nothing left to win.
What the right order is worth is answered by the negative control in a
moment.

The grid is row-major. Each row of stock prices is contiguous and the inner
loop walks along a row, so consecutive iterations touch consecutive doubles
and every 64-byte cache line is fully used.

## Level 4, loop splitting

The bottom row of the grid, where variance is zero, obeys a different and
simpler equation. Level 4 gives that row its own loop rather than testing
for it in every cell. Also zero, also by construction, because the baseline
already had the row peeled off. The control for this one puts the branch
back in on purpose, and that turns out to cost nothing, because the branch
goes the same way for an entire row and the predictor gets it right every
time.

The row is special because every diffusion term in the equation carries a
factor of the variance v, so at v equals zero they all vanish and only a
transport equation is left.

## Level 5, induction variables

Level 5 replaces the index arithmetic, row times width plus column
recomputed for all nine neighbours of every cell, with pointers that step
forward one place per cell. That folds about ten integer additions per cell
into the address arithmetic the hardware does anyway. It was worth about two
percent, which is real, and the only measurable effect after level 2.

It is small because the core does integer address arithmetic on separate
execution ports from floating point, so most of that work was already
happening alongside the floating-point math. The neighbour offsets that
remain are fixed displacements inside the load instructions, so the
address unit computes them as part of the load and no separate add
instruction is issued.

## Level 6, unrolling

Level 6 processes four cells per trip around the loop instead of one, giving
a quarter of the loop tests and more independent arithmetic in flight. Each
cell's update is a chain of dependent operations, so while one chain waits
on a result the core can issue from another, which hides instruction
latency. It was worth half a percent. The compiler already unrolls counted loops at -O2,
so doing it by hand mostly duplicated its work. That completes the ladder.
One rung out of six moved.

On the laptop the same change was a 23 percent regression. Same source,
opposite sign, which is why every number in this talk comes from the target
machine.

## The negative controls

Two rungs came back null because the baseline already did the right thing.
That claim can be tested by writing the wrong version deliberately and
timing it. Swapping the loop order, so the inner loop strides through memory
instead of streaming, makes the same arithmetic 2.1 times slower. That is
what traversal order was worth all along. Fusing the boundary row back in
costs nothing, because the branch predictor absorbs it. So one of the two
suspected optimisations was genuinely load-bearing and the other was not
needed at all.

The 2.1 comes from memory traffic. A modelled cache run shows the swapped
order moves about five times the memory per cell update, 40.8 bytes against
8.6, because each fetched cache line is used once and evicted before the
loop comes back for its neighbours.

## What the speedup cost

Every optimisation is a trade. This one cost fifteen units in the last
decimal place of the answer, a relative error of 3.3e-15, and the test suite
pins every level and both controls to that tolerance, on every build, before
any benchmark is allowed to run.
