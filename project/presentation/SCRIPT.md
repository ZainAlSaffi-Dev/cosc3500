# SCRIPT — the word-for-word script for the video

One section per take, in the same order as TIMELINE.md and the deck. Read
it as written or in your own words, the sentences are short on purpose so
they survive being spoken. At a normal pace the whole script runs about
nine and a half minutes, inside the ten-minute stop. If a take runs long,
the level 0, 1, 3 and 4 takes are the place to tighten.
SPEAKER_NOTES.md carries the longer explanations and the likely follow-up
answers behind these takes.

Numbers are written as digits so they are easy to read off. Say them
naturally, 52.6 is fifty-two point six.

## 0:00 Title

This is the price of a stock option, computed by solving a partial
differential equation on a grid. In the next ten minutes I will show you
how that works, how I made the solver two and a half times faster, and
the two things I was confident about that measurement proved wrong.

## 0:20 What an option is

A European call option is the right, but not the obligation, to buy a
share for a fixed strike price on a fixed day. At expiry its value is
known exactly, and the payoff has a kink at the strike. The problem is
the price today, ninety-one days before expiry.

## 0:40 Why not Black-Scholes

The standard model is Black-Scholes, and it assumes the share's
volatility is one constant number. Back the volatility out of real
market option prices and it is not constant, it bends across strikes.
The dashed line is the single Black-Scholes value, and the curve is what
my solver prices. So the volatility itself has to be allowed to move.

## 1:05 The Heston model

The Heston model does that with two equations. The share price moves
with volatility square root of v, and the variance v follows its own
random process. The two noises are negatively correlated, because
volatility tends to rise when prices fall.

## 1:20 The symbols

Every symbol in those equations is defined here. The two worth
remembering are kappa, the rate at which the variance returns to its
long-run level theta, and xi, the volatility of volatility, which sets
how large the variance fluctuations are.

## 1:35 Where the strike enters

One thing is missing. There is no strike price in the model, because
these equations describe the share, not the contract. The strike enters
later, through the payoff at expiry.

## 1:45 From model to equation

The fair price is an average, the expected payoff over every path the
market could take, discounted back to today. The Feynman-Kac theorem
says that average, the option price V, solves a deterministic equation.
So instead of simulating the randomness, the program solves the equation
the randomness leaves behind.

## 2:00 The notation

Capital V is the option price, as a function of share price, variance
and time. Lowercase v is the variance, a different quantity. A subscript
on V means a partial derivative.

## 2:10 Why the noise becomes second derivatives

Over one small step the noise moves the price up or down by about root v
times S, with both directions equally likely. Average the two and the
first derivatives cancel, because whatever the up move gains the down
move loses. The remaining term is half the move size squared times the
second derivative.

## 2:25 The equation, term by term

So every part of the model contributes exactly one term. The drifts give
first derivatives, the noises give second derivatives, the correlation
gives the mixed term, and discounting gives r times V.

## 2:35 The equation the program solves

This is the end point of the derivation. Everything in it except V is a
known number, the rates come from the market and the five Heston
parameters are chosen up front. The strike is not in the equation at
all, it sits in the expiry payoff. The one unknown is the whole function
V, and the answer is V read at today's share price and variance.

## 2:55 The solver

Share price runs across this grid, variance runs down it, and every cell
holds the option price for that pair. The spacing is the gap between
neighbouring values, about ten dollars between share prices here. Finer
spacing tracks the true curve better, and it multiplies the number of
cells. At expiry each cell equals the payoff, so the solver starts there
and steps backwards to today, rebuilding every cell from its nine
neighbours. A production solve is about a trillion cell updates, which
is why the speed of this one kernel matters.

## 3:25 Validation

Before benchmarking anything I test correctness. With xi at zero, Heston
reduces to Black-Scholes, which has an exact formula. The solver runs
unchanged and reproduces it with an error around ten to the minus three.
Put-call parity, a model-free identity, holds to five times ten to the
minus six.

## 3:45 Three independent methods

I also priced the same contract three independent ways. The grid gives
196.1684, a Fourier method that shares no code gives 196.1692, which
agrees to four parts in a million, and four hundred thousand Monte Carlo
paths bracket both.

## 4:00 Convergence

I also priced the option on a sequence of grids, each halving the
spacing. The error against the Fourier price falls about four times over
per halving, which is what a second-order scheme should do. It flattens
near one and a half cents, where the grid's cut-off and the reference's
own small error take over.

## 4:15 How the optimisation was structured

A single before and after gives one speedup number and no explanation of
where it came from. So there are seven complete solvers, each a copy of
the one below with exactly one technique added, and the difference
between neighbours measures one technique. Two negative controls remove
things the baseline already had.

## 4:30 Level 0

Level zero is the baseline arithmetic, line for line, behind the
ladder's kernel selection, a function pointer chosen once before the
time loop starts. It measured identical, so the dispatch itself costs
nothing.

## 4:45 Level 1

Level one hoists the invariant work, the first optimisation in every
textbook. Spacing products like two ds were recomputed in every cell and
now happen once per step, per-row coefficients like half v once per row.
No change, because GCC at O2 had already done all of it.

## 5:00 Level 2

Level two removes the five divisions per cell. V S was east minus west
divided by two ds. Each spacing constant is now inverted once per step,
and the loop multiplies by inv two ds instead. A divide costs up to
fourteen cycles and does not pipeline, and the compiler may not make
this change, because a rounded reciprocal changes the last bits.
Throughput went from 52.6 to 131.2 million cell updates per second.

## 5:25 Level 3

Level three walks memory in layout order, three row bases per row
instead of index calls, so the inner loop touches consecutive doubles.
No change, by construction, the baseline already walked this way. What
that order is really worth, the negative control will show.

## 5:40 Level 4

Level four splits the loops. The variance-zero row follows a simpler
equation, so it gets its own loop, and no cell ever tests which equation
it follows. Also no change, the baseline already peeled that row out.

## 5:55 Level 5

Level five replaces the row plus stock i arithmetic, about ten integer
additions per cell, with raw row pointers, row mid, row above and row
below. That folds the indexing into address arithmetic the hardware does
anyway. About two percent, the only measurable gain after level two.

## 6:10 Level 6

Level six unrolls the inner loop by four, level five's body copied with
numbered locals, four independent chains between one pair of branches.
Half a percent here. On my laptop the same change was a 23 percent
regression, same source, opposite sign, which is why every number in
this talk comes from the cluster.

## 6:30 The ladder, measured

Here is the whole ladder measured on the cluster, and only level two
moved. Everything the compiler was allowed to do it had already done, so
the hand optimisation that paid is the one it is forbidden to make.
Overall the solver ended up 2.56 times faster.

## 6:50 The negative controls

Two rungs came back null because the baseline already did the right
thing, and that claim can be tested by writing the wrong version on
purpose and timing it. Swapping the loop order makes identical
arithmetic 2.1 times slower, so traversal order was genuinely
load-bearing. Fusing the boundary branch back in costs nothing, the
predictor learns a branch that goes one way for a whole row, so that one
was never needed.

## 7:15 What the speedup cost

The speedup was not free. It moved the answer by fifteen units in the
last decimal place, a relative error of 3.3 times ten to the minus
fifteen, and an automated test holds every level and both controls to
that tolerance on every build.

## 7:30 The benchmark protocol

Every number comes from one protocol, a reserved compute node through
sbatch with exclusive access, a rebuild and full test pass before every
sweep, and the median of five repetitions, with every results file
stamped with its host, compiler and flags. One control still misled me.
The no-vectorise sweep at O2 proved nothing, because GCC 8.5 does not
vectorise at O2 at all. Re-run at O3, the vectoriser does nothing to my
kernels but nearly doubles the baseline.

## 7:55 Scaling with problem size

This sweep grows the working set from one mebibyte to sixty-four, and
the speedup holds between two and a half and three and a half times. The
first run of this sweep was itself misleading, it confounded page size
with cache effects, so I re-ran it with the allocator pinned. The
baseline staying flat across a sixty-four-fold growth in working set
backs up the division story.

## 8:20 The roofline

I had assumed this solver was memory-bound. It is not. The level six
kernel does thirty-one floating point operations per twenty-four bytes,
and a separate job measured the node's real scalar ceiling at 5.1
gigaflops rather than trusting a datasheet. Level six runs at 82 percent
of that ceiling, the baseline at 27, limited by its divisions.
Optimisation moved the kernel from division-bound to compute-bound.

## 8:45 The page-size effect

The largest single effect was not a code change at all. The same binary
on the same node runs at 134.8 million cell updates per second on
two-mebibyte pages and 81.1 on four-kibibyte pages, 1.66 times from page
size alone. Sixteen mebibytes of buffers is eight huge pages or four
thousand small ones, against a translation buffer with about a thousand
entries. It surfaced because my median of five was averaging two
regimes, huge pages on the first repetition, small on the rest.

## 9:15 Reflection

The main conclusion is that at O2 the compiler has already done the
textbook transformations, so the hand optimisations that pay are the
ones it is not allowed to make. The most surprising result was page
size, worth more than five of my six techniques combined. What I would
do differently is reserve the node exclusively from the first job. A
shared-node run said the loop-order control cost thirty times, and the
real answer was 2.1. For Milestone 2 the kernel is now compute-bound, so
SIMD first, then OpenMP across rows.

## 9:45 Close

I promised to show you how the pricer works, to make it faster, and to
name what I was wrong about. It is 2.56 times faster, it was never
memory-bound, and the compiler had already done five of my six
techniques. Thanks for watching.
