# SCRIPT — the word-for-word script for the video

One section per take, in the same order as TIMELINE.md and the deck. Read
it as written or in your own words, the sentences are short on purpose so
they survive being spoken. At a normal pace of about 150 words a minute
the whole script runs about nine and a half minutes, inside the
ten-minute stop. If a take runs long, the level 0, 1, 3 and 4 takes are
the place to tighten.
The time next to each heading is the start time from TIMELINE.md.

Numbers are written as digits so they are easy to read off. Say them
naturally, 52.6 is fifty-two point six.

## 0:00 Title

This is the price of a stock option, computed by solving a partial
differential equation on a grid. In the next nine and a half minutes I
will show you how that works, how I made the solver two and a half times
faster, and the two things I was confident about that measurement proved
wrong.

## 0:20 What an option is

A European call option is the right, but not the obligation, to buy a
share for a fixed strike price on a fixed day. At expiry its value is
known exactly, and the payoff has a kink at the strike. The problem is
the price today, ninety-one days before expiry.

## 0:40 Why not Black-Scholes

The standard model is Black-Scholes, and it assumes the share's
volatility is one constant number. If you take real market option prices
and compute the volatility they imply, it is not constant, it bends
across strikes. This curve comes from my own solver, and the dashed line
is the single Black-Scholes value. So I need a model where the
volatility itself can move.

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

The fair price is an average over every path the market could take,
discounted back to today. The Feynman-Kac theorem says that average,
which is the option price V, solves a deterministic equation. So my
program contains no random numbers.

## 1:55 The notation

Capital V is the option price, as a function of share price, variance
and time. Lowercase v is the variance, a different quantity. A subscript
on V means a partial derivative.

## 2:05 Why the noise becomes second derivatives

Over one time step the noise moves the price up or down by about root v
times S, with both directions equally likely. Averaging over the two
moves cancels the first derivatives, and the remaining term is half the
move size squared times the second derivative. The same argument gives
the variance term and the mixed term.

## 2:20 The equation, term by term

So every part of the model contributes exactly one term. The drifts give
first derivatives, the noises give second derivatives, the correlation
gives the mixed term, and discounting gives r times V.

## 2:30 The equation the program solves

This is the equation the program solves. Every quantity in it except V
is a known number before the solve starts. The unknown is the whole
option price function, and the number I report at the end is V at
today's share price and variance.

## 2:45 The solver

Share price runs across this grid, variance runs down it, and every cell
holds the option price for that pair. The spacing is the gap between
neighbouring values, about ten dollars between share prices here. Finer
spacing tracks the true curve better, and it multiplies the number of
cells. At expiry each cell equals the payoff, so the solver starts there
and steps backwards to today, rebuilding every cell from its nine
neighbours. A production solve is about a trillion cell updates, which
is why the speed of this one kernel matters.

## 3:15 Validation

Before benchmarking anything I test correctness. With xi at zero, Heston
reduces to Black-Scholes, which has an exact formula. The solver runs
unchanged and reproduces it with an error around ten to the minus three.
Put-call parity, a model-free identity, holds to five times ten to the
minus six.

## 3:35 Three independent methods

I also priced the same contract three independent ways. The grid gives
196.1684, a Fourier method that shares no code gives 196.1692, which
agrees to four parts in a million, and four hundred thousand Monte Carlo
paths bracket both.

## 3:50 Convergence

I also priced the option on a sequence of grids, each halving the
spacing. The error against the Fourier price falls about four times over
per halving, which is what a second-order scheme should do. It flattens
near one and a half cents, where the grid's cut-off and the reference's
own small error take over.

## 4:05 How the optimisation was structured

For optimisation I did not want a single before and after, because that
gives one number and no explanation. So there are seven complete
solvers, each a copy of the one below with exactly one technique added,
plus two negative controls that remove things the baseline already had.

## 4:20 Level 0

Level zero is the baseline arithmetic, line for line, behind the
ladder's kernel selection, a function pointer chosen once before the
time loop starts. It measured identical, so the dispatch itself costs
nothing.

## 4:35 Level 1

Level one hoists the invariant work. Spacing products like two ds and ds
squared were recomputed in every cell, and now happen once per step,
with per-row coefficients like half v and the mean reversion computed
once per row. No change, because GCC at O2 had already done all of it.

## 4:50 Level 2

Level two removes the five divisions per cell. V S was east minus west
divided by two ds. Each spacing constant is now inverted once per step,
and the loop multiplies by inv two ds instead. A divide costs up to
fourteen cycles and does not pipeline, and the compiler may not make
this change, because a rounded reciprocal changes the last bits.
Throughput went from 52.6 to 131.2 million cell updates per second.

## 5:15 Level 3

Level three walks memory in layout order. The g dot index calls become
three row bases, row, row above and row below, so the inner loop touches
consecutive doubles. No change, by construction, the baseline already
walked this way.

## 5:30 Level 4

Level four splits the loops. The variance-zero transport row gets its
own loop, the bounds become plain locals, and no cell ever tests which
equation it follows. Also no change, the baseline already peeled that
row out.

## 5:45 Level 5

Level five replaces the per-cell row plus stock i arithmetic, about ten
integer additions per cell, with raw row pointers, row mid, row above
and row below, so each access becomes plain pointer indexing. About two
percent, the only measurable gain after level two.

## 6:00 Level 6

Level six unrolls the inner loop by four. The body is level five's,
copied out with numbered locals V zero through V three, so four
independent dependency chains sit between one pair of loop-back
branches. Half a percent, the compiler already unrolls counted loops at
O2.

## 6:15 The ladder, measured

Here is the whole ladder measured on the cluster, and only level two
moved. Everything the compiler was allowed to do it had already done, so
the hand optimisation that paid is the one it is forbidden to make.
Overall the solver ended up 2.56 times faster.

## 6:35 The negative controls

The negative controls check that those null results were real. Swapping
the loop order makes identical arithmetic 2.1 times slower, so traversal
order genuinely matters, the baseline simply already had it right.
Fusing the boundary branch back in costs nothing, because the branch
goes the same way for an entire row and the predictor learns it. The
cache model on screen shows the swapped order moving five times the
memory per cell.

## 6:55 What the speedup cost

The speedup was not free. It moved the answer by fifteen units in the
last decimal place, a relative error of 3.3 times ten to the minus
fifteen, and an automated test holds every level and both controls to
that tolerance on every build.

## 7:10 The benchmark protocol

Every number so far comes from one protocol. Jobs run on a reserved
rangpur compute node through sbatch with exclusive access, and before
every sweep the binary rebuilds and the full test suite must pass. Each point is
the median of five repetitions, and every results file carries a header
naming the host, job id, compiler, flags, CPU and caches. The plot
scripts reject files without that header. One of my own controls still
misled me. The no-vectorise sweep at O2 proved nothing, because GCC 8.5
does not vectorise at O2 at all. Re-run at O3, the vectoriser does
nothing to my kernels but nearly doubles the baseline.

## 7:40 Scaling with problem size

This sweep grows the working set from one mebibyte to sixty-four, and
the speedup holds between two and a half and three and a half times. The
first run of this sweep was itself misleading, it confounded page size
with cache effects, so I re-ran it with the allocator pinned. The
baseline staying flat across a sixty-four-fold growth in working set
backs up the division story.

## 8:05 The roofline

I had assumed this solver was memory-bound. It is not. The level six
kernel does thirty-one floating point operations per twenty-four bytes,
and a separate job measured this node's actual scalar ceiling at 5.1
gigaflops rather than quoting a datasheet. Level six runs at 82 percent
of that ceiling while using about three gigabytes per second on a
four-channel machine. The baseline sits at 27 percent, limited by its
divisions. So optimisation moved the kernel from division-bound to
compute-bound.

## 8:35 The page-size effect

The largest single effect was not a code change at all. The same binary
on the same node runs at 134.8 million cell updates per second when its
buffers sit on two-mebibyte pages, and 81.1 on four-kibibyte pages. That
is 1.66 times from page size alone. Sixteen mebibytes of buffers needs
either eight huge pages or about four thousand small ones, against a
translation buffer with roughly a thousand entries. I caught it because
my median of five was averaging two different regimes, the allocator
gave the first repetition huge pages and the rest small ones.

## 9:10 Reflection

The main conclusion is that at O2 the compiler has already done the
textbook transformations, and the hand optimisations that pay are the
ones it is not allowed to make. The most
surprising result was page size, which was worth more than five of my
six techniques combined. What I would do differently is reserve the node
exclusively from the first job, because an early shared-node run told me
the loop-order control cost thirty times when the real answer was 2.1.
For Milestone 2, the kernel is now compute-bound, so SIMD comes first,
then OpenMP across rows, and the two-buffer design already makes every
cell independent within a step.

## 9:45 Close

I promised to show you how the pricer works, to make it faster, and to
name what I was wrong about. It is 2.56 times faster, it was never
memory-bound, and the compiler had already done five of my six
techniques. Thanks for watching.