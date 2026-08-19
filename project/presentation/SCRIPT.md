# SCRIPT — the word-for-word script for the video

One section per take, in the same order as TIMELINE.md and the deck. The
script is budgeted at 150 words per minute, which is a normal presentation
pace, and each take's word count fits its slot in TIMELINE.md. The whole
script is about 1440 words, which is 9 minutes 36 seconds of speech, so
there is no room for improvised elaboration. If a take does not fit its
slot when spoken, the fix is to re-record that take, not to talk faster.

Each take is written to stand on its own, so nothing needs to be added on
the day. SPEAKER_NOTES.md carries the longer versions and the likely
follow-up answers, for the interview rather than the recording.

Numbers are written as digits so they are easy to read off. Say them
naturally, 52.6 is fifty-two point six.

## 0:00 Title

This is the price of a stock option, computed by solving a partial
differential equation on a grid. In the next ten minutes I will show you
how that works, how I made the solver two and a half times faster, and
the two things measurement proved me wrong about.

## 0:20 What an option is

A European call option is the right, but not the obligation, to buy a
share for a fixed strike price on a fixed day. At expiry its value is
known exactly, and the payoff has a kink at the strike. The problem is
the price today, ninety-one days before expiry.

## 0:40 Why not Black-Scholes

The standard model is Black-Scholes, and it assumes the share's
volatility is one constant number. Back the volatility out of real
market option prices and it is not constant, it bends across strikes.
The dashed line is Black-Scholes, the curve is what my solver prices.
So the volatility itself has to move.

## 1:00 The Heston model

The Heston model does that with two equations. The share price moves
with volatility square root of v, and the variance v follows its own
random process. The two noises are negatively correlated, because
volatility tends to rise when prices fall.

## 1:15 The symbols

Every symbol in those equations is defined here. The two worth
remembering are kappa, the rate at which the variance returns to its
long-run level theta, and xi, the volatility of volatility, which sets
how large the variance fluctuations are.

## 1:30 From model to equation

The fair price is an average, the expected payoff over every path the
market could take, discounted back to today. The Feynman-Kac theorem
says that average, the option price V, solves a deterministic equation.
So instead of simulating the randomness, the program solves the equation
the randomness leaves behind.

## 1:50 The notation

Capital V is the option price, as a function of share price, variance
and time. Lowercase v is the variance, a different quantity. A subscript
on V means a partial derivative.

## 2:00 Why the noise becomes second derivatives

Over one small step the noise moves the price up or down by about root v
times S, both directions equally likely. Average the two and the first
derivatives cancel, whatever the up move gains the down move loses. The
remaining term is half the move size squared times the second
derivative.

## 2:20 The equation, term by term

So every part of the model contributes exactly one term. The drifts give
first derivatives, the noises give second derivatives, the correlation
gives the mixed term, and discounting gives r times V.

## 2:35 The equation the program solves

This is the end point of the derivation. Everything in it except V is a
known number, the rates come from the market and the five Heston
parameters are chosen up front. The strike is not in the equation, it
enters through the expiry payoff. The one unknown is the whole function V, and the
answer is V read at today's share price and variance.

## 3:05 The solver

Share price runs across this grid, variance runs down it, and every cell
holds the option price for that pair. The spacing is the gap between
neighbouring values, about ten dollars between share prices here. Finer
spacing tracks the true curve better, and it multiplies the number of
cells. At expiry each cell equals the payoff, so the solver starts there
and steps backwards to today, rebuilding every cell from its nine
neighbours. A production solve is about a trillion updates of this one
kernel.

## 3:35 Validation

Before benchmarking anything I test correctness. With xi at zero, Heston
reduces to Black-Scholes, which has an exact formula. The solver runs
unchanged and reproduces it with an error around ten to the minus three.
Put-call parity, a model-free identity, holds to five times ten to the
minus six.

## 3:55 Three independent methods

I also priced the same contract three independent ways. The grid gives
196.1684, a Fourier method that shares no code gives 196.1692, which
agrees to four parts in a million, and four hundred thousand Monte Carlo
paths bracket both.

## 4:10 How the optimisation was structured

A single before and after gives one speedup number and no explanation of
where it came from. So there are seven complete solvers, each a copy of
the one below it with exactly one technique added. The difference
between neighbouring levels measures one technique.

## 4:30 The levels that measured nothing

Four of the six techniques measured nothing. Level zero ran the baseline
arithmetic through the ladder's kernel selection, a function pointer
picked once before the time loop, and it measured identical, so the
harness is free. Level one hoisted work that never changes, like the
spacing product two ds, out of the loops, and GCC at O2 had already done
it. Levels three and four walked the grid in memory order and
gave the variance-zero row its own loop, and the baseline was already
written both ways.

## 5:05 Level 2, strength reduction

Level two removes the divisions, five per cell. V sub S was east minus
west, divided by two ds. Now each spacing constant is inverted once per
time step, and the loop multiplies by inv two ds instead. A divide costs
up to fourteen cycles and blocks the divide unit, a multiply costs about
one. The compiler is not allowed to make this change, because a
reciprocal rounds differently in the last bits. Throughput went from
52.6 to 131.2 million cell updates per second.

## 5:40 Level 5, induction variables

Level five replaces the index arithmetic. The code computed row plus
stock i for every neighbour, about ten integer additions per cell. Now
the loop keeps three raw pointers, row mid, row above and row below, and
steps them forward. Worth about two percent, the only measurable gain
after level two.

## 6:00 Level 6, unrolling

Level six unrolls the inner loop by four. The body is copied with
numbered locals, V zero to V three, so four independent chains of
arithmetic sit between each pair of loop tests. Worth half a percent,
the compiler already unrolls counted loops. On my laptop the same change
was a 23 percent regression, which is why every number in this talk
comes from the cluster.

## 6:25 The ladder, measured

Here is the whole ladder measured on the cluster, and only level two
moved. Everything the compiler was allowed to do, it had already done.
The one optimisation that paid is the one it is forbidden to make.
Swapping the loop order on purpose makes identical arithmetic 2.1 times
slower, so the traversal order the baseline already had was doing real
work. Overall the solver is 2.56 times faster, and an automated test
holds every level to the baseline answer, within fifteen units in the
last decimal place.

## 7:00 The benchmark protocol

Every number here comes from one protocol. The job reserves a compute
node exclusively, rebuilds, and passes the full test suite before it
measures anything. Each result is the median of five runs, stamped with
its host, compiler and flags. One control still misled me. GCC 8.5 does
not vectorise at O2 at all, so my no-vectorise sweep proved nothing.
Re-run at O3, the vectoriser does nothing to my kernels but nearly
doubles the baseline.

## 7:30 Scaling with problem size

This sweep grows the working set from one mebibyte to sixty-four, and
the speedup holds between two and a half and three and a half times. The
baseline stays flat the whole way, which is what you expect when the
limit is the divide unit rather than memory.

## 7:50 The roofline

I had assumed this solver was memory-bound. It is not. The level six
kernel does thirty-one floating point operations for every twenty-four
bytes it moves. A separate job measured the node's real scalar ceiling
at 5.1 gigaflops rather than trusting a datasheet. Level six runs at 82
percent of that ceiling, and the baseline at 27, held back by its
divisions.

## 8:15 The page-size effect

The largest single effect was not a code change at all. The same binary
on the same node runs at 134.8 million cell updates per second on
two-mebibyte pages, and 81.1 on four-kibibyte pages, 1.66 times from
page size alone. Sixteen mebibytes is eight huge pages or four thousand
small ones, against a translation buffer with about a thousand entries.
I found it because my median of five was averaging the two regimes.

## 8:45 Reflection

At O2 the compiler has already done the textbook transformations, so the
hand optimisations that pay are the ones it is not allowed to make. The
most surprising result was page size, worth more than five of my six
techniques combined. What I would do differently is reserve the node
exclusively from the first job. A shared node said my loop-order control
cost thirty times, when the real answer was 2.1. The kernel is now
compute-bound, so Milestone 2 is SIMD first, then OpenMP across rows.

## 9:20 Close

I promised to show you how the pricer works, to make it faster, and to
name what I was wrong about. It is 2.56 times faster, it was never
memory-bound, and the compiler had already done five of my six
techniques. Thanks for watching.
