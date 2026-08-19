# PRESENTATION.md — the Milestone 1 video plan

The 10-minute video is the entire assessed deliverable. The code is pass/fail
and already passes. This file is the brief: what the markers reward, how the
time is spent, what appears on screen for each beat, and what to say.

The recording side of this plan lives in the `presentation/` directory.
TIMELINE.md turns the beats below into a timestamped table, SCRIPT.md is
the word-for-word script with one section per take, SHOTLIST.md
says how each asset is shot and where the webcam sits, SPEAKER_NOTES.md
scripts the model-to-equation run and the optimisation slides in plain
English so they can be delivered and defended without re-deriving them,
and deck.tex is the
slide deck to record over (`make slides` builds deck.pdf via tectonic;
deck.html is the earlier HTML version, kept for its on-screen take timer).
Slides that carry an animation show a still frame from `make stills`, and
the MP4 itself is spliced in fullscreen during the edit, as SHOTLIST.md
prescribes. This file decides what to say. Those files decide how it gets
onto the screen.

---

## 1. What is actually being marked

From the milestone specification, page 2. These five are the graded criteria
and the three below them are pass/fail gates.

| criterion | weight | what the spec asks for |
|---|---|---|
| Introduction and Background | 20% | bring the assessor up to speed on the model being simulated. Must be **self-contained**, because the assessor will not remember Milestone 1 when they watch Milestone 2. |
| Optimisation | 25% | optimisations attempted and their performance, **techniques tried that failed or were considered and rejected**, how the problem was structured for speed, and the **trade-offs** made. |
| Benchmarking | 25% | a setup that avoids inaccurate or misleading results, **hotspot and bottleneck identification**, and **how performance scales with problem size**. Must run on a UQ cluster. (Amdahl and Gustafson are flagged for Milestone 2, not this one.) |
| Presentation | 20% | choice of topics and time given to each, verbal delivery, slide design, data visualisation. |
| Reflection / Conclusion | 10% | main conclusions, what was learnt or most surprising, what you would do differently, how you would extend it, and **plans for Milestone 2**. |

| gate | requirement |
|---|---|
| Identity verification | **your face visible throughout.** Splices are fine and abrupt head movement across them is tolerated. |
| Code submission | well-commented source plus Makefile and slurm scripts. Not graded directly, used to verify the model was really implemented. |
| Submission format | H.264, **under 100 MB** (Blackboard's limit, so about 1 Mbit/s), named `Milestone1_48008361.mp4`, and code plus video in one zip named `Milestone1_48008361.zip`. |

Two things follow directly from that table and are worth internalising.

**Optimisation and Benchmarking together are half the mark.** They are worth
more than twice the Introduction. Do not spend four minutes explaining the
Heston model.

**Null results are explicitly rewarded.** The spec asks for techniques that
failed or were rejected. Five of the six ladder rungs measured as nothing, and
two beliefs were disproved by measurement. That is not a weakness to hide. It is what the spec asks for.

---

## 2. Time budget

Target **9 minutes**, hard ceiling 10. Assessors stop watching at 10:00 and
the spec says a shorter high-quality presentation beats a longer low-quality
one. Milestone 1 in particular says shorter is acceptable.

| beat | time | running |
|---|---|---|
| Cold open and framing | 0:20 | 0:20 |
| Introduction and Background | 2:50 | 3:10 |
| Validation (bridges into benchmarking) | 0:55 | 4:05 |
| Optimisation | 3:05 | 7:10 |
| Benchmarking | 2:00 | 9:10 |
| Reflection and Milestone 2 | 0:35 | 9:45 |
| Close | 0:10 | 9:55 |

The Introduction is longer than the first draft of this plan allowed,
because the model-to-equation run now builds the pricing equation slide by
slide rather than asserting it. The spec permits this. It says some models
naturally need more background than others. The extra time comes out of the
cold open, the reflection and the slack, not out of Optimisation and
Benchmarking, which together still hold half the talk.

---

## 3. Beat sheet

### 0:00 — Cold open (face full screen, then the animation behind you)

Open on the weather map already moving. Say what the video will deliver, in
three promises, so the assessor knows the shape up front.

> "This is the price of a stock option, computed by solving a partial
> differential equation on a grid. In the next nine minutes I will show you how
> it works, how I made it two and a half times faster, and the two things I
> was confident about that turned out to be wrong when I measured them."

The screen shows `results/demo_timevalue.mp4`.

### 0:20 — Introduction and Background (20%)

Five short beats. Do not linger. The assessor needs enough to follow the rest,
not a finance course.

1. **The contract.** An option is the right, not the obligation, to buy a share
   at a fixed price on a fixed date. Show the payoff kink.
2. **Why the standard answer is not enough.** Black-Scholes assumes the
   share's volatility is one constant number forever. Markets disagree, and the
   evidence is the volatility smile. Play `results/xi_smile.mp4` and point at
   the curve. A flat line is what Black-Scholes predicts. The curve is what the
   market actually prices.
3. **Heston in one sentence.** The share price is volatile, the volatility
   is itself stochastic, and the two are negatively correlated because falling
   markets tend to be more volatile. The slide after the equations defines
   every symbol in words. Walk it in one pass so nothing on screen goes
   undefined. The slide after that says where the strike price enters. The
   equations describe the share and know nothing about the contract, and
   the strike enters through the payoff at expiry. Say this plainly,
   because a viewer used to the Black-Scholes formula will be looking
   for K in the equations and not finding it.
4. **Where the randomness went.** There is no random number generator
   anywhere in the program, and a five-slide run explains why, ending on the
   full equation. The Feynman-Kac slide states the idea, the fair price is an
   average over every path the market could take, and that average solves a
   deterministic equation. The notation slide then defines V and its
   subscripts before they are used. The curvature slide shows the one real step,
   averaging over an up move and a down move cancels the first derivatives
   and leaves half the move size squared times the second derivative. The term-by-term slide maps
   each part of the model to its term. The final slide assembles the whole
   equation and separates the known numbers from the one unknown function,
   which is the moment the viewer sees what the program actually solves.
   Walk the table slides at reading pace and keep the voice brief, because
   the rows carry the detail.
5. **The spreadsheet.** Columns are share prices, rows are variance levels, and
   the sheet is filled in backwards from expiry, where the answer is known
   exactly because there is no time left for anything to happen. Show the
   weather map again and let it run backwards.

**An optional ten seconds of extra depth, only if time allows.** Along the
bottom row the variance is zero, every diffusion term vanishes, and what is
left is a transport equation. Solving that is better than imposing a value
by hand, and the slope has to be one-sided because the drift points into
the grid.

The screen shows the payoff diagram, then `xi_smile.mp4`, then
`demo_timevalue.mp4`.

### 3:10 — Validation (0:55)

This beat is short but it buys credibility for everything after it, and it is
scored under Benchmarking because the spec asks for a setup that avoids
misleading results.

- **Turn the model off and it becomes the model we can check.** Set the
  volatility of volatility to zero and Heston collapses to Black-Scholes,
  which has an exact formula. Play `results/bs_collapse.mp4`.
- **Three independent methods, one number.** The grid says 196.1684. Fourier
  inversion says 196.1692, which agrees to four parts in a million. Four
  hundred thousand Monte Carlo paths say 196.1345 give or take 0.5524, which
  brackets both. Put all three on one slide.
- **Convergence.** `results/convergence.png`, and say the observed order out
  loud.
- **The test gate.** Show the four lines of `slurm/bench_common.sh` that refuse
  to benchmark a build whose tests have not just passed. Say plainly that no
  optimised version was ever timed until it had re-proved it gives the
  reference answer.

The screen shows `bs_collapse.mp4`, then the three-method table, then
`convergence.png`, then the `build_and_validate` snippet.

### 4:05 — Optimisation (25%)

Lead with the method, not the result.

1. **Why not a before and after.** One slow version and one fast version gives
   a ratio and tells you nothing about which change mattered. So there are
   seven complete solvers, each a copy of the one below with exactly one
   technique added and nothing else changed.
2. **Walk the ladder.** `results/bench_ladder_hugepage.png`, bar by bar.
   Baseline, then hoisting and common subexpression elimination, then strength
   reduction, then traversal order, loop splitting, induction variables,
   unrolling.
3. **The result.** Only one bar moved. Level 2 took it from 52.6 to 131.2
   million cell updates per second and every other rung came back inside the
   noise. Then say why. GCC at -O2 had already done all the others, and
   level 2 is the only transformation on the list the compiler is
   **forbidden** to make, because turning a division into a multiplication
   by a reciprocal changes the answer in the last bits.
   > "At -O2, the hand optimisations worth doing are the ones the compiler is
   > not allowed to do for you."
4. **The negative controls.** Levels 3 and 4 were null by construction, because
   the baseline already had the right loop order and already peeled the special
   row out. So instead of pretending, write the wrong version and time it.
   Swapping the loops costs 2.1 times. Fusing the special row back in is free,
   because the branch predictor handles it. Play `results/memory_cache.mp4`
   over the loop-order result.
5. **Trade-offs, said out loud.** The speedup cost 15 units in the last place
   of the answer, which is 3.3e-15 relative, and `test_opt_matches` holds every
   level to that. The kernel file is long and repetitive on purpose, because
   sharing setup between rungs would let one rung change what its neighbour
   measures.

The screen shows `bench_ladder_hugepage.png`, then `memory_cache.mp4`, then
the `test_opt_matches` output.

### 7:10 — Benchmarking (25%)

Hit the three things the spec names, in order.

1. **A setup that cannot mislead.** One slide, and it is the CSV header. Host,
   date, job id, compiler, exact flags, rep count, CPU model, all four cache
   sizes, stamped into every result file by the job itself. Then the hygiene:
   an sbatch job on a compute node, `--exclusive` so no co-tenant steals
   bandwidth, `--nodelist=r730-2` because the partition mixes 24-core Xeons
   with 8-core virtual machines, median of five with the min and max shown.
2. **Scaling with problem size.** `results/bench_scaling.png`, four grids
   from a working set of one mebibyte up to sixty-four. Two sentences
   earn real marks here. The first is that the first run of this sweep
   confounded page size with cache and was re-run with the allocator pinned.
   The second is that the baseline's flat line across a sixty-four-fold
   growth in working set confirms the division-bound story from another
   direction, because only the optimised kernel is fast enough to feel the
   memory system at all.
3. **Hotspot and bottleneck.** `results/roofline.png`. State the assumption first
   and then show the measurement. The solver was assumed to be memory-bound. It is not.
   Counting arithmetic and bytes straight off the level 6 kernel gives 31 flops
   per 24 bytes, and a separate job measured this node's real scalar ceiling at
   5.117 GFLOP/s rather than quoting a datasheet. Level 6 hits 82 percent of
   that while asking for 3 GB/s on a four-channel machine. The baseline sits at
   27 percent and is limited by five divisions per cell. So the optimisation
   moved the program from division-bound to compute-bound.
4. **The page-size discovery.** Same binary, same node, same flags. 134.8 or
   81.1 million cell updates per second depending on nothing but whether the
   two 8 MiB buffers sit on 2 MiB pages or 4 KiB pages, because 16 MiB needs
   either 8 address translation entries or 4096. Explain how it was caught:
   the median of five was silently averaging two different memory regimes,
   because the allocator happened to give the first repetition huge pages and
   the later ones small pages. It was established by pinning the allocator's
   mmap threshold by hand and watching the effect appear and disappear on
   demand.

The screen shows a CSV header, then `bench_scaling.png`, then
`roofline.png`, then the page-size table from RESULTS.md.

### 9:10 — Reflection and Conclusion (10%)

Four questions, one sentence each. The spec names all four.

- **Main conclusion.** At -O2 the compiler has already done the textbook
  transformations, so the hand optimisations that pay are the ones it is
  forbidden to do.
- **Most surprising.** Page size was worth more than five of the six techniques
  put together, and it was not a code change at all.
- **What I would do differently.** Run `--exclusive` from the very first job.
  An early measurement said the loop-order control cost 30 times, and it was
  actually 2.1, because the machine was shared with my own reference solve. It
  was a real measurement that answered a different question than the one I
  asked.
- **Milestone 2.** The roofline says compute-bound, so the next win has to come
  from doing more arithmetic per instruction, which is exactly what SIMD is.
  Then OpenMP on the outer loop, which needs zero synchronisation inside a
  timestep because the two-buffer design makes every cell independent of every
  other cell.

### 9:45 — Close

Name the three promises from the cold open and confirm you delivered them.

---

## 4. Things that must NOT be said

- **Do not quote 2295 seconds, or "about 40 minutes", as the reference solve
  runtime.** That is a laptop measurement and it sits in the same tables as
  cluster numbers. If the full-grid runtime is wanted, derive it from measured
  rangpur throughput and say you are doing so: about 2.1 hours at level 6 and
  about 5.4 hours for the baseline. See RESULTS.md.
- Do not present the benchmark configuration's price column as a price. Its
  maturity is four and a half hours, so the number is meaningless by design and
  only the throughput matters. `config/bench.cfg` says so in its first
  paragraph.
- Do not call a null rung a failure. The spec rewards it.
- **Do not claim sweep C ruled out the auto-vectoriser.** GCC 8.5 does not
  vectorise at -O2 at all, so that control disabled nothing, and the 1.00
  ratio it reported means nothing. Sweeps F and G measured the real answer,
  which is that the vectoriser does nothing to the optimised kernels and
  nearly doubles the baseline. If the control comes up at all, tell it as
  the control that was itself caught being misleading and then fixed. That
  version of the story earns Benchmarking marks, and the false version
  loses them to any assessor who knows GCC.

---

## 5. Production checklist

- [ ] Face visible throughout. Picture-in-picture webcam in a corner is the
      simplest way to satisfy the identity gate while slides carry the content.
- [ ] Record in takes and splice. Explicitly allowed.
- [ ] Export H.264. Target roughly 1 Mbit/s so 9 minutes lands well under
      100 MB.
- [ ] Name it `Milestone1_48008361.mp4`. Confirm that number is yours before
      exporting, since the spec gives it as the worked example.
- [ ] Zip video plus code as `Milestone1_48008361.zip`.
- [ ] Confirm `make clean && make && make test` works from a fresh clone before
      zipping, which is the last unticked box in PLAN.md §8.
- [ ] Every figure in the video carries its provenance caption. `bench_plot.py`
      already renders host, CPU, compiler, flags and date into the corner of
      each one, so do not crop it off.
- [ ] Place the webcam as SHOTLIST.md says for each shot, so the face never
      covers figure content. The spec warns about exactly this, in the words
      "be careful not to obscure important parts of your slides with your
      face."
- [ ] Keep a visible timer while recording. The target is nine minutes and
      the assessors stop watching at ten.
- [ ] Play the exported file back before zipping it, and check that it is
      H.264, that it is under one hundred megabytes, that the audio is
      audible, and that it is named exactly `Milestone1_48008361.mp4`.

---

## 6. Asset index

| asset | where it appears |
|---|---|
| `results/demo_timevalue.mp4` | the opening, and again for the backwards-from-expiry explanation |
| `results/payoff.png` (built by `make payoff`) | the introduction, showing what an option is |
| `results/xi_smile.mp4` | the introduction, showing why Black-Scholes is not enough |
| `results/bs_collapse.mp4` | validation |
| `results/convergence.png` | validation |
| `results/code_shots/*.png` (built by `make code-shots`) | the test gate, the rung slides, the loop-order control, the CSV provenance header and the test output |
| `results/bench_ladder_hugepage.png` | the optimisation ladder |
| `results/memory_cache.mp4` | the loop-order control |
| `results/bench_scaling.png` | scaling with problem size |
| `results/roofline.png` | finding the bottleneck |
| `results/demo_surface.mp4` | not used, because it duplicates what the flat map already shows |
| `results/blowup.mp4` | not used, kept as a spare in case a recording comes in short |
| `results/memory_buffers.mp4` | not used, because the point it makes fits in one spoken sentence |

TIMELINE.md records each of those cut decisions with its reason.
