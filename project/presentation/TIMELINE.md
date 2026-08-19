# TIMELINE — the minute-by-minute plan for the video

This file turns the beat sheet in PRESENTATION.md into timestamps and
concrete assets. The beat sheet decides what to say and why it earns marks.
This file is the thing to keep next to the recording software.

The slides are `presentation/deck.pdf`, built from `deck.tex` by
`make slides`. There is one slide for each row of the table below, in the
same order. The backup slides after the closing frame are interview
material and are not recorded. `deck.html` is the earlier HTML deck. Its
take timer is still useful, but its slide order no longer matches this
table, so record over the PDF.

The table ends at nine minutes thirty-five. The assessors stop watching at
ten, so takes must not drift. If the recording runs long, the four null
rungs, levels 0, 1, 3 and 4, are the place to tighten first.

| start | length | what happens | what is on screen | camera |
|---|---|---|---|---|
| 0:00 | 0:20 | The opening. Say the three promises of the video. | `results/demo_timevalue.mp4` looping behind the title | face large in the corner |
| 0:20 | 0:20 | What an option is. | `results/payoff.png` | small webcam |
| 0:40 | 0:25 | Why Black-Scholes is not enough. The market prices a smile. | `results/xi_smile.mp4`, then hold its final frame | small webcam |
| 1:05 | 0:15 | Heston in one sentence, and its two equations. | a text slide | small webcam |
| 1:20 | 0:15 | Every symbol in those equations, walked in one pass. | a table slide | small webcam |
| 1:35 | 0:10 | Where the strike price enters. The model and the contract are separate. | a text slide | small webcam |
| 1:45 | 0:10 | From stochastic model to deterministic equation. There is no random number generator in the program. | a text slide | small webcam |
| 1:55 | 0:10 | The derivative notation. Subscripts denote partial derivatives. | a table slide | small webcam |
| 2:05 | 0:15 | Why the noise terms become second derivatives. | a text slide | small webcam |
| 2:20 | 0:10 | The equation term by term. Each part of the model contributes one term. | a table slide | small webcam |
| 2:30 | 0:15 | The full equation, known numbers separated from the one unknown. | a text slide | small webcam |
| 2:45 | 0:25 | The grid, filled in backwards from expiry. | `results/demo_timevalue.mp4` again, this time explained | small webcam |
| 3:10 | 0:20 | Validation begins. Switch the model off and it becomes one we can check exactly. | `results/bs_collapse.mp4` | small webcam |
| 3:30 | 0:15 | Three methods that share no code agree on one number. | a table slide with the three prices | small webcam |
| 3:45 | 0:15 | Why the optimisation study is a ladder rather than a before and after. | a text slide | small webcam |
| 4:00 | 0:15 | Level 0, the harness control. It confirms the harness costs nothing. | `results/code_shots/shot_L0_anchor.png` | small webcam |
| 4:15 | 0:15 | Level 1, hoisting. The compiler had already done it. | `results/code_shots/shot_L1_hoisting.png` | small webcam |
| 4:30 | 0:25 | Level 2, strength reduction. The change the compiler is not allowed to make. | `results/code_shots/shot_L2_strength_reduction.png` | small webcam |
| 4:55 | 0:15 | Level 3, traversal order. Null by construction. | `results/code_shots/shot_L3_order.png` | small webcam |
| 5:10 | 0:15 | Level 4, loop splitting. Also null by construction. | `results/code_shots/shot_L4_split.png` | small webcam |
| 5:25 | 0:15 | Level 5, induction variables. The one small real win after level 2. | `results/code_shots/shot_L5_pointers.png` | small webcam |
| 5:40 | 0:15 | Level 6, unrolling. The compiler already unrolls counted loops. | `results/code_shots/shot_L6_unroll.png` | small webcam |
| 5:55 | 0:20 | The ladder, measured. Only level 2 moved, and why. | `results/bench_ladder_hugepage.png` | small webcam |
| 6:15 | 0:20 | The negative controls. Each removes a property the baseline already had. | `results/code_shots/shot_ctl_order.png`, then `results/memory_cache.mp4` | small webcam |
| 6:35 | 0:15 | What the speedup cost. | `results/code_shots/shot_optmatches.png` | small webcam |
| 6:50 | 0:30 | How the benchmarks were set up so they could not mislead. | the protocol table, with `results/code_shots/shot_csv_header.png` spliced fullscreen | small webcam |
| 7:20 | 0:25 | How performance scales with problem size, including the page-regime split. | `results/bench_scaling.png` | small webcam |
| 7:45 | 0:30 | The bottleneck. The solver was assumed to be memory-bound and is not. | `results/roofline.png` | small webcam |
| 8:15 | 0:35 | The page-size effect. | a table slide with the page-size numbers | small webcam |
| 8:50 | 0:35 | Reflection. The four questions the spec asks, one sentence each. | a text slide | face large |
| 9:25 | 0:10 | Close by returning to the promises from the opening. | face only | face full screen |

Somewhere inside the benchmarking section, spend about fifteen seconds on
the control that went wrong. One of my own
controls turned out to be vacuous, because GCC 8.5 does not vectorise at
-O2 at all, so my no-vectorise sweep disabled nothing and proved nothing. I
re-ran the control properly at -O3, and the answer split in two. The
vectoriser does nothing to my optimised kernels, but it nearly doubles the
baseline, because division-heavy code is exactly what it likes. My speedup
is real, and so was my broken control. That story is worth telling because
the spec explicitly rewards a benchmarking setup that avoids misleading
results, and catching your own control is the strongest possible evidence
you look for them. The numbers behind it are in RESULTS.md section 6.

## What was cut, and why

| asset | reason |
|---|---|
| `results/convergence.png` | The convergence slide moved to the unrecorded backup section. The result stays in RESULTS.md and the code ships the study. |
| `results/code_shots/shot_iron_rule.png` | The test-gate slide came out. Its content is one row of the protocol table and one spoken sentence in that take. |
| `results/demo_surface.mp4` | It does the same job as the flat weather map, and the three-dimensional view is slower for a viewer to parse in the time available. |
| `results/memory_buffers.mp4` | The two-buffer idea lands in one spoken sentence, so fifteen seconds of animation would not buy fifteen seconds of marks. |
| `results/blowup.mp4` | The instability demonstration is an aside. PRESENTATION.md keeps it listed as a spare in case a recording comes in short. |
| the older `weather_map` and `surface` renders | They are earlier versions of what `demo_timevalue.mp4` now does better. |
| every `.gif` file | The GIF files exist as twins of the MP4 files for embedding in documents. The video always uses the MP4. |

The rule behind every cut is that each second on screen has to serve one of
the five marked criteria, and the spec says plainly that a shorter
high-quality presentation beats a longer low-quality one.
