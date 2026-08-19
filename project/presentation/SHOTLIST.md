# SHOTLIST — how each asset is used on the day

TIMELINE.md says when each asset appears. This file says how to shoot it.
Two rules apply to every shot, so they are stated once here rather than
repeated in the table.

The first rule is that the provenance caption along the bottom of every
benchmark figure must never be cropped off. That little strip naming the
host, the CPU, the compiler, the flags and the date is the credibility of
the whole benchmarking section, and the spec's assessors are exactly the
people who will look for it.

The second rule is that the webcam must never cover the content of a slide.
The spec says this in as many words, so the webcam sits in the bottom right
corner by default, and the table below says where to move it on the few
slides whose content lives in that corner. The deck has a reserved webcam
area built in, and pressing the letter p while presenting shows its outline
so the recording software can be lined up against it.

| asset | how to play it | what to say over it | webcam |
|---|---|---|---|
| `demo_timevalue.mp4`, 11 seconds | loop it for the whole opening | the three promises of the video | bottom right, face may be large here |
| `payoff.png` | still | what an option is, pointing at the kink | bottom right |
| `xi_smile.mp4`, 7 seconds | play once, hold the final frame | the flat line is Black-Scholes, the curve is the market | bottom right |
| the Heston slide | still | the model in one sentence, then the two equations | bottom right |
| the randomness slide | still | there is no random number generator anywhere in the program | bottom right |
| `demo_timevalue.mp4` again | loop | shares across, variance down, colour is value, time runs backwards | bottom right |
| `bs_collapse.mp4`, 5 seconds | loop it twice | switch the volatility of volatility off and Heston becomes Black-Scholes | bottom right |
| the three-prices slide | still | three methods that share no code agree | bottom right |
| the ladder-idea slide | still | seven solvers, each one technique apart | bottom right |
| the null-levels table | still | levels 0, 1, 3 and 4 in one pass, one reason per zero | bottom right |
| `shot_L2_strength_reduction.png` | still | the divisions on the left became the multiplications on the right | bottom right |
| `shot_L5_pointers.png` | still | index arithmetic replaced with three row pointers | bottom right |
| `shot_L6_unroll.png` | still | the unroll by four, and the laptop regression | bottom right |
| `bench_ladder_hugepage.png` | still | walk the bars, stop on level 2, then the loop-order control and the test tolerance from the callout | top right, because the legend sits top left |
| `shot_csv_header.png` | still | the job stamps its own provenance into every result file | bottom right |
| `bench_scaling.png` | still | the shaded band is the page-regime split, and the first run of this sweep had to be redone | bottom right |
| `roofline.png` | still | the measured ceiling, and the memory-bound assumption disproved | bottom left, because the label column fills the right side |
| the page-size slide | still | the same binary runs two-thirds faster on nothing but page size | bottom right |
| the reflection slide | still | the four reflection questions, one sentence each | face large, centre right |
| the close | none | return to the three promises | face fills the screen |

## The recording day, in order

Record at 1080p or better, because the code screenshots are two thousand
pixels wide and only stay readable if the recording keeps that detail. Do a
microphone test with the actual first sentence and listen back to it before
recording anything long. Keep a visible countdown timer running, aim for
nine minutes, and remember the assessors stop watching at ten regardless of
what is mid-sentence.

Record one take per row of the timeline and number the takes to match, so
assembling the video is splicing rather than editing. The spec explicitly
allows splices and tolerates the face jumping between them. The face must be
visible in every take, including the ones that are mostly slides, because
identity verification is a pass or fail gate.

After exporting, play the file end to end before anything else. It must be
H.264, under one hundred megabytes, audible, and named exactly
`Milestone1_48008361.mp4`. Then `make package` builds the submission zip. It
refuses to run if the video is missing, and it exists because a plain zip of
the directory would ship the credentials file and half a gigabyte of
rendered frames. Before zipping, confirm that a fresh clone builds and
passes with `make clean && make && make test`, and that every figure used in
the video can be regenerated with `make bench-plot roofline code-shots
payoff convergence`.
