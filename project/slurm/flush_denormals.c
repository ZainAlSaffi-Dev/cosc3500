/* AI declaration
   I used Claude, an AI assistant made by Anthropic, while working on this file.
   It helped me debug the code, draft pseudocode and implement it. */

/* Switches on flush to zero and denormals are zero when a program starts.
   Loaded into heston with LD_PRELOAD by slurm/bench_common.sh, heston itself is unchanged. */

#include <xmmintrin.h>

/* Runs before main. Bit 15 (0x8000) is FTZ, results too small to be a normal double become 0.
   Bit 6 (0x0040) is DAZ, denormal inputs are read as 0. Threads made later copy this setting */
__attribute__((constructor)) static void flush_denormals(void) {
    _mm_setcsr(_mm_getcsr() | 0x8040);
}
