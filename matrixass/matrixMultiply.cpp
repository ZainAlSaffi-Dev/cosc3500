// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

// Provided header, brings in floatType, AVX and OpenMP
#include <matrixMultiply.h>
#define STUDENTID 48008361 //DO NOT REMOVE

// Defaults used when no tuning values are passed.
// k values per block
static const int kBlockSizeDefault = 384;
// Rows per block
static const int rowBlockSizeDefault = 96;
// 1 uses three real multiplies per complex multiply, 0 uses four
static const int useThreeMultDefault = 1;

// Buffers kept between calls. Inputs are split into real and imaginary arrays for the AVX code.
// N the buffers were sized for, -1 means none yet
static int scratchN = -1;
// k block size the pack buffers were sized for
static int scratchKBlockSize = 0;
// Real parts of the left matrix
static float* leftRe = 0;
// Imaginary parts of the left matrix
static float* leftIm = 0;
// Real plus imaginary
static float* leftSum = 0;
// Real parts of the right matrix
static float* rightRe = 0;
// Imaginary parts of the right matrix
static float* rightIm = 0;
// Real parts of the output
static float* outRe = 0;
// Imaginary parts of the output
static float* outIm = 0;
// One k block of the right matrix
static float* packRe = 0;
// Imaginary parts of the packed block
static float* packIm = 0;
// Real plus imaginary of the packed block
static float* packSum = 0;

// Rows are padded to a multiple of 4 and columns to a multiple of 16 so there are no edge cases.
// Only reallocates when N changes or the k block grows
static void ensureScratch(int N, int paddedRows, int paddedCols, int kBlockSize)
{
    // Already big enough, reuse them
    if (N == scratchN && kBlockSize <= scratchKBlockSize) { return; }
    // Free the old left buffers, freeing null is fine
    _mm_free(leftRe); _mm_free(leftIm); _mm_free(leftSum);
    // Free the old right buffers
    _mm_free(rightRe); _mm_free(rightIm);
    // Free the old output buffers
    _mm_free(outRe); _mm_free(outIm);
    // Free the old pack buffers
    _mm_free(packRe); _mm_free(packIm); _mm_free(packSum);
    // Never make the pack buffer smaller than the default
    const int packBlockSize = (kBlockSize > kBlockSizeDefault) ? kBlockSize : kBlockSizeDefault;
    // Floats in a padded left array
    const size_t leftFloats = (size_t)paddedRows * N;
    // Floats in a padded right array
    const size_t rightFloats = (size_t)N * paddedCols;
    // Floats in a padded output array
    const size_t outFloats = (size_t)paddedRows * paddedCols;
    // Floats in one packed k block
    const size_t packFloats = (size_t)packBlockSize * paddedCols;
    // Left real parts, 64 byte aligned so aligned AVX loads are safe
    leftRe = (float*)_mm_malloc(sizeof(float) * leftFloats, ALIGN);
    // Left imaginary parts
    leftIm = (float*)_mm_malloc(sizeof(float) * leftFloats, ALIGN);
    // Left sums
    leftSum = (float*)_mm_malloc(sizeof(float) * leftFloats, ALIGN);
    // Right real parts
    rightRe = (float*)_mm_malloc(sizeof(float) * rightFloats, ALIGN);
    // Right imaginary parts
    rightIm = (float*)_mm_malloc(sizeof(float) * rightFloats, ALIGN);
    // Output real parts
    outRe = (float*)_mm_malloc(sizeof(float) * outFloats, ALIGN);
    // Output imaginary parts
    outIm = (float*)_mm_malloc(sizeof(float) * outFloats, ALIGN);
    // Packed real parts
    packRe = (float*)_mm_malloc(sizeof(float) * packFloats, ALIGN);
    // Packed imaginary parts
    packIm = (float*)_mm_malloc(sizeof(float) * packFloats, ALIGN);
    // Packed sums
    packSum = (float*)_mm_malloc(sizeof(float) * packFloats, ALIGN);
    // The padding is zeroed once here and not written afterwards.
    // Zero left real, including the padding
    memset(leftRe, 0, sizeof(float) * leftFloats);
    // Zero left imaginary
    memset(leftIm, 0, sizeof(float) * leftFloats);
    // Zero left sums
    memset(leftSum, 0, sizeof(float) * leftFloats);
    // Zero right real, padded columns stay 0
    memset(rightRe, 0, sizeof(float) * rightFloats);
    // Zero right imaginary
    memset(rightIm, 0, sizeof(float) * rightFloats);
    // Remember the size we allocated for
    scratchN = N;
    // And the k block size
    scratchKBlockSize = packBlockSize;
}

// acc + x*y or acc - x*y, SIGN is fixed at compile time. There is no fused multiply add under -mavx.
// SIGN is known at compile time so there is no branch
template<int SIGN>
// __m256 holds 8 floats
static inline __m256 mulAddSigned(__m256 acc, __m256 x, __m256 y)
{
    // Add the product when SIGN is positive, otherwise subtract it
    return (SIGN > 0) ? _mm256_add_ps(acc, _mm256_mul_ps(x, y))
                      : _mm256_sub_ps(acc, _mm256_mul_ps(x, y));
}

// dst + t when SIGN is positive, dst - t otherwise.
// Same compile time sign trick
template<int SIGN>
// Used when adding a tile result into the output
static inline __m256 addSigned(__m256 dst, __m256 t)
{
    // Picked at compile time
    return (SIGN > 0) ? _mm256_add_ps(dst, t) : _mm256_sub_ps(dst, t);
}

// out += leftA * packA + SIGN * leftB * packB for a 4 by 16 block of the output.
// leftA and leftB point to four consecutive rows, leftStride apart. packA and packB hold 16 floats per k.
// SIGN is -1 for the real part and +1 for the imaginary part
template<int SIGN>
// Left inputs, packed right blocks, k count and the output tile
static void tile4x16TwoProducts(const float* leftA, const float* leftB, long leftStride,
                                const float* packA, const float* packB, int kCount,
                                float* out, long outStride)
{
    // accRH is output row R, half H of the 16 columns
    __m256 acc00 = _mm256_load_ps(out);
    // Row 0, columns 8 to 15
    __m256 acc01 = _mm256_load_ps(out + 8);
    // Row 1, columns 0 to 7
    __m256 acc10 = _mm256_load_ps(out + outStride);
    // Row 1, columns 8 to 15
    __m256 acc11 = _mm256_load_ps(out + outStride + 8);
    // Row 2, columns 0 to 7
    __m256 acc20 = _mm256_load_ps(out + 2 * outStride);
    // Row 2, columns 8 to 15
    __m256 acc21 = _mm256_load_ps(out + 2 * outStride + 8);
    // Row 3, columns 0 to 7
    __m256 acc30 = _mm256_load_ps(out + 3 * outStride);
    // Row 3, columns 8 to 15
    __m256 acc31 = _mm256_load_ps(out + 3 * outStride + 8);
    // Walk through this k block
    for (int k = 0; k < kCount; ++k) {
        // packA columns 0 to 7 at this k
        const __m256 a0 = _mm256_load_ps(packA + 16 * k);
        // packA columns 8 to 15
        const __m256 a1 = _mm256_load_ps(packA + 16 * k + 8);
        // packB columns 0 to 7
        const __m256 b0 = _mm256_load_ps(packB + 16 * k);
        // packB columns 8 to 15
        const __m256 b1 = _mm256_load_ps(packB + 16 * k + 8);

        // Row 0 leftA value copied into all 8 lanes
        __m256 l = _mm256_broadcast_ss(leftA + k);
        // Row 0 += leftA * packA
        acc00 = mulAddSigned<1>(acc00, l, a0);
        // Same for columns 8 to 15
        acc01 = mulAddSigned<1>(acc01, l, a1);
        // Row 0 leftB value
        l = _mm256_broadcast_ss(leftB + k);
        // Row 0 += SIGN * leftB * packB
        acc00 = mulAddSigned<SIGN>(acc00, l, b0);
        // Same for columns 8 to 15
        acc01 = mulAddSigned<SIGN>(acc01, l, b1);

        // Row 1 leftA value
        l = _mm256_broadcast_ss(leftA + leftStride + k);
        // Row 1 += leftA * packA
        acc10 = mulAddSigned<1>(acc10, l, a0);
        // Same for columns 8 to 15
        acc11 = mulAddSigned<1>(acc11, l, a1);
        // Row 1 leftB value
        l = _mm256_broadcast_ss(leftB + leftStride + k);
        // Row 1 += SIGN * leftB * packB
        acc10 = mulAddSigned<SIGN>(acc10, l, b0);
        // Same for columns 8 to 15
        acc11 = mulAddSigned<SIGN>(acc11, l, b1);

        // Row 2 leftA value
        l = _mm256_broadcast_ss(leftA + 2 * leftStride + k);
        // Row 2 += leftA * packA
        acc20 = mulAddSigned<1>(acc20, l, a0);
        // Same for columns 8 to 15
        acc21 = mulAddSigned<1>(acc21, l, a1);
        // Row 2 leftB value
        l = _mm256_broadcast_ss(leftB + 2 * leftStride + k);
        // Row 2 += SIGN * leftB * packB
        acc20 = mulAddSigned<SIGN>(acc20, l, b0);
        // Same for columns 8 to 15
        acc21 = mulAddSigned<SIGN>(acc21, l, b1);

        // Row 3 leftA value
        l = _mm256_broadcast_ss(leftA + 3 * leftStride + k);
        // Row 3 += leftA * packA
        acc30 = mulAddSigned<1>(acc30, l, a0);
        // Same for columns 8 to 15
        acc31 = mulAddSigned<1>(acc31, l, a1);
        // Row 3 leftB value
        l = _mm256_broadcast_ss(leftB + 3 * leftStride + k);
        // Row 3 += SIGN * leftB * packB
        acc30 = mulAddSigned<SIGN>(acc30, l, b0);
        // Same for columns 8 to 15
        acc31 = mulAddSigned<SIGN>(acc31, l, b1);
    }
    // Write row 0, columns 0 to 7 back
    _mm256_store_ps(out, acc00);
    // Row 0, columns 8 to 15
    _mm256_store_ps(out + 8, acc01);
    // Row 1, columns 0 to 7
    _mm256_store_ps(out + outStride, acc10);
    // Row 1, columns 8 to 15
    _mm256_store_ps(out + outStride + 8, acc11);
    // Row 2, columns 0 to 7
    _mm256_store_ps(out + 2 * outStride, acc20);
    // Row 2, columns 8 to 15
    _mm256_store_ps(out + 2 * outStride + 8, acc21);
    // Row 3, columns 0 to 7
    _mm256_store_ps(out + 3 * outStride, acc30);
    // Row 3, columns 8 to 15
    _mm256_store_ps(out + 3 * outStride + 8, acc31);
}

// t = left * pack for a 4 by 16 block, then out0 += SIGN0 * t and out1 += SIGN1 * t (skipped when SIGN1 is 0).
// Sign for each output, 0 means skip it
template<int SIGN0, int SIGN1>
// One left input, one packed right block and two outputs with the same stride
static void tile4x16OneProduct(const float* left, long leftStride,
                               const float* pack, int kCount,
                               float* out0, float* out1, long outStride)
{
    // Row 0, columns 0 to 7, starts at zero
    __m256 acc00 = _mm256_setzero_ps();
    // Row 0, columns 8 to 15
    __m256 acc01 = _mm256_setzero_ps();
    // Row 1, columns 0 to 7
    __m256 acc10 = _mm256_setzero_ps();
    // Row 1, columns 8 to 15
    __m256 acc11 = _mm256_setzero_ps();
    // Row 2, columns 0 to 7
    __m256 acc20 = _mm256_setzero_ps();
    // Row 2, columns 8 to 15
    __m256 acc21 = _mm256_setzero_ps();
    // Row 3, columns 0 to 7
    __m256 acc30 = _mm256_setzero_ps();
    // Row 3, columns 8 to 15
    __m256 acc31 = _mm256_setzero_ps();
    // Walk through this k block
    for (int k = 0; k < kCount; ++k) {
        // Packed columns 0 to 7 at this k
        const __m256 p0 = _mm256_load_ps(pack + 16 * k);
        // Packed columns 8 to 15
        const __m256 p1 = _mm256_load_ps(pack + 16 * k + 8);
        // Row 0 left value copied into all 8 lanes
        __m256 l = _mm256_broadcast_ss(left + k);
        // Row 0 += left * pack
        acc00 = mulAddSigned<1>(acc00, l, p0);
        // Same for columns 8 to 15
        acc01 = mulAddSigned<1>(acc01, l, p1);
        // Row 1 left value
        l = _mm256_broadcast_ss(left + leftStride + k);
        // Row 1 += left * pack
        acc10 = mulAddSigned<1>(acc10, l, p0);
        // Same for columns 8 to 15
        acc11 = mulAddSigned<1>(acc11, l, p1);
        // Row 2 left value
        l = _mm256_broadcast_ss(left + 2 * leftStride + k);
        // Row 2 += left * pack
        acc20 = mulAddSigned<1>(acc20, l, p0);
        // Same for columns 8 to 15
        acc21 = mulAddSigned<1>(acc21, l, p1);
        // Row 3 left value
        l = _mm256_broadcast_ss(left + 3 * leftStride + k);
        // Row 3 += left * pack
        acc30 = mulAddSigned<1>(acc30, l, p0);
        // Same for columns 8 to 15
        acc31 = mulAddSigned<1>(acc31, l, p1);
    }
    // out0 row 0 += SIGN0 * t
    _mm256_store_ps(out0, addSigned<SIGN0>(_mm256_load_ps(out0), acc00));
    // Row 0, columns 8 to 15
    _mm256_store_ps(out0 + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + 8), acc01));
    // Row 1, columns 0 to 7
    _mm256_store_ps(out0 + outStride, addSigned<SIGN0>(_mm256_load_ps(out0 + outStride), acc10));
    // Row 1, columns 8 to 15
    _mm256_store_ps(out0 + outStride + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + outStride + 8), acc11));
    // Row 2, columns 0 to 7
    _mm256_store_ps(out0 + 2 * outStride, addSigned<SIGN0>(_mm256_load_ps(out0 + 2 * outStride), acc20));
    // Row 2, columns 8 to 15
    _mm256_store_ps(out0 + 2 * outStride + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + 2 * outStride + 8), acc21));
    // Row 3, columns 0 to 7
    _mm256_store_ps(out0 + 3 * outStride, addSigned<SIGN0>(_mm256_load_ps(out0 + 3 * outStride), acc30));
    // Row 3, columns 8 to 15
    _mm256_store_ps(out0 + 3 * outStride + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + 3 * outStride + 8), acc31));
    // Compile time check, removed when SIGN1 is 0
    if (SIGN1 != 0) {
        // out1 row 0 += SIGN1 * t
        _mm256_store_ps(out1, addSigned<SIGN1>(_mm256_load_ps(out1), acc00));
        // Row 0, columns 8 to 15
        _mm256_store_ps(out1 + 8, addSigned<SIGN1>(_mm256_load_ps(out1 + 8), acc01));
        // Row 1, columns 0 to 7
        _mm256_store_ps(out1 + outStride, addSigned<SIGN1>(_mm256_load_ps(out1 + outStride), acc10));
        // Row 1, columns 8 to 15
        _mm256_store_ps(out1 + outStride + 8, addSigned<SIGN1>(_mm256_load_ps(out1 + outStride + 8), acc11));
        // Row 2, columns 0 to 7
        _mm256_store_ps(out1 + 2 * outStride, addSigned<SIGN1>(_mm256_load_ps(out1 + 2 * outStride), acc20));
        // Row 2, columns 8 to 15
        _mm256_store_ps(out1 + 2 * outStride + 8, addSigned<SIGN1>(_mm256_load_ps(out1 + 2 * outStride + 8), acc21));
        // Row 3, columns 0 to 7
        _mm256_store_ps(out1 + 3 * outStride, addSigned<SIGN1>(_mm256_load_ps(out1 + 3 * outStride), acc30));
        // Row 3, columns 8 to 15
        _mm256_store_ps(out1 + 3 * outStride + 8, addSigned<SIGN1>(_mm256_load_ps(out1 + 3 * outStride + 8), acc31));
    }
}

// Computes rows [rowBegin, rowEnd) of C and leaves the other rows alone.
// Not static because the MPI code calls it with a row range per rank.
// A, B and C are interleaved complex floats
void matrixMultiplyRowRange(int N, const floatType* A, const floatType* B, floatType* C, int rowBegin, int rowEnd, int* args, int argCount)
{
    // Nothing to do
    if (N <= 0) { return; }
    // Clamp the start into range
    if (rowBegin < 0) { rowBegin = 0; }
    // Clamp the end into range
    if (rowEnd > N) { rowEnd = N; }
    // Empty range
    if (rowBegin >= rowEnd) { return; }

    // The reference is column-major, which in row-major terms is C = B * A.
    // B is the left operand
    const float* left = reinterpret_cast<const float*>(B);
    // A is the right operand
    const float* right = reinterpret_cast<const float*>(A);
    // Treated as float pairs, std::complex multiply is slow
    float* out = reinterpret_cast<float*>(C);

    // Optional tuning values, missing or out of range values keep the default.
    // k values per cache block
    int kBlockSize = kBlockSizeDefault;
    // Rows per thread block
    int rowBlockSize = rowBlockSizeDefault;
    // Three or four real multiplies
    int useThreeMult = useThreeMultDefault;
    // Only read args if a pointer was passed
    if (args) {
        // Args[0] sets the k block size
        if (argCount > 0 && args[0] >= 16 && args[0] <= 4096) { kBlockSize = args[0]; }
        // Args[1] sets rows per block, rounded down to a multiple of 4
        if (argCount > 1 && args[1] >= 8 && args[1] <= 8192) { rowBlockSize = args[1] & ~3; }
        // Args[2] switches the three multiply trick on or off
        if (argCount > 2) { useThreeMult = (args[2] == 1) ? 1 : 0; }
    }

    // N rounded up to a multiple of 4
    const int paddedRows = (N + 3) & ~3;
    // N rounded up to a multiple of 16
    const int paddedCols = (N + 15) & ~15;
    // Number of 16 column strips
    const int colBlockCount = paddedCols / 16;
    // The kernels work on four rows at a time, so round the range out to multiples of 4.
    // The extra rows are computed but never copied out.
    // Round the start down to a multiple of 4
    const int alignedBegin = rowBegin & ~3;
    // Round the end up to a multiple of 4
    const int roundedEnd = (rowEnd + 3) & ~3;
    // But not past the padded rows
    const int alignedEnd = (roundedEnd < paddedRows) ? roundedEnd : paddedRows;
    // Allocate the buffers or reuse last call's
    ensureScratch(N, paddedRows, paddedCols, kBlockSize);

    // Start one team of threads for the whole function
    #pragma omp parallel
    {
        // Split the requested rows of the left input into real and imaginary arrays.
        // Static schedule, each thread gets an even share of rows
        // Each thread fills the rows it will read later
        #pragma omp for schedule(static)
        for (int i = rowBegin; i < rowEnd; ++i) {
            // Start of row i, two floats per complex value
            const float* src = left + (size_t)2 * i * N;
            // Row i of the real parts
            float* re = leftRe + (size_t)i * N;
            // Row i of the imaginary parts
            float* im = leftIm + (size_t)i * N;
            // Every column in the row
            for (int k = 0; k < N; ++k) {
                // Even floats are real parts
                re[k] = src[2 * k];
                // Odd floats are imaginary parts
                im[k] = src[2 * k + 1];
            }
            // Sums are only needed for the three multiply version
            if (useThreeMult) {
                // Row i of the sums
                float* sum = leftSum + (size_t)i * N;
                // Every column in the row
                for (int k = 0; k < N; ++k) {
                    // Real plus imaginary
                    sum[k] = re[k] + im[k];
                }
            }
        }

        // All rows of the right input are needed for any row range.
        // Split the right matrix rows across threads
        // Row k of the right matrix
        #pragma omp for schedule(static)
        for (int k = 0; k < N; ++k) {
            // Start of row k, interleaved
            const float* src = right + (size_t)2 * k * N;
            // Padded row k of the real parts
            float* re = rightRe + (size_t)k * paddedCols;
            // Padded row k of the imaginary parts
            float* im = rightIm + (size_t)k * paddedCols;
            // Every column in the row
            for (int j = 0; j < N; ++j) {
                // Real part
                re[j] = src[2 * j];
                // Imaginary part
                im[j] = src[2 * j + 1];
            }
        }

        // The output is accumulated over k blocks, so zero it first.
        // Each thread zeroes its share of the rows
        // Only the rows this call computes
        #pragma omp for schedule(static)
        for (int i = alignedBegin; i < alignedEnd; ++i) {
            // Zero the real row
            memset(outRe + (size_t)i * paddedCols, 0, sizeof(float) * paddedCols);
            // Zero the imaginary row
            memset(outIm + (size_t)i * paddedCols, 0, sizeof(float) * paddedCols);
        }

        // Every thread runs this loop, one k block at a time
        for (int kBegin = 0; kBegin < N; kBegin += kBlockSize) {
            // The last block can be shorter
            const int kCount = (N - kBegin < kBlockSize) ? (N - kBegin) : kBlockSize;

            // Repack this k block of the right matrix so each group of 16 columns is contiguous for every k.
            // Split the strips, the barrier after it means packing is done before use
            // One 16 column strip
            #pragma omp for schedule(static)
            for (int colBlock = 0; colBlock < colBlockCount; ++colBlock) {
                // Where this strip's real parts go
                float* re = packRe + (size_t)colBlock * kCount * 16;
                // Where its imaginary parts go
                float* im = packIm + (size_t)colBlock * kCount * 16;
                // Each k in the block
                for (int k = 0; k < kCount; ++k) {
                    // Copy 16 real floats, 64 bytes
                    memcpy(re + 16 * k, rightRe + (size_t)(kBegin + k) * paddedCols + 16 * colBlock, 64);
                    // Copy 16 imaginary floats
                    memcpy(im + 16 * k, rightIm + (size_t)(kBegin + k) * paddedCols + 16 * colBlock, 64);
                }
                // The three multiply version also needs the sums
                if (useThreeMult) {
                    // Where this strip's sums go
                    float* sum = packSum + (size_t)colBlock * kCount * 16;
                    // Each k in the block
                    for (int k = 0; k < kCount; ++k) {
                        // Sums for columns 0 to 7, real plus imaginary
                        _mm256_store_ps(sum + 16 * k,
                            _mm256_add_ps(_mm256_load_ps(re + 16 * k),
                                          _mm256_load_ps(im + 16 * k)));
                        // Sums for columns 8 to 15
                        _mm256_store_ps(sum + 16 * k + 8,
                            _mm256_add_ps(_mm256_load_ps(re + 16 * k + 8),
                                          _mm256_load_ps(im + 16 * k + 8)));
                    }
                }
            }

            // Each thread takes a row block and a column block. collapse(2) shares the pairs across threads.
            // Turns the two loops into one shared loop
            // One block of rows
            #pragma omp for collapse(2) schedule(static)
            for (int blockBegin = alignedBegin; blockBegin < alignedEnd; blockBegin += rowBlockSize) {
                // One 16 column strip
                for (int colBlock = 0; colBlock < colBlockCount; ++colBlock) {
                    // The last row block can be shorter
                    const int blockEnd = (blockBegin + rowBlockSize < alignedEnd) ? blockBegin + rowBlockSize : alignedEnd;
                    // This strip's packed real parts
                    const float* packedRe = packRe + (size_t)colBlock * kCount * 16;
                    // Packed imaginary parts
                    const float* packedIm = packIm + (size_t)colBlock * kCount * 16;
                    // Packed sums
                    const float* packedSum = packSum + (size_t)colBlock * kCount * 16;
                    // The left rows of the block stay in cache
                    for (int row = blockBegin; row < blockEnd; row += 4) {
                        // Left real parts for these 4 rows, from kBegin
                        const float* lre = leftRe + (size_t)row * N + kBegin;
                        // Left imaginary parts, from kBegin
                        const float* lim = leftIm + (size_t)row * N + kBegin;
                        // Output real tile
                        float* cre = outRe + (size_t)row * paddedCols + 16 * colBlock;
                        // Output imaginary tile
                        float* cim = outIm + (size_t)row * paddedCols + 16 * colBlock;
                        // Four multiply version
                        if (!useThreeMult) {
                            // re = lre*rre - lim*rim, SIGN -1 gives the minus
                            tile4x16TwoProducts<-1>(lre, lim, N, packedRe, packedIm, kCount, cre, paddedCols);
                            // im = lre*rim + lim*rre, SIGN +1 gives the plus
                            tile4x16TwoProducts<+1>(lre, lim, N, packedIm, packedRe, kCount, cim, paddedCols);
                        } else {
                            // Three multiply version
                            // Left sums, from kBegin
                            const float* lsum = leftSum + (size_t)row * N + kBegin;
                            // t1 = lre*rre, re += t1, im -= t1
                            tile4x16OneProduct<+1, -1>(lre, N, packedRe, kCount, cre, cim, paddedCols);
                            // t2 = lim*rim, re -= t2, im -= t2
                            tile4x16OneProduct<-1, -1>(lim, N, packedIm, kCount, cre, cim, paddedCols);
                            // t3 = (lre+lim)*(rre+rim), im += t3, out1 is skipped
                            tile4x16OneProduct<+1, 0>(lsum, N, packedSum, kCount, cim, cim, paddedCols);
                        }
                    }
                }
            }
        }

        // Write the computed rows back as interleaved complex values.
        // Split the rows across threads again
        // Only the requested rows, padding is dropped
        #pragma omp for schedule(static)
        for (int i = rowBegin; i < rowEnd; ++i) {
            // Row i of C, interleaved
            float* dst = out + (size_t)2 * i * N;
            // Computed real row
            const float* re = outRe + (size_t)i * paddedCols;
            // Computed imaginary row
            const float* im = outIm + (size_t)i * paddedCols;
            // Every column
            for (int j = 0; j < N; ++j) {
                // Real part into the even slot
                dst[2 * j] = re[j];
                // Imaginary part into the odd slot
                dst[2 * j + 1] = im[j];
            }
        }
        // End of the parallel region, threads join here
    }
}
/**
* @brief Implements an NxN matrix multiply C=A*B
*	  	 	 	  		   	  			 		 	      
* @param[in] N : dimension of square matrix (NxN)
* @param[in] A : pointer to input NxN matrix
* @param[in] B : pointer to input NxN matrix
* @param[out] C : pointer to output NxN matrix
* @param[in] args : pointer to array of integers which can be used for debugging and performance tweaks. Optional. If unused, set to zero
* @param[in] argCount : the length of the flags array
* @return : your student ID
*	  	 	 	  		   	  			 		 	      
* */
int matrixMultiply(int N, const floatType* A, const floatType* B, floatType* C, int* args, int argCount){	  	 	 	  		   	  			 		 	      
if (N<=0) { return STUDENTID;}//Your code must be able to deal with N=0 scenario without crashing.	  	 	 	  		   	  			 		 	      

// Whole matrix, rows 0 to N
matrixMultiplyRowRange(N, A, B, C, 0, N, args, argCount);

return STUDENTID;	  	 	 	  		   	  			 		 	      

}
	  	 	 	  		   	  			 		 	      
