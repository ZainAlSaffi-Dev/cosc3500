# Study guide for the matrix multiply assignment

This guide explains the three submitted files from the ground up, for someone who has never seen matrix multiply code, AVX, CUDA or MPI before. The aim is that you can explain every line in an interview.

Read Part 1 first. It covers the ideas that every later part relies on. Parts 2, 3 and 4 then walk through `matrixMultiply.cpp`, `matrixMultiplyGPU.cu` and `matrixMultiplyMPI.cpp` in the order the code appears, a few lines at a time. Part 5 is a set of practice interview questions with answers, and Part 6 is a glossary.

The code shown here leaves out the one-line comments that sit above each line in the real files, so the explanations have room. Where a stretch of code repeats the same pattern, part of it is replaced by `...`. Everything else is copied exactly.

---

## Part 1. Background

### 1.1 What a matrix multiply is

A matrix is a grid of numbers. An N by N matrix has N rows and N columns. Rows are numbered with i and columns with j, both starting at 0, so A(i, j) means the number in row i and column j of matrix A.

To multiply two matrices, C = A times B, each entry of C is built like this.

C(i, j) = A(i, 0) × B(0, j) + A(i, 1) × B(1, j) + ... + A(i, N-1) × B(N-1, j)

In words, take row i of A and column j of B, multiply the matching pairs, and add the results. The letter k is used for the position along that row and column, so C(i, j) is the sum over k of A(i, k) × B(k, j).

A 2 by 2 example.

```
A = | 1 2 |     B = | 5 6 |
    | 3 4 |         | 7 8 |

C(0,0) = 1×5 + 2×7 = 19      C(0,1) = 1×6 + 2×8 = 22
C(1,0) = 3×5 + 4×7 = 43      C(1,1) = 3×6 + 4×8 = 50
```

The simplest code is three nested loops.

```cpp
for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
        float sum = 0;
        for (int k = 0; k < N; ++k)
            sum += A[i][k] * B[k][j];
        C[i][j] = sum;
    }
```

There are N × N entries in C and each one needs N multiply-adds, so the total work is N × N × N. The graded size is N = 2048, which is about 8.6 billion multiply-adds. Everything in this assignment is about doing that amount of work quickly.

### 1.2 Complex numbers

This year the matrices hold complex numbers. The file `include/rubric.h` sets `floatType` to `std::complex<float>` on the CPU and `cuComplex` on the GPU.

A complex number has a real part and an imaginary part, written a + bi, where i is a special number with i × i = -1. Each part is a 32-bit float, so one complex number takes 8 bytes.

Adding two complex numbers adds the parts separately. Multiplying them follows from expanding the brackets and using i × i = -1.

(a + bi) × (c + di) = (ac - bd) + (ad + bc)i

So the real part of the product is ac - bd and the imaginary part is ad + bc. That is four real multiplies. For example, (1 + 2i) × (3 + 4i) = (3 - 8) + (4 + 6)i = -5 + 10i.

### 1.3 How a matrix sits in memory

Computer memory is one long line of bytes, so a grid has to be flattened into a line. The usual way in C and C++ is row-major order. Row 0 goes first, then row 1, and so on. Element (i, j) ends up at position i × N + j.

Each complex number is stored as two floats next to each other, real part first. This is called interleaved storage. A 2 by 2 complex matrix looks like this in memory.

```
[ re(0,0), im(0,0), re(0,1), im(0,1), re(1,0), im(1,0), re(1,1), im(1,1) ]
```

So if you treat the buffer as a plain array of floats, the real part of element (i, j) is float number 2 × (i × N + j) and the imaginary part is the float right after it.

### 1.4 The column-major twist, and why B comes first

The grader checks our answer against Intel MKL, and it calls MKL's `cgemm` routine. That routine comes from Fortran, which stores matrices column by column instead of row by row. This is called column-major order.

If you take a matrix stored column-major and read it as if it were row-major, you see its transpose, the same matrix with rows and columns swapped. So when our row-major code looks at the buffers, it sees Aᵀ and Bᵀ, and it has to write Cᵀ.

There is a rule for transposes of products, (A × B)ᵀ = Bᵀ × Aᵀ. The order flips. So the result our code must write is Cᵀ = Bᵀ × Aᵀ, which in the buffers we see means left matrix = the B buffer and right matrix = the A buffer.

That is why all three files multiply B times A. The first version did A times B and failed with an error of about 260 against a tolerance of 0.0000001 (job 568373). An error that large is never rounding. It means a completely different matrix, and swapping the operands fixed it (job 568381).

In the rest of this guide, left means the B buffer and right means the A buffer.

### 1.5 How the marks work

The harness program, GradeBot, runs each of our three functions and the reference on the same random matrices and times both. The score comes from the ratio of our time to the reference time, so lower is better. The references are MKL for the CPU and MPI versions and cuBLAS for the GPU version.

Full marks, 7 out of 7, need a ratio of 1.0 or less on the CPU, 1.5 or less on the GPU and 0.6 or less for MPI across two nodes. The answer also has to be correct. The summed squared difference from the reference over all N × N entries must stay under 1e-7.

The final run scored 21 out of 21. The CPU ratio was 0.81, the GPU ratio 1.40 and the MPI ratio 0.50.

GradeBot calls each function several times in a row. The first call is a warm-up and is not timed. That is why the code keeps its working buffers between calls instead of allocating them every time.

### 1.6 The hardware ideas

**Memory is slow and arithmetic is fast.** A CPU core can do many floating point operations in the time it takes to fetch one value from main memory. So fast code is mostly about getting data to the arithmetic units quickly.

**Caches.** A cache is a small, fast memory close to the core that keeps copies of recently used data. There are several levels. L1 is tiny and fastest, L2 is bigger and a bit slower, and L3 is bigger again. Data moves between memory and cache in 64-byte chunks called cache lines. Reading neighbouring addresses is cheap because they arrive in the same line. Jumping around memory is expensive because every jump pulls in a new line.

**Registers.** Registers are the few storage slots inside the core itself, and they are the fastest storage there is. The best inner loops keep their running totals in registers the whole time.

**SIMD and AVX.** SIMD means single instruction, multiple data. One AVX instruction works on 8 floats at once. The CPU has 16 AVX registers, each 256 bits wide, which is 8 floats. In C++ an AVX register is the type `__m256`, and operations on it are special functions called intrinsics, such as `_mm256_add_ps`. The grader compiles with `-mavx`, which gives the original AVX but not FMA (fused multiply-add). So a multiply and an add are always two separate instructions.

**Threads and OpenMP.** A CPU has several cores, and each can run its own thread of work. The harness gives us 4 threads. OpenMP is a library where you mark a loop with a `#pragma omp` line and it splits the loop's iterations across the threads for you.

**GPUs.** A GPU runs thousands of small threads at once. Threads are grouped into blocks, and the blocks together form a grid. Every 32 threads in a block run in lockstep as a warp. The GPU has large but slow global memory, a small fast shared memory that the threads of one block share, and registers private to each thread.

**MPI.** MPI runs several copies of the program, called ranks, often on different computers called nodes. Ranks do not share memory. They exchange data by sending messages. The MPI test uses two nodes.

### 1.7 The C++ ideas used

**Pointers.** A pointer is a variable that holds a memory address. `float* p` points at a float, and `p[k]` is the float k places after it. Pointer arithmetic counts in elements, not bytes, so `p + 8` is 8 floats (32 bytes) further on. `const float*` means you can read through the pointer but promise not to write.

**static at file level.** A variable declared `static` outside any function exists for the whole run of the program and keeps its value between function calls. It is also private to that file. A `static` function is likewise only visible inside its file.

**const.** A `const` variable cannot change after it is set. `static const int x = 384;` is a named constant.

**size_t.** An unsigned integer type big enough to hold any memory size. The code casts to `size_t` before multiplying indices, because something like 2 × i × N can go past the largest `int` (about 2.1 billion) for big matrices.

**Casts.** `(float*)p` or `reinterpret_cast<float*>(p)` tells the compiler to treat the same memory as a different type. The code uses this to view an array of complex numbers as an array of floats twice as long.

**Templates.** `template<int SIGN>` in front of a function makes the compiler build a separate copy of the function for each value of SIGN that the code uses. Because SIGN is fixed when compiling, a test like `SIGN > 0` is decided by the compiler and disappears from the running program, so the inner loop has no branch.

**inline.** A hint that the function body should be pasted where it is called, so there is no function call overhead.

**The ternary operator.** `cond ? a : b` means a if cond is true, otherwise b.

**Bit tricks for rounding.** `~3` is a number whose bits are all 1 except the lowest two. `x & ~3` clears the lowest two bits of x, which rounds it down to a multiple of 4. `(x + 3) & ~3` rounds up to a multiple of 4. The same idea with 15 rounds to multiples of 16. Also `x >> 2` divides by 4, `x & 3` is the remainder after dividing by 4, and `x ^ 1` flips 0 to 1 and 1 to 0.

**Pragmas.** A `#pragma` line is an instruction to the compiler rather than code. The files use `#pragma omp ...` for OpenMP and `#pragma unroll` on the GPU.

**#include and #define.** `#include <file>` pastes the contents of a header file in at that spot. `#define NAME value` makes the compiler replace NAME with value everywhere.

---

## Part 2. matrixMultiply.cpp, the CPU version

### 2.1 The plan before the lines

The first simple version, three plain loops, ran about 8 times slower than MKL (job 568383). The final version runs faster than MKL on the graded machine. It gets there with five ideas, and every line of the file belongs to one of them.

1. **Split the complex numbers into planes.** Before any arithmetic, each input is copied into two plain float arrays, one holding every real part and one holding every imaginary part. With interleaved storage, 8 neighbouring floats are a mix of real and imaginary parts, which need different arithmetic. In a plane, 8 neighbouring floats are 8 real parts, which is exactly what one AVX instruction wants.
2. **Use three real products instead of four (the 3M trick).** Section 2.6 explains it. It cuts the main arithmetic by a quarter.
3. **A small kernel that keeps a 4 by 16 patch of the answer in registers.** The innermost code computes 4 rows by 16 columns of the output at a time, holding all 64 running totals in 8 AVX registers while it loops over k.
4. **Cache blocking and packing.** The k direction is cut into blocks of 384, and each block of the right matrix is rearranged (packed) so the kernel reads it straight through memory. The rows are cut into blocks of 96 so the data being reused stays in cache.
5. **Four threads with OpenMP.** The work is shared across the 4 threads.

Each call does these steps in order.

1. Copy the left rows it needs into real, imaginary and sum planes.
2. Copy the whole right matrix into real and imaginary planes.
3. Zero the output planes.
4. For each block of 384 k values, pack that block of the right matrix into strips of 16 columns, then run the kernel over every (row block, strip) pair.
5. Copy the output planes back into C as interleaved complex numbers.

### 2.2 The top of the file

```cpp
#include <matrixMultiply.h>
#define STUDENTID 48008361 //DO NOT REMOVE
```

The first line pastes in the provided header. That header brings in `immintrin.h` (the AVX intrinsics), `omp.h` (OpenMP), `string.h` (`memset` and `memcpy`) and `rubric.h`, which defines `floatType` as `std::complex<float>` and `ALIGN` as 64. The second line was provided. GradeBot reads the student ID by calling each function with N = 0.

```cpp
static const int kBlockSizeDefault = 384;
static const int rowBlockSizeDefault = 96;
static const int useThreeMultDefault = 1;
```

These are the settings used on the graded run. `kBlockSizeDefault` is how many k values go in one cache block. `rowBlockSizeDefault` is how many rows go in one row block. `useThreeMultDefault = 1` switches on the 3M trick. The values came from timing runs rather than a formula. A sweep of block sizes (job 568658) found the timing nearly flat across the grid, and the run on the same kind of AMD EPYC chip as the graded nodes (job 568659) gave the best ratio, 0.806, with 384, 96 and 3M on.

### 2.3 The buffers kept between calls

```cpp
static int scratchN = -1;
static int scratchKBlockSize = 0;
static float* leftRe = 0;
static float* leftIm = 0;
static float* leftSum = 0;
static float* rightRe = 0;
static float* rightIm = 0;
static float* outRe = 0;
static float* outIm = 0;
static float* packRe = 0;
static float* packIm = 0;
static float* packSum = 0;
```

These are file-level `static` variables, so they survive from one call to the next. Allocating big buffers is slow, and GradeBot calls the function many times with the same N, so the buffers are made once and reused.

`scratchN` remembers which N the buffers were sized for. It starts at -1, meaning nothing has been allocated yet. `scratchKBlockSize` remembers how big a k block the pack buffers can hold.

Each `float*` starts at 0, which means a null pointer (pointing at nothing).

| Buffer | What it holds |
|---|---|
| `leftRe`, `leftIm` | The real and imaginary planes of the left matrix. Entry `i*N + k` is left(i, k). |
| `leftSum` | leftRe + leftIm for every entry. Only needed by the 3M trick. |
| `rightRe`, `rightIm` | The real and imaginary planes of the right matrix, with each row padded to a multiple of 16 columns. |
| `outRe`, `outIm` | The real and imaginary planes of the answer while it is being built. |
| `packRe`, `packIm`, `packSum` | One k block of the right matrix, rearranged into strips of 16 columns. |

### 2.4 ensureScratch, allocating the buffers

```cpp
static void ensureScratch(int N, int paddedRows, int paddedCols, int kBlockSize)
{
    if (N == scratchN && kBlockSize <= scratchKBlockSize) { return; }
```

If the buffers were already made for this N and are big enough for this block size, there is nothing to do, so the function returns straight away. This is the normal case on every call after the first.

```cpp
    _mm_free(leftRe); _mm_free(leftIm); _mm_free(leftSum);
    _mm_free(rightRe); _mm_free(rightIm);
    _mm_free(outRe); _mm_free(outIm);
    _mm_free(packRe); _mm_free(packIm); _mm_free(packSum);
```

Otherwise the old buffers are freed. `_mm_free` is the partner of `_mm_malloc` below. On the first call every pointer is still null, and freeing a null pointer does nothing, so this is safe.

```cpp
    const int packBlockSize = (kBlockSize > kBlockSizeDefault) ? kBlockSize : kBlockSizeDefault;
```

The pack buffers are sized for whichever is bigger, the requested block size or the default. This means a later call with the default size never has to reallocate.

```cpp
    const size_t leftFloats = (size_t)paddedRows * N;
    const size_t rightFloats = (size_t)N * paddedCols;
    const size_t outFloats = (size_t)paddedRows * paddedCols;
    const size_t packFloats = (size_t)packBlockSize * paddedCols;
```

These work out how many floats each buffer needs. The left planes have padded rows but exactly N columns, because the k direction never needs padding. The right planes have N rows (one per k) and padded columns. The output is padded both ways. A pack buffer holds one k block for every column. The `(size_t)` casts make the multiplication happen in the big integer type so it cannot overflow.

```cpp
    leftRe = (float*)_mm_malloc(sizeof(float) * leftFloats, ALIGN);
    leftIm = (float*)_mm_malloc(sizeof(float) * leftFloats, ALIGN);
    leftSum = (float*)_mm_malloc(sizeof(float) * leftFloats, ALIGN);
    rightRe = (float*)_mm_malloc(sizeof(float) * rightFloats, ALIGN);
    rightIm = (float*)_mm_malloc(sizeof(float) * rightFloats, ALIGN);
    outRe = (float*)_mm_malloc(sizeof(float) * outFloats, ALIGN);
    outIm = (float*)_mm_malloc(sizeof(float) * outFloats, ALIGN);
    packRe = (float*)_mm_malloc(sizeof(float) * packFloats, ALIGN);
    packIm = (float*)_mm_malloc(sizeof(float) * packFloats, ALIGN);
    packSum = (float*)_mm_malloc(sizeof(float) * packFloats, ALIGN);
```

Each line allocates one buffer. `sizeof(float)` is 4 bytes, so the size in bytes is 4 times the float count. `_mm_malloc` works like `malloc` but also takes an alignment, here `ALIGN` = 64, which guarantees the buffer starts at an address that is a multiple of 64. The fast AVX load `_mm256_load_ps` crashes unless its address is a multiple of 32, and 64 also lines the buffers up with cache lines. `_mm_malloc` returns a `void*` (an address with no type), so `(float*)` turns it into a float pointer.

```cpp
    memset(leftRe, 0, sizeof(float) * leftFloats);
    memset(leftIm, 0, sizeof(float) * leftFloats);
    memset(leftSum, 0, sizeof(float) * leftFloats);
    memset(rightRe, 0, sizeof(float) * rightFloats);
    memset(rightIm, 0, sizeof(float) * rightFloats);
```

`memset(p, 0, bytes)` fills memory with zero bytes, and a float made of zero bytes is 0.0. This zeroes the padding. The padding is never written again, so it stays zero forever, and zeros add nothing to a sum. That is what lets the kernel run over padded rows and columns without special edge cases. The output and pack buffers are not zeroed here because they are rewritten on every call.

```cpp
    scratchN = N;
    scratchKBlockSize = packBlockSize;
}
```

Finally the function records what it allocated, so the check at the top works next time.

### 2.5 Two small helpers for signed adds

```cpp
template<int SIGN>
static inline __m256 mulAddSigned(__m256 acc, __m256 x, __m256 y)
{
    return (SIGN > 0) ? _mm256_add_ps(acc, _mm256_mul_ps(x, y))
                      : _mm256_sub_ps(acc, _mm256_mul_ps(x, y));
}
```

This works on whole AVX registers of 8 floats. `_mm256_mul_ps(x, y)` multiplies the 8 pairs of floats lane by lane. `_mm256_add_ps` adds that product to the running total `acc`, and `_mm256_sub_ps` subtracts it instead. So the function returns acc + x × y when SIGN is positive and acc - x × y when it is not.

The sign matters because complex multiplication has a minus in it (ac - bd). Making SIGN a template parameter means the compiler builds one add version and one subtract version, and the `(SIGN > 0) ?` test is decided while compiling. The running code has no branch at all. The multiply and the add are two instructions because the grader's flags do not allow fused multiply-add. The CPU can run a multiply and an add in the same cycle on different units, so this costs less than it looks.

```cpp
template<int SIGN>
static inline __m256 addSigned(__m256 dst, __m256 t)
{
    return (SIGN > 0) ? _mm256_add_ps(dst, t) : _mm256_sub_ps(dst, t);
}
```

The same idea without the multiply. It returns dst + t or dst - t. It is used when a finished product is added into, or taken away from, the output.

### 2.6 The 3M trick

The complex product (a + bi) × (c + di) needs ac, bd, ad and bc, which is four real products. There is an old rearrangement, usually credited to Gauss, that needs only three.

```
t1 = a × c
t2 = b × d
t3 = (a + b) × (c + d)

real part      = t1 - t2
imaginary part = t3 - t1 - t2
```

Check it by expanding t3 = ac + ad + bc + bd. Taking away t1 = ac and t2 = bd leaves ad + bc, which is the imaginary part.

For whole matrices the same thing works with matrix products. With the left planes Lre and Lim and the right planes Rre and Rim, the code computes three big real products, T1 = Lre × Rre, T2 = Lim × Rim and T3 = (Lre + Lim) × (Rre + Rim). The sums Lre + Lim and Rre + Rim are computed once per entry, which is cheap next to the products. Three big products instead of four is 25 percent less of the dominant work.

The worry with this trick is accuracy, because the subtractions can cancel digits. The measured error went down, from 2.3e-08 with four products to about 5e-09 with three. The likely reason is that MKL itself uses a similar method, so our rounding matches the reference more closely. Either way it is far under the 1e-7 tolerance.

### 2.7 tile4x16OneProduct, the kernel used on the graded run

This function computes one of the three products, T1, T2 or T3, for a 4 row by 16 column patch of the output, over one k block. Then it adds the patch into one or two output planes with a chosen sign.

```cpp
template<int SIGN0, int SIGN1>
static void tile4x16OneProduct(const float* left, long leftStride,
                               const float* pack, int kCount,
                               float* out0, float* out1, long outStride)
{
```

`SIGN0` and `SIGN1` are the signs used when adding the finished patch into `out0` and `out1`. A `SIGN1` of 0 means skip the second output.

`left` points at row 0 of the four left rows, at the start of this k block. The next row is `leftStride` floats further on (N floats, since the left planes are N wide). `pack` points at the packed strip of 16 columns for this k block, and `kCount` is how many k values the block holds. `out0` and `out1` point at the top-left corner of the 4 by 16 output patch, and each output row is `outStride` floats after the last.

```cpp
    __m256 acc00 = _mm256_setzero_ps();
    __m256 acc01 = _mm256_setzero_ps();
    __m256 acc10 = _mm256_setzero_ps();
    __m256 acc11 = _mm256_setzero_ps();
    __m256 acc20 = _mm256_setzero_ps();
    __m256 acc21 = _mm256_setzero_ps();
    __m256 acc30 = _mm256_setzero_ps();
    __m256 acc31 = _mm256_setzero_ps();
```

Eight AVX registers, each holding 8 running totals, all set to zero. Together they hold the 4 × 16 = 64 totals of the patch. In the name `accRH`, R is the row (0 to 3) and H is which half of the 16 columns (0 means columns 0 to 7, 1 means columns 8 to 15). So `acc21` is row 2, columns 8 to 15. They start at zero because this product is built on its own before being added in with a sign.

```cpp
    for (int k = 0; k < kCount; ++k) {
        const __m256 p0 = _mm256_load_ps(pack + 16 * k);
        const __m256 p1 = _mm256_load_ps(pack + 16 * k + 8);
```

The loop walks through the k block. In the packed strip, the 16 values for step k sit together at `pack + 16*k`. `_mm256_load_ps` loads 8 floats into a register, so `p0` gets the right matrix values for columns 0 to 7 at this k and `p1` gets columns 8 to 15.

```cpp
        __m256 l = _mm256_broadcast_ss(left + k);
        acc00 = mulAddSigned<1>(acc00, l, p0);
        acc01 = mulAddSigned<1>(acc01, l, p1);
```

This is the heart of the whole CPU version. Output entry (row 0, column j) needs left(row 0, k) × right(k, j) added in for every k. For a fixed k, the left value is the same for all 16 columns. So `_mm256_broadcast_ss` copies that single float into all 8 lanes of a register `l`, and one multiply-add with `p0` updates 8 output totals at once, then one with `p1` updates the other 8. Two multiply and add pairs do 16 multiply-adds.

```cpp
        l = _mm256_broadcast_ss(left + leftStride + k);
        acc10 = mulAddSigned<1>(acc10, l, p0);
        acc11 = mulAddSigned<1>(acc11, l, p1);
        l = _mm256_broadcast_ss(left + 2 * leftStride + k);
        acc20 = mulAddSigned<1>(acc20, l, p0);
        acc21 = mulAddSigned<1>(acc21, l, p1);
        l = _mm256_broadcast_ss(left + 3 * leftStride + k);
        acc30 = mulAddSigned<1>(acc30, l, p0);
        acc31 = mulAddSigned<1>(acc31, l, p1);
    }
```

The same thing for rows 1, 2 and 3. Their left values are one, two and three `leftStride`s further on. Notice that `p0` and `p1` are loaded once and used by all four rows. That reuse is why the patch is 4 rows tall.

The patch shape comes from the register count. There are 16 AVX registers. The 8 accumulators, `p0`, `p1` and `l` use 11 of them. A bigger patch would not fit, and a smaller one would not have enough independent additions in flight to keep the add unit busy.

```cpp
    _mm256_store_ps(out0, addSigned<SIGN0>(_mm256_load_ps(out0), acc00));
    _mm256_store_ps(out0 + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + 8), acc01));
    _mm256_store_ps(out0 + outStride, addSigned<SIGN0>(_mm256_load_ps(out0 + outStride), acc10));
    _mm256_store_ps(out0 + outStride + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + outStride + 8), acc11));
    _mm256_store_ps(out0 + 2 * outStride, addSigned<SIGN0>(_mm256_load_ps(out0 + 2 * outStride), acc20));
    _mm256_store_ps(out0 + 2 * outStride + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + 2 * outStride + 8), acc21));
    _mm256_store_ps(out0 + 3 * outStride, addSigned<SIGN0>(_mm256_load_ps(out0 + 3 * outStride), acc30));
    _mm256_store_ps(out0 + 3 * outStride + 8, addSigned<SIGN0>(_mm256_load_ps(out0 + 3 * outStride + 8), acc31));
```

After the k loop the patch is finished. Each line reads 8 floats of the output plane, adds or subtracts the matching accumulator, and writes the 8 floats back with `_mm256_store_ps`. The eight lines cover rows 0 to 3, each in two halves.

```cpp
    if (SIGN1 != 0) {
        _mm256_store_ps(out1, addSigned<SIGN1>(_mm256_load_ps(out1), acc00));
        ...
    }
}
```

The same again for the second output, when there is one. SIGN1 is a template parameter, so when it is 0 the compiler deletes this whole block.

### 2.8 tile4x16TwoProducts, the four-product kernel

This kernel is only used when the optional args switch the 3M trick off. It was the main kernel before 3M was added and it stays so the two can be compared.

```cpp
template<int SIGN>
static void tile4x16TwoProducts(const float* leftA, const float* leftB, long leftStride,
                                const float* packA, const float* packB, int kCount,
                                float* out, long outStride)
```

It computes out += leftA × packA + SIGN × leftB × packB for a 4 by 16 patch. With four products, the real part is Lre × Rre - Lim × Rim and the imaginary part is Lre × Rim + Lim × Rre. Each of those is two products joined by a sign, so one call of this kernel builds one whole output plane.

The body is the same as 2.7 with three differences. The accumulators start by loading the current output (`_mm256_load_ps(out)`) instead of zero, because there is only one output and no separate sign to apply at the end. Each k step loads 16 values from both `packA` (as `a0`, `a1`) and `packB` (as `b0`, `b1`). And each row does two broadcasts, one from `leftA` multiplied into `a0`/`a1` and one from `leftB` multiplied into `b0`/`b1` with the sign SIGN. At the end the 8 accumulators are stored straight back to `out`.

### 2.9 matrixMultiplyRowRange, the main function

```cpp
void matrixMultiplyRowRange(int N, const floatType* A, const floatType* B, floatType* C, int rowBegin, int rowEnd, int* args, int argCount)
{
```

This function computes rows `rowBegin` up to but not including `rowEnd` of the answer, and leaves every other row of C alone. The CPU version calls it with the full range 0 to N. The MPI version calls it with one part of the rows on each node. It is not `static`, because the MPI file needs to call it.

```cpp
    if (N <= 0) { return; }
    if (rowBegin < 0) { rowBegin = 0; }
    if (rowEnd > N) { rowEnd = N; }
    if (rowBegin >= rowEnd) { return; }
```

Safety checks. An empty matrix means nothing to do. The row range is clamped into 0 to N, and an empty range returns.

```cpp
    const float* left = reinterpret_cast<const float*>(B);
    const float* right = reinterpret_cast<const float*>(A);
    float* out = reinterpret_cast<float*>(C);
```

Two things happen here. First, B becomes the left matrix and A the right matrix, which is the column-major swap from 1.4. Second, the complex arrays are viewed as plain float arrays, real and imaginary interleaved. Working on raw floats avoids the C++ complex multiply, which under the grader's flags calls a slow library function named `__mulsc3` on every multiply.

```cpp
    int kBlockSize = kBlockSizeDefault;
    int rowBlockSize = rowBlockSizeDefault;
    int useThreeMult = useThreeMultDefault;
    if (args) {
        if (argCount > 0 && args[0] >= 16 && args[0] <= 4096) { kBlockSize = args[0]; }
        if (argCount > 1 && args[1] >= 8 && args[1] <= 8192) { rowBlockSize = args[1] & ~3; }
        if (argCount > 2) { useThreeMult = (args[2] == 1) ? 1 : 0; }
    }
```

Start with the defaults. GradeBot can pass extra integers on its command line, and they arrive here as the `args` array with `argCount` entries. They were used for tuning runs so one job could try many settings without recompiling. Every read checks that the pointer is not null and the entry exists, and every value is range checked. `args[1] & ~3` rounds the row block size down to a multiple of 4, because the kernel always handles 4 rows at a time. On the graded run no args are passed, so the defaults are used.

```cpp
    const int paddedRows = (N + 3) & ~3;
    const int paddedCols = (N + 15) & ~15;
    const int colBlockCount = paddedCols / 16;
```

The rows are rounded up to a multiple of 4 and the columns to a multiple of 16, matching the 4 by 16 kernel. `colBlockCount` is the number of 16-column strips. For N = 2048 nothing changes and there are 128 strips. For N = 100 it gives 100 rows, 112 columns and 7 strips.

```cpp
    const int alignedBegin = rowBegin & ~3;
    const int roundedEnd = (rowEnd + 3) & ~3;
    const int alignedEnd = (roundedEnd < paddedRows) ? roundedEnd : paddedRows;
```

The kernel works on groups of 4 rows, so the requested range is widened to multiples of 4. For example rows 37 to 70 become 36 to 72. The end is capped at the padded row count. The extra rows are computed but never copied into C, so they do no harm. This only matters for MPI, where a node's range can start at any row.

```cpp
    ensureScratch(N, paddedRows, paddedCols, kBlockSize);
```

Make the buffers, or reuse the ones from last time.

#### The parallel region

```cpp
    #pragma omp parallel
    {
```

This starts the 4 OpenMP threads, and all of them run everything inside the braces. Starting the threads once for the whole function, instead of once per loop, saves the start-up cost each time.

Inside the region, `#pragma omp for` in front of a loop splits that loop's iterations among the threads. `schedule(static)` means the iterations are cut into equal chunks up front, which suits loops where every iteration costs the same. At the end of every `omp for` loop there is a hidden barrier. No thread goes past it until every thread has finished the loop. The steps below depend on those barriers to happen in order.

#### Step 1. Split the left rows into planes

```cpp
        #pragma omp for schedule(static)
        for (int i = rowBegin; i < rowEnd; ++i) {
            const float* src = left + (size_t)2 * i * N;
            float* re = leftRe + (size_t)i * N;
            float* im = leftIm + (size_t)i * N;
```

Each thread takes an equal share of the requested rows. `src` points at the start of row i in the input. Each complex number is two floats, so row i starts at float 2 × i × N. `re` and `im` point at row i of the two planes.

```cpp
            for (int k = 0; k < N; ++k) {
                re[k] = src[2 * k];
                im[k] = src[2 * k + 1];
            }
```

Walk along the row. The even floats are real parts and go into `re`. The odd floats are imaginary parts and go into `im`.

```cpp
            if (useThreeMult) {
                float* sum = leftSum + (size_t)i * N;
                for (int k = 0; k < N; ++k) {
                    sum[k] = re[k] + im[k];
                }
            }
        }
```

With the 3M trick on, also store real + imaginary for each entry. This is the Lre + Lim needed for T3.

#### Step 2. Split the right matrix into planes

```cpp
        #pragma omp for schedule(static)
        for (int k = 0; k < N; ++k) {
            const float* src = right + (size_t)2 * k * N;
            float* re = rightRe + (size_t)k * paddedCols;
            float* im = rightIm + (size_t)k * paddedCols;
            for (int j = 0; j < N; ++j) {
                re[j] = src[2 * j];
                im[j] = src[2 * j + 1];
            }
        }
```

The same split for the right matrix, but for all N rows, because every output row needs the whole right matrix. The only difference is that each plane row is `paddedCols` wide instead of N, and the padding columns stay at zero from `ensureScratch`.

#### Step 3. Zero the output

```cpp
        #pragma omp for schedule(static)
        for (int i = alignedBegin; i < alignedEnd; ++i) {
            memset(outRe + (size_t)i * paddedCols, 0, sizeof(float) * paddedCols);
            memset(outIm + (size_t)i * paddedCols, 0, sizeof(float) * paddedCols);
        }
```

The kernels add into the output planes a block at a time, so the rows being computed must start at zero. Each thread zeroes its share of the rows.

#### Step 4. Loop over the k blocks

```cpp
        for (int kBegin = 0; kBegin < N; kBegin += kBlockSize) {
            const int kCount = (N - kBegin < kBlockSize) ? (N - kBegin) : kBlockSize;
```

There is no `omp for` on this loop, so every thread runs it itself, all four in step. The k direction is processed 384 values at a time. `kCount` is the size of the current block, which is 384 except possibly for the last block, which is whatever is left.

Why block at all? A 2048 by 2048 plane is 16 MB, far bigger than any cache. Working on one k block at a time keeps the data the kernel is reusing small enough to stay in cache.

#### Step 4a. Pack this k block of the right matrix

```cpp
            #pragma omp for schedule(static)
            for (int colBlock = 0; colBlock < colBlockCount; ++colBlock) {
                float* re = packRe + (size_t)colBlock * kCount * 16;
                float* im = packIm + (size_t)colBlock * kCount * 16;
```

The strips are shared among the threads. Each strip gets `kCount × 16` floats in the pack buffers, one after another, and `re` and `im` point at where this strip's data goes.

```cpp
                for (int k = 0; k < kCount; ++k) {
                    memcpy(re + 16 * k, rightRe + (size_t)(kBegin + k) * paddedCols + 16 * colBlock, 64);
                    memcpy(im + 16 * k, rightIm + (size_t)(kBegin + k) * paddedCols + 16 * colBlock, 64);
                }
```

For each k in the block, copy the strip's 16 floats (64 bytes) from row kBegin + k of the right planes. `memcpy(dest, src, bytes)` copies bytes from one place to another.

This is what packing buys. In the right plane, the kernel would move from one k to the next by jumping a whole row, `paddedCols` floats, which is 8 KB at N = 2048. After packing, the 16 values for k + 1 sit right after the 16 values for k, so the kernel reads one straight line through memory, one 64-byte cache line per step.

```cpp
                if (useThreeMult) {
                    float* sum = packSum + (size_t)colBlock * kCount * 16;
                    for (int k = 0; k < kCount; ++k) {
                        _mm256_store_ps(sum + 16 * k,
                            _mm256_add_ps(_mm256_load_ps(re + 16 * k),
                                          _mm256_load_ps(im + 16 * k)));
                        _mm256_store_ps(sum + 16 * k + 8,
                            _mm256_add_ps(_mm256_load_ps(re + 16 * k + 8),
                                          _mm256_load_ps(im + 16 * k + 8)));
                    }
                }
            }
```

With 3M on, the packed sums Rre + Rim are also needed for T3. For each k, load 8 real and 8 imaginary values, add them lane by lane and store the result, then the same for the other 8 columns.

The hidden barrier at the end of this `omp for` matters. No thread starts the next step until every strip has been packed.

#### Step 4b. Run the kernels

```cpp
            #pragma omp for collapse(2) schedule(static)
            for (int blockBegin = alignedBegin; blockBegin < alignedEnd; blockBegin += rowBlockSize) {
                for (int colBlock = 0; colBlock < colBlockCount; ++colBlock) {
```

The work is every pair of (row block, column strip). `collapse(2)` tells OpenMP to treat the two loops as one long list of pairs and share that list among the threads. With N = 2048 there are 22 row blocks of 96 rows and 128 strips, so 2816 pairs for 4 threads, which spreads evenly. Without `collapse`, only the 22 row blocks would be shared, and 22 does not divide evenly by 4.

```cpp
                    const int blockEnd = (blockBegin + rowBlockSize < alignedEnd) ? blockBegin + rowBlockSize : alignedEnd;
                    const float* packedRe = packRe + (size_t)colBlock * kCount * 16;
                    const float* packedIm = packIm + (size_t)colBlock * kCount * 16;
                    const float* packedSum = packSum + (size_t)colBlock * kCount * 16;
```

`blockEnd` is where this row block stops, shorter for the last block. The three pointers find this strip's packed real, imaginary and sum data.

```cpp
                    for (int row = blockBegin; row < blockEnd; row += 4) {
                        const float* lre = leftRe + (size_t)row * N + kBegin;
                        const float* lim = leftIm + (size_t)row * N + kBegin;
                        float* cre = outRe + (size_t)row * paddedCols + 16 * colBlock;
                        float* cim = outIm + (size_t)row * paddedCols + 16 * colBlock;
```

Walk down the row block 4 rows at a time. `lre` and `lim` point at this group's left rows, starting at column kBegin, the start of the k block. `cre` and `cim` point at the top-left corner of this group's 4 by 16 patch in the output planes.

The same packed strip (384 × 16 floats, 24 KB per plane) is used by all 24 groups of 4 rows in the block, so it stays in cache while they use it. And because `collapse` hands each thread neighbouring pairs, a thread mostly works through strips of the same row block one after another, so those left rows stay in cache too.

```cpp
                        if (!useThreeMult) {
                            tile4x16TwoProducts<-1>(lre, lim, N, packedRe, packedIm, kCount, cre, paddedCols);
                            tile4x16TwoProducts<+1>(lre, lim, N, packedIm, packedRe, kCount, cim, paddedCols);
                        } else {
                            const float* lsum = leftSum + (size_t)row * N + kBegin;
                            tile4x16OneProduct<+1, -1>(lre, N, packedRe, kCount, cre, cim, paddedCols);
                            tile4x16OneProduct<-1, -1>(lim, N, packedIm, kCount, cre, cim, paddedCols);
                            tile4x16OneProduct<+1, 0>(lsum, N, packedSum, kCount, cim, cim, paddedCols);
                        }
```

The `if` branch is the four-product route, used only if args turned 3M off. The first call builds the real part, Lre × Rre - Lim × Rim, so SIGN is -1. The second builds the imaginary part, Lre × Rim + Lim × Rre, so SIGN is +1 and the packed planes are passed the other way round.

The `else` branch is the 3M route, used on the graded run. `lsum` points at the left sums for these rows.

1. T1 = Lre × Rre. It is added to the real part and taken away from the imaginary part, so the signs are +1 and -1.
2. T2 = Lim × Rim. It is taken away from both, so -1 and -1.
3. T3 = (Lre + Lim) × (Rre + Rim). It is added to the imaginary part only, so +1 and 0. The 0 means the second output is skipped, which is why `cim` is simply passed twice.

After all three, real = T1 - T2 and imaginary = T3 - T1 - T2, exactly the formulas from 2.6.

#### Step 5. Merge back into C

```cpp
        #pragma omp for schedule(static)
        for (int i = rowBegin; i < rowEnd; ++i) {
            float* dst = out + (size_t)2 * i * N;
            const float* re = outRe + (size_t)i * paddedCols;
            const float* im = outIm + (size_t)i * paddedCols;
            for (int j = 0; j < N; ++j) {
                dst[2 * j] = re[j];
                dst[2 * j + 1] = im[j];
            }
        }
    }
}
```

This reverses step 1. For each requested row, the real and imaginary output planes are woven back together, real parts into the even slots of C and imaginary parts into the odd slots. Only the requested rows and the first N columns are written, so the padding and the throwaway rows from the widened range never reach C. The closing brace of the parallel region is where the threads finish.

### 2.10 matrixMultiply, the function GradeBot calls

```cpp
int matrixMultiply(int N, const floatType* A, const floatType* B, floatType* C, int* args, int argCount){
if (N<=0) { return STUDENTID;}//Your code must be able to deal with N=0 scenario without crashing.

matrixMultiplyRowRange(N, A, B, C, 0, N, args, argCount);

return STUDENTID;
}
```

The signature, the N <= 0 guard and the `return STUDENTID;` line were provided and are left exactly as given. GradeBot calls with N = 0 to read the student ID, so the guard has to come first. The one line in the middle computes every row, 0 to N.

---

## Part 3. matrixMultiplyGPU.cu, the GPU version

### 3.1 GPU basics you need for this file

A `.cu` file mixes two kinds of code. Host code runs on the CPU and organises the work. Kernels run on the GPU. A kernel is marked `__global__`, and one call of it starts thousands of threads that all run the same kernel code at once. A helper function that only kernels call is marked `__device__`, and a host function can be marked `__host__`.

Threads are arranged in blocks, and blocks are arranged in a grid. Both can be 2D. Inside a kernel each thread can read where it is.

| Name | Meaning |
|---|---|
| `threadIdx.x`, `threadIdx.y` | This thread's position inside its block |
| `blockIdx.x`, `blockIdx.y` | This block's position inside the grid |
| `blockDim.x`, `blockDim.y` | How many threads a block has in each direction |

So `blockIdx.x * blockDim.x + threadIdx.x` gives a thread's overall x position across the whole grid. The code uses x for columns and y for rows.

The hardware runs threads in groups of 32 called warps. The 32 threads of a warp have consecutive `threadIdx.x` values and run each instruction together.

There are three kinds of memory that matter here.

1. **Global memory** is the GPU's main memory. It is large but each access takes hundreds of cycles. When the 32 threads of a warp read or write neighbouring addresses, the hardware merges them into a few wide transactions. This is called coalescing, and it is the first rule of GPU speed.
2. **Shared memory** is a small, fast scratchpad (48 KB here) shared by the threads of one block. It is split into 32 banks. Threads of a warp that hit different banks, or the exact same address, are served at once. Threads that hit different addresses in the same bank have to wait in turn, which is called a bank conflict.
3. **Registers** are private to each thread and fastest of all.

`__syncthreads()` is a barrier for one block. No thread in the block continues past it until every thread in the block has reached it. It is needed whenever threads read shared memory that other threads wrote.

On the GPU a complex number is a `cuComplex`, which is the same as `float2`, a pair of floats named `.x` (real part) and `.y` (imaginary part). `make_cuComplex(re, im)` builds one. A `float4` is four floats, `.x .y .z .w`, which is 16 bytes, so it can hold two complex numbers.

The grader compiles with `--gpu-architecture=sm_35`, a target for a very old GPU generation. The graded A100 GPUs run the program by compiling an intermediate code called PTX into their own machine code when the program starts. One result is that the newer hardware features, like the tensor cores, cannot be used, so everything runs on ordinary floating point units.

### 3.2 The top of the file and the device buffers

```cpp
#include <matrixMultiplyGPU.cuh>
#define STUDENTID 48008361 //DO NOT REMOVE

static const int kernelChoiceDefault = 0;
```

The provided header brings in `rubric.h` (which makes `floatTypeCUDA` mean `cuComplex`) and the declaration of the kernel we must define. The file has three kernels, and `kernelChoiceDefault` picks which one runs when no flags are passed. 0 is the register tiled kernel, the fastest.

```cpp
static int deviceN = -1;
static floatTypeCUDA* dLeft = 0;
static floatTypeCUDA* dRight = 0;
static floatTypeCUDA* dOut = 0;
```

Pointers to three buffers in GPU memory, one each for the left input, right input and output. The `d` at the front is a common way of marking device (GPU) memory. They are `static` for the same reason as on the CPU, so they are allocated once and reused across calls. `deviceN` remembers the size they were made for.

```cpp
static void ensureDevice(int N)
{
    if (N == deviceN) { return; }
    cudaFree(dLeft);
    cudaFree(dRight);
    cudaFree(dOut);
    const size_t bytes = sizeof(floatTypeCUDA) * (size_t)N * N;
    cudaMalloc((void**)&dLeft, bytes);
    cudaMalloc((void**)&dRight, bytes);
    cudaMalloc((void**)&dOut, bytes);
    deviceN = N;
}
```

The same pattern as `ensureScratch`. If the size is unchanged, return. Otherwise free the old GPU buffers (freeing a null pointer is allowed) and allocate new ones. One N by N complex matrix is 8 × N × N bytes. `cudaMalloc` allocates GPU memory and has to write the new address into our pointer variable, so we pass the address of the pointer, `&dLeft`, cast to the type `void**` that `cudaMalloc` expects.

### 3.3 The shared memory tiled kernel (choice 2)

This kernel is easier to follow than the default one, and the default one builds on its ideas, so it comes first.

The naive approach (3.6) has each thread read a whole row of the left matrix and a whole column of the right matrix from global memory. Neighbouring threads read many of the same values, so the same data is fetched from slow memory again and again. This kernel has each block of 32 × 32 threads load a 32 by 32 tile of each input into shared memory once, and then all 1024 threads read from the fast copy.

```cpp
__global__ void matrixMultiplyKernelTiled_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C)
{
    __shared__ floatTypeCUDA As[32][33];
    __shared__ floatTypeCUDA Bs[32][33];
```

Two tiles in shared memory, one for each input. They are 33 wide instead of 32. The extra column shifts each row onto different banks, which is a standard, cheap guard against bank conflicts when threads walk down a column. Inside the kernel, A and B are simply the kernel's left and right inputs. The host passes the B buffer as A here, because of the swap from 1.4.

```cpp
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int j = blockIdx.x * 32 + tx;
    const int i = blockIdx.y * 32 + ty;
    float re = 0.0f;
    float im = 0.0f;
```

Each thread computes one output entry, C(i, j). `tx` and `ty` are its column and row inside the block, and `j` and `i` are the column and row in the whole matrix. `re` and `im` are its running totals, kept in registers.

```cpp
    for (int t = 0; t < N; t += 32) {
        As[ty][tx] = (i < N && t + tx < N) ? A[(size_t)i * N + t + tx]
                                           : make_cuComplex(0.0f, 0.0f);
        Bs[ty][tx] = (t + ty < N && j < N) ? B[(size_t)(t + ty) * N + j]
                                           : make_cuComplex(0.0f, 0.0f);
        __syncthreads();
```

The loop steps along k 32 at a time. Each of the 1024 threads loads one element of each tile. Thread (tx, ty) loads A(i, t + tx) and B(t + ty, j). Anything past the edge of the matrix is loaded as zero, so the tiles always work at full size. The `__syncthreads()` waits until every thread has finished loading, because each thread is about to read values that other threads loaded.

```cpp
        for (int k = 0; k < 32; ++k) {
            const floatTypeCUDA l = As[ty][k];
            const floatTypeCUDA r = Bs[k][tx];
            re += l.x * r.x - l.y * r.y;
            im += l.x * r.y + l.y * r.x;
        }
        __syncthreads();
    }
```

Now each thread does 32 multiply-adds from shared memory, along row ty of the A tile and column tx of the B tile. The two lines inside are the complex multiply from 1.2 written out by hand, (ac - bd) and (ad + bc). The second `__syncthreads()` makes sure every thread has finished reading the tiles before the next loop iteration overwrites them.

```cpp
    if (i < N && j < N) {
        C[(size_t)i * N + j] = make_cuComplex(re, im);
    }
}
```

Threads that sit past the edge of the matrix still helped load tiles, but they do not write anything. Every other thread writes its answer.

Each value loaded from global memory is now used by 32 threads instead of 1. This kernel scored a ratio of 1.58 (job 568767).

### 3.4 loadOrZero, a small helper

```cpp
static __device__ __forceinline__ float2 loadOrZero(const float2* P, int N, int r, int c)
{
    if (r < N && c < N) { return __ldg(P + (size_t)r * N + c); }
    return make_cuComplex(0.0f, 0.0f);
}
```

A GPU helper used by the default kernel. It returns element (r, c) of a matrix, or zero if (r, c) is outside the matrix. `__forceinline__` makes the compiler paste it into the kernel. `__ldg` loads through the GPU's read-only data cache, which is allowed because the inputs never change during the kernel.

### 3.5 The register tiled kernel (choice 0, the default)

#### Why it is faster

In the tiled kernel, each thread does one complex multiply-add for every two values it reads from shared memory. Shared memory is fast but not free, and at that ratio it is shared memory, not arithmetic, that sets the speed.

The fix is for each thread to compute more outputs. Here each thread computes a 4 by 4 patch of 16 outputs, held in registers. For each k, it reads 4 values of the left tile and 4 of the right tile, and does 16 multiply-adds with them. That is 16 multiply-adds per 8 reads instead of 1 per 2, four times better. This took the ratio from 1.58 to 1.32 (job 568780), and the final graded run scored 1.40.

#### The layout

```
One block = 16 × 16 = 256 threads, computing a 64 × 64 tile of C.

Thread (tx, ty) owns these 16 outputs inside the tile
    rows    ty, ty+16, ty+32, ty+48
    columns tx, tx+16, tx+32, tx+48

k advances 16 at a time. For each step the block needs
    a 64 × 16 piece of the left matrix  (64 rows, 16 k values)
    a 16 × 64 piece of the right matrix (16 k values, 64 columns)
```

A thread's rows and columns are spaced 16 apart instead of being next to each other. That way, at every moment, the 16 threads in a row of the block touch 16 neighbouring columns, which keeps shared memory reads free of bank conflicts and the final global writes coalesced.

#### The code

```cpp
__global__ void matrixMultiplyKernelReg_GPU(int N, const floatTypeCUDA* __restrict__ A, const floatTypeCUDA* __restrict__ B, floatTypeCUDA* __restrict__ C)
{
```

`__restrict__` is a promise to the compiler that the three pointers never point at overlapping memory. That lets it keep loaded values around and reorder loads freely.

```cpp
    __shared__ float2 As[2][16][66];
    __shared__ float2 Bs[2][16][66];
```

The shared tiles. The first index, 2, means there are two copies of each tile (double buffering, explained below). The second index is k, 16 values. The third is the row (for As) or column (for Bs), 64 of them plus 2 of padding. The A tile is stored with k first, the transpose of how it sits in memory, so that for one k the 64 row values sit next to each other, just like the 64 column values of the B tile. The padding of 2 makes each row 66 × 8 = 528 bytes, which is a multiple of 16, so 16-byte `float4` stores into the tile stay aligned. Together the tiles use about 33 KB of the 48 KB of shared memory.

```cpp
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int tid = ty * 16 + tx;
    const int rowBase = blockIdx.y * 64;
    const int colBase = blockIdx.x * 64;
```

`tx` and `ty` both run 0 to 15. `tid` numbers the 256 threads 0 to 255 in one flat list, which is handy for dividing up the loading work. `rowBase` and `colBase` are the first row and column of this block's 64 by 64 tile.

```cpp
    const int aRow = tid >> 2;
    const int aK = (tid & 3) * 2;
    const int bK = tid >> 4;
    const int bCol = (tid & 15) * 2;
```

Each step needs 64 × 16 = 1024 left values and 16 × 64 = 1024 right values, which is 4 of each per thread. These lines give each thread its share.

For the left tile, `aRow = tid / 4` gives rows 0 to 63, four threads per row. `aK = (tid % 4) × 2` gives 0, 2, 4 or 6. Each thread loads k = aK and aK + 1 (one `float4`) and k = aK + 8 and aK + 9 (another). The four threads of a row between them cover all 16 k values.

For the right tile, `bK = tid / 16` gives k rows 0 to 15, sixteen threads per row. `bCol = (tid % 16) × 2` gives even columns 0 to 30. Each thread loads columns bCol and bCol + 1, then bCol + 32 and bCol + 33. The sixteen threads of a row cover all 64 columns, and in each load neighbouring threads read neighbouring addresses, so the loads are coalesced.

```cpp
    const float2* pA = (const float2*)A;
    const float2* pB = (const float2*)B;
```

The inputs viewed as `float2`, a pair of floats. For `cuComplex` this is the same memory under another name.

```cpp
    const bool vectorLoadA = (rowBase + 64 <= N) && ((N & 1) == 0);
    const bool vectorLoadB = (colBase + 64 <= N) && ((N & 1) == 0);
```

A `float4` load reads 16 bytes, two complex numbers, in one instruction, but only if the address is a multiple of 16. That holds when N is even (`N & 1` is the lowest bit, 0 for even numbers), because then every row starts on a 16-byte boundary. The tile must also sit fully inside the matrix so no element is read past the edge. When either fails, the kernel uses the slower one-at-a-time path with `loadOrZero`. For N = 2048 the fast path is always used.

```cpp
    float2 acc[4][4];
    #pragma unroll
    for (int m = 0; m < 4; ++m) {
        #pragma unroll
        for (int n = 0; n < 4; ++n) { acc[m][n] = make_cuComplex(0.0f, 0.0f); }
    }
```

The 16 running totals, all zeroed. `#pragma unroll` tells the compiler to write out the loop body once per iteration instead of looping. After unrolling, every `acc[m][n]` has a fixed index known at compile time, and that is what allows the array to live in registers. With a changing index it would be pushed out to slow memory.

```cpp
    As[0][aK][aRow] = loadOrZero(pA, N, rowBase + aRow, aK);
    As[0][aK + 1][aRow] = loadOrZero(pA, N, rowBase + aRow, aK + 1);
    As[0][aK + 8][aRow] = loadOrZero(pA, N, rowBase + aRow, aK + 8);
    As[0][aK + 9][aRow] = loadOrZero(pA, N, rowBase + aRow, aK + 9);
    Bs[0][bK][bCol] = loadOrZero(pB, N, bK, colBase + bCol);
    Bs[0][bK][bCol + 1] = loadOrZero(pB, N, bK, colBase + bCol + 1);
    Bs[0][bK][bCol + 32] = loadOrZero(pB, N, bK, colBase + bCol + 32);
    Bs[0][bK][bCol + 33] = loadOrZero(pB, N, bK, colBase + bCol + 33);
    __syncthreads();
```

Before the main loop, the first 16 k values (k = 0 to 15) are loaded straight into buffer 0. Each thread loads its four left values and its four right values. The left values go in with k first, `As[0][k][row]`, which is the transposed storage described above. The barrier makes sure the whole first tile is loaded before anyone reads it.

```cpp
    float4 nextA0, nextA1, nextB0, nextB1;
    int cur = 0;
    for (int t = 0; t < N; t += 16) {
        const int tNext = t + 16;
        const bool hasNext = tNext < N;
```

Four `float4` registers will hold this thread's share of the next tile while the current one is used. `cur` says which of the two shared buffers holds the current tile. The main loop steps along k 16 at a time. `tNext` is where the next tile starts, and `hasNext` is false on the last step.

#### Double buffering

Reading from global memory takes hundreds of cycles. If the threads loaded a tile and then waited for it before computing, they would sit idle for that whole time. Instead, each step starts the loads for the next tile into registers first, then computes on the current tile while those loads travel, and only then writes the new values into the other shared buffer.

```cpp
        if (hasNext) {
            if (vectorLoadA && tNext + 16 <= N) {
                const float2* base = pA + (size_t)(rowBase + aRow) * N + tNext + aK;
                nextA0 = __ldg((const float4*)base);
                nextA1 = __ldg((const float4*)(base + 8));
            } else {
                const float2 a0 = loadOrZero(pA, N, rowBase + aRow, tNext + aK);
                const float2 a1 = loadOrZero(pA, N, rowBase + aRow, tNext + aK + 1);
                const float2 a2 = loadOrZero(pA, N, rowBase + aRow, tNext + aK + 8);
                const float2 a3 = loadOrZero(pA, N, rowBase + aRow, tNext + aK + 9);
                nextA0.x = a0.x; nextA0.y = a0.y; nextA0.z = a1.x; nextA0.w = a1.y;
                nextA1.x = a2.x; nextA1.y = a2.y; nextA1.z = a3.x; nextA1.w = a3.y;
            }
```

Start loading the next left values. On the fast path, `base` points at this thread's first element, and two `float4` loads fetch k = aK, aK + 1 and then k = aK + 8, aK + 9. The cast `(const float4*)` reinterprets the address so 16 bytes are read at once. On the slow path, the four values are loaded one at a time, with zero for anything outside the matrix, and packed into the two `float4` registers by hand. `.x .y` hold the first complex number and `.z .w` the second.

```cpp
            if (vectorLoadB && tNext + 16 <= N) {
                const float2* base = pB + (size_t)(tNext + bK) * N + colBase + bCol;
                nextB0 = __ldg((const float4*)base);
                nextB1 = __ldg((const float4*)(base + 32));
            } else {
                const float2 b0 = loadOrZero(pB, N, tNext + bK, colBase + bCol);
                const float2 b1 = loadOrZero(pB, N, tNext + bK, colBase + bCol + 1);
                const float2 b2 = loadOrZero(pB, N, tNext + bK, colBase + bCol + 32);
                const float2 b3 = loadOrZero(pB, N, tNext + bK, colBase + bCol + 33);
                nextB0.x = b0.x; nextB0.y = b0.y; nextB0.z = b1.x; nextB0.w = b1.y;
                nextB1.x = b2.x; nextB1.y = b2.y; nextB1.z = b3.x; nextB1.w = b3.y;
            }
        }
```

The same for the right values. The fast path loads columns bCol, bCol + 1 and then bCol + 32, bCol + 33.

```cpp
        #pragma unroll
        for (int k = 0; k < 16; ++k) {
            float2 a[4];
            float2 b[4];
            #pragma unroll
            for (int m = 0; m < 4; ++m) { a[m] = As[cur][k][ty + 16 * m]; }
            #pragma unroll
            for (int n = 0; n < 4; ++n) { b[n] = Bs[cur][k][tx + 16 * n]; }
```

Now the compute on the current tile, for each of its 16 k values. The thread reads the 4 left values for its rows (ty, ty + 16, ty + 32, ty + 48) and the 4 right values for its columns (tx, tx + 16, tx + 32, tx + 48) from shared memory into registers. In a warp, the threads reading `As` share one or two `ty` values, so they read the same address and it is broadcast. The threads reading `Bs` have neighbouring `tx` values and read neighbouring entries. Neither pattern causes bank conflicts.

```cpp
            #pragma unroll
            for (int m = 0; m < 4; ++m) {
                #pragma unroll
                for (int n = 0; n < 4; ++n) {
                    acc[m][n].x += a[m].x * b[n].x - a[m].y * b[n].y;
                    acc[m][n].y += a[m].x * b[n].y + a[m].y * b[n].x;
                }
            }
        }
```

Every left value times every right value, 16 complex multiply-adds, each added into its own running total. This is where the kernel spends its time. After unrolling, this is 256 multiply-adds per tile step with no loop overhead.

```cpp
        if (hasNext) {
            const int next = cur ^ 1;
            As[next][aK][aRow] = make_cuComplex(nextA0.x, nextA0.y);
            As[next][aK + 1][aRow] = make_cuComplex(nextA0.z, nextA0.w);
            As[next][aK + 8][aRow] = make_cuComplex(nextA1.x, nextA1.y);
            As[next][aK + 9][aRow] = make_cuComplex(nextA1.z, nextA1.w);
            *(float4*)&Bs[next][bK][bCol] = nextB0;
            *(float4*)&Bs[next][bK][bCol + 32] = nextB1;
            __syncthreads();
            cur = next;
        }
    }
```

After the compute, the prefetched values have arrived, and they are written into the other buffer. `cur ^ 1` flips 0 to 1 or 1 to 0. The left values are unpacked one complex number at a time into the transposed layout. The right values go in as whole `float4`s. `&Bs[next][bK][bCol]` is the address of that slot, `(float4*)` treats it as the address of a `float4`, and the `*` in front stores through it, 16 bytes at once. Then one barrier, and the buffers swap roles.

One barrier per step is enough because of the two buffers. During a step, every thread reads only `cur` and writes only `next`, so there is no clash within a step. The buffer being written was last read in the previous step, and the barrier at the end of that step guarantees every thread had finished reading it.

```cpp
    float2* pC = (float2*)C;
    #pragma unroll
    for (int m = 0; m < 4; ++m) {
        #pragma unroll
        for (int n = 0; n < 4; ++n) {
            const int r = rowBase + ty + 16 * m;
            const int c = colBase + tx + 16 * n;
            if (r < N && c < N) { pC[(size_t)r * N + c] = acc[m][n]; }
        }
    }
}
```

When every k step is done, each thread writes its 16 results to global memory. `r` and `c` are the output's row and column in the whole matrix. For each (m, n), the 16 threads that share a `ty` write 16 neighbouring entries of one row, so the writes are coalesced. Outputs outside the matrix are skipped.

### 3.6 matrixMultiply_GPU, the host function GradeBot calls

```cpp
__host__ int matrixMultiply_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C, int* flags, int flagCount){
if (N<=0) { return STUDENTID;}//Your code must be able to deal with N=0 scenario without crashing.
```

The provided signature and guard. `__host__` means this function runs on the CPU.

```cpp
const size_t bytes = sizeof(floatTypeCUDA) * (size_t)N * N;
ensureDevice(N);
cudaMemcpy(dLeft, B, bytes, cudaMemcpyHostToDevice);
cudaMemcpy(dRight, A, bytes, cudaMemcpyHostToDevice);
```

Work out the size of one matrix in bytes, make sure the GPU buffers exist, and copy the inputs to the GPU. `cudaMemcpy(destination, source, bytes, direction)` copies memory. B goes in as the left input and A as the right, the swap from 1.4. Timing runs at several sizes showed that everything growing with N squared, the copies included, took only about 0.7 ms of a 5.9 ms call (job 568780). Copying this much over the normal CPU to GPU link would take longer than that, which suggests GradeBot's buffers already live on the GPU.

```cpp
int kernelChoice = kernelChoiceDefault;
if (flags && flagCount > 0 && flags[0] >= 0 && flags[0] <= 2) { kernelChoice = flags[0]; }
```

Use the default kernel unless the first optional flag asks for a different one. As on the CPU, the pointer, the count and the value are all checked. The flag was only used for comparison runs.

```cpp
if (kernelChoice == 1) {
    dim3 threads(32, 32);
    dim3 blocks((N + 31) / 32, (N + 31) / 32);
    matrixMultiplyKernel_GPU<<<blocks, threads>>>(N, dLeft, dRight, dOut, 0, 0, 0);
} else if (kernelChoice == 2) {
    dim3 threads(32, 32);
    dim3 blocks((N + 31) / 32, (N + 31) / 32);
    matrixMultiplyKernelTiled_GPU<<<blocks, threads>>>(N, dLeft, dRight, dOut);
} else {
    dim3 threads(16, 16);
    dim3 blocks((N + 63) / 64, (N + 63) / 64);
    matrixMultiplyKernelReg_GPU<<<blocks, threads>>>(N, dLeft, dRight, dOut);
}
```

Launch the chosen kernel. `dim3` is a small struct holding up to three sizes. `threads(32, 32)` asks for blocks of 32 × 32 = 1024 threads. `(N + 31) / 32` is N divided by 32, rounded up, so the blocks cover the whole matrix even when N is not a multiple of 32. The launch syntax `kernel<<<blocks, threads>>>(arguments)` starts the kernel with that grid and block size.

Choice 1 is the naive kernel, which takes three extra flag parameters that are not used, so it gets zeros. Choice 2 is the tiled kernel. The default uses 16 × 16 = 256 threads per block, each block covering a 64 by 64 tile, so the grid is N / 64 rounded up in each direction. For N = 2048 that is 32 × 32 = 1024 blocks.

```cpp
cudaMemcpy(C, dOut, bytes, cudaMemcpyDeviceToHost);

return STUDENTID;
}
```

Copy the answer back. A kernel launch returns straight away while the GPU works in the background. This copy waits for the kernel to finish before it runs, because both go through the same default queue of GPU work, so no separate wait is needed.

### 3.7 matrixMultiplyKernel_GPU, the naive kernel (choice 1)

The provided header declares this kernel, so it has to exist. It is kept as the simple baseline.

```cpp
__global__ void matrixMultiplyKernel_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C, int flag0, int flag1, int flag2){

const int j = blockIdx.x * blockDim.x + threadIdx.x;
const int i = blockIdx.y * blockDim.y + threadIdx.y;
if (i >= N || j >= N) { return; }
```

One thread per output entry. The fast-changing x index is used for the column j, so the 32 threads of a warp handle 32 neighbouring columns. Their reads of B and their writes to C then hit neighbouring addresses and are coalesced. Threads past the edge of the matrix stop immediately.

```cpp
float re = 0.0f;
float im = 0.0f;
for (int k = 0; k < N; ++k) {
    const floatTypeCUDA l = A[(size_t)i * N + k];
    const floatTypeCUDA r = B[(size_t)k * N + j];
    re += l.x * r.x - l.y * r.y;
    im += l.x * r.y + l.y * r.x;
}

C[(size_t)i * N + j] = make_cuComplex(re, im);
}
```

The textbook sum from 1.1, straight from global memory. Row i of the left input times column j of the right input, with the complex multiply written out by hand, then one write of the result.

---

## Part 4. matrixMultiplyMPI.cpp, the two-node version

### 4.1 The plan

MPI starts one copy of the program on each node. Each copy is a rank, numbered from 0, and together they form `MPI_COMM_WORLD`. GradeBot gives every rank its own full copy of A and B. The spec also says every node must end up with the full answer, because the benchmark feeds each result into the next multiply.

So the plan is short. Split the rows of the answer between the ranks. Each rank computes its rows with the CPU code from Part 2, using its own 4 threads. Then one MPI call swaps the pieces so every rank has the whole matrix. Using MPI between nodes and OpenMP inside each node is called a hybrid approach.

### 4.2 The top of the file

```cpp
#include <matrixMultiplyMPI.h>
#define STUDENTID 48008361 //DO NOT REMOVE

void matrixMultiplyRowRange(int N, const floatType* A, const floatType* B, floatType* C, int rowBegin, int rowEnd, int* args, int argCount);
```

The provided header brings in `matrixMultiply.h` and `mpi.h`. The third line is a declaration. It tells the compiler that a function with this name and these parameters exists somewhere, without giving its body. The body is in `matrixMultiply.cpp`. Both files are compiled separately and then joined into one program by the linker, which connects this call to that body. That is why `matrixMultiplyRowRange` is not `static` in the CPU file.

```cpp
static const int maxRanks = 256;
```

An upper limit on the number of ranks, so the arrays below can have a fixed size. The graded run uses 2.

### 4.3 matrixMultiply_MPI line by line

```cpp
int matrixMultiply_MPI(int N, const floatType* A, const floatType* B, floatType* C, int* flags, int flagCount){
if (N<=0) { return STUDENTID;}//Your code must be able to deal with N=0 scenario without crashing.
```

The provided signature and guard. The guard must come before any MPI call, because GradeBot calls with N = 0 to read the student ID, possibly before MPI has been started.

```cpp
int mpiReady = 0;
MPI_Initialized(&mpiReady);
```

`MPI_Initialized` is one of the few MPI functions that is safe to call before MPI is started. It sets `mpiReady` to 1 if `MPI_Init` has already run. It takes `&mpiReady`, the address of the variable, so it can write the answer into it.

```cpp
int worldSize = 1;
if (mpiReady) { MPI_Comm_size(MPI_COMM_WORLD, &worldSize); }
```

Ask MPI how many ranks there are. If MPI is not running, assume 1.

```cpp
if (!mpiReady || worldSize <= 1 || worldSize > maxRanks || (long long)N * N * 2 > 2147483647LL) {
    matrixMultiplyRowRange(N, A, B, C, 0, N, flags, flagCount);
    return STUDENTID;
}
```

A safe fallback. If MPI is not running, there is only one rank, there are too many ranks for the arrays, or the matrix is so big that its float count would not fit in an `int` (MPI counts are `int`s, which top out at 2,147,483,647), then just compute every row here and return. This is slower but always correct. The `(long long)` cast does the multiply in a 64-bit type so the check itself cannot overflow. None of these cases happen on the graded run.

```cpp
int rank = 0;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
```

Ask MPI which rank this copy of the program is, 0 or 1 on the graded run.

```cpp
const int rowsEach = N / worldSize;
const int extraRows = N % worldSize;
```

Split the N rows as evenly as possible. Every rank gets `N / worldSize` rows (integer division rounds down), and the leftover rows, `N % worldSize` of them, go one each to the first ranks. For N = 2048 and 2 ranks, each rank gets 1024 rows. For N = 5 and 3 ranks, the ranks get 2, 2 and 1 rows.

```cpp
static int recvCounts[maxRanks];
static int recvOffsets[maxRanks];
for (int p = 0; p < worldSize; ++p) {
    const int firstRow = p * rowsEach + (p < extraRows ? p : extraRows);
    const int rowCount = rowsEach + (p < extraRows ? 1 : 0);
    recvOffsets[p] = firstRow * N * 2;
    recvCounts[p] = rowCount * N * 2;
}
```

Work out every rank's piece, because the gather below needs to know where each piece goes. For rank p, `firstRow` is p times the base size, plus one for each earlier rank that got an extra row. `rowCount` is the base size, plus one if this rank gets an extra row.

The two arrays store each piece's size and starting position in C, counted in floats. A row has N complex numbers, which is 2 × N floats, hence the `* 2`. The data is sent as plain floats rather than an MPI complex type, because plain floats cannot go wrong. The arrays are `static` so they are not rebuilt in memory on every call.

```cpp
const int myFirstRow = rank * rowsEach + (rank < extraRows ? rank : extraRows);
const int myRowCount = rowsEach + (rank < extraRows ? 1 : 0);
matrixMultiplyRowRange(N, A, B, C, myFirstRow, myFirstRow + myRowCount, flags, flagCount);
```

The same formula for this rank's own rows, then compute them with the CPU code. The results land directly in C at the right rows. The other rows of C are not touched yet. On two nodes, each node does half the rows, so this step takes about half as long as the full CPU version. Each node still unpacks the whole right matrix, because every output row needs all of it.

```cpp
MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
               C, recvCounts, recvOffsets, MPI_FLOAT, MPI_COMM_WORLD);

return STUDENTID;
}
```

The exchange. `MPI_Allgatherv` is a collective call, meaning every rank calls it together. Every rank contributes its piece and receives every other rank's piece, so afterwards all ranks hold the whole of C.

- The first three arguments describe what this rank sends. `MPI_IN_PLACE` means this rank's piece is already sitting in C at the right place, so MPI takes it from there and nothing needs copying locally. The 0 and `MPI_DATATYPE_NULL` are ignored with `MPI_IN_PLACE`.
- `C` is where the pieces go.
- `recvCounts` and `recvOffsets` give each rank's piece size and position. The v at the end of `Allgatherv` means the pieces may have different sizes, which the leftover-row split needs.
- `MPI_FLOAT` says the data is floats, and `MPI_COMM_WORLD` says all ranks take part.

Each piece is a range of whole rows, and whole rows sit next to each other in memory, so each piece is one solid block. That keeps the exchange simple.

The measured two-node run took 0.149 s against 0.240 s for one node (job 568790). That means the exchange of about 17 MB costs only a few milliseconds, and the MPI ratio came out at 0.50 against the 0.6 needed for full marks. A more complex version that overlaps computing with sending was planned but not needed.

---

## Part 5. Practice interview questions

### Walk me through what happens when GradeBot calls your CPU function.

The N = 0 guard returns the student ID if the matrix is empty. Otherwise it calls the row range function for every row. That function treats B as the left matrix and A as the right because the reference is column-major. It makes sure the scratch buffers exist, then starts 4 OpenMP threads. The threads split both inputs into separate real and imaginary arrays, zero the output arrays, and then go through k in blocks of 384. For each block they pack the right matrix into 16-column strips, then share out every pair of row block and strip, and run an AVX kernel that computes 4 by 16 patches of the answer in registers. Using the 3M trick, each patch takes three real products instead of four. Finally they merge the real and imaginary output arrays back into C.

### Why is B the left operand?

GradeBot's reference calls MKL's `cgemm`, which uses Fortran's column-major layout. Reading a column-major matrix as row-major gives its transpose. Since (AB)ᵀ = BᵀAᵀ, row-major code has to compute B times A to produce the bytes the reference produces. The first version got this wrong and had an error of about 260, far too big to be rounding, which is what pointed to a layout problem.

### Why not just multiply std::complex values?

Under the grader's compile flags, `std::complex<float>` multiplication calls a library function, `__mulsc3`, that handles infinity and NaN cases with extra checks on every multiply. It is very slow and it stops the compiler from vectorising. Writing (ac - bd) and (ad + bc) by hand avoids it. Our inputs are ordinary finite numbers, so those checks buy nothing.

### What is the planar layout and why use it?

It stores all the real parts in one array and all the imaginary parts in another, instead of alternating them. AVX works on 8 neighbouring floats at once. With alternating storage those 8 floats are half real and half imaginary, which need different arithmetic and would need shuffle instructions to separate. In the planar layout, 8 neighbours are all real parts, so one instruction does 8 useful multiplies.

### What is the 3M trick and does it hurt accuracy?

A complex multiply normally needs four real products, ac, bd, ad and bc. The 3M trick computes t1 = ac, t2 = bd and t3 = (a + b)(c + d), then the real part is t1 - t2 and the imaginary part is t3 - t1 - t2. For matrices this means three big real products instead of four, 25 percent less arithmetic. In theory the subtractions can lose precision, but the measured error went down, from 2.3e-08 to about 5e-09, probably because MKL uses a similar method. Both are far under the 1e-7 tolerance.

### What does the 4 by 16 kernel do, and why that size?

It holds a 4 row by 16 column patch of the answer, 64 running totals, in 8 AVX registers. For each k it loads 16 values of the right matrix into 2 registers, then for each of the 4 rows broadcasts that row's single left value across a register and does two multiply-adds. The 16 loaded values are reused by all 4 rows. The size is set by the 16 AVX registers. The kernel uses 11, a bigger patch would not fit, and a smaller one would not keep enough additions in flight to hide the add instruction's delay.

### Why pack the right matrix?

In the planar array, moving from one k to the next means jumping a whole row, 8 KB at N = 2048, so every step would pull in a new cache line and use only a slice of it. Packing copies each 16-column strip of the current k block so the 16 values for each k come one after another. The kernel then reads one straight line through memory, and a strip of 384 k values is 24 KB, small enough to stay in cache while every row group uses it.

### Why pad the matrices?

The kernel always works on 4 rows and 16 columns. Rounding the scratch arrays up to those multiples and filling the padding with zeros means the kernel never needs special code for the edges. The zeros add nothing to any sum, and the padding is never copied into C.

### What do the OpenMP pragmas do?

`#pragma omp parallel` starts the 4 threads once for the whole function. Each `#pragma omp for` splits the next loop's iterations between them, and `schedule(static)` cuts them into equal chunks up front, which suits loops where every iteration costs the same. Each `omp for` ends with a hidden barrier, so for example no thread starts computing until every strip is packed. `collapse(2)` merges the row block loop and the strip loop into one list of 2816 pairs, which divides evenly between 4 threads. The row blocks alone, 22 of them, would not.

### Why are the buffers static, and is that safe?

Allocating tens of megabytes on every call would waste time, and GradeBot calls the function repeatedly with the same N. The buffers are allocated once and reused. It is safe because the allocation happens before the threads start, and GradeBot never calls the function from two places at once.

### On the GPU, what problem does shared memory solve?

Global memory is slow. In the naive kernel every thread reads a full row and a full column from it, and neighbouring threads read the same values over and over. In the tiled kernel a block loads a tile of each input into shared memory once, and all the threads in the block read from that fast copy. Each loaded value is used 32 times.

### Why is the register tiled kernel faster than the tiled kernel?

In the tiled kernel each thread does one complex multiply-add per two shared memory reads, so shared memory limits the speed. In the register kernel each thread owns 16 outputs in registers. For each k it reads 4 values from each tile and does 16 multiply-adds, which is four times more arithmetic per read. That moved the ratio from 1.58 to 1.32.

### What is double buffering and why is one barrier per step enough?

There are two copies of each shared tile. During a step, the threads start loading the next tile from global memory into registers, compute on the current tile while those loads are in flight, and then write the new values into the other copy. Every thread reads only the current copy and writes only the other one during a step, so they never clash within a step. The copy being written was last read in the step before, and the barrier at the end of that step guarantees all those reads had finished.

### What is coalescing and where does your code rely on it?

When the 32 threads of a warp access neighbouring addresses, the GPU merges them into a few wide memory transactions. The naive kernel maps the x thread index to columns so a warp reads 32 neighbouring B values and writes 32 neighbouring C values. The register kernel spaces each thread's columns 16 apart so, at each store, 16 threads write 16 neighbouring entries of a row, and its loads are arranged so neighbouring threads read neighbouring addresses.

### Why the float4 loads, and when can't you use them?

A `float4` load reads 16 bytes, two complex numbers, in one instruction, which halves the number of load instructions. It needs the address to be a multiple of 16, which holds when N is even, and the tile must be fully inside the matrix. Otherwise the code falls back to single loads that return zero past the edge.

### How is the MPI version split up?

Each rank computes a slice of the rows of the answer using the CPU code with its own 4 threads. The rows are split evenly and leftover rows go one each to the first ranks. Then `MPI_Allgatherv` with `MPI_IN_PLACE` gives every rank every other rank's rows, so all ranks end with the full answer, which the spec requires. The counts are in floats, two per complex number.

### Why does the N = 0 guard come before any MPI call?

GradeBot reads the student ID by calling each function with N = 0 and null pointers, possibly before MPI is running. Calling an MPI function before `MPI_Init` is an error, so the guard has to return first. After the guard, the code still checks `MPI_Initialized` and falls back to a single-node computation if MPI is not running.

### How did you check the answers were correct?

GradeBot compares every result with MKL or cuBLAS and fails anything with a summed squared error above 1e-7. Before using the cluster, the CPU code was also compiled for x86 and run on the Mac under Rosetta against a simple double precision reference at awkward sizes like 1, 3, 5, 33 and 257, to test the padding. The Mac has no NVIDIA GPU, so the GPU kernels were copied into plain C++ with loops standing in for the threads and checked the same way. No version was benchmarked until it passed.

### What would you improve next?

The GPU kernel runs at about 13 TFLOP/s against cuBLAS's 15. Closing that gap would take more tuning of the kernel itself. The compile target rules out the A100's tensor cores. On the CPU, the grader's flags rule out fused multiply-add, which the kernel would otherwise use. The MPI version could overlap the exchange with the computation, but the measured exchange cost only a few milliseconds, so it was not worth the added complexity.

### Did you use AI tools?

Yes, and each file declares it at the top. Claude helped with debugging, drafting pseudocode and implementing it. The design choices were checked with measurements on the cluster, and the changelog records what was tried and why.

---

## Part 6. Glossary

| Term | Meaning |
|---|---|
| AVX | The CPU's 256-bit vector instructions. One instruction works on 8 floats. |
| Bank conflict | Several GPU threads hitting different addresses in the same shared memory bank, which forces them to take turns. |
| Barrier | A point where threads wait until all of them arrive. `__syncthreads()` on the GPU, and the end of every `omp for` on the CPU. |
| Block (GPU) | A group of threads that can share shared memory and use `__syncthreads()`. |
| Broadcast | Copying one value into every lane of a vector register (`_mm256_broadcast_ss`), or one shared memory value being sent to many GPU threads at once. |
| Cache | Small fast memory near the core holding recently used data. |
| Cache blocking | Working on pieces of the data small enough to stay in cache while they are reused. |
| Cache line | The 64-byte chunk that memory moves in. |
| Coalescing | The GPU merging a warp's accesses to neighbouring addresses into a few wide transactions. |
| Column-major | Storing a matrix column by column. Fortran and MKL's `cgemm` use it. |
| FMA | Fused multiply-add, one instruction for a × b + c. Not allowed by the grader's `-mavx` flag. |
| Global memory | The GPU's main memory. Large and slow. |
| Interleaved | Real and imaginary parts stored alternately. |
| Kernel (GPU) | A function that runs on the GPU as many threads at once. |
| Kernel (CPU) | Here, the small innermost function that does most of the arithmetic, the 4 by 16 tile code. |
| MPI | A library for programs running on several machines to exchange messages. |
| OpenMP | A library for sharing loops across the threads of one machine using `#pragma omp`. |
| Packing | Copying data into the exact order the kernel will read it. |
| Padding | Extra zero rows or columns that make sizes a convenient multiple. |
| Planar layout | Real parts in one array, imaginary parts in another. |
| PTX | NVIDIA's intermediate GPU code, turned into machine code for the actual GPU when the program runs. |
| Rank | One copy of an MPI program. |
| Ratio | Our time divided by the reference time. Lower is better. |
| Register | Storage inside the core itself, the fastest there is. |
| Row-major | Storing a matrix row by row. The normal C and C++ order. |
| Shared memory | A small fast GPU memory shared by the threads of one block. |
| SIMD | Single instruction, multiple data. One instruction working on several values. |
| Template | C++ code the compiler copies for each compile-time value it is used with. |
| Thread | One independent stream of instructions. |
| Tile | A small square or rectangle of a matrix handled as one unit. |
| Warp | 32 GPU threads that run each instruction together. |
