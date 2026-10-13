// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

// Provided header, brings in cuComplex and the kernel declaration
#include <matrixMultiplyGPU.cuh>
#define STUDENTID 48008361 //DO NOT REMOVE

// Kernel used when no flags are passed. 0 register tiled, 1 one thread per element, 2 shared memory tiled.
// Register tiled was the fastest
static const int kernelChoiceDefault = 0;

// Device buffers, kept between calls and reallocated only when N changes.
// N the device buffers were sized for, -1 means none yet
static int deviceN = -1;
// Left matrix on the GPU
static floatTypeCUDA* dLeft = 0;
// Right matrix on the GPU
static floatTypeCUDA* dRight = 0;
// Result on the GPU
static floatTypeCUDA* dOut = 0;

// Allocate device memory for this N
static void ensureDevice(int N)
{
    // Same N as last time, reuse
    if (N == deviceN) { return; }
    // Free the old left buffer, freeing null is fine
    cudaFree(dLeft);
    // Free the old right buffer
    cudaFree(dRight);
    // Free the old output buffer
    cudaFree(dOut);
    // Bytes in one N by N complex matrix
    const size_t bytes = sizeof(floatTypeCUDA) * (size_t)N * N;
    // New left buffer
    cudaMalloc((void**)&dLeft, bytes);
    // New right buffer
    cudaMalloc((void**)&dRight, bytes);
    // New output buffer
    cudaMalloc((void**)&dOut, bytes);
    // Remember the size
    deviceN = N;
}

// Row-major C = A*B with 32 by 32 tiles loaded through shared memory.
// 32 by 32 threads, one output element each
__global__ void matrixMultiplyKernelTiled_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C)
{
    // Extra column avoids bank conflicts
    __shared__ floatTypeCUDA As[32][33];
    // Same padding for B
    __shared__ floatTypeCUDA Bs[32][33];
    // Column inside the tile
    const int tx = threadIdx.x;
    // Row inside the tile
    const int ty = threadIdx.y;
    // Global column of C
    const int j = blockIdx.x * 32 + tx;
    // Global row of C
    const int i = blockIdx.y * 32 + ty;
    // Running real sum
    float re = 0.0f;
    // Running imaginary sum
    float im = 0.0f;
    // Move along k one tile at a time
    for (int t = 0; t < N; t += 32) {
        // Load one tile of A and one of B, zero outside the matrix.
        // Each thread loads one element of A, zero if it is past the edge
        As[ty][tx] = (i < N && t + tx < N) ? A[(size_t)i * N + t + tx]
                                           : make_cuComplex(0.0f, 0.0f);
        // And one element of B, zero if it is past the edge
        Bs[ty][tx] = (t + ty < N && j < N) ? B[(size_t)(t + ty) * N + j]
                                           : make_cuComplex(0.0f, 0.0f);
        // Wait until the whole tile is loaded
        __syncthreads();
        // Complex multiply written out by hand
        for (int k = 0; k < 32; ++k) {
            // A value from shared memory
            const floatTypeCUDA l = As[ty][k];
            // B value from shared memory
            const floatTypeCUDA r = Bs[k][tx];
            // Real part, ac - bd
            re += l.x * r.x - l.y * r.y;
            // Imaginary part, ad + bc
            im += l.x * r.y + l.y * r.x;
        }
        // Wait before the tile is overwritten
        __syncthreads();
    }
    // Skip threads outside the matrix
    if (i < N && j < N) {
        // Write the complex result
        C[(size_t)i * N + j] = make_cuComplex(re, im);
    }
}

// Reads one element, or zero when the index is outside the matrix.
// Small device helper, inlined into the kernel
static __device__ __forceinline__ float2 loadOrZero(const float2* P, int N, int r, int c)
{
    // __ldg reads through the read only cache
    if (r < N && c < N) { return __ldg(P + (size_t)r * N + c); }
    // Outside the matrix gives zero
    return make_cuComplex(0.0f, 0.0f);
}

// The default kernel. Each block computes a 64 by 64 tile, each thread computes 4 by 4 of it in registers.
// __restrict__ means the pointers do not overlap
__global__ void matrixMultiplyKernelReg_GPU(int N, const floatTypeCUDA* __restrict__ A, const floatTypeCUDA* __restrict__ B, floatTypeCUDA* __restrict__ C)
{
    // Double buffered, 66 keeps float4 stores aligned and avoids bank conflicts
    __shared__ float2 As[2][16][66];
    // Same layout for B, k then column
    __shared__ float2 Bs[2][16][66];

    // Thread column in the block, 0 to 15
    const int tx = threadIdx.x;
    // Thread row in the block, 0 to 15
    const int ty = threadIdx.y;
    // Flat thread id, 0 to 255
    const int tid = ty * 16 + tx;
    // First row of this block's tile
    const int rowBase = blockIdx.y * 64;
    // First column of this block's tile
    const int colBase = blockIdx.x * 64;

    // Each thread copies two float4 of A and two of B per 16 k values.
    // A block is 64 rows by 16 k
    const int aRow = tid >> 2;
    // Which pair of k values this thread loads
    const int aK = (tid & 3) * 2;
    // B block is 16 k by 64 columns
    const int bK = tid >> 4;
    // Which pair of columns this thread loads
    const int bCol = (tid & 15) * 2;

    // View A as float2, real and imaginary
    const float2* pA = (const float2*)A;
    // Same for B
    const float2* pB = (const float2*)B;

    // float4 loads need the tile fully inside the matrix and an even N for 16 byte alignment.
    // float4 is safe for A
    const bool vectorLoadA = (rowBase + 64 <= N) && ((N & 1) == 0);
    // float4 is safe for B
    const bool vectorLoadB = (colBase + 64 <= N) && ((N & 1) == 0);

    // 4 by 4 complex results kept in registers
    float2 acc[4][4];
    // Fully unroll so acc stays in registers
    // Each row of acc
    #pragma unroll
    for (int m = 0; m < 4; ++m) {
        // Unroll the inner loop too
        // Clear each accumulator
        #pragma unroll
        for (int n = 0; n < 4; ++n) { acc[m][n] = make_cuComplex(0.0f, 0.0f); }
    }

    // The first 16 k values go straight into shared memory.
    // A is stored transposed, k then row
    As[0][aK][aRow] = loadOrZero(pA, N, rowBase + aRow, aK);
    // Same row, k + 1
    As[0][aK + 1][aRow] = loadOrZero(pA, N, rowBase + aRow, aK + 1);
    // Same row, k + 8
    As[0][aK + 8][aRow] = loadOrZero(pA, N, rowBase + aRow, aK + 8);
    // Same row, k + 9
    As[0][aK + 9][aRow] = loadOrZero(pA, N, rowBase + aRow, aK + 9);
    // B is stored k then column
    Bs[0][bK][bCol] = loadOrZero(pB, N, bK, colBase + bCol);
    // Same k, column + 1
    Bs[0][bK][bCol + 1] = loadOrZero(pB, N, bK, colBase + bCol + 1);
    // Same k, column + 32
    Bs[0][bK][bCol + 32] = loadOrZero(pB, N, bK, colBase + bCol + 32);
    // Same k, column + 33
    Bs[0][bK][bCol + 33] = loadOrZero(pB, N, bK, colBase + bCol + 33);
    // First tile is ready
    __syncthreads();

    // Registers holding the next tile
    float4 nextA0, nextA1, nextB0, nextB1;
    // Which shared buffer is being read
    int cur = 0;
    // Step along k 16 at a time
    for (int t = 0; t < N; t += 16) {
        // Start of the next tile
        const int tNext = t + 16;
        // False on the last tile
        const bool hasNext = tNext < N;
        // Only prefetch if there is a next tile
        if (hasNext) {
            // Load the next 16 k values into registers while the current ones are used.
            // Whole tile in range, use float4
            if (vectorLoadA && tNext + 16 <= N) {
                // This thread's first element
                const float2* base = pA + (size_t)(rowBase + aRow) * N + tNext + aK;
                // Two complex values in one load
                nextA0 = __ldg((const float4*)base);
                // Two more, 8 k values later
                nextA1 = __ldg((const float4*)(base + 8));
            } else {
                // Edge tile, load one at a time
                // Element at k
                const float2 a0 = loadOrZero(pA, N, rowBase + aRow, tNext + aK);
                // Element at k + 1
                const float2 a1 = loadOrZero(pA, N, rowBase + aRow, tNext + aK + 1);
                // Element at k + 8
                const float2 a2 = loadOrZero(pA, N, rowBase + aRow, tNext + aK + 8);
                // Element at k + 9
                const float2 a3 = loadOrZero(pA, N, rowBase + aRow, tNext + aK + 9);
                // Pack a0 and a1 into a float4
                nextA0.x = a0.x; nextA0.y = a0.y; nextA0.z = a1.x; nextA0.w = a1.y;
                // Pack a2 and a3
                nextA1.x = a2.x; nextA1.y = a2.y; nextA1.z = a3.x; nextA1.w = a3.y;
            }
            // Same idea for B
            if (vectorLoadB && tNext + 16 <= N) {
                // This thread's first element
                const float2* base = pB + (size_t)(tNext + bK) * N + colBase + bCol;
                // Two complex values in one load
                nextB0 = __ldg((const float4*)base);
                // Two more, 32 columns over
                nextB1 = __ldg((const float4*)(base + 32));
            } else {
                // Edge tile, load one at a time
                // Element at this column
                const float2 b0 = loadOrZero(pB, N, tNext + bK, colBase + bCol);
                // Element at column + 1
                const float2 b1 = loadOrZero(pB, N, tNext + bK, colBase + bCol + 1);
                // Element at column + 32
                const float2 b2 = loadOrZero(pB, N, tNext + bK, colBase + bCol + 32);
                // Element at column + 33
                const float2 b3 = loadOrZero(pB, N, tNext + bK, colBase + bCol + 33);
                // Pack b0 and b1 into a float4
                nextB0.x = b0.x; nextB0.y = b0.y; nextB0.z = b1.x; nextB0.w = b1.y;
                // Pack b2 and b3
                nextB1.x = b2.x; nextB1.y = b2.y; nextB1.z = b3.x; nextB1.w = b3.y;
            }
        }
        // Unroll the 16 k steps
        // Each k in the current tile
        #pragma unroll
        for (int k = 0; k < 16; ++k) {
            // One A value per output row
            float2 a[4];
            // One B value per output column
            float2 b[4];
            // Unroll this loop
            // Rows 16 apart, no bank conflicts
            #pragma unroll
            for (int m = 0; m < 4; ++m) { a[m] = As[cur][k][ty + 16 * m]; }
            // Unroll this loop
            // Columns 16 apart too
            #pragma unroll
            for (int n = 0; n < 4; ++n) { b[n] = Bs[cur][k][tx + 16 * n]; }
            // Unroll this loop
            // Each output row
            #pragma unroll
            for (int m = 0; m < 4; ++m) {
                // Unroll this loop
                // Each output column
                #pragma unroll
                for (int n = 0; n < 4; ++n) {
                    // Real part
                    acc[m][n].x += a[m].x * b[n].x - a[m].y * b[n].y;
                    // Imaginary part
                    acc[m][n].y += a[m].x * b[n].y + a[m].y * b[n].x;
                }
            }
        }
        // Store the prefetched tile
        if (hasNext) {
            // Write the next 16 k values into the other buffer. Reads of it finished before the last barrier.
            // The buffer nobody is reading
            const int next = cur ^ 1;
            // Unpack into A's transposed layout
            As[next][aK][aRow] = make_cuComplex(nextA0.x, nextA0.y);
            // Same row, k + 1
            As[next][aK + 1][aRow] = make_cuComplex(nextA0.z, nextA0.w);
            // Same row, k + 8
            As[next][aK + 8][aRow] = make_cuComplex(nextA1.x, nextA1.y);
            // Same row, k + 9
            As[next][aK + 9][aRow] = make_cuComplex(nextA1.z, nextA1.w);
            // B goes in as a whole float4
            *(float4*)&Bs[next][bK][bCol] = nextB0;
            // The other two columns, 32 over
            *(float4*)&Bs[next][bK][bCol + 32] = nextB1;
            // Next tile is ready for everyone
            __syncthreads();
            // Swap buffers
            cur = next;
        }
    }

    // Each thread's columns are 16 apart so the stores are coalesced.
    // View C as float2
    float2* pC = (float2*)C;
    // Unroll this loop
    // Each output row
    #pragma unroll
    for (int m = 0; m < 4; ++m) {
        // Unroll this loop
        // Each output column
        #pragma unroll
        for (int n = 0; n < 4; ++n) {
            // Global row
            const int r = rowBase + ty + 16 * m;
            // Global column
            const int c = colBase + tx + 16 * n;
            // Only write inside the matrix
            if (r < N && c < N) { pC[(size_t)r * N + c] = acc[m][n]; }
        }
    }
}
/**
* @brief Implements an NxN matrix multiply C=A*B
*	  	 	 	  		   	  			 		 	      
* @param[in] N : dimension of square matrix (NxN)
* @param[in] A : pointer to input NxN matrix
* @param[in] B : pointer to input NxN matrix
* @param[out] C : pointer to output NxN matrix
* @param[in] flags : pointer to array of integers which can be used for debugging and performance tweaks. Optional. If unused, set to zero
* @param[in] flagCount : the length of the flags array
* @return : your student ID
*	  	 	 	  		   	  			 		 	      
* */
__host__ int matrixMultiply_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C, int* flags, int flagCount){	  	 	 	  		   	  			 		 	      
if (N<=0) { return STUDENTID;}//Your code must be able to deal with N=0 scenario without crashing.	  	 	 	  		   	  			 		 	      

// The reference is column-major, which in row-major terms is C = B * A, so B is the left operand.
// Size of one matrix in bytes
const size_t bytes = sizeof(floatTypeCUDA) * (size_t)N * N;
// Allocate once, then reuse
ensureDevice(N);
// Copy B to the GPU as the left operand
cudaMemcpy(dLeft, B, bytes, cudaMemcpyHostToDevice);
// Copy A as the right operand
cudaMemcpy(dRight, A, bytes, cudaMemcpyHostToDevice);

// flags[0] picks the kernel if given, otherwise the default is used.
// Start with the default
int kernelChoice = kernelChoiceDefault;
// Only use flags[0] if it exists and is valid
if (flags && flagCount > 0 && flags[0] >= 0 && flags[0] <= 2) { kernelChoice = flags[0]; }

// One thread per element
if (kernelChoice == 1) {
    // 1024 threads per block
    dim3 threads(32, 32);
    // Enough blocks to cover N
    dim3 blocks((N + 31) / 32, (N + 31) / 32);
    // Flag arguments are unused
    matrixMultiplyKernel_GPU<<<blocks, threads>>>(N, dLeft, dRight, dOut, 0, 0, 0);
} else if (kernelChoice == 2) {
    // Shared memory tiled
    // 32 by 32 threads to match the tile
    dim3 threads(32, 32);
    // One block per 32 by 32 tile
    dim3 blocks((N + 31) / 32, (N + 31) / 32);
    // Launch the shared memory kernel
    matrixMultiplyKernelTiled_GPU<<<blocks, threads>>>(N, dLeft, dRight, dOut);
} else {
    // Register tiled
    // 256 threads, each does 4 by 4 outputs
    dim3 threads(16, 16);
    // One block per 64 by 64 tile
    dim3 blocks((N + 63) / 64, (N + 63) / 64);
    // Launch the register tiled kernel
    matrixMultiplyKernelReg_GPU<<<blocks, threads>>>(N, dLeft, dRight, dOut);
}

// The copy back waits for the kernel, so no separate synchronisation is needed.
// Copy the result back to the host
cudaMemcpy(C, dOut, bytes, cudaMemcpyDeviceToHost);

return STUDENTID;	  	 	 	  		   	  			 		 	      

}	  	 	 	  		   	  			 		 	      

//The kernel (device code) parameters have been setup almost the same as the host code, except the flags are passed in individually rather than as a pointer. This is done just so you don't have to copy the parameters to GPU memory first, you'll be able to pass in up to 3 on the function call.	  	 	 	  		   	  			 		 	      
__global__ void matrixMultiplyKernel_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C, int flag0, int flag1, int flag2){

// One thread per output element, kept as the baseline kernel.
// x along the row so reads of B and the write of C coalesce
const int j = blockIdx.x * blockDim.x + threadIdx.x;
// The y index gives the row
const int i = blockIdx.y * blockDim.y + threadIdx.y;
// Threads past the edge do nothing
if (i >= N || j >= N) { return; }
// Running real sum
float re = 0.0f;
// Running imaginary sum
float im = 0.0f;
// Row i of A times column j of B
for (int k = 0; k < N; ++k) {
    // Element A[i][k]
    const floatTypeCUDA l = A[(size_t)i * N + k];
    // Element B[k][j]
    const floatTypeCUDA r = B[(size_t)k * N + j];
    // Real part
    re += l.x * r.x - l.y * r.y;
    // Imaginary part
    im += l.x * r.y + l.y * r.x;
}
// write the result into C
C[(size_t)i * N + j] = make_cuComplex(re, im);	  	 	 	  		   	  			 		 	      

}	  	 	 	  		   	  			 		 	      
