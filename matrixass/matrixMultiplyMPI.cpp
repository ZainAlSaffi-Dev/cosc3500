// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

// Provided header, brings in mpi.h
#include <matrixMultiplyMPI.h>
#define STUDENTID 48008361 //DO NOT REMOVE

// Defined in matrixMultiply.cpp, computes rows [rowBegin, rowEnd) of C and leaves the other rows alone.
// Declared here so this file can call it
void matrixMultiplyRowRange(int N, const floatType* A, const floatType* B, floatType* C, int rowBegin, int rowEnd, int* args, int argCount);

// Upper bound on ranks so the count arrays can be static
static const int maxRanks = 256;
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
int matrixMultiply_MPI(int N, const floatType* A, const floatType* B, floatType* C, int* flags, int flagCount){	  	 	 	  		   	  			 		 	      
if (N<=0) { return STUDENTID;}//Your code must be able to deal with N=0 scenario without crashing.	  	 	 	  		   	  			 		 	      

// Each rank computes its own rows of C, then an allgather copies every rank's rows to every rank.
// The N<=0 return above happens before any MPI call.
// Becomes 1 if MPI_Init has run
int mpiReady = 0;
// Check this before any other MPI call
MPI_Initialized(&mpiReady);
// Assume one rank for now
int worldSize = 1;
// Number of ranks
if (mpiReady) { MPI_Comm_size(MPI_COMM_WORLD, &worldSize); }
// Fall back to one rank if MPI is not running or the counts would overflow an int.
// N*N*2 floats must fit in an int
if (!mpiReady || worldSize <= 1 || worldSize > maxRanks || (long long)N * N * 2 > 2147483647LL) {
    // Just compute every row here
    matrixMultiplyRowRange(N, A, B, C, 0, N, flags, flagCount);
    return STUDENTID;
}

// This process's rank
int rank = 0;
// Ask MPI for it
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// N / worldSize rows per rank, the first N % worldSize ranks get one more.
// Base rows per rank
const int rowsEach = N / worldSize;
// Leftover rows
const int extraRows = N % worldSize;
// Counted in floats, two per complex number
static int recvCounts[maxRanks];
// Where each rank's rows start in C, also in floats
static int recvOffsets[maxRanks];
// Work out every rank's share
for (int p = 0; p < worldSize; ++p) {
    // First row of rank p
    const int firstRow = p * rowsEach + (p < extraRows ? p : extraRows);
    // How many rows rank p gets
    const int rowCount = rowsEach + (p < extraRows ? 1 : 0);
    // Offset in floats
    recvOffsets[p] = firstRow * N * 2;
    // Count in floats
    recvCounts[p] = rowCount * N * 2;
}

// My first row, same formula
const int myFirstRow = rank * rowsEach + (rank < extraRows ? rank : extraRows);
// How many rows I get
const int myRowCount = rowsEach + (rank < extraRows ? 1 : 0);
// Compute only my rows using the OpenMP code
matrixMultiplyRowRange(N, A, B, C, myFirstRow, myFirstRow + myRowCount, flags, flagCount);

// In place since my rows are already in C, counts are in floats
MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
               C, recvCounts, recvOffsets, MPI_FLOAT, MPI_COMM_WORLD);

return STUDENTID;	  	 	 	  		   	  			 		 	      

}	  	 	 	  		   	  			 		 	      
