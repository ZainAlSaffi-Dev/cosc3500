#include <thread> //used for getting the number of cpu cores available (this is the hardware core count, not the slurm count)
#include <omp.h> //openMP library
#include <immintrin.h>//AVX library
#include <string.h> //This is included mostly so you can memset, e.g. for setting arrays to zero, memset(A,0,sizeof(floatType)*N*N);
#include <rubric.h>

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
* */
int matrixMultiply(int N, const floatType* A, const floatType* B, floatType* C, int* args, int argCount);