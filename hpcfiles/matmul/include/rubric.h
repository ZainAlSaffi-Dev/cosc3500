//This is the numeric type you'll be working with (complex64)
#define FLOATTYPE_COMPLEX64

#ifndef RUBRIC_CPU

//All GradeBot CPU memory allocations are aligned to this boundary (i.e. a 64-byte cacheline)
#define ALIGN 64

//enumerations that indicate which benchmarks to run from the command-line.
#define RUBRIC_CPU 0
#define RUBRIC_GPU 1
#define RUBRIC_MPI 2

//Default parameters (for when you don't specify anything on the command-line)
//The default matrix size
#define DEFAULT_N 2048

//The default thread count
#define DEFAULT_THREADCOUNT 4

#ifdef FLOATTYPE_FLOAT64
#define floatType double
#define floatTypeMKL double
#define floatTypeCUDA double
#define ERR_TOLERANCE 1e-20
#endif

#ifdef FLOATTYPE_COMPLEX128
#include <complex>
#include <cuComplex.h>
#define floatType std::complex<double>
#define floatTypeMKL MKL_Complex16
#define floatTypeCUDA cuDoubleComplex
#define ERR_TOLERANCE 1e-20
#endif

#ifdef FLOATTYPE_COMPLEX64
#include <complex>
#include <cuComplex.h>
#define floatType std::complex<float>
#define floatTypeMKL MKL_Complex8
#define floatTypeCUDA cuComplex
#define ERR_TOLERANCE 1e-7
#endif

#ifdef FLOATTYPE_FLOAT32
#define floatType float
#define floatTypeMKL float
#define floatTypeCUDA float
#define ERR_TOLERANCE 1e-7
#endif

//This is the main routine that runs the GradeBot
int GradeBot(int N, int threadCount, int *isEnabled, int* args, int argCount);

#endif