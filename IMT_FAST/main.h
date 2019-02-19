#ifndef MAIN_H
#define MAIN_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif


extern int main(int argc, char * argv[]);


#ifdef _OPENMP 
#include "omp.h"
//#define _PARALLEL_SEEDS
#define _PARALLEL_CONV
#endif


//#define _OLDCONV



// _GOFAST is for debugging purposes, it sets convergence criteria excessively loose and should not be used when you want useful results!
//#define _GOFAST

#define TOL_FUN 0.001
#define TOL_X 0.001




// This is based on a desired epsilon of 05
// ERROR_BOUND = min(log(1+epsilon), -(1-epsilon))
// see error analysis in the paper
#ifdef _GOFAST
#define _ERROR_BOUND 1
#define TOL_SIZE 1
#else
#define _ERROR_BOUND 0.4
#define TOL_SIZE 0.001
#endif

//#define _DIST_SINGLE

#ifdef _DIST_SINGLE
typedef float distType;
#define distMin FLT_MIN
#define distMax FLT_MAX

#else
typedef double distType;
#define distMin DBL_MIN
#define distMax DBL_MAX
#endif

#ifdef __INTEL_COMPILER
#define MALLOC(size) _mm_malloc(size, 32)
#define FREE(ptr) _mm_free(ptr)
#else
#define MALLOC(size) malloc(size)
#define FREE(ptr) free(ptr)
#endif

typedef struct  {
	distType * data;
	int data_size;

} configStruct;


#endif

//#define _ENABLE_FUNCTION_TRACING

// This stores the current function depth for printing function tracing information
int traceDepth;


#define _MULTI_ADAPT
