#ifndef MAIN_H
#define MAIN_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

extern int main(int argc, const char * const argv[]);

// This is based on a desired epsilon of 05
// ERROR_BOUND = min(log(1+epsilon), -(1-epsilon))
// see error analysis in the paper
#define _ERROR_BOUND 0.4

//#define _PARALLEL_SEEDS

#define TOL_FUN 0.0000001

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