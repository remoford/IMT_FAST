
 /* Include files */
#include "IMT_analysis_April2017.h"
#include "binned_conv.h"
#include "float.h"
#include "stdio.h"
#include "conv.h"
#include "window_conv.h"
#include "main.h"

//#define _VERBOSE
//#define _FASTIDXMETHOD
#define _SLOWIDXMETHOD
#define _WINDOW_MODE

/* Function Definitions */

void
binned_conv(const distType z[], const distType y[],
			  const double X[], const double x[], distType Y[],
			  double *logP0, int size_xyz, int size_XY,
			  double h)
{
    // function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )

    /*
       printf("X = \n");
       for (int i = 0; i < size_XY; i++) {
       if (i % 8 == 0)
       printf("\n");
       printf("%f ", X[i]);
       }
       printf("\n");
     */

    int allZero = 1;
    for (int i = 0; i < size_XY; i++) {
	if (X[i] > 0)
	    allZero = 0;
    }
    if (allZero)
	printf("X[] is not all positive!\n");

    // Initialize C big enough for the convolution
    int size_conv = 2 * size_xyz;
#ifdef __INTEL_COMPILER
    distType *C = (distType *) _mm_malloc(size_conv * sizeof(distType), 32);
#else
    distType *C = (distType *) malloc(size_conv * sizeof(distType));
#endif
    for (int i = 0; i < size_conv; i++) {
	C[i] = 0;
    }

#ifdef _WINDOW_MODE
    window_conv(y,z,C,h,size_xyz);
#else
	conv(y, z, C, h, size_xyz);
#endif

    // Display the convolution
#ifdef _VERBOSE
    printf("replacement_conv = \n");
    for (int i = 0; i < size_conv; i++) {
	if (i % 8 == 0)
	    printf("\n");
	printf("%.17f ", C[i]);
    }
    printf("\n\n");
#endif

	// Calculate the probability integral over the bin for each point in the data
	for (int i = 0; i < size_XY; i++) {
		// rightmost boundry of integration for this particular bin
		int rightBound = (int) (X[i] / h);

		// the bin width in terms of indices
		int goback = (int) (0.1 / h);

		// the leftmost boundry of integration
		int leftBound = rightBound - goback;

		// Calculate the right handed riemann sum
		Y[i] = 0;
		for (int j = leftBound + 1; j <= rightBound; j++) {
			Y[i] += C[j] * h;
		}

	}

	/*
	printf("Y = \n");
	for (int i = 0; i < size_XY; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%.17f ", Y[i]);
	}
	printf("\n\n");
	*/

	
	// Replace all zero values with DBL_MIN to abvoid unbounded numbers when taking a logarithm
	/*
	for (int i = 0; i < size_XY; i++) {
		if (Y[i] == 0){
			//Y[i] = 2.2250738585072014E-308;
			//printf("OMG OMG OMG OMG WE GOT A NUMBER THTS REALLY DAMN CLOSE TO ZERO!!!!!!!!!!!!!!!!!!!!!!\n");
			Y[i] = distMin;
		}
	}
	*/

	// FIX ME FIX ME FIX ME THIS IS BROKEN AND NEEDS FIXING ALWAYS RETURNS HTE SAME DMANED THING

	// Caluclate the log likelihood of the data given our computed distribution
	*logP0 = 0;
	for (int i = 0; i < size_XY; i++) {
		*logP0 += log(Y[i]);
	}
	//printf("OUR FOUND LOGP0 = %f", *logP0);


	

}
