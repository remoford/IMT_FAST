#include "imt_analysis.h"
#include "binned_conv.h"
#include "float.h"
#include "stdio.h"
#include "conv.h"
#include "window_conv.h"
#include "main.h"
#include "loglikelihood.h"

//#define _VERBOSE
//#define _FASTIDXMETHOD
#define _SLOWIDXMETHOD
#define _WINDOW_MODE
//#define _LAZY_MODE

void
binned_conv(const distType z[], const distType y[],
			  const double X[], const double x[], distType Y[],
			  double *logP0, int size_xyz, int size_XY,
			  double h)
{
	// Sanity checking
    int allZero = 1;
    for (int i = 0; i < size_XY; i++) {
	if (X[i] > 0)
	    allZero = 0;
    }
    if (allZero)
	printf("X[] is not all positive!\n");

    // Initialize C big enough for the convolution
    int size_conv = 2 * size_xyz;

	// C stores the result of the convolution
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


#ifdef _LAZY_MODE

	int maxX = 0;
	for (int i = 0; i < size_XY; i++) {
		if (X[i] > maxX)
			maxX = X[i];
	}

	int maxIdx = maxX * 10;

	int * weights = malloc(sizeof(int)*maxIdx);

	int * value = malloc(sizeof(int)*maxIdx);
	
	for (int i = 0; i <= maxIdx; i++)
		weights[i] = 0;
	
	for (int i = 0; i < size_XY; i++)
		weights[(int)(X[i] * 10)]++;

	/*
	printf("(");
	for (int i = 0; i <= maxIdx; i++)
		printf("%d ", weights[i]);
	printf(") ");
	*/

	for (int i = 0; i <= maxIdx; i++) {
		value[i] = 0;
		if (weights[i] > 0) {
			double dataPoint = ( (double) i ) / 10.0;

			// rightmost boundry of integration for this particular bin
			int rightBound = (int)(dataPoint / h);

			// the bin width in terms of indices
			int goback = (int)(0.1 / h);

			// the leftmost boundry of integration
			int leftBound = rightBound - goback;

			// Calculate the right handed riemann sum
			
			for (int j = leftBound + 1; j <= rightBound; j++) 
				value[i] += C[j] * h;
		}
	}

	for (int i = 0; i < size_XY; i++)
		Y[i] = value[(int)(X[i] * 10)];
	
	free(weights);
	free(value);

#else

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
	
#endif

	*logP0 = (double) loglikelihood(Y, size_XY);
}
