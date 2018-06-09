

 /* Include files */
#include "IMT_analysis_April2017.h"
#include "approxconvolv.h"
#include "float.h"
#include "stdio.h"

//#define _VERBOSE
//#define _FASTIDXMETHOD
#define _SLOWIDXMETHOD

/* Function Definitions */

void approxconvolv_replacement(const double z[], const double y[],
	const double X[], const double x[], double Y[],
	double *logP0, int size_xyz, int size_XY, double h) {
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
	double * C = _mm_malloc(size_conv * sizeof(double), 32);
	//double * C = malloc(size_conv * sizeof(double));
	for (int i = 0; i < size_conv; i++) {
		C[i] = 0;
	}

	// Do the convolution
	/*
	if (size_xyz > 1024*8) {

#pragma omp parallel for
		for (int i = 0; i < size_xyz; i++) {
			for (int j = 0; j < size_xyz; j++) {
				C[i + j] += z[i] * y[j] * h;
			}
		}
	}
	else {*/
		for (int i = 0; i < size_xyz; i++) {
#pragma vector aligned
			for (int j = 0; j < size_xyz; j++) {
				C[i + j] += z[i] * y[j] * h;
			}
		}

	//}

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

	// Find the largest value in X
	int MaxX = 0;
	for (int i = 0; i < size_XY; i++) {
		if (X[i] > MaxX)
			MaxX = X[i];
	}

#ifdef _VERBOSE
	printf("approxconvolv_replacement_indices = \n");
#endif
	double delta, minDelta;
	int nearestNeighborIdx, foundIdx, computedIdx, computedIdxLower, computedIdxUpper;
	for (int i = 0; i < size_XY; i++) {

#ifdef _FASTIDXMETHOD
		// Newer faster method
		computedIdxLower = (int)floor(X[i] / h);
		computedIdxUpper = (int)floor(X[i] / h) + 1;
		if (computedIdxUpper >= size_xyz) {
			//printf("WARNING OUT OF RANGE computedIdxUpper=%d size_xyz=%d\n", computedIdxUpper, size_xyz);
			computedIdxUpper = size_xyz - 1;
		}

		double deltaLower = fabs(X[i] - x[computedIdxLower]);
		double deltaUpper = fabs(X[i] - x[computedIdxUpper]);

		if (deltaLower <= deltaUpper)
			computedIdx = computedIdxLower;
		else
			computedIdx = computedIdxUpper;
#endif

#ifdef _SLOWIDXMETHOD
		// Older slower method
		//%[~, I(i)] = min((t(i) - x). ^ 2); % slow
		minDelta = DBL_MAX;
		foundIdx = 0;
		for (int k = 0; k <= size_xyz; k++) {
			delta = fabs(X[i] - x[k]);
			if (delta < minDelta) {
				foundIdx = k;
				minDelta = delta;
			}
		}
#endif

		/*
		printf("index check X[i]=%f x[%d]=%f x[%d]=%f X[i]/h=%f h=%f computedIdx=%d",
			X[i], computedIdxLower, x[computedIdxLower], computedIdxUpper, x[computedIdxUpper], X[i] / h, h, computedIdx);

		printf(" foundIdx=%d minDelta=%f x[%d]=%f\n", foundIdx, minDelta, foundIdx, x[foundIdx]);
		*/

#ifdef _FASTIDXMETHOD
#ifdef _SLOWIDXMETHOD
		if (foundIdx != computedIdx)
			printf("WARNING WARNING WARNING INDEX METHODS DIFFER WARNING WARNING WARNING\n");
#endif
#endif


#ifdef _SLOWIDXMETHOD
		nearestNeighborIdx = foundIdx;
#endif
#ifdef _FASTIDXMETHOD
		nearestNeighborIdx = computedIdx;
#endif

#ifdef _VERBOSE
		printf("%d ", nearestNeighborIdx);
		if (i % 8 == 0)
			printf("\n");
#endif


		if (X[i] > 0 && nearestNeighborIdx > 0 && nearestNeighborIdx < size_conv){
			// This is from matlab
			//P(i)=C(I(i)-1);

			// This is what I wrote the first time
			//Y[i] = C[computedIdx];

			// This mirrors what the matlab translated does
			//Y[i] = C[computedIdx-2];

			Y[i] = C[nearestNeighborIdx-1];

			// WHY????
		}
		else if (X[i] == 0){
			Y[i] = DBL_MIN;
		}
		else{
			printf("WARNING OUT OF RANGE X[%d]=%f nearestNeighborIdx=%d size_conv=%d size_XY=%d\n", i, X[i], nearestNeighborIdx, size_conv, size_XY);
			Y[i] = DBL_MIN;
		}
	}
#ifdef _VERBOSE
	printf("\n\n");
#endif

	*logP0 = 0;
	for (int i = 0; i < size_XY; i++) {
		*logP0 += log(Y[i]);
	}

	for (int i = 0; i < size_XY; i++) {
		if (Y[i] == 0)
			//Y[i] = 2.2250738585072014E-308;
			Y[i] = DBL_MIN;
	}

#ifdef _VERBOSE
	printf("approxconvolv_replacement = \n");
	for (int i = 0; i < size_XY; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%.17f ", Y[i]);
	}
	printf("\n\n");
#endif
	
	_mm_free(C);
	//free(C);
}




