#include "imt_analysis.h"
#include "nn_conv.h"
#include "float.h"
#include "stdio.h"
#include "conv.h"
#include "window_conv.h"
#include "loglikelihood.h"
#include "main.h"
#include "utility.h"

//#define _VERBOSE
#define _FASTIDXMETHOD
//#define _SLOWIDXMETHOD
#define _WINDOW_MODE

/*
Evaluating the algebraic convolution of vectors y and z of size size_xyz,
for each point in data[] of size size_XY,
find the value at the point, on the grid x with spacing h, nearest to the data point,
into Y[] and find the loglikelihood logP0

FIXME this really needs to be cleaned up...
*/
void nn_conv(const distType z[], const distType y[], const double X[], const double x[], distType Y[], int size_xyz, int size_XY, double h)
{
	beginTraceFun("nn_conv");

    int allZero = 1;
    for (int i = 0; i < size_XY; i++) {
	if (X[i] > 0)
	    allZero = 0;
    }
    if (allZero)
	printf("X[] is not all positive!\n");

    // Initialize C big enough for the convolution
    int size_conv = 2 * size_xyz;


	distType *C = (distType *)MALLOC(size_conv * sizeof(distType));


    for (int i = 0; i < size_conv; i++) {
		C[i] = 0;
    }

#ifdef _WINDOW_MODE
	window_conv(y, z, C, h, size_xyz);
#else
    conv(y,z,C,h,size_xyz);
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


#ifdef _VERBOSE
    printf("approxconvolv_replacement_indices = \n");
#endif
    double delta, minDelta;
    int nearestNeighborIdx, foundIdx;
    for (int i = 0; i < size_XY; i++) {

#ifdef _FASTIDXMETHOD
	// Newer faster method
	int computedIdx;
	int computedIdxLower = (int) floor(X[i] / h);
	int computedIdxUpper = (int) floor(X[i] / h) + 1;
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
	for (int k = 0; k < size_xyz; k++) {
	    delta = fabs(X[i] - x[k]);
	    if (delta < minDelta) {
		foundIdx = k;
		minDelta = delta;
	    }
	}
#endif

#ifdef _FASTIDXMETHOD
#ifdef _SLOWIDXMETHOD
	if (foundIdx != computedIdx)
	    printf
		("WARNING WARNING WARNING INDEX METHODS DIFFER WARNING WARNING WARNING\n");
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

	if (X[i] > 0 && nearestNeighborIdx > 0
	    && nearestNeighborIdx < size_conv) {
	    // This is from matlab
	    //P(i)=C(I(i)-1);

	    // This is what I wrote the first time
	    //Y[i] = C[computedIdx];

	    // This mirrors what the matlab translated does
	    //Y[i] = C[computedIdx-2];

	    Y[i] = C[nearestNeighborIdx - 1];

	    // WHY????
	} else if (X[i] == 0) {
	    //Y[i] = DBL_MIN;
		Y[i] = 0;
	} else {
	    printf
		("WARNING OUT OF RANGE X[%d]=%f nearestNeighborIdx=%d size_conv=%d size_XY=%d\n",
		 i, X[i], nearestNeighborIdx, size_conv, size_XY);
	    //Y[i] = DBL_MIN;
		Y[i] = 0;
	}
    }
#ifdef _VERBOSE
    printf("\n\n");
#endif

	for (int i = 0; i < size_XY; i++){
		if (Y[i] * h > 1) {
			printf("ERROR: nn_conv(): Y[%d]*%f = %f > 1 ; Y[%d] = %f\n", i, h, Y[i]*h, i, Y[i]);
			exit(1);
		}
	}


	//*logP0 = (double) loglikelihood(Y, size_XY);


#ifdef _VERBOSE
    printf("approxconvolv_replacement = \n");
    for (int i = 0; i < size_XY; i++) {
	if (i % 8 == 0)
	    printf("\n");
	printf("%.17f ", Y[i]);
    }
    printf("\n\n");
#endif

	FREE(C);

	endTraceFun("nn_conv");

}