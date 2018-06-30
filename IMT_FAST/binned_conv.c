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

void binned_conv(const distType z[], const distType y[],
			  const distType data[], const distType x[], distType Y[],
			  double *logP0, int size_xyz, int size_XY,
			  double h)
{
	// Sanity checking
    int allZero = 1;
    for (int i = 0; i < size_XY; i++) {
	if (data[i] > 0)
	    allZero = 0;
    }
    if (allZero)
	printf("data[] is not all positive!\n");

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
		if (data[i] > maxX)
			maxX = data[i];
	}

	int maxIdx = maxX * 10;

	int * weights = malloc(sizeof(int)*maxIdx);

	int * value = malloc(sizeof(int)*maxIdx);
	
	for (int i = 0; i <= maxIdx; i++)
		weights[i] = 0;
	
	for (int i = 0; i < size_XY; i++)
		weights[(int)(data[i] * 10)]++;

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
		Y[i] = value[(int)(data[i] * 10)];
	
	free(weights);
	free(value);

#else

	// Calculate the probability integral over the bin for each point in the data
	for (int i = 0; i < size_XY; i++) {
		// rightmost boundry of integration for this particular bin
		int rightBound = (int) (data[i] / h);

		if (rightBound >= size_conv)
			printf("ERROR: rightBound=%d>=size_conv=%d ", rightBound, size_conv);

		// the bin width in terms of indices
		int goback = (int) (0.1 / h);

		// the leftmost boundry of integration
		int leftBound = rightBound - goback;

		if (leftBound < 0)
			printf("ERROR: leftBound=%d<0 ", leftBound);

		// Calculate the right handed riemann sum
		Y[i] = 0;
		for (int j = leftBound + 1; j < rightBound; j++) {
			Y[i] += C[j] * h;
		}
	}
	
#endif

	*logP0 = (double) loglikelihood(Y, size_XY);

#ifdef __INTEL_COMPILER
	_mm_free(C);
#else
	free(C);
#endif
}



void threestage_binconv(const distType x[], const distType y[], const distType z[], const distType data[], distType Y[], double *logP0, int size_xyz, int dataSize, double h) {

	// Initialize C1 big enough for the first convolution
	int size_conv1 = 2 * size_xyz;

	// C stores the result of the convolution
#ifdef __INTEL_COMPILER
	distType *C1 = (distType *)_mm_malloc(size_conv1 * sizeof(distType), 32);
#else
	distType *C1 = (distType *)malloc(size_conv1 * sizeof(distType));
#endif

	for (int i = 0; i < size_conv1; i++)
		C1[i] = 0;

	window_conv(x, y, C1, h, size_xyz);

	// Copy z pdf into an expanded array
#ifdef __INTEL_COMPILER
	distType *expandedZ = (distType *)_mm_malloc(size_conv1 * sizeof(distType), 32);
#else
	distType *expandedZ = (distType *)malloc(size_conv1 * sizeof(distType));
#endif

	for (int i = 0; i < size_conv1; i++)
		expandedZ[i] = 0;

	for (int i = 0; i < size_xyz; i++)
		expandedZ[i] = z[i];

	// Initialize C2 big enough for the second convolution
	int size_conv2 = 2 * size_conv1;

#ifdef __INTEL_COMPILER
	distType *C2 = (distType *)_mm_malloc(size_conv2 * sizeof(distType), 32);
#else
	distType *C2 = (distType *)malloc(size_conv2 * sizeof(distType));
#endif

	window_conv(C1, expandedZ, C2, h, size_conv1);

	// Calculate the probability integral over the bin for each point in the data
	for (int i = 0; i < dataSize; i++) {
		// rightmost boundry of integration for this particular bin
		int rightBound = (int)(data[i] / h);

		if (rightBound >= size_conv2)
			printf("ERROR: rightBound=%d>=size_conv2=%d ", rightBound, size_conv2);

		// the bin width in terms of indices
		int goback = (int)(0.1 / h);

		// the leftmost boundry of integration
		int leftBound = rightBound - goback;

		if (leftBound < 0)
			printf("ERROR: leftBound=%d<0 ", leftBound);

		// Calculate the right handed riemann sum
		Y[i] = 0;
		for (int j = leftBound + 1; j < rightBound; j++) {
			Y[i] += C2[j] * h;
		}

		//printf("Y[%d]=%f ", i, Y[i]);
	}

	*logP0 = (double)loglikelihood(Y, dataSize);

	printf("logP0=%f ", *logP0);

#ifdef __INTEL_COMPILER
	_mm_free(C1);
	_mm_free(C2);
	_mm_free(expandedZ);
#else
	free(C1);
	free(C2);
	free(expandedZ);
#endif

	return;
}
