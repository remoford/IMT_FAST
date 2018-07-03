#include "imt_analysis.h"
#include "binned_conv.h"
#include "float.h"
#include "stdio.h"
#include "conv.h"
#include "window_conv.h"
#include "main.h"
#include "loglikelihood.h"
#include "utility.h"

//#define _VERBOSE
#define _WINDOW_MODE

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
	distType *C = (distType *)MALLOC(size_conv * sizeof(distType));

    for (int i = 0; i < size_conv; i++) {
		C[i] = 0;
    }

#ifdef _WINDOW_MODE
    window_conv(y,z,C,h,size_xyz);
#else
	conv(y, z, C, h, size_xyz);
#endif

	rightHandedRiemannSum((long)size_XY, data, h, (long)size_conv, 0.1, C, Y);

	*logP0 = (double) loglikelihood(Y, size_XY);

	FREE(C);
}

void threestage_binconv(const distType x[], const distType y[], const distType z[], const distType data[], distType Y[], double *logP0, int size_xyz, int dataSize, double h) {

	// Initialize C1 big enough for the first convolution
	int size_conv1 = 2 * size_xyz;

	distType *C1 = (distType *)MALLOC(size_conv1 * sizeof(distType));

	for (int i = 0; i < size_conv1; i++)
		C1[i] = 0;

	window_conv(x, y, C1, h, size_xyz);

	// Mask the new first convolved distribution beyond size_xyz
	// See error analysis section in paper for justification
	for (int i = size_xyz; i < size_conv1; i++)
		C1[i] = 0;

	distType *expandedZ = (distType *)MALLOC(size_conv1 * sizeof(distType));

	for (int i = 0; i < size_conv1; i++)
		expandedZ[i] = 0;

	for (int i = 0; i < size_xyz; i++)
		expandedZ[i] = z[i];

	// Initialize C2 big enough for the second convolution
	int size_conv2 = 2 * size_conv1;

	distType *C2 = (distType *)MALLOC(size_conv2 * sizeof(distType));

	window_conv(C1, expandedZ, C2, h, size_conv1);

	rightHandedRiemannSum((long)dataSize, data, h, (long)size_conv2, 0.1, C2, Y);

	*logP0 = (double)loglikelihood(Y, dataSize);

	FREE(C1);
	FREE(C2);
	FREE(expandedZ);

	return;
}
