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

/*
Evaluating the algebraic convolution of vectors y and z of size size_xyz,
for each point in data[] of size size_XY,
compute the integrated probability of the sampling bin containing the data point,
into Y[] using the given grid x with spacing h and find the compensated
loglikelihood logP0

FIXME do we even use the grid x anywhere in this function anyways?
*/
void binned_conv(const distType z[], const distType y[],
			  const distType data[], const distType x[], distType Y[],
			  double *logP0, int size_xyz, int size_XY,
			  double h)
{

	beginTraceFun("binned_conv");

	// Sanity checking, make sure all the data points are positive! Negative stopping times are bad.
    int allZero = 1;
    for (int i = 0; i < size_XY; i++) {
	if (data[i] > 0)
	    allZero = 0;
    }
    if (allZero)
	printf("data[] is not all positive!\n");

	// C stores the result of the convolution
	distType *C = (distType *)MALLOC(size_xyz * sizeof(distType));

	// Initialize C with zeros
    for (int i = 0; i < size_xyz; i++) {
		C[i] = 0;
    }

#ifdef _WINDOW_MODE
    window_conv(y,z,C,h,size_xyz);
#else
	conv(y, z, C, h, size_xyz);
#endif

	rightHandedRiemannSum((long)size_XY, data, h, (long)size_xyz, 0.1, C, Y);

	for (long i = 0; i < size_XY; i++) {
		if (Y[i] < 0)
			printf("ERROR: binned_conv(): L[%d] = %f < 0\n", i, Y[i]);
		if (Y[i] > 1)
			printf("ERROR: binned_conv(): L[%d] = %f > 1\n", i, Y[i]);
	}

	*logP0 = (double) loglikelihood(Y, size_XY);

	FREE(C);


	endTraceFun("binned_conv");
}

/*
Evaluating the algebraic convolution of vectors x, y and z of size size_xyz,
for each point in data[] of dataSize,
compute the integrated probability of the sampling bin containing the data point,
into Y[] using the grid spacing h and find the compensated
loglikelihood logP0
*/
void threestage_binconv(const distType x[], const distType y[], const distType z[], const distType data[], distType Y[], double *logP0, int size_xyz, int dataSize, double h) {

	distType *C1 = (distType *)MALLOC(size_xyz * sizeof(distType));

	for (int i = 0; i < size_xyz; i++)
		C1[i] = 0;

	window_conv(x, y, C1, h, size_xyz);

	distType *C2 = (distType *)MALLOC(size_xyz * sizeof(distType));

	for (int i = 0; i < size_xyz; i++)
		C2[i] = 0;

	window_conv(C1, z, C2, h, size_xyz);

	rightHandedRiemannSum((long)dataSize, data, h, (long)size_xyz, 0.1, C2, Y);

	*logP0 = (double)loglikelihood(Y, dataSize);

	FREE(C1);
	FREE(C2);

	return;
}
