#include "float.h"
#include "stdio.h"
#include "main.h"

#define _VERBOSE

void window_conv(const distType z[], const distType y[], distType C[], double h, unsigned long size_xyz)
{
    unsigned long size_conv = 2 * size_xyz;
	
	distType threshold = 0;

    for (unsigned long i = 0; i < size_conv; i++) {
		C[i] = 0.0;
    }

    /* Find the highest zero-valued index */
    unsigned long firstIdx = 0;
    for (unsigned long i = 0; i < size_xyz; i++) {
        if (z[i] > threshold)
			break;
        else
            firstIdx = i;

    }

	unsigned long lastIdx = size_xyz - 1;
	for (unsigned long i = size_xyz; i >= 0; i--) {
		if (z[i] > threshold)
			break;
		else
			lastIdx = i;
	}

	unsigned long firstYIdx = 0;
	for (unsigned long i = 0; i < size_xyz; i++) {
		if (y[i] > threshold)
			break;
		else
			firstYIdx = i;

	}

	unsigned long lastYIdx = size_xyz - 1;
	for (unsigned long i = size_xyz; i >= 0; i--) {
		if (y[i] > threshold)
			break;
		else
			lastYIdx = i;
	}

	unsigned long tripcount = 0;
    /* do the lopsided convolution */
    for (unsigned long i = firstIdx; i < lastIdx; i++) {
		for (unsigned long j = firstYIdx; j < lastYIdx; j++) {
			C[i + j] += z[i] * y[j] * h;
			tripcount++;
		}
    }

#ifdef _VERBOSE
	printf("[");

	
	if (firstIdx != 0)
		printf("fz=%lu ", firstIdx);
	if (lastIdx != size_xyz)
		printf("lz=%lu ", lastIdx);
	if (firstYIdx != 0)
		printf("fy=%lu ", firstYIdx);
	if (lastYIdx != size_xyz)
		printf("ly=%lu ", lastYIdx);
	if (tripcount != size_xyz * size_xyz)
		printf("sz=%lu ", size_xyz);

	unsigned long maximumTripcount = size_xyz * size_xyz;

	unsigned long skipCount = maximumTripcount - tripcount;

	printf("skp=%lu ", skipCount);
	
	distType skipPercentage = ((distType)skipCount) / ((distType)maximumTripcount);

	if (skipPercentage < 0)
		printf("\n\n\nERROR: NEGATIVE SKIP PERCENTAGE\n\n\n");

	printf("tripcount=%lu maxtripcount=%lu %g.3%%] ", tripcount, maximumTripcount, skipPercentage);
#endif

}
