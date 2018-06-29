#include "float.h"
#include "stdio.h"
#include "main.h"

#define _VERBOSE

void window_conv(const distType z[], const distType y[], distType C[], double h, long size_xyz)
{
    long size_conv = 2 * size_xyz;
	
	distType threshold = 0;

    for (long i = 0; i < size_conv; i++) {
		C[i] = 0.0;
    }

    /* Find the highest zero-valued index */
    long firstIdx = 0;
    for (long i = 0; i < size_xyz; i++) {
        if (z[i] > threshold)
			break;
        else
            firstIdx = i;

    }

	long lastIdx = size_xyz - 1;
	for (long i = size_xyz; i >= 0; i--) {
		if (z[i] > threshold)
			break;
		else
			lastIdx = i;
	}

	long firstYIdx = 0;
	for (long i = 0; i < size_xyz; i++) {
		if (y[i] > threshold)
			break;
		else
			firstYIdx = i;

	}

	long lastYIdx = size_xyz - 1;
	for (long i = size_xyz; i >= 0; i--) {
		if (y[i] > threshold)
			break;
		else
			lastYIdx = i;
	}

	long tripcount = 0;
    /* do the lopsided convolution */
    for (long i = firstIdx; i < lastIdx; i++) {
		for (long j = firstYIdx; j < lastYIdx; j++) {
			C[i + j] += z[i] * y[j] * h;
			tripcount++;
		}
    }

#ifdef _VERBOSE
	printf("[");

	
	if (firstIdx != 0)
		printf("fz=%lld ", firstIdx);
	if (lastIdx != size_xyz)
		printf("lz=%lld ", lastIdx);
	if (firstYIdx != 0)
		printf("fy=%lld ", firstYIdx);
	if (lastYIdx != size_xyz)
		printf("ly=%lld ", lastYIdx);
	if (tripcount != size_xyz * size_xyz)
		printf("sz=%lld ", size_xyz);

	printf("skp=%lld ", size_xyz*size_xyz - tripcount);
	

	printf("%g%%] ", 100 * ((double)(size_xyz*size_xyz) - (double)tripcount) / ((double)(size_xyz*size_xyz)));
#endif

}
