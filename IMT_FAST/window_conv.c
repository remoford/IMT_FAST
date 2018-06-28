#include "float.h"
#include "stdio.h"
#include "main.h"

#define _VERBOSE

void window_conv(const distType z[], const distType y[], distType C[], double h, long long int size_xyz)
{
    long long int size_conv = 2 * size_xyz;
	
	distType threshold = 0;

    for (long long int i = 0; i < size_conv; i++) {
		C[i] = 0.0;
    }

    /* Find the highest zero-valued index */
    long long int firstIdx = 0;
    for (long long int i = 0; i < size_xyz; i++) {
        if (z[i] > threshold)
			break;
        else
            firstIdx = i;

    }

	long long int lastIdx = size_xyz - 1;
	for (long long int i = size_xyz; i >= 0; i--) {
		if (z[i] > threshold)
			break;
		else
			lastIdx = i;
	}

	long long int firstYIdx = 0;
	for (long long int i = 0; i < size_xyz; i++) {
		if (y[i] > threshold)
			break;
		else
			firstYIdx = i;

	}

	long long int lastYIdx = size_xyz - 1;
	for (long long int i = size_xyz; i >= 0; i--) {
		if (y[i] > threshold)
			break;
		else
			lastYIdx = i;
	}

	long long int tripcount = 0;
    /* do the lopsided convolution */
    for (long long int i = firstIdx; i < lastIdx; i++) {
		for (long long int j = firstYIdx; j < lastYIdx; j++) {
			C[i + j] += z[i] * y[j] * h;
			tripcount++;
		}
    }

#ifdef _VERBOSE
	printf("[");

	/*
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
	*/

	printf("%f%%] ", 100 * ((double)(size_xyz*size_xyz) - (double)tripcount) / ((double)(size_xyz*size_xyz)));
#endif

}
