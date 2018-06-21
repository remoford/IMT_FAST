
 /* Include files */
#include "float.h"
#include "stdio.h"
#include "main.h"

#define _VERBOSE


void window_conv (const distType z[], const distType y[], distType C[], double h, int size_xyz)
{
    // Assume z is hte concentrated pdf
    //
    int size_conv = 2 * size_xyz;
	
	//double threshold = DBL_MIN;
	//double threshold = DBL_MIN * 2;
	distType threshold = 0;
	//distType threshold = distMin;

    for (int i = 0; i < size_conv; i++) {
		C[i] = 0.0;
    }


    /* Find the highest zero-valued index */
    int firstIdx = 0;
    for (int i = 0; i < size_xyz; i++) {
        if (z[i] > threshold)
			break;
        else
            firstIdx = i;

    }

	int lastIdx = size_xyz - 1;
	for (int i = size_xyz; i >= 0; i--) {
		if (z[i] > threshold)
			break;
		else
			lastIdx = i;
	}

	/* Find the highest zero-valued index */
	int firstYIdx = 0;
	for (int i = 0; i < size_xyz; i++) {
		if (y[i] > threshold)
			break;
		else
			firstYIdx = i;

	}

	int lastYIdx = size_xyz - 1;
	for (int i = size_xyz; i >= 0; i--) {
		if (y[i] > threshold)
			break;
		else
			lastYIdx = i;
	}

	int tripcount = 0;
    /* do the lopsided convolution */
    for (int i = firstIdx; i < lastIdx; i++) {
	for (int j = firstYIdx; j < lastYIdx; j++) {
	    C[i + j] += z[i] * y[j] * h;
		tripcount++;
	}
    }

#ifdef _VERBOSE
	printf("[");

	if (firstIdx != 0)
		printf("fz=%d ", firstIdx);
	if (lastIdx != size_xyz)
		printf("lz=%d ", lastIdx);
	if (firstYIdx != 0)
		printf("fy=%d ", firstYIdx);
	if (lastYIdx != size_xyz)
		printf("ly=%d ", lastYIdx);
	if (tripcount != size_xyz * size_xyz)
		printf("sz=%d ", size_xyz);

	printf("skp=%d] ", size_xyz*size_xyz - tripcount);
#endif

	/*
	for (int i = 0; i < size_conv; i++) {
		if (C[i] == 0)
			C[i] = distMin;

	}
	*/


}
