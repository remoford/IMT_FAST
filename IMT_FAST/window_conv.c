
 /* Include files */
#include "float.h"
#include "stdio.h"

//#define _VERBOSE

void window_conv (const double z[], const double y[], double C[], double h, int size_xyz)
{
    // Assume z is hte concentrated pdf
    //
    int size_conv = 2 * size_xyz;
	
#ifdef _VERBOSE
	printf("\nsize_xyz = %d ", size_xyz);
#endif

    for (int i = 0; i < size_conv; i++) {
	C[i] = 0;
    }


    /* Find the highest zero-valued index */
    int firstIdx = 0;
    for (int i = 0; i < size_xyz; i++) {
        if (z[i] > 0)
		break;
        else
            firstIdx = i;

    }
#ifdef _VERBOSE
	printf("First z idx = %d ", firstIdx);
#endif

	int lastIdx = size_xyz - 1;
	for (int i = size_xyz; i >= 0; i--) {
		if (z[i] > 0)
			break;
		else
			lastIdx = i;
	}
#ifdef _VERBOSE
	printf("Last z idx = %d ", lastIdx);
#endif

	/* Find the highest zero-valued index */
	int firstYIdx = 0;
	for (int i = 0; i < size_xyz; i++) {
		if (y[i] > 0)
			break;
		else
			firstYIdx = i;

	}
#ifdef _VERBOSE
	printf("First y idx = %d ", firstYIdx);
#endif

	int lastYIdx = size_xyz - 1;
	for (int i = size_xyz; i >= 0; i--) {
		if (z[i] > 0)
			break;
		else
			lastYIdx = i;
	}
#ifdef _VERBOSE
	printf("Last y idx = %d ", lastYIdx);
#endif

	int speedup = (lastIdx - firstIdx - size_xyz) + (lastYIdx - firstYIdx - size_xyz);

	if (speedup>1)
		printf(" Window Speedup=%d ", speedup);

#ifdef _VERBOSE
	printf("\n");
#endif

    /* do the lopsided convolution */
    for (int i = firstIdx; i < size_xyz; i++) {
	for (int j = 0; j < size_xyz; j++) {
	    C[i + j] += z[i] * y[j] * h;
	}
    }






}
