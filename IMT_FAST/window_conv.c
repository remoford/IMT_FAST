
 /* Include files */
#include "float.h"
#include "stdio.h"

//#define _VERBOSE

void window_conv (const double z[], const double y[], double C[], double h, int size_xyz)
{
    // Assume z is hte concentrated pdf
    //
    int size_conv = 2 * size_xyz;
	
	//double threshold = DBL_MIN;
	//double threshold = DBL_MIN * 2;
	double threshold = 0;

#ifdef _VERBOSE
	printf("\nsize_xyz = %d ", size_xyz);
#endif

    for (int i = 0; i < size_conv; i++) {
	C[i] = 0;
    }


    /* Find the highest zero-valued index */
    int firstIdx = 0;
    for (int i = 0; i < size_xyz; i++) {
        if (z[i] > threshold)
			break;
        else
            firstIdx = i;

    }
#ifdef _VERBOSE
	printf("First z idx = %d ", firstIdx);
#endif

	int lastIdx = size_xyz - 1;
	for (int i = size_xyz; i >= 0; i--) {
		if (z[i] > threshold)
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
		if (y[i] > threshold)
			break;
		else
			firstYIdx = i;

	}
#ifdef _VERBOSE
	printf("First y idx = %d ", firstYIdx);
#endif

	int lastYIdx = size_xyz - 1;
	for (int i = size_xyz; i >= 0; i--) {
		if (y[i] > threshold)
			break;
		else
			lastYIdx = i;
	}
#ifdef _VERBOSE
	printf("Last y idx = %d ", lastYIdx);


	int zspeedup = lastIdx - firstIdx - size_xyz;
	int yspeedup = lastYIdx - firstYIdx - size_xyz;

	//if (speedup>1)
		printf(" [%d %d]", zspeedup, yspeedup);

		if (abs(zspeedup) + abs(yspeedup) > 2)
			printf("!!!!!!!!!!! ");
		else
			printf(" ");
#endif

#ifdef _VERBOSE
	//printf("\n");
#endif
	int tripcount = 0;
    /* do the lopsided convolution */
    for (int i = firstIdx; i < lastIdx; i++) {
	for (int j = firstYIdx; j < lastYIdx; j++) {
	    C[i + j] += z[i] * y[j] * h;
		tripcount++;
	}
    }

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

	for (int i = 0; i < size_conv; i++) {
		if (C[i] == 0)
			C[i] == DBL_MIN;

	}


}
