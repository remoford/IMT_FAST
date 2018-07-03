#include "float.h"
#include "stdio.h"
#include "main.h"
#include "time.h"

#define _VERBOSE

void window_conv(const distType z[], const distType y[], distType C[], double h, unsigned long size_xyz)
{
	clock_t t;
	t = clock();

    unsigned long size_conv = 2 * size_xyz;

	printf("[sz=%luKB ", (sizeof(distType)*size_conv) / 1024);

	//if (size_xyz >= 65536)
	//	printf("ERROR: convolution steps too big to count with an unsigned long int!!!!\n");
	
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

	unsigned long newLastYIdx;
    /* do the lopsided convolution */
    for (unsigned long i = firstIdx; i < lastIdx; i++) {
		if ((size_xyz - i) < lastYIdx)
			newLastYIdx = size_xyz - i;
		else
			newLastYIdx = lastYIdx;

		for (unsigned long j = firstYIdx; j < newLastYIdx; j++)
			C[i + j] += z[i] * y[j] * h;
    }

	t = clock() - t;

	printf("%fs]\n", ((float)t) / CLOCKS_PER_SEC);
}
