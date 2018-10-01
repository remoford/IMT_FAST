#include "float.h"
#include "stdio.h"
#include "main.h"
#include "time.h"

/*
Compute the windowed algebraic convolution of vectors y and z of size size_xyz into C with a grid size of h

This skips inconsequential operations by picking a window for each vector and only convolve those windows. As
the runtime of convolutions looks like O(n*m), even small reductions in either n or m can give large speedups.

For a threshold of zero, inconsequential means that the result will be idential regardless of the use of
windowing. However, for nonzero thresholds error is introduced. This is difficult to characterize.
*/
void window_conv(const distType z[], const distType y[], distType C[], double h, unsigned long size_xyz)
{
	clock_t t;
	t = clock();

    //unsigned long size_conv = 2 * size_xyz;

	printf("[sz=%luKB ", ((unsigned long)sizeof(distType)*size_xyz) / 1024);

	if (size_xyz > 65535)
		printf(" WARNING ulong OVERFLOW!!! ");
	
	distType threshold = 0;

    for (unsigned long i = 0; i < size_xyz; i++) {
		C[i] = 0;
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
	for (unsigned long i = size_xyz - 1; i > 0; i--) {
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
	for (unsigned long i = size_xyz - 1; i > 0; i--) {
		if (y[i] > threshold)
			break;
		else
			lastYIdx = i;
	}

	double opsPerIteration = 3;
	unsigned long newLastYIdx;
    /* do the lopsided convolution */
	//double opCount = 0;

	unsigned long firstExceed = 0;
	unsigned long lastComputedIdx = lastIdx;

	if ((double)lastIdx - (double)firstIdx < 0 || (double)lastYIdx - (double)firstYIdx < 0) {
		printf("Skipping entire convolution as at least one pdf is all zeros! ");
	}
	else {

//#pragma omp parallel for
		for (unsigned long i = firstIdx; i < lastIdx; i++) {
			if ((size_xyz - i) < lastYIdx) {
				newLastYIdx = size_xyz - i;
				if (firstExceed == 0) {
					firstExceed = i;	
				}
			}
			else
				newLastYIdx = lastYIdx;

			if (newLastYIdx - firstYIdx < 1)
				lastComputedIdx = i;

			//opCount += opsPerIteration * (double)(newLastYIdx - firstYIdx);

#pragma omp parallel for
			for (unsigned long j = firstYIdx; j < newLastYIdx; j++)
				C[i + j] += z[i] * y[j] * h;
		}
	}

	double estimateOpCount = opsPerIteration * (((double)(firstExceed - firstIdx) * (double)(lastYIdx - firstYIdx)) + (0.5 * (double)(lastYIdx - firstYIdx) * (double)(lastComputedIdx)));

	t = clock() - t;

	double maxTripCount = opsPerIteration * (double)size_xyz * (double)size_xyz;

	double skipPercentage = 100*((maxTripCount - estimateOpCount) / maxTripCount);

	double runtime = ((float)t) / CLOCKS_PER_SEC;

	double megaFlopsPerSecond = (estimateOpCount / runtime) / 1000000;

	printf("%lu size_xyz (%lu %lu , %lu %lu ) %g fullOps %g estOps %f%% skipped %f Mflop/s %fs]\n",
		size_xyz, firstIdx, lastIdx, firstYIdx, lastYIdx,
		maxTripCount, estimateOpCount, skipPercentage, megaFlopsPerSecond, runtime);
}
