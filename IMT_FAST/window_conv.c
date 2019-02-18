#include "float.h"
#include "stdio.h"
#include "main.h"
#include "time.h"

#define _OLDCONV

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

    for (long i = 0; i < size_xyz; i++) {
		C[i] = 0;
    }

    /* Find the highest zero-valued index */
    long firstIdx = 0;
    for (unsigned long i = 0; i < size_xyz; i++) {
        if (z[i] > threshold)
			break;
        else
            firstIdx = i;
    }

	long lastIdx = size_xyz - 1;
	for (unsigned long i = size_xyz - 1; i > 0; i--) {
		if (z[i] > threshold)
			break;
		else
			lastIdx = i;
	}

	long firstYIdx = 0;
	for (unsigned long i = 0; i < size_xyz; i++) {
		if (y[i] > threshold)
			break;
		else
			firstYIdx = i;
	}

	long lastYIdx = size_xyz - 1;
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

	double estimateOpCount;
	double maxTripCount = opsPerIteration * (double)size_xyz * (double)size_xyz;

	if ((double)lastIdx - (double)firstIdx < 0 || (double)lastYIdx - (double)firstYIdx < 0) {
		printf("Skipping entire convolution as at least one pdf is all zeros! ");
	}
	else {


#ifdef _OLDCONV
		//unsigned long j;
		long j;
//#pragma omp parallel for schedule(static, 1) private(newLastYIdx, j)
		for (unsigned long i = firstIdx; i < lastIdx; i++) {
			if ((size_xyz - i) < lastYIdx) {
				newLastYIdx = size_xyz - i;
				if (firstExceed == 0) {
					firstExceed = i;	
				}
			}
			else
				newLastYIdx = lastYIdx;

			if (newLastYIdx - firstYIdx < 1 && lastComputedIdx == lastIdx)
				lastComputedIdx = i;

			//opCount += opsPerIteration * (double)(newLastYIdx - firstYIdx);
#ifdef _PARALLEL_CONV
#pragma omp parallel for private(j)
#endif
			for (j = firstYIdx; j < newLastYIdx; j++)
				C[i + j] += z[i] * y[j] * h;
		}
		estimateOpCount = opsPerIteration * (((double)(firstExceed - firstIdx) * (double)(lastYIdx - firstYIdx)) + (0.5 * (double)(lastYIdx - firstYIdx) * (double)(lastComputedIdx)));


#else
		//distType * CC = (distType *)MALLOC(sizeof(distType) * 2 * size_xyz);
		unsigned long int firstOutIdx = firstIdx + firstYIdx;
		unsigned long int lastOutIdx = lastIdx + lastYIdx;
		if (lastOutIdx >= size_xyz) {
			lastOutIdx = size_xyz;
			//printf("Clipping lastOutIdx\n");
		}
		if (firstOutIdx >= lastOutIdx)
			firstOutIdx = 0;
		distType acc;
		double skippedIndices = 0;
		unsigned long skippedIndicesInside;
		long j;
		long k;
		long firstKIdx = 0;
		long lastKIdx = 0;
#ifdef _PARALLEL_CONV
#pragma omp parallel for schedule(static, 16) private(acc, skippedIndicesInside, j, k, firstKIdx, lastKIdx)
#endif
		//for (unsigned long i = firstOutIdx; i < lastOutIdx; i++) {
		for (long i = 0 ; i < size_xyz ; i++){
			acc = 0;
			skippedIndicesInside = 0;



			/*
j < firstYIdx
i=169 k=161 j=8 firstKIdx=-2031 lastKIdx=-2031 firstIdx=0 lastIdx=2200 firstYIdx=12 lastYIdx=2200 !!!
i=202 k=198 j=4 firstKIdx=-1998 lastKIdx=-1998 firstIdx=0 lastIdx=2200 firstYIdx=12 lastYIdx=2200 !!!
j < firstYIdx
j < firstYIdx
i=148 k=142 j=6 firstKIdx=-2052 lastKIdx=-2052 firstIdx=0 lastIdx=2200 firstYIdx=12 lastYIdx=2200 !!!
i=202 k=199 j=3 firstKIdx=-1998 lastKIdx=-1998 firstIdx=0 lastIdx=2200 firstYIdx=12 lastYIdx=2200 !!!
j < firstYIdx
			*/



			firstKIdx = 0;
			lastKIdx = i;

			// Satisfy k >= firstIdx
			if(firstKIdx < firstIdx)
				firstKIdx = firstIdx;

			// Satisfy k <= lastIdx
			if (lastKIdx > lastIdx)
				lastKIdx = lastIdx;




			// Satisfy j <= lastYIdx

			long firstJIdx = i - firstKIdx;

			if (firstJIdx < i - firstYIdx)
				firstKIdx = firstYIdx - i;
			firstJIdx = i - firstKIdx;


			if (firstJIdx < i - lastYIdx)
				firstKIdx = lastYIdx - i;
			firstJIdx = i - firstKIdx;


			// Satisfy j >= firstYIdx

			long lastJIdx = i - lastKIdx;

			if (lastJIdx < i - firstYIdx) {

				printf("WOO");
				lastKIdx = i - firstYIdx;

			}
			lastJIdx = i - lastKIdx;





			// Careful about overlaps
			if (lastKIdx < firstKIdx)
				lastKIdx = firstKIdx;


			for (k = firstKIdx; k < lastKIdx; k++) {
			//for(k=firstKIdx ; k < lastKIdx ; k++){
				skippedIndicesInside++;
				//printf("i=%d k=%d j=%d firstKIdx=%d lastKIdx=%d firstIdx=%d lastIdx=%d firstYIdx=%d lastYIdx=%d !!!\n",
				//	i, k, j, firstKIdx, lastKIdx, firstIdx, lastIdx, firstYIdx, lastYIdx);

				j = i - k;
				//if (k >= firstIdx && k <= lastIdx && j >= firstYIdx && j <= lastYIdx) {
				if (j >= firstYIdx && j <= lastYIdx) {
					acc += z[k] * y[j];
				}
					
				
				else {
					
					printf("i=%ld k=%ld j=%ld firstKIdx=%ld lastKIdx=%ld firstJIdx=%ld lastJIdx=%ld firstIdx=%ld lastIdx=%ld firstYIdx=%ld lastYIdx=%ld !!!",
						i, k, j, firstKIdx, lastKIdx, firstJIdx, lastJIdx, firstIdx, lastIdx, firstYIdx, lastYIdx);

					if (k < firstIdx)
						printf("k < firstIdx ");
					if (k > lastIdx)
						printf("k > lastIdx ");
					if (j < firstYIdx)
						printf("j < firstYIdx ");
					if (j > lastYIdx)
						printf("j > lastYIdx ");
						
					printf("\n");
					/*
					if (k < firstKIdx || k > lastKIdx)
						printf("NEW BOUNDS NOT WORKING i=%d k=%d j=%d firstKIdx=%d lastKIdx=%d firstIdx=%d lastIdx=%d firstYIdx=%d lastYIdx=%d !!!\n",
							i, k, j, firstKIdx, lastKIdx, firstIdx, lastIdx, firstYIdx, lastYIdx);
					exit(1);
					*/

					skippedIndicesInside++;
					//	printf("Window skipping!\n");
				}
				
			}
			C[i] = acc * h;
			//skippedIndicesInside = i - skippedIndicesInside;
			//printf("skippedIndicesInside = %d", skippedIndicesInside);
			skippedIndices += skippedIndicesInside;
			
		}
		estimateOpCount = maxTripCount - (opsPerIteration * skippedIndices);
#endif
	}

	







	t = clock() - t;

	

	double skipPercentage = 100*((maxTripCount - estimateOpCount) / maxTripCount);

	double runtime = ((float)t) / CLOCKS_PER_SEC;

	double megaFlopsPerSecond = (estimateOpCount / runtime) / 1000000;

	printf("%lu size_xyz (%lu %lu , %lu %lu ) %g fullOps %g estOps %f%% skipped %f Mflop/s %fs]\n",
		size_xyz, firstIdx, lastIdx, firstYIdx, lastYIdx,
		maxTripCount, estimateOpCount, skipPercentage, megaFlopsPerSecond, runtime);
}
