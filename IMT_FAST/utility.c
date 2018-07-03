#include "float.h"
#include "main.h"
#include "stdio.h"

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

void rightHandedRiemannSum(long dataSize, distType * data, double gridSize, long partitionLength, double binSize, distType * C, distType * Y) {
	for (long i = 0; i < dataSize; i++) {

		// rightmost boundry of integration for this particular bin
		long rightBound = (long)(((double)data[i]) / gridSize);

		// the bin width in terms of indices
		long goback = (long)(binSize / gridSize);

		// the leftmost boundry of integration
		long leftBound = rightBound - goback;

		if (leftBound < 0)
			printf("ERROR: leftBound=%d<0 ", leftBound);

		if (rightBound > partitionLength)
			printf("ERROR: rightBound=%d>partitionLength=%d ", rightBound, partitionLength);

		//printf("leftBound=%ld ", leftBound);
		//printf("rightBound=%ld ", rightBound);

		// Calculate the right handed riemann sum
		Y[i] = 0;
		//printf("\n");
		for (long j = leftBound + 1; j <= rightBound; j++) {
			Y[i] += C[j] * gridSize;
			//printf("C[%d]=%.17f\n", j, C[j]);
		}

		//printf("Invg(%.17f)=Y[%d]=%.17f\n", data[i], i, Y[i]);
		//exit(1);
	}
}


distType * readfile(char * filename, int * arraySize)
{
	

	int readMax = 1000;

#ifdef __INTEL_COMPILER
	distType * readArray = (distType *)_mm_malloc(sizeof(distType) * readMax, 32);
#else
	distType * readArray = (distType *)malloc(sizeof(distType) * readMax);
#endif

	double currentNumber;

	int readCount = 0;

	FILE *fptr;
	fptr = fopen(filename, "r");

	if (fptr == NULL) {
		printf("Error: Unable to open file %s for reading!\n", filename);
	}
	else {
		while (fscanf(fptr, "%lf", &currentNumber) != EOF ) {
			readArray[readCount] = (distType) currentNumber;
			readCount++;
			if (readCount >= readMax) {
				printf("File had more than 1000 entries, truncating...\n");
				break;
			}
		}
	}

	printf("readCount = %d\n", readCount);

	fclose(fptr);

	distType * returnArray = NULL;

	if (readCount > 0) {
#ifdef __INTEL_COMPILER
		returnArray = (distType *) _mm_malloc(sizeof(distType)*readCount, 32);
#else
		returnArray = (distType *)malloc(sizeof(distType)*readCount);
#endif

		for (int i = 0; i < readCount; i++)
			returnArray[i] = (distType)readArray[i];
	}

#ifdef __INTEL_COMPILER
	_mm_free(readArray);
#else
	free(readArray);
#endif


	*arraySize = readCount;

	return returnArray;
}
