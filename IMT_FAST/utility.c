#include "float.h"
#include "main.h"
#include "stdio.h"

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

distType compare(const void * a, const void * b) {
	return (*(distType*)a - *(distType*)b);
}


void rightHandedRiemannSum(long dataSize, const distType data[], double gridSize, long partitionLength, double binSize, distType * C, distType * Y) {
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

			// Rectangular integration
			//Y[i] += C[j] * gridSize;

			// Trapezodial integration
			Y[i] += 0.5* gridSize * (C[j - 1] + C[j]);

			//printf("C[%d]=%.17f\n", j, C[j]);
		}

		//printf("Invg(%.17f)=Y[%d]=%.17f\n", data[i], i, Y[i]);
		//exit(1);
	}
}


distType * readfile(char * filename, int * arraySize)
{
	int readMax = 1000;

	distType * readArray = (distType *)MALLOC(sizeof(distType) * readMax);

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
		returnArray = (distType *)MALLOC(sizeof(distType)*readCount);

		for (int i = 0; i < readCount; i++)
			returnArray[i] = (distType)readArray[i];
	}

	FREE(readArray);

	*arraySize = readCount;

	return returnArray;
}
