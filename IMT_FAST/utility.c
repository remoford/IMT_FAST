#include "float.h"
#include "main.h"
#include "stdio.h"

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

distType * readfile(char * filename, int * arraySize)
{
	

	int readMax = 1000;

#ifdef __INTEL_COMPILER
	distType * readArray = (distType *)_mm_malloc(sizeof(distType) * readMax, 32);
#else
	distType * readArray = (distType *)malloc(sizeof(distType) * readMax);
#endif

	float currentNumber;

	int readCount = 0;

	FILE *fptr;
	fptr = fopen(filename, "r");

	if (fptr == NULL) {
		printf("Error: Unable to open file %s for reading!\n", filename);
	}
	else {
		while (fscanf(fptr, "%f", &currentNumber) != EOF ) {
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
