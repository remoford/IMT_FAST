#include "float.h"
#include "main.h"

#include "stdio.h"


distType * readfile(char * filename, int * arraySize)
{
	

	int readMax = 1000;

	float * readArray = (float *) malloc(sizeof(float) * readMax);

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
		returnArray = (distType *) malloc(sizeof(distType)*readCount);

		for (int i = 0; i < readCount; i++)
			returnArray[i] = (distType)readArray[i];
	}

	free(readArray);

	*arraySize = readCount;

	return returnArray;
}
