#include "float.h"
#include "main.h"
#include "stdio.h"

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

/*
Compare two distType numbers. This is used together with qsort elsewhere
*/
int compare(const void * a, const void * b) {

	distType A = *(distType*)a;
	distType B = *(distType*)b;

	return (A > B) - (A < B);
}

/*
Calculate the euclidian distance between point a and point b in the three stage parameter space

For use with two or one stage, set the unused coordinates to zero.
*/
double euclideanDistance(double a_m1, double a_s1, double a_m2, double a_s2, double a_m3, double a_s3, double b_m1, double b_s1, double b_m2, double b_s2, double b_m3, double b_s3) {

	double delta_m1 = a_m1 - b_m1;
	double delta_s1 = a_s1 - b_s1;
	double delta_m2 = a_m2 - b_m2;
	double delta_s2 = a_s2 - b_s2;
	double delta_m3 = a_m3 - b_m3;
	double delta_s3 = a_s3 - b_s3;

	printf("delta_m1=%f delta_s1=%f delta_m2=%f delta_s2=%f delta_m3=%f delta_s3=%f\n", delta_m1, delta_s1, delta_m2, delta_s3, delta_m3, delta_s3);

	double distance = sqrt(delta_m1*delta_m1 + delta_s1*delta_s1 + delta_m2*delta_m2 + delta_s2*delta_s2 + delta_m3*delta_m3 + delta_s3*delta_s3);

	printf("calculated distance = %.17f\n", distance);
	return distance;
}

/*
Calculate the integrated probability for the interval cooresponding to the sampling bin containing the
point on the provided probability density grid C of length partitionLength with spacing gridSize for
each data point in data[] of size dataSize and store that probability into Y[]

Despite the name, this is actually doing trapezoidal integration!
*/
void rightHandedRiemannSum(long dataSize, const distType data[], double gridSize, long partitionLength, double binSize, distType * C, distType * Y) {

	long skippedBins = 0;
	for (long i = 0; i < dataSize; i++) {

		// Don't waste time calculating things twice 
		if (i > 0) {
			if (data[i] == data[i - 1]) {
				Y[i] = Y[i - 1];
				skippedBins++;
				continue;
			}
		}

		// rightmost boundry of integration for this particular bin
		long rightBound = (long)(((double)data[i]) / gridSize);

		// the bin width in terms of indices
		long goback = (long)(binSize / gridSize);

		// the leftmost boundry of integration
		long leftBound = rightBound - goback;

		if (leftBound < 0)
			printf("ERROR(rightHandedRiemannSum): leftBound=%ld<0 ", leftBound);

		if (rightBound > partitionLength)
			printf("ERROR(rightHandedRiemannSum): rightBound=%ld>partitionLength=%ld ", rightBound, partitionLength);

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

	//printf("Skipped %ld bins %f%%\n", skippedBins, 100 * (double)skippedBins / (double)dataSize);

}

/*
Read a file in filename which is a list of numbers, one per line, into an array of size arraySize and
return the pointer to that array
*/
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


double checkNormalization(distType C[], long data_size, double gridSize)
{
	double totalProbability = 0;

	printf("data_size = %ld gridSize = %f ", data_size, gridSize);

	for (long i = 1; i < data_size; i++) {

		// Rectangular integration
		//totalProbability += C[i] * gridSize;

		// Trapezodial integration
		totalProbability += 0.5* gridSize * (C[i - 1] + C[i]);

		//printf("C[%d]=%.17f\n", j, C[j]);
	}

	printf("totalProbability = %f\n", totalProbability);

	return totalProbability;
}


void beginTraceFun(char * functionName) {

#ifdef _ENABLE_FUNCTION_TRACING

	for (int i = 0; i < traceDepth; i++)
		printf("\t");

	printf("BEGIN %s\n", functionName);


	traceDepth++;

#endif

}

void endTraceFun(char * functionName) {

#ifdef _ENABLE_FUNCTION_TRACING

	traceDepth--;

	for (int i = 0; i < traceDepth; i++)
		printf("\t");

	printf("END %s\n", functionName);

#endif

}