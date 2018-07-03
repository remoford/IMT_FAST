#include "main.h"

extern distType * readfile(char * filename, int * arraySize);

extern void rightHandedRiemannSum(long dataSize, distType * data, double gridSize, long partitionLength, double binSize, distType * C, distType * Y);