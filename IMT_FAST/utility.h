#include "main.h"

extern distType * readfile(char * filename, int * arraySize);

double euclideanDistance(double a_m1, double a_s1, double a_m2, double a_s2, double a_m3, double a_s3, double b_m1, double b_s1, double b_m2, double b_s2, double b_m3, double b_s3);

extern void rightHandedRiemannSum(long dataSize, const distType * data, double gridSize, long partitionLength, double binSize, distType * C, distType * Y);

extern int compare(const void * a, const void * b);

extern double checkNormalization(distType C[], long data_size, double gridSize);

extern void beginTraceFun(char * functionName);

extern void endTraceFun(char * functionName);
