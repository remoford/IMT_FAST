#ifndef MAIN_H
#define MAIN_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

extern int main(int argc, const char * const argv[]);

//#define _DIST_SINGLE

#ifdef _DIST_SINGLE
typedef float distType;
#define distMin FLT_MIN
#else
typedef double distType;
#define distMin DBL_MIN
#endif


typedef struct  {
	distType * data;
	int data_size;

} configStruct;

#endif