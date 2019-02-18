#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include "omp.h"
#endif
#include "main.h"

extern void binned_conv(const distType z[], const distType y[], const distType data[],
	const distType x[], distType Y[], double *logP0, int size_xyz, int size_XY, double h);

extern void threestage_binconv(const distType x[], const distType y[], const distType z[],
	const distType data[], distType Y[], double *logP0, int size_xyz, int dataSize, double h);