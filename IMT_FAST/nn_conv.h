


/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "main.h"

/* Function Declarations */



extern void nn_conv(const distType z[], const distType y[], const double X
	[], const double x[], distType Y[], double *logP0, int size_xyz, int size_XY, double h);


