

#ifndef CONVOLV_3INVG_NOV_H
#define CONVOLV_3INVG_NOV_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"


/* Function Declarations */


extern double convolv_3invG_nov_loglikelihood(const gsl_vector *v, void *params);

extern void convolv3waldpdf(double m1, double s1, double m2, double s2, double m3, double s3, const double X[266], double Y[266], int size_XY, double h);


#endif

