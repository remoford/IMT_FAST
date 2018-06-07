

#ifndef ONESTAGEPDF_LAG_H
#define ONESTAGEPDF_LAG_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"

/* Function Declarations */


extern double waldlag_loglikelihood(const gsl_vector *v, void *params);

extern void waldlagpdf(const double X[266], double mu, double s, double l, double Y[266], int size_XY);



#endif

