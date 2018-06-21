

#ifndef ONESTAGEPDF2_H
#define ONESTAGEPDF2_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"
#include "main.h"

/* Function Declarations */

extern void optimize_onestage(const double data[], int data_size);

extern void waldpdf(const double X[], double mu, double s, distType Y[], int size_XY);

extern double wald_loglikelihood(const gsl_vector *v, void *params);


#endif

