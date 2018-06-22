#ifndef ONESTAGEPDF_LAG_H
#define ONESTAGEPDF_LAG_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"
#include "main.h"

extern void optimize_onestagelag(const double data[], int data_size, configStruct config);

extern double waldlag_loglikelihood(const gsl_vector *v, void *params);

extern void waldlagpdf(const distType X[], double mu, double s, double l, distType Y[], int size_XY);

#endif