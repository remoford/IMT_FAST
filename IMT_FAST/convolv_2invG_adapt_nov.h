

#ifndef CONVOLV_2INVG_ADAPT_NOV_H
#define CONVOLV_2INVG_ADAPT_NOV_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"

/* Function Declarations */
extern double convolv_2invG_adapt_nov_loglikelihood(const gsl_vector *v, void *params);

extern void conv2waldpdf(const double X[], double m1, double s1, double m2, double s2, double Y[], double h, int adaptiveMode, int size_XY);




#endif
