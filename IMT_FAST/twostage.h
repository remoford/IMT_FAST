#ifndef CONVOLV_2INVG_ADAPT_NOV_H
#define CONVOLV_2INVG_ADAPT_NOV_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"
#include "main.h"

extern void optimize_twostage(int data_size, const double data[], int numseeds, double seeds[][4], configStruct config);

extern double convolv_2invG_adapt_nov_loglikelihood(const gsl_vector *v, void *params);

extern void conv2waldpdf(const distType data[], double m1, double s1, double m2, double s2, distType convolvedPDF[], double h, int adaptiveMode, int size_XY);

#endif