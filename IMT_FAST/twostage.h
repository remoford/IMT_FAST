#ifndef CONVOLV_2INVG_ADAPT_NOV_H
#define CONVOLV_2INVG_ADAPT_NOV_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "gsl/gsl_multimin.h"
#include "main.h"

enum adaptConvergenceMethod { LOGLIKELIHOOD, NORMALIZATION };

extern double ** twostage_seeds(double mean, double variance, int *numSeeds);

extern void optimize_twostage(int numseeds, double ** seeds, configStruct config);

extern double convolv_2invG_adapt_nov_loglikelihood(const gsl_vector *v, void *params);

extern void conv2waldpdf(const distType data[], double m1, double s1, double m2, double s2, distType convolvedPDF[], double h, int adaptiveMode, int size_XY, enum adaptConvergenceMethod convergenceMethod);

extern void twostage_bin(const distType data[], double m1, double s1, double m2, double s2, distType Y[], long dataSize, double gridSize, enum adaptConvergenceMethod convergenceMethod);

#endif