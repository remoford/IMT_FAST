

#ifndef EMGPDF_H
#define EMGPDF_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"

/* Function Declarations */
#ifdef _OLD_MATLAB_CODE
extern void emgpdf_old(const double X[266], double l, double m, double s, double Y[266]);
#endif

extern double emgpdf_loglikelihood(const gsl_vector *v, void *params);

extern void emgpdf(const double X[266], double l, double m, double s, double Y[266]);


#endif

