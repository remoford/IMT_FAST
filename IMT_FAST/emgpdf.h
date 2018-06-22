#ifndef EMGPDF_H
#define EMGPDF_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"
#include "main.h"

extern void optimize_emg(const double data[], int data_size, configStruct config);

extern double emgpdf_loglikelihood(const gsl_vector *v, void *params);

extern void emgpdf(const distType X[], double l, double m, double s, distType Y[], int data_size);

#endif