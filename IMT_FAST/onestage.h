#ifndef ONESTAGEPDF2_H
#define ONESTAGEPDF2_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"
#include "main.h"

extern void optimize_onestage(const distType data[], long data_size, configStruct config);

extern void wald_adapt(const distType data[], double mu, double s, distType Y[], long dataSize);

extern void wald_bin(const distType data[], double mu, double s, distType Y[], long dataSize, double gridSize);

extern void waldpdf(const distType data[], double mu, double s, distType Y[], long dataSize);

extern double wald_loglikelihood(const gsl_vector *v, void *params);

#endif