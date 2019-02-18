#ifndef CONVOLV_3INVG_NOV_H
#define CONVOLV_3INVG_NOV_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP 
#include "omp.h"
#endif
#include "gsl/gsl_multimin.h"
#include "main.h"

extern void optimize_threestage(const distType data[], int data_size, configStruct config);

extern double convolv_3invG_nov_loglikelihood(const gsl_vector *v, void *params);

extern void threestage_adapt(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], int dataSize);

extern void threestage_bin(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], long dataSize, double gridSize);

extern void threestage_double_adapt_bin(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], long dataSize, double gridSize);

#endif