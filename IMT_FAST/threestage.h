#ifndef CONVOLV_3INVG_NOV_H
#define CONVOLV_3INVG_NOV_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"
#include "main.h"

extern void optimize_threestage(const distType data[], int data_size, configStruct config);

extern double convolv_3invG_nov_loglikelihood(const gsl_vector *v, void *params);

#ifdef _OLD_THREESTAGE
extern void convolv3waldpdf(double m1, double s1, double m2, double s2, double m3, double s3, const double X[266], distType Y[266], int size_XY, double h);
#endif

extern void threestage_adapt(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], int dataSize);

extern void threestage_bin(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], long dataSize, double gridSize);

#endif