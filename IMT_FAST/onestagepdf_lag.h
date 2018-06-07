/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf_lag.h
 *
 * Code generation for function 'onestagepdf_lag'
 *
 */

#ifndef ONESTAGEPDF_LAG_H
#define ONESTAGEPDF_LAG_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"

/* Function Declarations */


extern double waldlag_loglikelihood(const gsl_vector *v, void *params);

extern void waldlagpdf(const double X[266], double mu, double s, double l, double Y[266], int size_XY);



#endif

/* End of code generation (onestagepdf_lag.h) */
