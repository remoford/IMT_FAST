/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf2.h
 *
 * Code generation for function 'onestagepdf2'
 *
 */

#ifndef ONESTAGEPDF2_H
#define ONESTAGEPDF2_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "gsl/gsl_multimin.h"

/* Function Declarations */



extern void waldpdf(const double X[], double mu, double s, double Y[], int size_XY);

extern double wald_loglikelihood(const gsl_vector *v, void *params);

#ifdef _OLD_MATLAB_CODE
extern void onestagepdf2(const double t[266], double mu, double s, double Y[266]);
#endif

#endif

/* End of code generation (onestagepdf2.h) */
