/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * approxconvolv.h
 *
 * Code generation for function 'approxconvolv'
 *
 */

#ifndef APPROXCONVOLV_H
#define APPROXCONVOLV_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"

/* Function Declarations */



extern void approxconvolv_replacement(const double z[], const double y[], const double X
	[], const double x[], double Y[], double *logP0, int size_xyz, int size_XY, double h);

#endif

/* End of code generation (approxconvolv.h) */
