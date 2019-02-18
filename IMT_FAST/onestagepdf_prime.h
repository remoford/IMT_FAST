#ifndef ONESTAGEPDF_PRIME_H
#define ONESTAGEPDF_PRIME_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP 
#include "omp.h"
#endif

extern void onestagepdf_prime_fixed(const double t[], int t_size, double m, double s, double Y[]);

#endif