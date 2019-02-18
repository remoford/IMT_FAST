#ifndef GP_MAX_H
#define GP_MAX_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef _ENABLE_OLDTAILMASS
extern double gp_max_fixed(double m, double s);
#endif

#endif