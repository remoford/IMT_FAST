

#ifndef TAILMASS_H
#define TAILMASS_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"

/* Function Declarations */
extern double tailmass(const double m[2], const double s[2], double T2, const double sd[2]);

extern int checktailmass(const double m_a, const double s_a, const double m_b, const double s_b, double T2, const double sd_a, const double sd_b);

#endif

