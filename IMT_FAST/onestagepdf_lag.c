/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf_lag.c
 *
 * Code generation for function 'onestagepdf_lag'
 *
 */

/* Include files */
#include "IMT_analysis_April2017.h"
#include "onestagepdf_lag.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "onestagepdf2.h"


/* Function Definitions */




double waldlag_loglikelihood(const gsl_vector *v, void *params)
{
	double *data = (double *)params;

	double m = gsl_vector_get(v, 0);
	double s = gsl_vector_get(v, 1);
	double l = gsl_vector_get(v, 2);

	double penalty = 0;
	if (m < 0 || s < 0 || l < 0)
		penalty = 1000;

	m = fabs(m);
	s = fabs(s);
	l = fabs(l);

	double Y[266];

	waldlagpdf(data, m, s, l, Y, 266);
	//onestagepdf_lag(data, m, s, l, Y);

	double loglikelihood = 0;
	for (int i = 0; i < 266; i++) {
		loglikelihood += log(Y[i]);
	}

	return penalty - loglikelihood;
}

void waldlagpdf(const double X[266], double mu, double s, double l, double Y[266], int size_XY)
{
	double a, b;
	for (int i = 0; i < size_XY; i++) {
		a = 1.0 / (s*pow(2 * M_PI  * pow(X[i] - l, 3.0), 0.5));
		b = (pow(mu*(X[i] - l) - 1, 2)) / (2.0 * s * s * (X[i] - l));
		Y[i] = a*exp(-b);
		if (Y[i] == 0)
			//Y[i] = 2.2250738585072014E-308;
			Y[i] = DBL_MIN;
	}
}


/* End of code generation (onestagepdf_lag.c) */
