

/* Include files */
#include "IMT_analysis_April2017.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sf_erf.h"

/* Function Definitions */

double emgpdf_loglikelihood(const gsl_vector *v, void *params)
{
	double *data = (double *)params;

	double l = gsl_vector_get(v, 0);
	double m = gsl_vector_get(v, 1);
	double s = gsl_vector_get(v, 2);

	double penalty = 0;
	if (m < 0 || s < 0 || l < 0)
		penalty = 1000;

	l = fabs(l);
	m = fabs(m);
	s = fabs(s);

	double Y[266];

	// emgpdf_replacement(data, l, m, s, Y);
	emgpdf(data, l, m, s, Y);
	

	double loglikelihood = 0;
	for (int i = 0; i < 266; i++) {
		loglikelihood += log(Y[i]);
	}

	return penalty - loglikelihood;
}

void emgpdf(const double X[266], double l, double m, double s, double Y[266])
{
	// Y=(l/2)*erfc((-X+m+l*s^2)/(s*2^(1/2))).*exp((l/2)*(-2*X+2*m+l*s^2));
	// https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	double a, b;
	for (int i = 0; i < 266; i++) {
		a = 0.5 * l * (2*m + l*pow(s, 2.0) - 2*X[i] );
		b = (m + l*pow(s, 2.0) - X[i]) / (s * pow(2, 0.5));
		Y[i] = 0.5 * l * exp(a) * gsl_sf_erfc(b);
	}
}

/* End of code generation (emgpdf.c) */
