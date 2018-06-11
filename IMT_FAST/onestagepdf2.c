
/* Include files */
#include "IMT_analysis_April2017.h"
#include "onestagepdf2.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "float.h"
#include "math.h"

//#define _VERBOSE

/* Function Definitions */

double wald_loglikelihood(const gsl_vector * v, void *params)
{
    double *data = (double *) params;

    double m = gsl_vector_get(v, 0);
    double s = gsl_vector_get(v, 1);

    double penalty = 0;
    if (m < 0 || s < 0)
	penalty = 1000;

    m = fabs(m);
    s = fabs(s);

    double Y[266];

    waldpdf(data, m, s, Y, 266);
    //onestagepdf2(data, m, s, Y);

    double loglikelihood = 0;
    for (int i = 0; i < 266; i++) {
	loglikelihood += log(Y[i]);
    }

    return penalty - loglikelihood;
}

void
waldpdf(const double X[], double mu, double s, double Y[], int size_XY)
{
    // Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t)));
    // https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution

    double a, b;
    for (int i = 0; i < size_XY; i++) {
	//a = 1.0 / (s*pow(2 * M_PI  * pow(X[i], 3.0), 0.5));
	a = 1.0 / (s * pow(6.2831853071795862 * pow(X[i], 3.0), 0.5));
	//b = (pow(mu*X[i] - 1, 2)) / (2.0 * s * s * X[i]);
	b = (pow(mu * X[i] - 1.0, 2.0)) / (2.0 * s * s * X[i]);
	Y[i] = a * exp(-b);
	if (Y[i] == 0)
	    Y[i] = 2.2250738585072014E-308;
	//Y[i] = DBL_MIN;

	if (isnan(Y[i]))
	    Y[i] = 2.2250738585072014E-308;
    }

#ifdef _VERBOSE
    printf("waldpdf = \n");
    for (int i = 0; i < size_XY; i++) {
	if (i % 8 == 0)
	    printf("\n");
	printf("%.17f ", Y[i]);
    }
    printf("\n\n");
#endif

}
