#include "imt_analysis.h"
#include "onestage.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "float.h"
#include "math.h"
#include "main.h"
#include "loglikelihood.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics_double.h"

//#define _VERBOSE

void optimize_onestage(const double data[], int data_size, configStruct config) {
	printf("onestagefitnomle\n\n");

	distType * l = (distType *)malloc(sizeof(distType)*data_size);

	double ld[16];

	/* prepare statistical variables */
	double mean = gsl_stats_mean(data, 1, data_size);
	double variance = gsl_stats_variance(data, 1, data_size);
	double vry[4] = { 0.5, 1.0, 1.5, 2.0 };

	double mu = 1.0 / mean;
	double sigma = pow((variance / pow(mean, 3)), 0.5);

	double m[4];
	for (int i = 0; i < 4; i++)
		m[i] = mu * vry[i];

	double s[4];
	for (int i = 0; i < 4; i++)
		s[i] = sigma * vry[i];

	/* prepare parameter seeds */
	int seedIdx = 0;
	double paramSeeds[16][2];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			paramSeeds[seedIdx][0] = m[i];
			paramSeeds[seedIdx][1] = s[j];
			seedIdx++;
		}
	}

	/* optimize parameters */
	double optimizedParams[27][3];
	//double le[27];

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (seedIdx = 0; seedIdx < 16; seedIdx++) {
		printf("i=%d\n", 1 + seedIdx);
		printf("  P=[%f %f]\n", paramSeeds[seedIdx][0],
			paramSeeds[seedIdx][1]);

		// https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives
		const gsl_multimin_fminimizer_type *T =
			gsl_multimin_fminimizer_nmsimplex2;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		gsl_multimin_function minex_func;

		size_t iter = 0;
		int status;
		double size;

		/* Starting point */
		x = gsl_vector_alloc(2);
		gsl_vector_set(x, 0, paramSeeds[seedIdx][0]);
		gsl_vector_set(x, 1, paramSeeds[seedIdx][1]);

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(2);
		gsl_vector_set_all(ss, 1.0);

		/* Initialize method and iterate */
		minex_func.n = 2;
		minex_func.f = wald_loglikelihood;
		minex_func.params = (void *) &config;

		s = gsl_multimin_fminimizer_alloc(T, 2);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

		do {
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);

			if (status)
				break;

			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, 1e-4);

			if (status == GSL_SUCCESS) {
				printf("converged to minimum at\n");
			}

			printf("%5d %.17f %.17f f() = %.17f size = %.3f\n",
				(int)iter,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1), s->fval, size);

		} while (status == GSL_CONTINUE && iter < 10000);

		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		printf("  p=[%f %f]\n", optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][1]);

		//onestagepdf2(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], l);
		waldpdf(data, optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][1], l, data_size);

		double ll = (double)loglikelihood(l, data_size);

		printf("  l=%.17f\n\n", ll);
		ld[seedIdx] = ll;
	}

	/* find the best optimized parameter set for all starting seeds tried */
	double maxLikelihood = 0;
	int ind_ld = 0;
	for (int i = 0; i < 16; i++) {
		if (ld[i] < maxLikelihood) {
			maxLikelihood = ld[i];
			ind_ld = i;
		}
	}

	printf("max_ld=%f row_ld=%d\n", maxLikelihood, ind_ld);
	printf("pd_max=[%f %f]\n\n", optimizedParams[ind_ld][0],
		optimizedParams[ind_ld][1]);

	free(l);
}

double wald_loglikelihood(const gsl_vector * v, void *params)
{
	configStruct config = *(configStruct *)params;

	distType * data = config.data;
	int data_size = config.data_size;

    double m = gsl_vector_get(v, 0);
    double s = gsl_vector_get(v, 1);

    double penalty = 0;
    if (m < 0 || s < 0)
	penalty = 1000;

    m = fabs(m);
    s = fabs(s);

	distType * Y = (distType *)malloc(sizeof(distType)*data_size);

    waldpdf(data, m, s, Y, data_size);
 
	double ll = (double) loglikelihood(Y, data_size);

	free(Y);

    return penalty - ll;
}

void
waldpdf(const double X[], double mu, double s, distType Y[], int size_XY)
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

		if (isnan(Y[i]))
			Y[i] = 0;
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
