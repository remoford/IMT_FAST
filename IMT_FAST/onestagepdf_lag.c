#include "imt_analysis.h"
#include "onestagepdf_lag.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "onestage.h"
#include "main.h"
#include "loglikelihood.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics_double.h"

/*
Optimize parameters for a one stage lag model using the data in the config struct
*/
void optimize_onestagelag(configStruct config) {

	distType * data = config.data;
	int data_size = config.data_size;

	distType * l = (distType *)malloc(sizeof(distType)*data_size);

	printf("onestagefitlagnomle\n\n");

	/* prepare statistical variables */
	double * doubleData = (double *)malloc(sizeof(double)*data_size);
	for (int i = 0; i < data_size; i++)
		doubleData[i] = (double)data[i];
	double mean = gsl_stats_mean(doubleData, 1, data_size);
	double variance = gsl_stats_variance(doubleData, 1, data_size);
	free(doubleData);

	// C3 = sum((data-C1).^3)/(length(data));
	double c3sum = 0;
	for (int i = 0; i < data_size; i++) {
		c3sum += pow(data[i] - mean, 3.0);
	}
	c3sum = c3sum / data_size;

	double mu = c3sum / (3 * pow(variance, 2.0));
	double sigma =
		pow((pow(c3sum, 3.0) / (27 * pow(variance, 5.0))), 0.5);
	double lag = mean - 3 * pow(variance, 2.0) / c3sum;

	double vryv[3] = { 0.5, 1, 2 };
	double vrym[3] = { 0.25, 0.5, 0.75 };

	double m[3];
	for (int i = 0; i < 3; i++)
		m[i] = mu * vrym[i];

	double s[3];
	for (int i = 0; i < 3; i++)
		s[i] = sigma * vryv[i];

	double L[3];
	for (int i = 0; i < 3; i++)
		L[i] = lag * vrym[i];

	/* prepare parameter seeds */
	int seedIdx = 0;
	double paramSeeds[64][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				paramSeeds[seedIdx][0] = m[i];
				paramSeeds[seedIdx][1] = s[j];
				paramSeeds[seedIdx][2] = L[k];
				seedIdx++;
			}
		}
	}

	/* optimize parameters */
	double optimizedParams[27][3];
	double le[27];

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (seedIdx = 0; seedIdx < 27; seedIdx++) {

		printf("i=%d\n", 1 + seedIdx);
		printf("  P=[%f %f %f]\n", paramSeeds[seedIdx][0],
			paramSeeds[seedIdx][1], paramSeeds[seedIdx][2]);

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
		x = gsl_vector_alloc(3);
		gsl_vector_set(x, 0, paramSeeds[seedIdx][0]);
		gsl_vector_set(x, 1, paramSeeds[seedIdx][1]);
		gsl_vector_set(x, 2, paramSeeds[seedIdx][2]);

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(3);
		gsl_vector_set_all(ss, 1.0);

		/* Initialize method and iterate */
		minex_func.n = 3;
		minex_func.f = waldlag_loglikelihood;
		minex_func.params = (void *) &config;

		s = gsl_multimin_fminimizer_alloc(T, 3);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

		do {
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);

			if (status)
				break;

			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
			if (status == GSL_SUCCESS) {
				printf("converged to minimum at\n");
			}

			printf("%5d %.17f %.17f %.17f f() = %.17f size = %.3f\n",
				(int)iter,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				gsl_vector_get(s->x, 2), s->fval, size);
#endif
		} while (status == GSL_CONTINUE && iter < 10000);

		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
		optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		printf("  p=[%f %f %f]\n", optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2]);

		waldlagpdf(data, optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2], l, data_size);

		/* calculate the log likelihood for our best fit */
		/*
		double loglikelihood = 0;
		for (int i = 0; i < data_size; i++) {
		loglikelihood += log(l[i]);
		}
		*/

		double ll = (double)loglikelihood(l, data_size);

		printf("  l=%.17f\n\n", ll);
		le[seedIdx] = ll;
	}

	/* find the best optimized parameter set for all starting seeds tried */
	double maxLikelihood = 0;
	int ind_ld = 0;
	for (int i = 0; i < 27; i++) {
		if (le[i] < maxLikelihood) {
			maxLikelihood = le[i];
			ind_ld = i;
		}
	}

	printf("max_ld=%f row_ld=%d\n", maxLikelihood, ind_ld + 1);
	printf("pd_max=[%f %f %f]\n\n", optimizedParams[ind_ld][0],
		optimizedParams[ind_ld][1], optimizedParams[ind_ld][2]);
	
	free(l);
}

double waldlag_loglikelihood(const gsl_vector * v, void *params)
{
	configStruct config = *(configStruct *)params;

	distType * data = config.data;
	int data_size = config.data_size;

    double m = gsl_vector_get(v, 0);
    double s = gsl_vector_get(v, 1);
    double l = gsl_vector_get(v, 2);

    double penalty = 0;
    if (m < 0 || s < 0 || l < 0)
	penalty = 1000;

    m = fabs(m);
    s = fabs(s);
    l = fabs(l);

	distType * Y = (distType *)malloc(sizeof(distType)*data_size);

    waldlagpdf(data, m, s, l, Y, data_size);

	double ll = (double) loglikelihood(Y, data_size);

	free(Y);

    return penalty - ll;
}

void waldlagpdf(const distType X[], double mu, double s, double l,
	   distType Y[], int size_XY)
{
    double a, b;
    for (int i = 0; i < size_XY; i++) {
		a = 1.0 / (s * pow(2 * M_PI * pow(X[i] - l, 3.0), 0.5));
		b = (pow(mu * (X[i] - l) - 1, 2)) / (2.0 * s * s * (X[i] - l));
		Y[i] = a * exp(-b);
    }
}
