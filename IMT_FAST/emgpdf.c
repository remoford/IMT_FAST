#include "imt_analysis.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sf_erf.h"
#include "main.h"
#include "loglikelihood.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics_double.h"


/*
Optimize parameters for the emg model using the data in the config struct
*/
void optimize_emg(configStruct config) {
	
	distType * data = config.data;
	int data_size = config.data_size;

	distType * l = (distType *)malloc(sizeof(distType)*data_size);

	printf("emgfitnomle\n\n");



	double * doubleData = (double *)malloc(sizeof(double)*data_size);
	for (int i = 0; i < data_size; i++)
		doubleData[i] = (double)data[i];


	/* prepare statistical variables */
	double mean = gsl_stats_mean(doubleData, 1, data_size);
	double variance = gsl_stats_variance(doubleData, 1, data_size);

	free(doubleData);


	double vry[3] = { 0.25, 0.5, 0.75 };

	/*we vary the parameters so the the Gaussian and exponential parts of
	the cell cycle are responsible for a fraction of the total mean and
	variance in the IMT.
	*/
	double lambda[3];
	for (int i = 0; i < 3; i++)
		lambda[i] = 1.0 / (mean * vry[i]);

	double mu[3];
	for (int i = 0; i < 3; i++)
		mu[i] = mean * vry[i];

	double sigma[3];
	for (int i = 0; i < 3; i++)
		sigma[i] = pow(variance * vry[i], 0.5);

	/* prepare parameter seeds */
	int seedIdx = 0;
	double paramSeeds[27][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				paramSeeds[seedIdx][0] = lambda[i];
				paramSeeds[seedIdx][1] = sigma[j];
				paramSeeds[seedIdx][2] = mu[k];
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

		printf("i = %d\n", 1 + seedIdx);
		printf("  P=[%.17f %.17f %.17f]\n", paramSeeds[seedIdx][0],
			paramSeeds[seedIdx][1], paramSeeds[seedIdx][2]);

		/* select our working seed */
		double paramSeed[3];
		for (int i = 0; i < 3; i++)
			paramSeed[i] = paramSeeds[seedIdx][i];

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
		gsl_vector_set(x, 0, paramSeed[0]);
		gsl_vector_set(x, 1, paramSeed[1]);
		gsl_vector_set(x, 2, paramSeed[2]);

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(3);
		gsl_vector_set_all(ss, 1.0);

		/* Initialize method and iterate */
		minex_func.n = 3;
		minex_func.f = emgpdf_loglikelihood;
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

			if (status == GSL_SUCCESS) {
				printf("converged to minimum at\n");
			}

			printf("%5d %.17f %.17f %.17f f() = %.17f size = %.3f\n",
				(int)iter,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				gsl_vector_get(s->x, 2), s->fval, size);

		} while (status == GSL_CONTINUE && iter < 10000);

		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
		optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		printf("  ep=[%.17f %.17f %.17f]\n",
			optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2]);

		/* recalculate our best fit */
		emgpdf(data, optimizedParams[seedIdx][0],
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
	int ind_le = 0;
	for (int i = 0; i < 27; i++) {
		if (le[i] < maxLikelihood) {
			maxLikelihood = le[i];
			ind_le = i;
		}
	}

	printf("max_le=%.17f ind_le=%d\n", maxLikelihood, ind_le + 1);
	printf("ep_max=[%.17f %.17f %.17f]\n\n",
		optimizedParams[ind_le][0], optimizedParams[ind_le][1],
		optimizedParams[ind_le][2]);
	
	free(l);
}

/*
Find and return the loglikelihood for an emg model with the parameters in v given the data in params
*/
double emgpdf_loglikelihood(const gsl_vector * v, void *params)
{
	configStruct config = *(configStruct *)params;

	distType * data = config.data;
	int data_size = config.data_size;

    double l = gsl_vector_get(v, 0);
    double m = gsl_vector_get(v, 1);
    double s = gsl_vector_get(v, 2);

    double penalty = 0;
    if (m < 0 || s < 0 || l < 0)
	penalty = 1000;

    l = fabs(l);
    m = fabs(m);
    s = fabs(s);

	distType * Y = (distType *)malloc(sizeof(distType)*data_size);

    emgpdf(data, l, m, s, Y, data_size);

	double ll = (double) loglikelihood(Y, data_size);

	free(Y);

    return penalty - ll;
}

/*
Evaluate the emg model with parameters l, m and s for each point in X[] of size data_size into Y[] by direct evaluation
*/
void emgpdf(const distType X[], double l, double m, double s, distType Y[], int data_size)
{
    // Y=(l/2)*erfc((-X+m+l*s^2)/(s*2^(1/2))).*exp((l/2)*(-2*X+2*m+l*s^2));
    // https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
    distType a, b;
    for (int i = 0; i < data_size; i++) {
		a = 0.5 * l * (2 * m + l * pow(s, 2.0) - 2 * X[i]);

		b = (m + l * pow(s, 2.0) - X[i]) / (s * pow(2, 0.5));

		Y[i] = 0.5 * l * exp(a) * gsl_sf_erfc(b);
    }
}