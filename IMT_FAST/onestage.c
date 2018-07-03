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
#include "stdio.h"
#include "time.h"
#include "utility.h"

#define _VERBOSE

void optimize_onestage(const distType data[], long data_size, configStruct config) {
	printf("onestagefitnomle\n\n");

	distType * l = (distType *)MALLOC(sizeof(distType)*data_size);

	double ld[16];

	/* prepare statistical variables */
	double * doubleData = (double *)malloc(sizeof(double)*data_size);
	for (int i = 0; i < data_size; i++) 
		doubleData[i] = (double)data[i];

	double mean = gsl_stats_mean(doubleData, 1, data_size);
	double variance = gsl_stats_variance(doubleData, 1, data_size);

	free(doubleData);

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

		printf("starting seedIdx=%d ", 1 + seedIdx);


		printf("P=[%f %f]\n", paramSeeds[seedIdx][0], paramSeeds[seedIdx][1]);

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
		printf("\n");

		do {
			iter++;

			clock_t t;
			t = clock();

			printf("iter=%d\n", (int)iter);
			status = gsl_multimin_fminimizer_iterate(s);

			t = clock() - t;

			if (status)
				break;

			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, 1e-4);

			if (status == GSL_SUCCESS) {
				printf("converged to minimum at\n");
			}

			printf("ll=%.17e [%.17f %.17f] size=%.3f %.3fs\n\n",
				s->fval,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				size,
				((float)t) / CLOCKS_PER_SEC);

		} while (status == GSL_CONTINUE && iter < 10000);

		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		printf("  p=[%f %f]\n", optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][1]);

		//onestagepdf2(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], l);
		wald_adapt(data, optimizedParams[seedIdx][0],
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

	FREE(l);
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

	distType * Y = (distType *)MALLOC(sizeof(distType)*data_size);

    //waldpdf(data, m, s, Y, data_size);
	wald_adapt(data, m, s, Y, data_size);
 
	double ll = (double) loglikelihood(Y, data_size);

	FREE(Y);

    return penalty - ll;
}

void wald_adapt(const distType data[], double mu, double s, distType Y[], long data_size) {

#ifdef _VERBOSE
	printf("[mu=%f s=%f\n", mu, s);
#endif

	distType gridSize = 0.01;
	distType E = 2000000000;
	distType ll_previous;
	distType ll_current;

	wald_bin(data, mu, s, Y, data_size, gridSize);

	ll_previous = loglikelihood(Y, data_size);

#ifdef _VERBOSE
	printf("  gridSize=%.17f E=%.17e eThr=%.17e ll=%.17e\n", gridSize, E, _ERROR_BOUND, ll_previous);
#endif

	while (E >= _ERROR_BOUND) {
		gridSize = gridSize * 0.5;

		wald_bin(data, mu, s, Y, data_size, gridSize);

		ll_current = loglikelihood(Y, data_size);

		E = fabs(ll_current - ll_previous);

#ifdef _VERBOSE
		printf("  gridSize=%.17f E=%.17e eThr=%.17e ll=%.17e\n", gridSize, E, _ERROR_BOUND, ll_current);
#endif

		ll_previous = ll_current;
	}

#ifdef _VERBOSE
	printf("]\n");
#endif
	return;
}

void wald_bin(const distType data[], double mu, double s, distType Y[], long dataSize, double gridSize) {
	
	//printf("dataSize=%d\n", dataSize);

	double binSize = 0.1;

	double maxData = 0;
	for (long i = 0; i < dataSize; i++) {
		if (data[i] > maxData)
			maxData = data[i];
	}

	long partitionLength = (long)(maxData / gridSize) + 1;

	distType * partition = (distType *)MALLOC(sizeof(distType)*partitionLength);

	for (long i = 0; i < partitionLength; i++) {
		partition[i] = i * gridSize;
	}

	distType * C = (distType *)MALLOC(sizeof(distType)*partitionLength);

	waldpdf(partition, mu, s, C, partitionLength);

	rightHandedRiemannSum((long)dataSize, data, gridSize, (long)partitionLength, binSize, C, Y);

	FREE(C);

}

void
waldpdf(const distType data[], double mu, double s, distType Y[], long dataSize)
{
    // Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t)));
    // https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution


	if (mu <= 0.0 || s <= 0.0)
		printf("ERROR: waldpdf:nonpositive mu or s!\n");

	distType a, b, bExp;

	//FILE *f;
	//f = fopen("onestage.log", "a");

	for (long i = 0; i < dataSize; i++) {

		errno = 0;

		int anyError = 0;

		a = 1.0 / (s * pow(6.2831853071795862 * pow(data[i], 3.0), 0.5));

		if (errno == ERANGE) {
			printf("ERROR: waldpdf:(a) range error! ");
			printf("data[%d]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}
		if (errno == EDOM) {
			printf("ERROR: waldpdf:(a) domain error! ");
			printf("data[%d]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}

		errno = 0;
		b = -((pow(mu * data[i] - 1.0, 2.0)) / (2.0 * s * s * data[i]));

		if (errno == EDOM) {
			printf("ERROR: waldpdf:(b) domain error! ");
			printf("data[%d]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}
		if (errno == ERANGE) {
			printf("ERROR: waldpdf:(b) domain error! ");
			printf("data[%d]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}

#ifdef _DIST_SINGLE
		if (b >= FLT_MAX_EXP)
#else
		if (b >= DBL_MAX_EXP)
#endif
			printf("ERROR: waldpdf:Predict exp(%f) overflow for type double!\n", b);

		errno = 0;
		bExp = exp(b);

		if (errno == ERANGE) {
			printf("ERROR: waldpdf:(bExp) range error! ");
			printf("data[%d]=%g mu=%g s=%g b=%g ", i, data[i], mu, s, b);
			anyError = 1;
		}
		if (errno == EDOM) {
			printf("ERROR: waldpdf:(bExp) domain error! ");
			printf("data[%d]=%g mu=%g s=%g b=%g ", i, data[i], mu, s, b);
			anyError = 1;
		}

		Y[i] = (distType)(a * bExp);


		if (anyError != 0)
			printf(" Y[%d]=%f\n", i, Y[i]);


		if (!(isfinite(a) && isfinite(b))) {
			/* Ok this fixup requires some explaination. Strictly wald(0) is indeterminate.
			We replace it here with zero as that is the practical value and this avoids
			blowing up a log likelihood calculation later! */
			//printf("{a=%f b=%f data[i]=%f Y[i]=%f} ", a, b, data[i], (double)Y[i]);
			Y[i] = 0;
		}

		if (!isfinite(Y[i])) {
			printf("ERROR: waldpdf:InvG(%f)=Y[%d]=inf mu=%f s=%f HOLY NORMALIZATION BATMAN\n", data[i], i, mu, s);
			Y[i] = 0;
		}

		/*
		if (Y[i] > 1) {
			printf("ERROR: waldpdf:InvG(%f)=Y[%d]>1=%f mu=%f s=%f HOLY NORMALIZATION BATMAN\n", data[i], i, Y[i], mu, s);
			Y[i] = 0;
		}
		*/

		if (Y[i] < 0) {
			printf("ERROR: waldpdf:InvG(%f)=Y[%d]<0=%f mu=%f s=%f HOLY NORMALIZATION BATMAN\n", data[i], i, Y[i], mu, s);
			Y[i] = 0;
		}

		//if (i >= dataSize - 10)
		//	printf("InvG(%f)=Y[%d]=%f\n", data[i], i, Y[i]);


		//fprintf(f, "InvG(%f)=Y[%d]=%f\n", data[i], i, Y[i]);
		
    }

	//fclose(f);
}
