#include "imt_analysis.h"
#include "twostage.h"
#include "emgpdf.h"
#include "onestagepdf_lag.h"
#include "onestage.h"
#include "tailmass.h"
#include "gsl/gsl_multimin.h"
#include "math.h"
#include "binned_conv.h"
#include "nn_conv.h"
#include "main.h"
#include "loglikelihood.h"
#include "stdio.h"
#include "gsl/gsl_statistics_double.h"
#include "time.h"


//#define _CONV2INVG
#define _CONV2WALD
#define _VERBOSE
//#define _PARALLEL_PDF
#define _BINNED_MODE
#define _PARALLEL_SEEDS

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3
#endif

typedef struct {
	double f1[2];
} cell_wrap_3;


void optimize_twostage(int data_size, const double data[], int numseeds, double seeds[][4], configStruct config) {
	printf("twostagefitnomle\n");

	//numseeds = 5;

	double optimizedParams[45][4];

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {

		distType * seedll = (distType *)malloc(sizeof(distType)*data_size);
		conv2waldpdf(data, seeds[seedIdx][0],
			seeds[seedIdx][1], seeds[seedIdx][2],
			seeds[seedIdx][3], seedll, 0.01, 1,
			data_size);

		double seedll_sum = (double)loglikelihood(seedll, data_size);

		free(seedll);

		printf("\nstarting seedIdx=%d p=[%f %f %f %f] ll=%f\n", seedIdx, seeds[seedIdx][0],
			seeds[seedIdx][1], seeds[seedIdx][2],
			seeds[seedIdx][3], seedll_sum);

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
		x = gsl_vector_alloc(4);
		gsl_vector_set(x, 0, seeds[seedIdx][0]);
		gsl_vector_set(x, 1, seeds[seedIdx][1]);
		gsl_vector_set(x, 2, seeds[seedIdx][2]);
		gsl_vector_set(x, 3, seeds[seedIdx][3]);

		/*
		// Store previous state
		gsl_vector *prev;
		prev = gsl_vector_alloc(4);
		gsl_vector_set_all(prev, 0);
		*/

		//int repeated = 0;

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(4);
		gsl_vector_set_all(ss, 1.0);

		/* Initialize method and iterate */
		minex_func.n = 4;
		minex_func.f = convolv_2invG_adapt_nov_loglikelihood;
		minex_func.params = (void *) &config;

		s = gsl_multimin_fminimizer_alloc(T, 4);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

		do {
			iter++;

			clock_t t;
			t = clock();
			status = gsl_multimin_fminimizer_iterate(s);
			t = clock() - t;
			//printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);
			printf("%.3f ", ((float)t) / CLOCKS_PER_SEC);

			if (iter % 1 == 0)
				printf("\n\n");

			if (status)
				break;

			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
			if (status == GSL_SUCCESS) {
				printf("converged to minimum at\n");
			}

			/*
			printf("%5d %.17f %.17f %.17f %.17f f() = %.17f size = %.3f\n",
			(int)iter,
			gsl_vector_get(s->x, 0),
			gsl_vector_get(s->x, 1),
			gsl_vector_get(s->x, 2),
			gsl_vector_get(s->x, 3),
			s->fval, size);
			*/
#endif
		} while (status == GSL_CONTINUE && iter < 10000);
		printf("\n");

		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
		optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));
		optimizedParams[seedIdx][3] = fabs(gsl_vector_get(s->x, 3));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		distType * ll = (distType *)malloc(sizeof(distType)*data_size);
		conv2waldpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], ll, 0.01, 1,
			data_size);

		double l_sum = (double)loglikelihood(ll, data_size);

		free(ll);

		printf("\nfinished seedIdx=%d p=[%f %f %f %f] ll=%f\n", seedIdx, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], l_sum);
	}

	/*  we previously optimized with a larger step size, recalculate with */
	/*  a smaller stepsize after the fact */
	printf("\nrecalculating canidate solutions with smaller stepsize\n");

	distType * loglikelihoods = (distType *)malloc(sizeof(distType)*numseeds);
	distType * likelihoods = (distType *)malloc(sizeof(distType)*data_size);

	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		conv2waldpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2], optimizedParams[seedIdx][3],
			likelihoods, 0.001, 1, data_size);

		loglikelihoods[seedIdx] = loglikelihood(likelihoods, data_size);

		printf("id=%d p=[%f %f %f %f] ll=%f\n", seedIdx, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], loglikelihoods[seedIdx]);
	}
	free(likelihoods);

	// Find the best log likelihood
	double max_ld = 0.0;
	int row_id = 0;
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		if (loglikelihoods[seedIdx] > max_ld) {
			max_ld = (double)loglikelihoods[seedIdx];
			row_id = seedIdx;
		}
	}

	free(loglikelihoods);

	printf("\nBest fit: row_id=%d\n", row_id);
	printf("loglikelihood=%f ", max_ld);
	printf("p=[ %f %f %f %f ]\n", optimizedParams[row_id][0], optimizedParams[row_id][1], optimizedParams[row_id][2], optimizedParams[row_id][3]);

	return;
}

double
convolv_2invG_adapt_nov_loglikelihood(const gsl_vector * v, void * params)
{
    configStruct config = *(configStruct *)params;

	distType * data = config.data;
	int data_size = config.data_size;

    double m1 = gsl_vector_get(v, 0);
    double s1 = gsl_vector_get(v, 1);
    double m2 = gsl_vector_get(v, 2);
    double s2 = gsl_vector_get(v, 3);

    double penalty = 0;
    if (m1 < 0 || s1 < 0 || m2 < 0 || s2 < 0)
	penalty = 1000;

    m1 = fabs(m1);
    s1 = fabs(s1);
    m2 = fabs(m2);
    s2 = fabs(s2);

	distType * Y = (double *)malloc(sizeof(double)*data_size);

    conv2waldpdf(data, m1, s1, m2, s2, Y, 0.01, 1, data_size);

	double ll = (double) loglikelihood(Y, data_size);

	free(Y);

    return penalty - ll;
}

void
conv2waldpdf(const double X[], double m1, double s1, double m2, double s2,
	distType Y[], double h, int adaptiveMode, int size_XY) {

	int flag = 0;		// remember if we applied the approximation
	double eps = 0.01;		// a constant that is used in determining if the Dirac approximation should be applied.

	double E;
	if (adaptiveMode){
		E = 2000000000;
		//printf("Error E: %f\n", E);
	}
    else
		E = 0;

    double m_a = m1;
    double m_b = m2;
    double s_a = s1;
    double s_b = s2;

    double m[2];
    double s[2];
    m[0] = m1;
    m[1] = m2;
    s[0] = s1;
    s[1] = s2;

    // find the variance for both sub-distributions
    double v_a = (s_a * s_a) / (m_a * m_a * m_a);
    double v_b = (s_b * s_b) / (m_b * m_b * m_b);

    // find the standard deviation for each sub-distribution
    double sd_a = pow(v_a, 0.5);
    double sd_b = pow(v_b, 0.5);

    /* reorder m and s so that the sub-distribution with the smallest
       variance comes first.  So the fist part might be approximated as a Dirac delta. */
    if (sd_a > sd_b) {
	double tmp;
	tmp = m_a;
	m_a = m_b;
	m_b = tmp;
	tmp = s_a;
	s_a = s_b;
	s_b = tmp;
	tmp = v_a;
	v_a = v_b;
	v_b = tmp;
	tmp = sd_a;
	sd_a = sd_b;
	sd_b = tmp;
    }

    // find the largest point in t
    double maxX = 0;
    for (int i = 0; i < size_XY; i++) {
	if (X[i] > maxX)
	    maxX = X[i];
    }

    // These represent our loglikelihoods
    double logP0;
    double logP1;

    /* if the first pdf is very concentrated, check to see if it can be
       approximated as a point-mass distribution */
    if (sd_a < 0.01 && checktailmass(m_a, s_a, m_b, s_b, eps, sd_a, sd_b)) {
#ifdef _VERBOSE
	printf("using dirac delta approximation\n");
#endif

	/* If there is not much probability in the tails of the convolution,
	   we use a Dirac Delta for part 1.
	   Flag is set to 1 to indicate this. */
	double lag = 1 / m_a;

	waldlagpdf(X, m_b, s_b, lag, Y, size_XY);
	flag = 1;
    } else {

	/* produce a range of evenly spaced points to evaluate at between
	   0 and Maxt with step size h. the even spacing is important when
	   calculating the convolution later */
	int partitionLength = (int) (maxX / h);

	// This represents the partition
#ifdef __INTEL_COMPILER
	double *x =
	    (double *) _mm_malloc(partitionLength * sizeof(double), 32);
#else
	double *x = (double *) malloc(partitionLength * sizeof(double));
#endif

	// There are our two wald distributions
#ifdef __INTEL_COMPILER
	distType *y =
	    (distType *) _mm_malloc(partitionLength * sizeof(distType), 32);
#else
	distType *y = (distType *) malloc(partitionLength * sizeof(distType));
#endif

#ifdef __INTEL_COMPILER
	distType *z =
	    (distType *) _mm_malloc(partitionLength * sizeof(distType), 32);
#else
	distType *z = (distType *) malloc(partitionLength * sizeof(distType));
#endif

	double tally = 0;
	for (int i = 0; i < partitionLength; i++) {
	    x[i] = tally;
	    tally += h;
	}

	// evaluate the sub-distributions at each point in x
	waldpdf(x, m[0], s[0], y, partitionLength);
	waldpdf(x, m[1], s[1], z, partitionLength);

	//approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);

#ifdef _BINNED_MODE
	binned_conv(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);
#else
	nn_conv(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);
#endif

#ifdef __INTEL_COMPILER
	_mm_free(x);
	_mm_free(y);
	_mm_free(z);
#else
	free(x);
	free(y);
	free(z);
#endif

	while (E >= 0.001 * fabs(logP0)) {

	    h = h * 0.5;	// Shrink the step size

#ifdef _VERBOSE
	    printf("h=%f ", h);
#endif

	    int partitionLength = (int) (maxX / h);
#ifdef __INTEL_COMPILER
	    double *x =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
	    distType *y =
		(distType *) _mm_malloc(partitionLength * sizeof(distType),
				      32);
	    distType *z =
		(distType *) _mm_malloc(partitionLength * sizeof(distType),
				      32);
#else
	    double *x =
		(double *) malloc(partitionLength * sizeof(double));
	    distType *y =
		(distType *) malloc(partitionLength * sizeof(distType));
	    distType *z =
		(distType *) malloc(partitionLength * sizeof(distType));
#endif

	    // fill the partition
	    double tally = 0;
	    for (int i = 0; i < partitionLength; i++) {
			x[i] = tally;
			tally += h;
	    }

	    // evaluate the sub-distributions at each point in x
		waldpdf(x, m[0], s[0], y, partitionLength);
		waldpdf(x, m[1], s[1], z, partitionLength);

#ifdef _BINNED_MODE
		binned_conv(z, y, X, x, Y, &logP1, partitionLength, size_XY, h);
#else
		nn_conv(z, y, X, x, Y, &logP1, partitionLength, size_XY, h);
#endif

#ifdef __INTEL_COMPILER
	    _mm_free(x);
	    _mm_free(y);
	    _mm_free(z);
#else
	    free(x);
	    free(y);
	    free(z);
#endif

	    E = fabs(logP1 - logP0);

	    logP0 = logP1;

	}
    }
#ifdef _VERBOSE
	printf("\n");
#endif
}
