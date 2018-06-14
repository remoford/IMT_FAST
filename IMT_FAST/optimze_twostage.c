



#include "stdio.h"
#include "convolv_2invG_adapt_nov.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics_double.h"
#include "time.h"

/* Type Definitions */
#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3

typedef struct {
	double f1[2];
} cell_wrap_3;

#endif				/*typedef_cell_wrap_3 */

void optimize_twostage(int data_size, const double data[])
{
	printf("twostagefitnomle\n");

	int numseeds = 45;
	double seeds[45][4] = {
		{ 0.339700, 0.235900, 0.339700, 0.235900 } ,
		{ 0.339700, 0.235900, 0.339700, 0.118000 } ,
		{ 0.339700, 0.235900, 0.339700, 0.078600 } ,
		{ 0.339700, 0.235900, 0.169900, 0.235900 } ,
		{ 0.339700, 0.235900, 0.169900, 0.118000 } ,
		{ 0.339700, 0.235900, 0.169900, 0.078600 } ,
		{ 0.339700, 0.235900, 0.113200, 0.235900 } ,
		{ 0.339700, 0.235900, 0.113200, 0.118000 } ,
		{ 0.339700, 0.235900, 0.113200, 0.078600 } ,
		{ 0.339700, 0.118000, 0.339700, 0.118000 } ,
		{ 0.339700, 0.118000, 0.339700, 0.078600 } ,
		{ 0.339700, 0.118000, 0.169900, 0.235900 } ,
		{ 0.339700, 0.118000, 0.169900, 0.118000 } ,
		{ 0.339700, 0.118000, 0.169900, 0.078600 } ,
		{ 0.339700, 0.118000, 0.113200, 0.235900 } ,
		{ 0.339700, 0.118000, 0.113200, 0.118000 } ,
		{ 0.339700, 0.118000, 0.113200, 0.078600 } ,
		{ 0.339700, 0.078600, 0.339700, 0.078600 } ,
		{ 0.339700, 0.078600, 0.169900, 0.235900 } ,
		{ 0.339700, 0.078600, 0.169900, 0.118000 } ,
		{ 0.339700, 0.078600, 0.169900, 0.078600 } ,
		{ 0.339700, 0.078600, 0.113200, 0.235900 } ,
		{ 0.339700, 0.078600, 0.113200, 0.118000 } ,
		{ 0.339700, 0.078600, 0.113200, 0.078600 } ,
		{ 0.169900, 0.235900, 0.169900, 0.235900 } ,
		{ 0.169900, 0.235900, 0.169900, 0.118000 } ,
		{ 0.169900, 0.235900, 0.169900, 0.078600 } ,
		{ 0.169900, 0.235900, 0.113200, 0.235900 } ,
		{ 0.169900, 0.235900, 0.113200, 0.118000 } ,
		{ 0.169900, 0.235900, 0.113200, 0.078600 } ,
		{ 0.169900, 0.118000, 0.169900, 0.118000 } ,
		{ 0.169900, 0.118000, 0.169900, 0.078600 } ,
		{ 0.169900, 0.118000, 0.113200, 0.235900 } ,
		{ 0.169900, 0.118000, 0.113200, 0.118000 } ,
		{ 0.169900, 0.118000, 0.113200, 0.078600 } ,
		{ 0.169900, 0.078600, 0.169900, 0.078600 } ,
		{ 0.169900, 0.078600, 0.113200, 0.235900 } ,
		{ 0.169900, 0.078600, 0.113200, 0.118000 } ,
		{ 0.169900, 0.078600, 0.113200, 0.078600 } ,
		{ 0.113200, 0.235900, 0.113200, 0.235900 } ,
		{ 0.113200, 0.235900, 0.113200, 0.118000 } ,
		{ 0.113200, 0.235900, 0.113200, 0.078600 } ,
		{ 0.113200, 0.118000, 0.113200, 0.118000 } ,
		{ 0.113200, 0.118000, 0.113200, 0.078600 } ,
		{ 0.113200, 0.078600, 0.113200, 0.078600 }
	};

	double optimizedParams[45][4];

	/* prepare statistical variables */
	double mean = gsl_stats_mean(data, 1, data_size);
	double variance = gsl_stats_variance(data, 1, data_size);

	double vry[3] = { 0.25, 0.5, 0.75 };

	double m[3];
	double s[3];
	for (int i = 0; i < 3; i++) {
		m[i] = 1 / (mean * vry[i]);
		s[i] =
			pow((variance * vry[i]) / (pow(mean * vry[i], 3.0)), 0.5);
	}

	printf("\n");
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		printf("i=%d ", seedIdx);
		printf("  P=[%f %f %f %f]\n", seeds[seedIdx][0],
			seeds[seedIdx][1], seeds[seedIdx][2],
			seeds[seedIdx][3]);
	}
	printf("\n");

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++){

		printf("i=%d ", seedIdx);
		printf("  P=[%f %f %f %f]\n", seeds[seedIdx][0],
			seeds[seedIdx][1], seeds[seedIdx][2],
			seeds[seedIdx][3]);
		
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

		// Store previous state
		gsl_vector *prev;
		prev = gsl_vector_alloc(4);
		gsl_vector_set_all(prev, 0);
		int repeated = 0;

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(4);
		gsl_vector_set_all(ss, 1.0);

		/* Initialize method and iterate */
		minex_func.n = 4;
		minex_func.f = convolv_2invG_adapt_nov_loglikelihood;
		minex_func.params = (void *)data;

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

		printf("  p=[%f %f %f %f]\n", optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], optimizedParams[seedIdx][3]);

		double l[266];
		conv2waldpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], l, 0.01, 1,
			data_size);

		double l_sum = 0;
		for (int i = 0; i < data_size; i++) {
			l_sum += log(l[i]);
		}
		/* FIXME FIXME FIXME we are not actually reporting hp, flag and E here!*/
		printf("  l=%f hp=%f flag=%f E=%f WARNING HP FLAG AND E ARE BOGUS NUMBERS\n\n", l_sum, 9999.0, 9999.0, 9999.0);
	}

	/*  we previously optimized with a larger step size, recalculate with */
	/*  a smaller stepsize after the fact */
	printf("recalculating canidate solutions with smaller stepsize\n");
	double loglikelihoods[45];
	double l[266];
#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (int seedIdx = 0; seedIdx < 45; seedIdx++) {
		conv2waldpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][0], optimizedParams[seedIdx][0],
			l, 0.001, 1, data_size);
		double l_sum = 0;
		for (int i = 0; i < data_size; i++)
			l_sum += log(l[i]);

		loglikelihoods[seedIdx] = l_sum;
	}

	// Find the best log likelihood
	double max_ld = DBL_MIN;
	int row_id = 0;
	for (int seedIdx = 0; seedIdx < 45; seedIdx++) {
		if (loglikelihoods[seedIdx] > max_ld) {
			max_ld = loglikelihoods[seedIdx];
			row_id = seedIdx;
		}
	}

	printf("Best fit: row_id=%d\n", row_id + 1);
	printf("loglikelihood=%f ", max_ld);
	printf("p=[ %f %f %f %f ]\n", optimizedParams[row_id][0], optimizedParams[row_id][1], optimizedParams[row_id][2], optimizedParams[row_id][3]);

	return;
}