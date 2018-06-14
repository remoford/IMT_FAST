



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

double optimize_twostage(int data_size, const double data[])
{

	printf("twostagefitnomle\n");


	static const signed char id[90] =
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
		2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6,
		6, 6, 6,
		7,
		7, 7, 8, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 2, 3, 4, 5, 6, 7, 8, 9,
		3, 4, 5,
		6,
		7, 8, 9, 4, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9, 6, 7, 8, 9, 7, 8, 9, 8,
		9, 9
	};


	

	/* prepare statistical variables */
	double mean = gsl_stats_mean(data, 1, data_size);
	double variance = gsl_stats_variance(data, 1, data_size);

	double vry[3] = { 0.25, 0.5, 0.75 };
	//double vrys[3] = { 0.01, 1, 10 };

	double m[3];
	double s[3];
	for (int i = 0; i < 3; i++) {
		m[i] = 1 / (mean * vry[i]);
		s[i] =
			pow((variance * vry[i]) / (pow(mean * vry[i], 3.0)), 0.5);
	}

	double pcomb[9][2];
	int pcomb_idx = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			pcomb[pcomb_idx][0] = m[i];
			pcomb[pcomb_idx][1] = s[j];
		}
	}

	/*  prepare parameter seeds */
	/* get all pairs of the form [m(i),s(j)] */
	/* these pairs represent all possible unique  */
	/* parameter choices for each part of the cell */
	/* cycle.   */
	/* pcomb = allcomb(m,s) */
	/* 'IMT_analysis_April2017:572' pcomb = [0.3397, 0.2359; 0.3397, 0.1180; 0.3397, 0.0786; 0.1699, 0.2359; 0.1699, 0.1180; 0.1699, 0.0786; 0.1132, 0.2359; 0.1132, 0.1180; 0.1132, 0.0786]; */
	/* place paramter pairs into a cell.  The parameters choices for each part */
	/* are now indexed */

	int itmp;
	cell_wrap_3 pcell[9];
	for (itmp = 0; itmp < 2; itmp++) {
		pcell[0].f1[itmp] = 0.3397 + -0.1038 * (double)itmp;
		pcell[1].f1[itmp] = 0.3397 + -0.2217 * (double)itmp;
		pcell[2].f1[itmp] = 0.3397 + -0.2611 * (double)itmp;
		pcell[3].f1[itmp] = 0.1699 + 0.066 * (double)itmp;
		pcell[4].f1[itmp] = 0.1699 + -0.0519 * (double)itmp;
		pcell[5].f1[itmp] =
			0.1699 + -0.091299999999999992 * (double)itmp;
		pcell[6].f1[itmp] = 0.1132 + 0.1227 * (double)itmp;
		pcell[7].f1[itmp] =
			0.1132 + 0.0047999999999999987 * (double)itmp;
		pcell[8].f1[itmp] =
			0.1132 + -0.034599999999999992 * (double)itmp;
	}

	/* get all pairs of indices for the parameter  */
	/* choices for each part of the cycle to get all  */
	/* parameter choices for the entire cycle */



	int emgfitnomle;
	double c_P[180];
	for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {
		/* 'IMT_analysis_April2017:596' P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}]; */
		for (itmp = 0; itmp < 2; itmp++) {
			c_P[emgfitnomle + 45 * itmp] =
				pcell[id[emgfitnomle] - 1].f1[itmp];
			c_P[emgfitnomle + 45 * (itmp + 2)] =
				pcell[id[45 + emgfitnomle] - 1].f1[itmp];
		}
	}

	double c_p[4];
#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (emgfitnomle = 0; emgfitnomle < 10; emgfitnomle++) {
		//for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {

		printf("i=%f\n", 1.0 + (double)emgfitnomle);

		printf("  P=[%f %f %f %f]\n", c_P[emgfitnomle],
			c_P[45 + emgfitnomle], c_P[90 + emgfitnomle],
			c_P[135 + emgfitnomle]);

		for (itmp = 0; itmp < 4; itmp++) {
			c_p[itmp] = c_P[emgfitnomle + 45 * itmp];
		}
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
		gsl_vector_set(x, 0, c_p[0]);
		gsl_vector_set(x, 1, c_p[1]);
		gsl_vector_set(x, 2, c_p[2]);
		gsl_vector_set(x, 3, c_p[3]);

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

			/*
			// Detect repeated states
			gsl_vector_sub(prev, s->x);
			if (gsl_vector_isnull(prev))
			repeated++;
			else
			repeated = 0;
			if (repeated > 10){
			printf("ERROR - QUITTING BECAUSE EXCESSIVE REPEATS!!!\n");
			status = GSL_SUCCESS;
			}
			gsl_vector_set(prev, 0, gsl_vector_get(s->x, 0));
			gsl_vector_set(prev, 1, gsl_vector_get(s->x, 1));
			gsl_vector_set(prev, 2, gsl_vector_get(s->x, 2));
			gsl_vector_set(prev, 3, gsl_vector_get(s->x, 3));
			*/

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
		c_p[0] = fabs(gsl_vector_get(s->x, 0));
		c_p[1] = fabs(gsl_vector_get(s->x, 1));
		c_p[2] = fabs(gsl_vector_get(s->x, 2));
		c_p[3] = fabs(gsl_vector_get(s->x, 3));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		printf("  p=[%f %f %f %f]\n", c_p[0], c_p[1], c_p[2], c_p[3]);

		double b_pd[180];
		for (itmp = 0; itmp < 4; itmp++) {
			b_pd[emgfitnomle + 45 * itmp] = c_p[itmp];
		}

		double l[266];
		conv2waldpdf(data, c_p[0], c_p[1], c_p[2], c_p[3], l, 0.01, 1,
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

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif

	double b_pd[180];
	double l[266];
	double b_flag[45];
	for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {

		conv2waldpdf(data, b_pd[emgfitnomle], b_pd[45 + emgfitnomle],
			b_pd[90 + emgfitnomle], b_pd[135 + emgfitnomle],
			l, 0.001, 1, data_size);

		double l_sum = 0;
		for (int i = 0; i < data_size; i++)
			l_sum += log(l[i]);

		b_flag[emgfitnomle] = l_sum;
	}

	double max_ld = DBL_MIN;
	int row_id = 0;
	for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {
		if (b_flag[emgfitnomle] > max_ld) {
			max_ld = b_flag[emgfitnomle];
			row_id = emgfitnomle + 1;
		}
	}

	printf("max_ld=%f row_id=%d\n", max_ld, row_id);
	printf("pd_max=[%f %f %f %f]\n\n", b_pd[row_id - 1],
		b_pd[45 + row_id - 1], b_pd[90 + row_id - 1],
		b_pd[135 + row_id - 1]);

	/*  common to each fit, consider factoring out */
	emgfitnomle = 1;
	double mtmp;
	mtmp = b_flag[0];
	itmp = 0;

	if (emgfitnomle < 45) {
		while (emgfitnomle + 1 < 46) {
			if (b_flag[emgfitnomle] > mtmp) {
				mtmp = b_flag[emgfitnomle];
				itmp = emgfitnomle;
			}

			emgfitnomle++;
		}
	}

	printf("max_ld=%f row_ld=%f\n", mtmp, (double)itmp + 1);
	printf("pd_max=[%f %f %f %f]\n\n", b_pd[itmp], b_pd[45 + itmp],
		b_pd[90 + itmp], b_pd[135 + itmp]);

	/*  END FUNCTION FIT_TWOSTAGE */
	



}