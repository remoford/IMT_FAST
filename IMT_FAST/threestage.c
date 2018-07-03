#include "imt_analysis.h"
#include "threestage.h"
#include "onestage.h"
#include "twostage.h"
#include "gp_max.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sort.h"
#include "binned_conv.h"
#include "main.h"
#include "loglikelihood.h"
#include "time.h"
#include "binned_conv.h"

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3

typedef struct {
	double f1[2];
} cell_wrap_3;

#endif

#define _GSL_GP_MAX_FIXED
#define _BINNED_MODE

void optimize_threestage(const distType data[], int data_size, configStruct config) {

	int itmp;
	cell_wrap_3 b_pcell[4];
	static const double dv61[2] = { 0.8494, 0.0843 };
	int emgfitnomle;
	double d_P[120];
	static const signed char b_id[60] =
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 1, 1, 1, 1, 2, 2, 2,
		3, 3, 4, 2, 2, 2, 3, 3, 4, 3, 3, 4, 4, 1, 2, 3, 4, 2, 3, 4, 3, 4, 4, 2, 3, 4,
		3, 4, 4, 3, 4, 4, 4
	};
	double c_pd[120];
	double d_p[6];

	distType * l = (distType *)MALLOC(sizeof(distType)*data_size);

	double c_flag[20];
	double mtmp;

	printf("threestagefitnomle\n");

	for (itmp = 0; itmp < 2; itmp++) {
		b_pcell[0].f1[itmp] =
			0.8494 + -0.25960000000000005 * (double)itmp;
		b_pcell[1].f1[itmp] = dv61[itmp];
		b_pcell[2].f1[itmp] =
			0.1213 + 0.46849999999999997 * (double)itmp;
		b_pcell[3].f1[itmp] =
			0.1213 + -0.037000000000000005 * (double)itmp;
	}

	for (emgfitnomle = 0; emgfitnomle < 20; emgfitnomle++) {
		for (itmp = 0; itmp < 2; itmp++) {
			d_P[emgfitnomle + 20 * itmp] =
				b_pcell[b_id[emgfitnomle] - 1].f1[itmp];
			d_P[emgfitnomle + 20 * (itmp + 2)] =
				b_pcell[b_id[20 + emgfitnomle] - 1].f1[itmp];
			d_P[emgfitnomle + 20 * (itmp + 4)] =
				b_pcell[b_id[40 + emgfitnomle] - 1].f1[itmp];
		}
	}

	memset(&c_pd[0], 0, 120U * sizeof(double));

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (emgfitnomle = 0; emgfitnomle < 17; emgfitnomle++) {

		printf("i=%f ", 1.0 + (double)emgfitnomle);

		printf("P=[%f %f %f %f %f %f]\n", d_P[emgfitnomle],
			d_P[20 + emgfitnomle], d_P[40 + emgfitnomle],
			d_P[60 + emgfitnomle], d_P[80 + emgfitnomle],
			d_P[100 + emgfitnomle]);

		for (itmp = 0; itmp < 6; itmp++) {
			d_p[itmp] = d_P[emgfitnomle + 20 * itmp];
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
		x = gsl_vector_alloc(6);
		gsl_vector_set(x, 0, d_p[0]);
		gsl_vector_set(x, 1, d_p[1]);
		gsl_vector_set(x, 2, d_p[2]);
		gsl_vector_set(x, 3, d_p[3]);
		gsl_vector_set(x, 4, d_p[4]);
		gsl_vector_set(x, 5, d_p[5]);

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(6);
		gsl_vector_set_all(ss, 1.0);

		/* Initialize method and iterate */
		minex_func.n = 6;
		minex_func.f = convolv_3invG_nov_loglikelihood;
		minex_func.params = (void *) &config;

		s = gsl_multimin_fminimizer_alloc(T, 6);
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

			printf("ll=%g [%.8f %.8f %.8f %.8f %.8f %.8f] size=%.3f %.3fs\n\n",
				s->fval,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				gsl_vector_get(s->x, 2),
				gsl_vector_get(s->x, 3),
				gsl_vector_get(s->x, 4),
				gsl_vector_get(s->x, 5),
				size,
				((float)t) / CLOCKS_PER_SEC);
		} while (status == GSL_CONTINUE && iter < 10000);

		d_p[0] = fabs(gsl_vector_get(s->x, 0));
		d_p[1] = fabs(gsl_vector_get(s->x, 1));
		d_p[2] = fabs(gsl_vector_get(s->x, 2));
		d_p[3] = fabs(gsl_vector_get(s->x, 3));
		d_p[4] = fabs(gsl_vector_get(s->x, 4));
		d_p[5] = fabs(gsl_vector_get(s->x, 5));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		printf("  p=[%f %f %f %f %f %f]\n", d_p[0], d_p[1], d_p[2],
			d_p[3], d_p[4], d_p[5]);

		for (itmp = 0; itmp < 6; itmp++) {
			c_pd[emgfitnomle + 20 * itmp] = d_p[itmp];
		}


		threestage_adapt(data, d_p[0], d_p[1], d_p[2], d_p[3], d_p[4], d_p[5], l, data_size);

		double l_sum = (double)loglikelihood(l, data_size);

		/* FIXME FIXME FIXME we are not reporting the real values here */
		printf("  l=%f hp=%f flag=%f E=%f WARNING HP FLAG AND E ARE BOGUS\n\n", l_sum, 9999.0, 9999.0, 9999.0);

	}

	printf("Recalculating canidate solutions with smaller stepsize\n");

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (emgfitnomle = 0; emgfitnomle < 20; emgfitnomle++) {


		threestage_adapt(data, c_pd[emgfitnomle], c_pd[20 + emgfitnomle],
			c_pd[40 + emgfitnomle], c_pd[60 + emgfitnomle],
			c_pd[80 + emgfitnomle],
			c_pd[100 + emgfitnomle], l, data_size);

		double l_sum = (double)loglikelihood(l, data_size);

		c_flag[emgfitnomle] = l_sum;
	}

	emgfitnomle = 1;
	mtmp = c_flag[0];
	itmp = 0;

	if (emgfitnomle < 20) {
		while (emgfitnomle + 1 < 21) {
			if (c_flag[emgfitnomle] > mtmp) {
				mtmp = c_flag[emgfitnomle];
				itmp = emgfitnomle;
			}

			emgfitnomle++;
		}
	}

	printf("max_ld=%f row_ld=%f\n", mtmp, (double)itmp + 1);

	FREE(l);
}

double convolv_3invG_nov_loglikelihood(const gsl_vector * v, void *params)
{
	configStruct config = *(configStruct *)params;

	distType * data = config.data;
	int data_size = config.data_size;

    double m1 = gsl_vector_get(v, 0);
    double s1 = gsl_vector_get(v, 1);
    double m2 = gsl_vector_get(v, 2);
    double s2 = gsl_vector_get(v, 3);
    double m3 = gsl_vector_get(v, 4);
    double s3 = gsl_vector_get(v, 5);

    double penalty = 0;
    if (m1 < 0 || s1 < 0 || m2 < 0 || s2 < 0 || m3 < 0 || s3 < 0)
		penalty = 1000;

    m1 = fabs(m1);
    s1 = fabs(s1);
    m2 = fabs(m2);
    s2 = fabs(s2);
    m3 = fabs(m3);
    s3 = fabs(s3);

	distType * Y = (distType *)MALLOC(sizeof(distType)*data_size);

    for (int i = 0; i < data_size; i++) {
		Y[i] = 0;
    }

	threestage_adapt(data, m1, s1, m2, s2, m3, s3, Y, data_size);

	double ll = (double) loglikelihood(Y, data_size);

    double objective = penalty - ll;

	FREE(Y);

    return objective;
}

void threestage_adapt(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], int dataSize) {

	printf("starting (%g %g %g %g %g %g)\n", m1, s1, m2, s2, m3, s3);

	distType gridSize = 0.01;

	// find the largest point in t
	double maxData = 0;
	for (int i = 0; i < dataSize; i++) {
		if (data[i] > maxData)
			maxData = data[i];
	}

	// These represent our loglikelihoods
	double logP0;
	double logP1;

	double E = DBL_MAX;

	threestage_bin(data, m1, s1, m2, s2, m3, s3, Y, dataSize, gridSize);

	logP0 = (double)loglikelihood(Y, dataSize);

	printf("ll=%f gridSize=%.17f E=%.17e EB=%g\n", logP0, gridSize, E, _ERROR_BOUND);

	while (E >= _ERROR_BOUND) {

		gridSize = gridSize * 0.5;	// Shrink the step size

		threestage_bin(data, m1, s1, m2, s2, m3, s3, Y, dataSize, gridSize);

		logP1 = (double)loglikelihood(Y, dataSize);

		E = fabs(logP1 - logP0);

		logP0 = logP1;

		printf("ll=%f gridSize=%.17f E=%.17e EB=%g\n", logP1, gridSize, E, _ERROR_BOUND);
	}
	
	printf("\n");
}

void threestage_bin(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], long dataSize, double gridSize) {

	double maxData = 0;
	for (long i = 0; i < dataSize; i++) {
		if (data[i] > maxData)
			maxData = data[i];
	}

	int partitionLength = (int)(maxData / gridSize) + 1;

	distType *partition = (distType *)MALLOC(partitionLength * sizeof(double));
	distType *x = (distType *)MALLOC(partitionLength * sizeof(distType));
	distType *y = (distType *)MALLOC(partitionLength * sizeof(distType));
	distType *z = (distType *)MALLOC(partitionLength * sizeof(distType));

	// fill the partition
	distType tally = 0;
	for (int i = 0; i < partitionLength; i++) {
		partition[i] = tally;
		tally += gridSize;
	}

	// evaluate the sub-distributions at each point in x
	waldpdf(partition, m1, s1, x, partitionLength);
	waldpdf(partition, m2, s2, y, partitionLength);
	waldpdf(partition, m3, s3, z, partitionLength);

	double logP1;

	threestage_binconv(x, y, z, data, Y, &logP1, partitionLength, dataSize, gridSize);

	FREE(partition);
	FREE(x);
	FREE(y);
	FREE(z);

}