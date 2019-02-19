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
#include "gsl/gsl_statistics.h"
#include "utility.h"

#ifdef _WIN32
#include "windows.h"
#endif

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3

typedef struct {
	double f1[2];
} cell_wrap_3;

#endif

#define _GSL_GP_MAX_FIXED
#define _BINNED_MODE

/*
Generate seeds for a threestage convolutional model using the partial method of cumulants.
*/
double ** threestage_seeds(double mean, double variance, int *numSeeds) {
	/*
	This works similarly to twostage_seeds except that we then further split the right mean and right variance again for the second two distributions.
	*/


	int numRatios = 2;

	double ratios[2] = { 0.1, 0.75 };

	int numseeds_twostage;

	// Note that we only want the number of seeds that twostage generated so we throw away the actual seeds afterwards and generate our own
	double ** seeds_twostage = twostage_seeds(mean, variance, &numseeds_twostage);

	*numSeeds = numRatios * numRatios * numseeds_twostage;

	// throw away the existing seeds
	for (int seedIdx = 0; seedIdx < numseeds_twostage; seedIdx++)
		FREE(seeds_twostage[seedIdx]);
	FREE(seeds_twostage);


	double ** seeds = (double **)MALLOC(sizeof(double *) * *numSeeds);

	for (int seedIdx = 0; seedIdx < *numSeeds; seedIdx++) {

		seeds[seedIdx] = (double *)MALLOC(sizeof(double) * 6);

		int meanIdx = ((seedIdx / numseeds_twostage) / numRatios);
		int varIdx = ((seedIdx / numseeds_twostage) % numRatios);
		int twostageIdx = seedIdx % numseeds_twostage;


		double leftMean = mean * ratios[meanIdx];
		double rightMean = mean * fabs(1 - ratios[meanIdx]);

		double leftVariance = variance * ratios[varIdx];
		double rightVariance = variance * fabs(1 - ratios[varIdx]);

		//printf("meanIdx = %d  varIdx = %d  towstageIdx = %d      leftMean = %f  rightMean = %f  leftVariance = %f  rightVariance = %f\n",  meanIdx, varIdx, twostageIdx, leftMean, rightMean, leftVariance, rightVariance);

		double ** seeds_twostage = twostage_seeds(rightMean, rightVariance, &numseeds_twostage);

		seeds[seedIdx][0] = 1 / leftMean;
		seeds[seedIdx][1] = pow(leftVariance / pow(leftMean, 3), 0.5);

		seeds[seedIdx][2] = seeds_twostage[twostageIdx][0];
		seeds[seedIdx][3] = seeds_twostage[twostageIdx][1];

		seeds[seedIdx][4] = seeds_twostage[twostageIdx][2];
		seeds[seedIdx][5] = seeds_twostage[twostageIdx][3];

		for (int seedIdx = 0; seedIdx < numseeds_twostage; seedIdx++)
			FREE(seeds_twostage[seedIdx]);
		FREE(seeds_twostage);
	}
	return seeds;
}

/*
Optimize parameters for a three stage model using the data in the config struct for the given seeds
*/
void optimize_threestage(const distType data[], int data_size, configStruct config) {

	distType * l = (distType *)MALLOC(sizeof(distType)*data_size);

	/*
	Get some seeds
	*/

	double * doubleData = (double *)MALLOC(sizeof(double)*data_size);
	for (int i = 0; i < data_size; i++)
		doubleData[i] = (double)data[i];
	double mean = gsl_stats_mean(doubleData, 1, data_size);
	double variance = gsl_stats_variance(doubleData, 1, data_size);
	FREE(doubleData);
	int numseeds;
	double ** seeds = threestage_seeds(mean, variance, &numseeds);
	printf("\nseeds[%d]:\n", numseeds);
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++)
		printf("%f %f %f %f %f %f\n", seeds[seedIdx][0], seeds[seedIdx][1], seeds[seedIdx][2], seeds[seedIdx][3], seeds[seedIdx][4], seeds[seedIdx][5]);
	printf("\n");

	// Need to store optimized parameters for each seed
	double ** optimizedParams = (double **)MALLOC(sizeof(double *) * numseeds);
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		optimizedParams[seedIdx] = (double *)MALLOC(sizeof(double) * 6);
	}

	//Store iteration count for each seed
	size_t * numIters = (size_t *)MALLOC(sizeof(size_t) * numseeds);

	//Store the runtimes for each seed
	float * runtimes = (float *)MALLOC(sizeof(float) * numseeds);

	printf("threestagefitnomle\n\n");

	distType * loglikelihoods = (distType *)MALLOC(sizeof(distType)*numseeds);

	/*
	For each seed
	*/
#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {

		
		clock_t seed_t = clock();

		/*
		Optimize this seed
		*/

		printf("seed[%d]=[%f %f %f %f %f %f]\n\n", seedIdx, seeds[seedIdx][0], seeds[seedIdx][1], seeds[seedIdx][2], seeds[seedIdx][3], seeds[seedIdx][4], seeds[seedIdx][5]);

		// https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives

		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;	// We are using the nelder-mead implementation
		gsl_multimin_fminimizer *s = NULL;	// s is the minimizer itself
		gsl_vector *ss, *x;					// ss is the step and x is the current position
		gsl_multimin_function minex_func;	// this is the function that the minimizer will minimize

		size_t iter = 0;					// Optimizer iteration count
		int status;							// This lets us know when to stop optimizing or continue
		double size;						// This is the size of the optimizers current simplex

		/* Starting point */
		x = gsl_vector_alloc(6);
		gsl_vector_set(x, 0, seeds[seedIdx][0]);
		gsl_vector_set(x, 1, seeds[seedIdx][1]);
		gsl_vector_set(x, 2, seeds[seedIdx][2]);
		gsl_vector_set(x, 3, seeds[seedIdx][3]);
		gsl_vector_set(x, 4, seeds[seedIdx][4]);
		gsl_vector_set(x, 5, seeds[seedIdx][5]);

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(6);
		gsl_vector_set_all(ss, 1.0);

		/* Initialize method and iterate */
		minex_func.n = 6;
		minex_func.f = convolv_3invG_nov_loglikelihood;
		minex_func.params = (void *)&config;

		// Allocate the minimizer and set it up
		s = gsl_multimin_fminimizer_alloc(T, 6);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
		printf("\n");

		distType prevll = 0;				// We need to store the loglikelihood of prior iterations for comparison
		distType ll_delta = 0;				// difference in these loglikelihoods

		/*
		Iterate the optimizer until we satisfy convergence criteria or we exceed the iteration limit
		*/
		do {
			iter++;


			clock_t t;
			t = clock();

			printf("seed=%d iter=%d\n", (int) seedIdx, (int)iter);
			status = gsl_multimin_fminimizer_iterate(s);

			ll_delta = prevll - s->fval;

			prevll = s->fval;

			t = clock() - t;

			if (status)
				break;

			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, TOL_SIZE);

			if (status == GSL_SUCCESS) {
				printf("converged to minimum at\n");
			}

			if (fabs(ll_delta) < TOL_FUN) {
				printf("declaring victory!\n");
				//status = GSL_SUCCESS;
			}

#ifdef _WIN32
			HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
			CONSOLE_SCREEN_BUFFER_INFO consoleInfo;
			WORD saved_attributes;

			/* Save current attributes */
			GetConsoleScreenBufferInfo(hConsole, &consoleInfo);
			saved_attributes = consoleInfo.wAttributes;

			SetConsoleTextAttribute(hConsole, BACKGROUND_INTENSITY);
#endif
			printf("ll=%g [%.8f %.8f %.8f %.8f %.8f %.8f] size=%.3f %.3fs",
				s->fval,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				gsl_vector_get(s->x, 2),
				gsl_vector_get(s->x, 3),
				gsl_vector_get(s->x, 4),
				gsl_vector_get(s->x, 5),
				size,
				((float)t) / CLOCKS_PER_SEC);
#ifdef _WIN32
			/* Restore original attributes */
			SetConsoleTextAttribute(hConsole, saved_attributes);
#endif

			printf("\n\n");

		} while (status == GSL_CONTINUE && iter < 10000);

		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
		optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));
		optimizedParams[seedIdx][3] = fabs(gsl_vector_get(s->x, 3));
		optimizedParams[seedIdx][4] = fabs(gsl_vector_get(s->x, 4));
		optimizedParams[seedIdx][5] = fabs(gsl_vector_get(s->x, 5));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		threestage_adapt(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], optimizedParams[seedIdx][4], optimizedParams[seedIdx][5], l, data_size);

		double l_sum = (double)loglikelihood(l, data_size);

		loglikelihoods[seedIdx] = l_sum;

		seed_t = clock() - seed_t;

		runtimes[seedIdx] = ((float)seed_t) / CLOCKS_PER_SEC;

#ifdef _WIN32
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
		CONSOLE_SCREEN_BUFFER_INFO consoleInfo;
		WORD saved_attributes;

		/* Save current attributes */
		GetConsoleScreenBufferInfo(hConsole, &consoleInfo);
		saved_attributes = consoleInfo.wAttributes;

		SetConsoleTextAttribute(hConsole, BACKGROUND_RED);
#endif
		printf("finished seedIdx=%d p=[%f %f %f %f %f %f] ll=%f time=%fs",
			seedIdx, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2], optimizedParams[seedIdx][3],
			optimizedParams[seedIdx][4], optimizedParams[seedIdx][5],
			l_sum, ((float)seed_t) / CLOCKS_PER_SEC);

#ifdef _WIN32
		/* Restore original attributes */
		SetConsoleTextAttribute(hConsole, saved_attributes);
#endif

		numIters[seedIdx] = iter;

		printf("\n\n");
		
	}

	// Find the best log likelihood
	distType max_ld = -distMax;
	int row_id = 0;
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		printf("row_id=%d [%f %f %f %f %f %f] -> [%f %f %f %f %f %f] iters=%d runtime=%fs ll=%f\n", seedIdx,
			seeds[seedIdx][0], seeds[seedIdx][1], seeds[seedIdx][2], seeds[seedIdx][3], seeds[seedIdx][4], seeds[seedIdx][5],
			optimizedParams[seedIdx][0], optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2], optimizedParams[seedIdx][3],
			optimizedParams[seedIdx][4], optimizedParams[seedIdx][5],
			(int)numIters[seedIdx], runtimes[seedIdx], loglikelihoods[seedIdx]);
		if (loglikelihoods[seedIdx] > max_ld) {
			max_ld = loglikelihoods[seedIdx];
			row_id = seedIdx;
		}
	}

	printf("max_ld=%f row_ld=%d [%f %f %f %f %f %f]\n", max_ld, row_id, optimizedParams[row_id][0], optimizedParams[row_id][1], optimizedParams[row_id][2], optimizedParams[row_id][3], optimizedParams[row_id][4], optimizedParams[row_id][5]);

	FREE(l);

	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++)
		FREE(seeds[seedIdx]);
	FREE(seeds);

	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		FREE(optimizedParams[seedIdx]);
	}
	FREE(optimizedParams);


}

/*
Find and return the loglikelihood for a threestage model with the parameters in v given the data in params
*/
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

	double ll;

	
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

		ll = (double)loglikelihood(Y, data_size);

		FREE(Y);
	
	
    double objective =  ll;

	

    return objective;
}

/*
Evaluate the threestage model with parameters m1, s1, m2, s2, m3, s3 for each point in data[] into Y[] indirectly
using adapatation to adjust the gridSize to ensure acceptable error
*/
void threestage_adapt(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], int dataSize) {

	beginTraceFun("threestage_adapt");

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


	// find the variance for each sub-distributions
	double v1 = (s1 * s1) / (m1 * m1 * m1);
	double v2 = (s2 * s2) / (m2 * m2 * m2);
	double v3 = (s3 * s3) / (m3 * m3 * m3);

	// Ensure that the lowest variance distribution is never the last one to prevent extra work
	double m_tmp, s_tmp;
	if (v3 < v2) {
		m_tmp = m2;
		s_tmp = s2;
		m2 = m3;
		s2 = s3;
		m3 = m_tmp;
		s3 = s_tmp;
	}

	printf("starting (%g %g %g %g %g %g)\n", m1, s1, m2, s2, m3, s3);

	//threestage_bin(data, m1, s1, m2, s2, m3, s3, Y, dataSize, gridSize);
	threestage_double_adapt_bin(data, m1, s1, m2, s2, m3, s3, Y, dataSize, gridSize);

	logP0 = (double)loglikelihood(Y, dataSize);

	printf("ll=%f gridSize=%.17f E=%.17e EB=%g threestage\n", logP0, gridSize, E, (double)_ERROR_BOUND);

	while (E >= _ERROR_BOUND) {

		gridSize = gridSize * 0.5;	// Shrink the step size

		printf("Shrinking gridSize... gridSize = %f\n", gridSize);

		//threestage_bin(data, m1, s1, m2, s2, m3, s3, Y, dataSize, gridSize);
		threestage_double_adapt_bin(data, m1, s1, m2, s2, m3, s3, Y, dataSize, gridSize);

		logP1 = (double)loglikelihood(Y, dataSize);

		E = fabs(logP1 - logP0);

		logP0 = logP1;

		printf("ll=%f gridSize=%.17f E=%.17e EB=%g threestage\n", logP1, gridSize, E, (double)_ERROR_BOUND);
	}
	
	endTraceFun("threestage_adapt");
}

/*
Evaluate the threestage model with parameters m1, s1, m2, s2, m3, s3 for each point in data[], the integrated
probability of the sampling bin containing the point, into Y[] using the given gridSize
*/
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

/*
Evaluate the threestage model with parameters m1, s1, m2, s2, m3, s3 for each point in data[], the integrated
probability of the sampling bin containing the point, into Y[] using the given gridSize where the first two
distributions are convolved using a seperate adapatation than the last convolution
*/
void threestage_double_adapt_bin(const distType data[], double m1, double s1, double m2, double s2, double m3, double s3, distType Y[], long dataSize, double gridSize) {


	beginTraceFun("threestage_double_adapt_bin");

	double maxData = 0;
	for (long i = 0; i < dataSize; i++) {
		if (data[i] > maxData)
			maxData = data[i];
	}

	int partitionLength = (int)(maxData / gridSize) + 1;

	distType *partition = (distType *)MALLOC(partitionLength * sizeof(double));
	distType *x = (distType *)MALLOC(partitionLength * sizeof(distType));
	distType *y = (distType *)MALLOC(partitionLength * sizeof(distType));
#ifndef _MULTI_ADAPT
	distType *z = (distType *)MALLOC(partitionLength * sizeof(distType));
#endif

	// fill the partition
	distType tally = 0;
	for (int i = 0; i < partitionLength; i++) {
		partition[i] = tally;
		tally += gridSize;
	}

	// evaluate the sub-distributions at each point in x
	waldpdf(partition, m3, s3, y, partitionLength);
#ifdef _MULTI_ADAPT
	conv2waldpdf(partition, m1, s1, m2, s2, x, gridSize, 1, partitionLength, NORMALIZATION);
#else
	waldpdf(partition, m1, s1, x, partitionLength);
	waldpdf(partition, m3, s3, z, partitionLength);
#endif

	double logP1;

#ifdef _MULTI_ADAPT
	binned_conv(x, y, data, partition, Y, &logP1, partitionLength, dataSize, gridSize);
#else
	threestage_binconv(x, y, z, data, Y, &logP1, partitionLength, dataSize, gridSize);
#endif

	FREE(partition);
	FREE(x);
	FREE(y);
#ifndef _MULTI_ADAPT
	FREE(z);
#endif

	endTraceFun("threestage_double_adapt_bin");

}
