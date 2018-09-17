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
#include "utility.h"
#include "gsl/gsl_statistics.h"


#define _VERBOSE
//#define _TWOSTAGE_BINNED_MODE

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3
#endif

typedef struct {
	double f1[2];
} cell_wrap_3;


/*
Generate seeds for a twostage convolutional model using the partial method of cumulants.
*/
double ** twostage_seeds(double mean, double variance, int *numSeeds) {
	/*
	 It is a property of convolution that for any nth cumulant, the nth cumulant of the convolution is the sum of the nth cumulants of
	the input distributions. With a regular method of cumulants we would match these cumulants to the data directly. However, because
	the noise in the computed cumulants for the data gets large for higher order cumulants and small data set sizes, we really only want
	to work with the first few where we have a reasonable confidence. As such, we simply generate a mesh of "ratios" and we split the
	cumulants value between the left and right distribution according to these ratios.
	*/


	int numRatios = 2;

	double ratios[2] = { 0.2, 0.4 };

	*numSeeds = numRatios * numRatios;

	double ** seeds = (double **) MALLOC(sizeof(double *) * *numSeeds);

	for (int seedIdx = 0; seedIdx < *numSeeds; seedIdx++) {

		seeds[seedIdx] = (double *)MALLOC(sizeof(double) * 4);

		int meanIdx = seedIdx % numRatios;
		int varIdx = seedIdx / numRatios;

		double leftMean = mean * ratios[meanIdx];
		double rightMean = mean * (1 - ratios[meanIdx]);

		double leftVariance = variance * ratios[varIdx];
		double rightVariance = variance * (1 - ratios[varIdx]);


		seeds[seedIdx][0] = 1 / leftMean;
		seeds[seedIdx][1] = pow( leftVariance / pow(leftMean, 3), 0.5);

		seeds[seedIdx][2] = 1 / rightMean;
		seeds[seedIdx][3] = pow(rightVariance / pow(rightMean, 3), 0.5);

	}
	return seeds;
}

/*
Optimize parameters for a one stage model using the data in the config struct for the given seeds
*/
void optimize_twostage(int numseeds, double ** seeds, configStruct config) {
	printf("twostagefitnomle\n");

	int data_size = config.data_size;
	distType * data = config.data;

	/*
	We need somewhere to store the optimized parameters corresponding to the given seeds
	*/
	double ** optimizedParams = (double **)MALLOC(sizeof(double *) * numseeds);
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		optimizedParams[seedIdx] = (double *)MALLOC(sizeof(double) * 4);
	}


	/*
	For each seed
	*/
#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		/*
		Optimize this seed
		*/

		// This stores the likelihoods of the distribution with our starting seed parameters
		distType * seedll = (distType *)MALLOC(sizeof(distType)*data_size);

		// Calculate the likelihoods of observing the specified data points
		conv2waldpdf(data, seeds[seedIdx][0],
			seeds[seedIdx][1], seeds[seedIdx][2],
			seeds[seedIdx][3], seedll, 0.01, 1,
			data_size);

		// Find the loglikelihood of the set of observations taken together
		double seedll_sum = (double)loglikelihood(seedll, data_size);

		FREE(seedll);

		printf("\nstarting seedIdx=%d p=[%f %f %f %f] ll=%f\n", seedIdx, seeds[seedIdx][0],
			seeds[seedIdx][1], seeds[seedIdx][2],
			seeds[seedIdx][3], seedll_sum);

		/*
		Set up the GSL derivative free optimizer
		https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives
		*/

		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;	// We are using the nelder-mead implementation
		gsl_multimin_fminimizer *s = NULL;	// s is the minimizer itself
		gsl_vector *ss, *x;					// ss is the step and x is the current position
		gsl_multimin_function minex_func;	// this is the function that the minimizer will minimize
		size_t iter = 0;					// Optimizer iteration count
		int status;							// This lets us know when to stop optimizing or continue
		double size;						// This is the size of the optimizers current simplex

		/* Set the starting seed parameters */
		x = gsl_vector_alloc(4);
		gsl_vector_set(x, 0, seeds[seedIdx][0]);
		gsl_vector_set(x, 1, seeds[seedIdx][1]);
		gsl_vector_set(x, 2, seeds[seedIdx][2]);
		gsl_vector_set(x, 3, seeds[seedIdx][3]);

		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc(4);
		gsl_vector_set_all(ss, 1.0);

		/* Tell minimizer what function to minimize and with what fixed data */
		minex_func.n = 4;
		minex_func.f = convolv_2invG_adapt_nov_loglikelihood;
		minex_func.params = (void *) &config;

		// Allocate the minimizer and set it up
		s = gsl_multimin_fminimizer_alloc(T, 4);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
		printf("\n");


		distType prevll = 0;	// We need to store the loglikelihood of prior iterations for comparison
		distType ll_delta = 0;	// difference in these loglikelihoods
		double prevm1 = seeds[seedIdx][0];	// We also store the prior iterations parameters
		double prevs1 = seeds[seedIdx][1];
		double prevm2 = seeds[seedIdx][2];
		double prevs2 = seeds[seedIdx][3];
		int victory_TOL_X = 0;		// Have we satisfied our change in parameters convergence criterion?
		int victory_TOL_FUN = 0;	// Have we satisfied our loglikelihood convergence criterion?

		/*
		Iterate the optimizer until we satisfy convergence criteria or we exceed the iteration limit
		*/
		do {
			iter++;

			clock_t t;
			t = clock();
			
			printf("iter=%d\n", (int)iter);

			// Iterate the optimizer one step
			status = gsl_multimin_fminimizer_iterate(s);

			// Look at the new loglikelihood and compare it with the old one
			ll_delta = prevll - s->fval;
			prevll = s->fval;

			t = clock() - t;

			// Get the "size" of the current simplex
			size = gsl_multimin_fminimizer_size(s);


			if (iter % 1 == 0)
				printf("\n\n");

			// If the optimizer wants us to stop, stop!
			if (status)
				break;

			// Check to see if simplex is small enough to stop iterating
			status = gsl_multimin_test_size(size, TOL_SIZE);


			if (status == GSL_SUCCESS) {
				printf("converged on default convergence criteria at\n");
			}

			// Check the distance between this point int the parameter space and the prior one
			if (euclideanDistance(gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3), 0, 0, prevm1, prevs1, prevm2, prevs2, 0, 0) < TOL_X) {
				printf("declaring victory TOL_X!\n");
				victory_TOL_X = 1;
			}
			prevm1 = gsl_vector_get(s->x, 0);
			prevs1 = gsl_vector_get(s->x, 1);
			prevm2 = gsl_vector_get(s->x, 2);
			prevs2 = gsl_vector_get(s->x, 3);

			// Check the difference in the loglikelihood between this point and the prior one
			if (fabs(ll_delta) < TOL_FUN) {
				printf("declaring victory TOL_FUN!\n");
				victory_TOL_FUN = 1;
				//status = GSL_SUCCESS;
			}

			/*
			if (victory_TOL_X && victory_TOL_FUN) {
				printf("converged on custom convergence critera at\n");
				//status = GSL_SUCCESS;
			}
			else {
				victory_TOL_FUN = 0;
				victory_TOL_X = 0;
			}
			*/

			printf("ll=%g [%f %f %f %f] size=%.3f %.3fs \n\n\n",
				s->fval,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				gsl_vector_get(s->x, 2),
				gsl_vector_get(s->x, 3),
				size,
				((float)t) / CLOCKS_PER_SEC);

		} while (status == GSL_CONTINUE && iter < 10000);
		printf("\n");


		// Store the optimized parameters
		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
		optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));
		optimizedParams[seedIdx][3] = fabs(gsl_vector_get(s->x, 3));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		// Calculate the loglikelihood of our optimized parameters
		distType * ll = (distType *)MALLOC(sizeof(distType)*data_size);
		conv2waldpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], ll, 0.01, 1, data_size);
		double l_sum = (double)loglikelihood(ll, data_size);

		FREE(ll);

		printf("\nfinished seedIdx=%d p=[%f %f %f %f] ll=%f\n", seedIdx, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], l_sum);
	}

	/*  we previously optimized with a larger step size, recalculate with */
	/*  a smaller stepsize after the fact */
	printf("\nrecalculating canidate solutions with smaller stepsize\n");
	
	distType * loglikelihoods = (distType *)MALLOC(sizeof(distType)*numseeds);
	distType * likelihoods = (distType *)MALLOC(sizeof(distType)*data_size);

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		conv2waldpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2], optimizedParams[seedIdx][3],
			likelihoods, 0.001, 1, data_size);

		loglikelihoods[seedIdx] = loglikelihood(likelihoods, data_size);

		printf("id=%d p=[%f %f %f %f] ll=%f\n", seedIdx, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1],
			optimizedParams[seedIdx][2], optimizedParams[seedIdx][3], loglikelihoods[seedIdx]);
	}

	FREE(likelihoods);

	// Find the best log likelihood
	distType max_ld = -distMax;
	int row_id = 0;
	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		if (loglikelihoods[seedIdx] > max_ld) {
			max_ld = loglikelihoods[seedIdx];
			row_id = seedIdx;
		}
	}

	FREE(loglikelihoods);

	printf("\nBest fit: row_id=%d\n", row_id);
	printf("loglikelihood=%f ", max_ld);
	printf("p=[ %f %f %f %f ]\n", optimizedParams[row_id][0], optimizedParams[row_id][1], optimizedParams[row_id][2], optimizedParams[row_id][3]);

	for (int seedIdx = 0; seedIdx < numseeds; seedIdx++) {
		FREE(optimizedParams[seedIdx]);
	}
	FREE(optimizedParams);

	return;
}

/*
Find and return the loglikelihood for a twostage model with the parameters in v given the data in params
*/
double convolv_2invG_adapt_nov_loglikelihood(const gsl_vector * v, void * params)
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

	distType * Y = (distType *)MALLOC(sizeof(distType)*data_size);

    conv2waldpdf(data, m1, s1, m2, s2, Y, 0.01, 1, data_size);

	double ll = (double) loglikelihood(Y, data_size);

	FREE(Y);

    return penalty - ll;
}

/*
Evaluate the twostage model with parameters m1, s1, m2, s2 for each point in data[] into convolvedPDF[] indirectly
using adapatation to adjust the gridSize to ensure acceptable error with a starting grid spacing of hand a flag,
adapativeMode ,for disabling adapation
*/
void conv2waldpdf(const distType data[], double m1, double s1, double m2, double s2, distType convolvedPDF[], double h, int adaptiveMode, int size_XY) {



	beginTraceFun("conv2waldpdf");

	double initialH = h;

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

    /*
	Reorder m and s so that the sub-distribution with the smallest
    variance comes first.  So the fist part might be approximated as a Dirac delta.
	*/
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
    double maxData = 0;
    for (int i = 0; i < size_XY; i++) {
	if (data[i] > maxData)
	    maxData = data[i];
    }

    // These represent our loglikelihoods
    double logP0 = 0;
    double logP1 = 0;

	double totalProbability = 0;
	double priorTotalProbability = 0;
	double deltaTotalProbability = 0;

    /*
	If the first pdf is very concentrated, check to see if it can be approximated as a
	point-mass distribution
	 */
    if (sd_a < 0.01 && checktailmass(m_a, s_a, m_b, s_b, eps, sd_a, sd_b)) {
#ifdef _VERBOSE
		printf("using dirac delta approximation\n");
#endif

		/*
		If there is not much probability in the tails of the convolution,
		we use a Dirac Delta for part 1.
		Flag is set to 1 to indicate this.
		*/
		double lag = 1 / m_a;

		// Use the log model approximation rather than do the convolution
		waldlagpdf(data, m_b, s_b, lag, convolvedPDF, size_XY);

		flag = 1;
    } else {

		/*
		Calculuate the binned convolution
		*/
		twostage_bin(data, m[0], s[0], m[1], s[1], convolvedPDF, size_XY, h);

		//logP0 = (double)loglikelihood(convolvedPDF, size_XY);

		//printf("ll=%g h=%.17f E=%g EB=%g twostage\n", logP0, h, E, _ERROR_BOUND);

		priorTotalProbability = checkNormalization(convolvedPDF, size_XY, initialH);

		deltaTotalProbability = 1;

		//while (E >= _ERROR_BOUND) {
		while( deltaTotalProbability > 0.000001 ){
			/*
			Try again with a smaller grid size until the error is acceptable
			*/

			h = h * 0.5;	// Shrink the step size

			printf("Shrinking h... h = %f\n", h);

			twostage_bin(data, m[0], s[0], m[1], s[1], convolvedPDF, size_XY, h);

			//logP1 = (double)loglikelihood(convolvedPDF, size_XY);

			//E = fabs(logP1 - logP0);

			//logP0 = logP1;

			totalProbability = checkNormalization(convolvedPDF, size_XY, initialH);

			deltaTotalProbability = fabs(priorTotalProbability - totalProbability);

			printf("totalProbability = %.17f deltaTotalProbability = %.17f\n", totalProbability, deltaTotalProbability);

			priorTotalProbability = totalProbability;


#ifdef _VERBOSE
			printf("ll=%g h=%.17f E=%g EB=%g twostage\n", logP1, h, E, (double)_ERROR_BOUND);
#endif
		}
    }
#ifdef _VERBOSE
	//printf("\n");
#endif

	endTraceFun("conv2waldpdf");
}

/*
Evaluate the twostage model with parameters m1, s1, m2, s2 for each point in data[], the integrated probability of the sampling bin containing the point, into Y[] using the given gridSize
*/
void twostage_bin(const distType data[], double m1, double s1, double m2, double s2, distType Y[], long dataSize, double gridSize) {

	beginTraceFun("twostage_bin");


	// Find the largest value in data
	double maxData = 0;
	for (long i = 0; i < dataSize; i++) {
		if (data[i] > maxData)
			maxData = data[i];
	}

	// The grids need to span the observed range of the data
	int partitionLength = (int)(maxData / gridSize) + 1;

	distType *partition = (distType *)MALLOC(partitionLength * sizeof(distType));
	distType *y = (distType *)MALLOC(partitionLength * sizeof(distType));
	distType *z = (distType *)MALLOC(partitionLength * sizeof(distType));

	// fill the partition
	distType tally = 0;
	for (int i = 0; i < partitionLength; i++) {
		partition[i] = tally;
		tally += gridSize;
	}

	// evaluate the sub-distributions at each point in x
	waldpdf(partition, m1, s1, y, partitionLength);
	waldpdf(partition, m2, s2, z, partitionLength);

	double logP1;

#ifdef _TWOSTAGE_BINNED_MODE
	binned_conv(z, y, data, partition, Y, &logP1, partitionLength, dataSize, gridSize);
#else
	nn_conv(z, y, data, partition, Y, &logP1, partitionLength, dataSize, gridSize);
#endif

	FREE(partition);
	FREE(y);
	FREE(z);


	endTraceFun("twostage_bin");
}