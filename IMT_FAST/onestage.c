/*
onestage.c

Functions pertaining to one stage models

*/



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

/*
Optimize parameters for a one stage model using the data in the config struct
*/
void optimize_onestage(configStruct config) {
	printf("onestagefitnomle\n\n");


	distType * data = config.data;
	long data_size = config.data_size;


	distType * l = (distType *)MALLOC(sizeof(distType)*data_size);

	double ld[16];

	/*
	Prepare seeds. Generate a rectangular cloud of points centered on the method of moments parameters
	*/
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

	/*
	We need somewhere to store the optimized parameters corresponding to the given seeds
	*/
	double optimizedParams[27][3];
	//double le[27];

	/*
	For each seed
	*/
#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (seedIdx = 0; seedIdx < 16; seedIdx++) {
		/*
		Optimize this seed
		*/


		printf("starting seedIdx=%d ", 1 + seedIdx);


		printf("P=[%f %f]\n", paramSeeds[seedIdx][0], paramSeeds[seedIdx][1]);


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

		// Allocate the minimizer and set it up
		s = gsl_multimin_fminimizer_alloc(T, 2);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
		printf("\n");

		distType prevll = 0;		// We need to store the loglikelihood of prior iterations for comparison
		distType ll_delta = 0;		// difference in these loglikelihoods
		double prevm = paramSeeds[seedIdx][0];	// We also store the prior iterations parameters
		double prevs = paramSeeds[seedIdx][1];
		int victory_TOL_X = 0;			// Have we satisfied our change in parameters convergence criterion?
		int victory_TOL_FUN = 0;		// Have we satisfied our loglikelihood convergence criterion?
		double delta_X;					// Change in parameters
		double delta_FUN;				// Change in loglikelihood

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

			if (status)
				break;

			// Get the "size" of the current simplex
			size = gsl_multimin_fminimizer_size(s);

			// Check to see if simplex is small enough to stop iterating
			status = gsl_multimin_test_size(size, TOL_SIZE);

			if (status == GSL_SUCCESS) {
				printf("converged to minimum at\n");
			}

			// Check the distance between this point int the parameter space and the prior one
			delta_X = euclideanDistance(gsl_vector_get(gsl_multimin_fminimizer_x(s), 0), gsl_vector_get(gsl_multimin_fminimizer_x(s), 1), 0, 0, 0, 0, prevm, prevs, 0, 0, 0, 0);
			prevm = gsl_vector_get(gsl_multimin_fminimizer_x(s), 0);
			prevs = gsl_vector_get(gsl_multimin_fminimizer_x(s), 1);

			// Probably obselete safety check
			if (ll_delta < 0)
				delta_FUN = DBL_MAX;
			else
				delta_FUN = ll_delta;

			printf("delta_X=%f delta_FUN=%f\n", delta_X, delta_FUN);

			/*
			Apply convergence criteria
			*/
			if (delta_X < TOL_X) {
				printf("declaring victory TOL_X!\n");
				victory_TOL_X = 1;
			}

			if (delta_FUN < TOL_FUN) {
				printf("declaring victory TOL_FUN!\n");
				victory_TOL_FUN = 1;
			}
			
			if (victory_TOL_X && victory_TOL_FUN) {
				printf("converged on custom convergence critera at\n");
				//status = GSL_SUCCESS;
			} else {
				victory_TOL_FUN = 0;
				victory_TOL_X = 0;
			}

			printf("ll=%.17e [%.17f %.17f] size=%.17f %.3fs\n\n",
				s->fval,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				size,
				((float)t) / CLOCKS_PER_SEC);

		} while (status == GSL_CONTINUE && iter < 10000);

		// Store the optimized parameters
		optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
		optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		printf("  p=[%f %f]\n", optimizedParams[seedIdx][0],
			optimizedParams[seedIdx][1]);

		// Calculate the loglikelihood of our optimized parameters
		//onestagepdf2(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], l);
		wald_adapt(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], l, data_size);

		double ll = (double)loglikelihood(l, data_size);

		printf("  l=%.17f\n\n", ll);
		ld[seedIdx] = ll;
	}

	/* find the best optimized parameter set for all starting seeds tried */
	double maxLikelihood = 0;
	int ind_ld = 0;
	printf("\n");
	for (int i = 0; i < 16; i++) {
		printf("ll=%f ms=%f s=%f\n", ld[i], optimizedParams[i][0], optimizedParams[i][1]);
		if (ld[i] < maxLikelihood) {
			maxLikelihood = ld[i];
			ind_ld = i;
		}
	}

	printf("\nmax_ld=%f row_ld=%d\n", maxLikelihood, ind_ld);
	printf("pd_max=[%f %f]\n\n", optimizedParams[ind_ld][0],
		optimizedParams[ind_ld][1]);

	FREE(l);
}

/*
Find and return the loglikelihood for a onestage model with the parameters in v given the data in params (this is unfortunate naming)
*/
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

/*
Evaluate the onestage model with parameters mu and s for each point in data[] into Y[] indirectly using adapatation to adjust the gridSize to ensure acceptable error
*/
void wald_adapt(const distType data[], double mu, double s, distType Y[], long data_size) {

	clock_t t;
	t = clock();

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
	printf("  gridSize=%.17f E=%.17e eThr=%.17e ll=%.17e\n", gridSize, E, (double)_ERROR_BOUND, ll_previous);
#endif

	while (E >= _ERROR_BOUND) {
		gridSize = gridSize * 0.5;

		wald_bin(data, mu, s, Y, data_size, gridSize);

		ll_current = loglikelihood(Y, data_size);

		E = fabs(ll_current - ll_previous);

#ifdef _VERBOSE
		printf("  gridSize=%.17f E=%.17e eThr=%.17e ll=%.17e\n", gridSize, E, (double)_ERROR_BOUND, ll_current);
#endif

		ll_previous = ll_current;
	}

#ifdef _VERBOSE
	printf("%fs]\n", ((float)t) / CLOCKS_PER_SEC);
#endif
	return;
}

/*
Evaluate the onestage model with parameters mu and s for each point in data[], the integrated probability of the sampling bin containing the point, into Y[] using the given gridSize
*/
void wald_bin(const distType data[], double mu, double s, distType Y[], long dataSize, double gridSize) {
	
	printf("dataSize=%ld\n", dataSize);

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

/*
Evaluate the onestage model with parameters mu and s for each point in data[] into Y[] by direct evaluation
*/
void waldpdf(const distType data[], double mu, double s, distType Y[], long dataSize)
{
	
	beginTraceFun("waldpdf");


    // Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t)));
    // https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution

	int badParams = 0;
	if (mu <= 0.0 || s <= 0.0) {
		printf("ERROR: waldpdf:nonpositive mu or s!\n");
		badParams = 1;
	}

	/*
	if (s == 0)
		s = distMin;
	if (mu == 0)
		mu = distMin;
	*/

	distType a, b, bExp;

	int underflow = 0;
	int underflowBegin;
	double percentUnderflow;

	for (long i = 0; i < dataSize; i++) {


		/*
		printf("LOG OF GSL_LOG_DBL_MIN %f %f", log(GSL_LOG_DBL_MIN), GSL_LOG_DBL_MIN);

		
		distType logwald = -(1 / 3)*log(2 * M_PI*s*s*data[i]) - pow(1 - mu * data[i], 2) / (2 * s*s*data[i]);

		if (logwald < distMin) {
			//printf("UNDERFLOW ON LOGWALD, replacing with zero\n");
			//Y[i] = 0;
			Y[i] = exp(logwald);
			continue;
		} else {
			Y[i] = exp(logwald);
			continue;
		}


		*/

		int anyError = 0;

		errno = 0;
		a = 1.0 / (s * pow(6.2831853071795862 * pow(data[i], 3.0), 0.5));

		if (errno == ERANGE) {
			printf("ERROR: waldpdf:(a) range error! ");
			printf("data[%ld]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}
		if (errno == EDOM) {
			printf("ERROR: waldpdf:(a) domain error! ");
			printf("data[%ld]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}

		if (a == HUGE_VAL) {
			a = distMax;

		}


		errno = 0;
		b = -((pow(mu * data[i] - 1.0, 2.0)) / (2.0 * s * s * data[i]));

		if (errno == EDOM) {
			printf("ERROR: waldpdf:(b) domain error! ");
			printf("data[%ld]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}
		if (errno == ERANGE) {
			printf("ERROR: waldpdf:(b) domain error! ");
			printf("data[%ld]=%g mu=%g s=%g ", i, data[i], mu, s);
			anyError = 1;
		}

#ifdef _DIST_SINGLE
		if (b >= FLT_MAX_EXP) {
			printf("ERROR: waldpdf:Predict exp(%f) overflow for type float! data[%d]=%g mu=%g s=%g\n", b, i, data[i], mu, s);
			b = FLT_MAX_EXP - 1;
		}
		if (b <= FLT_MIN_EXP) {
			b = FLT_MIN_EXP + 1;
			if (!underflow) {
				underflow = 1;
				underflowBegin = i;
				//printf("{UNDERFLOW BEGIN data[%d]=%g ", i, data[i]);
			}
		} else {
			if (underflow) {
				percentUnderflow = 100 * ((double)i - (double)underflowBegin - 1) / ((double)dataSize - 1);
				//printf("END data[%d]=%g dataSize=%d %f%% mu=%g s=%g}\n", i - 1, data[i - 1], dataSize, percentUnderflow, mu, s);
				underflow = 0;
			}
		}
#else
		if (b >= DBL_MAX_EXP) {
			printf("ERROR: waldpdf:Predict exp(%lf) overflow for type double! data[%ld]=%g mu=%g s=%g\n", b, i, data[i], mu, s);
			b = DBL_MAX_EXP - 1;
		}
		if (b <= DBL_MIN_EXP) {
			b = DBL_MIN_EXP + 1;
			if (!underflow) {
				underflow = 1;
				underflowBegin = i;
				//printf("{UNDERFLOW BEGIN data[%d]=%g ", i, data[i]);
			}
		} else {
			if (underflow) {
				percentUnderflow = 100 * ((double)i - (double)underflowBegin - 1) / ((double)dataSize - 1);
				//printf("END data[%d]=%g dataSize=%d %f%% mu=%g s=%g}\n", i - 1, data[i - 1], dataSize, percentUnderflow, mu, s);
				underflow = 0;
			}
		}
#endif	

		errno = 0;
		bExp = exp(b);

		if (errno == EDOM) {
			printf("ERROR: waldpdf:(bExp) domain error! ");
			printf("data[%ld]=%g mu=%g s=%g b=%g ", i, data[i], mu, s, b);
			anyError = 1;
		}

		// This is where the bad things happen, bExp has underflowed and then gets multiplied by a giving a non-underflowed value but with a loss of accuracy
		Y[i] = (distType)(a * bExp);

		if (anyError != 0)
			printf(" Y[%ld]=%f\n", i, Y[i]);

		/* Ok this fixup requires some explaination. Strictly wald(0) is indeterminate.
		We replace it here with zero as that is the practical value and this avoids
		blowing up a log likelihood calculation later! */
		if (isnan(Y[i]) || !isfinite(Y[i])) {
			if (data[i] == 0)
				Y[i] = 0;
			else
				printf("ERROR: NAN or Inf a=%f b=%f data[%ld]=%f Y[%ld]=%f\n", a, b, i, data[i], i, (double)Y[i]);
		}

		if (Y[i] != 0) {
			if (Y[i] >= distMax) {
				printf("ERROR: Denormal value Y[%ld]=%g >= DBL_MAX\n", i, Y[i]);
			} else if (Y[i] <= distMin ){
				//printf("ERROR: Denormal value Y[%d]=%g <= DBL_MIN\n", i, Y[i]);
				Y[i] = 0;
			}
		}

		if (Y[i] < 0) {
			printf("ERROR: waldpdf:InvG(%f)=Y[%ld]<0=%f mu=%f s=%f HOLY NORMALIZATION BATMAN\n", data[i], i, Y[i], mu, s);
			Y[i] = 0;
		}
    }
	if (underflow) {
		percentUnderflow = 100* ((double)dataSize - (double)underflowBegin - 1) / ((double)dataSize -1);
		printf("END data[%ld]=%g dataSize=%ld %f%% mu=%g s=%g}\n", dataSize-1, data[dataSize-1], dataSize, percentUnderflow, mu, s);
		underflow = 0;
	}
	/*
	if (badParams)
		for(int i = 0; i < dataSize; i++)
			printf("Y[%d]=%f ", i, Y[i]);
	*/


	endTraceFun("waldpdf");
}
