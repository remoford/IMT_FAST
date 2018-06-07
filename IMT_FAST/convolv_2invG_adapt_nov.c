/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * convolv_2invG_adapt_nov.c
 *
 * Code generation for function 'convolv_2invG_adapt_nov'
 *
 */

/* Include files */
#include "IMT_analysis_April2017.h"
#include "convolv_2invG_adapt_nov.h"
#include "emgpdf.h"
#include "onestagepdf_lag.h"
#include "approxconvolv.h"
#include "onestagepdf2.h"
#include "tailmass.h"
#include "gsl/gsl_multimin.h"
#include "math.h"


//#define _CONV2INVG
#define _CONV2WALD
//#define _VERBOSE
//#define _PARALLEL_PDF


/* Function Definitions */

double convolv_2invG_adapt_nov_loglikelihood(const gsl_vector *v, void *params)
{
	double *data = (double *)params;

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

	double Y[266];

	double Y_WALD[266];

	double Y_INVG[266];

#ifdef _CONV2INVG
#ifdef _VERBOSE
	printf("starting convolv_2invG_adapt_nov\n");
#endif
	convolv_2invG_adapt_nov(m1, s1, m2, s2, Y_INVG);
#ifdef _VERBOSE
	printf("convolv_2invG_adapt_nov = \n");
	for (int i = 0; i < 266; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%f ", Y_INVG[i]);
	}
	printf("\n\n");
#endif
#endif

#ifdef _CONV2WALD
#ifdef _VERBOSE
	printf("starting conv2waldpdf\n");
#endif
	conv2waldpdf(data, m1, s1, m2, s2, Y_WALD, 0.01, 1, 266);
#ifdef _VERBOSE
	printf("conv2waldpdf = \n");
	for (int i = 0; i < 266; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%f ", Y_WALD[i]);
	}
	printf("\n\n");
#endif
#endif

#ifdef _CONV2INVG
#ifdef _CONV2WALD
	printf("difference between convolv_2invG_adapt_nov and conv2wald = \n");
	for (int i = 0; i < 266; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%.17f ", Y_INVG[i] - Y_WALD[i]);
	}
	printf("\n\n");
#endif
#endif

#ifdef _CONV2WALD
	for (int i = 0; i < 266; i++)
		Y[i] = Y_WALD[i];
#endif

#ifdef _CONV2INVG
	for (int i = 0; i < 266; i++)
		Y[i] = Y_INVG[i];
#endif


	double loglikelihood = 0;
	for (int i = 0; i < 266; i++) {
		loglikelihood += log(Y[i]);
		//printf("data[%d]=%.17f Y[%d]=%.17f log(%.17f)=%.17f\n", i, data[i], i, Y[i], Y[i], log(Y[i]));

	}
	//printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nloglikelihood = %f\n\n\n", loglikelihood);

	//exit(1);

	return penalty - loglikelihood;
}

void conv2waldpdf(const double X[], double m1, double s1, double m2, double s2, double Y[], double h, int adaptiveMode, int size_XY)
{

	int flag = 0; // remember if we applied the approximation
	double eps = 0.01; // a constant that is used in determining if the Dirac approximation should be applied.

	double E;
	if (adaptiveMode)
		E = FP_INFINITE;
	else
		E = 0;

	//int N = 266; // number of points to be evaluated
	int N = size_XY;

	double m_a = m1;
	double m_b = m2;
	double s_a = s1;
	double s_b = s2;

	double m[2];
	double s[2];
	m[0] = m1; m[1] = m2;
	s[0] = s1; s[1] = s2;

	// find the variance for both sub-distributions
	double v_a = (s_a * s_a) / (m_a * m_a * m_a);
	double v_b = (s_b * s_b) / (m_b * m_b * m_b);

	// find the standard deviation for each sub-distribution
	double sd_a = pow(v_a, 0.5);
	double sd_b = pow(v_b, 0.5);

	/* reorder m and s so that the sub-distribution with the smallest
	variance comes first.  So the fist part might be approximated as a Dirac delta.*/
	if (sd_a > sd_b) {
		double tmp;
		tmp = m_a; m_a = m_b; m_b = tmp;
		tmp = s_a; s_a = s_b; s_b = tmp;
		tmp = v_a; v_a = v_b; v_b = tmp;
		tmp = sd_a; sd_a = sd_b; sd_b = tmp;
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
		//onestagepdf_lag(X, m_b, s_b, lag, Y);
		waldlagpdf(X, m_b, s_b, lag, Y, size_XY);
		flag = 1;
	}
	else {

		/* produce a range of evenly spaced points to evaluate at between
		0 and Maxt with step size h. the even spacing is important when
		calculating the convolution later */
		int partitionLength = (int)(maxX / h);

		// This represents the partition
		double * x = _mm_malloc(partitionLength * sizeof(double), 32);
		//double * x = malloc(partitionLength * sizeof(double));

		// There are our two wald distributions
		double * y = _mm_malloc(partitionLength * sizeof(double), 32);
		//double * y = malloc(partitionLength * sizeof(double));


		double * z = _mm_malloc(partitionLength * sizeof(double), 32);
		//double * z = malloc(partitionLength * sizeof(double));


		// fill the partition
		/*
		for (int i = 0; i < partitionLength; i++)
		x[i] = i*h;
		*/
		double tally = 0;
		for (int i = 0; i < partitionLength; i++) {
			x[i] = tally;
			tally += h;
		}

		// evaluate the sub-distributions at each point in x

		/*
		waldpdf(x, m_a, s_a, y, partitionLength);
		waldpdf(x, m_b, s_b, z, partitionLength);
		*/
		double * w[2];

		w[0] = y;
		w[1] = z;

#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
		for (int i = 0; i < 2; i++) {
			waldpdf(x, m[i], s[i], w[i], partitionLength);
		}

		//printf("\n\ncalling approxconv_rep from conv2waldpdf partitionLength=%d size_XY=%d h=%f\n", partitionLength, size_XY, h);
		//approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength, 266, h);
		approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);

		_mm_free(x); _mm_free(y); _mm_free(z);
		//free(x); free(y); free(z);

		while (E >= 0.001*fabs(logP0)) {
			h = h * 0.5; // Shrink the step size
#ifdef _VERBOSE
			printf("h=%f logP0=%f ", h, logP0);
#endif
			int partitionLength = (int)(maxX / h);
			double * x = _mm_malloc(partitionLength * sizeof(double), 32);
			double * y = _mm_malloc(partitionLength * sizeof(double), 32);
			double * z = _mm_malloc(partitionLength * sizeof(double), 32);
			//double * x = malloc(partitionLength * sizeof(double));
			//double * y = malloc(partitionLength * sizeof(double));
			//double * z = malloc(partitionLength * sizeof(double));

			// fill the partition
			/* CAN I REPLACE THIS WITH AN ADDITIVE FILLER? x[i+1] = x[i]+h
			HOTSPOT ANALYSIS SUGGESTS THIS IS 18% OF RUNTIME */

			/*
			for (int i = 0; i < partitionLength; i++)
			x[i] = i*h;
			*/
			double tally = 0;
			for (int i = 0; i < partitionLength; i++) {
				x[i] = tally;
				tally += h;
			}


			// evaluate the sub-distributions at each point in x
			/*
			waldpdf(x, m_a, s_a, y, partitionLength);
			waldpdf(x, m_b, s_b, z, partitionLength);
			*/
			w[0] = y;
			w[1] = z;

#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
			for (int i = 0; i < 2; i++) {
				waldpdf(x, m[i], s[i], w[i], partitionLength);
			}


			//printf("\n\ncalling approxconv_rep from conv2waldpdf partitionLength=%d size_XY=%d h=%f\n", partitionLength, size_XY, h);
			//approxconvolv_replacement(z, y, X, x, Y, &logP1, partitionLength, 266, h);
			approxconvolv_replacement(z, y, X, x, Y, &logP1, partitionLength, size_XY, h);

			_mm_free(x); _mm_free(y); _mm_free(z);
			//free(x); free(y); free(z);

			E = fabs(logP1 - logP0);

			logP0 = logP1;

			//printf("logP1=%f E=%f\n", logP1, E);

			//exit(1);
		}
	}

}


