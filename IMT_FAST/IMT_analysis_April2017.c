
/* Include files */
#include "IMT_analysis_April2017.h"
#include "emgpdf.h"
#include "onestagepdf2.h"
#include "onestagepdf_lag.h"
#include "convolv_2invG_adapt_nov.h"
#include "convolv_3invG_nov.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics_double.h"
#include "time.h"
#include "optimze_twostage.h"
#include "main.h"
#include "loglikelihood.h"

#define _VERBOSE
#define _CONV2WALD
//#define _PARALLEL_SEEDS
#define _ENABLE_EMG
#define _ENABLE_ONESTAGE
#define _ENABLE_ONESTAGELAG
#define _ENABLE_TWOSTAGE
#define _ENABLE_THREESTAGE

/* Type Definitions */
#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3

typedef struct {
    double f1[2];
} cell_wrap_3;

#endif				/*typedef_cell_wrap_3 */

/* Function Definitions */

void IMT_analysis_April2017(const char *model)
{

    printf("Using GSL minimizer\n");

	int data_size = 266;
	static const double data[266] =
	{ 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
		12.4, 13.7, 12.4, 11.9, 10.3, 12.9, 14.7, 11.6, 13.4, 13.4, 11.6,
		10.3,
		9.3,
		13.7, 9.6, 10.1, 9.8, 10.9, 16.0, 9.3, 9.6, 10.3, 11.4, 10.6, 8.5,
		10.3,
		11.1, 8.0, 10.6, 7.5, 12.9, 9.0, 8.5, 12.4, 11.6, 9.6, 9.6, 14.7,
		9.8,
		10.3,
		12.1, 8.8, 10.6, 12.1, 13.4, 12.4, 8.8, 13.2, 10.1, 11.6, 11.1,
		15.8,
		12.1,
		12.7, 12.7, 11.1, 13.2, 11.9, 12.4, 13.2, 14.0, 8.0, 8.8, 9.3,
		16.5, 14.5,
		10.1, 14.2, 7.8, 13.2, 8.8, 8.8, 10.1, 11.9, 12.9, 14.5, 10.9,
		10.6, 14.0,
		8.8, 8.8, 9.0, 10.9, 14.5, 9.6, 12.4, 11.9, 12.4, 11.1, 14.5, 10.3,
		12.4,
		12.7, 11.9, 10.3, 13.7, 15.5, 14.5, 11.6, 10.6, 15.5, 14.7, 8.8,
		11.6,
		8.3,
		17.6, 12.4, 11.6, 15.0, 13.7, 12.7, 10.9, 7.2, 8.5, 8.3, 9.6, 11.4,
		12.9,
		11.6, 13.4, 10.1, 11.6, 8.8, 12.4, 10.3, 16.3, 10.9, 10.1, 8.8,
		9.3, 15.2,
		8.5, 11.1, 8.3, 11.4, 11.9, 9.3, 9.8, 16.3, 12.7, 9.0, 11.9, 9.3,
		10.3,
		13.4,
		11.4, 12.9, 12.4, 9.6, 10.3, 13.2, 10.6, 9.8, 11.9, 14.2, 13.4,
		9.3, 9.6,
		12.1, 11.9, 10.1, 14.0, 12.9, 21.7, 11.6, 12.1, 10.3, 9.8, 14.2,
		13.7,
		7.2,
		10.9, 10.1, 9.6, 13.4, 13.2, 16.3, 11.6, 14.0, 10.9, 14.2, 12.4,
		12.4,
		13.4,
		17.6, 10.1, 10.9, 14.0, 12.9, 9.0, 13.4, 15.0, 16.0, 8.0, 9.8,
		12.4, 8.5,
		9.6, 12.7, 12.1, 15.0, 16.0, 10.9, 14.2, 13.7, 11.9, 16.8, 11.4,
		13.4,
		12.4,
		22.0, 12.4, 16.8, 12.1, 10.3, 13.4, 11.6, 10.1, 14.5, 10.6, 11.9,
		15.5,
		9.8,
		12.4, 10.1, 8.0, 9.0, 9.3, 13.2, 11.1, 12.7, 12.1, 10.1, 13.2,
		14.5, 10.1,
		12.7, 12.9, 11.9, 12.4, 11.1, 8.5, 14.5, 16.5, 12.4, 9.0, 11.1,
		9.8, 11.1,
		11.1, 8.8, 13.2, 17.6, 16.8, 10.9, 12.4, 8.5, 14.7
	};

    int twostagefitnomle;
    int onestagelagnomle;
    int onestagefitnomle;
    int emgfitnomle;
    int threestagefitnomle;
    int itmp;
    //double P[81];
    double mtmp;
    //int ix;
    //double p[3];
    //double ep[81];
    
    //double b_P[32];
    double ld[16];
    //double b_p[2];
    //double pd[32];

    cell_wrap_3 pcell[9];
    double c_P[180];
     /**/


    double b_pd[180];
    //double flag;
    //double E;
    double b_flag[45];
    double c_p[4];
    cell_wrap_3 b_pcell[4];
    static const double dv61[2] = { 0.8494, 0.0843 };

    double d_P[120];
    static const signed char b_id[60] =
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
	2, 2, 3, 3, 3, 4, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 2, 2, 2, 3, 3, 4,
	    3, 3, 4,
	4,
	1, 2, 3, 4, 2, 3, 4, 3, 4, 4, 2, 3, 4, 3, 4, 4, 3, 4, 4, 4
    };

    double c_pd[120];
    double c_flag[20];
    double d_p[6];

	distType l[266];
	

    printf("FUCCI Data\n\n");

    /* choose model to fit */
    twostagefitnomle = 0;
    const char *twostagefitnomle_str = "twostage";
    onestagelagnomle = 0;
    const char *onestagelagnomle_str = "onestagelag";
    onestagefitnomle = 0;
    const char *onestagefitnomle_str = "onestage";
    emgfitnomle = 0;
    const char *emgfitnomle_str = "emg";
    threestagefitnomle = 0;
    const char *threestagefitnomle_str = "threestage";
    const char *all_str = "all";

    if (strcmp(model, twostagefitnomle_str) == 0) {
	twostagefitnomle = 1;
    } else if (strcmp(model, onestagelagnomle_str) == 0) {
	onestagelagnomle = 1;
    } else if (strcmp(model, onestagefitnomle_str) == 0) {
	onestagefitnomle = 1;
    } else if (strcmp(model, emgfitnomle_str) == 0) {
	emgfitnomle = 1;
    } else if (strcmp(model, threestagefitnomle_str) == 0) {
	threestagefitnomle = 1;
    } else if (strcmp(model, all_str) == 0) {
	emgfitnomle = 1;
	onestagefitnomle = 1;
	onestagelagnomle = 1;
	twostagefitnomle = 1;
	threestagefitnomle = 1;
    } else {

	printf("No valid model was selected, quitting.\n");
    }

	if (emgfitnomle == 1) {
		optimize_emg(data, data_size);
	}

	if (onestagefitnomle == 1) {
		optimize_onestage(data, data_size);
	}

#ifdef _ENABLE_ONESTAGELAG
    /* Fit one-stage model with lag with fminsearch */
    if (onestagelagnomle == 1) {

	printf("onestagefitlagnomle\n\n");

	/* prepare statistical variables */
	double mean = gsl_stats_mean(data, 1, data_size);
	double variance = gsl_stats_variance(data, 1, data_size);

	// C3 = sum((data-C1).^3)/(length(data));
	double c3sum = 0;
	for (int i = 0; i < data_size; i++) {
	    c3sum += pow(data[i] - mean, 3.0);
	}
	c3sum = c3sum / data_size;

	double mu = c3sum / (3 * pow(variance, 2.0));
	double sigma =
	    pow((pow(c3sum, 3.0) / (27 * pow(variance, 5.0))), 0.5);
	double lag = mean - 3 * pow(variance, 2.0) / c3sum;

	double vryv[3] = { 0.5, 1, 2 };
	double vrym[3] = { 0.25, 0.5, 0.75 };

	double m[3];
	for (int i = 0; i < 3; i++)
	    m[i] = mu * vrym[i];

	double s[3];
	for (int i = 0; i < 3; i++)
	    s[i] = sigma * vryv[i];

	double L[3];
	for (int i = 0; i < 3; i++)
	    L[i] = lag * vrym[i];

	/* prepare parameter seeds */
	int seedIdx = 0;
	double paramSeeds[64][3];
	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
		    paramSeeds[seedIdx][0] = m[i];
		    paramSeeds[seedIdx][1] = s[j];
		    paramSeeds[seedIdx][2] = L[k];
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

	    printf("i=%d\n", 1 + seedIdx);
	    printf("  P=[%f %f %f]\n", paramSeeds[seedIdx][0],
		   paramSeeds[seedIdx][1], paramSeeds[seedIdx][2]);

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
	    gsl_vector_set(x, 0, paramSeeds[seedIdx][0]);
	    gsl_vector_set(x, 1, paramSeeds[seedIdx][1]);
	    gsl_vector_set(x, 2, paramSeeds[seedIdx][2]);

	    /* Set initial step sizes to 1 */
	    ss = gsl_vector_alloc(3);
	    gsl_vector_set_all(ss, 1.0);

	    /* Initialize method and iterate */
	    minex_func.n = 3;
	    minex_func.f = waldlag_loglikelihood;
	    minex_func.params = (void *) data;

	    s = gsl_multimin_fminimizer_alloc(T, 3);
	    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	    do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
		    break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
		if (status == GSL_SUCCESS) {
		    printf("converged to minimum at\n");
		}

		printf("%5d %.17f %.17f %.17f f() = %.17f size = %.3f\n",
		       (int) iter,
		       gsl_vector_get(s->x, 0),
		       gsl_vector_get(s->x, 1),
		       gsl_vector_get(s->x, 2), s->fval, size);
#endif
	    }
	    while (status == GSL_CONTINUE && iter < 10000);

	    optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
	    optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
	    optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));

	    gsl_vector_free(x);
	    gsl_vector_free(ss);
	    gsl_multimin_fminimizer_free(s);

	    printf("  p=[%f %f %f]\n", optimizedParams[seedIdx][0],
		   optimizedParams[seedIdx][1],
		   optimizedParams[seedIdx][2]);

	    waldlagpdf(data, optimizedParams[seedIdx][0],
		       optimizedParams[seedIdx][1],
		       optimizedParams[seedIdx][2], l, data_size);

	    /* calculate the log likelihood for our best fit */
		/*
	    double loglikelihood = 0;
	    for (int i = 0; i < data_size; i++) {
		loglikelihood += log(l[i]);
	    }
		*/

		double ll = (double) loglikelihood(l, data_size);

	    printf("  l=%.17f\n\n", ll);
	    le[seedIdx] = ll;
	}

	/* find the best optimized parameter set for all starting seeds tried */
	double maxLikelihood = 0;
	int ind_ld = 0;
	for (int i = 0; i < 27; i++) {
	    if (le[i] < maxLikelihood) {
		maxLikelihood = le[i];
		ind_ld = i;
	    }
	}

	printf("max_ld=%f row_ld=%d\n", maxLikelihood, ind_ld + 1);
	printf("pd_max=[%f %f %f]\n\n", optimizedParams[ind_ld][0],
	       optimizedParams[ind_ld][1], optimizedParams[ind_ld][2]);
    }
#endif

#ifdef _ENABLE_TWOSTAGE
	if (twostagefitnomle) {
		int numseeds_twostage = 45;
		double seeds_twostage[45][4] = {
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

		optimize_twostage(data_size, data, numseeds_twostage, seeds_twostage);
	}
#endif

#ifdef _ENABLE_THREESTAGE

    if (threestagefitnomle == 1) {

	printf("threestagefitnomle\n");

	for (itmp = 0; itmp < 2; itmp++) {
	    b_pcell[0].f1[itmp] =
		0.8494 + -0.25960000000000005 * (double) itmp;
	    b_pcell[1].f1[itmp] = dv61[itmp];
	    b_pcell[2].f1[itmp] =
		0.1213 + 0.46849999999999997 * (double) itmp;
	    b_pcell[3].f1[itmp] =
		0.1213 + -0.037000000000000005 * (double) itmp;
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

	/* options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs','Display','iter'); */
	/* 'IMT_analysis_April2017:858' options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs'); */
	/*  WARNING WARNING WARNING */
	/*  I've skipped over P(18),P(19) and P(20) for the moment b/c P(18) */
	/*  is pathological. this needs fixed!!!! */
	/*  WARNING WARNING WARNING */
	/* 'IMT_analysis_April2017:863' for i=1:length(P)-3 */
#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (emgfitnomle = 0; emgfitnomle < 17; emgfitnomle++) {

	    printf("i=%f\n", 1.0 + (double) emgfitnomle);

	    printf("  P=[%f %f %f %f %f %f]\n", d_P[emgfitnomle],
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
	    minex_func.params = (void *) data;

	    s = gsl_multimin_fminimizer_alloc(T, 6);
	    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	    do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
		    break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
		if (status == GSL_SUCCESS) {
		    printf("converged to minimum at\n");
		}

		printf
		    ("%5d %.8f %.8f %.8f %.8f %.8f %.8f f() = %.17f size = %.3f\n",
		     (int) iter, gsl_vector_get(s->x, 0),
		     gsl_vector_get(s->x, 1), gsl_vector_get(s->x, 2),
		     gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4),
		     gsl_vector_get(s->x, 5), s->fval, size);
#endif
	    }
	    while (status == GSL_CONTINUE && iter < 10000);

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

	    convolv3waldpdf(d_p[0], d_p[1], d_p[2], d_p[3], d_p[4], d_p[5],
			    data, l, data_size, 0.01);

		/*
	    double l_sum = 0;
	    for (int i = 0; i < data_size; i++) {
		l_sum += log(l[i]);
	    }
		*/

		double l_sum = (double) loglikelihood(l, data_size);

		/* FIXME FIXME FIXME we are not reporting the real values here */
	    printf("  l=%f hp=%f flag=%f E=%f WARNING HP FLAG AND E ARE BOGUS\n\n", l_sum, 9999.0, 9999.0, 9999.0);

	}

	printf("Recalculating canidate solutions with smaller stepsize\n");

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
	for (emgfitnomle = 0; emgfitnomle < 20; emgfitnomle++) {

	    convolv3waldpdf(c_pd[emgfitnomle], c_pd[20 + emgfitnomle],
			    c_pd[40 + emgfitnomle], c_pd[60 + emgfitnomle],
			    c_pd[80 + emgfitnomle],
			    c_pd[100 + emgfitnomle], data, l, data_size, 0.001);

		/*
	    double l_sum = 0;
	    for (int i = 0; i < data_size; i++) {
		l_sum += log(l[i]);
	    }
		*/

		double l_sum = (double) loglikelihood(l, data_size);

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

	printf("max_ld=%f row_ld=%f\n", mtmp, (double) itmp + 1);

    }
#endif

}
