

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

typedef struct
{
  double f1[2];
} cell_wrap_3;

#endif /*typedef_cell_wrap_3 */

/* Function Definitions */


void
IMT_analysis_April2017 (const char *model)
{


  printf ("Using GSL minimizer\n");


  int twostagefitnomle;
  int onestagelagnomle;
  int onestagefitnomle;
  int emgfitnomle;
  int threestagefitnomle;
  int itmp;
  double P[81];
  double mtmp;
  int ix;
  double p[3];
  double ep[81];
  double l[266];
  double b_P[32];
  double ld[16];
  double b_p[2];
  double pd[32];


  cell_wrap_3 pcell[9];
  double c_P[180];
   /**/
    static const signed char id[90] =
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6,
    7,
    7, 7, 8, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5,
    6,
    7, 8, 9, 4, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9, 6, 7, 8, 9, 7, 8, 9, 8, 9, 9
  };

  double b_pd[180];
  double flag;
  double E;
  double b_flag[45];
  double c_p[4];
  cell_wrap_3 b_pcell[4];
  static const double dv61[2] = { 0.8494, 0.0843 };

  double d_P[120];
  static const signed char b_id[60] =
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
    2, 2, 3, 3, 3, 4, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 2, 2, 2, 3, 3, 4, 3, 3, 4,
    4,
    1, 2, 3, 4, 2, 3, 4, 3, 4, 4, 2, 3, 4, 3, 4, 4, 3, 4, 4, 4
  };

  double c_pd[120];
  double c_flag[20];
  double d_p[6];


  static const double data[266] =
    { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
    12.4, 13.7, 12.4, 11.9, 10.3, 12.9, 14.7, 11.6, 13.4, 13.4, 11.6, 10.3,
    9.3,
    13.7, 9.6, 10.1, 9.8, 10.9, 16.0, 9.3, 9.6, 10.3, 11.4, 10.6, 8.5, 10.3,
    11.1, 8.0, 10.6, 7.5, 12.9, 9.0, 8.5, 12.4, 11.6, 9.6, 9.6, 14.7, 9.8,
    10.3,
    12.1, 8.8, 10.6, 12.1, 13.4, 12.4, 8.8, 13.2, 10.1, 11.6, 11.1, 15.8,
    12.1,
    12.7, 12.7, 11.1, 13.2, 11.9, 12.4, 13.2, 14.0, 8.0, 8.8, 9.3, 16.5, 14.5,
    10.1, 14.2, 7.8, 13.2, 8.8, 8.8, 10.1, 11.9, 12.9, 14.5, 10.9, 10.6, 14.0,
    8.8, 8.8, 9.0, 10.9, 14.5, 9.6, 12.4, 11.9, 12.4, 11.1, 14.5, 10.3, 12.4,
    12.7, 11.9, 10.3, 13.7, 15.5, 14.5, 11.6, 10.6, 15.5, 14.7, 8.8, 11.6,
    8.3,
    17.6, 12.4, 11.6, 15.0, 13.7, 12.7, 10.9, 7.2, 8.5, 8.3, 9.6, 11.4, 12.9,
    11.6, 13.4, 10.1, 11.6, 8.8, 12.4, 10.3, 16.3, 10.9, 10.1, 8.8, 9.3, 15.2,
    8.5, 11.1, 8.3, 11.4, 11.9, 9.3, 9.8, 16.3, 12.7, 9.0, 11.9, 9.3, 10.3,
    13.4,
    11.4, 12.9, 12.4, 9.6, 10.3, 13.2, 10.6, 9.8, 11.9, 14.2, 13.4, 9.3, 9.6,
    12.1, 11.9, 10.1, 14.0, 12.9, 21.7, 11.6, 12.1, 10.3, 9.8, 14.2, 13.7,
    7.2,
    10.9, 10.1, 9.6, 13.4, 13.2, 16.3, 11.6, 14.0, 10.9, 14.2, 12.4, 12.4,
    13.4,
    17.6, 10.1, 10.9, 14.0, 12.9, 9.0, 13.4, 15.0, 16.0, 8.0, 9.8, 12.4, 8.5,
    9.6, 12.7, 12.1, 15.0, 16.0, 10.9, 14.2, 13.7, 11.9, 16.8, 11.4, 13.4,
    12.4,
    22.0, 12.4, 16.8, 12.1, 10.3, 13.4, 11.6, 10.1, 14.5, 10.6, 11.9, 15.5,
    9.8,
    12.4, 10.1, 8.0, 9.0, 9.3, 13.2, 11.1, 12.7, 12.1, 10.1, 13.2, 14.5, 10.1,
    12.7, 12.9, 11.9, 12.4, 11.1, 8.5, 14.5, 16.5, 12.4, 9.0, 11.1, 9.8, 11.1,
    11.1, 8.8, 13.2, 17.6, 16.8, 10.9, 12.4, 8.5, 14.7
  };






  printf ("FUCCI Data\n\n");


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

  if (strcmp (model, twostagefitnomle_str) == 0)
    {
      twostagefitnomle = 1;
    }
  else if (strcmp (model, onestagelagnomle_str) == 0)
    {
      onestagelagnomle = 1;
    }
  else if (strcmp (model, onestagefitnomle_str) == 0)
    {
      onestagefitnomle = 1;
    }
  else if (strcmp (model, emgfitnomle_str) == 0)
    {
      emgfitnomle = 1;
    }
  else if (strcmp (model, threestagefitnomle_str) == 0)
    {
      threestagefitnomle = 1;
    }
  else if (strcmp (model, all_str) == 0)
    {
      emgfitnomle = 1;
      onestagefitnomle = 1;
      onestagelagnomle = 1;
      twostagefitnomle = 1;
      threestagefitnomle = 1;
    }
  else
    {

      printf ("No valid model was selected, quitting.\n");
    }


#ifdef _ENABLE_EMG
  /* Fit the EMG model */
  if (emgfitnomle == 1)
    {

      printf ("emgfitnomle\n\n");

      /* prepare statistical variables */
      double mean = gsl_stats_mean (data, 1, 266);
      double variance = gsl_stats_variance (data, 1, 266);
      double vry[3] = { 0.25, 0.5, 0.75 };

      /*we vary the parameters so the the Gaussian and exponential parts of
         the cell cycle are responsible for a fraction of the total mean and
         variance in the IMT.
       */
      double lambda[3];
      for (int i = 0; i < 3; i++)
	lambda[i] = 1.0 / (mean * vry[i]);

      double mu[3];
      for (int i = 0; i < 3; i++)
	mu[i] = mean * vry[i];

      double sigma[3];
      for (int i = 0; i < 3; i++)
	sigma[i] = pow (variance * vry[i], 0.5);

      /* prepare parameter seeds */
      int seedIdx = 0;
      double paramSeeds[27][3];
      for (int i = 0; i < 3; i++)
	{
	  for (int j = 0; j < 3; j++)
	    {
	      for (int k = 0; k < 3; k++)
		{
		  paramSeeds[seedIdx][0] = lambda[i];
		  paramSeeds[seedIdx][1] = sigma[j];
		  paramSeeds[seedIdx][2] = mu[k];
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
      for (seedIdx = 0; seedIdx < 27; seedIdx++)
	{

	  printf ("i = %d\n", 1 + seedIdx);
	  printf ("  P=[%.17f %.17f %.17f]\n", paramSeeds[seedIdx][0],
		  paramSeeds[seedIdx][1], paramSeeds[seedIdx][2]);

	  /* select our working seed */
	  double paramSeed[3];
	  for (int i = 0; i < 3; i++)
	    paramSeed[i] = paramSeeds[seedIdx][i];

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
	  x = gsl_vector_alloc (3);
	  gsl_vector_set (x, 0, paramSeed[0]);
	  gsl_vector_set (x, 1, paramSeed[1]);
	  gsl_vector_set (x, 2, paramSeed[2]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc (3);
	  gsl_vector_set_all (ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 3;
	  minex_func.f = emgpdf_loglikelihood;
	  minex_func.params = (void *) data;

	  s = gsl_multimin_fminimizer_alloc (T, 3);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  do
	    {
	      iter++;
	      status = gsl_multimin_fminimizer_iterate (s);

	      if (status)
		break;

	      size = gsl_multimin_fminimizer_size (s);
	      status = gsl_multimin_test_size (size, 1e-4);

#ifdef _VERBOSE
	      if (status == GSL_SUCCESS)
		{
		  printf ("converged to minimum at\n");
		}

	      printf ("%5d %.17f %.17f %.17f f() = %.17f size = %.3f\n",
		      (int) iter,
		      gsl_vector_get (s->x, 0),
		      gsl_vector_get (s->x, 1),
		      gsl_vector_get (s->x, 2), s->fval, size);
#endif
	    }
	  while (status == GSL_CONTINUE && iter < 10000);

	  optimizedParams[seedIdx][0] = fabs (gsl_vector_get (s->x, 0));
	  optimizedParams[seedIdx][1] = fabs (gsl_vector_get (s->x, 1));
	  optimizedParams[seedIdx][2] = fabs (gsl_vector_get (s->x, 2));

	  gsl_vector_free (x);
	  gsl_vector_free (ss);
	  gsl_multimin_fminimizer_free (s);




	  printf ("  ep=[%.17f %.17f %.17f]\n", optimizedParams[seedIdx][0],
		  optimizedParams[seedIdx][1], optimizedParams[seedIdx][2]);

	  /* recalculate our best fit */
	  emgpdf (data, optimizedParams[seedIdx][0],
		  optimizedParams[seedIdx][1], optimizedParams[seedIdx][2],
		  l);

	  /* calculate the log likelihood for our best fit */
	  double loglikelihood = 0;
	  for (int i = 0; i < 266; i++)
	    {
	      loglikelihood += log (l[i]);
	    }

	  printf ("  l=%.17f\n\n", loglikelihood);
	  le[seedIdx] = loglikelihood;
	}

      /* find the best optimized parameter set for all starting seeds tried */
      double maxLikelihood = 0;
      int ind_le = 0;
      for (int i = 0; i < 27; i++)
	{
	  if (le[i] < maxLikelihood)
	    {
	      maxLikelihood = le[i];
	      ind_le = i;
	    }
	}

      printf ("max_le=%.17f ind_le=%d\n", maxLikelihood, ind_le + 1);
      printf ("ep_max=[%.17f %.17f %.17f]\n\n", optimizedParams[ind_le][0],
	      optimizedParams[ind_le][1], optimizedParams[ind_le][2]);
    }
#endif

#ifdef _ENABLE_ONESTAGE


  if (onestagefitnomle == 1)
    {

      printf ("onestagefitnomle\n\n");

      /* prepare statistical variables */
      double mean = gsl_stats_mean (data, 1, 266);
      double variance = gsl_stats_variance (data, 1, 266);
      double vry[4] = { 0.5, 1.0, 1.5, 2.0 };

      double mu = 1.0 / mean;
      double sigma = pow ((variance / pow (mean, 3)), 0.5);

      double m[4];
      for (int i = 0; i < 4; i++)
	m[i] = mu * vry[i];

      double s[4];
      for (int i = 0; i < 4; i++)
	s[i] = sigma * vry[i];

      /* prepare parameter seeds */
      int seedIdx = 0;
      double paramSeeds[16][2];
      for (int i = 0; i < 4; i++)
	{
	  for (int j = 0; j < 4; j++)
	    {
	      paramSeeds[seedIdx][0] = m[i];
	      paramSeeds[seedIdx][1] = s[j];
	      seedIdx++;
	    }
	}

      /* optimize parameters */
      double optimizedParams[27][3];
      double le[27];
#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
      for (seedIdx = 0; seedIdx < 16; seedIdx++)
	{
	  printf ("i=%d\n", 1 + seedIdx);
	  printf ("  P=[%f %f]\n", paramSeeds[seedIdx][0],
		  paramSeeds[seedIdx][1]);



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
	  x = gsl_vector_alloc (2);
	  gsl_vector_set (x, 0, paramSeeds[seedIdx][0]);
	  gsl_vector_set (x, 1, paramSeeds[seedIdx][1]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc (2);
	  gsl_vector_set_all (ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 2;
	  minex_func.f = wald_loglikelihood;
	  minex_func.params = (void *) data;

	  s = gsl_multimin_fminimizer_alloc (T, 2);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  do
	    {
	      iter++;
	      status = gsl_multimin_fminimizer_iterate (s);

	      if (status)
		break;

	      size = gsl_multimin_fminimizer_size (s);
	      status = gsl_multimin_test_size (size, 1e-4);

#ifdef _VERBOSE
	      if (status == GSL_SUCCESS)
		{
		  printf ("converged to minimum at\n");
		}

	      printf ("%5d %.17f %.17f f() = %.17f size = %.3f\n",
		      (int) iter,
		      gsl_vector_get (s->x, 0),
		      gsl_vector_get (s->x, 1), s->fval, size);
#endif
	    }
	  while (status == GSL_CONTINUE && iter < 10000);

	  optimizedParams[seedIdx][0] = fabs (gsl_vector_get (s->x, 0));
	  optimizedParams[seedIdx][1] = fabs (gsl_vector_get (s->x, 1));

	  gsl_vector_free (x);
	  gsl_vector_free (ss);
	  gsl_multimin_fminimizer_free (s);



	  printf ("  p=[%f %f]\n", optimizedParams[seedIdx][0],
		  optimizedParams[seedIdx][1]);

	  //onestagepdf2(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], l);
	  waldpdf (data, optimizedParams[seedIdx][0],
		   optimizedParams[seedIdx][1], l, 266);


	  /* calculate the log likelihood for our best fit */
	  double loglikelihood = 0;
	  for (int i = 0; i < 266; i++)
	    {
	      loglikelihood += log (l[i]);
	    }

	  printf ("  l=%.17f\n\n", loglikelihood);
	  ld[seedIdx] = loglikelihood;
	}

      /* find the best optimized parameter set for all starting seeds tried */
      double maxLikelihood = 0;
      int ind_ld = 0;
      for (int i = 0; i < 16; i++)
	{
	  if (ld[i] < maxLikelihood)
	    {
	      maxLikelihood = ld[i];
	      ind_ld = i;
	    }
	}

      printf ("max_ld=%f row_ld=%d\n", maxLikelihood, ind_ld);
      printf ("pd_max=[%f %f]\n\n", optimizedParams[ind_ld][0],
	      optimizedParams[ind_ld][1]);
    }
#endif

#ifdef _ENABLE_ONESTAGELAG
  /* Fit one-stage model with lag with fminsearch */
  if (onestagelagnomle == 1)
    {

      printf ("onestagefitlagnomle\n\n");

      /* prepare statistical variables */
      double mean = gsl_stats_mean (data, 1, 266);
      double variance = gsl_stats_variance (data, 1, 266);

      // C3 = sum((data-C1).^3)/(length(data));
      double c3sum = 0;
      for (int i = 0; i < 266; i++)
	{
	  c3sum += pow (data[i] - mean, 3.0);
	}
      c3sum = c3sum / 266;

      double mu = c3sum / (3 * pow (variance, 2.0));
      double sigma =
	pow ((pow (c3sum, 3.0) / (27 * pow (variance, 5.0))), 0.5);
      double lag = mean - 3 * pow (variance, 2.0) / c3sum;

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
      for (int i = 0; i < 3; i++)
	{
	  for (int j = 0; j < 3; j++)
	    {
	      for (int k = 0; k < 3; k++)
		{
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
      for (seedIdx = 0; seedIdx < 27; seedIdx++)
	{

	  printf ("i=%d\n", 1 + seedIdx);
	  printf ("  P=[%f %f %f]\n", paramSeeds[seedIdx][0],
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
	  x = gsl_vector_alloc (3);
	  gsl_vector_set (x, 0, paramSeeds[seedIdx][0]);
	  gsl_vector_set (x, 1, paramSeeds[seedIdx][1]);
	  gsl_vector_set (x, 2, paramSeeds[seedIdx][2]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc (3);
	  gsl_vector_set_all (ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 3;
	  minex_func.f = waldlag_loglikelihood;
	  minex_func.params = (void *) data;

	  s = gsl_multimin_fminimizer_alloc (T, 3);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  do
	    {
	      iter++;
	      status = gsl_multimin_fminimizer_iterate (s);

	      if (status)
		break;

	      size = gsl_multimin_fminimizer_size (s);
	      status = gsl_multimin_test_size (size, 1e-4);

#ifdef _VERBOSE
	      if (status == GSL_SUCCESS)
		{
		  printf ("converged to minimum at\n");
		}

	      printf ("%5d %.17f %.17f %.17f f() = %.17f size = %.3f\n",
		      (int) iter,
		      gsl_vector_get (s->x, 0),
		      gsl_vector_get (s->x, 1),
		      gsl_vector_get (s->x, 2), s->fval, size);
#endif
	    }
	  while (status == GSL_CONTINUE && iter < 10000);

	  optimizedParams[seedIdx][0] = fabs (gsl_vector_get (s->x, 0));
	  optimizedParams[seedIdx][1] = fabs (gsl_vector_get (s->x, 1));
	  optimizedParams[seedIdx][2] = fabs (gsl_vector_get (s->x, 2));

	  gsl_vector_free (x);
	  gsl_vector_free (ss);
	  gsl_multimin_fminimizer_free (s);



	  printf ("  p=[%f %f %f]\n", optimizedParams[seedIdx][0],
		  optimizedParams[seedIdx][1], optimizedParams[seedIdx][2]);

	  waldlagpdf (data, optimizedParams[seedIdx][0],
		      optimizedParams[seedIdx][1],
		      optimizedParams[seedIdx][2], l, 266);


	  /* calculate the log likelihood for our best fit */
	  double loglikelihood = 0;
	  for (int i = 0; i < 266; i++)
	    {
	      loglikelihood += log (l[i]);
	    }

	  printf ("  l=%.17f\n\n", loglikelihood);
	  le[seedIdx] = loglikelihood;
	}

      /* find the best optimized parameter set for all starting seeds tried */
      double maxLikelihood = 0;
      int ind_ld = 0;
      for (int i = 0; i < 27; i++)
	{
	  if (le[i] < maxLikelihood)
	    {
	      maxLikelihood = le[i];
	      ind_ld = i;
	    }
	}

      printf ("max_ld=%f row_ld=%d\n", maxLikelihood, ind_ld + 1);
      printf ("pd_max=[%f %f %f]\n\n", optimizedParams[ind_ld][0],
	      optimizedParams[ind_ld][1], optimizedParams[ind_ld][2]);
    }
#endif

#ifdef _ENABLE_TWOSTAGE

  if (twostagefitnomle == 1)
    {

      printf ("twostagefitnomle\n");

      /* prepare statistical variables */
      double mean = gsl_stats_mean (data, 1, 266);
      double variance = gsl_stats_variance (data, 1, 266);

      double vry[3] = { 0.25, 0.5, 0.75 };
      double vrys[3] = { 0.01, 1, 10 };

      double m[3];
      double s[3];
      for (int i = 0; i < 3; i++)
	{
	  m[i] = 1 / (mean * vry[i]);
	  s[i] = pow ((variance * vry[i]) / (pow (mean * vry[i], 3.0)), 0.5);
	}

      double pcomb[9][2];
      int pcomb_idx = 0;
      for (int i = 0; i < 3; i++)
	{
	  for (int j = 0; j < 3; j++)
	    {
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

      for (itmp = 0; itmp < 2; itmp++)
	{
	  pcell[0].f1[itmp] = 0.3397 + -0.1038 * (double) itmp;
	  pcell[1].f1[itmp] = 0.3397 + -0.2217 * (double) itmp;
	  pcell[2].f1[itmp] = 0.3397 + -0.2611 * (double) itmp;
	  pcell[3].f1[itmp] = 0.1699 + 0.066 * (double) itmp;
	  pcell[4].f1[itmp] = 0.1699 + -0.0519 * (double) itmp;
	  pcell[5].f1[itmp] = 0.1699 + -0.091299999999999992 * (double) itmp;
	  pcell[6].f1[itmp] = 0.1132 + 0.1227 * (double) itmp;
	  pcell[7].f1[itmp] = 0.1132 + 0.0047999999999999987 * (double) itmp;
	  pcell[8].f1[itmp] = 0.1132 + -0.034599999999999992 * (double) itmp;
	}

      /* get all pairs of indices for the parameter  */
      /* choices for each part of the cycle to get all  */
      /* parameter choices for the entire cycle */

#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
      for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++)
	{
	  /* 'IMT_analysis_April2017:596' P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}]; */
	  for (itmp = 0; itmp < 2; itmp++)
	    {
	      c_P[emgfitnomle + 45 * itmp] =
		pcell[id[emgfitnomle] - 1].f1[itmp];
	      c_P[emgfitnomle + 45 * (itmp + 2)] =
		pcell[id[45 + emgfitnomle] - 1].f1[itmp];
	    }
	}


#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
      for (emgfitnomle = 0; emgfitnomle < 10; emgfitnomle++)
	{
	  //for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {


	  printf ("i=%f\n", 1.0 + (double) emgfitnomle);

	  printf ("  P=[%f %f %f %f]\n", c_P[emgfitnomle],
		  c_P[45 + emgfitnomle], c_P[90 + emgfitnomle],
		  c_P[135 + emgfitnomle]);


	  for (itmp = 0; itmp < 4; itmp++)
	    {
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
	  x = gsl_vector_alloc (4);
	  gsl_vector_set (x, 0, c_p[0]);
	  gsl_vector_set (x, 1, c_p[1]);
	  gsl_vector_set (x, 2, c_p[2]);
	  gsl_vector_set (x, 3, c_p[3]);


	  // Store previous state
	  gsl_vector *prev;
	  prev = gsl_vector_alloc (4);
	  gsl_vector_set_all (prev, 0);
	  int repeated = 0;


	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc (4);
	  gsl_vector_set_all (ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 4;
	  minex_func.f = convolv_2invG_adapt_nov_loglikelihood;
	  minex_func.params = (void *) data;

	  s = gsl_multimin_fminimizer_alloc (T, 4);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  do
	    {
	      iter++;

	      clock_t t;
	      t = clock ();
	      status = gsl_multimin_fminimizer_iterate (s);
	      t = clock () - t;
	      //printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);
	      printf ("%.3f ", ((float) t) / CLOCKS_PER_SEC);

	      if (status)
		break;

	      size = gsl_multimin_fminimizer_size (s);
	      status = gsl_multimin_test_size (size, 1e-4);

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
	      if (status == GSL_SUCCESS)
		{
		  printf ("converged to minimum at\n");
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
	    }
	  while (status == GSL_CONTINUE && iter < 10000);
	  printf ("\n");
	  c_p[0] = fabs (gsl_vector_get (s->x, 0));
	  c_p[1] = fabs (gsl_vector_get (s->x, 1));
	  c_p[2] = fabs (gsl_vector_get (s->x, 2));
	  c_p[3] = fabs (gsl_vector_get (s->x, 3));

	  gsl_vector_free (x);
	  gsl_vector_free (ss);
	  gsl_multimin_fminimizer_free (s);



	  printf ("  p=[%f %f %f %f]\n", c_p[0], c_p[1], c_p[2], c_p[3]);


	  for (itmp = 0; itmp < 4; itmp++)
	    {
	      b_pd[emgfitnomle + 45 * itmp] = c_p[itmp];
	    }




	  conv2waldpdf (data, c_p[0], c_p[1], c_p[2], c_p[3], l, 0.01, 1,
			266);



	  int l_sum = 0;
	  for (int i = 0; i < 266; i++)
	    {
	      l_sum += log (l[i]);
	    }

	  printf ("  l=%f hp=%f flag=%f E=%f\n\n", l_sum, 42, 42, 42);
	}

      /*  we previously optimized with a larger step size, recalculate with */
      /*  a smaller stepsize after the fact */
      printf ("recalculating canidate solutions with smaller stepsize\n");


#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
      for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++)
	{


	  conv2waldpdf (data, b_pd[emgfitnomle], b_pd[45 + emgfitnomle],
			b_pd[90 + emgfitnomle], b_pd[135 + emgfitnomle], l,
			0.001, 1, 266);



	  int l_sum = 0;
	  for (int i = 0; i < 266; i++)
	    l_sum += log (l[i]);

	  b_flag[emgfitnomle] = l_sum;
	}

      double max_ld = DBL_MIN;
      int row_id = 0;
      for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++)
	{
	  if (b_flag[emgfitnomle] > max_ld)
	    {
	      max_ld = b_flag[emgfitnomle];
	      row_id = emgfitnomle + 1;
	    }
	}

      printf ("max_ld=%f row_ld=%f\n", max_ld, row_id);
      printf ("pd_max=[%f %f %f %f]\n\n", b_pd[row_id - 1],
	      b_pd[45 + row_id - 1], b_pd[90 + row_id - 1],
	      b_pd[135 + row_id - 1]);


      /*  common to each fit, consider factoring out */
      emgfitnomle = 1;
      mtmp = b_flag[0];
      itmp = 0;


      if (emgfitnomle < 45)
	{
	  while (emgfitnomle + 1 < 46)
	    {
	      if (b_flag[emgfitnomle] > mtmp)
		{
		  mtmp = b_flag[emgfitnomle];
		  itmp = emgfitnomle;
		}

	      emgfitnomle++;
	    }
	}


      printf ("max_ld=%f row_ld=%f\n", mtmp, (double) itmp + 1);
      printf ("pd_max=[%f %f %f %f]\n\n", b_pd[itmp], b_pd[45 + itmp],
	      b_pd[90 + itmp], b_pd[135 + itmp]);

      /*  END FUNCTION FIT_TWOSTAGE */
    }
#endif

#ifdef _ENABLE_THREESTAGE

  if (threestagefitnomle == 1)
    {

      printf ("threestagefitnomle\n");


      for (itmp = 0; itmp < 2; itmp++)
	{
	  b_pcell[0].f1[itmp] = 0.8494 + -0.25960000000000005 * (double) itmp;
	  b_pcell[1].f1[itmp] = dv61[itmp];
	  b_pcell[2].f1[itmp] = 0.1213 + 0.46849999999999997 * (double) itmp;
	  b_pcell[3].f1[itmp] =
	    0.1213 + -0.037000000000000005 * (double) itmp;
	}


      for (emgfitnomle = 0; emgfitnomle < 20; emgfitnomle++)
	{
	  for (itmp = 0; itmp < 2; itmp++)
	    {
	      d_P[emgfitnomle + 20 * itmp] =
		b_pcell[b_id[emgfitnomle] - 1].f1[itmp];
	      d_P[emgfitnomle + 20 * (itmp + 2)] =
		b_pcell[b_id[20 + emgfitnomle] - 1].f1[itmp];
	      d_P[emgfitnomle + 20 * (itmp + 4)] =
		b_pcell[b_id[40 + emgfitnomle] - 1].f1[itmp];
	    }
	}


      memset (&c_pd[0], 0, 120U * sizeof (double));


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
      for (emgfitnomle = 0; emgfitnomle < 17; emgfitnomle++)
	{

	  printf ("i=%f\n", 1.0 + (double) emgfitnomle);


	  printf ("  P=[%f %f %f %f %f %f]\n", d_P[emgfitnomle],
		  d_P[20 + emgfitnomle], d_P[40 + emgfitnomle],
		  d_P[60 + emgfitnomle], d_P[80 + emgfitnomle],
		  d_P[100 + emgfitnomle]);


	  for (itmp = 0; itmp < 6; itmp++)
	    {
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
	  x = gsl_vector_alloc (6);
	  gsl_vector_set (x, 0, d_p[0]);
	  gsl_vector_set (x, 1, d_p[1]);
	  gsl_vector_set (x, 2, d_p[2]);
	  gsl_vector_set (x, 3, d_p[3]);
	  gsl_vector_set (x, 4, d_p[4]);
	  gsl_vector_set (x, 5, d_p[5]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc (6);
	  gsl_vector_set_all (ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 6;
	  minex_func.f = convolv_3invG_nov_loglikelihood;
	  minex_func.params = (void *) data;

	  s = gsl_multimin_fminimizer_alloc (T, 6);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  do
	    {
	      iter++;
	      status = gsl_multimin_fminimizer_iterate (s);

	      if (status)
		break;

	      size = gsl_multimin_fminimizer_size (s);
	      status = gsl_multimin_test_size (size, 1e-4);

#ifdef _VERBOSE
	      if (status == GSL_SUCCESS)
		{
		  printf ("converged to minimum at\n");
		}

	      printf
		("%5d %.8f %.8f %.8f %.8f %.8f %.8f f() = %.17f size = %.3f\n",
		 (int) iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x,
								       1),
		 gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3),
		 gsl_vector_get (s->x, 4), gsl_vector_get (s->x, 5), s->fval,
		 size);
#endif
	    }
	  while (status == GSL_CONTINUE && iter < 10000);

	  d_p[0] = fabs (gsl_vector_get (s->x, 0));
	  d_p[1] = fabs (gsl_vector_get (s->x, 1));
	  d_p[2] = fabs (gsl_vector_get (s->x, 2));
	  d_p[3] = fabs (gsl_vector_get (s->x, 3));
	  d_p[4] = fabs (gsl_vector_get (s->x, 4));
	  d_p[5] = fabs (gsl_vector_get (s->x, 5));

	  gsl_vector_free (x);
	  gsl_vector_free (ss);
	  gsl_multimin_fminimizer_free (s);



	  printf ("  p=[%f %f %f %f %f %f]\n", d_p[0], d_p[1], d_p[2], d_p[3],
		  d_p[4], d_p[5]);


	  for (itmp = 0; itmp < 6; itmp++)
	    {
	      c_pd[emgfitnomle + 20 * itmp] = d_p[itmp];
	    }

	  convolv3waldpdf (d_p[0], d_p[1], d_p[2], d_p[3], d_p[4], d_p[5],
			   data, l, 266, 0.01);



	  int l_sum = 0;
	  for (int i = 0; i < 266; i++)
	    {
	      l_sum += log (l[i]);
	    }

	  printf ("  l=%f hp=%f flag=%f E=%f\n\n", l_sum, 42, 42, 42);



	}


      printf ("Recalculating canidate solutions with smaller stepsize\n");


#ifdef _PARALLEL_SEEDS
#pragma omp parallel for
#endif
      for (emgfitnomle = 0; emgfitnomle < 20; emgfitnomle++)
	{

	  convolv3waldpdf (c_pd[emgfitnomle], c_pd[20 + emgfitnomle],
			   c_pd[40 + emgfitnomle], c_pd[60 + emgfitnomle],
			   c_pd[80 + emgfitnomle], c_pd[100 + emgfitnomle],
			   data, l, 266, 0.001);


	  int l_sum = 0;
	  for (int i = 0; i < 266; i++)
	    {
	      l_sum += log (l[i]);
	    }

	  c_flag[emgfitnomle] = l_sum;
	}


      emgfitnomle = 1;
      mtmp = c_flag[0];
      itmp = 0;


      if (emgfitnomle < 20)
	{
	  while (emgfitnomle + 1 < 21)
	    {
	      if (c_flag[emgfitnomle] > mtmp)
		{
		  mtmp = c_flag[emgfitnomle];
		  itmp = emgfitnomle;
		}

	      emgfitnomle++;
	    }
	}


      printf ("max_ld=%f row_ld=%f\n", mtmp, (double) itmp + 1);


    }

#endif

}
