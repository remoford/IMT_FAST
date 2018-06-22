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

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3

typedef struct {
	double f1[2];
} cell_wrap_3;

#endif

#define _GSL_GP_MAX_FIXED
//#define _PARALLEL_PDF
#define _BINNED_MODE

void optimize_threestage(const double data[], int data_size, configStruct config) {

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

	distType * l = (distType *)malloc(sizeof(distType)*data_size);
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

		printf("i=%f\n", 1.0 + (double)emgfitnomle);

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
		minex_func.params = (void *) &config;

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
				(int)iter, gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1), gsl_vector_get(s->x, 2),
				gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4),
				gsl_vector_get(s->x, 5), s->fval, size);
#endif
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

		convolv3waldpdf(d_p[0], d_p[1], d_p[2], d_p[3], d_p[4], d_p[5],
			data, l, data_size, 0.01);

		/*
		double l_sum = 0;
		for (int i = 0; i < data_size; i++) {
		l_sum += log(l[i]);
		}
		*/

		double l_sum = (double)loglikelihood(l, data_size);

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

	
	free(l);
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

	distType * Y = (distType *)malloc(sizeof(distType)*data_size);

    for (int i = 0; i < data_size; i++) {
		Y[i] = 0;
    }

    convolv3waldpdf(m1, s1, m2, s2, m3, s3, data, Y, data_size, 0.1);

	double ll = (double) loglikelihood(Y, data_size);

    double objective = penalty - ll;

	free(Y);

    return objective;
}

void
convolv3waldpdf(double m1, double s1, double m2, double s2, double m3,
		double s3, const double X[], distType Y[], int size_XY,
		double h)
{

    int flag = 0;

    double maxX = 0;
    for (int i = 0; i < size_XY; i++) {
	if (X[i] > maxX)
	    maxX = X[i];
    }

    int partitionLength = (int) (maxX / h);

    // This represents the partition
#ifdef __INTEL_COMPILER
    double *x =
	(double *) _mm_malloc(partitionLength * sizeof(double), 32);
#else
    double *x = malloc(partitionLength * sizeof(double));
#endif

    // fill the partition
    double tally = 0;
    for (int i = 0; i < partitionLength; i++) {
		x[i] = tally;
		tally += h;
    }

    double E = DBL_MAX;
    double eps = 0.01;

    double m[3];
    m[0] = m1;
    m[1] = m2;
    m[2] = m3;

    double s[3];
    s[0] = s1;
    s[1] = s2;
    s[2] = s3;

    m[0] = fabs(m[0]);
    m[1] = fabs(m[1]);
    m[2] = fabs(m[2]);

    s[0] = fabs(s[0]);
    s[1] = fabs(s[1]);
    s[2] = fabs(s[2]);

    // find the variance for both sub-distributions
    double v[3];
    v[0] = (s[0] * s[0]) / (m[0] * m[0] * m[0]);
    v[1] = (s[1] * s[1]) / (m[1] * m[1] * m[1]);
    v[2] = (s[2] * s[2]) / (m[2] * m[2] * m[2]);

    size_t sortedIndices[3];
    gsl_sort_index(sortedIndices, v, 1, 3);

    double m_copy[3];
    double s_copy[3];
    double v_copy[3];

    for (int i = 0; i < 3; i++) {
		m_copy[i] = m[i];
		s_copy[i] = s[i];
		v_copy[i] = v[i];
    }

    for (int i = 0; i < 3; i++) {
		m[sortedIndices[i]] = m_copy[i];
		s[sortedIndices[i]] = s_copy[i];
		v[sortedIndices[i]] = v_copy[i];
    }

    double sd[3];
    sd[0] = pow(v[0], 0.5);
    sd[1] = pow(v[1], 0.5);
    sd[2] = pow(v[2], 0.5);

    double T1 =
	(1.0 / m[0]) *
	(pow(1 + (9.0 / 4.0) * (pow(s[0], 4) / pow(m[0], 2)), 0.5) -
	 (3.0 / 2.0) * (pow(s[0], 2) / m[0]));
    double T2 =
	(1.0 / m[1]) *
	(pow(1 + (9.0 / 4.0) * (pow(s[1], 4) / pow(m[1], 2)), 0.5) -
	 (3.0 / 2.0) * (pow(s[1], 2) / m[1]));
    double T3 =
	(1.0 / m[2]) *
	(pow(1 + (9.0 / 4.0) * (pow(s[2], 4) / pow(m[2], 2)), 0.5) -
	 (3.0 / 2.0) * (pow(s[2], 2) / m[2]));

    distType *w[2];

    if (sd[0] < 0.01) {

		/*
		   %to estimate the error (in calculating probabilities the of the data),
		   %that results from approximating the first pdf as a point-mass
		   %distribtion, find the maximum of the absolute value of the
		   %derivative of the second pdf.
		 */
		double gp = DBL_MAX;
		double gp_tmp;
		for (int i = 1; i < 3; i++) {
			gp_tmp = gp_max_fixed(m[i], s[i]);
			if (gp_tmp < gp)
			gp = gp_tmp;
		}

		/*
		   %define the radius, r, of a small interval over which the second pdf is
		   %nearly constant and close to zero for t<r.
		   T=max(T2,T3);
		 */
		double T = DBL_MIN;
		if (T2 > T)
			T = T2;
		if (T3 > T)
			T = T3;

		/*
		   r=min(eps/(3*gp),T);
		 */
		double r = DBL_MAX;
		if (eps / (3 * gp) < r)
			r = eps / (3 * gp);
		if (T < r)
			r = T;

		/*
		   checkval=min(onestagepdf2(r,m(2),s(2)),onestagepdf2(r,m(3),s(3)));
		 */
		double checkval = DBL_MAX;
		double check_domain[1];
		check_domain[0] = r;
		distType check_Y1[1];
		distType check_Y2[1];
		waldpdf(check_domain, m[1], s[1], check_Y1, 1);
		waldpdf(check_domain, m[2], s[2], check_Y2, 1);
		if (check_Y1[0] < checkval)
			checkval = check_Y1[0];
		if (check_Y2[0] < checkval)
			checkval = check_Y2[0];

		/*
		   while checkval>=eps/2
		 */
		while (checkval >= eps / 2) {
			/*
			   r=r/2;
			   checkval=min(onestagepdf2(r,m(2),s(2)),onestagepdf2(r,m(3),s(3))) ;
			 */
			r = r / 2.0;
			waldpdf(check_domain, m[1], s[1], check_Y1, 1);
			waldpdf(check_domain, m[2], s[2], check_Y2, 1);
			if (check_Y1[0] < checkval)
			checkval = check_Y1[0];
			if (check_Y2[0] < checkval)
			checkval = check_Y2[0];
		}			/* end */

		/*
		   %estimate the maximum value of the convolution of the second two pdfs
		   gm=min(onestagepdf2(T2,m(2),s(2)), onestagepdf2(T3,m(3),s(3)));
		 */
		double gm = DBL_MAX;
		double T2_domain[1];
		T2_domain[0] = T2;
		distType T2_Y[1];
		waldpdf(T2_domain, m[1], s[1], T2_Y, 1);
		double T3_domain[1];
		T3_domain[0] = T3;
		distType T3_Y[1];
		waldpdf(T3_domain, m[2], s[2], T3_Y, 1);
		if (T2_Y[0] < gm)
			gm = T2_Y[0];
		if (T3_Y[0] < gm)
			gm = T3_Y[0];

		/*
		   %get the average value of the first pdf.  This is the point at which
		   %its mass is concentrated.
		 */
		double nu = 1.0 / m[0];

		//   %get upper limit of integral for approximating
		double Tu = DBL_MIN;
		if (100 > Tu)
			Tu = 100;
		if (nu + r + 1000.0 * sd[0] > Tu)
			Tu = nu + r + 1000.0 * sd[0];

		/*
		   checkerror=100;
		   hh=.001;
		 */
		double checkerror = 100;
		double hh = 0.001;

		/*
		   check1=.001*sum(onestagepdf2((0:.001:nu-r),m(1),s(1)))+.001*sum(onestagepdf2((nu+r:.001:Tu),m(1),s(1)));
		 */
		double check1 = 0;
		int partitionLength1 = (int) ((nu - r) / hh);
		int partitionLength2 = (int) ((Tu - (nu + r)) / hh);
		double *check1_X1 =
			(double *) malloc(partitionLength1 * sizeof(double));
		double *check1_X2 =
			(double *) malloc(partitionLength2 * sizeof(double));
		distType *check1_Y1 =
			(distType *) malloc(partitionLength1 * sizeof(distType));
		distType *check1_Y2 =
			(distType *) malloc(partitionLength2 * sizeof(distType));
		double tally = 0;
		for (int i = 0; i < partitionLength1; i++) {
			check1_X1[i] = tally;
			tally += hh;
		}
		tally = 0;

		/* FIXME FIXME FIXME */
		for (int i = (int) nu + r; i < partitionLength2; i++) {
			check1_X2[i] = tally;
			tally += hh;
		}
		waldpdf(check1_X1, m[0], s[0], check1_Y1, partitionLength1);
		waldpdf(check1_X2, m[0], s[0], check1_Y2, partitionLength2);
		for (int i = 0; i < partitionLength1; i++)
			check1 += check1_Y1[i];
		for (int i = 0; i < partitionLength2; i++)
			check1 += check1_Y2[i];
		check1 = check1 * hh;
		free(check1_X1);
		free(check1_Y1);
		free(check1_X2);
		free(check1_Y2);

		/*
		   while checkerror>10^(-4)
		 */
		while (checkerror > 0.0001) {
			/*
			   hh=.5*hh;
			   ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1)));
			   checkerror=abs(check1-ck1);
			   check1=ck1;
			 */
			hh = 0.5 * hh;
			double ck1 = 0;
			partitionLength1 = (int) ((nu - r) / hh);
			partitionLength2 = (int) ((Tu - (nu + r)) / hh);
			double *ck1_X1 =
			(double *) malloc(partitionLength1 * sizeof(double));
			double *ck1_X2 =
			(double *) malloc(partitionLength2 * sizeof(double));
			distType *ck1_Y1 =
			(distType *) malloc(partitionLength1 * sizeof(distType));
			distType *ck1_Y2 =
			(distType *) malloc(partitionLength2 * sizeof(distType));
			tally = 0;
			for (int i = 0; i < partitionLength1; i++) {
			ck1_X1[i] = tally;
			tally += hh;
			}
			tally = 0;
			for (int i = nu + r; i < partitionLength2; i++) {
			ck1_X2[i] = tally;
			tally += hh;
			}
			waldpdf(ck1_X1, m[0], s[0], ck1_Y1, partitionLength1);
			waldpdf(ck1_X2, m[0], s[0], ck1_Y2, partitionLength2);
			for (int i = 0; i < partitionLength1; i++)
			ck1 += ck1_Y1[i];
			for (int i = 0; i < partitionLength2; i++)
			ck1 += ck1_Y2[i];
			ck1 = ck1 * hh;
			free(ck1_X1);
			free(ck1_Y1);
			free(ck1_X2);
			free(ck1_Y2);
			check1 = ck1;

		}			/* end */

		/*
		   check2=gm*check1;
		 */
		double check2 = gm * check1;

		double logP0;

		/* if check2<=eps/3 */
		double l;

		if (check2 <= eps / 3.0) {
			/*
			   flag=1;
			   l=1/m(1);
			   P=convolv_2invG_adapt_nov(t-l,m(2),s(2),m(3),s(3),h);
			 */
			flag = 1;
			l = 1.0 / m[0];
#ifdef __INTEL_COMPILER
			double *X_tminusone =
			(double *) _mm_malloc(size_XY * sizeof(double), 32);
#else
			double *X_tminusone =
			(double *) malloc(size_XY * sizeof(double));
#endif
			for (int i = 0; i < size_XY; i++)
			X_tminusone[i] = X[i] - 1;
			conv2waldpdf(X_tminusone, m[1], s[1], m[2], s[2], Y, h, 0,
				 size_XY);
			free(X_tminusone);
		}
		/* else */
		else {
			/*
			   [y,flag]=convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h);
			 */
#ifdef __INTEL_COMPILER
			distType *y =
			(distType *) _mm_malloc(partitionLength * sizeof(distType),
						  32);
#else
			distType *y = (distType *) malloc(partitionLength * sizeof(distType));
#endif
			w[0] = y;
			//conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h, 0, partitionLength);

			/*
			   z=onestagepdf2(x',m(1),s(1));
			 */
#ifdef __INTEL_COMPILER
			distType *z =
			(distType *) _mm_malloc(partitionLength * sizeof(distType),
						  32);
#else
			distType *z = (distType *) malloc(partitionLength * sizeof(distType));
#endif

			w[1] = z;
			//waldpdf(x, m[0], s[0], z, partitionLength);
#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
			for (int i = 0; i < 2; i++) {
			if (i == 0)
				conv2waldpdf(x, m[1], s[1], m[2], s[2], w[0], h, 0,
					 partitionLength);
			else
				waldpdf(x, m[0], s[0], w[1], partitionLength);
			}

			/*
			   % approximate the convolution of the two pdfs
			   % the ith element of v gives approximates convolution over [0, x(i+1)]
			   % as a left-hand Riemann sum
			 */
			//printf("\n\ncalling approxconv_rep from conv3waldpdf partitionLength=%d size_XY=%d h=%f\n", partitionLength, size_XY, h);
			//approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);

#ifdef _BINNED_MODE
			binned_conv(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);
#else
			nn_conv(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);
#endif

#ifdef __INTEL_COMPILER
			_mm_free(y);
			_mm_free(z);
#else
			free(y);
			free(z);
#endif
		}			/* end */

		free(x);

		/*
		   P=P';
		   P0=max(realmin,P);
		   logP0=sum(log(P0));
		 */
		/*
		for (int i = 0; i < size_XY; i++) {
			//if (Y[i] < DBL_MIN)
			//Y[i] = DBL_MIN;
			logP0 += log(Y[i]);;
		}
		*/

		logP0 = loglikelihood(Y, size_XY);

		/*
		   while E>=.001*abs(logP0)

		 */
		double h1;
		double logP1;
		while (E >= 0.001 * fabs(logP0)) {
			/*
			   h1=.5*h;
			 */
			h1 = 0.5 * h;

			/* x=0:h1:Maxt; */
			/* produce a range of evenly spaced points to evaluate at between
			   0 and Maxt with step size h. the even spacing is important when
			   calculating the convolution later */
			partitionLength = (int) (maxX / h1);

			// This represents the partition
#ifdef __INTEL_COMPILER
			double *x =
			(double *) _mm_malloc(partitionLength * sizeof(double),
						  32);
#else
			double *x =
			(double *) malloc(partitionLength * sizeof(double));
#endif

			// fill the partition
			double tally = 0;
			for (int i = 0; i < partitionLength; i++) {
			x[i] = tally;
			tally += h1;
			}

			/*
			   [y,flag]=convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h1);
			 */
#ifdef __INTEL_COMPILER
			distType *y =
			(distType *) _mm_malloc(partitionLength * sizeof(distType),
						  32);
#else
			distType *y = (distType *) malloc(partitionLength * sizeof(distType));
#endif
			w[0] = y;
			//conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h1, 0, partitionLength);

			/*
			   z=onestagepdf2(x',m(1),s(1));
			 */
#ifdef __INTEL_COMPILER
			distType *z =
			(distType *) _mm_malloc(partitionLength * sizeof(distType),
						  32);
#else
			distType *z = (distType *) malloc(partitionLength * sizeof(distType));
#endif
			w[1] = z;
			//waldpdf(x, m[0], s[0], z, partitionLength);

#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
			for (int i = 0; i < 2; i++) {
			if (i == 0)
				conv2waldpdf(x, m[1], s[1], m[2], s[2], w[0], h1, 0,
					 partitionLength);
			else
				waldpdf(x, m[0], s[0], w[1], partitionLength);
			}

			/*
			   % approximate the convolution of the two pdfs
			   % the ith element of v gives approximates convolution over [.001, x(i+1)]
			   % as a left-hand Riemann sum
			 */
			//printf("\n\ncalling approxconv_rep from conv3waldpdf partitionLength=%d size_XY=%d h1=%f\n", partitionLength, size_XY, h1);
			//approxconvolv_replacement(z, y, X, x, Y, &logP1, partitionLength, size_XY, h1);

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

			// THIS IS WHERE I SHOULDBE UPDATING E
			E = fabs(logP1 - logP0);
			logP0 = logP1;
			h = h1;

			//printf("logP1=%f E=%f\n", logP1, E);

		}			/* end */

    }

    /*
       %otherwise, perform the convolution
       else
     */
    else {
	/*
	   flag = 0;
	   [y, flag] = convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h);

	 */
#ifdef __INTEL_COMPILER
	distType *y =
	    (distType *) _mm_malloc(partitionLength * sizeof(distType), 32);
#else
	distType *y = (distType *) malloc(partitionLength * sizeof(distType));
#endif
	w[0] = y;
	//conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h, 0, partitionLength);

	/*
	   z = onestagepdf2(x',m(1),s(1));

	 */
#ifdef __INTEL_COMPILER
	distType *z =
	    (distType *) _mm_malloc(partitionLength * sizeof(distType), 32);
#else
	distType *z = (distType *) malloc(partitionLength * sizeof(distType));
#endif
	w[1] = z;
	//waldpdf(x, m[0], s[0], z, partitionLength);

#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
	for (int i = 0; i < 2; i++) {
	    if (i == 0)
		conv2waldpdf(x, m[1], s[1], m[2], s[2], w[0], h, 0,
			     partitionLength);
	    else
		waldpdf(x, m[0], s[0], w[1], partitionLength);
	}

	/*
	   % approximate the convolution of the two pdfs
	   % the ith element of v gives approximates convolution over[.001, x(i + 1)]
	   % as a left - hand Riemann sum
	 */
	double logP0;
	//printf("\n\ncalling approxconv_rep from conv3waldpdf partitionLength=%d size_XY=%d h=%f\n", partitionLength, size_XY, h);
	//approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);
#ifdef _BINNED_MODE
				  binned_conv(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);
#else
				  nn_conv(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);
#endif

#ifdef __INTEL_COMPILER
	_mm_free(y);
	_mm_free(z);
#else
	free(y);
	free(z);
#endif
	/*
	   while E >= .001*abs(logP0)
	 */
	double h1;
	double logP1;
	while (E >= 0.001 * fabs(logP0)) {
	    /*
	       h1 = .5*h;
	     */
	    h1 = 0.5 * h;

	    /*
	       x = 0:h1:Maxt;
	     */
	    partitionLength = (int) (maxX / h1);
#ifdef __INTEL_COMPILER
	    double *x =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
#else
	    double *x = malloc(partitionLength * sizeof(double));
#endif
	    // fill the partition
	    double tally = 0;
	    for (int i = 0; i < partitionLength; i++) {
		x[i] = tally;
		tally += h1;
	    }

	    /*
	       [y, flag] = convolv_2invG_nov(x',m(2),s(2),m(3),s(3),h1);
	     */
#ifdef __INTEL_COMPILER
	    distType *y =
		(distType *) _mm_malloc(partitionLength * sizeof(distType),
				      32);
#else
	    distType *y = (distType *) malloc(partitionLength * sizeof(distType));
#endif
	    w[0] = y;
	    //conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h1, 0, partitionLength);

	    /*
	       z = onestagepdf2(x',m(1),s(1));
	     */
#ifdef __INTEL_COMPILER
	    distType *z =
		(distType *) _mm_malloc(partitionLength * sizeof(distType),
				      32);
#else
	    distType *z = (distType *) malloc(partitionLength * sizeof(distType));
#endif
	    w[1] = z;
	    //waldpdf(x, m[0], s[0], z, partitionLength);

#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
	    for (int i = 0; i < 2; i++) {
		if (i == 0)
		    conv2waldpdf(x, m[1], s[1], m[2], s[2], w[0], h1, 0,
				 partitionLength);
		else
		    waldpdf(x, m[0], s[0], w[1], partitionLength);

	    }

	    /*
	       % approximate the convolution of the two pdfs
	       % the ith element of v gives approximates convolution over[.001, x(i + 1)]
	       % as a left - hand Riemann sum
	     */
	    //printf("\n\ncalling approxconv_rep from conv3waldpdf partitionLength=%d size_XY=%d h1=%f\n", partitionLength, size_XY, h1);
	    //approxconvolv_replacement(z, y, X, x, Y, &logP1, partitionLength, size_XY, h1);
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

	    /*
	       E = abs(logP1 - logP0);

	       P0 = P1;
	       logP0 = logP1;
	       h = h1;
	     */
	    E = fabs(logP1 - logP0);
	    logP0 = logP1;
	    h = h1;

	    //printf("logP1=%f E=%f\n", logP1, E);

	}			/* end */
	/* P = P0; */

    }				/* end */

    // free(x);
}

/* End of code generation (convolv_3invG_nov.c) */
