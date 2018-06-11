
/* Include files */
#include "IMT_analysis_April2017.h"
#include "convolv_3invG_nov.h"
#include "onestagepdf2.h"
#include "convolv_2invG_adapt_nov.h"
#include "convolv_2invG_nov.h"
#include "gp_max.h"
#include "emgpdf.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sort.h"
#include "approxconvolv.h"

#define _GSL_GP_MAX_FIXED
//#define _PARALLEL_PDF

/* Function Definitions */

double convolv_3invG_nov_loglikelihood(const gsl_vector * v, void *params)
{
    double *data = (double *) params;

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

    double Y[266];
    double Y_replacement[266];
    for (int i = 0; i < 266; i++) {
	Y[i] = 0;
	Y_replacement[i] = 0;
    }

    convolv3waldpdf(m1, s1, m2, s2, m3, s3, data, Y, 266, 0.1);
    /*
       printf("Y = \n");
       for (int i = 0; i < 266; i++) {
       printf("%.17f ", Y[i]);
       if (i % 8 == 7)
       printf("\n");
       }
       printf("\n");
     */

    /*
       convolv_3invG_nov(m1, s1, m2, s2, m3, s3, Y);
       printf("Y = \n");
       for (int i = 0; i < 266; i++) {
       printf("%.17f ", Y[i]);
       if (i % 8 == 7)
       printf("\n");
       }
       printf("\n");
     */

    double loglikelihood = 0;
    double elementLogLikelihood;
    for (int i = 0; i < 266; i++) {
	elementLogLikelihood = log(Y[i]);
	//if (elementLogLikelihood > DBL_MAX)
	//printf("elementLogLikelihood=inf Y[%d]=%f\n", i, Y[i]);
	loglikelihood += elementLogLikelihood;
    }

    //exit(1);

    double objective = penalty - loglikelihood;
    //printf("\nobjective = %f\n\n", objective);

    return objective;
}

/*
%this function evaluates the convolution of three inverse gaussian
%distributions at vector t
function [P,h,flag,E]=convolv_3invG_nov(t,m1,s1,m2,s2,m3,s3,h)
*/
//void convolv3waldpdf(double m1, double s1, double m2, double s2, double m3, double s3, const double X[266], double Y[266], int size_XY, double h)
void
convolv3waldpdf(double m1, double s1, double m2, double s2, double m3,
		double s3, const double X[], double Y[], int size_XY,
		double h)
{
    /*
       flag=0;
     */
    int flag = 0;

    /*
       n=length(t);
       Maxt=max(t);
     */
    // find the largest point in t
    double maxX = 0;
    for (int i = 0; i < size_XY; i++) {
	if (X[i] > maxX)
	    maxX = X[i];
    }

    /*
       x=0:h:Maxt;
     */
    /* produce a range of evenly spaced points to evaluate at between
       0 and Maxt with step size h. the even spacing is important when
       calculating the convolution later */
    int partitionLength = (int) (maxX / h);

    // This represents the partition
#ifdef __INTEL_COMPILER
    double *x =
	(double *) _mm_malloc(partitionLength * sizeof(double), 32);
#else
    double *x = malloc(partitionLength * sizeof(double));
#endif

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

    /*
       E=Inf;
       eps=.01;
     */
    double E = DBL_MAX;
    double eps = 0.01;

    /*
       m=[m1 m2 m3];
       s=[s1 s2 s3];
     */
    double m[3];
    m[0] = m1;
    m[1] = m2;
    m[2] = m3;
    double s[3];
    s[0] = s1;
    s[1] = s2;
    s[2] = s3;

    /*
       m=abs(m);
       s=abs(s);
     */
    m[0] = fabs(m[0]);
    m[1] = fabs(m[1]);
    m[2] = fabs(m[2]);
    s[0] = fabs(s[0]);
    s[1] = fabs(s[1]);
    s[2] = fabs(s[2]);

    /*
       v=(s.^2)./(m.^3);
     */
    // find the variance for both sub-distributions
    double v[3];
    v[0] = (s[0] * s[0]) / (m[0] * m[0] * m[0]);
    v[1] = (s[1] * s[1]) / (m[1] * m[1] * m[1]);
    v[2] = (s[2] * s[2]) / (m[2] * m[2] * m[2]);

    /*
       [v,I]=sort(v);
       m=m(I);
       s=s(I);
     */
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

    /*
       sd=v.^.5;
     */
    double sd[3];
    sd[0] = pow(v[0], 0.5);
    sd[1] = pow(v[1], 0.5);
    sd[2] = pow(v[2], 0.5);

    /*
       T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
       T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));
       T3=(1/m(3))*((1+(9/4)*(s(3)^4/m(3)^2))^.5-(3/2)*(s(3)^2/m(3)));
     */
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

    double *w[2];

    /*
       if sd(1)<.01
     */
    if (sd[0] < 0.01) {

	/*
	   %to estimate the error (in calculating probabilities the of the data),
	   %that results from approximating the first pdf as a point-mass
	   %distribtion, find the maximum of the absolute value of the
	   %derivative of the second pdf.

	   gp=min(gp_max(m(2),s(2)),gp_max(m(3),s(3)));
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
	double check_Y1[1];
	double check_Y2[1];
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
	double T2_Y[1];
	waldpdf(T2_domain, m[1], s[1], T2_Y, 1);
	double T3_domain[1];
	T3_domain[0] = T3;
	double T3_Y[1];
	waldpdf(T3_domain, m[2], s[2], T3_Y, 1);
	if (T2_Y[0] < gm)
	    gm = T2_Y[0];
	if (T3_Y[0] < gm)
	    gm = T3_Y[0];

	/*
	   %get the average value of the first pdf.  This is the point at which
	   %its mass is concentrated.
	   nu=1/m(1);
	 */
	double nu = 1.0 / m[0];

	/*
	   %get upper limit of integral for approximating
	   %int_{r+nu}^{infty}f(s)ds.
	   Tu=max(100,nu+r+1000*sd(1));
	 */
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
	double *check1_Y1 =
	    (double *) malloc(partitionLength1 * sizeof(double));
	double *check1_Y2 =
	    (double *) malloc(partitionLength2 * sizeof(double));
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
	    double *ck1_Y1 =
		(double *) malloc(partitionLength1 * sizeof(double));
	    double *ck1_Y2 =
		(double *) malloc(partitionLength2 * sizeof(double));
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
	    double *y =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
#else
	    double *y = malloc(partitionLength * sizeof(double));
#endif
	    w[0] = y;
	    //conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h, 0, partitionLength);

	    /*
	       z=onestagepdf2(x',m(1),s(1));
	     */
#ifdef __INTEL_COMPILER
	    double *z =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
#else
	    double *z = malloc(partitionLength * sizeof(double));
#endif

	    w[1] = z;
	    //waldpdf(x, m[0], s[0], z, partitionLength);
#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
	    for (int i = 0; i < 2; i++) {
		if (i = 0)
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
	    approxconvolv_replacement(z, y, X, x, Y, &logP0,
				      partitionLength, size_XY, h);

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
	for (int i = 0; i < size_XY; i++) {
	    if (Y[i] < DBL_MIN)
		Y[i] = DBL_MIN;
	    logP0 += log(Y[i]);;
	}

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
	    double *y =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
#else
	    double *y = malloc(partitionLength * sizeof(double));
#endif
	    w[0] = y;
	    //conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h1, 0, partitionLength);

	    /*
	       z=onestagepdf2(x',m(1),s(1));
	     */
#ifdef __INTEL_COMPILER
	    double *z =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
#else
	    double *z = malloc(partitionLength * sizeof(double));
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
	    approxconvolv_replacement(z, y, X, x, Y, &logP1,
				      partitionLength, size_XY, h1);

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
	double *y =
	    (double *) _mm_malloc(partitionLength * sizeof(double), 32);
#else
	double *y = malloc(partitionLength * sizeof(double));
#endif
	w[0] = y;
	//conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h, 0, partitionLength);

	/*
	   z = onestagepdf2(x',m(1),s(1));

	 */
#ifdef __INTEL_COMPILER
	double *z =
	    (double *) _mm_malloc(partitionLength * sizeof(double), 32);
#else
	double *z = malloc(partitionLength * sizeof(double));
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
	approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength,
				  size_XY, h);

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
	    double *y =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
#else
	    double *y = malloc(partitionLength * sizeof(double));
#endif
	    w[0] = y;
	    //conv2waldpdf(x, m[1], s[1], m[2], s[2], y, h1, 0, partitionLength);

	    /*
	       z = onestagepdf2(x',m(1),s(1));
	     */
#ifdef __INTEL_COMPILER
	    double *z =
		(double *) _mm_malloc(partitionLength * sizeof(double),
				      32);
#else
	    double *z = malloc(partitionLength * sizeof(double));
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
	    approxconvolv_replacement(z, y, X, x, Y, &logP1,
				      partitionLength, size_XY, h1);
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
