#include "imt_analysis.h"
#include "main.h"
#include "stdio.h"
#include "math.h"       /* math_errhandling */
#include <errno.h>      /* errno, EDOM */
#include <fenv.h>       /* feclearexcept, fetestexcept, FE_ALL_EXCEPT, FE_INVALID */
#include "float.h"
#include "onestage.h"
#include "utility.h"
#include "loglikelihood.h"
#include "onestage.h"
#include "twostage.h"
#include "threestage.h"
#include "time.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

traceDepth = 0;


int main(int argc, char * argv[])
{
    (void) argc;
    (void) argv;

	printf("argc = %d", argc);
	for(int i = 0; i < argc ; i++){
		printf("argv[%d] = %s", i, argv[i]);
	}

	/*
	Precision can vary from machine to machine, show the ground truth
	*/
	printf("FLT_MIN=%g\n", (double)FLT_MIN);
	printf("FLT_MAX=%g\n", (double)FLT_MAX);
	printf("FLT_MIN_EXP=%g\n", (double)FLT_MIN_EXP);
	printf("FLT_MAX_EXP=%g\n", (double)FLT_MAX_EXP);

	printf("DBL_MIN=%g\n", DBL_MIN);
	printf("DBL_MAX=%g\n", DBL_MAX);
	printf("DBL_MIN_EXP=%g\n", (double)DBL_MIN_EXP);
	printf("DBL_MAX_EXP=%g\n", (double)DBL_MAX_EXP);

	printf("LDBL_MIN=%f\n", (double)LDBL_MIN);
	printf("LDBL_MAX=%f\n", (double)LDBL_MAX);
	printf("LDBL_MIN_EXP=%f\n", (double)LDBL_MIN_EXP);
	printf("LDBL_MAX_EXP=%f\n", (double)LDBL_MAX_EXP);


	/*
	Error checking depends on compiler options, show the ground truth
	*/
	errno = 0;
	if (math_errhandling & MATH_ERREXCEPT) feclearexcept(FE_ALL_EXCEPT);

	printf("Error handling: %d\n", math_errhandling);

	sqrt(-1);
	if (math_errhandling & MATH_ERRNO) {
		if (errno == EDOM) printf("errno set to EDOM\n");
	}
	if (math_errhandling & MATH_ERREXCEPT) {
		if (fetestexcept(FE_INVALID)) printf("FE_INVALID raised\n");
	}
	errno = 0;
	printf("\n");


	/*
	Parse command line arguments
	*/
	if (argc < 2) {
		printf("No arguments given! Usage: <usage goes here>\n");
	}
	else if( strcmp(argv[1], "analysis") == 0 ) {
		/*
		Run a maximum likelihood analysis with model argv[2] and an optional data file argv[3]
		*/

		if (argc == 3)
			IMT_analysis_April2017(argv[2], "");

		else if (argc == 4)
			IMT_analysis_April2017(argv[2], argv[3]);

		else
			printf("Invalid number of arguments for optimize mode. Please provide the model and optionally a data file like \"optimize twostage data/DMSO.txt\"\n");
	}

	else if (strcmp(argv[1], "single") == 0) {
		/*
		Run a single model evaluation
		*/

		clock_t t;
		t = clock();

		if (argc == 2)
			printf("Invalid number of arguments for single mode. Please provide the model and it's requisite parameters and optionally a data file like \"single twostage 0.03 0.2 0.5 0.45 data/DMSO.txt\"\n");

		else if (strcmp(argv[2], "onestage") == 0) {
			/*
			Evaluate a one stage model with parameters m=argv[3] and s=argv[4] and an optional data file argv[5] and print its loglikelihood
			*/
			if (argc == 5 || argc == 6) {

				int data_size;
				distType * data;

				if (argc == 5)
					data = get_default_data(&data_size);
				else
					data = readfile((char *)argv[5], &data_size);

				distType * Y = (distType *)MALLOC(sizeof(distType)*data_size);

				double m = atof(argv[3]);
				double s = atof(argv[4]);

				wald_adapt(data, m, s, Y, data_size);

				double ll = (double)loglikelihood(Y, data_size);

				printf("loglikelihood = %f\n", ll);

				FREE(Y);
				FREE(data);
			}
			else
				printf("Invalid number of parameters\n");
		}

		else if (strcmp(argv[2], "twostage") == 0) {
			/*
			Evaluate a two stage model with parameters argv[3] through argv[6] and an optional data file argv[7] and print its loglikelihood
			*/

			if (argc == 7 || argc == 8) {

				int data_size;
				distType * data;

				if (argc == 7)
					data = get_default_data(&data_size);
				else
					data = readfile((char *)argv[7], &data_size);

				distType * Y = (distType *)MALLOC(sizeof(distType)*data_size);

				double m1 = atof(argv[3]);
				double s1 = atof(argv[4]);
				double m2 = atof(argv[5]);
				double s2 = atof(argv[6]);

				conv2waldpdf(data, m1, s1, m2, s2, Y, 0.01, 1, data_size, LOGLIKELIHOOD);

				double ll = (double)loglikelihood(Y, data_size);

				printf("loglikelihood = %f\n", ll);

				FREE(Y);
				FREE(data);

			} else
				printf("Invalid number of parameters\n");
		}

		else if (strcmp(argv[2], "threestage") == 0) {
			/*
			Evaluate a three stage model with parameters argv[3] through argv[8] and an optional data file argv[9] and print its loglikelihood
			*/

			if (argc == 9 || argc == 10) {

				int data_size;
				distType * data;

				if (argc == 9)
					data = get_default_data(&data_size);
				else
					data = readfile((char *)argv[9], &data_size);

				distType * Y = (distType *)MALLOC(sizeof(distType)*data_size);

				double m1 = atof(argv[3]);
				double s1 = atof(argv[4]);
				double m2 = atof(argv[5]);
				double s2 = atof(argv[6]);
				double m3 = atof(argv[7]);
				double s3 = atof(argv[8]);

				threestage_adapt(data, m1, s1, m2, s2, m3, s3, Y, data_size);

				double ll = (double)loglikelihood(Y, data_size);

				printf("loglikelihood = %f\n", ll);

				FREE(Y);
				FREE(data);

			} else
				printf("Invalid number of parameters\n");
		}

		printf("single evaluation took %fs\n", ((float)t) / CLOCKS_PER_SEC);
	}
	else
		printf("Invalid mode %s\n", argv[1]);
	
    return 0;
}
