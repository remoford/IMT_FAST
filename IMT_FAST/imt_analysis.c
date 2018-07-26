#include "imt_analysis.h"
#include "emgpdf.h"
#include "onestage.h"
#include "onestagepdf_lag.h"
#include "twostage.h"
#include "threestage.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics_double.h"
#include "time.h"
#include "main.h"
#include "loglikelihood.h"
#include "utility.h"
#include "time.h"
#include "gsl/gsl_statistics.h"

#define _VERBOSE
#define _CONV2WALD


distType * get_default_data(int * data_size) {
	printf("Using default FUCCI data\n");

	*data_size = 266;

	distType default_data[266] =
	{ 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
		12.4, 13.7, 12.4, 11.9, 10.3, 12.9, 14.7, 11.6, 13.4, 13.4, 11.6, 10.3, 9.3,
		13.7, 9.6, 10.1, 9.8, 10.9, 16.0, 9.3, 9.6, 10.3, 11.4, 10.6, 8.5, 10.3,
		11.1, 8.0, 10.6, 7.5, 12.9, 9.0, 8.5, 12.4, 11.6, 9.6, 9.6, 14.7, 9.8, 10.3,
		12.1, 8.8, 10.6, 12.1, 13.4, 12.4, 8.8, 13.2, 10.1, 11.6, 11.1, 15.8, 12.1,
		12.7, 12.7, 11.1, 13.2, 11.9, 12.4, 13.2, 14.0, 8.0, 8.8, 9.3, 16.5, 14.5,
		10.1, 14.2, 7.8, 13.2, 8.8, 8.8, 10.1, 11.9, 12.9, 14.5, 10.9, 10.6, 14.0,
		8.8, 8.8, 9.0, 10.9, 14.5, 9.6, 12.4, 11.9, 12.4, 11.1, 14.5, 10.3, 12.4,
		12.7, 11.9, 10.3, 13.7, 15.5, 14.5, 11.6, 10.6, 15.5, 14.7, 8.8, 11.6, 8.3,
		17.6, 12.4, 11.6, 15.0, 13.7, 12.7, 10.9, 7.2, 8.5, 8.3, 9.6, 11.4, 12.9,
		11.6, 13.4, 10.1, 11.6, 8.8, 12.4, 10.3, 16.3, 10.9, 10.1, 8.8, 9.3, 15.2,
		8.5, 11.1, 8.3, 11.4, 11.9, 9.3, 9.8, 16.3, 12.7, 9.0, 11.9, 9.3, 10.3, 13.4,
		11.4, 12.9, 12.4, 9.6, 10.3, 13.2, 10.6, 9.8, 11.9, 14.2, 13.4, 9.3, 9.6,
		12.1, 11.9, 10.1, 14.0, 12.9, 21.7, 11.6, 12.1, 10.3, 9.8, 14.2, 13.7, 7.2,
		10.9, 10.1, 9.6, 13.4, 13.2, 16.3, 11.6, 14.0, 10.9, 14.2, 12.4, 12.4, 13.4,
		17.6, 10.1, 10.9, 14.0, 12.9, 9.0, 13.4, 15.0, 16.0, 8.0, 9.8, 12.4, 8.5,
		9.6, 12.7, 12.1, 15.0, 16.0, 10.9, 14.2, 13.7, 11.9, 16.8, 11.4, 13.4, 12.4,
		22.0, 12.4, 16.8, 12.1, 10.3, 13.4, 11.6, 10.1, 14.5, 10.6, 11.9, 15.5, 9.8,
		12.4, 10.1, 8.0, 9.0, 9.3, 13.2, 11.1, 12.7, 12.1, 10.1, 13.2, 14.5, 10.1,
		12.7, 12.9, 11.9, 12.4, 11.1, 8.5, 14.5, 16.5, 12.4, 9.0, 11.1, 9.8, 11.1,
		11.1, 8.8, 13.2, 17.6, 16.8, 10.9, 12.4, 8.5, 14.7 };

	distType * data;
	data = (distType *)MALLOC(sizeof(distType)* *data_size);

	for (int i = 0; i < *data_size; i++)
		data[i] = default_data[i];

	return data;

}

void IMT_analysis_April2017(const char *model, char * data_filename) {
	printf("Using GSL minimizer\n");

	clock_t t;
	t = clock();

	int data_size;
	distType * data;
	if  (strcmp(data_filename, "") == 0) {
		data = get_default_data(&data_size);

	}
	else {
		printf("Reading data from file %s\n", data_filename);
		data = readfile(data_filename, &data_size);
	}

	if (data_size <= 0) {
		printf("ERROR: returned data_size <= 0 !!!\n");
		return;
	}

	if (data == NULL) {
		printf("ERROR: returned data == NULL\n");
		return;
	}

	qsort(data, data_size, sizeof(distType), compare);

	printf("\ndata[%d]={\n", data_size);
	for (int i = 0; i < data_size; i++) {
		printf("%.17f ", data[i]);
		if (i % 8 == 7)
			printf("\n");
	}
	printf("}\n");

	configStruct config;

	config.data = data;
	config.data_size = data_size;

	int emgfitnomle = 0;
    int twostagefitnomle = 0;
    int onestagelagnomle = 0;
    int onestagefitnomle = 0;
    int threestagefitnomle = 0;

    const char *twostagefitnomle_str = "twostage";
    const char *onestagelagnomle_str = "onestagelag";
    const char *onestagefitnomle_str = "onestage";
    const char *emgfitnomle_str = "emg";
    const char *threestagefitnomle_str = "threestage";
    const char *all_str = "all";

	/* choose model to fit */
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
		optimize_emg(data, data_size, config);
	}

	if (onestagefitnomle == 1) {
		optimize_onestage(data, data_size, config);
	}

	if (onestagelagnomle == 1) {
		optimize_onestagelag(data, data_size, config);
	}

	if (twostagefitnomle) {

		int numseeds_twostage;

		double mean = gsl_stats_mean(data, 1, data_size);

		double variance = gsl_stats_variance(data, 1, data_size);

		double ** seeds_twostage = twostage_seeds(mean, variance, &numseeds_twostage);

		printf("\nseeds[%d]:\n", numseeds_twostage);
		for (int seedIdx = 0; seedIdx < numseeds_twostage; seedIdx++)
			printf("%f %f %f %f\n", seeds_twostage[seedIdx][0], seeds_twostage[seedIdx][1], seeds_twostage[seedIdx][2], seeds_twostage[seedIdx][3]);
		printf("\n");


		optimize_twostage(data_size, data, numseeds_twostage, seeds_twostage, config);

		for (int seedIdx = 0; seedIdx < numseeds_twostage; seedIdx++)
			FREE(seeds_twostage[seedIdx]);
		FREE(seeds_twostage);
	}


	if (threestagefitnomle == 1) {
		optimize_threestage(data, data_size, config);
	}

	FREE(data);

	t = clock() - t;
	printf("total runtime = %fs\n", ((float)t) / CLOCKS_PER_SEC);
}
