#include "main.h"
#include <float.h>
#include "stdio.h"

distType loglikelihood(distType likelihood[], long data_size){

	distType llSum = 0;
	distType ll = 0;
	for (long i = 0; i < data_size; i++) {

		if (!isfinite(likelihood[i]))
			printf("OH NOES NONFINITE LIKELIHOOD!!!");

		
		if (likelihood[i]+0.0 == 0.0)
			ll = (distType) log((double)distMin);
			//ll = distMin;
		else
			ll = (distType) log((double)likelihood[i]);

		if (!isfinite(ll)) {
			printf("[OH NOES NONFINITE LOGLIKELIHOOD!!! likelihood=%f  ---  ", (double)likelihood[i]);
			
			
			for (long k = 0; k < data_size; k++)
				printf("%f ", likelihood[k]);

				

			printf("]");
			ll = distMin;
		}
		
		llSum += ll;
	}

	/*
	if (!isfinite(llSum)) {
		printf("OH NOES NONFINITE LOGLIKELIHOODSUM!!!");
		for (long i = 0; i < data_size; i++)
			printf("%f ", likelihood[i]);
	}
	*/

	return llSum;
}