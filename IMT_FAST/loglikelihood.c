#include "main.h"
#include <float.h>

distType loglikelihood(distType likelihood[], long long int data_size){

	distType ll = 0;
	for (long long int i = 0; i < data_size; i++) {
		if (likelihood[i] == 0)
			ll += (distType) log((long double)distMin);
		else
			ll += (distType) log((long double)likelihood[i]);
	}
	return ll;
}