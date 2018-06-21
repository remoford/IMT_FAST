


#include "main.h"
#include <float.h>

distType loglikelihood(distType likelihood[], int data_size){

	double ll = 0;
	for (int i = 0; i < data_size; i++) {
		if (likelihood[i] == 0)
			ll += log(distMin);
		else
			ll += log(likelihood[i]);
	}
	return ll;
}