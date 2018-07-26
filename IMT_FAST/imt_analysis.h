#ifndef IMT_ANALYSIS_APRIL2017_H
#define IMT_ANALYSIS_APRIL2017_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include "main.h"

extern distType * get_default_data(int * data_size);

extern void IMT_analysis_April2017(const char *model, char * data_filename);

#endif