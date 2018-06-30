#!/bin/bash
gcc -lm -lgsl -lgslcblas -fopenmp -std=c11 -O3 *.c -o gcc_fast
