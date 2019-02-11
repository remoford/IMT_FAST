#!/bin/bash
gcc -lm -lgsl -lgslcblas -fopenmp -std=c11 -O3 -march=bdver1 -ftree-vectorizer-verbose=2 *.c -o gcc_fast
