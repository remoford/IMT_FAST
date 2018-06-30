#!/bin/bash
icc -lgsl -lgslcblas -fopenmp -std=c11 -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -no-fast-transcendentals *.c -o icc_fast
