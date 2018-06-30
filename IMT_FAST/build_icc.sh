#!/bin/bash
icc -lgsl -lgslcblas -std=c11 -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -qopenmp-link=static -qopenmp *.c -o icc_fast
