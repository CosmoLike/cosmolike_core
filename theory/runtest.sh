#!/bin/bash 

rm tests

gcc -I/usr/local/include -L/usr/local/lib -o ./tests tests.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -L../class -lclass

./tests
