#!/bin/bash 

gcc -I/usr/local/include -L/usr/local/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm

./compute_covariances_fourier 100000

# setenv LD_LIBRARY_PATH /home/teifler/lib
# gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm 

#  ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov_parallel2/ NCoyote_fid_LSST_conti 1 0 1 20 30 5000 1 2 2 2 2 1 1 0 1 4 4 4 4 1 8
#./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../emulator/cov_playground2/ s-lhs.100_6_CosmicEmuOnly_range 1 0 1 20 30 5000 1 2 2 2 2 1 1 0 1 4 4 4 4 1 8

 #  ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov_parallel2/ NCoyote_fid_LSST_conti 1 0 1 20 30 5000 1 2 2 2 2 1 1 0 1 4 4 4 4 1 8

# 
# ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov/ Coyote_fid_DES_conti 1 0 1 20 30 5000 0 2 1 1 1 1 1 0 0 0 0 0 0 1
# 
# ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov/ Coyote_fid_LSST_conti 1 0 1 20 30 5000 0 2 3 3 3 1 1 0 0 0 0 0 0 1
