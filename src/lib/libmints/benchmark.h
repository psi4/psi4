#ifndef _psi_src_lib_libmints_bench_h
#define _psi_src_lib_libmints_bench_h

namespace psi {

/**
* Perform a benchmark traverse of BLAS 1 routines on 
* the current hardware
* \param N maximum dimension exponent (requires ~ 3 (2^N x 2^N) 
* double matrices
* \param min_time minimum amount of time to run each routine [s]
* \param nthread maximum number of threads to use
**/
void benchmark_blas1(int N, double min_time);
/**
* Perform a benchmark traverse of BLAS 2 routines on 
* the current hardware
* \param N maximum dimension exponent (requires ~ 3 (2^N x 2^N) 
* double matrices
* \param min_time minimum amount of time to run each routine [s]
* \param nthread maximum number of threads to use
**/
void benchmark_blas2(int N, double min_time);
/**
* Perform a benchmark traverse of BLAS 3 and LAPACK routines on 
* the current hardware
* \param N maximum dimension exponent (requires ~ 4 (2^N x 2^N) 
* double matrices
* \param min_time minimum amount of time to run each routine [s]
* \param nthread maximum number of threads to use
**/
void benchmark_blas3(int N, double min_time, int nthread = 1);
/**
* Perform a benchmark of PSIO disk performance on
* the current hardware
* \param N maximum dimension exponent (requires 1 (2^N x 2^N) 
* double matrices
* \param min_time minimum amount of time to run each routine [s]
**/
void benchmark_disk(int N, double min_time);
/**
* Perform a benchmark of psi integrals (of libmints type)
* on the current hardware
* All integrals will be called from different centers
* \param max_am maximum am to consider 
* \param min_time minimum time to run each shell combination of
* each integral type
**/
void benchmark_integrals(int max_am, double min_time);
/**
* Perform a benchmark of common double floating
* point operations, including most of cmath
* \param min_time minimum amount of time to run each routine [s]
**/
void benchmark_math(double min_time);

}

#endif
