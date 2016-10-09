#ifndef __MIC_DGEMM_H__
#define __MIC_DGEMM_H__


#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <mkl.h>
#include <sys/time.h>


#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)
#define D(x) ((double)(x))

void init_pcl_dgemm (int matrix_size, int usemic);

void deinit_pcl_dgemm (void);

void pcl_dgemm (CBLAS_ORDER order, CBLAS_TRANSPOSE transA,
                CBLAS_TRANSPOSE transB, int M, int N, int K, double alpha,
                double *A, int lda, double *B, int ldb, double beta,
                double *C, int ldc);


#endif /* __MIC_DGEMM_H__ */