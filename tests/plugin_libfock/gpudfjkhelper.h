#ifndef GPUDFJKHELPER_H
#define GPUDFJKHELPER_H

/**
 * fortran-ordered blas routines
 */
#include"blas.h"

namespace psi{
/**
 * Class GPUDFJKHelper
 *
 * Helper class that must be compiled with nvcc.  
 * Meant to be an extension of the GPUDFJK class
 * defined in jk.h
 *
 */
class GPUDFJKHelper{ 
 
public:
    /**
     * Constructor
     */
    GPUDFJKHelper();

    /**
     * Destructor
     */
    ~GPUDFJKHelper();

    /**
     * cublasDgemm wrapper. This function assumes that the input and output
     * to cublasDgemm can fit in device memory.
     */
    void GPU_DGEMM(char transa,char transb,int m,int n,int k,double alpha,double*A,int lda,double*B,int ldb,double beta,double*C,int ldc);

    /**
     * cublasDgemm wrapper that uses a 2-dimensional tile. The function
     * determines optimal block sizes such that input and output arrays can
     * fit in device memory.  For now, this is a threaded call, where
     * all threads < num_gpus execute dgemm on the gpu while all threads
     * >= num_gpus execute on the cpu.  TODO: make a separate function call
     * for the threaded version.  Currently, this function only works for
     * transa='n' and transb='t'.
     */
    void GPU_DGEMM_2DTile(char transa,char transb,int m,int n,int k,double alpha,double*A,int lda,double*B,int ldb,double beta,double*C,int ldc,int thread);

    /**
     * free cpu and gpu memory.
     */
    void Finalize();

    /**
     * check for errors from cuda
     */
    inline void Check_CUDA_Error(FILE*fp,const char *message);

    /**
     * variables:
     */
    double **tmp;
    int num_gpus,nthreads;

    /**
     * initialize temporary cpu memory
     */
    void  Initialize(int max_rows,int max_nocc,int nbf);

};
} // end of namespace

#endif
