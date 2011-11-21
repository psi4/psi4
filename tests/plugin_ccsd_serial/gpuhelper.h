#ifndef GPUHELPER_H
#define GPUHELPER_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// psi headers
#include <libplugin/plugin.h>
#include"psi4-dec.h"
#include<boost/shared_ptr.hpp>
#include<liboptions/liboptions.h>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>
#include<libchkpt/chkpt.h>
#include<libiwl/iwl.h>
#include <libpsio/psio.hpp>

#define NUMTHREADS 8
#define MAXBLOCKS 65535

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{


class GPUHelper{
  public:
    GPUHelper();
    ~GPUHelper();

    // tiled dgemm for gpu
    void GPUTiledDGEMM(char transa,char transb,long int m, long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc);
    void GPU_DGEMM(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc);

    // threaded tiled dgemm for gpu
    void GPU_DGEMM_2DTile_nn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc);
    void GPU_DGEMM_2DTile_nt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc);
    void GPU_DGEMM_2DTile_tt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc);
    void GPU_DGEMM_2DTile_tn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc);

    void CudaInit(Options&options);
    inline void Check_CUDA_Error(FILE*fp,const char *message);

    /**
      * the maximum amount of cpu memory dedicated to mapped memory.
      * the default value is num_gpus * (gpumemory-extraroom), which
      * can be quite large.
      */
    ULI max_mapped_memory;
    // available gpu memory
    ULI gpumemory;
    // wasted gpu memory
    ULI extraroom;
    // how large must a gemm be to go on the gpu?
    ULI gputhresh;
    // pointers to gpu and mapped cpu memory
    double**gpuarray,**tmp;

    // tiling
    void Tiling(double mem1,double mem1,long int m,long int n,long int k);
    long int ntilesN,ntilesM,ntilesK;
    long int tilesizeK,tilesizeN,tilesizeM;
    long int lasttileK,lasttileN,lasttileM;
    long int *tilesizesM,*tilesizesN,*tilesizesK;

    int ndoccact,nvirt,nmo,num_gpus;

};
};

#endif
