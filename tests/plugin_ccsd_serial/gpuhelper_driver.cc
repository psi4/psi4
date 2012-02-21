#include"blas.h"
#include"gpuhelper.h"

using namespace psi;
using namespace boost;

namespace psi{

GPUHelper::GPUHelper()
{}
GPUHelper::~GPUHelper()
{}

/**
  *  initialize cuda (if we have it)
  */
void GPUHelper::CudaInit(Options&options){
  max_mapped_memory=0;
  num_gpus=gpumemory=extraroom=0;
  num_cpus=0;
  num_gpus=0;
  #ifdef CUDA
    CudaInitGPU(options);
  #endif
}
/**
  *  clean up 
  */
void GPUHelper::Cleanup(Options&options){
  #ifdef CUDA
    CleanupGPU(options);
  #endif 
}

/**
  *  wrappers to gpu dgemm
  */
void GPUHelper::GPUTiledDGEMM(char transa,char transb,long int m, long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){
  if (num_gpus<1){
     F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     return;
  }
  if (transa=='n'){
     if (transb=='n'){
        //GPU_DGEMM_2DTile_nn_threaded_WithCpuStealing(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        GPU_DGEMM_2DTile_nn_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
     else{
        //GPU_DGEMM_2DTile_nt_threaded_WithCpuStealing(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        GPU_DGEMM_2DTile_nt_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
  }
  else{
     if (transb=='n'){
        //GPU_DGEMM_2DTile_tn_threaded_WithCpuStealing(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        GPU_DGEMM_2DTile_tn_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
     else{
        //GPU_DGEMM_2DTile_tt_threaded_WithCpuStealing(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        GPU_DGEMM_2DTile_tt_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
  }
}

void GPUHelper::GPUTiledDGEMM_NoThread(char transa,char transb,long int m, long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc,int thread){
  if (num_gpus<1){
     F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     return;
  }
  if (thread>=num_gpus){
     F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     return;
  }
  if (transa=='n'){
     if (transb=='n'){
        GPU_DGEMM_2DTile_nn(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,thread);
     }
     else{
        GPU_DGEMM_2DTile_nt(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,thread);
     }
  }
  else{
     if (transb=='n'){
        GPU_DGEMM_2DTile_tn(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,thread);
     }
     else{
        GPU_DGEMM_2DTile_tt(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,thread);
     }
  }
}




}
