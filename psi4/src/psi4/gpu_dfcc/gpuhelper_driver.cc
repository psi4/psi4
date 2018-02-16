/*
 *@BEGIN LICENSE
 *
 * GPU-accelerated density-fitted coupled-cluster, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include"blas.h"
#include"gpuhelper.h"
#define CUDA

#include <omp.h>

using namespace psi;
using namespace std;

namespace psi{ namespace fnocc{

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
  * finalize cuda (if we have it)
  */
void GPUHelper::CudaFinalize(Options&options){
  #ifdef CUDA
    CudaFinalizeGPU(options);
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
  if (num_cpus>0){
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
  }else{
     if (transa=='n'){
        if (transb=='n'){
           GPU_DGEMM_2DTile_nn_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        }
        else{
           GPU_DGEMM_2DTile_nt_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        }
     }
     else{
        if (transb=='n'){
           GPU_DGEMM_2DTile_tn_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        }
        else{
           GPU_DGEMM_2DTile_tt_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        }
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




}}
