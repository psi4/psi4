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

#include"globals.h"
#include"blas.h"
#include"gpuhelper.h"
#include"gpuonly.h"
#ifdef _OPENMP
  #include<omp.h>
#endif


using namespace psi;
using namespace boost;

namespace psi{

GPUHelper::GPUHelper()
{}
GPUHelper::~GPUHelper()
{}

inline void GPUHelper::Check_CUDA_Error(FILE*fp,const char *message){
  cudaError_t error = cudaGetLastError();
  if (error!=cudaSuccess) {
     fprintf(fp,"\n  ERROR: %s: %s\n\n", message, cudaGetErrorString(error) );
     fflush(fp);
     exit(-1);
  }
}

/*===================================================================

  initialize cublas and get device properties

===================================================================*/
void GPUHelper::CudaInit(Options&options){

  num_gpus=gpumemory=extraroom=0;
  cudaGetDeviceCount(&num_gpus);
  if (num_gpus>0){
     cublasInit();
     struct cudaDeviceProp cudaProp;
     int gpu_id;
     cudaGetDevice(&gpu_id);
     cudaGetDeviceProperties( &cudaProp,gpu_id );
     fprintf(outfile,
       "\n  _________________________________________________________\n");
     fprintf(outfile,"  CUDA device properties:\n");
     fprintf(outfile,"  name:                 %20s\n",cudaProp.name);
     fprintf(outfile,"  major version:        %20d\n",cudaProp.major);
     fprintf(outfile,"  minor version:        %20d\n",cudaProp.minor);
     fprintf(outfile,"  canMapHostMemory:     %20d\n",cudaProp.canMapHostMemory);
     fprintf(outfile,"  totalGlobalMem:       %20lu mb\n",
       cudaProp.totalGlobalMem/(1024*1024));
     fprintf(outfile,"  sharedMemPerBlock:    %20lu\n",cudaProp.sharedMemPerBlock);
     fprintf(outfile,"  clockRate:            %20.3f ghz\n",
       cudaProp.clockRate/1.0e6);
     fprintf(outfile,"  regsPerBlock:         %20d\n",cudaProp.regsPerBlock);
     fprintf(outfile,"  warpSize:             %20d\n",cudaProp.warpSize);
     fprintf(outfile,"  maxThreadsPerBlock:   %20d\n",cudaProp.maxThreadsPerBlock);
     fprintf(outfile,
       "  _________________________________________________________\n\n");
     fflush(outfile);

     gpumemory = cudaProp.totalGlobalMem;

     extraroom = 350*1024*1024;
     
     cudaThreadExit();

     // default memory for mapped cpu memory is the sum of all gpu memory
     max_mapped_memory = num_gpus * (gpumemory-extraroom);
     if (options["MAX_MAPPED_MEMORY"].has_changed()){
        ULI temp_mem = options.get_int("MAX_MAPPED_MEMORY");
        temp_mem *= 1024*1024;
        if (temp_mem<max_mapped_memory)
           max_mapped_memory = options.get_int("MAX_MAPPED_MEMORY");
     }

     fprintf(outfile,"\n");
     fprintf(outfile,"  allocating gpu memory...");
     fflush(outfile);
     tmp = (double**)malloc(num_gpus*sizeof(double*));
     gpuarray = (double**)malloc(num_gpus*sizeof(double*));
     #pragma omp parallel for schedule (static) num_threads(num_gpus)
     for (int i=0; i<num_gpus; i++){
         int thread = 0;
         #ifdef _OPENMP
           thread = omp_get_thread_num();
         #endif
         cudaSetDevice(thread);
         Check_CUDA_Error(stdout,"cudaSetDevice");
         cudaMallocHost((void**)&tmp[thread],max_mapped_memory/num_gpus);  
         Check_CUDA_Error(outfile,"cpu tmp");
         cudaMalloc((void**)&gpuarray[thread],gpumemory-extraroom);
         Check_CUDA_Error(outfile,"gpu memory");
     }
     fprintf(outfile,"done.\n");
     fprintf(outfile,"\n");
     fflush(outfile);
  }
}

/**
 * dgemm assuming no tiling is necessary
 */
void GPUHelper::GPU_DGEMM(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){
  double*gpuA,*gpuB,*gpuC;
  cudaMalloc((void**)&gpuA,m*k*sizeof(double));
  cudaMalloc((void**)&gpuB,n*k*sizeof(double));
  cudaMalloc((void**)&gpuC,m*n*sizeof(double));
  cudaMemcpy(gpuA,A,m*k*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(gpuB,B,n*k*sizeof(double),cudaMemcpyHostToDevice);
  cublasDgemm(transa,transb,m,n,k,alpha,gpuA,lda,gpuB,ldb,beta,gpuC,ldc);
  cudaMemcpy(C,gpuC,m*n*sizeof(double),cudaMemcpyDeviceToHost);
  cudaFree(gpuA);
  cudaFree(gpuB);
  cudaFree(gpuC);
}
/**
 * dgemm using a 2-dimensional tile.
 */
void GPUHelper::GPUTiledDGEMM(char transa,char transb,long int m, long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){
  if (num_gpus<1){
     F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     return;
  }
  /* if (thread>=num_gpus){
     F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
  }*/
  if (transa=='n'){
     if (transb=='n'){
        GPU_DGEMM_2DTile_nn_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        //F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
     else{
        GPU_DGEMM_2DTile_nt_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        //F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
  }
  else{
     if (transb=='n'){
        GPU_DGEMM_2DTile_tn_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        //F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
     else{
        GPU_DGEMM_2DTile_tt_threaded(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
        //F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     }
  }
}
/**
 * dgemm using a 2-dimensional tile - threaded versions for multiple gpus
 */
void GPUHelper::GPU_DGEMM_2DTile_nn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8.,max_mapped_memory/num_gpus/8.,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      for (long int tk=0; tk<ntilesK; tk++){

          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
          cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
          for (long int j=0; j<tilesizesN[tn]; j++){
              F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
          }
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_nt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8.,max_mapped_memory/num_gpus/8.,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
          cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
          for (long int j=0; j<tilesizesN[tn]; j++){
              F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
          }
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8.,max_mapped_memory/num_gpus/8.,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
          cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
          for (long int j=0; j<tilesizesN[tn]; j++){
              F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
          }
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8.,max_mapped_memory/num_gpus/8.,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
          cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
          for (long int j=0; j<tilesizesN[tn]; j++){
              F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
          }
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}

void GPUHelper::Tiling(double mem1,double mem2,long int m,long int n,long int k){

  // first tile according to how much space is on gpu
  tilesizeN = n;
  tilesizeM = m;
  tilesizeK = k;
  ntilesM=ntilesN=ntilesK=1;
  while(tilesizeN*tilesizeM+tilesizeK*(tilesizeN+tilesizeM)>mem1){
     if (tilesizeN>tilesizeM){
        if (tilesizeN>tilesizeK){
           ntilesN++;
           tilesizeN = n/ntilesN;
           if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
        }
        else{
           ntilesK++;
           tilesizeK = k/ntilesK;
           if (k/ntilesK<(double)k/ntilesK) tilesizeK++;
        }
     }
     else{
        if (tilesizeM>tilesizeK){
           ntilesM++;
           tilesizeM = m/ntilesM;
           if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
        }
        else{
           ntilesK++;
           tilesizeK = k/ntilesK;
           if (k/ntilesK<(double)k/ntilesK) tilesizeK++;
        }
     }
  }

  // ensure each block of A, B, and C will fit in the temporary CPU buffer
  while(tilesizeN*tilesizeM>mem2){
     if (tilesizeN>tilesizeM){
        ntilesN++;
        tilesizeN = n/ntilesN;
        if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
     }
     else{
        ntilesM++;
        tilesizeM = m/ntilesM;
        if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
     }
  }
  while(tilesizeN*tilesizeK>mem2){
     if (tilesizeN>tilesizeK){
        ntilesN++;
        tilesizeN = n/ntilesN;
        if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
     }
     else{
        ntilesK++;
        tilesizeK = k/ntilesK;
        if (k/ntilesK<(double)k/ntilesK) tilesizeK++;
     }
  }
  while(tilesizeK*tilesizeM>mem2){
     if (tilesizeK>tilesizeM){
        ntilesK++;
        tilesizeK = k/ntilesK;
        if (k/ntilesK<(double)k/ntilesK) tilesizeK++;
     }
     else{
        ntilesM++;
        tilesizeM = m/ntilesM;
        if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
     }
  }

  // finally make sure that we've tiled enough so each gpu has something to do
  // TODO:  i need a way to make sure these end up 3:1, not 2:2
  //        ...should probably be more general than this...
  while(ntilesN*ntilesM<num_gpus){
     if (tilesizeN>tilesizeM){
        ntilesN++;
        tilesizeN = n/ntilesN;
        if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
     }
     else{
        ntilesM++;
        tilesizeM = m/ntilesM;
        if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
     }
  }

  lasttileN = n - (ntilesN-1)*tilesizeN;
  lasttileM = m - (ntilesM-1)*tilesizeM;
  lasttileK = k - (ntilesK-1)*tilesizeK;

  tilesizesM = (long int*)malloc(ntilesM*sizeof(long int));
  tilesizesN = (long int*)malloc(ntilesN*sizeof(long int));
  tilesizesK = (long int*)malloc(ntilesK*sizeof(long int));
  for (long int i=0; i<ntilesM-1; i++) tilesizesM[i] = tilesizeM;
  for (long int i=0; i<ntilesN-1; i++) tilesizesN[i] = tilesizeN;
  for (long int i=0; i<ntilesK-1; i++) tilesizesK[i] = tilesizeK;
  tilesizesM[ntilesM-1] = lasttileM;
  tilesizesN[ntilesN-1] = lasttileN;
  tilesizesK[ntilesK-1] = lasttileK;

}

}//end of namespace psi


