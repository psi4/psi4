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
#ifdef _OPENMP
  #include<omp.h>
#endif

// cuda libraries
#include<cuda.h>
#include<cublas.h>
#include<cuda_runtime.h>

#include"blas.h"
#include"gpuhelper.h"

using namespace psi;
using namespace boost;

namespace psi{

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
void GPUHelper::CudaInitGPU(Options&options){

  max_mapped_memory=0;
  num_gpus=gpumemory=extraroom=0;
  int n;
  cudaGetDeviceCount(&n);
  num_gpus = n;
  num_cpus=1;
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

     extraroom = 350L*1024L*1024L;
     
     cudaThreadExit();

     // default memory for mapped cpu memory is the sum of all gpu memory
     max_mapped_memory = (num_gpus+num_cpus) * (gpumemory-extraroom);
     if (options["MAX_MAPPED_MEMORY"].has_changed()){
        long int temp_mem = options.get_int("MAX_MAPPED_MEMORY");
        temp_mem *= 1024L*1024L;
        if (temp_mem<max_mapped_memory)
           max_mapped_memory = temp_mem;
     }
     max_mapped_memory_per_thread = max_mapped_memory/(num_gpus+num_cpus);

     fprintf(outfile,"\n");
     fprintf(outfile,"  allocating gpu memory...");
     fflush(outfile);
     tmp = (double**)malloc(num_gpus*sizeof(double*));
     gpuarray = (double**)malloc(num_gpus*sizeof(double*));
     #pragma omp parallel for schedule (static) num_threads(num_gpus)
     for (long int i=0; i<num_gpus; i++){
         long int thread = 0;
         #ifdef _OPENMP
           thread = omp_get_thread_num();
         #endif
         cudaSetDevice(thread);
         Check_CUDA_Error(stdout,"cudaSetDevice");
         cudaMallocHost((void**)&tmp[thread],max_mapped_memory_per_thread);  
         Check_CUDA_Error(outfile,"cpu tmp");
         cudaMalloc((void**)&gpuarray[thread],gpumemory-extraroom);
         Check_CUDA_Error(outfile,"gpu memory");
     }
     fprintf(outfile,"done.\n");
     fprintf(outfile,"\n");
     fflush(outfile);

     // some cpu memory for cores to use when stealing gpu work 
     cpuarray = (double**)malloc(num_cpus*sizeof(double*));
     for (long int i=0; i<num_cpus; i++){
         // TODO: need to be more intelligent about this...
         cpuarray[i] = (double*)malloc(3*max_mapped_memory_per_thread+20*max_mapped_memory_per_thread/30);
     }
  }
}
/*===================================================================

  free gpu and mapped cpu memory

===================================================================*/
void GPUHelper::CudaFinalizeGPU(Options&options){
  if (num_gpus>0){
     #pragma omp parallel for schedule (static) num_threads(num_gpus)
     for (long int i=0; i<num_gpus; i++){
         long int thread = 0;
         #ifdef _OPENMP
           thread = omp_get_thread_num();
         #endif
         cudaSetDevice(thread);
         Check_CUDA_Error(stdout,"cudaSetDevice (free)");
         cudaFreeHost(tmp[thread]);
         Check_CUDA_Error(outfile,"cpu tmp (free)");
         cudaFree(gpuarray[thread]);
         Check_CUDA_Error(outfile,"gpu memory (free)");
     }
     free(tmp);
     free(gpuarray);
     for (long int i=0; i<num_cpus; i++){
         free(cpuarray[i]);
     }
     free(cpuarray);
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
 * dgemm using a 2-dimensional tile - threaded versions for multiple gpus
 */
void GPUHelper::GPU_DGEMM_2DTile_nn_threaded_WithCpuStealing(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  TilingWithCpuStealing((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);
  //Tiling((gpumemory-extraroom)/8L,max_mapped_memory/num_gpus/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel num_threads(num_gpus+num_cpus)
  {

  long int thread = 0;
  #ifdef _OPENMP
    thread = omp_get_thread_num();
  #endif

  double*gpuA,*gpuB,*gpuC;
  // pointers to gpu memory
  if (thread<num_gpus){
     gpuA = gpuarray[thread];
     gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
     gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
  }
  // pointers to cpu memory
  else {
     gpuA = cpuarray[thread-num_gpus];
     gpuB = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK;
     gpuC = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK+tilesizeNprime*tilesizeK;
  }

  // cpu takes some of the 'N' tile
  if (StolenDimension=='N'){
     for (long int tm=0; tm<ntilesM; tm++){
         for (long int tk=0; tk<ntilesK; tk++){

             // this is for the gpus:
             if (thread<num_gpus){
                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
             
                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

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
             // this if for any cpu cores that might be helping:
             else{
                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,gpuA+i*tilesizesM[tm],1);
                }

                for (long int tn=0; tn<ntilesNprime; tn++){
                    if ((tm*ntilesNprime+tn)%num_cpus + num_gpus!=thread) continue;
                    for (long int i=0; i<tilesizesNprime[tn]; i++){
                        F_DCOPY(tilesizesK[tk],B+(NprimeOffSet+i+tn*tilesizeNprime)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(NprimeOffSet+j+tn*tilesizeNprime)*ldc+tm*tilesizeM,1);
                    }
                }
             }
         }
     }   
  }   
  // cpu takes some of the 'M' tile
  else if (StolenDimension=='M'){
     for (long int tn=0; tn<ntilesN; tn++){
         for (long int tk=0; tk<ntilesK; tk++){

             // this is for the gpus:
             if (thread<num_gpus){

                for (long int i=0; i<tilesizesN[tn]; i++){
                    F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
             
                for (long int tm=0; tm<ntilesM; tm++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             // this if for any cpu cores that might be helping:
             else{
                for (long int i=0; i<tilesizesN[tn]; i++){
                    F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                }
             
                for (long int tm=0; tm<ntilesMprime; tm++){
                    if ((tm*ntilesN+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesMprime[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeMprime+MprimeOffSet,1,gpuA+i*tilesizesMprime[tm],1);
                    }
                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesMprime[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
                    }
                }
             }
         }
     }   
  }
  else{
     if (thread<num_gpus){
        for (long int tm=0; tm<ntilesM; tm++){
            for (long int tk=0; tk<ntilesK; tk++){

                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                
                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

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
        }
     }
  }

  }
  free(tilesizesMprime);
  free(tilesizesNprime);
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
/**
 * dgemm using a 2-dimensional tile - threaded versions for multiple gpus
 */
void GPUHelper::GPU_DGEMM_2DTile_nn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){

          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_nn(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc,int thread){

  TilingNoThread((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  for (long int mn=0; mn<ntilesM*ntilesN; mn++){

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){

          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_nt_threaded_WithCpuStealing(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  //TilingWithCpuStealing((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);
  TilingWithCpuStealing((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;


  #pragma omp parallel num_threads(num_gpus+num_cpus)
  {

  long int thread = 0;
  #ifdef _OPENMP
    thread = omp_get_thread_num();
  #endif

  double*gpuA,*gpuB,*gpuC;

  // pointers to gpu memory
  if (thread<num_gpus){
     gpuA = gpuarray[thread];
     gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
     gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
  }
  // pointers to cpu memory
  else {
     gpuA = cpuarray[thread-num_gpus];
     gpuB = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK;
     gpuC = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK+tilesizeNprime*tilesizeK;
  }

  // cpu takes some of the 'N' tile
  if (StolenDimension=='N'){
     for (long int tm=0; tm<ntilesM; tm++){
         for (long int tk=0; tk<ntilesK; tk++){
             if (thread<num_gpus){

                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

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
             else{

                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,gpuA+i*tilesizesM[tm],1);
                }

                for (long int tn=0; tn<ntilesNprime; tn++){
                    if ((tm*ntilesNprime+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesNprime[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeNprime+NprimeOffSet,1,gpuB+i*tilesizesNprime[tn],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesNprime[tn],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(j+tn*tilesizeNprime+NprimeOffSet)*ldc+tm*tilesizeM,1);
                    }
                }

             }
         }
     }
  }
  else if (StolenDimension=='M'){
     for (long int tn=0; tn<ntilesN; tn++){
         for (long int tk=0; tk<ntilesK; tk++){
             if (thread<num_gpus){

                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tm=0; tm<ntilesM; tm++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }

             }
             else{

                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,gpuB+i*tilesizesN[tn],1);
                }

                for (long int tm=0; tm<ntilesMprime; tm++){
                    if ((tm*ntilesN+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesMprime[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeMprime+MprimeOffSet,1,gpuA+i*tilesizesMprime[tm],1);
                    }

                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesMprime[tm],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
                    }
                }

             }
         }
     }
  }
  else{
     for (long int tm=0; tm<ntilesM; tm++){
         for (long int tn=0; tn<ntilesN; tn++){
             if (thread<num_gpus){
                if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
                for (long int tk=0; tk<ntilesK; tk++){
                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
                }
                cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                for (long int j=0; j<tilesizesN[tn]; j++){
                    F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                }

             }
         }
     }
  }
  }

  free(tilesizesNprime);
  free(tilesizesMprime);
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_nt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_nt(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc,int thread){

  TilingNoThread((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  for (long int mn=0; mn<ntilesM*ntilesN; mn++){

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tn_threaded_WithCpuStealing(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  TilingWithCpuStealing((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel num_threads(num_gpus+num_cpus)
  {

  long int thread = 0;
  #ifdef _OPENMP
    thread = omp_get_thread_num();
  #endif

  double*gpuA,*gpuB,*gpuC;

  // pointers to gpu memory
  if (thread<num_gpus){
     gpuA = gpuarray[thread];
     gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
     gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
  }
  // pointers to cpu memory
  else {
     gpuA = cpuarray[thread-num_gpus];
     gpuB = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK;
     gpuC = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK+tilesizeNprime*tilesizeK;
  }

  // cpu takes some of the 'N' tile
  StolenDimension=' ';
  if (StolenDimension=='N'){
     for (long int tm=0; tm<ntilesM; tm++){
         for (long int tk=0; tk<ntilesK; tk++){
             if (thread<num_gpus){
                for (long int i=0; i<tilesizesM[tm]; i++){
                    F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tn=0; tn<ntilesN; tn++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

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
             else{
                for (long int i=0; i<tilesizesM[tm]; i++){
                    F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                }

                for (long int tn=0; tn<ntilesNprime; tn++){

                    if ((tm*ntilesNprime+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesNprime[tn]; i++){
                        F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeNprime+NprimeOffSet)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(j+tn*tilesizeNprime+NprimeOffSet)*ldc+tm*tilesizeM,1);
                    }
                }
             }
         }
     }
  }
  else if (StolenDimension=='M'){
     for (long int tn=0; tn<ntilesN; tn++){
         for (long int tk=0; tk<ntilesK; tk++){
             if (thread<num_gpus){
                for (long int i=0; i<tilesizesN[tn]; i++){
                    F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tm=0; tm<ntilesM; tm++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesM[tm]; i++){
                        F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             else{
                for (long int i=0; i<tilesizesN[tn]; i++){
                    F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                }

                for (long int tm=0; tm<ntilesMprime; tm++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesMprime[tm]; i++){
                        F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeMprime+MprimeOffSet)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
                    }
                }
             }
         }
     }
  }
  else{
     if (thread<num_gpus){
        for (long int tm=0; tm<ntilesM; tm++){
            for (long int tn=0; tn<ntilesN; tn++){
                cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
                for (long int tk=0; tk<ntilesK; tk++){
                    for (long int i=0; i<tilesizesM[tm]; i++){
                        F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    for (long int i=0; i<tilesizesN[tn]; i++){
                        F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
                }
                cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                for (long int j=0; j<tilesizesN[tn]; j++){
                    F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                }

            }
        }
     }
  }


  }
  free(tilesizesMprime);
  free(tilesizesNprime);
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tn(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc,int thread){

  TilingNoThread((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  for (long int mn=0; mn<ntilesM*ntilesN; mn++){

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              F_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tt_threaded_WithCpuStealing(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  TilingWithCpuStealing((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;


  #pragma omp parallel num_threads(num_gpus+num_cpus)
  {

  long int thread = 0;
  #ifdef _OPENMP
    thread = omp_get_thread_num();
  #endif

  double*gpuA,*gpuB,*gpuC;

  // pointers to gpu memory
  if (thread<num_gpus){
     gpuA = gpuarray[thread];
     gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
     gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
  }
  // pointers to cpu memory
  else {
     gpuA = cpuarray[thread-num_gpus];
     gpuB = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK;
     gpuC = cpuarray[thread-num_gpus]+tilesizeMprime*tilesizeK+tilesizeNprime*tilesizeK;
  }

  // cpu takes some of the 'N' tile
  StolenDimension=' ';
  if (StolenDimension=='N'){
     for (long int tm=0; tm<ntilesM; tm++){
         for (long int tk=0; tk<ntilesK; tk++){
             if (thread<num_gpus){
                for (long int i=0; i<tilesizesM[tm]; i++){
                    F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

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
             else{
                for (long int i=0; i<tilesizesM[tm]; i++){
                    F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                }
                for (long int tn=0; tn<ntilesNprime; tn++){
                    if ((tm*ntilesNprime+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesNprime[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeNprime+NprimeOffSet,1,gpuB+i*tilesizesN[tn],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesNprime[tn],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(j+tn*tilesizeNprime+NprimeOffSet)*ldc+tm*tilesizeM,1);
                    }
                }
             }
         }
     }
  }
  else if (StolenDimension=='M'){
     for (long int tn=0; tn<ntilesN; tn++){
         for (long int tk=0; tk<ntilesK; tk++){
             if (thread<num_gpus){
                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                for (long int tm=0; tm<ntilesM; tm++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesM[tm]; i++){
                        F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             else{
                for (long int i=0; i<tilesizesK[tk]; i++){
                    F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,gpuB+i*tilesizesN[tn],1);
                }
                for (long int tm=0; tm<ntilesMprime; tm++){

                    if ((tm*ntilesN+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesMprime[tm]; i++){
                        F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeMprime+MprimeOffSet)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        F_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
                    }
                }
             }
         }
     }
  }
  else{
     if (thread<num_gpus){
        for (long int tm=0; tm<ntilesM; tm++){
            for (long int tn=0; tn<ntilesN; tn++){
                cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
                for (long int tk=0; tk<ntilesK; tk++){
                    for (long int i=0; i<tilesizesM[tm]; i++){
                        F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    for (long int i=0; i<tilesizesK[tk]; i++){
                        F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
                }
                cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                for (long int j=0; j<tilesizesN[tn]; j++){
                    F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                }
            }
        }
     }
  }

  }

  free(tilesizesMprime);
  free(tilesizesNprime);
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}
void GPUHelper::GPU_DGEMM_2DTile_tt(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc,int thread){

  TilingNoThread((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  for (long int mn=0; mn<ntilesM*ntilesN; mn++){

      // pointers to gpu memory
      double*gpuA = gpuarray[thread];
      double*gpuB = gpuarray[thread]+tilesizeM*tilesizeK;
      double*gpuC = gpuarray[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              F_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              F_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          F_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  }
  free(tilesizesM);
  free(tilesizesN);
  free(tilesizesK);
}

void GPUHelper::TilingNoThread(long int mem1,long int mem2,long int m,long int n,long int k){

  // first tile according to how much space is on gpu
  tilesizeN = n;
  tilesizeM = m;
  tilesizeK = k;
  ntilesM=ntilesN=ntilesK=1L;
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

  lasttileN = n - (ntilesN-1L)*tilesizeN;
  lasttileM = m - (ntilesM-1L)*tilesizeM;
  lasttileK = k - (ntilesK-1L)*tilesizeK;

  tilesizesM = (long int*)malloc(ntilesM*sizeof(long int));
  tilesizesN = (long int*)malloc(ntilesN*sizeof(long int));
  tilesizesK = (long int*)malloc(ntilesK*sizeof(long int));
  for (long int i=0; i<ntilesM-1L; i++) tilesizesM[i] = tilesizeM;
  for (long int i=0; i<ntilesN-1L; i++) tilesizesN[i] = tilesizeN;
  for (long int i=0; i<ntilesK-1L; i++) tilesizesK[i] = tilesizeK;
  tilesizesM[ntilesM-1L] = lasttileM;
  tilesizesN[ntilesN-1L] = lasttileN;
  tilesizesK[ntilesK-1L] = lasttileK;

  //printf("%5li %5li %5li (n^5)\n",ntilesM,ntilesN,ntilesK);fflush(stdout);

}
void GPUHelper::Tiling(long int mem1,long int mem2,long int m,long int n,long int k){

  // first tile according to how much space is on gpu
  tilesizeN = n;
  tilesizeM = m;
  tilesizeK = k;
  ntilesM=ntilesN=ntilesK=1L;
  while(tilesizeN*tilesizeM+tilesizeK*(tilesizeN+tilesizeM)>mem1){
     if (ntilesN*ntilesM<num_gpus){
        if (tilesizeN>tilesizeM){
           ntilesN++;
           ntilesN++;
           tilesizeN = n/ntilesN;
           if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
        }
        else{
           ntilesM++;
           ntilesM++;
           tilesizeM = m/ntilesM;
           if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
        }
     }
     else{
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
  }

  // ensure each block of A, B, and C will fit in the temporary CPU buffer
  while(tilesizeN*tilesizeM>mem2){
     if (ntilesN*ntilesM<num_gpus){
        if (tilesizeN>tilesizeM){
           ntilesN++;
           ntilesN++;
           tilesizeN = n/ntilesN;
           if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
        }
        else{
           ntilesM++;
           ntilesM++;
           tilesizeM = m/ntilesM;
           if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
        }
     }
     else{
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
  }

  while(tilesizeN*tilesizeK>mem2){
     if (ntilesN*ntilesM<num_gpus){
        ntilesN++;
        ntilesN++;
        tilesizeN = n/ntilesN;
        if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
     }
     else{
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
  }
  while(tilesizeK*tilesizeM>mem2){
     if (ntilesN*ntilesM<num_gpus){
        ntilesM++;
        ntilesM++;
        tilesizeM = m/ntilesM;
        if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
     }
     else{
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
  }

  // finally make sure that we've tiled enough so each gpu has something to do
  // TODO:  i need a way to make sure these end up 3:1, not 2:2
  //        ...should probably be more general than this...
  while(ntilesN*ntilesM<num_gpus){
     if (tilesizeN>tilesizeM){
        ntilesN++;
        ntilesN++;
        tilesizeN = n/ntilesN;
        if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
     }
     else{
        ntilesM++;
        ntilesM++;
        tilesizeM = m/ntilesM;
        if (m/ntilesM<(double)m/ntilesM) tilesizeM++;
     }
  }

  lasttileN = n - (ntilesN-1L)*tilesizeN;
  lasttileM = m - (ntilesM-1L)*tilesizeM;
  lasttileK = k - (ntilesK-1L)*tilesizeK;

  tilesizesM = (long int*)malloc(ntilesM*sizeof(long int));
  tilesizesN = (long int*)malloc(ntilesN*sizeof(long int));
  tilesizesK = (long int*)malloc(ntilesK*sizeof(long int));
  for (long int i=0; i<ntilesM-1L; i++) tilesizesM[i] = tilesizeM;
  for (long int i=0; i<ntilesN-1L; i++) tilesizesN[i] = tilesizeN;
  for (long int i=0; i<ntilesK-1L; i++) tilesizesK[i] = tilesizeK;
  tilesizesM[ntilesM-1L] = lasttileM;
  tilesizesN[ntilesN-1L] = lasttileN;
  tilesizesK[ntilesK-1L] = lasttileK;

  //printf("%5li %5li %5li\n",ntilesM,ntilesN,ntilesK);fflush(stdout);

}
void GPUHelper::TilingWithCpuStealing(long int mem1,long int mem2,long int m,long int n,long int k){
  // compute normal tiling
  Tiling(mem1,mem2,m,n,k);
  // take a slice of the larger of m or n for the cpu
  
  // first let's just try taking a sliver of the last tile of N...
  ntilesNprime = num_cpus;
  ntilesMprime = num_cpus;
  tilesizesNprime = (long int*)malloc(ntilesNprime*sizeof(long int));
  tilesizesMprime = (long int*)malloc(ntilesMprime*sizeof(long int));

  // which dimension will cpu work on?
  if (tilesizeN>tilesizeM){
     // assume the gpu is ~30x faster than a single core:
     if (tilesizeN<30){
        StolenDimension = ' ';
        return;
     }
     tilesizeNprime = tilesizeN/30;

     StolenDimension = 'N';

     // need to figure out new tiles in N (and might as well make them even)
     // TODO: should make these multiples of the warp size, too
     long int newn = n-num_cpus*tilesizeNprime;
     tilesizeN = newn/ntilesN-1;
     lasttileN = tilesizeN;
     for (long int i=0; i<ntilesN; i++)
         tilesizesN[i] = tilesizeN;

     // redo Nprime's numbers
     ntilesNprime   = num_cpus;
     NprimeOffSet   = ntilesN*tilesizeN;
     tilesizeNprime = (n - NprimeOffSet)/ntilesNprime;
     if (tilesizeNprime*ntilesNprime<(n-NprimeOffSet)) tilesizeNprime++;
     lasttileNprime = (n - NprimeOffSet)-(ntilesNprime-1)*tilesizeNprime;

     for (long int i=0; i<ntilesNprime-1; i++) tilesizesNprime[i] = tilesizeNprime;
     tilesizesNprime[ntilesNprime-1] = lasttileNprime;

     // set this just for memory mapping
     lasttileMprime = 0;
     tilesizesMprime[0] = lasttileMprime;
     tilesizeMprime = tilesizeM;
  }
  // do M instead:
  else{
     // assume the gpu is ~30x faster than a single core:
     if (tilesizeM<30){
        StolenDimension = ' ';
        return;
     }
     tilesizeMprime = tilesizeM/30;

     StolenDimension = 'M';

     // need to figure out new tiles in N (and might as well make them even)
     // TODO: should make these multiples of the warp size, too
     long int newm = m-num_cpus*tilesizeMprime;
     tilesizeM = newm/ntilesM-1;
     lasttileM = tilesizeM;
     for (long int i=0; i<ntilesM; i++)
         tilesizesM[i] = tilesizeM;

     // redo Mprime's numbers
     ntilesMprime   = num_cpus;
     MprimeOffSet   = ntilesM*tilesizeM;
     tilesizeMprime = (m - MprimeOffSet)/ntilesMprime;
     if (tilesizeMprime*ntilesMprime<(m-MprimeOffSet)) tilesizeMprime++;
     lasttileMprime = (m - MprimeOffSet)-(ntilesMprime-1)*tilesizeMprime;

     for (long int i=0; i<ntilesMprime-1; i++) tilesizesMprime[i] = tilesizeMprime;
     tilesizesMprime[ntilesMprime-1] = lasttileMprime;

     // set this just for memory mapping
     lasttileNprime = 0;
     tilesizesNprime[0] = lasttileNprime;
     tilesizeNprime = tilesizeN;

     //printf("hey the tile is %5li (%5li) %5li %5li\n",tilesizeMprime,m,ntilesM,tilesizeM);fflush(stdout);
  }
}

}//end of namespace psi


