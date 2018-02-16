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

// TODO: interleaved dgemm seems to be broken
bool interleaved_dgemm = true;

#include<psi4/libplugin/plugin.h>
#include<psi4/psi4-dec.h>
#include"psi4/liboptions/liboptions.h"
#include<psi4/libqt/qt.h>
#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>
#include<psi4/libmints/matrix.h>
#include<psi4/libmints/vector.h>
#include<psi4/libiwl/iwl.h>
#include<psi4/libpsio/psio.hpp>
#include<psi4/libparallel/process.h>
#include"blas.h"

#include"gpuhelper.h"
#include"gpuonly.h"
#include<omp.h>
using namespace psi;
using namespace std;

namespace psi{namespace fnocc{

void GPUHelper::Check_CUDA_Error(FILE*fp,const char *message){
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
  size_t free;
  size_t total;
  cudaGetDeviceCount(&n);
  num_gpus = n;
  num_cpus=0;
  if (options["NUM_GPUS"].has_changed())
     num_gpus = options.get_int("NUM_GPUS");

  if (num_gpus>0){
     cublasInit();
     struct cudaDeviceProp cudaProp;
     int gpu_id;
     cudaGetDevice(&gpu_id);
     cudaGetDeviceProperties( &cudaProp,gpu_id );
     outfile->Printf(
       "\n  _________________________________________________________\n");
     outfile->Printf("  CUDA device properties:\n");
     outfile->Printf("  name:                 %20s\n",cudaProp.name);
     outfile->Printf("  major version:        %20d\n",cudaProp.major);
     outfile->Printf("  minor version:        %20d\n",cudaProp.minor);
     outfile->Printf("  canMapHostMemory:     %20d\n",cudaProp.canMapHostMemory);
     outfile->Printf("  totalGlobalMem:       %20lu mb\n",
       cudaProp.totalGlobalMem/(1024*1024));
     outfile->Printf("  sharedMemPerBlock:    %20lu\n",cudaProp.sharedMemPerBlock);
     outfile->Printf("  clockRate:            %20.3f ghz\n",
       cudaProp.clockRate/1.0e6);
     outfile->Printf("  regsPerBlock:         %20d\n",cudaProp.regsPerBlock);
     outfile->Printf("  warpSize:             %20d\n",cudaProp.warpSize);
     outfile->Printf("  maxThreadsPerBlock:   %20d\n",cudaProp.maxThreadsPerBlock);
     outfile->Printf(
       "  _________________________________________________________\n\n");
     //fflush(outfile);

     //gpumemory = cudaProp.totalGlobalMem;
     
     
     cudaMemGetInfo(&free,&total);
     gpumemory = free; 
     extraroom = 200L*1024L*1024L;
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

     outfile->Printf("\n");
     outfile->Printf("  allocating gpu memory...");
     //fflush(outfile);
     tmp = (double**)malloc(num_gpus*sizeof(double*));   
     gpubuffer = (double**)malloc(num_gpus*sizeof(double*));
     #pragma omp parallel for schedule (static) num_threads(num_gpus)
     for (long int i=0; i<num_gpus; i++){
         long int thread = 0;
         #ifdef _OPENMP
           thread = omp_get_thread_num();
         #endif
         cudaSetDevice(thread);
         Check_CUDA_Error(stdout,"cudaSetDevice");
         cudaMallocHost((void**)&tmp[thread],max_mapped_memory_per_thread);  
         //tmp[thread] = (double*)malloc(max_mapped_memory_per_thread*sizeof(double));
         Check_CUDA_Error(stdout,"cpu tmp");
         //cudaMemGetInfo(&free,&total);
         cudaMalloc((void**)&gpubuffer[thread],gpumemory-extraroom);
    //     cudaMalloc((void**)&gpubuffer[thread],gpumemory-extraroom);   
         Check_CUDA_Error(stdout,"gpu memory");

     }
     // thread-safe tiling info: TODO: these are never free'd at the end
     myntilesM = (long int*)malloc(num_gpus*sizeof(long int));
     myntilesN = (long int*)malloc(num_gpus*sizeof(long int));
     myntilesK = (long int*)malloc(num_gpus*sizeof(long int));
     mytilesizeM = (long int*)malloc(num_gpus*sizeof(long int));
     mytilesizeN = (long int*)malloc(num_gpus*sizeof(long int));
     mytilesizeK = (long int*)malloc(num_gpus*sizeof(long int));
     mylasttileM = (long int*)malloc(num_gpus*sizeof(long int));
     mylasttileN = (long int*)malloc(num_gpus*sizeof(long int));
     mylasttileK = (long int*)malloc(num_gpus*sizeof(long int));
     mytilesizesM = (long int**)malloc(num_gpus*sizeof(long int*));
     mytilesizesN = (long int**)malloc(num_gpus*sizeof(long int*));
     mytilesizesK = (long int**)malloc(num_gpus*sizeof(long int*));

     //fflush(outfile);

     // some cpu memory for cores to use when stealing gpu work 
     //cpuarray = (double**)malloc(num_cpus*sizeof(double*));
     //for (long int i=0; i<num_cpus; i++){
     //    // TODO: need to be more intelligent about this...
     //    cpuarray[i] = (double*)malloc(3*max_mapped_memory_per_thread+20*max_mapped_memory_per_thread/30);
     //}
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
         Check_CUDA_Error(stdout,"cpu tmp (free)");
         cudaFree(gpubuffer[thread]);
         Check_CUDA_Error(stdout,"gpu memory (free)");
     }
     free(tmp);
     free(gpubuffer);
     //for (long int i=0; i<num_cpus; i++){
     //    free(cpuarray[i]);
     //}
     //free(cpuarray);
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

  throw PsiException("GPU_DGEMM_2DTile_nn_threaded_WithCpuStealing: not implemented",__FILE__,__LINE__);
//DPG commented out to remove statement unreachable warning
/*
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
     cudaSetDevice(thread);
     gpuA = gpubuffer[thread];
     gpuB = gpubuffer[thread]+tilesizeM*tilesizeK;
     gpuC = gpubuffer[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
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
                    C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
             
                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesN[tn]; i++){
                        C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             // this if for any cpu cores that might be helping:
             else{
                for (long int i=0; i<tilesizesK[tk]; i++){
                    C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,gpuA+i*tilesizesM[tm],1);
                }

                for (long int tn=0; tn<ntilesNprime; tn++){
                    if ((tm*ntilesNprime+tn)%num_cpus + num_gpus!=thread) continue;
                    for (long int i=0; i<tilesizesNprime[tn]; i++){
                        C_DCOPY(tilesizesK[tk],B+(NprimeOffSet+i+tn*tilesizeNprime)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(NprimeOffSet+j+tn*tilesizeNprime)*ldc+tm*tilesizeM,1);
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
                    C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
             
                for (long int tm=0; tm<ntilesM; tm++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             // this if for any cpu cores that might be helping:
             else{
                for (long int i=0; i<tilesizesN[tn]; i++){
                    C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                }
             
                for (long int tm=0; tm<ntilesMprime; tm++){
                    if ((tm*ntilesN+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesMprime[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeMprime+MprimeOffSet,1,gpuA+i*tilesizesMprime[tm],1);
                    }
                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesMprime[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
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
                    C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                
                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesN[tn]; i++){
                        C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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
*/
}

/**
 * dgemm using a 2-dimensional tile - threaded versions for multiple gpus
 */
void report_num_threads(int level)
{
    #pragma omp single
    {
        printf("Level %d: number of threads in the team - %d\n",
                  level, omp_get_num_threads());
    }
 }
void GPUHelper::GPU_DGEMM_2DTile_nn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  omp_set_nested(1);
  omp_set_dynamic(0);
  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      cudaSetDevice(thread);

      // pointers to gpu memory
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+tilesizeM*tilesizeK*2;
      double*gpuC = gpubuffer[thread]+tilesizeM*tilesizeK*2+tilesizeN*tilesizeK*2;

      long int offsetA = tilesizeM * tilesizeK;
      long int offsetB = tilesizeN * tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;
      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
      omp_set_nested(1);
      omp_set_dynamic(0);
if (interleaved_dgemm) {
      // create streams:
      cudaStream_t stream1;
      cudaStreamCreate(&stream1);
      cudaEvent_t estart1,estop1;
      cudaEventCreate(&estart1);
      cudaEventCreate(&estop1);
      cublasSetKernelStream(stream1);

      cudaStream_t stream2;
      cudaStreamCreate(&stream2);
      cudaEvent_t estart2,estop2;
      cudaEventCreate(&estart2);
      cudaEventCreate(&estop2);

      double start = omp_get_wtime();

      // need to transfer data for first tile
      for (long int i=0; i<tilesizesK[0]; i++){
          C_DCOPY(tilesizesM[tm],A+(i+0*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
      }
      cudaMemcpyAsync(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice,stream1);
      cudaStreamSynchronize(stream1);
      for (long int i=0; i<tilesizesN[tn]; i++){
          C_DCOPY(tilesizesK[0],B+(i+tn*tilesizeN)*ldb+0*tilesizeK,1,tmp[thread]+i*tilesizesK[0],1);
      }
      cudaMemcpyAsync(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice,stream1);
      cudaStreamSynchronize(stream1);

      for (long int tk=0; tk<ntilesK; tk++){

          #pragma omp parallel num_threads(2)
          {

              long int thread2 = omp_get_thread_num();
              if (thread2 == 0) {

                  double * A_curr = ( tk % 2 == 0 ) ? gpuA : gpuA + offsetA;
                  double * B_curr = ( tk % 2 == 0 ) ? gpuB : gpuB + offsetB;

                  cudaEventRecord(estart1,stream1);
                      cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,A_curr,tilesizesM[tm],B_curr,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
                      cudaStreamSynchronize(stream1);
                  cudaEventRecord(estop1,stream1);

              } else {
                  // only copy next tiles if we need them:
                  if ( tk < ntilesK - 1 ) {
                      double * A_next = ( tk % 2 == 0 ) ? gpuA + offsetA : gpuA;
                      double * B_next = ( tk % 2 == 0 ) ? gpuB + offsetB : gpuB;
                      cudaEventRecord(estart2,stream2);
                          for (long int i=0; i<tilesizesK[tk+1]; i++){
                              C_DCOPY(tilesizesM[tm],A+(i+(tk+1)*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                          }
                          cudaMemcpyAsync(A_next,tmp[thread],tilesizesM[tm]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                          for (long int i=0; i<tilesizesN[tn]; i++){
                              C_DCOPY(tilesizesK[tk+1],B+(i+tn*tilesizeN)*ldb+(tk+1)*tilesizeK,1,tmp[thread]+i*tilesizesK[tk+1],1);
                          }
                          cudaMemcpyAsync(B_next,tmp[thread],tilesizesN[tn]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                      cudaEventRecord(estop2,stream2);
                  }
              }
          }
          cudaThreadSynchronize();
      }
      cublasSetKernelStream(NULL);
      cudaEventDestroy(estart2);
      cudaEventDestroy(estart1);
      cudaEventDestroy(estop1);
      cudaEventDestroy(estop2);
      cudaStreamDestroy(stream1);
      cudaStreamDestroy(stream2);
}else {
      // original version:
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesK[tk]; i++){
              C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
      }
}
      omp_set_nested(0);
      omp_set_dynamic(1);
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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

  cudaSetDevice(thread);

  for (long int mn=0; mn<myntilesM[thread]*myntilesN[thread]; mn++){

      // pointers to gpu memory
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread];
      double*gpuC = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread]+mytilesizeN[thread]*mytilesizeK[thread];

      long int tn = mn%myntilesN[thread];
      long int tm = (mn-tn)/myntilesN[thread];

      cudaMemset((void*)gpuC,'\0',mytilesizesM[thread][tm]*mytilesizesN[thread][tn]*sizeof(double));
      for (long int tk=0; tk<myntilesK[thread]; tk++){

          for (long int i=0; i<mytilesizesK[thread][tk]; i++){
              C_DCOPY(mytilesizesM[thread][tm],A+(i+tk*mytilesizeK[thread])*lda+tm*mytilesizeM[thread],1,tmp[thread]+i*mytilesizesM[thread][tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],mytilesizesM[thread][tm]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<mytilesizesN[thread][tn]; i++){
              C_DCOPY(mytilesizesK[thread][tk],B+(i+tn*mytilesizeN[thread])*ldb+tk*mytilesizeK[thread],1,tmp[thread]+i*mytilesizesK[thread][tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],mytilesizesN[thread][tn]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,mytilesizesM[thread][tm],mytilesizesN[thread][tn],mytilesizesK[thread][tk],alpha,gpuA,mytilesizesM[thread][tm],gpuB,mytilesizesK[thread][tk],1.0,gpuC,mytilesizesM[thread][tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,mytilesizesN[thread][tn]*mytilesizesM[thread][tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<mytilesizesN[thread][tn]; j++){
          C_DAXPY(mytilesizesM[thread][tm],1.0,tmp[thread]+j*mytilesizesM[thread][tm],1,C+(j+tn*mytilesizeN[thread])*ldc+tm*mytilesizeM[thread],1);
      }
  }
  free(mytilesizesM[thread]);
  free(mytilesizesN[thread]);
  free(mytilesizesK[thread]);
}
void GPUHelper::GPU_DGEMM_2DTile_nt_threaded_WithCpuStealing(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){


  throw PsiException("GPU_DGEMM_2DTile_nt_threaded_WithCpuStealing: not implemented",__FILE__,__LINE__);
//DPG commented out to remove statement unreachable warning
/*
  //Tiling((gpumemory-extraroom)/8L,max_mapped_memory/num_gpus/8L,m,n,k);
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
     cudaSetDevice(thread);
     gpuA = gpubuffer[thread];
     gpuB = gpubuffer[thread]+tilesizeM*tilesizeK;
     gpuC = gpubuffer[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
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
                    C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
                   cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                   for (long int j=0; j<tilesizesN[tn]; j++){
                       C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                   }
                }

             }
             else{

                for (long int i=0; i<tilesizesK[tk]; i++){
                    C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,gpuA+i*tilesizesM[tm],1);
                }

                for (long int tn=0; tn<ntilesNprime; tn++){
                    if ((tm*ntilesNprime+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesNprime[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeNprime+NprimeOffSet,1,gpuB+i*tilesizesNprime[tn],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesNprime[tn],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(j+tn*tilesizeNprime+NprimeOffSet)*ldc+tm*tilesizeM,1);
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
                    C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tm=0; tm<ntilesM; tm++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }

             }
             else{

                for (long int i=0; i<tilesizesK[tk]; i++){
                    C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,gpuB+i*tilesizesN[tn],1);
                }

                for (long int tm=0; tm<ntilesMprime; tm++){
                    if ((tm*ntilesN+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesMprime[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeMprime+MprimeOffSet,1,gpuA+i*tilesizesMprime[tm],1);
                    }

                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesMprime[tm],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
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
                        C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
                }
                cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                for (long int j=0; j<tilesizesN[tn]; j++){
                    C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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
*/
}
void GPUHelper::GPU_DGEMM_2DTile_nt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;


  omp_set_nested(1);
  omp_set_dynamic(0);
  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      cudaSetDevice(thread);

      // pointers to gpu memory ... keep in mind that tilesizeK has been reduced by at least a factor of 2.
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+tilesizeM*tilesizeK*2;
      double*gpuC = gpubuffer[thread]+tilesizeM*tilesizeK*2+tilesizeN*tilesizeK*2;

      long int offsetA = tilesizeM * tilesizeK;
      long int offsetB = tilesizeN * tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;
      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));

      omp_set_nested(1);
      omp_set_dynamic(0);
if (interleaved_dgemm) {
      // create streams:
      cudaStream_t stream1;
      cudaStreamCreate(&stream1);
      cudaEvent_t estart1,estop1;
      cudaEventCreate(&estart1);
      cudaEventCreate(&estop1);
      cublasSetKernelStream(stream1);

      cudaStream_t stream2;
      cudaStreamCreate(&stream2);
      cudaEvent_t estart2,estop2;
      cudaEventCreate(&estart2);
      cudaEventCreate(&estop2);

      double start = omp_get_wtime();

      // need to transfer data for first tile
      for (long int i=0; i<tilesizesK[0]; i++){
          C_DCOPY(tilesizesM[tm],A+(i+0*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
      }
      cudaMemcpyAsync(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice,stream1);
      cudaStreamSynchronize(stream1);
     for (long int i=0; i<tilesizesK[0]; i++){
          C_DCOPY(tilesizesN[tn],B+(i+0*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
      }
      cudaMemcpyAsync(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice,stream1);
      cudaStreamSynchronize(stream1);
      for (long int tk=0; tk<ntilesK; tk++){
	#pragma omp parallel num_threads(2)
	{
              int thread2 = omp_get_thread_num();
              if (thread2 == 0) {

                  double * A_curr = ( tk % 2 == 0 ) ? gpuA : gpuA + offsetA;
                  double * B_curr = ( tk % 2 == 0 ) ? gpuB : gpuB + offsetB;

                  cudaEventRecord(estart1,stream1);
                      cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,A_curr,tilesizesM[tm],B_curr,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
                      cudaStreamSynchronize(stream1);
                  cudaEventRecord(estop1,stream1);

              } else {
                  // only copy next tiles if we need them:
                  if ( tk < ntilesK - 1) {
                      double * A_next = ( tk % 2 == 0 ) ? gpuA + offsetA : gpuA;
                      double * B_next = ( tk % 2 == 0 ) ? gpuB + offsetB : gpuB;
                      cudaEventRecord(estart2,stream2);
                          for (long int i=0; i<tilesizesK[tk+1]; i++){
                              C_DCOPY(tilesizesM[tm],A+(i+(tk+1)*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
                          }
                          cudaMemcpyAsync(A_next,tmp[thread],tilesizesM[tm]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                          for (long int i=0; i<tilesizesK[(tk+1)]; i++){
                              C_DCOPY(tilesizesN[tn],B+(i+(tk+1)*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                          }
                          cudaMemcpyAsync(B_next,tmp[thread],tilesizesN[tn]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                      cudaEventRecord(estop2,stream2);
                  }
              }
          }
          cudaThreadSynchronize();
      }
      cublasSetKernelStream(NULL);
      cudaEventDestroy(estart2);
      cudaEventDestroy(estart1);
      cudaEventDestroy(estop1);
      cudaEventDestroy(estop2);
      cudaStreamDestroy(stream1);
      cudaStreamDestroy(stream2);
}else {
      // original version:
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesK[tk]; i++){
              C_DCOPY(tilesizesM[tm],A+(i+tk*tilesizeK)*lda+tm*tilesizeM,1,tmp[thread]+i*tilesizesM[tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesM[tm],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
      }
}
      omp_set_nested(0);
      omp_set_dynamic(1);
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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

  cudaSetDevice(thread);

  for (long int mn=0; mn<myntilesM[thread]*myntilesN[thread]; mn++){

      // pointers to gpu memory
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread];
      double*gpuC = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread]+mytilesizeN[thread]*mytilesizeK[thread];

      long int tn = mn%myntilesN[thread];
      long int tm = (mn-tn)/myntilesN[thread];

      cudaMemset((void*)gpuC,'\0',mytilesizesM[thread][tm]*mytilesizesN[thread][tn]*sizeof(double));
      for (long int tk=0; tk<myntilesK[thread]; tk++){
          for (long int i=0; i<mytilesizesK[thread][tk]; i++){
              C_DCOPY(mytilesizesM[thread][tm],A+(i+tk*mytilesizeK[thread])*lda+tm*mytilesizeM[thread],1,tmp[thread]+i*mytilesizesM[thread][tm],1);
          }
          cudaMemcpy(gpuA,tmp[thread],mytilesizesM[thread][tm]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<mytilesizesK[thread][tk]; i++){
              C_DCOPY(mytilesizesN[thread][tn],B+(i+tk*mytilesizeK[thread])*ldb+tn*mytilesizeN[thread],1,tmp[thread]+i*mytilesizesN[thread][tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],mytilesizesN[thread][tn]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,mytilesizesM[thread][tm],mytilesizesN[thread][tn],mytilesizesK[thread][tk],alpha,gpuA,mytilesizesM[thread][tm],gpuB,mytilesizesN[thread][tn],1.0,gpuC,mytilesizesM[thread][tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,mytilesizesN[thread][tn]*mytilesizesM[thread][tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<mytilesizesN[thread][tn]; j++){
          C_DAXPY(mytilesizesM[thread][tm],1.0,tmp[thread]+j*mytilesizesM[thread][tm],1,C+(j+tn*mytilesizeN[thread])*ldc+tm*mytilesizeM[thread],1);
      }
  }
  free(mytilesizesM[thread]);
  free(mytilesizesN[thread]);
  free(mytilesizesK[thread]);
}
void GPUHelper::GPU_DGEMM_2DTile_tn_threaded_WithCpuStealing(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  throw PsiException("GPU_DGEMM_2DTile_tn_threaded_WithCpuStealing: not implemented",__FILE__,__LINE__);
//DPG commented out to remove statement unreachable warning
/*
  //Tiling((gpumemory-extraroom)/8L,max_mapped_memory/num_gpus/8L,m,n,k);
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
     cudaSetDevice(thread);
     gpuA = gpubuffer[thread];
     gpuB = gpubuffer[thread]+tilesizeM*tilesizeK;
     gpuC = gpubuffer[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
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
                    C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tn=0; tn<ntilesN; tn++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesN[tn]; i++){
                        C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             else{
                for (long int i=0; i<tilesizesM[tm]; i++){
                    C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                }

                for (long int tn=0; tn<ntilesNprime; tn++){

                    if ((tm*ntilesNprime+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesNprime[tn]; i++){
                        C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeNprime+NprimeOffSet)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(j+tn*tilesizeNprime+NprimeOffSet)*ldc+tm*tilesizeM,1);
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
                    C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);

                for (long int tm=0; tm<ntilesM; tm++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesM[tm]; i++){
                        C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             else{
                for (long int i=0; i<tilesizesN[tn]; i++){
                    C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,gpuB+i*tilesizesK[tk],1);
                }

                for (long int tm=0; tm<ntilesMprime; tm++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesMprime[tm]; i++){
                        C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeMprime+MprimeOffSet)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
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
                        C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    for (long int i=0; i<tilesizesN[tn]; i++){
                        C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
                }
                cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                for (long int j=0; j<tilesizesN[tn]; j++){
                    C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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
*/
}
void GPUHelper::GPU_DGEMM_2DTile_tn_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      cudaSetDevice(thread);

      // pointers to gpu memory
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+tilesizeM*tilesizeK*2;
      double*gpuC = gpubuffer[thread]+tilesizeM*tilesizeK*2+tilesizeN*tilesizeK*2;

      long int offsetA = tilesizeM * tilesizeK;
      long int offsetB = tilesizeN * tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));
    
      omp_set_nested(1);
      omp_set_dynamic(0);
if (interleaved_dgemm) {
      // create streams:
      cudaStream_t stream1;
      cudaStreamCreate(&stream1);
      cudaEvent_t estart1,estop1;
      cudaEventCreate(&estart1);
      cudaEventCreate(&estop1);
      cublasSetKernelStream(stream1);

      cudaStream_t stream2;
      cudaStreamCreate(&stream2);
      cudaEvent_t estart2,estop2;
      cudaEventCreate(&estart2);
      cudaEventCreate(&estop2);

      double start = omp_get_wtime();

      // need to transfer data for first tile
      for (long int i=0; i<tilesizesM[tm]; i++){
          C_DCOPY(tilesizesK[0],A+(i+tm*tilesizeM)*lda+0*tilesizeK,1,tmp[thread]+i*tilesizesK[0],1);
      }
      cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice);
      cudaStreamSynchronize(stream1);
      for (long int i=0; i<tilesizesN[tn]; i++){
          C_DCOPY(tilesizesK[0],B+(i+tn*tilesizeN)*ldb+0*tilesizeK,1,tmp[thread]+i*tilesizesK[0],1);
      }
      cudaMemcpyAsync(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice,stream1);
      cudaStreamSynchronize(stream1);

      for (long int tk=0; tk<ntilesK; tk++){

          #pragma omp parallel num_threads(2)
          {

              long int thread2 = omp_get_thread_num();
              if (thread2 == 0) {

                  double * A_curr = ( tk % 2 == 0 ) ? gpuA : gpuA + offsetA;
                  double * B_curr = ( tk % 2 == 0 ) ? gpuB : gpuB + offsetB;

                  cudaEventRecord(estart1,stream1);
                      cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,A_curr,tilesizesK[tk],B_curr,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
                      cudaStreamSynchronize(stream1);
                  cudaEventRecord(estop1,stream1);

              } else {
                  // only copy next tiles if we need them:
                  if ( tk < ntilesK - 1 ) {
                      double * A_next = ( tk % 2 == 0 ) ? gpuA + offsetA : gpuA;
                      double * B_next = ( tk % 2 == 0 ) ? gpuB + offsetB : gpuB;
                      cudaEventRecord(estart2,stream2);
                          for (long int i=0; i<tilesizesM[tm]; i++){
                              C_DCOPY(tilesizesK[tk+1],A+(i+tm*tilesizeM)*lda+(tk+1)*tilesizeK,1,tmp[thread]+i*tilesizesK[tk+1],1);
                          }
                          cudaMemcpyAsync(A_next,tmp[thread],tilesizesM[tm]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                          for (long int i=0; i<tilesizesN[tn]; i++){
                              C_DCOPY(tilesizesK[tk+1],B+(i+tn*tilesizeN)*ldb+(tk+1)*tilesizeK,1,tmp[thread]+i*tilesizesK[tk+1],1);
                          }
                          cudaMemcpyAsync(B_next,tmp[thread],tilesizesN[tn]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                      cudaEventRecord(estop2,stream2);
                      //while( cudaEventQuery(estop) == cudaErrorNotReady );


                  }
              }
              cudaThreadSynchronize();
// TODO: something is wrong with this one ... how to fix ... 
          }
      }
      cublasSetKernelStream(NULL);
      cudaEventDestroy(estart2);
      cudaEventDestroy(estart1);
      cudaEventDestroy(estop1);
      cudaEventDestroy(estop2);
      cudaStreamDestroy(stream1);
      cudaStreamDestroy(stream2);
}else {
      // original version:
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesN[tn]; i++){
              C_DCOPY(tilesizesK[tk],B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesK[tk],1.0,gpuC,tilesizesM[tm]);
      }
}
      omp_set_nested(0);
      omp_set_dynamic(1);
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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

  cudaSetDevice(thread);

  for (long int mn=0; mn<myntilesM[thread]*myntilesN[thread]; mn++){

      // pointers to gpu memory
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread];
      double*gpuC = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread]+mytilesizeN[thread]*mytilesizeK[thread];

      long int tn = mn%myntilesN[thread];
      long int tm = (mn-tn)/myntilesN[thread];

      cudaMemset((void*)gpuC,'\0',mytilesizesM[thread][tm]*mytilesizesN[thread][tn]*sizeof(double));
      for (long int tk=0; tk<myntilesK[thread]; tk++){
          for (long int i=0; i<mytilesizesM[thread][tm]; i++){
              C_DCOPY(mytilesizesK[thread][tk],A+(i+tm*mytilesizeM[thread])*lda+tk*mytilesizeK[thread],1,tmp[thread]+i*mytilesizesK[thread][tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],mytilesizesM[thread][tm]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<mytilesizesN[thread][tn]; i++){
              C_DCOPY(mytilesizesK[thread][tk],B+(i+tn*mytilesizeN[thread])*ldb+tk*mytilesizeK[thread],1,tmp[thread]+i*mytilesizesK[thread][tk],1);
          }
          cudaMemcpy(gpuB,tmp[thread],mytilesizesN[thread][tn]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,mytilesizesM[thread][tm],mytilesizesN[thread][tn],mytilesizesK[thread][tk],alpha,gpuA,mytilesizesK[thread][tk],gpuB,mytilesizesK[thread][tk],1.0,gpuC,mytilesizesM[thread][tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,mytilesizesN[thread][tn]*mytilesizesM[thread][tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<mytilesizesN[thread][tn]; j++){
          C_DAXPY(mytilesizesM[thread][tm],1.0,tmp[thread]+j*mytilesizesM[thread][tm],1,C+(j+tn*mytilesizeN[thread])*ldc+tm*mytilesizeM[thread],1);
      }
  }
  free(mytilesizesM[thread]);
  free(mytilesizesN[thread]);
  free(mytilesizesK[thread]);
}
void GPUHelper::GPU_DGEMM_2DTile_tt_threaded_WithCpuStealing(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  throw PsiException("GPU_DGEMM_2DTile_tt_threaded_WithCpuStealing: not implemented",__FILE__,__LINE__);
//DPG commented out to remove statement unreachable warning
/*
  //Tiling((gpumemory-extraroom)/8L,max_mapped_memory/num_gpus/8L,m,n,k);
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
     cudaSetDevice(thread);
     gpuA = gpubuffer[thread];
     gpuB = gpubuffer[thread]+tilesizeM*tilesizeK;
     gpuC = gpubuffer[thread]+tilesizeM*tilesizeK+tilesizeN*tilesizeK;
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
                    C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                }
                cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                for (long int tn=0; tn<ntilesN; tn++){
                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             else{
                for (long int i=0; i<tilesizesM[tm]; i++){
                    C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                }
                for (long int tn=0; tn<ntilesNprime; tn++){
                    if ((tm*ntilesNprime+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesNprime[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeNprime+NprimeOffSet,1,gpuB+i*tilesizesN[tn],1);
                    }
                    F_DGEMM(transa,transb,tilesizesM[tm],tilesizesNprime[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesNprime[tn],0.0,gpuC,tilesizesM[tm]);
                    for (long int j=0; j<tilesizesNprime[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,gpuC+j*tilesizesM[tm],1,C+(j+tn*tilesizeNprime+NprimeOffSet)*ldc+tm*tilesizeM,1);
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
                    C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                }
                cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                for (long int tm=0; tm<ntilesM; tm++){

                    if ((tm*ntilesN+tn)%num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesM[tm]; i++){
                        C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesM[tm]);
                    cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
                    }
                }
             }
             else{
                for (long int i=0; i<tilesizesK[tk]; i++){
                    C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,gpuB+i*tilesizesN[tn],1);
                }
                for (long int tm=0; tm<ntilesMprime; tm++){

                    if ((tm*ntilesN+tn)%num_cpus+num_gpus!=thread) continue;

                    for (long int i=0; i<tilesizesMprime[tm]; i++){
                        C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeMprime+MprimeOffSet)*lda+tk*tilesizeK,1,gpuA+i*tilesizesK[tk],1);
                    }
                    F_DGEMM(transa,transb,tilesizesMprime[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],0.0,gpuC,tilesizesMprime[tm]);
                    for (long int j=0; j<tilesizesN[tn]; j++){
                        C_DAXPY(tilesizesMprime[tm],1.0,gpuC+j*tilesizesMprime[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeMprime+MprimeOffSet,1);
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
                        C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
                    }
                    cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    for (long int i=0; i<tilesizesK[tk]; i++){
                        C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                    }
                    cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
                    cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
                }
                cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
                for (long int j=0; j<tilesizesN[tn]; j++){
                    C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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
*/
}
// TODO: not thoroughly tested yet.
void GPUHelper::GPU_DGEMM_2DTile_tt_threaded(char transa,char transb,long int m,long int n,long int k,double alpha,double*A,long int lda,double*B,long int ldb,double beta,double*C,long int ldc){

  Tiling((gpumemory-extraroom)/8L,max_mapped_memory_per_thread/8L,m,n,k);

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (long int i=0; i<n*ldc; i++) C[i] *= beta;

  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (long int mn=0; mn<ntilesM*ntilesN; mn++){
      long int thread = 0;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      cudaSetDevice(thread);

      // pointers to gpu memory
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+tilesizeM*tilesizeK*2;
      double*gpuC = gpubuffer[thread]+tilesizeM*tilesizeK*2+tilesizeN*tilesizeK*2;

      long int offsetA = tilesizeM * tilesizeK;
      long int offsetB = tilesizeN * tilesizeK;

      long int tn = mn%ntilesN;
      long int tm = (mn-tn)/ntilesN;

      cudaMemset((void*)gpuC,'\0',tilesizesM[tm]*tilesizesN[tn]*sizeof(double));

      // create streams:
      omp_set_nested(1);
      omp_set_dynamic(0);
if (interleaved_dgemm) {
      // create streams:
      cudaStream_t stream1;
      cudaStreamCreate(&stream1);
      cudaEvent_t estart1,estop1;
      cudaEventCreate(&estart1);
      cudaEventCreate(&estop1);
      cublasSetKernelStream(stream1);

      cudaStream_t stream2;
      cudaStreamCreate(&stream2);
      cudaEvent_t estart2,estop2;
      cudaEventCreate(&estart2);
      cudaEventCreate(&estop2);

      double start = omp_get_wtime();

      // need to transfer data for first tile
      for (long int i=0; i<tilesizesM[tm]; i++){
          C_DCOPY(tilesizesK[0],A+(i+tm*tilesizeM)*lda+0*tilesizeK,1,tmp[thread]+i*tilesizesK[0],1);
      }
      cudaMemcpyAsync(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice,stream1);
      cudaStreamSynchronize(stream1);
      for (long int i=0; i<tilesizesK[0]; i++){
          C_DCOPY(tilesizesN[tn],B+(i+0*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
      }
      cudaMemcpyAsync(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[0]*sizeof(double),cudaMemcpyHostToDevice,stream1);
      cudaStreamSynchronize(stream1);

      for (long int tk=0; tk<ntilesK; tk++){

          #pragma omp parallel num_threads(2)
          {

              long int thread2 = omp_get_thread_num();
              if (thread2 == 0) {

                  double * A_curr = ( tk % 2 == 0 ) ? gpuA : gpuA + offsetA;
                  double * B_curr = ( tk % 2 == 0 ) ? gpuB : gpuB + offsetB;

                  cudaEventRecord(estart1,stream1);
                      cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,A_curr,tilesizesK[tk],B_curr,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
                      cudaStreamSynchronize(stream1);
                  cudaEventRecord(estop1,stream1);

              } else {
                  // only copy next tiles if we need them:
                  if ( tk < ntilesK - 1 ) {
                      double * A_next = ( tk % 2 == 0 ) ? gpuA + offsetA : gpuA;
                      double * B_next = ( tk % 2 == 0 ) ? gpuB + offsetB : gpuB;
                      cudaEventRecord(estart2,stream2);
                          for (long int i=0; i<tilesizesM[tm]; i++){
                              C_DCOPY(tilesizesK[tk+1],A+(i+tm*tilesizeM)*lda+(tk+1)*tilesizeK,1,tmp[thread]+i*tilesizesK[tk+1],1);
                          }
                          cudaMemcpyAsync(A_next,tmp[thread],tilesizesM[tm]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                          for (long int i=0; i<tilesizesK[tk+1]; i++){
                              C_DCOPY(tilesizesN[tn],B+(i+(tk+1)*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
                          }
                          cudaMemcpyAsync(B_next,tmp[thread],tilesizesN[tn]*tilesizesK[tk+1]*sizeof(double),cudaMemcpyHostToDevice,stream2);
                          cudaStreamSynchronize(stream2);
                      cudaEventRecord(estop2,stream2);
                  }
              }
          }
          cudaThreadSynchronize();
      }
      cublasSetKernelStream(NULL);
      cudaEventDestroy(estart2);
      cudaEventDestroy(estart1);
      cudaEventDestroy(estop1);
      cudaEventDestroy(estop2);
      cudaStreamDestroy(stream1);
      cudaStreamDestroy(stream2);
}else {
      // original version:
      for (long int tk=0; tk<ntilesK; tk++){
          for (long int i=0; i<tilesizesM[tm]; i++){
              C_DCOPY(tilesizesK[tk],A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizesK[tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizesM[tm]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<tilesizesK[tk]; i++){
              C_DCOPY(tilesizesN[tn],B+(i+tk*tilesizeK)*ldb+tn*tilesizeN,1,tmp[thread]+i*tilesizesN[tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizesN[tn]*tilesizesK[tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizesM[tm],tilesizesN[tn],tilesizesK[tk],alpha,gpuA,tilesizesK[tk],gpuB,tilesizesN[tn],1.0,gpuC,tilesizesM[tm]);
      }
}
      omp_set_nested(0);
      omp_set_dynamic(1);
      cudaMemcpy(tmp[thread],gpuC,tilesizesN[tn]*tilesizesM[tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<tilesizesN[tn]; j++){
          C_DAXPY(tilesizesM[tm],1.0,tmp[thread]+j*tilesizesM[tm],1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
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

  cudaSetDevice(thread);

  for (long int mn=0; mn<myntilesM[thread]*myntilesN[thread]; mn++){

      // pointers to gpu memory
      double*gpuA = gpubuffer[thread];
      double*gpuB = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread];
      double*gpuC = gpubuffer[thread]+mytilesizeM[thread]*mytilesizeK[thread]+mytilesizeN[thread]*mytilesizeK[thread];

      long int tn = mn%myntilesN[thread];
      long int tm = (mn-tn)/myntilesN[thread];

      cudaMemset((void*)gpuC,'\0',mytilesizesM[thread][tm]*mytilesizesN[thread][tn]*sizeof(double));
      for (long int tk=0; tk<myntilesK[thread]; tk++){
          for (long int i=0; i<mytilesizesM[thread][tm]; i++){
              C_DCOPY(mytilesizesK[thread][tk],A+(i+tm*mytilesizeM[thread])*lda+tk*mytilesizeK[thread],1,tmp[thread]+i*mytilesizesK[thread][tk],1);
          }
          cudaMemcpy(gpuA,tmp[thread],mytilesizesM[thread][tm]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          for (long int i=0; i<mytilesizesK[thread][tk]; i++){
              C_DCOPY(mytilesizesN[thread][tn],B+(i+tk*mytilesizeK[thread])*ldb+tn*mytilesizeN[thread],1,tmp[thread]+i*mytilesizesN[thread][tn],1);
          }
          cudaMemcpy(gpuB,tmp[thread],mytilesizesN[thread][tn]*mytilesizesK[thread][tk]*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,mytilesizesM[thread][tm],mytilesizesN[thread][tn],mytilesizesK[thread][tk],alpha,gpuA,mytilesizesK[thread][tk],gpuB,mytilesizesN[thread][tn],1.0,gpuC,mytilesizesM[thread][tm]);
      }
      cudaMemcpy(tmp[thread],gpuC,mytilesizesN[thread][tn]*mytilesizesM[thread][tm]*sizeof(double),cudaMemcpyDeviceToHost);
      for (long int j=0; j<mytilesizesN[thread][tn]; j++){
          C_DAXPY(mytilesizesM[thread][tm],1.0,tmp[thread]+j*mytilesizesM[thread][tm],1,C+(j+tn*mytilesizeN[thread])*ldc+tm*mytilesizeM[thread],1);
      }
  }
  free(mytilesizesM[thread]);
  free(mytilesizesN[thread]);
  free(mytilesizesK[thread]);
}

void GPUHelper::TilingNoThread(long int mem1,long int mem2,long int m,long int n,long int k){

  long int thread = 0;
  #ifdef _OPENMP
    thread = omp_get_thread_num();
  #endif

  // first tile according to how much space is on gpu
  mytilesizeN[thread] = n;
  mytilesizeM[thread] = m;
  mytilesizeK[thread] = k;
  myntilesM[thread]=myntilesN[thread]=myntilesK[thread]=1L;
  while(mytilesizeN[thread]*mytilesizeM[thread]+mytilesizeK[thread]*(mytilesizeN[thread]+mytilesizeM[thread])>mem1){
     if (mytilesizeN[thread]>mytilesizeM[thread]){
        if (mytilesizeN[thread]>mytilesizeK[thread]){
           myntilesN[thread]++;
           mytilesizeN[thread] = n/myntilesN[thread];
           if (n/myntilesN[thread]<(double)n/myntilesN[thread]) mytilesizeN[thread]++;
        }
        else{
           myntilesK[thread]++;
           mytilesizeK[thread] = k/myntilesK[thread];
           if (k/myntilesK[thread]<(double)k/myntilesK[thread]) mytilesizeK[thread]++;
        }
     }
     else{
        if (mytilesizeM[thread]>mytilesizeK[thread]){
           myntilesM[thread]++;
           mytilesizeM[thread] = m/myntilesM[thread];
           if (m/myntilesM[thread]<(double)m/myntilesM[thread]) mytilesizeM[thread]++;
        }
        else{
           myntilesK[thread]++;
           mytilesizeK[thread] = k/myntilesK[thread];
           if (k/myntilesK[thread]<(double)k/myntilesK[thread]) mytilesizeK[thread]++;
        }
     }
  }

  // ensure each block of A, B, and C will fit in the temporary CPU buffer
  while(mytilesizeN[thread]*mytilesizeM[thread]>mem2){
     if (mytilesizeN[thread]>mytilesizeM[thread]){
        myntilesN[thread]++;
        mytilesizeN[thread] = n/myntilesN[thread];
        if (n/myntilesN[thread]<(double)n/myntilesN[thread]) mytilesizeN[thread]++;
     }
     else{
        myntilesM[thread]++;
        mytilesizeM[thread] = m/myntilesM[thread];
        if (m/myntilesM[thread]<(double)m/myntilesM[thread]) mytilesizeM[thread]++;
     }
  }

  while(mytilesizeN[thread]*mytilesizeK[thread]>mem2){
     if (mytilesizeN[thread]>mytilesizeK[thread]){
        myntilesN[thread]++;
        mytilesizeN[thread] = n/myntilesN[thread];
        if (n/myntilesN[thread]<(double)n/myntilesN[thread]) mytilesizeN[thread]++;
     }
     else{
        myntilesK[thread]++;
        mytilesizeK[thread] = k/myntilesK[thread];
        if (k/myntilesK[thread]<(double)k/myntilesK[thread]) mytilesizeK[thread]++;
     }
  }
  while(mytilesizeK[thread]*mytilesizeM[thread]>mem2){
     if (mytilesizeK[thread]>mytilesizeM[thread]){
        myntilesK[thread]++;
        mytilesizeK[thread] = k/myntilesK[thread];
        if (k/myntilesK[thread]<(double)k/myntilesK[thread]) mytilesizeK[thread]++;
     }
     else{
        myntilesM[thread]++;
        mytilesizeM[thread] = m/myntilesM[thread];
        if (m/myntilesM[thread]<(double)m/myntilesM[thread]) mytilesizeM[thread]++;
     }
  }

  mylasttileN[thread] = n - (myntilesN[thread]-1L)*mytilesizeN[thread];
  mylasttileM[thread] = m - (myntilesM[thread]-1L)*mytilesizeM[thread];
  mylasttileK[thread] = k - (myntilesK[thread]-1L)*mytilesizeK[thread];

  mytilesizesM[thread] = (long int*)malloc(myntilesM[thread]*sizeof(long int));
  mytilesizesN[thread] = (long int*)malloc(myntilesN[thread]*sizeof(long int));
  mytilesizesK[thread] = (long int*)malloc(myntilesK[thread]*sizeof(long int));
  for (long int i=0; i<myntilesM[thread]-1L; i++) mytilesizesM[thread][i] = mytilesizeM[thread];
  for (long int i=0; i<myntilesN[thread]-1L; i++) mytilesizesN[thread][i] = mytilesizeN[thread];
  for (long int i=0; i<myntilesK[thread]-1L; i++) mytilesizesK[thread][i] = mytilesizeK[thread];
  mytilesizesM[thread][myntilesM[thread]-1L] = mylasttileM[thread];
  mytilesizesN[thread][myntilesN[thread]-1L] = mylasttileN[thread];
  mytilesizesK[thread][myntilesK[thread]-1L] = mylasttileK[thread];


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
           tilesizeN = n/ntilesN;
           if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
        }
        else{
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
           //ntilesN++;
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
        //ntilesN++;
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
        //ntilesM++;
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
  // also, make sure we're load balanced - each GPU has the same work to do
  while(ntilesN*ntilesM<num_gpus && (num_gpus % (ntilesN*ntilesM)) == 0){
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

  // double tiling in K so we can pipeline communication/computation
  //ntilesN *= 2;
  //tilesizeN = n/ntilesN;
  //if (n/ntilesN<(double)n/ntilesN) tilesizeN++;
  //ntilesM *= 2;
  //tilesizeM = m/ntilesM;
  //if (m/ntilesM<(double)m/ntilesM) tilesizeM++;


//AED - something is wrong with the tiling ... 
  if (ntilesK < 4)  {
      ntilesK = 8;
      tilesizeK = k/ntilesK;
      if (k/ntilesK<(double)k/ntilesK) tilesizeK++;
  }
    else{
      ntilesK *= 2;
      tilesizeK = k/ntilesK;
      if (k/ntilesK<(double)k/ntilesK) tilesizeK++;
  }
//AED

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

void GPUHelper::DGEMM_Timings() {
    long int m,n,k;
    m = n = k = 20000;
    double * A = (double*)malloc(m*k*sizeof(double));
    double * B = (double*)malloc(n*k*sizeof(double));
    double * C = (double*)malloc(m*n*sizeof(double));
    memset((void*)A,'\0',m*k*sizeof(double));
    memset((void*)B,'\0',n*k*sizeof(double));
    memset((void*)C,'\0',m*n*sizeof(double));

    printf("begin tn:\n");
    m = n = k = 10;
    double start = omp_get_wtime();
    GPUTiledDGEMM('t','n',m,n,k,1.0,A,m,B,n,0.0,C,m);
    cudaThreadSynchronize();
    double end = omp_get_wtime();
    printf("%5i %20.12lf\n",10,m*n*k*2.0/(end-start)/1024./1024./1024.);
    fflush(stdout);
    for (long int i = 1; i < 81; i++) {
        m = n = k = 250 * i;
        double start = omp_get_wtime();
        GPUTiledDGEMM('t','n',m,n,k,1.0,A,m,B,n,0.0,C,m);
        cudaThreadSynchronize();
        double end = omp_get_wtime();
        printf("%5li %20.12lf\n",250*i,m*n*k*2.0/(end-start)/1024./1024./1024.);
        fflush(stdout);
    }
    printf("begin nt:\n");
    m = n = k = 10;
    start = omp_get_wtime();
    GPUTiledDGEMM('n','t',m,n,k,1.0,A,m,B,n,0.0,C,m);
    cudaThreadSynchronize();
    end = omp_get_wtime();
    printf("%5i %20.12lf\n",10,m*n*k*2.0/(end-start)/1024./1024./1024.);
    fflush(stdout);
    for (long int i = 1; i < 81; i++) {
        m = n = k = 250 * i;
        double start = omp_get_wtime();
        GPUTiledDGEMM('n','t',m,n,k,1.0,A,m,B,n,0.0,C,m);
        cudaThreadSynchronize();
        double end = omp_get_wtime();
        printf("%5li %20.12lf\n",250*i,m*n*k*2.0/(end-start)/1024./1024./1024.);
        fflush(stdout);
    }
    printf("begin tt:\n");
    m = n = k = 10;
    start = omp_get_wtime();
    GPUTiledDGEMM('t','t',m,n,k,1.0,A,m,B,n,0.0,C,m);
    cudaThreadSynchronize();
    end = omp_get_wtime();
    printf("%5i %20.12lf\n",10,m*n*k*2.0/(end-start)/1024./1024./1024.);
    fflush(stdout);
    for (long int i = 1; i < 81; i++) {
        m = n = k = 250 * i;
        double start = omp_get_wtime();
        GPUTiledDGEMM('t','t',m,n,k,1.0,A,m,B,n,0.0,C,m);
        cudaThreadSynchronize();
        double end = omp_get_wtime();
        printf("%5li %20.12lf\n",250*i,m*n*k*2.0/(end-start)/1024./1024./1024.);
        fflush(stdout);
    }
    printf("begin nn:\n");
    m = n = k = 10;
    start = omp_get_wtime();
    GPUTiledDGEMM('n','n',m,n,k,1.0,A,m,B,n,0.0,C,m);
    cudaThreadSynchronize();
    end = omp_get_wtime();
    printf("%5i %20.12lf\n",10,m*n*k*2.0/(end-start)/1024./1024./1024.);
    fflush(stdout);
    for (long int i = 1; i < 81; i++) {
        m = n = k = 250 * i;
        double start = omp_get_wtime();
        GPUTiledDGEMM('n','n',m,n,k,1.0,A,m,B,n,0.0,C,m);
        cudaThreadSynchronize();
        double end = omp_get_wtime();
        printf("%5li %20.12lf\n",250*i,m*n*k*2.0/(end-start)/1024./1024./1024.);
        fflush(stdout);
    }
    free(A);
    free(B);
    free(C);
    exit(0);
}

}}//end of namespaces


