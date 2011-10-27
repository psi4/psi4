#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cuda.h>
#include<cublas.h>
#include<cuda_runtime.h>
#include"gpudfjkhelper.h"
#include"blas.h"
#ifdef _OPENMP
   #include<omp.h>
#endif

using namespace psi;

namespace psi{

/**
 * constructor
 */
GPUDFJKHelper::GPUDFJKHelper(){}

/**
 * destructor
 */
GPUDFJKHelper::~GPUDFJKHelper(){}

/**
 * check for errors following cuda calls
 */
inline void GPUDFJKHelper::Check_CUDA_Error(FILE*fp,const char *message){
  cudaError_t error = cudaGetLastError();
  if (error!=cudaSuccess) {
     fprintf(fp,"\n  ERROR: %s: %s\n\n", message, cudaGetErrorString(error) );
     fflush(fp);
     exit(-1);
  }
}

/**
 * free host memory
 */
void GPUDFJKHelper::Finalize(){
  // free cpu memory and reset each device
  #pragma omp parallel for schedule (static)
  for (int i=0; i<num_gpus; i++){
      int nthreads = omp_get_num_threads();
      int thread = 1;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      if (i<nthreads){
         cudaFreeHost(tmp[thread]);
         Check_CUDA_Error(stdout,"cudaFreeHost");
         cudaThreadExit();
         Check_CUDA_Error(stdout,"cudaThreadExit");
      }
  }
}

/**
 * initialize temporary array for mapping cpu memory before transfer
 * to gpu.  also, set number of gpus. also, initialize cublas
 */
void GPUDFJKHelper::Initialize(int max_rows,int max_nocc,int nbf){

  // allocate cpu memory
  int gpu_id; 
  cudaGetDevice(&gpu_id);
  struct cudaDeviceProp cudaProp;
  cudaGetDeviceProperties( &cudaProp,gpu_id );

  // memory available to the device (less a little extra)
  double memory = cudaProp.totalGlobalMem/8. - 200.*1024*1024/8.;
  int dim = (int)memory;

  cublasInit();
  Check_CUDA_Error(stdout,"cublasInit");
  cudaGetDeviceCount(&num_gpus);
  Check_CUDA_Error(stdout,"cudaGetDeviceCount");
  cudaThreadExit();
  Check_CUDA_Error(stdout,"cudaThreadExit");
  tmp = (double**)malloc(num_gpus*sizeof(double*));

  // initialize each device and device memory:
  #pragma omp parallel for schedule (static)
  for (int i=0; i<num_gpus; i++){
      int nthreads = omp_get_num_threads();
      int thread = 1;
      #ifdef _OPENMP
        thread = omp_get_thread_num();
      #endif
      if (i<nthreads){
         cudaSetDevice(thread);
         Check_CUDA_Error(stdout,"cudaSetDevice");
         cudaMallocHost((void**)&tmp[thread],dim*sizeof(double));
         Check_CUDA_Error(stdout,"cudaMallocHost");
     }
  }

}


/**
 * dgemm assuming no tiling is necessary
 */
void GPUDFJKHelper::GPU_DGEMM(char transa,char transb,int m,int n,int k,double alpha,double*A,int lda,double*B,int ldb,double beta,double*C,int ldc){
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
 * dgemm using a 2-dimensional tile.  this guy is also threaded.
 */
void GPUDFJKHelper::GPU_DGEMM_2DTile(char transa,char transb,int m,int n,int k,double alpha,double*A,int lda,double*B,int ldb,double beta,double*C,int ldc,int thread){

  // cpu threads:
  if (thread>=num_gpus){
     F_DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
     return;
  } 

  int nthreads = omp_get_num_threads();

  struct cudaDeviceProp cudaProp;
  int gpu_id;
  cudaGetDevice(&gpu_id);
  cudaGetDeviceProperties( &cudaProp,gpu_id );

  // memory available to the device (less a little extra)
  double memory = cudaProp.totalGlobalMem/8. - 200.*1024*1024/8.;

  // test not having enough memory:
  //memory = (m*n+k*(m+n))/7.;
     
  // determine what tiling should be:
  int ntilesN,ntilesM,ntilesK,tilesizeK,tilesizeN,tilesizeM,lasttileK,lasttileN,lasttileM;
  tilesizeN = n;
  tilesizeM = m;
  tilesizeK = k;
  ntilesM=ntilesN=ntilesK=1;
  while(tilesizeN*tilesizeM+tilesizeK*(tilesizeN+tilesizeM)>(int)memory){
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
  lasttileN = n - (ntilesN-1)*tilesizeN;
  lasttileM = m - (ntilesM-1)*tilesizeM;
  lasttileK = k - (ntilesK-1)*tilesizeK;

  // allocate gpu memory
  double*gpuA,*gpuB,*gpuC;
  cudaMalloc((void**)&gpuA,tilesizeM*tilesizeK*sizeof(double));
  cudaMalloc((void**)&gpuB,tilesizeN*tilesizeK*sizeof(double));
  cudaMalloc((void**)&gpuC,tilesizeM*tilesizeN*sizeof(double));

  int tm,tn,tk,i,j;

  // initialize result
  if (beta==0.0) 
     memset((void*)C,'\0',n*ldc*sizeof(double));
  else           
     for (i=0; i<n*ldc; i++) C[i] *= beta;

  for (tm=0; tm<ntilesM-1; tm++){
      for (tn=0; tn<ntilesN-1; tn++){
          for (tk=0; tk<ntilesK-1; tk++){

              for (i=0; i<tilesizeM; i++){
                  F_DCOPY(tilesizeK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
              }
              cudaMemcpy(gpuA,tmp[thread],tilesizeM*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
              for (i=0; i<tilesizeN; i++){
                  F_DCOPY(tilesizeK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
              }
              cudaMemcpy(gpuB,tmp[thread],tilesizeN*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
              cublasDgemm(transa,transb,tilesizeM,tilesizeN,tilesizeK,alpha,gpuA,tilesizeK,gpuB,tilesizeK,0.0,gpuC,tilesizeM);
              cudaMemcpy(tmp[thread],gpuC,tilesizeN*tilesizeM*sizeof(double),cudaMemcpyDeviceToHost);
              for (j=0; j<tilesizeN; j++){
                  F_DAXPY(tilesizeM,1.0,tmp[thread]+j*tilesizeM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
              }
          } // end of tiles over k

          // last tile of k
          tk = ntilesK-1;

          for (i=0; i<tilesizeM; i++){
              F_DCOPY(lasttileK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizeM*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
          for (i=0; i<tilesizeN; i++){
              F_DCOPY(lasttileK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizeN*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizeM,tilesizeN,lasttileK,alpha,gpuA,lasttileK,gpuB,lasttileK,0.0,gpuC,tilesizeM);
          cudaMemcpy(tmp[thread],gpuC,tilesizeN*tilesizeM*sizeof(double),cudaMemcpyDeviceToHost);
          for (j=0; j<tilesizeN; j++){
              F_DAXPY(tilesizeM,1.0,tmp[thread]+j*tilesizeM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
          }

      } // end of tiles over n

  } // end of tiles over m

  // last tiles of m and n:
  tm = ntilesM-1;
  for (tn=0; tn<ntilesN-1; tn++){
      for (tk=0; tk<ntilesK-1; tk++){

          for (i=0; i<lasttileM; i++){
              F_DCOPY(tilesizeK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
          }
          cudaMemcpy(gpuA,tmp[thread],lasttileM*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
          for (i=0; i<tilesizeN; i++){
              F_DCOPY(tilesizeK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
          }
          cudaMemcpy(gpuB,tmp[thread],tilesizeN*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,lasttileM,tilesizeN,tilesizeK,alpha,gpuA,tilesizeK,gpuB,tilesizeK,0.0,gpuC,lasttileM);
          cudaMemcpy(tmp[thread],gpuC,tilesizeN*lasttileM*sizeof(double),cudaMemcpyDeviceToHost);
          for (j=0; j<tilesizeN; j++){
              F_DAXPY(lasttileM,1.0,tmp[thread]+j*lasttileM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
          }
      } // end of tiles over k

      // last tile of k
      tk = ntilesK-1;
 
      for (i=0; i<lasttileM; i++){
          F_DCOPY(lasttileK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
      }
      cudaMemcpy(gpuA,tmp[thread],lasttileM*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
      for (i=0; i<tilesizeN; i++){
          F_DCOPY(lasttileK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
      }
      cudaMemcpy(gpuB,tmp[thread],tilesizeN*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
      cublasDgemm(transa,transb,lasttileM,tilesizeN,lasttileK,alpha,gpuA,lasttileK,gpuB,lasttileK,0.0,gpuC,lasttileM);
      cudaMemcpy(tmp[thread],gpuC,tilesizeN*lasttileM*sizeof(double),cudaMemcpyDeviceToHost);
      for (j=0; j<tilesizeN; j++){
          F_DAXPY(lasttileM,1.0,tmp[thread]+j*lasttileM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  } // end of tiles over n


  tn = ntilesN-1; 
  for (tm=0; tm<ntilesM-1; tm++){
      for (tk=0; tk<ntilesK-1; tk++){

          for (i=0; i<tilesizeM; i++){
              F_DCOPY(tilesizeK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
          }
          cudaMemcpy(gpuA,tmp[thread],tilesizeM*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
          for (i=0; i<lasttileN; i++){
              F_DCOPY(tilesizeK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
          }
          cudaMemcpy(gpuB,tmp[thread],lasttileN*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
          cublasDgemm(transa,transb,tilesizeM,lasttileN,tilesizeK,alpha,gpuA,tilesizeK,gpuB,tilesizeK,0.0,gpuC,tilesizeM);
          cudaMemcpy(tmp[thread],gpuC,lasttileN*tilesizeM*sizeof(double),cudaMemcpyDeviceToHost);
          for (j=0; j<lasttileN; j++){
              F_DAXPY(tilesizeM,1.0,tmp[thread]+j*tilesizeM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
          }
      } // end of tiles over k

      // last tile of k
      tk = ntilesK-1;

      for (i=0; i<tilesizeM; i++){
          F_DCOPY(lasttileK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
      }
      cudaMemcpy(gpuA,tmp[thread],tilesizeM*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
      for (i=0; i<lasttileN; i++){
          F_DCOPY(lasttileK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
      }
      cudaMemcpy(gpuB,tmp[thread],lasttileN*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
      cublasDgemm(transa,transb,tilesizeM,lasttileN,lasttileK,alpha,gpuA,lasttileK,gpuB,lasttileK,0.0,gpuC,tilesizeM);
      cudaMemcpy(tmp[thread],gpuC,lasttileN*tilesizeM*sizeof(double),cudaMemcpyDeviceToHost);
      for (j=0; j<lasttileN; j++){
          F_DAXPY(tilesizeM,1.0,tmp[thread]+j*tilesizeM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  } // end of tiles over m

  tn = ntilesN-1; 
  tm = ntilesM-1; 
  for (tk=0; tk<ntilesK-1; tk++){

      for (i=0; i<lasttileM; i++){
          F_DCOPY(tilesizeK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
      }
      cudaMemcpy(gpuA,tmp[thread],lasttileM*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
      for (i=0; i<lasttileN; i++){
          F_DCOPY(tilesizeK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*tilesizeK,1);
      }
      cudaMemcpy(gpuB,tmp[thread],lasttileN*tilesizeK*sizeof(double),cudaMemcpyHostToDevice);
      cublasDgemm(transa,transb,lasttileM,lasttileN,tilesizeK,alpha,gpuA,tilesizeK,gpuB,tilesizeK,0.0,gpuC,lasttileM);
      cudaMemcpy(tmp[thread],gpuC,lasttileN*lasttileM*sizeof(double),cudaMemcpyDeviceToHost);
      for (j=0; j<lasttileN; j++){
          F_DAXPY(lasttileM,1.0,tmp[thread]+j*lasttileM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
      }
  } // end of tiles over k

  // last tile of k
  tk = ntilesK-1;

  for (i=0; i<lasttileM; i++){
      F_DCOPY(lasttileK,A+(i+tm*tilesizeM)*lda+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
  }
  cudaMemcpy(gpuA,tmp[thread],lasttileM*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
  for (i=0; i<lasttileN; i++){
      F_DCOPY(lasttileK,B+(i+tn*tilesizeN)*ldb+tk*tilesizeK,1,tmp[thread]+i*lasttileK,1);
  }
  cudaMemcpy(gpuB,tmp[thread],lasttileN*lasttileK*sizeof(double),cudaMemcpyHostToDevice);
  cublasDgemm(transa,transb,lasttileM,lasttileN,lasttileK,alpha,gpuA,lasttileK,gpuB,lasttileK,0.0,gpuC,lasttileM);
  cudaMemcpy(tmp[thread],gpuC,lasttileN*lasttileM*sizeof(double),cudaMemcpyDeviceToHost);
  for (j=0; j<lasttileN; j++){
      F_DAXPY(lasttileM,1.0,tmp[thread]+j*lasttileM,1,C+(j+tn*tilesizeN)*ldc+tm*tilesizeM,1);
  }
  cudaFree(gpuA);
  cudaFree(gpuB);
  cudaFree(gpuC);

}
}//end of namespace
