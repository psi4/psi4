#ifndef GPU_KERNELS_H
#define GPU_KERNELS_H

#ifdef _OPENMP
  #include<omp.h>
#endif

// cuda libraries
#include<cuda.h>
#include<cublas.h>
#include<cuda_runtime.h>

#define NUMTHREADS 8
#define MAXBLOCKS 65535

namespace psi{

// gpu functions that nvcc only needs to see
__global__ void GPUt2Plust1(int ns,int h,double*t2,double*t1);
__global__ void GPUt2Plust1_and_E2klcd3(int ns,int h,double*t2,double*t1,double*E2klcd1,double*E2klcd3);
__global__ void AddPermutedOnGPU(int ns,int h,double*in,double*out);
__global__ void GPUFill_I2iajk_add(int ns,int h,double*in,double*out);
__global__ void GPUc2Sym1(int ns,int h,double*in,double*out);
__global__ void GPUc2Sym2(int ns,int h,double*in,double*out);
__global__ void GPUSymmAdd1(int ns,int h,double*in,double*out);
__global__ void GPUSymmAdd2(int ns,int h,double*in,double*out);
__global__ void gpuCopyArray(double*out,double*in,int dim);
__device__ int GPUPosition(int i,int j);

__global__ void GPUFill_I2iajk_and_c2Sym1(int ns,int h,double*gpuv,
                   double*gpuw,double*gput2,double*gputempw);
__global__ void GPUc2Sym2_onefunction(int ns,int h,double*gput2,
                   double*gput1,double*gputempw);

__global__ void GPUt2Plus2t1(int ns,int h,double*out,double*t2,double*t1,double*newt2);
__global__ void GPUv2MinusHalfv2(int ns,int h,double*out,double*v2);
__global__ void GPUPermute_iabj_to_aijb(int ns,int h,double*in,double*out);
__global__ void GPUFill_I2iabj(int ns,int h,double*in1,double*in2,double*out);
__global__ void GPU2t2Minust2(int ns,int h,double*out,double*t2);
__global__ void GPUFill_t2_I2iajb(int ns,int h,double*in,double*out);

__global__ void GPUt2Plus2t1_2(int ns,int h,double*out,double*t2,double*t1);
__global__ void GPUt2Plus2t1_and_E2klcd2(int ns,int h,double*out,double*t2,double*t1,double*E2klcd1,double*E2klcd2);
__global__ void GPUFillI2iajb1(int ns,int h,double*in,double*out);
__global__ void GPUFillI2iajb2(int ns,int h,double*in,double*out);
__global__ void GPUFill_t2_I2iabj1(int ns,int h,double*in,double*out);
__global__ void GPUFill_t2_I2iabj2(int ns,int h,double*in,double*out);
__global__ void GPUPermute_tikbc(double*in,double*out,int ns,int h);
};

#endif
