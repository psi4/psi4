#ifdef USING_cuEST

#include "mp2_kernels.h"

// this kernel is adapted from the CUDALibrarySamples mp2 code
__global__
void make_pair_amplitudes_ij(
    long long i,
    long long j,
    long long nvir,
    const double* iajb,
    const double* eps_occ,
    const double* eps_vir,
    double* iajb_d
) {
    // for maximum memory coalescing, idx x handles b and idx y handles a
    long long rel_a = (long long) blockIdx.y * blockDim.y + threadIdx.y;
    long long rel_b = (long long) blockIdx.x * blockDim.x + threadIdx.x;

    if (rel_a < nvir && rel_b < nvir) {
        double denominator = eps_occ[i] + eps_occ[j] - eps_vir[rel_a] - eps_vir[rel_b];
        double numerator = iajb[rel_a * nvir + rel_b];
        iajb_d[rel_a * nvir + rel_b] = numerator * (fabs(denominator) < 1e-20 ? 1e20 : 1.0 / denominator);
    }
}

void pair_amplitudes_ij(
    int i,
    int j,
    int nvir,
    const double* iajb,
    const double* eps_occ,
    const double* eps_vir,
    double* iajb_d,
    cudaStream_t stream
) {
    // schedule thread blocks of 256 elements each (16x16)
    dim3 DimGrid((nvir + 16 - 1) / 16, (nvir + 16 - 1) / 16, 1);
    dim3 DimBlock(16, 16, 1);
    make_pair_amplitudes_ij<<<DimGrid, DimBlock, 0, stream>>>(i, j, nvir, iajb, eps_occ, eps_vir, iajb_d);
}

#endif // USING_cuEST