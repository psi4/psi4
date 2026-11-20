#include "HRR_F_D_SAMECENTER.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_F_D_SAMECENTER(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        T* base_integrals,         //!< work array for the base [LA+LB|0] integrals of the HRR (not used by all kernels)
        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)
        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively
        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)
        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair
        int gb_size)               //!< number of gridpoints in this grid batch 
{
    if (nprim < 1) return;

    // The gridpoints are processed in blocks.  The integrals of a block are contracted in
    // a buffer that stays in the L1 cache and are then copied to the output arrays.
    enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
    alignas(64) T base_integral[21*g_block];
    T* __restrict base_integral_0 = &base_integral[0*g_block];
    T* __restrict base_integral_1 = &base_integral[1*g_block];
    T* __restrict base_integral_2 = &base_integral[2*g_block];
    T* __restrict base_integral_3 = &base_integral[3*g_block];
    T* __restrict base_integral_4 = &base_integral[4*g_block];
    T* __restrict base_integral_5 = &base_integral[5*g_block];
    T* __restrict base_integral_6 = &base_integral[6*g_block];
    T* __restrict base_integral_7 = &base_integral[7*g_block];
    T* __restrict base_integral_8 = &base_integral[8*g_block];
    T* __restrict base_integral_9 = &base_integral[9*g_block];
    T* __restrict base_integral_10 = &base_integral[10*g_block];
    T* __restrict base_integral_11 = &base_integral[11*g_block];
    T* __restrict base_integral_12 = &base_integral[12*g_block];
    T* __restrict base_integral_13 = &base_integral[13*g_block];
    T* __restrict base_integral_14 = &base_integral[14*g_block];
    T* __restrict base_integral_15 = &base_integral[15*g_block];
    T* __restrict base_integral_16 = &base_integral[16*g_block];
    T* __restrict base_integral_17 = &base_integral[17*g_block];
    T* __restrict base_integral_18 = &base_integral[18*g_block];
    T* __restrict base_integral_19 = &base_integral[19*g_block];
    T* __restrict base_integral_20 = &base_integral[20*g_block];

    // fewer gridpoints than one block: work on a zero padded copy of the inputs
    const bool padded = (gb_size < g_block);
    alignas(64) T pad[9*g_block];
    if (padded) {
        for (int k = 0; k < 9*g_block; ++k) pad[k] = 0;
    }

    for (int g_lo = 0; g_lo < gb_size; g_lo += g_block) {

        // all blocks are complete: the last one overlaps its predecessor
        if (g_lo + g_block > gb_size && !padded) g_lo = gb_size - g_block;

        for (int prim = 0; prim < nprim; ++prim) {

            if (padded) {
                for (int k = 0; k < gb_size; ++k) {
                    pad[0*g_block+k] = v_PGx[prim*gb_size+k];
                    pad[1*g_block+k] = v_PGy[prim*gb_size+k];
                    pad[2*g_block+k] = v_PGz[prim*gb_size+k];
                    pad[3*g_block+k] = v_scaled_boys[prim*6*gb_size+0*gb_size+k];
                    pad[4*g_block+k] = v_scaled_boys[prim*6*gb_size+1*gb_size+k];
                    pad[5*g_block+k] = v_scaled_boys[prim*6*gb_size+2*gb_size+k];
                    pad[6*g_block+k] = v_scaled_boys[prim*6*gb_size+3*gb_size+k];
                    pad[7*g_block+k] = v_scaled_boys[prim*6*gb_size+4*gb_size+k];
                    pad[8*g_block+k] = v_scaled_boys[prim*6*gb_size+5*gb_size+k];
                }
            }
            const T* __restrict PGx = padded ? &pad[0*g_block] : &v_PGx[prim*gb_size + g_lo];
            const T* __restrict PGy = padded ? &pad[1*g_block] : &v_PGy[prim*gb_size + g_lo];
            const T* __restrict PGz = padded ? &pad[2*g_block] : &v_PGz[prim*gb_size + g_lo];
            const T* __restrict scaled_boys_0 = padded ? &pad[3*g_block] : &v_scaled_boys[prim*6*gb_size+0*gb_size + g_lo];
            const T* __restrict scaled_boys_1 = padded ? &pad[4*g_block] : &v_scaled_boys[prim*6*gb_size+1*gb_size + g_lo];
            const T* __restrict scaled_boys_2 = padded ? &pad[5*g_block] : &v_scaled_boys[prim*6*gb_size+2*gb_size + g_lo];
            const T* __restrict scaled_boys_3 = padded ? &pad[6*g_block] : &v_scaled_boys[prim*6*gb_size+3*gb_size + g_lo];
            const T* __restrict scaled_boys_4 = padded ? &pad[7*g_block] : &v_scaled_boys[prim*6*gb_size+4*gb_size + g_lo];
            const T* __restrict scaled_boys_5 = padded ? &pad[8*g_block] : &v_scaled_boys[prim*6*gb_size+5*gb_size + g_lo];

            const T tpi = 0.5 / v_p[prim];
            const T tpi_2 = 2*tpi;
            const T tpi_3 = 3*tpi;
            const T tpi_4 = 4*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x1 = -PGy[g]*scaled_boys_5[g];
                    const T x2 = tpi*x0 - PGy[g]*x1;
                    const T x3 = -PGy[g]*x0;
                    const T x4 = tpi_2*x3 - PGy[g]*x2;
                    const T x5 = -PGx[g]*x4;
                    base_integral_11[g] = -PGz[g]*x5;
                    const T x6 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x7 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x8 = x6 - x7;
                    const T x9 = tpi*x8 - PGy[g]*x3;
                    const T x10 = tpi_3*x9 - PGy[g]*x4;
                    base_integral_10[g] = -PGx[g]*x10;
                    const T x11 = -tpi_2*x8 - x9;
                    const T x12 = PGy[g]*x11;
                    base_integral_6[g] = tpi*x12 - PGx[g]*x5;
                    base_integral_15[g] = tpi_4*x12 - PGy[g]*x10;
                    base_integral_16[g] = -PGz[g]*x10;
                    const T x13 = -PGx[g]*x0;
                    const T x14 = tpi*x8 - PGx[g]*x13;
                    const T x15 = -PGx[g]*scaled_boys_5[g];
                    const T x16 = tpi*x0 - PGx[g]*x15;
                    const T x17 = -PGy[g]*x16;
                    const T x18 = tpi*x14 - PGy[g]*x17;
                    base_integral_7[g] = -PGz[g]*x18;
                    const T x19 = tpi_2*x13 - PGx[g]*x16;
                    const T x20 = tpi_3*x14 - PGx[g]*x19;
                    base_integral_1[g] = -PGy[g]*x20;
                    base_integral_2[g] = -PGz[g]*x20;
                    const T x21 = -tpi_2*x8 - x14;
                    const T x22 = PGx[g]*x21;
                    base_integral_0[g] = tpi_4*x22 - PGx[g]*x20;
                    const T x23 = -PGy[g]*x19;
                    base_integral_3[g] = tpi*x22 - PGy[g]*x23;
                    const T x24 = -PGz[g]*x19;
                    base_integral_5[g] = tpi*x22 - PGz[g]*x24;
                    base_integral_4[g] = -PGy[g]*x24;
                    const T x25 = -PGz[g]*scaled_boys_5[g];
                    const T x26 = tpi*x0 - PGz[g]*x25;
                    const T x27 = -PGz[g]*x0;
                    const T x28 = tpi_2*x27 - PGz[g]*x26;
                    const T x29 = -PGx[g]*x28;
                    const T x30 = tpi*x8 - PGz[g]*x27;
                    const T x31 = -tpi_2*x8 - x30;
                    const T x32 = PGz[g]*x31;
                    base_integral_9[g] = tpi*x32 - PGx[g]*x29;
                    const T x33 = -PGy[g]*x28;
                    base_integral_13[g] = -PGx[g]*x33;
                    base_integral_18[g] = tpi*x32 - PGy[g]*x33;
                    const T x34 = tpi_3*x30 - PGz[g]*x28;
                    base_integral_20[g] = tpi_4*x32 - PGz[g]*x34;
                    base_integral_14[g] = -PGx[g]*x34;
                    base_integral_19[g] = -PGy[g]*x34;
                    const T x35 = -PGy[g]*x26;
                    const T x36 = tpi*x30 - PGy[g]*x35;
                    base_integral_12[g] = -PGx[g]*x36;
                    const T x37 = -PGy[g]*x30;
                    base_integral_17[g] = tpi_2*x37 - PGy[g]*x36;
                    const T x38 = -PGx[g]*x35;
                    base_integral_8[g] = tpi*x37 - PGx[g]*x38;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x1 = -PGy[g]*scaled_boys_5[g];
                    const T x2 = tpi*x0 - PGy[g]*x1;
                    const T x3 = -PGy[g]*x0;
                    const T x4 = tpi_2*x3 - PGy[g]*x2;
                    const T x5 = -PGx[g]*x4;
                    base_integral_11[g] += -PGz[g]*x5;
                    const T x6 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x7 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x8 = x6 - x7;
                    const T x9 = tpi*x8 - PGy[g]*x3;
                    const T x10 = tpi_3*x9 - PGy[g]*x4;
                    base_integral_10[g] += -PGx[g]*x10;
                    const T x11 = -tpi_2*x8 - x9;
                    const T x12 = PGy[g]*x11;
                    base_integral_6[g] += tpi*x12 - PGx[g]*x5;
                    base_integral_15[g] += tpi_4*x12 - PGy[g]*x10;
                    base_integral_16[g] += -PGz[g]*x10;
                    const T x13 = -PGx[g]*x0;
                    const T x14 = tpi*x8 - PGx[g]*x13;
                    const T x15 = -PGx[g]*scaled_boys_5[g];
                    const T x16 = tpi*x0 - PGx[g]*x15;
                    const T x17 = -PGy[g]*x16;
                    const T x18 = tpi*x14 - PGy[g]*x17;
                    base_integral_7[g] += -PGz[g]*x18;
                    const T x19 = tpi_2*x13 - PGx[g]*x16;
                    const T x20 = tpi_3*x14 - PGx[g]*x19;
                    base_integral_1[g] += -PGy[g]*x20;
                    base_integral_2[g] += -PGz[g]*x20;
                    const T x21 = -tpi_2*x8 - x14;
                    const T x22 = PGx[g]*x21;
                    base_integral_0[g] += tpi_4*x22 - PGx[g]*x20;
                    const T x23 = -PGy[g]*x19;
                    base_integral_3[g] += tpi*x22 - PGy[g]*x23;
                    const T x24 = -PGz[g]*x19;
                    base_integral_5[g] += tpi*x22 - PGz[g]*x24;
                    base_integral_4[g] += -PGy[g]*x24;
                    const T x25 = -PGz[g]*scaled_boys_5[g];
                    const T x26 = tpi*x0 - PGz[g]*x25;
                    const T x27 = -PGz[g]*x0;
                    const T x28 = tpi_2*x27 - PGz[g]*x26;
                    const T x29 = -PGx[g]*x28;
                    const T x30 = tpi*x8 - PGz[g]*x27;
                    const T x31 = -tpi_2*x8 - x30;
                    const T x32 = PGz[g]*x31;
                    base_integral_9[g] += tpi*x32 - PGx[g]*x29;
                    const T x33 = -PGy[g]*x28;
                    base_integral_13[g] += -PGx[g]*x33;
                    base_integral_18[g] += tpi*x32 - PGy[g]*x33;
                    const T x34 = tpi_3*x30 - PGz[g]*x28;
                    base_integral_20[g] += tpi_4*x32 - PGz[g]*x34;
                    base_integral_14[g] += -PGx[g]*x34;
                    base_integral_19[g] += -PGy[g]*x34;
                    const T x35 = -PGy[g]*x26;
                    const T x36 = tpi*x30 - PGy[g]*x35;
                    base_integral_12[g] += -PGx[g]*x36;
                    const T x37 = -PGy[g]*x30;
                    base_integral_17[g] += tpi_2*x37 - PGy[g]*x36;
                    const T x38 = -PGx[g]*x35;
                    base_integral_8[g] += tpi*x37 - PGx[g]*x38;

                } // for g
            } // if/else prim
        } // for prim

        // copy the block to the output arrays
        static const short row[60] = {
            0, 1, 2, 3, 4, 5, 1, 3, 4, 6, 7, 8, 2, 4, 5, 7, 8, 9, 3, 6,
            7, 10, 11, 12, 4, 7, 8, 11, 12, 13, 5, 8, 9, 12, 13, 14, 6, 10, 11, 15,
            16, 17, 7, 11, 12, 16, 17, 18, 8, 12, 13, 17, 18, 19, 9, 13, 14, 18, 19, 20,
        };
        for (int i = 0; i < 60; ++i) {
            const T* __restrict src = &base_integral[row[i]*g_block];
            T* __restrict dst = &integrals[i*gb_size + g_lo];
            if (padded) {
                SNK_COPY(dst, src, gb_size);
            } else {
                SNK_COPY(dst, src, g_block);
            }
        }
    } // for g_lo
} // HRR_F_D_SAMECENTER

template void HRR_F_D_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_F_D_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

