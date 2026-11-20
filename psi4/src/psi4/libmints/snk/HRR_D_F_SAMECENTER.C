#include "HRR_D_F_SAMECENTER.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_D_F_SAMECENTER(
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
            const T tpi_p2 = tpi*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x1 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x2 = x0 - x1;
                    const T x3 = tpi_p2*x2;
                    const T x4 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x5 = tpi*x4;
                    const T x6 = PGy[g]*x5;
                    const T x7 = x3 + PGy[g]*x6;
                    const T x8 = 2*x3 + x7;
                    const T x9 = PGy[g]*x8;
                    const T x10 = PGy[g]*scaled_boys_5[g];
                    const T x11 = x5 + PGy[g]*x10;
                    const T x12 = PGy[g]*x11 + 2*x6;
                    const T x13 = PGx[g]*x12;
                    const T x14 = PGy[g]*x12 + 3*x7;
                    base_integral_6[g] = -PGx[g]*x13 - x9;
                    base_integral_11[g] = -PGz[g]*x13;
                    base_integral_15[g] = -4*x9 - PGy[g]*x14;
                    base_integral_10[g] = -PGx[g]*x14;
                    base_integral_16[g] = -PGz[g]*x14;
                    const T x15 = PGx[g]*x5;
                    const T x16 = PGx[g]*x15 + x3;
                    const T x17 = PGx[g]*scaled_boys_5[g];
                    const T x18 = PGx[g]*x17 + x5;
                    const T x19 = PGy[g]*x18;
                    const T x20 = PGx[g]*x18 + 2*x15;
                    const T x21 = x16 + PGy[g]*x19;
                    base_integral_7[g] = -PGz[g]*x21;
                    const T x22 = PGx[g]*x20 + 3*x16;
                    const T x23 = 2*x3 + x16;
                    const T x24 = PGx[g]*x23;
                    base_integral_0[g] = -PGx[g]*x22 - 4*x24;
                    const T x25 = PGy[g]*x20;
                    base_integral_3[g] = -x24 - PGy[g]*x25;
                    const T x26 = PGz[g]*x20;
                    base_integral_5[g] = -x24 - PGz[g]*x26;
                    base_integral_4[g] = -PGy[g]*x26;
                    base_integral_1[g] = -PGy[g]*x22;
                    base_integral_2[g] = -PGz[g]*x22;
                    const T x27 = PGz[g]*x5;
                    const T x28 = x3 + PGz[g]*x27;
                    const T x29 = 2*x3 + x28;
                    const T x30 = PGz[g]*x29;
                    const T x31 = PGz[g]*scaled_boys_5[g];
                    const T x32 = x5 + PGz[g]*x31;
                    const T x33 = PGz[g]*x32 + 2*x27;
                    const T x34 = PGx[g]*x33;
                    base_integral_9[g] = -PGx[g]*x34 - x30;
                    const T x35 = 3*x28 + PGz[g]*x33;
                    base_integral_20[g] = -4*x30 - PGz[g]*x35;
                    const T x36 = PGy[g]*x33;
                    base_integral_18[g] = -x30 - PGy[g]*x36;
                    base_integral_13[g] = -PGx[g]*x36;
                    base_integral_14[g] = -PGx[g]*x35;
                    base_integral_19[g] = -PGy[g]*x35;
                    const T x37 = PGy[g]*x32;
                    const T x38 = PGy[g]*x37 + x28;
                    base_integral_12[g] = -PGx[g]*x38;
                    const T x39 = PGy[g]*x28;
                    base_integral_17[g] = -2*x39 - PGy[g]*x38;
                    const T x40 = PGx[g]*x37;
                    base_integral_8[g] = -PGx[g]*x40 - x39;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x1 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x2 = x0 - x1;
                    const T x3 = tpi_p2*x2;
                    const T x4 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x5 = tpi*x4;
                    const T x6 = PGy[g]*x5;
                    const T x7 = x3 + PGy[g]*x6;
                    const T x8 = 2*x3 + x7;
                    const T x9 = PGy[g]*x8;
                    const T x10 = PGy[g]*scaled_boys_5[g];
                    const T x11 = x5 + PGy[g]*x10;
                    const T x12 = PGy[g]*x11 + 2*x6;
                    const T x13 = PGx[g]*x12;
                    const T x14 = PGy[g]*x12 + 3*x7;
                    base_integral_6[g] += -PGx[g]*x13 - x9;
                    base_integral_11[g] += -PGz[g]*x13;
                    base_integral_15[g] += -4*x9 - PGy[g]*x14;
                    base_integral_10[g] += -PGx[g]*x14;
                    base_integral_16[g] += -PGz[g]*x14;
                    const T x15 = PGx[g]*x5;
                    const T x16 = PGx[g]*x15 + x3;
                    const T x17 = PGx[g]*scaled_boys_5[g];
                    const T x18 = PGx[g]*x17 + x5;
                    const T x19 = PGy[g]*x18;
                    const T x20 = PGx[g]*x18 + 2*x15;
                    const T x21 = x16 + PGy[g]*x19;
                    base_integral_7[g] += -PGz[g]*x21;
                    const T x22 = PGx[g]*x20 + 3*x16;
                    const T x23 = 2*x3 + x16;
                    const T x24 = PGx[g]*x23;
                    base_integral_0[g] += -PGx[g]*x22 - 4*x24;
                    const T x25 = PGy[g]*x20;
                    base_integral_3[g] += -x24 - PGy[g]*x25;
                    const T x26 = PGz[g]*x20;
                    base_integral_5[g] += -x24 - PGz[g]*x26;
                    base_integral_4[g] += -PGy[g]*x26;
                    base_integral_1[g] += -PGy[g]*x22;
                    base_integral_2[g] += -PGz[g]*x22;
                    const T x27 = PGz[g]*x5;
                    const T x28 = x3 + PGz[g]*x27;
                    const T x29 = 2*x3 + x28;
                    const T x30 = PGz[g]*x29;
                    const T x31 = PGz[g]*scaled_boys_5[g];
                    const T x32 = x5 + PGz[g]*x31;
                    const T x33 = PGz[g]*x32 + 2*x27;
                    const T x34 = PGx[g]*x33;
                    base_integral_9[g] += -PGx[g]*x34 - x30;
                    const T x35 = 3*x28 + PGz[g]*x33;
                    base_integral_20[g] += -4*x30 - PGz[g]*x35;
                    const T x36 = PGy[g]*x33;
                    base_integral_18[g] += -x30 - PGy[g]*x36;
                    base_integral_13[g] += -PGx[g]*x36;
                    base_integral_14[g] += -PGx[g]*x35;
                    base_integral_19[g] += -PGy[g]*x35;
                    const T x37 = PGy[g]*x32;
                    const T x38 = PGy[g]*x37 + x28;
                    base_integral_12[g] += -PGx[g]*x38;
                    const T x39 = PGy[g]*x28;
                    base_integral_17[g] += -2*x39 - PGy[g]*x38;
                    const T x40 = PGx[g]*x37;
                    base_integral_8[g] += -PGx[g]*x40 - x39;

                } // for g
            } // if/else prim
        } // for prim

        // copy the block to the output arrays
        static const short row[60] = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 3, 4, 6, 7, 8, 10, 11, 12, 13,
            2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 3, 6, 7, 10, 11, 12, 15, 16, 17, 18,
            4, 7, 8, 11, 12, 13, 16, 17, 18, 19, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20,
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
} // HRR_D_F_SAMECENTER

template void HRR_D_F_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_D_F_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

