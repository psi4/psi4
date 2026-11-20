#include "HRR_P_F_SAMECENTER.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_P_F_SAMECENTER(
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
    if (53*sizeof(T) >= 280) {

        if (nprim < 1) return;

        // The gridpoints are processed in blocks.  The integrals of a block are contracted in
        // a buffer that stays in the L1 cache and are then copied to the output arrays.
        enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
        alignas(64) T base_integral[15*g_block];
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

        // fewer gridpoints than one block: work on a zero padded copy of the inputs
        const bool padded = (gb_size < g_block);
        alignas(64) T pad[8*g_block];
        if (padded) {
            for (int k = 0; k < 8*g_block; ++k) pad[k] = 0;
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
                        pad[3*g_block+k] = v_scaled_boys[prim*5*gb_size+0*gb_size+k];
                        pad[4*g_block+k] = v_scaled_boys[prim*5*gb_size+1*gb_size+k];
                        pad[5*g_block+k] = v_scaled_boys[prim*5*gb_size+2*gb_size+k];
                        pad[6*g_block+k] = v_scaled_boys[prim*5*gb_size+3*gb_size+k];
                        pad[7*g_block+k] = v_scaled_boys[prim*5*gb_size+4*gb_size+k];
                    }
                }
                const T* __restrict PGx = padded ? &pad[0*g_block] : &v_PGx[prim*gb_size + g_lo];
                const T* __restrict PGy = padded ? &pad[1*g_block] : &v_PGy[prim*gb_size + g_lo];
                const T* __restrict PGz = padded ? &pad[2*g_block] : &v_PGz[prim*gb_size + g_lo];
                const T* __restrict scaled_boys_0 = padded ? &pad[3*g_block] : &v_scaled_boys[prim*5*gb_size+0*gb_size + g_lo];
                const T* __restrict scaled_boys_1 = padded ? &pad[4*g_block] : &v_scaled_boys[prim*5*gb_size+1*gb_size + g_lo];
                const T* __restrict scaled_boys_2 = padded ? &pad[5*g_block] : &v_scaled_boys[prim*5*gb_size+2*gb_size + g_lo];
                const T* __restrict scaled_boys_3 = padded ? &pad[6*g_block] : &v_scaled_boys[prim*5*gb_size+3*gb_size + g_lo];
                const T* __restrict scaled_boys_4 = padded ? &pad[7*g_block] : &v_scaled_boys[prim*5*gb_size+4*gb_size + g_lo];

                const T tpi = 0.5 / v_p[prim];
                const T tpi_p2 = tpi*tpi;

                // overwrite integrals for first primitive, accumulate otherwise
                if (prim==0) {
                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x1 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x2 = x0 - x1;
                        const T x3 = tpi_p2*x2;
                        const T x4 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x5 = tpi*x4;
                        const T x6 = PGy[g]*x5;
                        const T x7 = x3 + PGy[g]*x6;
                        const T x8 = PGy[g]*scaled_boys_4[g];
                        const T x9 = x5 + PGy[g]*x8;
                        const T x10 = PGy[g]*x9 + 2*x6;
                        base_integral_10[g] = 3*x7 + PGy[g]*x10;
                        const T x11 = PGx[g]*x9;
                        base_integral_3[g] = PGx[g]*x11 + x7;
                        base_integral_7[g] = PGz[g]*x11;
                        base_integral_6[g] = PGx[g]*x10;
                        base_integral_11[g] = PGz[g]*x10;
                        const T x12 = PGx[g]*scaled_boys_4[g];
                        const T x13 = PGx[g]*x12 + x5;
                        const T x14 = PGz[g]*x13;
                        base_integral_4[g] = PGy[g]*x14;
                        const T x15 = PGx[g]*x5;
                        const T x16 = PGx[g]*x15 + x3;
                        const T x17 = PGx[g]*x13 + 2*x15;
                        base_integral_5[g] = x16 + PGz[g]*x14;
                        base_integral_0[g] = PGx[g]*x17 + 3*x16;
                        base_integral_1[g] = PGy[g]*x17;
                        base_integral_2[g] = PGz[g]*x17;
                        const T x18 = PGz[g]*x5;
                        const T x19 = x3 + PGz[g]*x18;
                        const T x20 = PGz[g]*scaled_boys_4[g];
                        const T x21 = x5 + PGz[g]*x20;
                        const T x22 = PGz[g]*x21 + 2*x18;
                        base_integral_14[g] = PGz[g]*x22 + 3*x19;
                        const T x23 = PGy[g]*x21;
                        base_integral_12[g] = PGy[g]*x23 + x19;
                        base_integral_8[g] = PGx[g]*x23;
                        base_integral_9[g] = PGx[g]*x22;
                        base_integral_13[g] = PGy[g]*x22;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x1 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x2 = x0 - x1;
                        const T x3 = tpi_p2*x2;
                        const T x4 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x5 = tpi*x4;
                        const T x6 = PGy[g]*x5;
                        const T x7 = x3 + PGy[g]*x6;
                        const T x8 = PGy[g]*scaled_boys_4[g];
                        const T x9 = x5 + PGy[g]*x8;
                        const T x10 = PGy[g]*x9 + 2*x6;
                        base_integral_10[g] += 3*x7 + PGy[g]*x10;
                        const T x11 = PGx[g]*x9;
                        base_integral_3[g] += PGx[g]*x11 + x7;
                        base_integral_7[g] += PGz[g]*x11;
                        base_integral_6[g] += PGx[g]*x10;
                        base_integral_11[g] += PGz[g]*x10;
                        const T x12 = PGx[g]*scaled_boys_4[g];
                        const T x13 = PGx[g]*x12 + x5;
                        const T x14 = PGz[g]*x13;
                        base_integral_4[g] += PGy[g]*x14;
                        const T x15 = PGx[g]*x5;
                        const T x16 = PGx[g]*x15 + x3;
                        const T x17 = PGx[g]*x13 + 2*x15;
                        base_integral_5[g] += x16 + PGz[g]*x14;
                        base_integral_0[g] += PGx[g]*x17 + 3*x16;
                        base_integral_1[g] += PGy[g]*x17;
                        base_integral_2[g] += PGz[g]*x17;
                        const T x18 = PGz[g]*x5;
                        const T x19 = x3 + PGz[g]*x18;
                        const T x20 = PGz[g]*scaled_boys_4[g];
                        const T x21 = x5 + PGz[g]*x20;
                        const T x22 = PGz[g]*x21 + 2*x18;
                        base_integral_14[g] += PGz[g]*x22 + 3*x19;
                        const T x23 = PGy[g]*x21;
                        base_integral_12[g] += PGy[g]*x23 + x19;
                        base_integral_8[g] += PGx[g]*x23;
                        base_integral_9[g] += PGx[g]*x22;
                        base_integral_13[g] += PGy[g]*x22;

                    } // for g
                } // if/else prim
            } // for prim

            // copy the block to the output arrays
            static const short row[30] = {
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 3, 4, 6, 7, 8, 10, 11, 12, 13,
                2, 4, 5, 7, 8, 9, 11, 12, 13, 14,
            };
            for (int i = 0; i < 30; ++i) {
                const T* __restrict src = &base_integral[row[i]*g_block];
                T* __restrict dst = &integrals[i*gb_size + g_lo];
                if (padded) {
                    SNK_COPY(dst, src, gb_size);
                } else {
                    SNK_COPY(dst, src, g_block);
                }
            }
        } // for g_lo

    } else {

        T* __restrict base_integral_0 = &base_integrals[0*gb_size];
        T* __restrict base_integral_1 = &base_integrals[1*gb_size];
        T* __restrict base_integral_2 = &base_integrals[2*gb_size];
        T* __restrict base_integral_3 = &base_integrals[3*gb_size];
        T* __restrict base_integral_4 = &base_integrals[4*gb_size];
        T* __restrict base_integral_5 = &base_integrals[5*gb_size];
        T* __restrict base_integral_6 = &base_integrals[6*gb_size];
        T* __restrict base_integral_7 = &base_integrals[7*gb_size];
        T* __restrict base_integral_8 = &base_integrals[8*gb_size];
        T* __restrict base_integral_9 = &base_integrals[9*gb_size];
        T* __restrict base_integral_10 = &base_integrals[10*gb_size];
        T* __restrict base_integral_11 = &base_integrals[11*gb_size];
        T* __restrict base_integral_12 = &base_integrals[12*gb_size];
        T* __restrict base_integral_13 = &base_integrals[13*gb_size];
        T* __restrict base_integral_14 = &base_integrals[14*gb_size];

        for (int prim = 0; prim < nprim; ++prim) {

            const T* __restrict PGx = &v_PGx[prim*gb_size];
            const T* __restrict PGy = &v_PGy[prim*gb_size];
            const T* __restrict PGz = &v_PGz[prim*gb_size];
            const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*5*gb_size+0*gb_size];
            const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*5*gb_size+1*gb_size];
            const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*5*gb_size+2*gb_size];
            const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*5*gb_size+3*gb_size];
            const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*5*gb_size+4*gb_size];

            const T tpi = 0.5 / v_p[prim];
            const T tpi_p2 = tpi*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x1 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x2 = x0 - x1;
                    const T x3 = tpi_p2*x2;
                    const T x4 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x5 = tpi*x4;
                    const T x6 = PGy[g]*x5;
                    const T x7 = x3 + PGy[g]*x6;
                    const T x8 = PGy[g]*scaled_boys_4[g];
                    const T x9 = x5 + PGy[g]*x8;
                    const T x10 = PGy[g]*x9 + 2*x6;
                    const T x11 = PGx[g]*x9;
                    const T x12 = PGx[g]*scaled_boys_4[g];
                    const T x13 = PGx[g]*x12 + x5;
                    const T x14 = PGz[g]*x13;
                    const T x15 = PGx[g]*x5;
                    const T x16 = PGx[g]*x15 + x3;
                    const T x17 = PGx[g]*x13 + 2*x15;
                    const T x18 = PGz[g]*x5;
                    const T x19 = x3 + PGz[g]*x18;
                    const T x20 = PGz[g]*scaled_boys_4[g];
                    const T x21 = x5 + PGz[g]*x20;
                    const T x22 = PGz[g]*x21 + 2*x18;
                    const T x23 = PGy[g]*x21;

                    base_integral_0[g] = PGx[g]*x17 + 3*x16;
                    base_integral_1[g] = PGy[g]*x17;
                    base_integral_2[g] = PGz[g]*x17;
                    base_integral_3[g] = PGx[g]*x11 + x7;
                    base_integral_4[g] = PGy[g]*x14;
                    base_integral_5[g] = x16 + PGz[g]*x14;
                    base_integral_6[g] = PGx[g]*x10;
                    base_integral_7[g] = PGz[g]*x11;
                    base_integral_8[g] = PGx[g]*x23;
                    base_integral_9[g] = PGx[g]*x22;
                    base_integral_10[g] = 3*x7 + PGy[g]*x10;
                    base_integral_11[g] = PGz[g]*x10;
                    base_integral_12[g] = PGy[g]*x23 + x19;
                    base_integral_13[g] = PGy[g]*x22;
                    base_integral_14[g] = PGz[g]*x22 + 3*x19;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x1 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x2 = x0 - x1;
                    const T x3 = tpi_p2*x2;
                    const T x4 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x5 = tpi*x4;
                    const T x6 = PGy[g]*x5;
                    const T x7 = x3 + PGy[g]*x6;
                    const T x8 = PGy[g]*scaled_boys_4[g];
                    const T x9 = x5 + PGy[g]*x8;
                    const T x10 = PGy[g]*x9 + 2*x6;
                    const T x11 = PGx[g]*x9;
                    const T x12 = PGx[g]*scaled_boys_4[g];
                    const T x13 = PGx[g]*x12 + x5;
                    const T x14 = PGz[g]*x13;
                    const T x15 = PGx[g]*x5;
                    const T x16 = PGx[g]*x15 + x3;
                    const T x17 = PGx[g]*x13 + 2*x15;
                    const T x18 = PGz[g]*x5;
                    const T x19 = x3 + PGz[g]*x18;
                    const T x20 = PGz[g]*scaled_boys_4[g];
                    const T x21 = x5 + PGz[g]*x20;
                    const T x22 = PGz[g]*x21 + 2*x18;
                    const T x23 = PGy[g]*x21;

                    base_integral_0[g] += PGx[g]*x17 + 3*x16;
                    base_integral_1[g] += PGy[g]*x17;
                    base_integral_2[g] += PGz[g]*x17;
                    base_integral_3[g] += PGx[g]*x11 + x7;
                    base_integral_4[g] += PGy[g]*x14;
                    base_integral_5[g] += x16 + PGz[g]*x14;
                    base_integral_6[g] += PGx[g]*x10;
                    base_integral_7[g] += PGz[g]*x11;
                    base_integral_8[g] += PGx[g]*x23;
                    base_integral_9[g] += PGx[g]*x22;
                    base_integral_10[g] += 3*x7 + PGy[g]*x10;
                    base_integral_11[g] += PGz[g]*x10;
                    base_integral_12[g] += PGy[g]*x23 + x19;
                    base_integral_13[g] += PGy[g]*x22;
                    base_integral_14[g] += PGz[g]*x22 + 3*x19;

                } // for g
            } // if/else prim
        } // for prim

        T* __restrict integral_0 = &integrals[0*gb_size];
        T* __restrict integral_1 = &integrals[1*gb_size];
        T* __restrict integral_2 = &integrals[2*gb_size];
        T* __restrict integral_3 = &integrals[3*gb_size];
        T* __restrict integral_4 = &integrals[4*gb_size];
        T* __restrict integral_5 = &integrals[5*gb_size];
        T* __restrict integral_6 = &integrals[6*gb_size];
        T* __restrict integral_7 = &integrals[7*gb_size];
        T* __restrict integral_8 = &integrals[8*gb_size];
        T* __restrict integral_9 = &integrals[9*gb_size];
        T* __restrict integral_10 = &integrals[10*gb_size];
        T* __restrict integral_11 = &integrals[11*gb_size];
        T* __restrict integral_12 = &integrals[12*gb_size];
        T* __restrict integral_13 = &integrals[13*gb_size];
        T* __restrict integral_14 = &integrals[14*gb_size];
        T* __restrict integral_15 = &integrals[15*gb_size];
        T* __restrict integral_16 = &integrals[16*gb_size];
        T* __restrict integral_17 = &integrals[17*gb_size];
        T* __restrict integral_18 = &integrals[18*gb_size];
        T* __restrict integral_19 = &integrals[19*gb_size];
        T* __restrict integral_20 = &integrals[20*gb_size];
        T* __restrict integral_21 = &integrals[21*gb_size];
        T* __restrict integral_22 = &integrals[22*gb_size];
        T* __restrict integral_23 = &integrals[23*gb_size];
        T* __restrict integral_24 = &integrals[24*gb_size];
        T* __restrict integral_25 = &integrals[25*gb_size];
        T* __restrict integral_26 = &integrals[26*gb_size];
        T* __restrict integral_27 = &integrals[27*gb_size];
        T* __restrict integral_28 = &integrals[28*gb_size];
        T* __restrict integral_29 = &integrals[29*gb_size];

        #pragma omp simd
        for (int g = 0; g < gb_size; ++g) {


            integral_0[g] = base_integral_0[g];
            integral_1[g] = base_integral_1[g];
            integral_2[g] = base_integral_2[g];
            integral_3[g] = base_integral_3[g];
            integral_4[g] = base_integral_4[g];
            integral_5[g] = base_integral_5[g];
            integral_6[g] = base_integral_6[g];
            integral_7[g] = base_integral_7[g];
            integral_8[g] = base_integral_8[g];
            integral_9[g] = base_integral_9[g];
            integral_10[g] = base_integral_1[g];
            integral_11[g] = base_integral_3[g];
            integral_12[g] = base_integral_4[g];
            integral_13[g] = base_integral_6[g];
            integral_14[g] = base_integral_7[g];
            integral_15[g] = base_integral_8[g];
            integral_16[g] = base_integral_10[g];
            integral_17[g] = base_integral_11[g];
            integral_18[g] = base_integral_12[g];
            integral_19[g] = base_integral_13[g];
            integral_20[g] = base_integral_2[g];
            integral_21[g] = base_integral_4[g];
            integral_22[g] = base_integral_5[g];
            integral_23[g] = base_integral_7[g];
            integral_24[g] = base_integral_8[g];
            integral_25[g] = base_integral_9[g];
            integral_26[g] = base_integral_11[g];
            integral_27[g] = base_integral_12[g];
            integral_28[g] = base_integral_13[g];
            integral_29[g] = base_integral_14[g];

        } // for g

    }
} // HRR_P_F_SAMECENTER

template void HRR_P_F_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_P_F_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

