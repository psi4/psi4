#include "HRR_P_D.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_P_D(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        T* base_integrals,         //!< work array for the base [LA+LB|0] integrals of the HRR (not used by all kernels)
        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)
        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively
        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PAx,            //!< (Px - Ax); P is the GPT combined center; A is the center of multishell A
        const T* v_PAy,            //!< (Px - Ax); P is the GPT combined center; A is the center of multishell A
        const T* v_PAz,            //!< (Px - Ax); P is the GPT combined center; A is the center of multishell A
        T ABx,                     //!< (Ax - Bx); A/B is the center of multishell A/B, respectively
        T ABy,                     //!< (Ax - Bx); A/B is the center of multishell A/B, respectively
        T ABz,                     //!< (Ax - Bx); A/B is the center of multishell A/B, respectively
        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair
        int gb_size)               //!< number of gridpoints in this grid batch 
{
    if (41*sizeof(T) >= 280) {

        if (nprim < 1) return;

        // The gridpoints are processed in blocks.  The integrals of a block are contracted in
        // a buffer that stays in the L1 cache and are then copied to the output arrays.
        enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
        alignas(64) T base_integral[16*g_block];
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
        alignas(64) T block[18*g_block];
        T* __restrict block_0 = &block[0*g_block];
        T* __restrict block_1 = &block[1*g_block];
        T* __restrict block_2 = &block[2*g_block];
        T* __restrict block_3 = &block[3*g_block];
        T* __restrict block_4 = &block[4*g_block];
        T* __restrict block_5 = &block[5*g_block];
        T* __restrict block_6 = &block[6*g_block];
        T* __restrict block_7 = &block[7*g_block];
        T* __restrict block_8 = &block[8*g_block];
        T* __restrict block_9 = &block[9*g_block];
        T* __restrict block_10 = &block[10*g_block];
        T* __restrict block_11 = &block[11*g_block];
        T* __restrict block_12 = &block[12*g_block];
        T* __restrict block_13 = &block[13*g_block];
        T* __restrict block_14 = &block[14*g_block];
        T* __restrict block_15 = &block[15*g_block];
        T* __restrict block_16 = &block[16*g_block];
        T* __restrict block_17 = &block[17*g_block];

        // fewer gridpoints than one block: work on a zero padded copy of the inputs
        const bool padded = (gb_size < g_block);
        alignas(64) T pad[7*g_block];
        if (padded) {
            for (int k = 0; k < 7*g_block; ++k) pad[k] = 0;
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
                        pad[3*g_block+k] = v_scaled_boys[prim*4*gb_size+0*gb_size+k];
                        pad[4*g_block+k] = v_scaled_boys[prim*4*gb_size+1*gb_size+k];
                        pad[5*g_block+k] = v_scaled_boys[prim*4*gb_size+2*gb_size+k];
                        pad[6*g_block+k] = v_scaled_boys[prim*4*gb_size+3*gb_size+k];
                    }
                }
                const T* __restrict PGx = padded ? &pad[0*g_block] : &v_PGx[prim*gb_size + g_lo];
                const T* __restrict PGy = padded ? &pad[1*g_block] : &v_PGy[prim*gb_size + g_lo];
                const T* __restrict PGz = padded ? &pad[2*g_block] : &v_PGz[prim*gb_size + g_lo];
                const T* __restrict scaled_boys_0 = padded ? &pad[3*g_block] : &v_scaled_boys[prim*4*gb_size+0*gb_size + g_lo];
                const T* __restrict scaled_boys_1 = padded ? &pad[4*g_block] : &v_scaled_boys[prim*4*gb_size+1*gb_size + g_lo];
                const T* __restrict scaled_boys_2 = padded ? &pad[5*g_block] : &v_scaled_boys[prim*4*gb_size+2*gb_size + g_lo];
                const T* __restrict scaled_boys_3 = padded ? &pad[6*g_block] : &v_scaled_boys[prim*4*gb_size+3*gb_size + g_lo];

                const T PBx = v_PAx[prim] + ABx;
                const T PBy = v_PAy[prim] + ABy;
                const T PBz = v_PAz[prim] + ABz;
                const T tpi = 0.5 / v_p[prim];
                const T tpi_2 = 2*tpi;

                // overwrite integrals for first primitive, accumulate otherwise
                if (prim==0) {
                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        base_integral_1[g] = PBx*x0 - PGx[g]*x1;
                        const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x3 = tpi*x2 + PBy*x0 - PGy[g]*x1;
                        base_integral_3[g] = x3;
                        const T x4 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x5 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x6 = PBy*x1 + tpi*x4 - PGy[g]*x5;
                        const T x7 = x0 - x1;
                        base_integral_12[g] = PBy*x3 + tpi_2*x7 - PGy[g]*x6;
                        base_integral_9[g] = PBx*x3 - PGx[g]*x6;
                        base_integral_13[g] = PBz*x3 - PGz[g]*x6;
                        const T x8 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x9 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        base_integral_4[g] = PBy*x8 - PGy[g]*x9;
                        const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x12 = PBx*x10 + tpi*x2 - PGx[g]*x11;
                        base_integral_0[g] = x12;
                        const T x13 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x14 = PBx*x11 + tpi*x4 - PGx[g]*x13;
                        const T x15 = x10 - x11;
                        base_integral_6[g] = PBx*x12 + tpi_2*x15 - PGx[g]*x14;
                        const T x16 = PBx*x8 - PGx[g]*x9;
                        base_integral_2[g] = x16;
                        const T x17 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x18 = PBx*x9 - PGx[g]*x17;
                        base_integral_10[g] = PBy*x16 - PGy[g]*x18;
                        const T x19 = tpi*x2 + PBz*x8 - PGz[g]*x9;
                        base_integral_5[g] = x19;
                        const T x20 = PBz*x9 + tpi*x4 - PGz[g]*x17;
                        const T x21 = x8 - x9;
                        base_integral_15[g] = PBz*x19 + tpi_2*x21 - PGz[g]*x20;
                        base_integral_7[g] = PBy*x12 - PGy[g]*x14;
                        base_integral_8[g] = PBz*x12 - PGz[g]*x14;
                        base_integral_11[g] = PBx*x19 - PGx[g]*x20;
                        base_integral_14[g] = PBy*x19 - PGy[g]*x20;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        base_integral_1[g] += PBx*x0 - PGx[g]*x1;
                        const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x3 = tpi*x2 + PBy*x0 - PGy[g]*x1;
                        base_integral_3[g] += x3;
                        const T x4 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x5 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x6 = PBy*x1 + tpi*x4 - PGy[g]*x5;
                        const T x7 = x0 - x1;
                        base_integral_12[g] += PBy*x3 + tpi_2*x7 - PGy[g]*x6;
                        base_integral_9[g] += PBx*x3 - PGx[g]*x6;
                        base_integral_13[g] += PBz*x3 - PGz[g]*x6;
                        const T x8 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x9 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        base_integral_4[g] += PBy*x8 - PGy[g]*x9;
                        const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x12 = PBx*x10 + tpi*x2 - PGx[g]*x11;
                        base_integral_0[g] += x12;
                        const T x13 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x14 = PBx*x11 + tpi*x4 - PGx[g]*x13;
                        const T x15 = x10 - x11;
                        base_integral_6[g] += PBx*x12 + tpi_2*x15 - PGx[g]*x14;
                        const T x16 = PBx*x8 - PGx[g]*x9;
                        base_integral_2[g] += x16;
                        const T x17 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x18 = PBx*x9 - PGx[g]*x17;
                        base_integral_10[g] += PBy*x16 - PGy[g]*x18;
                        const T x19 = tpi*x2 + PBz*x8 - PGz[g]*x9;
                        base_integral_5[g] += x19;
                        const T x20 = PBz*x9 + tpi*x4 - PGz[g]*x17;
                        const T x21 = x8 - x9;
                        base_integral_15[g] += PBz*x19 + tpi_2*x21 - PGz[g]*x20;
                        base_integral_7[g] += PBy*x12 - PGy[g]*x14;
                        base_integral_8[g] += PBz*x12 - PGz[g]*x14;
                        base_integral_11[g] += PBx*x19 - PGx[g]*x20;
                        base_integral_14[g] += PBy*x19 - PGy[g]*x20;

                    } // for g
                } // if/else prim
            } // for prim

            #pragma omp simd
            for (int g = 0; g < g_block; ++g) {

                block_0[g] = base_integral_6[g] - ABx*base_integral_0[g];
                block_6[g] = base_integral_7[g] - ABy*base_integral_0[g];
                block_1[g] = base_integral_7[g] - ABx*base_integral_1[g];
                block_12[g] = base_integral_8[g] - ABz*base_integral_0[g];
                block_2[g] = base_integral_8[g] - ABx*base_integral_2[g];
                block_7[g] = base_integral_9[g] - ABy*base_integral_1[g];
                block_3[g] = base_integral_9[g] - ABx*base_integral_3[g];
                block_9[g] = base_integral_12[g] - ABy*base_integral_3[g];
                block_13[g] = base_integral_10[g] - ABz*base_integral_1[g];
                block_8[g] = base_integral_10[g] - ABy*base_integral_2[g];
                block_4[g] = base_integral_10[g] - ABx*base_integral_4[g];
                block_14[g] = base_integral_11[g] - ABz*base_integral_2[g];
                block_5[g] = base_integral_11[g] - ABx*base_integral_5[g];
                block_15[g] = base_integral_13[g] - ABz*base_integral_3[g];
                block_10[g] = base_integral_13[g] - ABy*base_integral_4[g];
                block_16[g] = base_integral_14[g] - ABz*base_integral_4[g];
                block_11[g] = base_integral_14[g] - ABy*base_integral_5[g];
                block_17[g] = base_integral_15[g] - ABz*base_integral_5[g];

            } // for g

            // copy the block to the output arrays
            for (int i = 0; i < 18; ++i) {
                const T* __restrict src = &block[i*g_block];
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
        T* __restrict base_integral_15 = &base_integrals[15*gb_size];

        for (int prim = 0; prim < nprim; ++prim) {

            const T* __restrict PGx = &v_PGx[prim*gb_size];
            const T* __restrict PGy = &v_PGy[prim*gb_size];
            const T* __restrict PGz = &v_PGz[prim*gb_size];
            const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*4*gb_size+0*gb_size];
            const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*4*gb_size+1*gb_size];
            const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*4*gb_size+2*gb_size];
            const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*4*gb_size+3*gb_size];

            const T PBx = v_PAx[prim] + ABx;
            const T PBy = v_PAy[prim] + ABy;
            const T PBz = v_PAz[prim] + ABz;
            const T tpi = 0.5 / v_p[prim];
            const T tpi_2 = 2*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x3 = tpi*x2 + PBy*x0 - PGy[g]*x1;
                    const T x4 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x5 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x6 = PBy*x1 + tpi*x4 - PGy[g]*x5;
                    const T x7 = x0 - x1;
                    const T x8 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x9 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x12 = PBx*x10 + tpi*x2 - PGx[g]*x11;
                    const T x13 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x14 = PBx*x11 + tpi*x4 - PGx[g]*x13;
                    const T x15 = x10 - x11;
                    const T x16 = PBx*x8 - PGx[g]*x9;
                    const T x17 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x18 = PBx*x9 - PGx[g]*x17;
                    const T x19 = tpi*x2 + PBz*x8 - PGz[g]*x9;
                    const T x20 = PBz*x9 + tpi*x4 - PGz[g]*x17;
                    const T x21 = x8 - x9;

                    base_integral_0[g] = x12;
                    base_integral_1[g] = PBx*x0 - PGx[g]*x1;
                    base_integral_2[g] = x16;
                    base_integral_3[g] = x3;
                    base_integral_4[g] = PBy*x8 - PGy[g]*x9;
                    base_integral_5[g] = x19;
                    base_integral_6[g] = PBx*x12 + tpi_2*x15 - PGx[g]*x14;
                    base_integral_7[g] = PBy*x12 - PGy[g]*x14;
                    base_integral_8[g] = PBz*x12 - PGz[g]*x14;
                    base_integral_9[g] = PBx*x3 - PGx[g]*x6;
                    base_integral_10[g] = PBy*x16 - PGy[g]*x18;
                    base_integral_11[g] = PBx*x19 - PGx[g]*x20;
                    base_integral_12[g] = PBy*x3 + tpi_2*x7 - PGy[g]*x6;
                    base_integral_13[g] = PBz*x3 - PGz[g]*x6;
                    base_integral_14[g] = PBy*x19 - PGy[g]*x20;
                    base_integral_15[g] = PBz*x19 + tpi_2*x21 - PGz[g]*x20;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x3 = tpi*x2 + PBy*x0 - PGy[g]*x1;
                    const T x4 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x5 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x6 = PBy*x1 + tpi*x4 - PGy[g]*x5;
                    const T x7 = x0 - x1;
                    const T x8 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x9 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x12 = PBx*x10 + tpi*x2 - PGx[g]*x11;
                    const T x13 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x14 = PBx*x11 + tpi*x4 - PGx[g]*x13;
                    const T x15 = x10 - x11;
                    const T x16 = PBx*x8 - PGx[g]*x9;
                    const T x17 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x18 = PBx*x9 - PGx[g]*x17;
                    const T x19 = tpi*x2 + PBz*x8 - PGz[g]*x9;
                    const T x20 = PBz*x9 + tpi*x4 - PGz[g]*x17;
                    const T x21 = x8 - x9;

                    base_integral_0[g] += x12;
                    base_integral_1[g] += PBx*x0 - PGx[g]*x1;
                    base_integral_2[g] += x16;
                    base_integral_3[g] += x3;
                    base_integral_4[g] += PBy*x8 - PGy[g]*x9;
                    base_integral_5[g] += x19;
                    base_integral_6[g] += PBx*x12 + tpi_2*x15 - PGx[g]*x14;
                    base_integral_7[g] += PBy*x12 - PGy[g]*x14;
                    base_integral_8[g] += PBz*x12 - PGz[g]*x14;
                    base_integral_9[g] += PBx*x3 - PGx[g]*x6;
                    base_integral_10[g] += PBy*x16 - PGy[g]*x18;
                    base_integral_11[g] += PBx*x19 - PGx[g]*x20;
                    base_integral_12[g] += PBy*x3 + tpi_2*x7 - PGy[g]*x6;
                    base_integral_13[g] += PBz*x3 - PGz[g]*x6;
                    base_integral_14[g] += PBy*x19 - PGy[g]*x20;
                    base_integral_15[g] += PBz*x19 + tpi_2*x21 - PGz[g]*x20;

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

        #pragma omp simd
        for (int g = 0; g < gb_size; ++g) {


            integral_0[g] = base_integral_6[g] - ABx*base_integral_0[g];
            integral_1[g] = base_integral_7[g] - ABx*base_integral_1[g];
            integral_2[g] = base_integral_8[g] - ABx*base_integral_2[g];
            integral_3[g] = base_integral_9[g] - ABx*base_integral_3[g];
            integral_4[g] = base_integral_10[g] - ABx*base_integral_4[g];
            integral_5[g] = base_integral_11[g] - ABx*base_integral_5[g];
            integral_6[g] = base_integral_7[g] - ABy*base_integral_0[g];
            integral_7[g] = base_integral_9[g] - ABy*base_integral_1[g];
            integral_8[g] = base_integral_10[g] - ABy*base_integral_2[g];
            integral_9[g] = base_integral_12[g] - ABy*base_integral_3[g];
            integral_10[g] = base_integral_13[g] - ABy*base_integral_4[g];
            integral_11[g] = base_integral_14[g] - ABy*base_integral_5[g];
            integral_12[g] = base_integral_8[g] - ABz*base_integral_0[g];
            integral_13[g] = base_integral_10[g] - ABz*base_integral_1[g];
            integral_14[g] = base_integral_11[g] - ABz*base_integral_2[g];
            integral_15[g] = base_integral_13[g] - ABz*base_integral_3[g];
            integral_16[g] = base_integral_14[g] - ABz*base_integral_4[g];
            integral_17[g] = base_integral_15[g] - ABz*base_integral_5[g];

        } // for g

    }
} // HRR_P_D

template void HRR_P_D<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        const float* v_PAx, 
        const float* v_PAy, 
        const float* v_PAz, 
        float ABx, 
        float ABy, 
        float ABz, 
        int nprim, 
        int gb_size);

template void HRR_P_D<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        const double* v_PAx, 
        const double* v_PAy, 
        const double* v_PAz, 
        double ABx, 
        double ABy, 
        double ABz, 
        int nprim, 
        int gb_size);

