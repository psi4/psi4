#include "HRR_F_P.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_F_P(
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
    if (63*sizeof(T) >= 280) {

        if (nprim < 1) return;

        // The gridpoints are processed in blocks.  The integrals of a block are contracted in
        // a buffer that stays in the L1 cache and are then copied to the output arrays.
        enum { g_block = (sizeof(T) <= 4) ? 32 : 16 };
        alignas(64) T base_integral[25*g_block];
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
        T* __restrict base_integral_21 = &base_integral[21*g_block];
        T* __restrict base_integral_22 = &base_integral[22*g_block];
        T* __restrict base_integral_23 = &base_integral[23*g_block];
        T* __restrict base_integral_24 = &base_integral[24*g_block];
        alignas(64) T block[30*g_block];
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
        T* __restrict block_18 = &block[18*g_block];
        T* __restrict block_19 = &block[19*g_block];
        T* __restrict block_20 = &block[20*g_block];
        T* __restrict block_21 = &block[21*g_block];
        T* __restrict block_22 = &block[22*g_block];
        T* __restrict block_23 = &block[23*g_block];
        T* __restrict block_24 = &block[24*g_block];
        T* __restrict block_25 = &block[25*g_block];
        T* __restrict block_26 = &block[26*g_block];
        T* __restrict block_27 = &block[27*g_block];
        T* __restrict block_28 = &block[28*g_block];
        T* __restrict block_29 = &block[29*g_block];

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

                const T PAx = v_PAx[prim];
                const T PAy = v_PAy[prim];
                const T PAz = v_PAz[prim];
                const T tpi = 0.5 / v_p[prim];
                const T tpi_2 = 2*tpi;
                const T tpi_3 = 3*tpi;

                // overwrite integrals for first primitive, accumulate otherwise
                if (prim==0) {
                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x2 = PAx*x0 - PGx[g]*x1;
                        const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x4 = PAx*x1 - PGx[g]*x3;
                        base_integral_4[g] = PAy*x2 - PGy[g]*x4;
                        const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x6 = tpi*x5 + PAz*x0 - PGz[g]*x1;
                        const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x8 = tpi*x7 + PAz*x1 - PGz[g]*x3;
                        base_integral_5[g] = PAx*x6 - PGx[g]*x8;
                        const T x9 = x0 - x1;
                        const T x10 = PAz*x6 + tpi_2*x9 - PGz[g]*x8;
                        base_integral_9[g] = x10;
                        const T x11 = PAy*x6 - PGy[g]*x8;
                        base_integral_8[g] = x11;
                        const T x12 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x13 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x14 = PAz*x3 + tpi*x12 - PGz[g]*x13;
                        const T x15 = x1 - x3;
                        const T x16 = PAz*x8 + tpi_2*x15 - PGz[g]*x14;
                        base_integral_19[g] = PAx*x10 - PGx[g]*x16;
                        const T x17 = PAy*x8 - PGy[g]*x14;
                        const T x18 = x6 - x8;
                        base_integral_22[g] = PAy*x11 + tpi*x18 - PGy[g]*x17;
                        base_integral_18[g] = PAx*x11 - PGx[g]*x17;
                        base_integral_24[g] = PAz*x10 + tpi_3*x18 - PGz[g]*x16;
                        base_integral_23[g] = PAy*x10 - PGy[g]*x16;
                        const T x19 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x20 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x21 = PAx*x19 + tpi*x5 - PGx[g]*x20;
                        const T x22 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x23 = PAx*x20 + tpi*x7 - PGx[g]*x22;
                        base_integral_1[g] = PAy*x21 - PGy[g]*x23;
                        const T x24 = x19 - x20;
                        const T x25 = PAx*x21 + tpi_2*x24 - PGx[g]*x23;
                        base_integral_0[g] = x25;
                        const T x26 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x27 = PAx*x22 + tpi*x12 - PGx[g]*x26;
                        const T x28 = x20 - x22;
                        const T x29 = PAx*x23 + tpi_2*x28 - PGx[g]*x27;
                        base_integral_11[g] = PAy*x25 - PGy[g]*x29;
                        const T x30 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x31 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x32 = tpi*x5 + PAy*x30 - PGy[g]*x31;
                        const T x33 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x34 = tpi*x7 + PAy*x31 - PGy[g]*x33;
                        const T x35 = x30 - x31;
                        const T x36 = PAy*x32 + tpi_2*x35 - PGy[g]*x34;
                        base_integral_6[g] = x36;
                        const T x37 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x38 = PAy*x33 + tpi*x12 - PGy[g]*x37;
                        const T x39 = x31 - x33;
                        const T x40 = PAy*x34 + tpi_2*x39 - PGy[g]*x38;
                        base_integral_16[g] = PAx*x36 - PGx[g]*x40;
                        base_integral_7[g] = PAz*x32 - PGz[g]*x34;
                        base_integral_12[g] = PAz*x25 - PGz[g]*x29;
                        const T x41 = x21 - x23;
                        base_integral_10[g] = PAx*x25 + tpi_3*x41 - PGx[g]*x29;
                        const T x42 = PAz*x21 - PGz[g]*x23;
                        base_integral_2[g] = x42;
                        const T x43 = PAz*x23 - PGz[g]*x27;
                        base_integral_15[g] = PAz*x42 + tpi*x41 - PGz[g]*x43;
                        base_integral_14[g] = PAy*x42 - PGy[g]*x43;
                        base_integral_21[g] = PAz*x36 - PGz[g]*x40;
                        const T x44 = x32 - x34;
                        base_integral_20[g] = PAy*x36 + tpi_3*x44 - PGy[g]*x40;
                        const T x45 = PAx*x32 - PGx[g]*x34;
                        base_integral_3[g] = x45;
                        const T x46 = PAx*x34 - PGx[g]*x38;
                        base_integral_13[g] = PAx*x45 + tpi*x44 - PGx[g]*x46;
                        base_integral_17[g] = PAz*x45 - PGz[g]*x46;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x2 = PAx*x0 - PGx[g]*x1;
                        const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x4 = PAx*x1 - PGx[g]*x3;
                        base_integral_4[g] += PAy*x2 - PGy[g]*x4;
                        const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x6 = tpi*x5 + PAz*x0 - PGz[g]*x1;
                        const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x8 = tpi*x7 + PAz*x1 - PGz[g]*x3;
                        base_integral_5[g] += PAx*x6 - PGx[g]*x8;
                        const T x9 = x0 - x1;
                        const T x10 = PAz*x6 + tpi_2*x9 - PGz[g]*x8;
                        base_integral_9[g] += x10;
                        const T x11 = PAy*x6 - PGy[g]*x8;
                        base_integral_8[g] += x11;
                        const T x12 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x13 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x14 = PAz*x3 + tpi*x12 - PGz[g]*x13;
                        const T x15 = x1 - x3;
                        const T x16 = PAz*x8 + tpi_2*x15 - PGz[g]*x14;
                        base_integral_19[g] += PAx*x10 - PGx[g]*x16;
                        const T x17 = PAy*x8 - PGy[g]*x14;
                        const T x18 = x6 - x8;
                        base_integral_22[g] += PAy*x11 + tpi*x18 - PGy[g]*x17;
                        base_integral_18[g] += PAx*x11 - PGx[g]*x17;
                        base_integral_24[g] += PAz*x10 + tpi_3*x18 - PGz[g]*x16;
                        base_integral_23[g] += PAy*x10 - PGy[g]*x16;
                        const T x19 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x20 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x21 = PAx*x19 + tpi*x5 - PGx[g]*x20;
                        const T x22 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x23 = PAx*x20 + tpi*x7 - PGx[g]*x22;
                        base_integral_1[g] += PAy*x21 - PGy[g]*x23;
                        const T x24 = x19 - x20;
                        const T x25 = PAx*x21 + tpi_2*x24 - PGx[g]*x23;
                        base_integral_0[g] += x25;
                        const T x26 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x27 = PAx*x22 + tpi*x12 - PGx[g]*x26;
                        const T x28 = x20 - x22;
                        const T x29 = PAx*x23 + tpi_2*x28 - PGx[g]*x27;
                        base_integral_11[g] += PAy*x25 - PGy[g]*x29;
                        const T x30 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x31 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x32 = tpi*x5 + PAy*x30 - PGy[g]*x31;
                        const T x33 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x34 = tpi*x7 + PAy*x31 - PGy[g]*x33;
                        const T x35 = x30 - x31;
                        const T x36 = PAy*x32 + tpi_2*x35 - PGy[g]*x34;
                        base_integral_6[g] += x36;
                        const T x37 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x38 = PAy*x33 + tpi*x12 - PGy[g]*x37;
                        const T x39 = x31 - x33;
                        const T x40 = PAy*x34 + tpi_2*x39 - PGy[g]*x38;
                        base_integral_16[g] += PAx*x36 - PGx[g]*x40;
                        base_integral_7[g] += PAz*x32 - PGz[g]*x34;
                        base_integral_12[g] += PAz*x25 - PGz[g]*x29;
                        const T x41 = x21 - x23;
                        base_integral_10[g] += PAx*x25 + tpi_3*x41 - PGx[g]*x29;
                        const T x42 = PAz*x21 - PGz[g]*x23;
                        base_integral_2[g] += x42;
                        const T x43 = PAz*x23 - PGz[g]*x27;
                        base_integral_15[g] += PAz*x42 + tpi*x41 - PGz[g]*x43;
                        base_integral_14[g] += PAy*x42 - PGy[g]*x43;
                        base_integral_21[g] += PAz*x36 - PGz[g]*x40;
                        const T x44 = x32 - x34;
                        base_integral_20[g] += PAy*x36 + tpi_3*x44 - PGy[g]*x40;
                        const T x45 = PAx*x32 - PGx[g]*x34;
                        base_integral_3[g] += x45;
                        const T x46 = PAx*x34 - PGx[g]*x38;
                        base_integral_13[g] += PAx*x45 + tpi*x44 - PGx[g]*x46;
                        base_integral_17[g] += PAz*x45 - PGz[g]*x46;

                    } // for g
                } // if/else prim
            } // for prim

            #pragma omp simd
            for (int g = 0; g < g_block; ++g) {

                block_10[g] = ABy*base_integral_3[g] + base_integral_16[g];
                block_18[g] = base_integral_16[g] + ABx*base_integral_6[g];
                block_19[g] = ABy*base_integral_6[g] + base_integral_20[g];
                block_20[g] = ABz*base_integral_6[g] + base_integral_21[g];
                block_22[g] = base_integral_21[g] + ABy*base_integral_7[g];
                block_9[g] = base_integral_13[g] + ABx*base_integral_3[g];
                block_4[g] = ABy*base_integral_1[g] + base_integral_13[g];
                block_11[g] = ABz*base_integral_3[g] + base_integral_17[g];
                block_21[g] = base_integral_17[g] + ABx*base_integral_7[g];
                block_13[g] = base_integral_17[g] + ABy*base_integral_4[g];
                block_23[g] = ABz*base_integral_7[g] + base_integral_22[g];
                block_25[g] = base_integral_22[g] + ABy*base_integral_8[g];
                block_3[g] = base_integral_11[g] + ABx*base_integral_1[g];
                block_1[g] = ABy*base_integral_0[g] + base_integral_11[g];
                block_0[g] = base_integral_10[g] + ABx*base_integral_0[g];
                block_2[g] = ABz*base_integral_0[g] + base_integral_12[g];
                block_5[g] = ABz*base_integral_1[g] + base_integral_14[g];
                block_12[g] = base_integral_14[g] + ABx*base_integral_4[g];
                block_6[g] = base_integral_12[g] + ABx*base_integral_2[g];
                block_7[g] = base_integral_14[g] + ABy*base_integral_2[g];
                block_8[g] = ABz*base_integral_2[g] + base_integral_15[g];
                block_14[g] = ABz*base_integral_4[g] + base_integral_18[g];
                block_24[g] = base_integral_18[g] + ABx*base_integral_8[g];
                block_15[g] = base_integral_15[g] + ABx*base_integral_5[g];
                block_16[g] = base_integral_18[g] + ABy*base_integral_5[g];
                block_17[g] = ABz*base_integral_5[g] + base_integral_19[g];
                block_26[g] = ABz*base_integral_8[g] + base_integral_23[g];
                block_27[g] = base_integral_19[g] + ABx*base_integral_9[g];
                block_28[g] = base_integral_23[g] + ABy*base_integral_9[g];
                block_29[g] = ABz*base_integral_9[g] + base_integral_24[g];

            } // for g

            // copy the block to the output arrays
            for (int i = 0; i < 30; ++i) {
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
        T* __restrict base_integral_16 = &base_integrals[16*gb_size];
        T* __restrict base_integral_17 = &base_integrals[17*gb_size];
        T* __restrict base_integral_18 = &base_integrals[18*gb_size];
        T* __restrict base_integral_19 = &base_integrals[19*gb_size];
        T* __restrict base_integral_20 = &base_integrals[20*gb_size];
        T* __restrict base_integral_21 = &base_integrals[21*gb_size];
        T* __restrict base_integral_22 = &base_integrals[22*gb_size];
        T* __restrict base_integral_23 = &base_integrals[23*gb_size];
        T* __restrict base_integral_24 = &base_integrals[24*gb_size];

        for (int prim = 0; prim < nprim; ++prim) {

            const T* __restrict PGx = &v_PGx[prim*gb_size];
            const T* __restrict PGy = &v_PGy[prim*gb_size];
            const T* __restrict PGz = &v_PGz[prim*gb_size];
            const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*5*gb_size+0*gb_size];
            const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*5*gb_size+1*gb_size];
            const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*5*gb_size+2*gb_size];
            const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*5*gb_size+3*gb_size];
            const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*5*gb_size+4*gb_size];

            const T PAx = v_PAx[prim];
            const T PAy = v_PAy[prim];
            const T PAz = v_PAz[prim];
            const T tpi = 0.5 / v_p[prim];
            const T tpi_2 = 2*tpi;
            const T tpi_3 = 3*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PAx*x0 - PGx[g]*x1;
                    const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PAx*x1 - PGx[g]*x3;
                    const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x6 = tpi*x5 + PAz*x0 - PGz[g]*x1;
                    const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x8 = tpi*x7 + PAz*x1 - PGz[g]*x3;
                    const T x9 = x0 - x1;
                    const T x10 = PAz*x6 + tpi_2*x9 - PGz[g]*x8;
                    const T x11 = PAy*x6 - PGy[g]*x8;
                    const T x12 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x13 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x14 = PAz*x3 + tpi*x12 - PGz[g]*x13;
                    const T x15 = x1 - x3;
                    const T x16 = PAz*x8 + tpi_2*x15 - PGz[g]*x14;
                    const T x17 = PAy*x8 - PGy[g]*x14;
                    const T x18 = x6 - x8;
                    const T x19 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x20 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x21 = PAx*x19 + tpi*x5 - PGx[g]*x20;
                    const T x22 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x23 = PAx*x20 + tpi*x7 - PGx[g]*x22;
                    const T x24 = x19 - x20;
                    const T x25 = PAx*x21 + tpi_2*x24 - PGx[g]*x23;
                    const T x26 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x27 = PAx*x22 + tpi*x12 - PGx[g]*x26;
                    const T x28 = x20 - x22;
                    const T x29 = PAx*x23 + tpi_2*x28 - PGx[g]*x27;
                    const T x30 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x31 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x32 = tpi*x5 + PAy*x30 - PGy[g]*x31;
                    const T x33 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x34 = tpi*x7 + PAy*x31 - PGy[g]*x33;
                    const T x35 = x30 - x31;
                    const T x36 = PAy*x32 + tpi_2*x35 - PGy[g]*x34;
                    const T x37 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x38 = PAy*x33 + tpi*x12 - PGy[g]*x37;
                    const T x39 = x31 - x33;
                    const T x40 = PAy*x34 + tpi_2*x39 - PGy[g]*x38;
                    const T x41 = x21 - x23;
                    const T x42 = PAz*x21 - PGz[g]*x23;
                    const T x43 = PAz*x23 - PGz[g]*x27;
                    const T x44 = x32 - x34;
                    const T x45 = PAx*x32 - PGx[g]*x34;
                    const T x46 = PAx*x34 - PGx[g]*x38;

                    base_integral_0[g] = x25;
                    base_integral_1[g] = PAy*x21 - PGy[g]*x23;
                    base_integral_2[g] = x42;
                    base_integral_3[g] = x45;
                    base_integral_4[g] = PAy*x2 - PGy[g]*x4;
                    base_integral_5[g] = PAx*x6 - PGx[g]*x8;
                    base_integral_6[g] = x36;
                    base_integral_7[g] = PAz*x32 - PGz[g]*x34;
                    base_integral_8[g] = x11;
                    base_integral_9[g] = x10;
                    base_integral_10[g] = PAx*x25 + tpi_3*x41 - PGx[g]*x29;
                    base_integral_11[g] = PAy*x25 - PGy[g]*x29;
                    base_integral_12[g] = PAz*x25 - PGz[g]*x29;
                    base_integral_13[g] = PAx*x45 + tpi*x44 - PGx[g]*x46;
                    base_integral_14[g] = PAy*x42 - PGy[g]*x43;
                    base_integral_15[g] = PAz*x42 + tpi*x41 - PGz[g]*x43;
                    base_integral_16[g] = PAx*x36 - PGx[g]*x40;
                    base_integral_17[g] = PAz*x45 - PGz[g]*x46;
                    base_integral_18[g] = PAx*x11 - PGx[g]*x17;
                    base_integral_19[g] = PAx*x10 - PGx[g]*x16;
                    base_integral_20[g] = PAy*x36 + tpi_3*x44 - PGy[g]*x40;
                    base_integral_21[g] = PAz*x36 - PGz[g]*x40;
                    base_integral_22[g] = PAy*x11 + tpi*x18 - PGy[g]*x17;
                    base_integral_23[g] = PAy*x10 - PGy[g]*x16;
                    base_integral_24[g] = PAz*x10 + tpi_3*x18 - PGz[g]*x16;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PAx*x0 - PGx[g]*x1;
                    const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PAx*x1 - PGx[g]*x3;
                    const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x6 = tpi*x5 + PAz*x0 - PGz[g]*x1;
                    const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x8 = tpi*x7 + PAz*x1 - PGz[g]*x3;
                    const T x9 = x0 - x1;
                    const T x10 = PAz*x6 + tpi_2*x9 - PGz[g]*x8;
                    const T x11 = PAy*x6 - PGy[g]*x8;
                    const T x12 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x13 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x14 = PAz*x3 + tpi*x12 - PGz[g]*x13;
                    const T x15 = x1 - x3;
                    const T x16 = PAz*x8 + tpi_2*x15 - PGz[g]*x14;
                    const T x17 = PAy*x8 - PGy[g]*x14;
                    const T x18 = x6 - x8;
                    const T x19 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x20 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x21 = PAx*x19 + tpi*x5 - PGx[g]*x20;
                    const T x22 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x23 = PAx*x20 + tpi*x7 - PGx[g]*x22;
                    const T x24 = x19 - x20;
                    const T x25 = PAx*x21 + tpi_2*x24 - PGx[g]*x23;
                    const T x26 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x27 = PAx*x22 + tpi*x12 - PGx[g]*x26;
                    const T x28 = x20 - x22;
                    const T x29 = PAx*x23 + tpi_2*x28 - PGx[g]*x27;
                    const T x30 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x31 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x32 = tpi*x5 + PAy*x30 - PGy[g]*x31;
                    const T x33 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x34 = tpi*x7 + PAy*x31 - PGy[g]*x33;
                    const T x35 = x30 - x31;
                    const T x36 = PAy*x32 + tpi_2*x35 - PGy[g]*x34;
                    const T x37 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x38 = PAy*x33 + tpi*x12 - PGy[g]*x37;
                    const T x39 = x31 - x33;
                    const T x40 = PAy*x34 + tpi_2*x39 - PGy[g]*x38;
                    const T x41 = x21 - x23;
                    const T x42 = PAz*x21 - PGz[g]*x23;
                    const T x43 = PAz*x23 - PGz[g]*x27;
                    const T x44 = x32 - x34;
                    const T x45 = PAx*x32 - PGx[g]*x34;
                    const T x46 = PAx*x34 - PGx[g]*x38;

                    base_integral_0[g] += x25;
                    base_integral_1[g] += PAy*x21 - PGy[g]*x23;
                    base_integral_2[g] += x42;
                    base_integral_3[g] += x45;
                    base_integral_4[g] += PAy*x2 - PGy[g]*x4;
                    base_integral_5[g] += PAx*x6 - PGx[g]*x8;
                    base_integral_6[g] += x36;
                    base_integral_7[g] += PAz*x32 - PGz[g]*x34;
                    base_integral_8[g] += x11;
                    base_integral_9[g] += x10;
                    base_integral_10[g] += PAx*x25 + tpi_3*x41 - PGx[g]*x29;
                    base_integral_11[g] += PAy*x25 - PGy[g]*x29;
                    base_integral_12[g] += PAz*x25 - PGz[g]*x29;
                    base_integral_13[g] += PAx*x45 + tpi*x44 - PGx[g]*x46;
                    base_integral_14[g] += PAy*x42 - PGy[g]*x43;
                    base_integral_15[g] += PAz*x42 + tpi*x41 - PGz[g]*x43;
                    base_integral_16[g] += PAx*x36 - PGx[g]*x40;
                    base_integral_17[g] += PAz*x45 - PGz[g]*x46;
                    base_integral_18[g] += PAx*x11 - PGx[g]*x17;
                    base_integral_19[g] += PAx*x10 - PGx[g]*x16;
                    base_integral_20[g] += PAy*x36 + tpi_3*x44 - PGy[g]*x40;
                    base_integral_21[g] += PAz*x36 - PGz[g]*x40;
                    base_integral_22[g] += PAy*x11 + tpi*x18 - PGy[g]*x17;
                    base_integral_23[g] += PAy*x10 - PGy[g]*x16;
                    base_integral_24[g] += PAz*x10 + tpi_3*x18 - PGz[g]*x16;

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


            integral_0[g] = base_integral_10[g] + ABx*base_integral_0[g];
            integral_1[g] = ABy*base_integral_0[g] + base_integral_11[g];
            integral_2[g] = ABz*base_integral_0[g] + base_integral_12[g];
            integral_3[g] = base_integral_11[g] + ABx*base_integral_1[g];
            integral_4[g] = ABy*base_integral_1[g] + base_integral_13[g];
            integral_5[g] = ABz*base_integral_1[g] + base_integral_14[g];
            integral_6[g] = base_integral_12[g] + ABx*base_integral_2[g];
            integral_7[g] = base_integral_14[g] + ABy*base_integral_2[g];
            integral_8[g] = ABz*base_integral_2[g] + base_integral_15[g];
            integral_9[g] = base_integral_13[g] + ABx*base_integral_3[g];
            integral_10[g] = ABy*base_integral_3[g] + base_integral_16[g];
            integral_11[g] = ABz*base_integral_3[g] + base_integral_17[g];
            integral_12[g] = base_integral_14[g] + ABx*base_integral_4[g];
            integral_13[g] = base_integral_17[g] + ABy*base_integral_4[g];
            integral_14[g] = ABz*base_integral_4[g] + base_integral_18[g];
            integral_15[g] = base_integral_15[g] + ABx*base_integral_5[g];
            integral_16[g] = base_integral_18[g] + ABy*base_integral_5[g];
            integral_17[g] = ABz*base_integral_5[g] + base_integral_19[g];
            integral_18[g] = base_integral_16[g] + ABx*base_integral_6[g];
            integral_19[g] = ABy*base_integral_6[g] + base_integral_20[g];
            integral_20[g] = ABz*base_integral_6[g] + base_integral_21[g];
            integral_21[g] = base_integral_17[g] + ABx*base_integral_7[g];
            integral_22[g] = base_integral_21[g] + ABy*base_integral_7[g];
            integral_23[g] = ABz*base_integral_7[g] + base_integral_22[g];
            integral_24[g] = base_integral_18[g] + ABx*base_integral_8[g];
            integral_25[g] = base_integral_22[g] + ABy*base_integral_8[g];
            integral_26[g] = ABz*base_integral_8[g] + base_integral_23[g];
            integral_27[g] = base_integral_19[g] + ABx*base_integral_9[g];
            integral_28[g] = base_integral_23[g] + ABy*base_integral_9[g];
            integral_29[g] = ABz*base_integral_9[g] + base_integral_24[g];

        } // for g

    }
} // HRR_F_P

template void HRR_F_P<float>(
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

template void HRR_F_P<double>(
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

