#include "HRR_D_D.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_D_D(
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
    if (nprim < 1) return;

    // The gridpoints are processed in blocks.  The integrals of a block are contracted in
    // a buffer that stays in the L1 cache and are then copied to the output arrays.
    enum { g_block = (sizeof(T) <= 4) ? 32 : 16 };
    alignas(64) T base_integral[31*g_block];
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
    T* __restrict base_integral_25 = &base_integral[25*g_block];
    T* __restrict base_integral_26 = &base_integral[26*g_block];
    T* __restrict base_integral_27 = &base_integral[27*g_block];
    T* __restrict base_integral_28 = &base_integral[28*g_block];
    T* __restrict base_integral_29 = &base_integral[29*g_block];
    T* __restrict base_integral_30 = &base_integral[30*g_block];
    alignas(64) T block[36*g_block];
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
    T* __restrict block_30 = &block[30*g_block];
    T* __restrict block_31 = &block[31*g_block];
    T* __restrict block_32 = &block[32*g_block];
    T* __restrict block_33 = &block[33*g_block];
    T* __restrict block_34 = &block[34*g_block];
    T* __restrict block_35 = &block[35*g_block];

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

                    const T x0 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    base_integral_1[g] = PAx*x0 - PGx[g]*x1;
                    const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x3 = tpi*x2 + PAy*x0 - PGy[g]*x1;
                    base_integral_3[g] = x3;
                    const T x4 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x5 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    base_integral_4[g] = PAy*x4 - PGy[g]*x5;
                    const T x6 = PAx*x4 - PGx[g]*x5;
                    base_integral_2[g] = x6;
                    const T x7 = tpi*x2 + PAz*x4 - PGz[g]*x5;
                    base_integral_5[g] = x7;
                    const T x8 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x9 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x10 = PAx*x8 + tpi*x2 - PGx[g]*x9;
                    base_integral_0[g] = x10;
                    const T x11 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x12 = PAx*x5 - PGx[g]*x11;
                    base_integral_10[g] = PAy*x6 - PGy[g]*x12;
                    const T x13 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x14 = PAz*x5 + tpi*x13 - PGz[g]*x11;
                    base_integral_11[g] = PAx*x7 - PGx[g]*x14;
                    const T x15 = x4 - x5;
                    const T x16 = PAz*x7 + tpi_2*x15 - PGz[g]*x14;
                    base_integral_15[g] = x16;
                    const T x17 = PAy*x7 - PGy[g]*x14;
                    base_integral_14[g] = x17;
                    const T x18 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x19 = PAx*x9 + tpi*x13 - PGx[g]*x18;
                    base_integral_7[g] = PAy*x10 - PGy[g]*x19;
                    const T x20 = x8 - x9;
                    const T x21 = PAx*x10 + tpi_2*x20 - PGx[g]*x19;
                    base_integral_6[g] = x21;
                    const T x22 = PAz*x10 - PGz[g]*x19;
                    base_integral_8[g] = x22;
                    const T x23 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x24 = PAy*x1 + tpi*x13 - PGy[g]*x23;
                    base_integral_13[g] = PAz*x3 - PGz[g]*x24;
                    const T x25 = x0 - x1;
                    const T x26 = PAy*x3 + tpi_2*x25 - PGy[g]*x24;
                    base_integral_12[g] = x26;
                    const T x27 = PAx*x3 - PGx[g]*x24;
                    base_integral_9[g] = x27;
                    const T x28 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x29 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x30 = PAx*x18 + tpi*x29 - PGx[g]*x28;
                    const T x31 = x9 - x18;
                    const T x32 = PAx*x19 + tpi_2*x31 - PGx[g]*x30;
                    base_integral_17[g] = PAy*x21 - PGy[g]*x32;
                    const T x33 = x10 - x19;
                    const T x34 = PAz*x19 - PGz[g]*x30;
                    base_integral_21[g] = PAz*x22 + tpi*x33 - PGz[g]*x34;
                    base_integral_20[g] = PAy*x22 - PGy[g]*x34;
                    base_integral_16[g] = PAx*x21 + tpi_3*x33 - PGx[g]*x32;
                    base_integral_18[g] = PAz*x21 - PGz[g]*x32;
                    const T x35 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x36 = PAy*x23 + tpi*x29 - PGy[g]*x35;
                    const T x37 = x1 - x23;
                    const T x38 = PAy*x24 + tpi_2*x37 - PGy[g]*x36;
                    base_integral_22[g] = PAx*x26 - PGx[g]*x38;
                    const T x39 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x40 = PAz*x11 + tpi*x29 - PGz[g]*x39;
                    const T x41 = x5 - x11;
                    const T x42 = PAz*x14 + tpi_2*x41 - PGz[g]*x40;
                    base_integral_25[g] = PAx*x16 - PGx[g]*x42;
                    const T x43 = PAx*x24 - PGx[g]*x36;
                    const T x44 = x3 - x24;
                    base_integral_19[g] = PAx*x27 + tpi*x44 - PGx[g]*x43;
                    base_integral_23[g] = PAz*x27 - PGz[g]*x43;
                    base_integral_26[g] = PAy*x26 + tpi_3*x44 - PGy[g]*x38;
                    base_integral_27[g] = PAz*x26 - PGz[g]*x38;
                    const T x45 = PAy*x14 - PGy[g]*x40;
                    const T x46 = x7 - x14;
                    base_integral_28[g] = PAy*x17 + tpi*x46 - PGy[g]*x45;
                    base_integral_24[g] = PAx*x17 - PGx[g]*x45;
                    base_integral_30[g] = PAz*x16 + tpi_3*x46 - PGz[g]*x42;
                    base_integral_29[g] = PAy*x16 - PGy[g]*x42;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    base_integral_1[g] += PAx*x0 - PGx[g]*x1;
                    const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x3 = tpi*x2 + PAy*x0 - PGy[g]*x1;
                    base_integral_3[g] += x3;
                    const T x4 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x5 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    base_integral_4[g] += PAy*x4 - PGy[g]*x5;
                    const T x6 = PAx*x4 - PGx[g]*x5;
                    base_integral_2[g] += x6;
                    const T x7 = tpi*x2 + PAz*x4 - PGz[g]*x5;
                    base_integral_5[g] += x7;
                    const T x8 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x9 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x10 = PAx*x8 + tpi*x2 - PGx[g]*x9;
                    base_integral_0[g] += x10;
                    const T x11 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x12 = PAx*x5 - PGx[g]*x11;
                    base_integral_10[g] += PAy*x6 - PGy[g]*x12;
                    const T x13 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x14 = PAz*x5 + tpi*x13 - PGz[g]*x11;
                    base_integral_11[g] += PAx*x7 - PGx[g]*x14;
                    const T x15 = x4 - x5;
                    const T x16 = PAz*x7 + tpi_2*x15 - PGz[g]*x14;
                    base_integral_15[g] += x16;
                    const T x17 = PAy*x7 - PGy[g]*x14;
                    base_integral_14[g] += x17;
                    const T x18 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x19 = PAx*x9 + tpi*x13 - PGx[g]*x18;
                    base_integral_7[g] += PAy*x10 - PGy[g]*x19;
                    const T x20 = x8 - x9;
                    const T x21 = PAx*x10 + tpi_2*x20 - PGx[g]*x19;
                    base_integral_6[g] += x21;
                    const T x22 = PAz*x10 - PGz[g]*x19;
                    base_integral_8[g] += x22;
                    const T x23 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x24 = PAy*x1 + tpi*x13 - PGy[g]*x23;
                    base_integral_13[g] += PAz*x3 - PGz[g]*x24;
                    const T x25 = x0 - x1;
                    const T x26 = PAy*x3 + tpi_2*x25 - PGy[g]*x24;
                    base_integral_12[g] += x26;
                    const T x27 = PAx*x3 - PGx[g]*x24;
                    base_integral_9[g] += x27;
                    const T x28 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x29 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x30 = PAx*x18 + tpi*x29 - PGx[g]*x28;
                    const T x31 = x9 - x18;
                    const T x32 = PAx*x19 + tpi_2*x31 - PGx[g]*x30;
                    base_integral_17[g] += PAy*x21 - PGy[g]*x32;
                    const T x33 = x10 - x19;
                    const T x34 = PAz*x19 - PGz[g]*x30;
                    base_integral_21[g] += PAz*x22 + tpi*x33 - PGz[g]*x34;
                    base_integral_20[g] += PAy*x22 - PGy[g]*x34;
                    base_integral_16[g] += PAx*x21 + tpi_3*x33 - PGx[g]*x32;
                    base_integral_18[g] += PAz*x21 - PGz[g]*x32;
                    const T x35 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x36 = PAy*x23 + tpi*x29 - PGy[g]*x35;
                    const T x37 = x1 - x23;
                    const T x38 = PAy*x24 + tpi_2*x37 - PGy[g]*x36;
                    base_integral_22[g] += PAx*x26 - PGx[g]*x38;
                    const T x39 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x40 = PAz*x11 + tpi*x29 - PGz[g]*x39;
                    const T x41 = x5 - x11;
                    const T x42 = PAz*x14 + tpi_2*x41 - PGz[g]*x40;
                    base_integral_25[g] += PAx*x16 - PGx[g]*x42;
                    const T x43 = PAx*x24 - PGx[g]*x36;
                    const T x44 = x3 - x24;
                    base_integral_19[g] += PAx*x27 + tpi*x44 - PGx[g]*x43;
                    base_integral_23[g] += PAz*x27 - PGz[g]*x43;
                    base_integral_26[g] += PAy*x26 + tpi_3*x44 - PGy[g]*x38;
                    base_integral_27[g] += PAz*x26 - PGz[g]*x38;
                    const T x45 = PAy*x14 - PGy[g]*x40;
                    const T x46 = x7 - x14;
                    base_integral_28[g] += PAy*x17 + tpi*x46 - PGy[g]*x45;
                    base_integral_24[g] += PAx*x17 - PGx[g]*x45;
                    base_integral_30[g] += PAz*x16 + tpi_3*x46 - PGz[g]*x42;
                    base_integral_29[g] += PAy*x16 - PGy[g]*x42;

                } // for g
            } // if/else prim
        } // for prim

        #pragma omp simd
        for (int g = 0; g < g_block; ++g) {

            const T x0 = base_integral_19[g] + ABx*base_integral_9[g];
            const T x1 = base_integral_9[g] + ABx*base_integral_3[g];
            block_18[g] = x0 + ABx*x1;
            const T x2 = ABy*base_integral_3[g] + base_integral_12[g];
            const T x3 = ABy*base_integral_12[g] + base_integral_26[g];
            block_21[g] = ABy*x2 + x3;
            const T x4 = ABy*base_integral_9[g] + base_integral_22[g];
            block_19[g] = x4 + ABx*x2;
            const T x5 = ABy*base_integral_1[g] + base_integral_9[g];
            block_9[g] = ABy*x5 + x4;
            const T x6 = ABy*base_integral_7[g] + base_integral_19[g];
            block_7[g] = x6 + ABx*x5;
            const T x7 = ABy*base_integral_0[g] + base_integral_7[g];
            block_3[g] = x6 + ABy*x7;
            const T x8 = base_integral_7[g] + ABx*base_integral_1[g];
            const T x9 = base_integral_17[g] + ABx*base_integral_7[g];
            block_6[g] = x9 + ABx*x8;
            const T x10 = base_integral_6[g] + ABx*base_integral_0[g];
            block_1[g] = ABy*x10 + x9;
            const T x11 = base_integral_16[g] + ABx*base_integral_6[g];
            block_0[g] = x11 + ABx*x10;
            const T x12 = base_integral_18[g] + ABx*base_integral_8[g];
            block_2[g] = ABz*x10 + x12;
            const T x13 = ABz*base_integral_0[g] + base_integral_8[g];
            const T x14 = ABz*base_integral_7[g] + base_integral_20[g];
            block_4[g] = x14 + ABy*x13;
            const T x15 = ABz*base_integral_1[g] + base_integral_10[g];
            block_8[g] = x14 + ABx*x15;
            const T x16 = ABz*base_integral_8[g] + base_integral_21[g];
            block_5[g] = ABz*x13 + x16;
            const T x17 = base_integral_8[g] + ABx*base_integral_2[g];
            block_12[g] = x12 + ABx*x17;
            const T x18 = base_integral_20[g] + ABx*base_integral_10[g];
            block_13[g] = ABy*x17 + x18;
            const T x19 = base_integral_10[g] + ABx*base_integral_4[g];
            block_24[g] = x18 + ABx*x19;
            const T x20 = ABz*base_integral_9[g] + base_integral_23[g];
            block_10[g] = ABy*x15 + x20;
            const T x21 = ABy*base_integral_10[g] + base_integral_23[g];
            const T x22 = base_integral_10[g] + ABy*base_integral_2[g];
            block_15[g] = x21 + ABy*x22;
            const T x23 = ABz*base_integral_10[g] + base_integral_24[g];
            block_11[g] = ABz*x15 + x23;
            const T x24 = ABz*base_integral_2[g] + base_integral_11[g];
            block_14[g] = x16 + ABx*x24;
            block_16[g] = x23 + ABy*x24;
            const T x25 = ABz*base_integral_3[g] + base_integral_13[g];
            block_20[g] = x20 + ABx*x25;
            const T x26 = base_integral_13[g] + ABy*base_integral_4[g];
            block_25[g] = x21 + ABx*x26;
            const T x27 = ABy*base_integral_13[g] + base_integral_27[g];
            block_22[g] = ABz*x2 + x27;
            block_27[g] = x27 + ABy*x26;
            const T x28 = ABz*base_integral_13[g] + base_integral_28[g];
            block_23[g] = ABz*x25 + x28;
            const T x29 = ABz*base_integral_4[g] + base_integral_14[g];
            block_26[g] = x23 + ABx*x29;
            block_28[g] = x28 + ABy*x29;
            const T x30 = base_integral_21[g] + ABx*base_integral_11[g];
            const T x31 = base_integral_11[g] + ABx*base_integral_5[g];
            block_30[g] = x30 + ABx*x31;
            const T x32 = base_integral_24[g] + ABy*base_integral_11[g];
            const T x33 = base_integral_14[g] + ABy*base_integral_5[g];
            block_31[g] = x32 + ABx*x33;
            const T x34 = base_integral_28[g] + ABy*base_integral_14[g];
            block_33[g] = ABy*x33 + x34;
            const T x35 = ABz*base_integral_11[g] + base_integral_25[g];
            block_17[g] = ABz*x24 + x35;
            const T x36 = ABz*base_integral_14[g] + base_integral_29[g];
            block_29[g] = ABz*x29 + x36;
            const T x37 = ABz*base_integral_5[g] + base_integral_15[g];
            block_32[g] = x35 + ABx*x37;
            block_34[g] = x36 + ABy*x37;
            const T x38 = ABz*base_integral_15[g] + base_integral_30[g];
            block_35[g] = ABz*x37 + x38;

        } // for g

        // copy the block to the output arrays
        for (int i = 0; i < 36; ++i) {
            const T* __restrict src = &block[i*g_block];
            T* __restrict dst = &integrals[i*gb_size + g_lo];
            if (padded) {
                SNK_COPY(dst, src, gb_size);
            } else {
                SNK_COPY(dst, src, g_block);
            }
        }
    } // for g_lo
} // HRR_D_D

template void HRR_D_D<float>(
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

template void HRR_D_D<double>(
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

