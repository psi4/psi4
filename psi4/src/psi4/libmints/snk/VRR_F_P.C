#include "VRR_F_P.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void VRR_F_P(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
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
    if (38*sizeof(T) >= 280) {

        if (nprim < 1) return;

        // The gridpoints are processed in blocks.  The integrals of a block are contracted in
        // a buffer that stays in the L1 cache and are then copied to the output arrays.
        enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
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
                const T PBx = v_PAx[prim] + ABx;
                const T PBy = v_PAy[prim] + ABy;
                const T PBz = v_PAz[prim] + ABz;
                const T tpi = 0.5 / v_p[prim];
                const T tpi_2 = 2*tpi;
                const T tpi_3 = 3*tpi;

                // overwrite integrals for first primitive, accumulate otherwise
                if (prim==0) {
                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x2 = PAy*x0 - PGy[g]*x1;
                        const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x4 = PAy*x1 - PGy[g]*x3;
                        const T x5 = PAx*x2 - PGx[g]*x4;
                        const T x6 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x7 = PAy*x3 - PGy[g]*x6;
                        const T x8 = PAx*x4 - PGx[g]*x7;
                        const T x9 = x2 - x4;
                        block_12[g] = PBx*x5 + tpi*x9 - PGx[g]*x8;
                        const T x10 = x0 - x1;
                        const T x11 = x1 - x3;
                        const T x12 = PAx*x10 - PGx[g]*x11;
                        block_13[g] = tpi*x12 + PBy*x5 - PGy[g]*x8;
                        const T x13 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x14 = tpi*x13 + PAz*x0 - PGz[g]*x1;
                        const T x15 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x16 = tpi*x15 + PAz*x1 - PGz[g]*x3;
                        const T x17 = tpi_2*x10 + PAz*x14 - PGz[g]*x16;
                        const T x18 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x19 = tpi*x18 + PAz*x3 - PGz[g]*x6;
                        const T x20 = tpi_2*x11 + PAz*x16 - PGz[g]*x19;
                        block_27[g] = PBx*x17 - PGx[g]*x20;
                        block_28[g] = PBy*x17 - PGy[g]*x20;
                        const T x21 = x14 - x16;
                        block_29[g] = tpi_3*x21 + PBz*x17 - PGz[g]*x20;
                        const T x22 = PAx*x14 - PGx[g]*x16;
                        const T x23 = PAx*x16 - PGx[g]*x19;
                        block_15[g] = PBx*x22 + tpi*x21 - PGx[g]*x23;
                        const T x24 = PAy*x14 - PGy[g]*x16;
                        const T x25 = PAy*x16 - PGy[g]*x19;
                        block_25[g] = tpi*x21 + PBy*x24 - PGy[g]*x25;
                        block_16[g] = PBy*x22 - PGy[g]*x23;
                        block_17[g] = tpi_2*x12 + PBz*x22 - PGz[g]*x23;
                        block_24[g] = PBx*x24 - PGx[g]*x25;
                        block_26[g] = tpi_2*x9 + PBz*x24 - PGz[g]*x25;
                        const T x26 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x27 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x28 = x26 - x27;
                        const T x29 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x30 = x27 - x29;
                        const T x31 = PAx*x28 - PGx[g]*x30;
                        block_14[g] = tpi*x31 + PBz*x5 - PGz[g]*x8;
                        const T x32 = tpi*x13 + PAy*x26 - PGy[g]*x27;
                        const T x33 = tpi*x15 + PAy*x27 - PGy[g]*x29;
                        const T x34 = tpi_2*x28 + PAy*x32 - PGy[g]*x33;
                        const T x35 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x36 = tpi*x18 + PAy*x29 - PGy[g]*x35;
                        const T x37 = tpi_2*x30 + PAy*x33 - PGy[g]*x36;
                        block_18[g] = PBx*x34 - PGx[g]*x37;
                        const T x38 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x39 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x40 = PAx*x38 + tpi*x13 - PGx[g]*x39;
                        const T x41 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x42 = PAx*x39 + tpi*x15 - PGx[g]*x41;
                        const T x43 = x38 - x39;
                        const T x44 = PAx*x40 + tpi_2*x43 - PGx[g]*x42;
                        const T x45 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x46 = PAx*x41 + tpi*x18 - PGx[g]*x45;
                        const T x47 = x39 - x41;
                        const T x48 = PAx*x42 + tpi_2*x47 - PGx[g]*x46;
                        block_1[g] = PBy*x44 - PGy[g]*x48;
                        block_2[g] = PBz*x44 - PGz[g]*x48;
                        const T x49 = x40 - x42;
                        block_0[g] = PBx*x44 + tpi_3*x49 - PGx[g]*x48;
                        block_20[g] = PBz*x34 - PGz[g]*x37;
                        const T x50 = x32 - x33;
                        block_19[g] = tpi_3*x50 + PBy*x34 - PGy[g]*x37;
                        const T x51 = PAz*x40 - PGz[g]*x42;
                        const T x52 = PAz*x42 - PGz[g]*x46;
                        block_6[g] = PBx*x51 + tpi_2*x12 - PGx[g]*x52;
                        const T x53 = PAy*x40 - PGy[g]*x42;
                        const T x54 = PAy*x42 - PGy[g]*x46;
                        block_3[g] = PBx*x53 + tpi_2*x31 - PGx[g]*x54;
                        block_4[g] = tpi*x49 + PBy*x53 - PGy[g]*x54;
                        block_5[g] = PBz*x53 - PGz[g]*x54;
                        block_8[g] = tpi*x49 + PBz*x51 - PGz[g]*x52;
                        block_7[g] = PBy*x51 - PGy[g]*x52;
                        const T x55 = PAx*x32 - PGx[g]*x33;
                        const T x56 = PAx*x33 - PGx[g]*x36;
                        block_10[g] = tpi_2*x31 + PBy*x55 - PGy[g]*x56;
                        const T x57 = PAz*x32 - PGz[g]*x33;
                        const T x58 = PAz*x33 - PGz[g]*x36;
                        block_22[g] = tpi_2*x9 + PBy*x57 - PGy[g]*x58;
                        block_9[g] = PBx*x55 + tpi*x50 - PGx[g]*x56;
                        block_11[g] = PBz*x55 - PGz[g]*x56;
                        block_23[g] = tpi*x50 + PBz*x57 - PGz[g]*x58;
                        block_21[g] = PBx*x57 - PGx[g]*x58;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x2 = PAy*x0 - PGy[g]*x1;
                        const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x4 = PAy*x1 - PGy[g]*x3;
                        const T x5 = PAx*x2 - PGx[g]*x4;
                        const T x6 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x7 = PAy*x3 - PGy[g]*x6;
                        const T x8 = PAx*x4 - PGx[g]*x7;
                        const T x9 = x2 - x4;
                        block_12[g] += PBx*x5 + tpi*x9 - PGx[g]*x8;
                        const T x10 = x0 - x1;
                        const T x11 = x1 - x3;
                        const T x12 = PAx*x10 - PGx[g]*x11;
                        block_13[g] += tpi*x12 + PBy*x5 - PGy[g]*x8;
                        const T x13 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x14 = tpi*x13 + PAz*x0 - PGz[g]*x1;
                        const T x15 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x16 = tpi*x15 + PAz*x1 - PGz[g]*x3;
                        const T x17 = tpi_2*x10 + PAz*x14 - PGz[g]*x16;
                        const T x18 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x19 = tpi*x18 + PAz*x3 - PGz[g]*x6;
                        const T x20 = tpi_2*x11 + PAz*x16 - PGz[g]*x19;
                        block_27[g] += PBx*x17 - PGx[g]*x20;
                        block_28[g] += PBy*x17 - PGy[g]*x20;
                        const T x21 = x14 - x16;
                        block_29[g] += tpi_3*x21 + PBz*x17 - PGz[g]*x20;
                        const T x22 = PAx*x14 - PGx[g]*x16;
                        const T x23 = PAx*x16 - PGx[g]*x19;
                        block_15[g] += PBx*x22 + tpi*x21 - PGx[g]*x23;
                        const T x24 = PAy*x14 - PGy[g]*x16;
                        const T x25 = PAy*x16 - PGy[g]*x19;
                        block_25[g] += tpi*x21 + PBy*x24 - PGy[g]*x25;
                        block_16[g] += PBy*x22 - PGy[g]*x23;
                        block_17[g] += tpi_2*x12 + PBz*x22 - PGz[g]*x23;
                        block_24[g] += PBx*x24 - PGx[g]*x25;
                        block_26[g] += tpi_2*x9 + PBz*x24 - PGz[g]*x25;
                        const T x26 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x27 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x28 = x26 - x27;
                        const T x29 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x30 = x27 - x29;
                        const T x31 = PAx*x28 - PGx[g]*x30;
                        block_14[g] += tpi*x31 + PBz*x5 - PGz[g]*x8;
                        const T x32 = tpi*x13 + PAy*x26 - PGy[g]*x27;
                        const T x33 = tpi*x15 + PAy*x27 - PGy[g]*x29;
                        const T x34 = tpi_2*x28 + PAy*x32 - PGy[g]*x33;
                        const T x35 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x36 = tpi*x18 + PAy*x29 - PGy[g]*x35;
                        const T x37 = tpi_2*x30 + PAy*x33 - PGy[g]*x36;
                        block_18[g] += PBx*x34 - PGx[g]*x37;
                        const T x38 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x39 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x40 = PAx*x38 + tpi*x13 - PGx[g]*x39;
                        const T x41 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x42 = PAx*x39 + tpi*x15 - PGx[g]*x41;
                        const T x43 = x38 - x39;
                        const T x44 = PAx*x40 + tpi_2*x43 - PGx[g]*x42;
                        const T x45 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x46 = PAx*x41 + tpi*x18 - PGx[g]*x45;
                        const T x47 = x39 - x41;
                        const T x48 = PAx*x42 + tpi_2*x47 - PGx[g]*x46;
                        block_1[g] += PBy*x44 - PGy[g]*x48;
                        block_2[g] += PBz*x44 - PGz[g]*x48;
                        const T x49 = x40 - x42;
                        block_0[g] += PBx*x44 + tpi_3*x49 - PGx[g]*x48;
                        block_20[g] += PBz*x34 - PGz[g]*x37;
                        const T x50 = x32 - x33;
                        block_19[g] += tpi_3*x50 + PBy*x34 - PGy[g]*x37;
                        const T x51 = PAz*x40 - PGz[g]*x42;
                        const T x52 = PAz*x42 - PGz[g]*x46;
                        block_6[g] += PBx*x51 + tpi_2*x12 - PGx[g]*x52;
                        const T x53 = PAy*x40 - PGy[g]*x42;
                        const T x54 = PAy*x42 - PGy[g]*x46;
                        block_3[g] += PBx*x53 + tpi_2*x31 - PGx[g]*x54;
                        block_4[g] += tpi*x49 + PBy*x53 - PGy[g]*x54;
                        block_5[g] += PBz*x53 - PGz[g]*x54;
                        block_8[g] += tpi*x49 + PBz*x51 - PGz[g]*x52;
                        block_7[g] += PBy*x51 - PGy[g]*x52;
                        const T x55 = PAx*x32 - PGx[g]*x33;
                        const T x56 = PAx*x33 - PGx[g]*x36;
                        block_10[g] += tpi_2*x31 + PBy*x55 - PGy[g]*x56;
                        const T x57 = PAz*x32 - PGz[g]*x33;
                        const T x58 = PAz*x33 - PGz[g]*x36;
                        block_22[g] += tpi_2*x9 + PBy*x57 - PGy[g]*x58;
                        block_9[g] += PBx*x55 + tpi*x50 - PGx[g]*x56;
                        block_11[g] += PBz*x55 - PGz[g]*x56;
                        block_23[g] += tpi*x50 + PBz*x57 - PGz[g]*x58;
                        block_21[g] += PBx*x57 - PGx[g]*x58;

                    } // for g
                } // if/else prim
            } // for prim

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
            const T PBx = v_PAx[prim] + ABx;
            const T PBy = v_PAy[prim] + ABy;
            const T PBz = v_PAz[prim] + ABz;
            const T tpi = 0.5 / v_p[prim];
            const T tpi_2 = 2*tpi;
            const T tpi_3 = 3*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PAy*x0 - PGy[g]*x1;
                    const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PAy*x1 - PGy[g]*x3;
                    const T x5 = PAx*x2 - PGx[g]*x4;
                    const T x6 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x7 = PAy*x3 - PGy[g]*x6;
                    const T x8 = PAx*x4 - PGx[g]*x7;
                    const T x9 = x2 - x4;
                    const T x10 = x0 - x1;
                    const T x11 = x1 - x3;
                    const T x12 = PAx*x10 - PGx[g]*x11;
                    const T x13 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x14 = tpi*x13 + PAz*x0 - PGz[g]*x1;
                    const T x15 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x16 = tpi*x15 + PAz*x1 - PGz[g]*x3;
                    const T x17 = tpi_2*x10 + PAz*x14 - PGz[g]*x16;
                    const T x18 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x19 = tpi*x18 + PAz*x3 - PGz[g]*x6;
                    const T x20 = tpi_2*x11 + PAz*x16 - PGz[g]*x19;
                    const T x21 = x14 - x16;
                    const T x22 = PAx*x14 - PGx[g]*x16;
                    const T x23 = PAx*x16 - PGx[g]*x19;
                    const T x24 = PAy*x14 - PGy[g]*x16;
                    const T x25 = PAy*x16 - PGy[g]*x19;
                    const T x26 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x27 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x28 = x26 - x27;
                    const T x29 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x30 = x27 - x29;
                    const T x31 = PAx*x28 - PGx[g]*x30;
                    const T x32 = tpi*x13 + PAy*x26 - PGy[g]*x27;
                    const T x33 = tpi*x15 + PAy*x27 - PGy[g]*x29;
                    const T x34 = tpi_2*x28 + PAy*x32 - PGy[g]*x33;
                    const T x35 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x36 = tpi*x18 + PAy*x29 - PGy[g]*x35;
                    const T x37 = tpi_2*x30 + PAy*x33 - PGy[g]*x36;
                    const T x38 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x39 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x40 = PAx*x38 + tpi*x13 - PGx[g]*x39;
                    const T x41 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x42 = PAx*x39 + tpi*x15 - PGx[g]*x41;
                    const T x43 = x38 - x39;
                    const T x44 = PAx*x40 + tpi_2*x43 - PGx[g]*x42;
                    const T x45 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x46 = PAx*x41 + tpi*x18 - PGx[g]*x45;
                    const T x47 = x39 - x41;
                    const T x48 = PAx*x42 + tpi_2*x47 - PGx[g]*x46;
                    const T x49 = x40 - x42;
                    const T x50 = x32 - x33;
                    const T x51 = PAz*x40 - PGz[g]*x42;
                    const T x52 = PAz*x42 - PGz[g]*x46;
                    const T x53 = PAy*x40 - PGy[g]*x42;
                    const T x54 = PAy*x42 - PGy[g]*x46;
                    const T x55 = PAx*x32 - PGx[g]*x33;
                    const T x56 = PAx*x33 - PGx[g]*x36;
                    const T x57 = PAz*x32 - PGz[g]*x33;
                    const T x58 = PAz*x33 - PGz[g]*x36;

                    integral_0[g] = PBx*x44 + tpi_3*x49 - PGx[g]*x48;
                    integral_1[g] = PBy*x44 - PGy[g]*x48;
                    integral_2[g] = PBz*x44 - PGz[g]*x48;
                    integral_3[g] = PBx*x53 + tpi_2*x31 - PGx[g]*x54;
                    integral_4[g] = tpi*x49 + PBy*x53 - PGy[g]*x54;
                    integral_5[g] = PBz*x53 - PGz[g]*x54;
                    integral_6[g] = PBx*x51 + tpi_2*x12 - PGx[g]*x52;
                    integral_7[g] = PBy*x51 - PGy[g]*x52;
                    integral_8[g] = tpi*x49 + PBz*x51 - PGz[g]*x52;
                    integral_9[g] = PBx*x55 + tpi*x50 - PGx[g]*x56;
                    integral_10[g] = tpi_2*x31 + PBy*x55 - PGy[g]*x56;
                    integral_11[g] = PBz*x55 - PGz[g]*x56;
                    integral_12[g] = PBx*x5 + tpi*x9 - PGx[g]*x8;
                    integral_13[g] = tpi*x12 + PBy*x5 - PGy[g]*x8;
                    integral_14[g] = tpi*x31 + PBz*x5 - PGz[g]*x8;
                    integral_15[g] = PBx*x22 + tpi*x21 - PGx[g]*x23;
                    integral_16[g] = PBy*x22 - PGy[g]*x23;
                    integral_17[g] = tpi_2*x12 + PBz*x22 - PGz[g]*x23;
                    integral_18[g] = PBx*x34 - PGx[g]*x37;
                    integral_19[g] = tpi_3*x50 + PBy*x34 - PGy[g]*x37;
                    integral_20[g] = PBz*x34 - PGz[g]*x37;
                    integral_21[g] = PBx*x57 - PGx[g]*x58;
                    integral_22[g] = tpi_2*x9 + PBy*x57 - PGy[g]*x58;
                    integral_23[g] = tpi*x50 + PBz*x57 - PGz[g]*x58;
                    integral_24[g] = PBx*x24 - PGx[g]*x25;
                    integral_25[g] = tpi*x21 + PBy*x24 - PGy[g]*x25;
                    integral_26[g] = tpi_2*x9 + PBz*x24 - PGz[g]*x25;
                    integral_27[g] = PBx*x17 - PGx[g]*x20;
                    integral_28[g] = PBy*x17 - PGy[g]*x20;
                    integral_29[g] = tpi_3*x21 + PBz*x17 - PGz[g]*x20;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PAy*x0 - PGy[g]*x1;
                    const T x3 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PAy*x1 - PGy[g]*x3;
                    const T x5 = PAx*x2 - PGx[g]*x4;
                    const T x6 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x7 = PAy*x3 - PGy[g]*x6;
                    const T x8 = PAx*x4 - PGx[g]*x7;
                    const T x9 = x2 - x4;
                    const T x10 = x0 - x1;
                    const T x11 = x1 - x3;
                    const T x12 = PAx*x10 - PGx[g]*x11;
                    const T x13 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x14 = tpi*x13 + PAz*x0 - PGz[g]*x1;
                    const T x15 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x16 = tpi*x15 + PAz*x1 - PGz[g]*x3;
                    const T x17 = tpi_2*x10 + PAz*x14 - PGz[g]*x16;
                    const T x18 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x19 = tpi*x18 + PAz*x3 - PGz[g]*x6;
                    const T x20 = tpi_2*x11 + PAz*x16 - PGz[g]*x19;
                    const T x21 = x14 - x16;
                    const T x22 = PAx*x14 - PGx[g]*x16;
                    const T x23 = PAx*x16 - PGx[g]*x19;
                    const T x24 = PAy*x14 - PGy[g]*x16;
                    const T x25 = PAy*x16 - PGy[g]*x19;
                    const T x26 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x27 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x28 = x26 - x27;
                    const T x29 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x30 = x27 - x29;
                    const T x31 = PAx*x28 - PGx[g]*x30;
                    const T x32 = tpi*x13 + PAy*x26 - PGy[g]*x27;
                    const T x33 = tpi*x15 + PAy*x27 - PGy[g]*x29;
                    const T x34 = tpi_2*x28 + PAy*x32 - PGy[g]*x33;
                    const T x35 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x36 = tpi*x18 + PAy*x29 - PGy[g]*x35;
                    const T x37 = tpi_2*x30 + PAy*x33 - PGy[g]*x36;
                    const T x38 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x39 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x40 = PAx*x38 + tpi*x13 - PGx[g]*x39;
                    const T x41 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x42 = PAx*x39 + tpi*x15 - PGx[g]*x41;
                    const T x43 = x38 - x39;
                    const T x44 = PAx*x40 + tpi_2*x43 - PGx[g]*x42;
                    const T x45 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x46 = PAx*x41 + tpi*x18 - PGx[g]*x45;
                    const T x47 = x39 - x41;
                    const T x48 = PAx*x42 + tpi_2*x47 - PGx[g]*x46;
                    const T x49 = x40 - x42;
                    const T x50 = x32 - x33;
                    const T x51 = PAz*x40 - PGz[g]*x42;
                    const T x52 = PAz*x42 - PGz[g]*x46;
                    const T x53 = PAy*x40 - PGy[g]*x42;
                    const T x54 = PAy*x42 - PGy[g]*x46;
                    const T x55 = PAx*x32 - PGx[g]*x33;
                    const T x56 = PAx*x33 - PGx[g]*x36;
                    const T x57 = PAz*x32 - PGz[g]*x33;
                    const T x58 = PAz*x33 - PGz[g]*x36;

                    integral_0[g] += PBx*x44 + tpi_3*x49 - PGx[g]*x48;
                    integral_1[g] += PBy*x44 - PGy[g]*x48;
                    integral_2[g] += PBz*x44 - PGz[g]*x48;
                    integral_3[g] += PBx*x53 + tpi_2*x31 - PGx[g]*x54;
                    integral_4[g] += tpi*x49 + PBy*x53 - PGy[g]*x54;
                    integral_5[g] += PBz*x53 - PGz[g]*x54;
                    integral_6[g] += PBx*x51 + tpi_2*x12 - PGx[g]*x52;
                    integral_7[g] += PBy*x51 - PGy[g]*x52;
                    integral_8[g] += tpi*x49 + PBz*x51 - PGz[g]*x52;
                    integral_9[g] += PBx*x55 + tpi*x50 - PGx[g]*x56;
                    integral_10[g] += tpi_2*x31 + PBy*x55 - PGy[g]*x56;
                    integral_11[g] += PBz*x55 - PGz[g]*x56;
                    integral_12[g] += PBx*x5 + tpi*x9 - PGx[g]*x8;
                    integral_13[g] += tpi*x12 + PBy*x5 - PGy[g]*x8;
                    integral_14[g] += tpi*x31 + PBz*x5 - PGz[g]*x8;
                    integral_15[g] += PBx*x22 + tpi*x21 - PGx[g]*x23;
                    integral_16[g] += PBy*x22 - PGy[g]*x23;
                    integral_17[g] += tpi_2*x12 + PBz*x22 - PGz[g]*x23;
                    integral_18[g] += PBx*x34 - PGx[g]*x37;
                    integral_19[g] += tpi_3*x50 + PBy*x34 - PGy[g]*x37;
                    integral_20[g] += PBz*x34 - PGz[g]*x37;
                    integral_21[g] += PBx*x57 - PGx[g]*x58;
                    integral_22[g] += tpi_2*x9 + PBy*x57 - PGy[g]*x58;
                    integral_23[g] += tpi*x50 + PBz*x57 - PGz[g]*x58;
                    integral_24[g] += PBx*x24 - PGx[g]*x25;
                    integral_25[g] += tpi*x21 + PBy*x24 - PGy[g]*x25;
                    integral_26[g] += tpi_2*x9 + PBz*x24 - PGz[g]*x25;
                    integral_27[g] += PBx*x17 - PGx[g]*x20;
                    integral_28[g] += PBy*x17 - PGy[g]*x20;
                    integral_29[g] += tpi_3*x21 + PBz*x17 - PGz[g]*x20;

                } // for g
            } // if/else prim
        } // for prim

    }
} // VRR_F_P

template void VRR_F_P<float>(
        float* integrals, 
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

template void VRR_F_P<double>(
        double* integrals, 
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

