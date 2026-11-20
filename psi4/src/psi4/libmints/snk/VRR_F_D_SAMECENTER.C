#include "VRR_F_D_SAMECENTER.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void VRR_F_D_SAMECENTER(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)
        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively
        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)
        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair
        int gb_size)               //!< number of gridpoints in this grid batch 
{
    if (69*sizeof(T) >= 280) {

        if (nprim < 1) return;

        // The gridpoints are processed in blocks.  The integrals of a block are contracted in
        // a buffer that stays in the L1 cache and are then copied to the output arrays.
        enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
        alignas(64) T block[21*g_block];
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
                        block_6[g] = -PGx[g]*x13 - x9;
                        block_11[g] = -PGz[g]*x13;
                        block_15[g] = -4*x9 - PGy[g]*x14;
                        block_10[g] = -PGx[g]*x14;
                        block_16[g] = -PGz[g]*x14;
                        const T x15 = PGx[g]*x5;
                        const T x16 = PGx[g]*x15 + x3;
                        const T x17 = PGx[g]*scaled_boys_5[g];
                        const T x18 = PGx[g]*x17 + x5;
                        const T x19 = PGy[g]*x18;
                        const T x20 = PGx[g]*x18 + 2*x15;
                        const T x21 = x16 + PGy[g]*x19;
                        block_7[g] = -PGz[g]*x21;
                        const T x22 = PGx[g]*x20 + 3*x16;
                        const T x23 = 2*x3 + x16;
                        const T x24 = PGx[g]*x23;
                        block_0[g] = -PGx[g]*x22 - 4*x24;
                        const T x25 = PGy[g]*x20;
                        block_3[g] = -x24 - PGy[g]*x25;
                        const T x26 = PGz[g]*x20;
                        block_5[g] = -x24 - PGz[g]*x26;
                        block_4[g] = -PGy[g]*x26;
                        block_1[g] = -PGy[g]*x22;
                        block_2[g] = -PGz[g]*x22;
                        const T x27 = PGz[g]*x5;
                        const T x28 = x3 + PGz[g]*x27;
                        const T x29 = 2*x3 + x28;
                        const T x30 = PGz[g]*x29;
                        const T x31 = PGz[g]*scaled_boys_5[g];
                        const T x32 = x5 + PGz[g]*x31;
                        const T x33 = PGz[g]*x32 + 2*x27;
                        const T x34 = PGx[g]*x33;
                        block_9[g] = -PGx[g]*x34 - x30;
                        const T x35 = 3*x28 + PGz[g]*x33;
                        block_20[g] = -4*x30 - PGz[g]*x35;
                        const T x36 = PGy[g]*x33;
                        block_18[g] = -x30 - PGy[g]*x36;
                        block_13[g] = -PGx[g]*x36;
                        block_14[g] = -PGx[g]*x35;
                        block_19[g] = -PGy[g]*x35;
                        const T x37 = PGy[g]*x32;
                        const T x38 = PGy[g]*x37 + x28;
                        block_12[g] = -PGx[g]*x38;
                        const T x39 = PGy[g]*x28;
                        block_17[g] = -2*x39 - PGy[g]*x38;
                        const T x40 = PGx[g]*x37;
                        block_8[g] = -PGx[g]*x40 - x39;

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
                        block_6[g] += -PGx[g]*x13 - x9;
                        block_11[g] += -PGz[g]*x13;
                        block_15[g] += -4*x9 - PGy[g]*x14;
                        block_10[g] += -PGx[g]*x14;
                        block_16[g] += -PGz[g]*x14;
                        const T x15 = PGx[g]*x5;
                        const T x16 = PGx[g]*x15 + x3;
                        const T x17 = PGx[g]*scaled_boys_5[g];
                        const T x18 = PGx[g]*x17 + x5;
                        const T x19 = PGy[g]*x18;
                        const T x20 = PGx[g]*x18 + 2*x15;
                        const T x21 = x16 + PGy[g]*x19;
                        block_7[g] += -PGz[g]*x21;
                        const T x22 = PGx[g]*x20 + 3*x16;
                        const T x23 = 2*x3 + x16;
                        const T x24 = PGx[g]*x23;
                        block_0[g] += -PGx[g]*x22 - 4*x24;
                        const T x25 = PGy[g]*x20;
                        block_3[g] += -x24 - PGy[g]*x25;
                        const T x26 = PGz[g]*x20;
                        block_5[g] += -x24 - PGz[g]*x26;
                        block_4[g] += -PGy[g]*x26;
                        block_1[g] += -PGy[g]*x22;
                        block_2[g] += -PGz[g]*x22;
                        const T x27 = PGz[g]*x5;
                        const T x28 = x3 + PGz[g]*x27;
                        const T x29 = 2*x3 + x28;
                        const T x30 = PGz[g]*x29;
                        const T x31 = PGz[g]*scaled_boys_5[g];
                        const T x32 = x5 + PGz[g]*x31;
                        const T x33 = PGz[g]*x32 + 2*x27;
                        const T x34 = PGx[g]*x33;
                        block_9[g] += -PGx[g]*x34 - x30;
                        const T x35 = 3*x28 + PGz[g]*x33;
                        block_20[g] += -4*x30 - PGz[g]*x35;
                        const T x36 = PGy[g]*x33;
                        block_18[g] += -x30 - PGy[g]*x36;
                        block_13[g] += -PGx[g]*x36;
                        block_14[g] += -PGx[g]*x35;
                        block_19[g] += -PGy[g]*x35;
                        const T x37 = PGy[g]*x32;
                        const T x38 = PGy[g]*x37 + x28;
                        block_12[g] += -PGx[g]*x38;
                        const T x39 = PGy[g]*x28;
                        block_17[g] += -2*x39 - PGy[g]*x38;
                        const T x40 = PGx[g]*x37;
                        block_8[g] += -PGx[g]*x40 - x39;

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
                const T* __restrict src = &block[row[i]*g_block];
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
        T* __restrict integral_30 = &integrals[30*gb_size];
        T* __restrict integral_31 = &integrals[31*gb_size];
        T* __restrict integral_32 = &integrals[32*gb_size];
        T* __restrict integral_33 = &integrals[33*gb_size];
        T* __restrict integral_34 = &integrals[34*gb_size];
        T* __restrict integral_35 = &integrals[35*gb_size];
        T* __restrict integral_36 = &integrals[36*gb_size];
        T* __restrict integral_37 = &integrals[37*gb_size];
        T* __restrict integral_38 = &integrals[38*gb_size];
        T* __restrict integral_39 = &integrals[39*gb_size];
        T* __restrict integral_40 = &integrals[40*gb_size];
        T* __restrict integral_41 = &integrals[41*gb_size];
        T* __restrict integral_42 = &integrals[42*gb_size];
        T* __restrict integral_43 = &integrals[43*gb_size];
        T* __restrict integral_44 = &integrals[44*gb_size];
        T* __restrict integral_45 = &integrals[45*gb_size];
        T* __restrict integral_46 = &integrals[46*gb_size];
        T* __restrict integral_47 = &integrals[47*gb_size];
        T* __restrict integral_48 = &integrals[48*gb_size];
        T* __restrict integral_49 = &integrals[49*gb_size];
        T* __restrict integral_50 = &integrals[50*gb_size];
        T* __restrict integral_51 = &integrals[51*gb_size];
        T* __restrict integral_52 = &integrals[52*gb_size];
        T* __restrict integral_53 = &integrals[53*gb_size];
        T* __restrict integral_54 = &integrals[54*gb_size];
        T* __restrict integral_55 = &integrals[55*gb_size];
        T* __restrict integral_56 = &integrals[56*gb_size];
        T* __restrict integral_57 = &integrals[57*gb_size];
        T* __restrict integral_58 = &integrals[58*gb_size];
        T* __restrict integral_59 = &integrals[59*gb_size];

        for (int prim = 0; prim < nprim; ++prim) {

            const T* __restrict PGx = &v_PGx[prim*gb_size];
            const T* __restrict PGy = &v_PGy[prim*gb_size];
            const T* __restrict PGz = &v_PGz[prim*gb_size];
            const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*6*gb_size+0*gb_size];
            const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*6*gb_size+1*gb_size];
            const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*6*gb_size+2*gb_size];
            const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*6*gb_size+3*gb_size];
            const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*6*gb_size+4*gb_size];
            const T* __restrict scaled_boys_5 = &v_scaled_boys[prim*6*gb_size+5*gb_size];

            const T tpi = 0.5 / v_p[prim];
            const T tpi_p2 = tpi*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

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
                    const T x15 = -PGx[g]*x13 - x9;
                    const T x16 = -PGz[g]*x13;
                    const T x17 = -PGx[g]*x14;
                    const T x18 = -PGz[g]*x14;
                    const T x19 = PGx[g]*x5;
                    const T x20 = PGx[g]*x19 + x3;
                    const T x21 = PGx[g]*scaled_boys_5[g];
                    const T x22 = PGx[g]*x21 + x5;
                    const T x23 = PGy[g]*x22;
                    const T x24 = PGx[g]*x22 + 2*x19;
                    const T x25 = x20 + PGy[g]*x23;
                    const T x26 = -PGz[g]*x25;
                    const T x27 = PGx[g]*x24 + 3*x20;
                    const T x28 = 2*x3 + x20;
                    const T x29 = PGx[g]*x28;
                    const T x30 = PGy[g]*x24;
                    const T x31 = -x29 - PGy[g]*x30;
                    const T x32 = PGz[g]*x24;
                    const T x33 = -x29 - PGz[g]*x32;
                    const T x34 = -PGy[g]*x32;
                    const T x35 = -PGy[g]*x27;
                    const T x36 = -PGz[g]*x27;
                    const T x37 = PGz[g]*x5;
                    const T x38 = x3 + PGz[g]*x37;
                    const T x39 = 2*x3 + x38;
                    const T x40 = PGz[g]*x39;
                    const T x41 = PGz[g]*scaled_boys_5[g];
                    const T x42 = x5 + PGz[g]*x41;
                    const T x43 = PGz[g]*x42 + 2*x37;
                    const T x44 = PGx[g]*x43;
                    const T x45 = -PGx[g]*x44 - x40;
                    const T x46 = 3*x38 + PGz[g]*x43;
                    const T x47 = PGy[g]*x43;
                    const T x48 = -x40 - PGy[g]*x47;
                    const T x49 = -PGx[g]*x47;
                    const T x50 = -PGx[g]*x46;
                    const T x51 = -PGy[g]*x46;
                    const T x52 = PGy[g]*x42;
                    const T x53 = PGy[g]*x52 + x38;
                    const T x54 = -PGx[g]*x53;
                    const T x55 = PGy[g]*x38;
                    const T x56 = -2*x55 - PGy[g]*x53;
                    const T x57 = PGx[g]*x52;
                    const T x58 = -PGx[g]*x57 - x55;

                    integral_0[g] = -PGx[g]*x27 - 4*x29;
                    integral_1[g] = x35;
                    integral_2[g] = x36;
                    integral_3[g] = x31;
                    integral_4[g] = x34;
                    integral_5[g] = x33;
                    integral_6[g] = x35;
                    integral_7[g] = x31;
                    integral_8[g] = x34;
                    integral_9[g] = x15;
                    integral_10[g] = x26;
                    integral_11[g] = x58;
                    integral_12[g] = x36;
                    integral_13[g] = x34;
                    integral_14[g] = x33;
                    integral_15[g] = x26;
                    integral_16[g] = x58;
                    integral_17[g] = x45;
                    integral_18[g] = x31;
                    integral_19[g] = x15;
                    integral_20[g] = x26;
                    integral_21[g] = x17;
                    integral_22[g] = x16;
                    integral_23[g] = x54;
                    integral_24[g] = x34;
                    integral_25[g] = x26;
                    integral_26[g] = x58;
                    integral_27[g] = x16;
                    integral_28[g] = x54;
                    integral_29[g] = x49;
                    integral_30[g] = x33;
                    integral_31[g] = x58;
                    integral_32[g] = x45;
                    integral_33[g] = x54;
                    integral_34[g] = x49;
                    integral_35[g] = x50;
                    integral_36[g] = x15;
                    integral_37[g] = x17;
                    integral_38[g] = x16;
                    integral_39[g] = -4*x9 - PGy[g]*x14;
                    integral_40[g] = x18;
                    integral_41[g] = x56;
                    integral_42[g] = x26;
                    integral_43[g] = x16;
                    integral_44[g] = x54;
                    integral_45[g] = x18;
                    integral_46[g] = x56;
                    integral_47[g] = x48;
                    integral_48[g] = x58;
                    integral_49[g] = x54;
                    integral_50[g] = x49;
                    integral_51[g] = x56;
                    integral_52[g] = x48;
                    integral_53[g] = x51;
                    integral_54[g] = x45;
                    integral_55[g] = x49;
                    integral_56[g] = x50;
                    integral_57[g] = x48;
                    integral_58[g] = x51;
                    integral_59[g] = -4*x40 - PGz[g]*x46;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

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
                    const T x15 = -PGx[g]*x13 - x9;
                    const T x16 = -PGz[g]*x13;
                    const T x17 = -PGx[g]*x14;
                    const T x18 = -PGz[g]*x14;
                    const T x19 = PGx[g]*x5;
                    const T x20 = PGx[g]*x19 + x3;
                    const T x21 = PGx[g]*scaled_boys_5[g];
                    const T x22 = PGx[g]*x21 + x5;
                    const T x23 = PGy[g]*x22;
                    const T x24 = PGx[g]*x22 + 2*x19;
                    const T x25 = x20 + PGy[g]*x23;
                    const T x26 = -PGz[g]*x25;
                    const T x27 = PGx[g]*x24 + 3*x20;
                    const T x28 = 2*x3 + x20;
                    const T x29 = PGx[g]*x28;
                    const T x30 = PGy[g]*x24;
                    const T x31 = -x29 - PGy[g]*x30;
                    const T x32 = PGz[g]*x24;
                    const T x33 = -x29 - PGz[g]*x32;
                    const T x34 = -PGy[g]*x32;
                    const T x35 = -PGy[g]*x27;
                    const T x36 = -PGz[g]*x27;
                    const T x37 = PGz[g]*x5;
                    const T x38 = x3 + PGz[g]*x37;
                    const T x39 = 2*x3 + x38;
                    const T x40 = PGz[g]*x39;
                    const T x41 = PGz[g]*scaled_boys_5[g];
                    const T x42 = x5 + PGz[g]*x41;
                    const T x43 = PGz[g]*x42 + 2*x37;
                    const T x44 = PGx[g]*x43;
                    const T x45 = -PGx[g]*x44 - x40;
                    const T x46 = 3*x38 + PGz[g]*x43;
                    const T x47 = PGy[g]*x43;
                    const T x48 = -x40 - PGy[g]*x47;
                    const T x49 = -PGx[g]*x47;
                    const T x50 = -PGx[g]*x46;
                    const T x51 = -PGy[g]*x46;
                    const T x52 = PGy[g]*x42;
                    const T x53 = PGy[g]*x52 + x38;
                    const T x54 = -PGx[g]*x53;
                    const T x55 = PGy[g]*x38;
                    const T x56 = -2*x55 - PGy[g]*x53;
                    const T x57 = PGx[g]*x52;
                    const T x58 = -PGx[g]*x57 - x55;

                    integral_0[g] += -PGx[g]*x27 - 4*x29;
                    integral_1[g] += x35;
                    integral_2[g] += x36;
                    integral_3[g] += x31;
                    integral_4[g] += x34;
                    integral_5[g] += x33;
                    integral_6[g] += x35;
                    integral_7[g] += x31;
                    integral_8[g] += x34;
                    integral_9[g] += x15;
                    integral_10[g] += x26;
                    integral_11[g] += x58;
                    integral_12[g] += x36;
                    integral_13[g] += x34;
                    integral_14[g] += x33;
                    integral_15[g] += x26;
                    integral_16[g] += x58;
                    integral_17[g] += x45;
                    integral_18[g] += x31;
                    integral_19[g] += x15;
                    integral_20[g] += x26;
                    integral_21[g] += x17;
                    integral_22[g] += x16;
                    integral_23[g] += x54;
                    integral_24[g] += x34;
                    integral_25[g] += x26;
                    integral_26[g] += x58;
                    integral_27[g] += x16;
                    integral_28[g] += x54;
                    integral_29[g] += x49;
                    integral_30[g] += x33;
                    integral_31[g] += x58;
                    integral_32[g] += x45;
                    integral_33[g] += x54;
                    integral_34[g] += x49;
                    integral_35[g] += x50;
                    integral_36[g] += x15;
                    integral_37[g] += x17;
                    integral_38[g] += x16;
                    integral_39[g] += -4*x9 - PGy[g]*x14;
                    integral_40[g] += x18;
                    integral_41[g] += x56;
                    integral_42[g] += x26;
                    integral_43[g] += x16;
                    integral_44[g] += x54;
                    integral_45[g] += x18;
                    integral_46[g] += x56;
                    integral_47[g] += x48;
                    integral_48[g] += x58;
                    integral_49[g] += x54;
                    integral_50[g] += x49;
                    integral_51[g] += x56;
                    integral_52[g] += x48;
                    integral_53[g] += x51;
                    integral_54[g] += x45;
                    integral_55[g] += x49;
                    integral_56[g] += x50;
                    integral_57[g] += x48;
                    integral_58[g] += x51;
                    integral_59[g] += -4*x40 - PGz[g]*x46;

                } // for g
            } // if/else prim
        } // for prim

    }
} // VRR_F_D_SAMECENTER

template void VRR_F_D_SAMECENTER<float>(
        float* integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void VRR_F_D_SAMECENTER<double>(
        double* integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

