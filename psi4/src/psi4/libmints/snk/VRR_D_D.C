#include "VRR_D_D.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void VRR_D_D(
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
    if (44*sizeof(T) >= 280) {

        if (nprim < 1) return;

        // The gridpoints are processed in blocks.  The integrals of a block are contracted in
        // a buffer that stays in the L1 cache and are then copied to the output arrays.
        enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
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
                const T PBx = v_PAx[prim] + ABx;
                const T PBy = v_PAy[prim] + ABy;
                const T PBz = v_PAz[prim] + ABz;
                const T tpi = 0.5 / v_p[prim];
                const T tpi_2 = 2*tpi;
                const T tpi_p2 = tpi*tpi;

                // overwrite integrals for first primitive, accumulate otherwise
                if (prim==0) {
                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x1 = tpi*x0;
                        const T x2 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x3 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x4 = PAx*x2 + x1 - PGx[g]*x3;
                        const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x6 = tpi*x5;
                        const T x7 = x0 - x5;
                        const T x8 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x9 = PAx*x3 + x6 - PGx[g]*x8;
                        const T x10 = x2 - x3;
                        const T x11 = tpi*x10;
                        const T x12 = PBx*x4 + 2*x11 - PGx[g]*x9;
                        const T x13 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x14 = tpi*x13;
                        const T x15 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x16 = PAx*x8 + x14 - PGx[g]*x15;
                        const T x17 = x3 - x8;
                        const T x18 = tpi*x17;
                        const T x19 = PBx*x9 + 2*x18 - PGx[g]*x16;
                        const T x20 = PBx*x11 + tpi_p2*x7 - PGx[g]*x18;
                        const T x21 = x4 - x9;
                        block_0[g] = PBx*x12 + 2*x20 + tpi*x21 - PGx[g]*x19;
                        const T x22 = PBy*x4 - PGy[g]*x9;
                        const T x23 = PBy*x9 - PGy[g]*x16;
                        block_3[g] = tpi*x21 + PBy*x22 - PGy[g]*x23;
                        const T x24 = PBz*x4 - PGz[g]*x9;
                        const T x25 = PBz*x9 - PGz[g]*x16;
                        block_5[g] = tpi*x21 + PBz*x24 - PGz[g]*x25;
                        block_4[g] = PBy*x24 - PGy[g]*x25;
                        block_1[g] = PBy*x12 - PGy[g]*x19;
                        block_2[g] = PBz*x12 - PGz[g]*x19;
                        const T x26 = PBz*x11 - PGz[g]*x18;
                        const T x27 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x28 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x29 = PAx*x27 - PGx[g]*x28;
                        const T x30 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x31 = PAx*x28 - PGx[g]*x30;
                        const T x32 = x11 + PBz*x29 - PGz[g]*x31;
                        const T x33 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x34 = PAx*x30 - PGx[g]*x33;
                        const T x35 = x18 + PBz*x31 - PGz[g]*x34;
                        const T x36 = x29 - x31;
                        block_17[g] = x26 + tpi*x36 + PBz*x32 - PGz[g]*x35;
                        const T x37 = PBy*x29 - PGy[g]*x31;
                        const T x38 = PBy*x31 - PGy[g]*x34;
                        block_15[g] = tpi*x36 + PBy*x37 - PGy[g]*x38;
                        const T x39 = x27 - x28;
                        const T x40 = tpi*x39;
                        const T x41 = PBx*x29 + x40 - PGx[g]*x31;
                        const T x42 = x28 - x30;
                        const T x43 = tpi*x42;
                        const T x44 = PBx*x31 + x43 - PGx[g]*x34;
                        const T x45 = PBx*x40 - PGx[g]*x43;
                        block_12[g] = PBx*x41 + x45 + tpi*x36 - PGx[g]*x44;
                        block_13[g] = PBy*x41 - PGy[g]*x44;
                        const T x46 = tpi_p2*x7 + PBz*x40 - PGz[g]*x43;
                        block_14[g] = PBx*x32 + x46 - PGx[g]*x35;
                        block_16[g] = PBy*x32 - PGy[g]*x35;
                        const T x47 = x1 + PAz*x27 - PGz[g]*x28;
                        const T x48 = x6 + PAz*x28 - PGz[g]*x30;
                        const T x49 = 2*x40 + PBz*x47 - PGz[g]*x48;
                        const T x50 = x14 + PAz*x30 - PGz[g]*x33;
                        const T x51 = 2*x43 + PBz*x48 - PGz[g]*x50;
                        block_32[g] = PBx*x49 - PGx[g]*x51;
                        const T x52 = PAy*x27 - PGy[g]*x28;
                        const T x53 = PAy*x28 - PGy[g]*x30;
                        const T x54 = PAy*x30 - PGy[g]*x33;
                        const T x55 = PBx*x52 - PGx[g]*x53;
                        const T x56 = PBx*x53 - PGx[g]*x54;
                        block_25[g] = x45 + PBy*x55 - PGy[g]*x56;
                        const T x57 = x52 - x53;
                        block_24[g] = PBx*x55 + tpi*x57 - PGx[g]*x56;
                        block_34[g] = PBy*x49 - PGy[g]*x51;
                        const T x58 = x47 - x48;
                        block_35[g] = 2*x46 + tpi*x58 + PBz*x49 - PGz[g]*x51;
                        const T x59 = PBx*x47 - PGx[g]*x48;
                        const T x60 = PBx*x48 - PGx[g]*x50;
                        block_30[g] = PBx*x59 + tpi*x58 - PGx[g]*x60;
                        const T x61 = PBy*x47 - PGy[g]*x48;
                        const T x62 = PBy*x48 - PGy[g]*x50;
                        block_33[g] = tpi*x58 + PBy*x61 - PGy[g]*x62;
                        block_31[g] = PBx*x61 - PGx[g]*x62;
                        const T x63 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x64 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x65 = x63 - x64;
                        const T x66 = tpi*x65;
                        const T x67 = x66 + PBz*x52 - PGz[g]*x53;
                        const T x68 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x69 = x64 - x68;
                        const T x70 = tpi*x69;
                        const T x71 = x70 + PBz*x53 - PGz[g]*x54;
                        block_28[g] = x46 + PBy*x67 - PGy[g]*x71;
                        const T x72 = PBz*x66 - PGz[g]*x70;
                        block_29[g] = x72 + tpi*x57 + PBz*x67 - PGz[g]*x71;
                        block_26[g] = PBx*x67 - PGx[g]*x71;
                        const T x73 = PAx*x63 - PGx[g]*x64;
                        const T x74 = PAx*x64 - PGx[g]*x68;
                        const T x75 = PBx*x73 + x66 - PGx[g]*x74;
                        const T x76 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x77 = PAx*x68 - PGx[g]*x76;
                        const T x78 = PBx*x74 + x70 - PGx[g]*x77;
                        block_7[g] = x20 + PBy*x75 - PGy[g]*x78;
                        const T x79 = x73 - x74;
                        const T x80 = PBx*x66 - PGx[g]*x70;
                        block_6[g] = PBx*x75 + x80 + tpi*x79 - PGx[g]*x78;
                        const T x81 = PBz*x73 - PGz[g]*x74;
                        const T x82 = PBz*x74 - PGz[g]*x77;
                        block_8[g] = PBx*x81 + x72 - PGx[g]*x82;
                        block_10[g] = PBy*x81 + x26 - PGy[g]*x82;
                        block_11[g] = tpi*x79 + PBz*x81 - PGz[g]*x82;
                        const T x83 = x1 + PBy*x63 - PGy[g]*x64;
                        const T x84 = x6 + PBy*x64 - PGy[g]*x68;
                        const T x85 = PBy*x1 - PGy[g]*x6;
                        const T x86 = x1 + PAy*x63 - PGy[g]*x64;
                        const T x87 = x66 + PAy*x83 + x85 - PGy[g]*x84;
                        const T x88 = x66 + PBy*x83 + x85 - PGy[g]*x84;
                        const T x89 = x14 + PBy*x68 - PGy[g]*x76;
                        const T x90 = PBy*x6 - PGy[g]*x14;
                        const T x91 = x6 + PAy*x64 - PGy[g]*x68;
                        const T x92 = x14 + PAy*x68 - PGy[g]*x76;
                        const T x93 = x70 + PAy*x84 + x90 - PGy[g]*x89;
                        const T x94 = x70 + PBy*x84 + x90 - PGy[g]*x89;
                        const T x95 = x83 - x84;
                        const T x96 = x86 - x91;
                        block_21[g] = tpi*x96 + PBy*x87 + tpi_2*x95 - PGy[g]*x93;
                        const T x97 = PBx*x86 - PGx[g]*x91;
                        const T x98 = PBx*x91 - PGx[g]*x92;
                        block_18[g] = PBx*x97 + tpi*x96 - PGx[g]*x98;
                        const T x99 = PBz*x86 - PGz[g]*x91;
                        const T x100 = PBz*x91 - PGz[g]*x92;
                        block_23[g] = tpi*x96 + PBz*x99 - PGz[g]*x100;
                        block_20[g] = PBx*x99 - PGx[g]*x100;
                        block_9[g] = PAx*x88 - PGx[g]*x94;
                        block_27[g] = PAz*x88 - PGz[g]*x94;
                        block_19[g] = PBx*x87 - PGx[g]*x93;
                        block_22[g] = PBz*x87 - PGz[g]*x93;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x1 = tpi*x0;
                        const T x2 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x3 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x4 = PAx*x2 + x1 - PGx[g]*x3;
                        const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x6 = tpi*x5;
                        const T x7 = x0 - x5;
                        const T x8 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x9 = PAx*x3 + x6 - PGx[g]*x8;
                        const T x10 = x2 - x3;
                        const T x11 = tpi*x10;
                        const T x12 = PBx*x4 + 2*x11 - PGx[g]*x9;
                        const T x13 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x14 = tpi*x13;
                        const T x15 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x16 = PAx*x8 + x14 - PGx[g]*x15;
                        const T x17 = x3 - x8;
                        const T x18 = tpi*x17;
                        const T x19 = PBx*x9 + 2*x18 - PGx[g]*x16;
                        const T x20 = PBx*x11 + tpi_p2*x7 - PGx[g]*x18;
                        const T x21 = x4 - x9;
                        block_0[g] += PBx*x12 + 2*x20 + tpi*x21 - PGx[g]*x19;
                        const T x22 = PBy*x4 - PGy[g]*x9;
                        const T x23 = PBy*x9 - PGy[g]*x16;
                        block_3[g] += tpi*x21 + PBy*x22 - PGy[g]*x23;
                        const T x24 = PBz*x4 - PGz[g]*x9;
                        const T x25 = PBz*x9 - PGz[g]*x16;
                        block_5[g] += tpi*x21 + PBz*x24 - PGz[g]*x25;
                        block_4[g] += PBy*x24 - PGy[g]*x25;
                        block_1[g] += PBy*x12 - PGy[g]*x19;
                        block_2[g] += PBz*x12 - PGz[g]*x19;
                        const T x26 = PBz*x11 - PGz[g]*x18;
                        const T x27 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x28 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x29 = PAx*x27 - PGx[g]*x28;
                        const T x30 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x31 = PAx*x28 - PGx[g]*x30;
                        const T x32 = x11 + PBz*x29 - PGz[g]*x31;
                        const T x33 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x34 = PAx*x30 - PGx[g]*x33;
                        const T x35 = x18 + PBz*x31 - PGz[g]*x34;
                        const T x36 = x29 - x31;
                        block_17[g] += x26 + tpi*x36 + PBz*x32 - PGz[g]*x35;
                        const T x37 = PBy*x29 - PGy[g]*x31;
                        const T x38 = PBy*x31 - PGy[g]*x34;
                        block_15[g] += tpi*x36 + PBy*x37 - PGy[g]*x38;
                        const T x39 = x27 - x28;
                        const T x40 = tpi*x39;
                        const T x41 = PBx*x29 + x40 - PGx[g]*x31;
                        const T x42 = x28 - x30;
                        const T x43 = tpi*x42;
                        const T x44 = PBx*x31 + x43 - PGx[g]*x34;
                        const T x45 = PBx*x40 - PGx[g]*x43;
                        block_12[g] += PBx*x41 + x45 + tpi*x36 - PGx[g]*x44;
                        block_13[g] += PBy*x41 - PGy[g]*x44;
                        const T x46 = tpi_p2*x7 + PBz*x40 - PGz[g]*x43;
                        block_14[g] += PBx*x32 + x46 - PGx[g]*x35;
                        block_16[g] += PBy*x32 - PGy[g]*x35;
                        const T x47 = x1 + PAz*x27 - PGz[g]*x28;
                        const T x48 = x6 + PAz*x28 - PGz[g]*x30;
                        const T x49 = 2*x40 + PBz*x47 - PGz[g]*x48;
                        const T x50 = x14 + PAz*x30 - PGz[g]*x33;
                        const T x51 = 2*x43 + PBz*x48 - PGz[g]*x50;
                        block_32[g] += PBx*x49 - PGx[g]*x51;
                        const T x52 = PAy*x27 - PGy[g]*x28;
                        const T x53 = PAy*x28 - PGy[g]*x30;
                        const T x54 = PAy*x30 - PGy[g]*x33;
                        const T x55 = PBx*x52 - PGx[g]*x53;
                        const T x56 = PBx*x53 - PGx[g]*x54;
                        block_25[g] += x45 + PBy*x55 - PGy[g]*x56;
                        const T x57 = x52 - x53;
                        block_24[g] += PBx*x55 + tpi*x57 - PGx[g]*x56;
                        block_34[g] += PBy*x49 - PGy[g]*x51;
                        const T x58 = x47 - x48;
                        block_35[g] += 2*x46 + tpi*x58 + PBz*x49 - PGz[g]*x51;
                        const T x59 = PBx*x47 - PGx[g]*x48;
                        const T x60 = PBx*x48 - PGx[g]*x50;
                        block_30[g] += PBx*x59 + tpi*x58 - PGx[g]*x60;
                        const T x61 = PBy*x47 - PGy[g]*x48;
                        const T x62 = PBy*x48 - PGy[g]*x50;
                        block_33[g] += tpi*x58 + PBy*x61 - PGy[g]*x62;
                        block_31[g] += PBx*x61 - PGx[g]*x62;
                        const T x63 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x64 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x65 = x63 - x64;
                        const T x66 = tpi*x65;
                        const T x67 = x66 + PBz*x52 - PGz[g]*x53;
                        const T x68 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x69 = x64 - x68;
                        const T x70 = tpi*x69;
                        const T x71 = x70 + PBz*x53 - PGz[g]*x54;
                        block_28[g] += x46 + PBy*x67 - PGy[g]*x71;
                        const T x72 = PBz*x66 - PGz[g]*x70;
                        block_29[g] += x72 + tpi*x57 + PBz*x67 - PGz[g]*x71;
                        block_26[g] += PBx*x67 - PGx[g]*x71;
                        const T x73 = PAx*x63 - PGx[g]*x64;
                        const T x74 = PAx*x64 - PGx[g]*x68;
                        const T x75 = PBx*x73 + x66 - PGx[g]*x74;
                        const T x76 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x77 = PAx*x68 - PGx[g]*x76;
                        const T x78 = PBx*x74 + x70 - PGx[g]*x77;
                        block_7[g] += x20 + PBy*x75 - PGy[g]*x78;
                        const T x79 = x73 - x74;
                        const T x80 = PBx*x66 - PGx[g]*x70;
                        block_6[g] += PBx*x75 + x80 + tpi*x79 - PGx[g]*x78;
                        const T x81 = PBz*x73 - PGz[g]*x74;
                        const T x82 = PBz*x74 - PGz[g]*x77;
                        block_8[g] += PBx*x81 + x72 - PGx[g]*x82;
                        block_10[g] += PBy*x81 + x26 - PGy[g]*x82;
                        block_11[g] += tpi*x79 + PBz*x81 - PGz[g]*x82;
                        const T x83 = x1 + PBy*x63 - PGy[g]*x64;
                        const T x84 = x6 + PBy*x64 - PGy[g]*x68;
                        const T x85 = PBy*x1 - PGy[g]*x6;
                        const T x86 = x1 + PAy*x63 - PGy[g]*x64;
                        const T x87 = x66 + PAy*x83 + x85 - PGy[g]*x84;
                        const T x88 = x66 + PBy*x83 + x85 - PGy[g]*x84;
                        const T x89 = x14 + PBy*x68 - PGy[g]*x76;
                        const T x90 = PBy*x6 - PGy[g]*x14;
                        const T x91 = x6 + PAy*x64 - PGy[g]*x68;
                        const T x92 = x14 + PAy*x68 - PGy[g]*x76;
                        const T x93 = x70 + PAy*x84 + x90 - PGy[g]*x89;
                        const T x94 = x70 + PBy*x84 + x90 - PGy[g]*x89;
                        const T x95 = x83 - x84;
                        const T x96 = x86 - x91;
                        block_21[g] += tpi*x96 + PBy*x87 + tpi_2*x95 - PGy[g]*x93;
                        const T x97 = PBx*x86 - PGx[g]*x91;
                        const T x98 = PBx*x91 - PGx[g]*x92;
                        block_18[g] += PBx*x97 + tpi*x96 - PGx[g]*x98;
                        const T x99 = PBz*x86 - PGz[g]*x91;
                        const T x100 = PBz*x91 - PGz[g]*x92;
                        block_23[g] += tpi*x96 + PBz*x99 - PGz[g]*x100;
                        block_20[g] += PBx*x99 - PGx[g]*x100;
                        block_9[g] += PAx*x88 - PGx[g]*x94;
                        block_27[g] += PAz*x88 - PGz[g]*x94;
                        block_19[g] += PBx*x87 - PGx[g]*x93;
                        block_22[g] += PBz*x87 - PGz[g]*x93;

                    } // for g
                } // if/else prim
            } // for prim

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
            const T tpi_p2 = tpi*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x1 = tpi*x0;
                    const T x2 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x3 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x4 = PAx*x2 + x1 - PGx[g]*x3;
                    const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x6 = tpi*x5;
                    const T x7 = x0 - x5;
                    const T x8 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x9 = PAx*x3 + x6 - PGx[g]*x8;
                    const T x10 = x2 - x3;
                    const T x11 = tpi*x10;
                    const T x12 = PBx*x4 + 2*x11 - PGx[g]*x9;
                    const T x13 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x14 = tpi*x13;
                    const T x15 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x16 = PAx*x8 + x14 - PGx[g]*x15;
                    const T x17 = x3 - x8;
                    const T x18 = tpi*x17;
                    const T x19 = PBx*x9 + 2*x18 - PGx[g]*x16;
                    const T x20 = PBx*x11 + tpi_p2*x7 - PGx[g]*x18;
                    const T x21 = x4 - x9;
                    const T x22 = PBy*x4 - PGy[g]*x9;
                    const T x23 = PBy*x9 - PGy[g]*x16;
                    const T x24 = PBz*x4 - PGz[g]*x9;
                    const T x25 = PBz*x9 - PGz[g]*x16;
                    const T x26 = PBz*x11 - PGz[g]*x18;
                    const T x27 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x28 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x29 = PAx*x27 - PGx[g]*x28;
                    const T x30 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x31 = PAx*x28 - PGx[g]*x30;
                    const T x32 = x11 + PBz*x29 - PGz[g]*x31;
                    const T x33 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x34 = PAx*x30 - PGx[g]*x33;
                    const T x35 = x18 + PBz*x31 - PGz[g]*x34;
                    const T x36 = x29 - x31;
                    const T x37 = PBy*x29 - PGy[g]*x31;
                    const T x38 = PBy*x31 - PGy[g]*x34;
                    const T x39 = x27 - x28;
                    const T x40 = tpi*x39;
                    const T x41 = PBx*x29 + x40 - PGx[g]*x31;
                    const T x42 = x28 - x30;
                    const T x43 = tpi*x42;
                    const T x44 = PBx*x31 + x43 - PGx[g]*x34;
                    const T x45 = PBx*x40 - PGx[g]*x43;
                    const T x46 = tpi_p2*x7 + PBz*x40 - PGz[g]*x43;
                    const T x47 = x1 + PAz*x27 - PGz[g]*x28;
                    const T x48 = x6 + PAz*x28 - PGz[g]*x30;
                    const T x49 = 2*x40 + PBz*x47 - PGz[g]*x48;
                    const T x50 = x14 + PAz*x30 - PGz[g]*x33;
                    const T x51 = 2*x43 + PBz*x48 - PGz[g]*x50;
                    const T x52 = PAy*x27 - PGy[g]*x28;
                    const T x53 = PAy*x28 - PGy[g]*x30;
                    const T x54 = PAy*x30 - PGy[g]*x33;
                    const T x55 = PBx*x52 - PGx[g]*x53;
                    const T x56 = PBx*x53 - PGx[g]*x54;
                    const T x57 = x52 - x53;
                    const T x58 = x47 - x48;
                    const T x59 = PBx*x47 - PGx[g]*x48;
                    const T x60 = PBx*x48 - PGx[g]*x50;
                    const T x61 = PBy*x47 - PGy[g]*x48;
                    const T x62 = PBy*x48 - PGy[g]*x50;
                    const T x63 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x64 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x65 = x63 - x64;
                    const T x66 = tpi*x65;
                    const T x67 = x66 + PBz*x52 - PGz[g]*x53;
                    const T x68 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x69 = x64 - x68;
                    const T x70 = tpi*x69;
                    const T x71 = x70 + PBz*x53 - PGz[g]*x54;
                    const T x72 = PBz*x66 - PGz[g]*x70;
                    const T x73 = PAx*x63 - PGx[g]*x64;
                    const T x74 = PAx*x64 - PGx[g]*x68;
                    const T x75 = PBx*x73 + x66 - PGx[g]*x74;
                    const T x76 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x77 = PAx*x68 - PGx[g]*x76;
                    const T x78 = PBx*x74 + x70 - PGx[g]*x77;
                    const T x79 = x73 - x74;
                    const T x80 = PBx*x66 - PGx[g]*x70;
                    const T x81 = PBz*x73 - PGz[g]*x74;
                    const T x82 = PBz*x74 - PGz[g]*x77;
                    const T x83 = x1 + PBy*x63 - PGy[g]*x64;
                    const T x84 = x6 + PBy*x64 - PGy[g]*x68;
                    const T x85 = PBy*x1 - PGy[g]*x6;
                    const T x86 = x1 + PAy*x63 - PGy[g]*x64;
                    const T x87 = x66 + PAy*x83 + x85 - PGy[g]*x84;
                    const T x88 = x66 + PBy*x83 + x85 - PGy[g]*x84;
                    const T x89 = x14 + PBy*x68 - PGy[g]*x76;
                    const T x90 = PBy*x6 - PGy[g]*x14;
                    const T x91 = x6 + PAy*x64 - PGy[g]*x68;
                    const T x92 = x14 + PAy*x68 - PGy[g]*x76;
                    const T x93 = x70 + PAy*x84 + x90 - PGy[g]*x89;
                    const T x94 = x70 + PBy*x84 + x90 - PGy[g]*x89;
                    const T x95 = x83 - x84;
                    const T x96 = x86 - x91;
                    const T x97 = PBx*x86 - PGx[g]*x91;
                    const T x98 = PBx*x91 - PGx[g]*x92;
                    const T x99 = PBz*x86 - PGz[g]*x91;
                    const T x100 = PBz*x91 - PGz[g]*x92;

                    integral_0[g] = PBx*x12 + 2*x20 + tpi*x21 - PGx[g]*x19;
                    integral_1[g] = PBy*x12 - PGy[g]*x19;
                    integral_2[g] = PBz*x12 - PGz[g]*x19;
                    integral_3[g] = tpi*x21 + PBy*x22 - PGy[g]*x23;
                    integral_4[g] = PBy*x24 - PGy[g]*x25;
                    integral_5[g] = tpi*x21 + PBz*x24 - PGz[g]*x25;
                    integral_6[g] = PBx*x75 + x80 + tpi*x79 - PGx[g]*x78;
                    integral_7[g] = x20 + PBy*x75 - PGy[g]*x78;
                    integral_8[g] = PBx*x81 + x72 - PGx[g]*x82;
                    integral_9[g] = PAx*x88 - PGx[g]*x94;
                    integral_10[g] = PBy*x81 + x26 - PGy[g]*x82;
                    integral_11[g] = tpi*x79 + PBz*x81 - PGz[g]*x82;
                    integral_12[g] = PBx*x41 + x45 + tpi*x36 - PGx[g]*x44;
                    integral_13[g] = PBy*x41 - PGy[g]*x44;
                    integral_14[g] = PBx*x32 + x46 - PGx[g]*x35;
                    integral_15[g] = tpi*x36 + PBy*x37 - PGy[g]*x38;
                    integral_16[g] = PBy*x32 - PGy[g]*x35;
                    integral_17[g] = x26 + tpi*x36 + PBz*x32 - PGz[g]*x35;
                    integral_18[g] = PBx*x97 + tpi*x96 - PGx[g]*x98;
                    integral_19[g] = PBx*x87 - PGx[g]*x93;
                    integral_20[g] = PBx*x99 - PGx[g]*x100;
                    integral_21[g] = tpi*x96 + PBy*x87 + tpi_2*x95 - PGy[g]*x93;
                    integral_22[g] = PBz*x87 - PGz[g]*x93;
                    integral_23[g] = tpi*x96 + PBz*x99 - PGz[g]*x100;
                    integral_24[g] = PBx*x55 + tpi*x57 - PGx[g]*x56;
                    integral_25[g] = x45 + PBy*x55 - PGy[g]*x56;
                    integral_26[g] = PBx*x67 - PGx[g]*x71;
                    integral_27[g] = PAz*x88 - PGz[g]*x94;
                    integral_28[g] = x46 + PBy*x67 - PGy[g]*x71;
                    integral_29[g] = x72 + tpi*x57 + PBz*x67 - PGz[g]*x71;
                    integral_30[g] = PBx*x59 + tpi*x58 - PGx[g]*x60;
                    integral_31[g] = PBx*x61 - PGx[g]*x62;
                    integral_32[g] = PBx*x49 - PGx[g]*x51;
                    integral_33[g] = tpi*x58 + PBy*x61 - PGy[g]*x62;
                    integral_34[g] = PBy*x49 - PGy[g]*x51;
                    integral_35[g] = 2*x46 + tpi*x58 + PBz*x49 - PGz[g]*x51;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x1 = tpi*x0;
                    const T x2 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x3 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x4 = PAx*x2 + x1 - PGx[g]*x3;
                    const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x6 = tpi*x5;
                    const T x7 = x0 - x5;
                    const T x8 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x9 = PAx*x3 + x6 - PGx[g]*x8;
                    const T x10 = x2 - x3;
                    const T x11 = tpi*x10;
                    const T x12 = PBx*x4 + 2*x11 - PGx[g]*x9;
                    const T x13 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x14 = tpi*x13;
                    const T x15 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x16 = PAx*x8 + x14 - PGx[g]*x15;
                    const T x17 = x3 - x8;
                    const T x18 = tpi*x17;
                    const T x19 = PBx*x9 + 2*x18 - PGx[g]*x16;
                    const T x20 = PBx*x11 + tpi_p2*x7 - PGx[g]*x18;
                    const T x21 = x4 - x9;
                    const T x22 = PBy*x4 - PGy[g]*x9;
                    const T x23 = PBy*x9 - PGy[g]*x16;
                    const T x24 = PBz*x4 - PGz[g]*x9;
                    const T x25 = PBz*x9 - PGz[g]*x16;
                    const T x26 = PBz*x11 - PGz[g]*x18;
                    const T x27 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x28 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x29 = PAx*x27 - PGx[g]*x28;
                    const T x30 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x31 = PAx*x28 - PGx[g]*x30;
                    const T x32 = x11 + PBz*x29 - PGz[g]*x31;
                    const T x33 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x34 = PAx*x30 - PGx[g]*x33;
                    const T x35 = x18 + PBz*x31 - PGz[g]*x34;
                    const T x36 = x29 - x31;
                    const T x37 = PBy*x29 - PGy[g]*x31;
                    const T x38 = PBy*x31 - PGy[g]*x34;
                    const T x39 = x27 - x28;
                    const T x40 = tpi*x39;
                    const T x41 = PBx*x29 + x40 - PGx[g]*x31;
                    const T x42 = x28 - x30;
                    const T x43 = tpi*x42;
                    const T x44 = PBx*x31 + x43 - PGx[g]*x34;
                    const T x45 = PBx*x40 - PGx[g]*x43;
                    const T x46 = tpi_p2*x7 + PBz*x40 - PGz[g]*x43;
                    const T x47 = x1 + PAz*x27 - PGz[g]*x28;
                    const T x48 = x6 + PAz*x28 - PGz[g]*x30;
                    const T x49 = 2*x40 + PBz*x47 - PGz[g]*x48;
                    const T x50 = x14 + PAz*x30 - PGz[g]*x33;
                    const T x51 = 2*x43 + PBz*x48 - PGz[g]*x50;
                    const T x52 = PAy*x27 - PGy[g]*x28;
                    const T x53 = PAy*x28 - PGy[g]*x30;
                    const T x54 = PAy*x30 - PGy[g]*x33;
                    const T x55 = PBx*x52 - PGx[g]*x53;
                    const T x56 = PBx*x53 - PGx[g]*x54;
                    const T x57 = x52 - x53;
                    const T x58 = x47 - x48;
                    const T x59 = PBx*x47 - PGx[g]*x48;
                    const T x60 = PBx*x48 - PGx[g]*x50;
                    const T x61 = PBy*x47 - PGy[g]*x48;
                    const T x62 = PBy*x48 - PGy[g]*x50;
                    const T x63 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x64 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x65 = x63 - x64;
                    const T x66 = tpi*x65;
                    const T x67 = x66 + PBz*x52 - PGz[g]*x53;
                    const T x68 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x69 = x64 - x68;
                    const T x70 = tpi*x69;
                    const T x71 = x70 + PBz*x53 - PGz[g]*x54;
                    const T x72 = PBz*x66 - PGz[g]*x70;
                    const T x73 = PAx*x63 - PGx[g]*x64;
                    const T x74 = PAx*x64 - PGx[g]*x68;
                    const T x75 = PBx*x73 + x66 - PGx[g]*x74;
                    const T x76 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x77 = PAx*x68 - PGx[g]*x76;
                    const T x78 = PBx*x74 + x70 - PGx[g]*x77;
                    const T x79 = x73 - x74;
                    const T x80 = PBx*x66 - PGx[g]*x70;
                    const T x81 = PBz*x73 - PGz[g]*x74;
                    const T x82 = PBz*x74 - PGz[g]*x77;
                    const T x83 = x1 + PBy*x63 - PGy[g]*x64;
                    const T x84 = x6 + PBy*x64 - PGy[g]*x68;
                    const T x85 = PBy*x1 - PGy[g]*x6;
                    const T x86 = x1 + PAy*x63 - PGy[g]*x64;
                    const T x87 = x66 + PAy*x83 + x85 - PGy[g]*x84;
                    const T x88 = x66 + PBy*x83 + x85 - PGy[g]*x84;
                    const T x89 = x14 + PBy*x68 - PGy[g]*x76;
                    const T x90 = PBy*x6 - PGy[g]*x14;
                    const T x91 = x6 + PAy*x64 - PGy[g]*x68;
                    const T x92 = x14 + PAy*x68 - PGy[g]*x76;
                    const T x93 = x70 + PAy*x84 + x90 - PGy[g]*x89;
                    const T x94 = x70 + PBy*x84 + x90 - PGy[g]*x89;
                    const T x95 = x83 - x84;
                    const T x96 = x86 - x91;
                    const T x97 = PBx*x86 - PGx[g]*x91;
                    const T x98 = PBx*x91 - PGx[g]*x92;
                    const T x99 = PBz*x86 - PGz[g]*x91;
                    const T x100 = PBz*x91 - PGz[g]*x92;

                    integral_0[g] += PBx*x12 + 2*x20 + tpi*x21 - PGx[g]*x19;
                    integral_1[g] += PBy*x12 - PGy[g]*x19;
                    integral_2[g] += PBz*x12 - PGz[g]*x19;
                    integral_3[g] += tpi*x21 + PBy*x22 - PGy[g]*x23;
                    integral_4[g] += PBy*x24 - PGy[g]*x25;
                    integral_5[g] += tpi*x21 + PBz*x24 - PGz[g]*x25;
                    integral_6[g] += PBx*x75 + x80 + tpi*x79 - PGx[g]*x78;
                    integral_7[g] += x20 + PBy*x75 - PGy[g]*x78;
                    integral_8[g] += PBx*x81 + x72 - PGx[g]*x82;
                    integral_9[g] += PAx*x88 - PGx[g]*x94;
                    integral_10[g] += PBy*x81 + x26 - PGy[g]*x82;
                    integral_11[g] += tpi*x79 + PBz*x81 - PGz[g]*x82;
                    integral_12[g] += PBx*x41 + x45 + tpi*x36 - PGx[g]*x44;
                    integral_13[g] += PBy*x41 - PGy[g]*x44;
                    integral_14[g] += PBx*x32 + x46 - PGx[g]*x35;
                    integral_15[g] += tpi*x36 + PBy*x37 - PGy[g]*x38;
                    integral_16[g] += PBy*x32 - PGy[g]*x35;
                    integral_17[g] += x26 + tpi*x36 + PBz*x32 - PGz[g]*x35;
                    integral_18[g] += PBx*x97 + tpi*x96 - PGx[g]*x98;
                    integral_19[g] += PBx*x87 - PGx[g]*x93;
                    integral_20[g] += PBx*x99 - PGx[g]*x100;
                    integral_21[g] += tpi*x96 + PBy*x87 + tpi_2*x95 - PGy[g]*x93;
                    integral_22[g] += PBz*x87 - PGz[g]*x93;
                    integral_23[g] += tpi*x96 + PBz*x99 - PGz[g]*x100;
                    integral_24[g] += PBx*x55 + tpi*x57 - PGx[g]*x56;
                    integral_25[g] += x45 + PBy*x55 - PGy[g]*x56;
                    integral_26[g] += PBx*x67 - PGx[g]*x71;
                    integral_27[g] += PAz*x88 - PGz[g]*x94;
                    integral_28[g] += x46 + PBy*x67 - PGy[g]*x71;
                    integral_29[g] += x72 + tpi*x57 + PBz*x67 - PGz[g]*x71;
                    integral_30[g] += PBx*x59 + tpi*x58 - PGx[g]*x60;
                    integral_31[g] += PBx*x61 - PGx[g]*x62;
                    integral_32[g] += PBx*x49 - PGx[g]*x51;
                    integral_33[g] += tpi*x58 + PBy*x61 - PGy[g]*x62;
                    integral_34[g] += PBy*x49 - PGy[g]*x51;
                    integral_35[g] += 2*x46 + tpi*x58 + PBz*x49 - PGz[g]*x51;

                } // for g
            } // if/else prim
        } // for prim

    }
} // VRR_D_D

template void VRR_D_D<float>(
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

template void VRR_D_D<double>(
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

