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

                // overwrite integrals for first primitive, accumulate otherwise
                if (prim==0) {
                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x2 = PBx*x0 - PGx[g]*x1;
                        const T x3 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x4 = PBx*x1 - PGx[g]*x3;
                        const T x5 = x2 - x4;
                        const T x6 = PAz*x2 - PGz[g]*x4;
                        const T x7 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x8 = PBx*x3 - PGx[g]*x7;
                        const T x9 = PAz*x4 - PGz[g]*x8;
                        block_31[g] = tpi*x5 + PAz*x6 - PGz[g]*x9;
                        const T x10 = x0 - x1;
                        const T x11 = PAx*x2 + tpi*x10 - PGx[g]*x4;
                        const T x12 = x1 - x3;
                        const T x13 = PAx*x4 + tpi*x12 - PGx[g]*x8;
                        block_13[g] = PAz*x11 - PGz[g]*x13;
                        const T x14 = PAx*x10 - PGx[g]*x12;
                        block_1[g] = PAx*x11 + tpi*x5 + tpi*x14 - PGx[g]*x13;
                        const T x15 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x16 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x17 = x15 - x16;
                        const T x18 = tpi*x17 + PAy*x2 - PGy[g]*x4;
                        const T x19 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x20 = x16 - x19;
                        const T x21 = tpi*x20 + PAy*x4 - PGy[g]*x8;
                        const T x22 = PAy*x17 - PGy[g]*x20;
                        block_19[g] = tpi*x5 + PAy*x18 + tpi*x22 - PGy[g]*x21;
                        block_25[g] = PAz*x18 - PGz[g]*x21;
                        const T x23 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x24 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x25 = x23 - x24;
                        const T x26 = tpi*x25 + PAy*x10 - PGy[g]*x12;
                        block_7[g] = PAx*x18 + tpi*x26 - PGx[g]*x21;
                        const T x27 = tpi*x23 + PBy*x0 - PGy[g]*x1;
                        const T x28 = tpi*x24 + PBy*x1 - PGy[g]*x3;
                        const T x29 = PAx*x27 - PGx[g]*x28;
                        const T x30 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x31 = tpi*x30 + PBy*x3 - PGy[g]*x7;
                        const T x32 = PAx*x28 - PGx[g]*x31;
                        const T x33 = x27 - x28;
                        block_3[g] = PAx*x29 + tpi*x33 - PGx[g]*x32;
                        const T x34 = tpi_2*x10 + PAy*x27 - PGy[g]*x28;
                        const T x35 = tpi_2*x12 + PAy*x28 - PGy[g]*x31;
                        block_21[g] = tpi*x33 + tpi_2*x26 + PAy*x34 - PGy[g]*x35;
                        const T x36 = PAz*x27 - PGz[g]*x28;
                        const T x37 = PAz*x28 - PGz[g]*x31;
                        block_33[g] = tpi*x33 + PAz*x36 - PGz[g]*x37;
                        block_15[g] = PAx*x36 - PGx[g]*x37;
                        block_9[g] = PAx*x34 - PGx[g]*x35;
                        block_27[g] = PAz*x34 - PGz[g]*x35;
                        const T x38 = PBx*x15 + tpi*x23 - PGx[g]*x16;
                        const T x39 = PBx*x16 + tpi*x24 - PGx[g]*x19;
                        const T x40 = x38 - x39;
                        const T x41 = PAy*x38 - PGy[g]*x39;
                        const T x42 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x43 = PBx*x19 + tpi*x30 - PGx[g]*x42;
                        const T x44 = PAy*x39 - PGy[g]*x43;
                        block_18[g] = tpi*x40 + PAy*x41 - PGy[g]*x44;
                        const T x45 = PAx*x38 + tpi_2*x17 - PGx[g]*x39;
                        const T x46 = PAx*x39 + tpi_2*x20 - PGx[g]*x43;
                        const T x47 = PAx*x17 + tpi*x25 - PGx[g]*x20;
                        block_0[g] = PAx*x45 + tpi*x40 + tpi_2*x47 - PGx[g]*x46;
                        const T x48 = PAz*x38 - PGz[g]*x39;
                        const T x49 = PAz*x39 - PGz[g]*x43;
                        block_30[g] = tpi*x40 + PAz*x48 - PGz[g]*x49;
                        block_24[g] = PAy*x48 - PGy[g]*x49;
                        block_6[g] = PAy*x45 - PGy[g]*x46;
                        block_12[g] = PAz*x45 - PGz[g]*x46;
                        const T x50 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x51 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x52 = tpi*x23 + PBz*x50 - PGz[g]*x51;
                        const T x53 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x54 = tpi*x24 + PBz*x51 - PGz[g]*x53;
                        const T x55 = PAx*x52 - PGx[g]*x54;
                        const T x56 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x57 = tpi*x30 + PBz*x53 - PGz[g]*x56;
                        const T x58 = PAx*x54 - PGx[g]*x57;
                        const T x59 = x52 - x54;
                        block_5[g] = PAx*x55 + tpi*x59 - PGx[g]*x58;
                        const T x60 = PAy*x52 - PGy[g]*x54;
                        const T x61 = PAy*x54 - PGy[g]*x57;
                        block_11[g] = PAx*x60 - PGx[g]*x61;
                        block_23[g] = tpi*x59 + PAy*x60 - PGy[g]*x61;
                        const T x62 = x50 - x51;
                        const T x63 = tpi_2*x62 + PAz*x52 - PGz[g]*x54;
                        const T x64 = x51 - x53;
                        const T x65 = tpi_2*x64 + PAz*x54 - PGz[g]*x57;
                        block_17[g] = PAx*x63 - PGx[g]*x65;
                        block_29[g] = PAy*x63 - PGy[g]*x65;
                        const T x66 = tpi*x23 + PAz*x50 - PGz[g]*x51;
                        const T x67 = tpi*x24 + PAz*x51 - PGz[g]*x53;
                        const T x68 = x66 - x67;
                        block_35[g] = tpi*x59 + tpi_2*x68 + PAz*x63 - PGz[g]*x65;
                        const T x69 = PBx*x66 - PGx[g]*x67;
                        const T x70 = tpi*x30 + PAz*x53 - PGz[g]*x56;
                        const T x71 = PBx*x67 - PGx[g]*x70;
                        const T x72 = PBx*x62 - PGx[g]*x64;
                        const T x73 = PAz*x17 - PGz[g]*x20;
                        block_32[g] = PAz*x69 + tpi*x72 + tpi*x73 - PGz[g]*x71;
                        const T x74 = PBy*x66 - PGy[g]*x67;
                        const T x75 = PBy*x67 - PGy[g]*x70;
                        const T x76 = PBy*x62 - PGy[g]*x64;
                        const T x77 = PAz*x10 - PGz[g]*x12;
                        block_34[g] = PAz*x74 + tpi*x76 + tpi*x77 - PGz[g]*x75;
                        block_14[g] = PAx*x69 + tpi*x68 - PGx[g]*x71;
                        block_26[g] = PAy*x69 - PGy[g]*x71;
                        block_28[g] = tpi*x68 + PAy*x74 - PGy[g]*x75;
                        block_16[g] = PAx*x74 - PGx[g]*x75;
                        const T x78 = PAx*x50 - PGx[g]*x51;
                        const T x79 = PAx*x51 - PGx[g]*x53;
                        const T x80 = PAx*x78 + tpi*x62 - PGx[g]*x79;
                        const T x81 = PAx*x53 - PGx[g]*x56;
                        const T x82 = PAx*x79 + tpi*x64 - PGx[g]*x81;
                        const T x83 = x78 - x79;
                        block_2[g] = PBx*x80 + tpi_2*x83 - PGx[g]*x82;
                        block_4[g] = PBy*x80 - PGy[g]*x82;
                        const T x84 = PAy*x50 - PGy[g]*x51;
                        const T x85 = PAy*x51 - PGy[g]*x53;
                        const T x86 = tpi*x62 + PAy*x84 - PGy[g]*x85;
                        const T x87 = PAy*x53 - PGy[g]*x56;
                        const T x88 = tpi*x64 + PAy*x85 - PGy[g]*x87;
                        block_20[g] = PBx*x86 - PGx[g]*x88;
                        const T x89 = x84 - x85;
                        block_22[g] = tpi_2*x89 + PBy*x86 - PGy[g]*x88;
                        const T x90 = PAx*x84 - PGx[g]*x85;
                        const T x91 = PAx*x85 - PGx[g]*x87;
                        block_8[g] = PBx*x90 + tpi*x89 - PGx[g]*x91;
                        block_10[g] = tpi*x83 + PBy*x90 - PGy[g]*x91;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x2 = PBx*x0 - PGx[g]*x1;
                        const T x3 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x4 = PBx*x1 - PGx[g]*x3;
                        const T x5 = x2 - x4;
                        const T x6 = PAz*x2 - PGz[g]*x4;
                        const T x7 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x8 = PBx*x3 - PGx[g]*x7;
                        const T x9 = PAz*x4 - PGz[g]*x8;
                        block_31[g] += tpi*x5 + PAz*x6 - PGz[g]*x9;
                        const T x10 = x0 - x1;
                        const T x11 = PAx*x2 + tpi*x10 - PGx[g]*x4;
                        const T x12 = x1 - x3;
                        const T x13 = PAx*x4 + tpi*x12 - PGx[g]*x8;
                        block_13[g] += PAz*x11 - PGz[g]*x13;
                        const T x14 = PAx*x10 - PGx[g]*x12;
                        block_1[g] += PAx*x11 + tpi*x5 + tpi*x14 - PGx[g]*x13;
                        const T x15 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x16 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x17 = x15 - x16;
                        const T x18 = tpi*x17 + PAy*x2 - PGy[g]*x4;
                        const T x19 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x20 = x16 - x19;
                        const T x21 = tpi*x20 + PAy*x4 - PGy[g]*x8;
                        const T x22 = PAy*x17 - PGy[g]*x20;
                        block_19[g] += tpi*x5 + PAy*x18 + tpi*x22 - PGy[g]*x21;
                        block_25[g] += PAz*x18 - PGz[g]*x21;
                        const T x23 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x24 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x25 = x23 - x24;
                        const T x26 = tpi*x25 + PAy*x10 - PGy[g]*x12;
                        block_7[g] += PAx*x18 + tpi*x26 - PGx[g]*x21;
                        const T x27 = tpi*x23 + PBy*x0 - PGy[g]*x1;
                        const T x28 = tpi*x24 + PBy*x1 - PGy[g]*x3;
                        const T x29 = PAx*x27 - PGx[g]*x28;
                        const T x30 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x31 = tpi*x30 + PBy*x3 - PGy[g]*x7;
                        const T x32 = PAx*x28 - PGx[g]*x31;
                        const T x33 = x27 - x28;
                        block_3[g] += PAx*x29 + tpi*x33 - PGx[g]*x32;
                        const T x34 = tpi_2*x10 + PAy*x27 - PGy[g]*x28;
                        const T x35 = tpi_2*x12 + PAy*x28 - PGy[g]*x31;
                        block_21[g] += tpi*x33 + tpi_2*x26 + PAy*x34 - PGy[g]*x35;
                        const T x36 = PAz*x27 - PGz[g]*x28;
                        const T x37 = PAz*x28 - PGz[g]*x31;
                        block_33[g] += tpi*x33 + PAz*x36 - PGz[g]*x37;
                        block_15[g] += PAx*x36 - PGx[g]*x37;
                        block_9[g] += PAx*x34 - PGx[g]*x35;
                        block_27[g] += PAz*x34 - PGz[g]*x35;
                        const T x38 = PBx*x15 + tpi*x23 - PGx[g]*x16;
                        const T x39 = PBx*x16 + tpi*x24 - PGx[g]*x19;
                        const T x40 = x38 - x39;
                        const T x41 = PAy*x38 - PGy[g]*x39;
                        const T x42 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x43 = PBx*x19 + tpi*x30 - PGx[g]*x42;
                        const T x44 = PAy*x39 - PGy[g]*x43;
                        block_18[g] += tpi*x40 + PAy*x41 - PGy[g]*x44;
                        const T x45 = PAx*x38 + tpi_2*x17 - PGx[g]*x39;
                        const T x46 = PAx*x39 + tpi_2*x20 - PGx[g]*x43;
                        const T x47 = PAx*x17 + tpi*x25 - PGx[g]*x20;
                        block_0[g] += PAx*x45 + tpi*x40 + tpi_2*x47 - PGx[g]*x46;
                        const T x48 = PAz*x38 - PGz[g]*x39;
                        const T x49 = PAz*x39 - PGz[g]*x43;
                        block_30[g] += tpi*x40 + PAz*x48 - PGz[g]*x49;
                        block_24[g] += PAy*x48 - PGy[g]*x49;
                        block_6[g] += PAy*x45 - PGy[g]*x46;
                        block_12[g] += PAz*x45 - PGz[g]*x46;
                        const T x50 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x51 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x52 = tpi*x23 + PBz*x50 - PGz[g]*x51;
                        const T x53 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x54 = tpi*x24 + PBz*x51 - PGz[g]*x53;
                        const T x55 = PAx*x52 - PGx[g]*x54;
                        const T x56 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x57 = tpi*x30 + PBz*x53 - PGz[g]*x56;
                        const T x58 = PAx*x54 - PGx[g]*x57;
                        const T x59 = x52 - x54;
                        block_5[g] += PAx*x55 + tpi*x59 - PGx[g]*x58;
                        const T x60 = PAy*x52 - PGy[g]*x54;
                        const T x61 = PAy*x54 - PGy[g]*x57;
                        block_11[g] += PAx*x60 - PGx[g]*x61;
                        block_23[g] += tpi*x59 + PAy*x60 - PGy[g]*x61;
                        const T x62 = x50 - x51;
                        const T x63 = tpi_2*x62 + PAz*x52 - PGz[g]*x54;
                        const T x64 = x51 - x53;
                        const T x65 = tpi_2*x64 + PAz*x54 - PGz[g]*x57;
                        block_17[g] += PAx*x63 - PGx[g]*x65;
                        block_29[g] += PAy*x63 - PGy[g]*x65;
                        const T x66 = tpi*x23 + PAz*x50 - PGz[g]*x51;
                        const T x67 = tpi*x24 + PAz*x51 - PGz[g]*x53;
                        const T x68 = x66 - x67;
                        block_35[g] += tpi*x59 + tpi_2*x68 + PAz*x63 - PGz[g]*x65;
                        const T x69 = PBx*x66 - PGx[g]*x67;
                        const T x70 = tpi*x30 + PAz*x53 - PGz[g]*x56;
                        const T x71 = PBx*x67 - PGx[g]*x70;
                        const T x72 = PBx*x62 - PGx[g]*x64;
                        const T x73 = PAz*x17 - PGz[g]*x20;
                        block_32[g] += PAz*x69 + tpi*x72 + tpi*x73 - PGz[g]*x71;
                        const T x74 = PBy*x66 - PGy[g]*x67;
                        const T x75 = PBy*x67 - PGy[g]*x70;
                        const T x76 = PBy*x62 - PGy[g]*x64;
                        const T x77 = PAz*x10 - PGz[g]*x12;
                        block_34[g] += PAz*x74 + tpi*x76 + tpi*x77 - PGz[g]*x75;
                        block_14[g] += PAx*x69 + tpi*x68 - PGx[g]*x71;
                        block_26[g] += PAy*x69 - PGy[g]*x71;
                        block_28[g] += tpi*x68 + PAy*x74 - PGy[g]*x75;
                        block_16[g] += PAx*x74 - PGx[g]*x75;
                        const T x78 = PAx*x50 - PGx[g]*x51;
                        const T x79 = PAx*x51 - PGx[g]*x53;
                        const T x80 = PAx*x78 + tpi*x62 - PGx[g]*x79;
                        const T x81 = PAx*x53 - PGx[g]*x56;
                        const T x82 = PAx*x79 + tpi*x64 - PGx[g]*x81;
                        const T x83 = x78 - x79;
                        block_2[g] += PBx*x80 + tpi_2*x83 - PGx[g]*x82;
                        block_4[g] += PBy*x80 - PGy[g]*x82;
                        const T x84 = PAy*x50 - PGy[g]*x51;
                        const T x85 = PAy*x51 - PGy[g]*x53;
                        const T x86 = tpi*x62 + PAy*x84 - PGy[g]*x85;
                        const T x87 = PAy*x53 - PGy[g]*x56;
                        const T x88 = tpi*x64 + PAy*x85 - PGy[g]*x87;
                        block_20[g] += PBx*x86 - PGx[g]*x88;
                        const T x89 = x84 - x85;
                        block_22[g] += tpi_2*x89 + PBy*x86 - PGy[g]*x88;
                        const T x90 = PAx*x84 - PGx[g]*x85;
                        const T x91 = PAx*x85 - PGx[g]*x87;
                        block_8[g] += PBx*x90 + tpi*x89 - PGx[g]*x91;
                        block_10[g] += tpi*x83 + PBy*x90 - PGy[g]*x91;

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

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x2 = PBx*x0 - PGx[g]*x1;
                    const T x3 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x4 = PBx*x1 - PGx[g]*x3;
                    const T x5 = x2 - x4;
                    const T x6 = PAz*x2 - PGz[g]*x4;
                    const T x7 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x8 = PBx*x3 - PGx[g]*x7;
                    const T x9 = PAz*x4 - PGz[g]*x8;
                    const T x10 = x0 - x1;
                    const T x11 = PAx*x2 + tpi*x10 - PGx[g]*x4;
                    const T x12 = x1 - x3;
                    const T x13 = PAx*x4 + tpi*x12 - PGx[g]*x8;
                    const T x14 = PAx*x10 - PGx[g]*x12;
                    const T x15 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x16 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x17 = x15 - x16;
                    const T x18 = tpi*x17 + PAy*x2 - PGy[g]*x4;
                    const T x19 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x20 = x16 - x19;
                    const T x21 = tpi*x20 + PAy*x4 - PGy[g]*x8;
                    const T x22 = PAy*x17 - PGy[g]*x20;
                    const T x23 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x24 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x25 = x23 - x24;
                    const T x26 = tpi*x25 + PAy*x10 - PGy[g]*x12;
                    const T x27 = tpi*x23 + PBy*x0 - PGy[g]*x1;
                    const T x28 = tpi*x24 + PBy*x1 - PGy[g]*x3;
                    const T x29 = PAx*x27 - PGx[g]*x28;
                    const T x30 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x31 = tpi*x30 + PBy*x3 - PGy[g]*x7;
                    const T x32 = PAx*x28 - PGx[g]*x31;
                    const T x33 = x27 - x28;
                    const T x34 = tpi_2*x10 + PAy*x27 - PGy[g]*x28;
                    const T x35 = tpi_2*x12 + PAy*x28 - PGy[g]*x31;
                    const T x36 = PAz*x27 - PGz[g]*x28;
                    const T x37 = PAz*x28 - PGz[g]*x31;
                    const T x38 = PBx*x15 + tpi*x23 - PGx[g]*x16;
                    const T x39 = PBx*x16 + tpi*x24 - PGx[g]*x19;
                    const T x40 = x38 - x39;
                    const T x41 = PAy*x38 - PGy[g]*x39;
                    const T x42 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x43 = PBx*x19 + tpi*x30 - PGx[g]*x42;
                    const T x44 = PAy*x39 - PGy[g]*x43;
                    const T x45 = PAx*x38 + tpi_2*x17 - PGx[g]*x39;
                    const T x46 = PAx*x39 + tpi_2*x20 - PGx[g]*x43;
                    const T x47 = PAx*x17 + tpi*x25 - PGx[g]*x20;
                    const T x48 = PAz*x38 - PGz[g]*x39;
                    const T x49 = PAz*x39 - PGz[g]*x43;
                    const T x50 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x51 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x52 = tpi*x23 + PBz*x50 - PGz[g]*x51;
                    const T x53 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x54 = tpi*x24 + PBz*x51 - PGz[g]*x53;
                    const T x55 = PAx*x52 - PGx[g]*x54;
                    const T x56 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x57 = tpi*x30 + PBz*x53 - PGz[g]*x56;
                    const T x58 = PAx*x54 - PGx[g]*x57;
                    const T x59 = x52 - x54;
                    const T x60 = PAy*x52 - PGy[g]*x54;
                    const T x61 = PAy*x54 - PGy[g]*x57;
                    const T x62 = x50 - x51;
                    const T x63 = tpi_2*x62 + PAz*x52 - PGz[g]*x54;
                    const T x64 = x51 - x53;
                    const T x65 = tpi_2*x64 + PAz*x54 - PGz[g]*x57;
                    const T x66 = tpi*x23 + PAz*x50 - PGz[g]*x51;
                    const T x67 = tpi*x24 + PAz*x51 - PGz[g]*x53;
                    const T x68 = x66 - x67;
                    const T x69 = PBx*x66 - PGx[g]*x67;
                    const T x70 = tpi*x30 + PAz*x53 - PGz[g]*x56;
                    const T x71 = PBx*x67 - PGx[g]*x70;
                    const T x72 = PBx*x62 - PGx[g]*x64;
                    const T x73 = PAz*x17 - PGz[g]*x20;
                    const T x74 = PBy*x66 - PGy[g]*x67;
                    const T x75 = PBy*x67 - PGy[g]*x70;
                    const T x76 = PBy*x62 - PGy[g]*x64;
                    const T x77 = PAz*x10 - PGz[g]*x12;
                    const T x78 = PAx*x50 - PGx[g]*x51;
                    const T x79 = PAx*x51 - PGx[g]*x53;
                    const T x80 = PAx*x78 + tpi*x62 - PGx[g]*x79;
                    const T x81 = PAx*x53 - PGx[g]*x56;
                    const T x82 = PAx*x79 + tpi*x64 - PGx[g]*x81;
                    const T x83 = x78 - x79;
                    const T x84 = PAy*x50 - PGy[g]*x51;
                    const T x85 = PAy*x51 - PGy[g]*x53;
                    const T x86 = tpi*x62 + PAy*x84 - PGy[g]*x85;
                    const T x87 = PAy*x53 - PGy[g]*x56;
                    const T x88 = tpi*x64 + PAy*x85 - PGy[g]*x87;
                    const T x89 = x84 - x85;
                    const T x90 = PAx*x84 - PGx[g]*x85;
                    const T x91 = PAx*x85 - PGx[g]*x87;

                    integral_0[g] = PAx*x45 + tpi*x40 + tpi_2*x47 - PGx[g]*x46;
                    integral_1[g] = PAx*x11 + tpi*x5 + tpi*x14 - PGx[g]*x13;
                    integral_2[g] = PBx*x80 + tpi_2*x83 - PGx[g]*x82;
                    integral_3[g] = PAx*x29 + tpi*x33 - PGx[g]*x32;
                    integral_4[g] = PBy*x80 - PGy[g]*x82;
                    integral_5[g] = PAx*x55 + tpi*x59 - PGx[g]*x58;
                    integral_6[g] = PAy*x45 - PGy[g]*x46;
                    integral_7[g] = PAx*x18 + tpi*x26 - PGx[g]*x21;
                    integral_8[g] = PBx*x90 + tpi*x89 - PGx[g]*x91;
                    integral_9[g] = PAx*x34 - PGx[g]*x35;
                    integral_10[g] = tpi*x83 + PBy*x90 - PGy[g]*x91;
                    integral_11[g] = PAx*x60 - PGx[g]*x61;
                    integral_12[g] = PAz*x45 - PGz[g]*x46;
                    integral_13[g] = PAz*x11 - PGz[g]*x13;
                    integral_14[g] = PAx*x69 + tpi*x68 - PGx[g]*x71;
                    integral_15[g] = PAx*x36 - PGx[g]*x37;
                    integral_16[g] = PAx*x74 - PGx[g]*x75;
                    integral_17[g] = PAx*x63 - PGx[g]*x65;
                    integral_18[g] = tpi*x40 + PAy*x41 - PGy[g]*x44;
                    integral_19[g] = tpi*x5 + PAy*x18 + tpi*x22 - PGy[g]*x21;
                    integral_20[g] = PBx*x86 - PGx[g]*x88;
                    integral_21[g] = tpi*x33 + tpi_2*x26 + PAy*x34 - PGy[g]*x35;
                    integral_22[g] = tpi_2*x89 + PBy*x86 - PGy[g]*x88;
                    integral_23[g] = tpi*x59 + PAy*x60 - PGy[g]*x61;
                    integral_24[g] = PAy*x48 - PGy[g]*x49;
                    integral_25[g] = PAz*x18 - PGz[g]*x21;
                    integral_26[g] = PAy*x69 - PGy[g]*x71;
                    integral_27[g] = PAz*x34 - PGz[g]*x35;
                    integral_28[g] = tpi*x68 + PAy*x74 - PGy[g]*x75;
                    integral_29[g] = PAy*x63 - PGy[g]*x65;
                    integral_30[g] = tpi*x40 + PAz*x48 - PGz[g]*x49;
                    integral_31[g] = tpi*x5 + PAz*x6 - PGz[g]*x9;
                    integral_32[g] = PAz*x69 + tpi*x72 + tpi*x73 - PGz[g]*x71;
                    integral_33[g] = tpi*x33 + PAz*x36 - PGz[g]*x37;
                    integral_34[g] = PAz*x74 + tpi*x76 + tpi*x77 - PGz[g]*x75;
                    integral_35[g] = tpi*x59 + tpi_2*x68 + PAz*x63 - PGz[g]*x65;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x2 = PBx*x0 - PGx[g]*x1;
                    const T x3 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x4 = PBx*x1 - PGx[g]*x3;
                    const T x5 = x2 - x4;
                    const T x6 = PAz*x2 - PGz[g]*x4;
                    const T x7 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x8 = PBx*x3 - PGx[g]*x7;
                    const T x9 = PAz*x4 - PGz[g]*x8;
                    const T x10 = x0 - x1;
                    const T x11 = PAx*x2 + tpi*x10 - PGx[g]*x4;
                    const T x12 = x1 - x3;
                    const T x13 = PAx*x4 + tpi*x12 - PGx[g]*x8;
                    const T x14 = PAx*x10 - PGx[g]*x12;
                    const T x15 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x16 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x17 = x15 - x16;
                    const T x18 = tpi*x17 + PAy*x2 - PGy[g]*x4;
                    const T x19 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x20 = x16 - x19;
                    const T x21 = tpi*x20 + PAy*x4 - PGy[g]*x8;
                    const T x22 = PAy*x17 - PGy[g]*x20;
                    const T x23 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x24 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x25 = x23 - x24;
                    const T x26 = tpi*x25 + PAy*x10 - PGy[g]*x12;
                    const T x27 = tpi*x23 + PBy*x0 - PGy[g]*x1;
                    const T x28 = tpi*x24 + PBy*x1 - PGy[g]*x3;
                    const T x29 = PAx*x27 - PGx[g]*x28;
                    const T x30 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x31 = tpi*x30 + PBy*x3 - PGy[g]*x7;
                    const T x32 = PAx*x28 - PGx[g]*x31;
                    const T x33 = x27 - x28;
                    const T x34 = tpi_2*x10 + PAy*x27 - PGy[g]*x28;
                    const T x35 = tpi_2*x12 + PAy*x28 - PGy[g]*x31;
                    const T x36 = PAz*x27 - PGz[g]*x28;
                    const T x37 = PAz*x28 - PGz[g]*x31;
                    const T x38 = PBx*x15 + tpi*x23 - PGx[g]*x16;
                    const T x39 = PBx*x16 + tpi*x24 - PGx[g]*x19;
                    const T x40 = x38 - x39;
                    const T x41 = PAy*x38 - PGy[g]*x39;
                    const T x42 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x43 = PBx*x19 + tpi*x30 - PGx[g]*x42;
                    const T x44 = PAy*x39 - PGy[g]*x43;
                    const T x45 = PAx*x38 + tpi_2*x17 - PGx[g]*x39;
                    const T x46 = PAx*x39 + tpi_2*x20 - PGx[g]*x43;
                    const T x47 = PAx*x17 + tpi*x25 - PGx[g]*x20;
                    const T x48 = PAz*x38 - PGz[g]*x39;
                    const T x49 = PAz*x39 - PGz[g]*x43;
                    const T x50 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x51 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x52 = tpi*x23 + PBz*x50 - PGz[g]*x51;
                    const T x53 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x54 = tpi*x24 + PBz*x51 - PGz[g]*x53;
                    const T x55 = PAx*x52 - PGx[g]*x54;
                    const T x56 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x57 = tpi*x30 + PBz*x53 - PGz[g]*x56;
                    const T x58 = PAx*x54 - PGx[g]*x57;
                    const T x59 = x52 - x54;
                    const T x60 = PAy*x52 - PGy[g]*x54;
                    const T x61 = PAy*x54 - PGy[g]*x57;
                    const T x62 = x50 - x51;
                    const T x63 = tpi_2*x62 + PAz*x52 - PGz[g]*x54;
                    const T x64 = x51 - x53;
                    const T x65 = tpi_2*x64 + PAz*x54 - PGz[g]*x57;
                    const T x66 = tpi*x23 + PAz*x50 - PGz[g]*x51;
                    const T x67 = tpi*x24 + PAz*x51 - PGz[g]*x53;
                    const T x68 = x66 - x67;
                    const T x69 = PBx*x66 - PGx[g]*x67;
                    const T x70 = tpi*x30 + PAz*x53 - PGz[g]*x56;
                    const T x71 = PBx*x67 - PGx[g]*x70;
                    const T x72 = PBx*x62 - PGx[g]*x64;
                    const T x73 = PAz*x17 - PGz[g]*x20;
                    const T x74 = PBy*x66 - PGy[g]*x67;
                    const T x75 = PBy*x67 - PGy[g]*x70;
                    const T x76 = PBy*x62 - PGy[g]*x64;
                    const T x77 = PAz*x10 - PGz[g]*x12;
                    const T x78 = PAx*x50 - PGx[g]*x51;
                    const T x79 = PAx*x51 - PGx[g]*x53;
                    const T x80 = PAx*x78 + tpi*x62 - PGx[g]*x79;
                    const T x81 = PAx*x53 - PGx[g]*x56;
                    const T x82 = PAx*x79 + tpi*x64 - PGx[g]*x81;
                    const T x83 = x78 - x79;
                    const T x84 = PAy*x50 - PGy[g]*x51;
                    const T x85 = PAy*x51 - PGy[g]*x53;
                    const T x86 = tpi*x62 + PAy*x84 - PGy[g]*x85;
                    const T x87 = PAy*x53 - PGy[g]*x56;
                    const T x88 = tpi*x64 + PAy*x85 - PGy[g]*x87;
                    const T x89 = x84 - x85;
                    const T x90 = PAx*x84 - PGx[g]*x85;
                    const T x91 = PAx*x85 - PGx[g]*x87;

                    integral_0[g] += PAx*x45 + tpi*x40 + tpi_2*x47 - PGx[g]*x46;
                    integral_1[g] += PAx*x11 + tpi*x5 + tpi*x14 - PGx[g]*x13;
                    integral_2[g] += PBx*x80 + tpi_2*x83 - PGx[g]*x82;
                    integral_3[g] += PAx*x29 + tpi*x33 - PGx[g]*x32;
                    integral_4[g] += PBy*x80 - PGy[g]*x82;
                    integral_5[g] += PAx*x55 + tpi*x59 - PGx[g]*x58;
                    integral_6[g] += PAy*x45 - PGy[g]*x46;
                    integral_7[g] += PAx*x18 + tpi*x26 - PGx[g]*x21;
                    integral_8[g] += PBx*x90 + tpi*x89 - PGx[g]*x91;
                    integral_9[g] += PAx*x34 - PGx[g]*x35;
                    integral_10[g] += tpi*x83 + PBy*x90 - PGy[g]*x91;
                    integral_11[g] += PAx*x60 - PGx[g]*x61;
                    integral_12[g] += PAz*x45 - PGz[g]*x46;
                    integral_13[g] += PAz*x11 - PGz[g]*x13;
                    integral_14[g] += PAx*x69 + tpi*x68 - PGx[g]*x71;
                    integral_15[g] += PAx*x36 - PGx[g]*x37;
                    integral_16[g] += PAx*x74 - PGx[g]*x75;
                    integral_17[g] += PAx*x63 - PGx[g]*x65;
                    integral_18[g] += tpi*x40 + PAy*x41 - PGy[g]*x44;
                    integral_19[g] += tpi*x5 + PAy*x18 + tpi*x22 - PGy[g]*x21;
                    integral_20[g] += PBx*x86 - PGx[g]*x88;
                    integral_21[g] += tpi*x33 + tpi_2*x26 + PAy*x34 - PGy[g]*x35;
                    integral_22[g] += tpi_2*x89 + PBy*x86 - PGy[g]*x88;
                    integral_23[g] += tpi*x59 + PAy*x60 - PGy[g]*x61;
                    integral_24[g] += PAy*x48 - PGy[g]*x49;
                    integral_25[g] += PAz*x18 - PGz[g]*x21;
                    integral_26[g] += PAy*x69 - PGy[g]*x71;
                    integral_27[g] += PAz*x34 - PGz[g]*x35;
                    integral_28[g] += tpi*x68 + PAy*x74 - PGy[g]*x75;
                    integral_29[g] += PAy*x63 - PGy[g]*x65;
                    integral_30[g] += tpi*x40 + PAz*x48 - PGz[g]*x49;
                    integral_31[g] += tpi*x5 + PAz*x6 - PGz[g]*x9;
                    integral_32[g] += PAz*x69 + tpi*x72 + tpi*x73 - PGz[g]*x71;
                    integral_33[g] += tpi*x33 + PAz*x36 - PGz[g]*x37;
                    integral_34[g] += PAz*x74 + tpi*x76 + tpi*x77 - PGz[g]*x75;
                    integral_35[g] += tpi*x59 + tpi_2*x68 + PAz*x63 - PGz[g]*x65;

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

