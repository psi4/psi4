#include "VRR_D_F.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void VRR_D_F(
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
    if (69*sizeof(T) >= 280) {

        if (nprim < 1) return;

        // The gridpoints are processed in blocks.  The integrals of a block are contracted in
        // a buffer that stays in the L1 cache and are then copied to the output arrays.
        enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
        alignas(64) T block[60*g_block];
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
        T* __restrict block_36 = &block[36*g_block];
        T* __restrict block_37 = &block[37*g_block];
        T* __restrict block_38 = &block[38*g_block];
        T* __restrict block_39 = &block[39*g_block];
        T* __restrict block_40 = &block[40*g_block];
        T* __restrict block_41 = &block[41*g_block];
        T* __restrict block_42 = &block[42*g_block];
        T* __restrict block_43 = &block[43*g_block];
        T* __restrict block_44 = &block[44*g_block];
        T* __restrict block_45 = &block[45*g_block];
        T* __restrict block_46 = &block[46*g_block];
        T* __restrict block_47 = &block[47*g_block];
        T* __restrict block_48 = &block[48*g_block];
        T* __restrict block_49 = &block[49*g_block];
        T* __restrict block_50 = &block[50*g_block];
        T* __restrict block_51 = &block[51*g_block];
        T* __restrict block_52 = &block[52*g_block];
        T* __restrict block_53 = &block[53*g_block];
        T* __restrict block_54 = &block[54*g_block];
        T* __restrict block_55 = &block[55*g_block];
        T* __restrict block_56 = &block[56*g_block];
        T* __restrict block_57 = &block[57*g_block];
        T* __restrict block_58 = &block[58*g_block];
        T* __restrict block_59 = &block[59*g_block];

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

                        const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x2 = PBy*x0 - PGy[g]*x1;
                        const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x4 = PBy*x1 - PGy[g]*x3;
                        const T x5 = PBx*x2 - PGx[g]*x4;
                        const T x6 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x7 = PBy*x3 - PGy[g]*x6;
                        const T x8 = PBx*x4 - PGx[g]*x7;
                        const T x9 = x2 - x4;
                        const T x10 = PAx*x5 + tpi*x9 - PGx[g]*x8;
                        const T x11 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                        const T x12 = PBy*x6 - PGy[g]*x11;
                        const T x13 = PBx*x7 - PGx[g]*x12;
                        const T x14 = x4 - x7;
                        const T x15 = PAx*x8 + tpi*x14 - PGx[g]*x13;
                        const T x16 = x5 - x8;
                        const T x17 = PAx*x9 - PGx[g]*x14;
                        block_4[g] = PAx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                        const T x18 = x0 - x1;
                        const T x19 = x1 - x3;
                        const T x20 = PBx*x18 - PGx[g]*x19;
                        const T x21 = tpi*x20 + PAy*x5 - PGy[g]*x8;
                        const T x22 = x3 - x6;
                        const T x23 = PBx*x19 - PGx[g]*x22;
                        const T x24 = tpi*x23 + PAy*x8 - PGy[g]*x13;
                        const T x25 = PAy*x20 - PGy[g]*x23;
                        block_34[g] = tpi*x16 + PAy*x21 + tpi*x25 - PGy[g]*x24;
                        const T x26 = x18 - x19;
                        const T x27 = tpi*x26 + PAy*x9 - PGy[g]*x14;
                        block_14[g] = PAx*x21 + tpi*x27 - PGx[g]*x24;
                        const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x29 = tpi*x28 + PBz*x0 - PGz[g]*x1;
                        const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x31 = tpi*x30 + PBz*x1 - PGz[g]*x3;
                        const T x32 = tpi_2*x18 + PBz*x29 - PGz[g]*x31;
                        const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x34 = tpi*x33 + PBz*x3 - PGz[g]*x6;
                        const T x35 = tpi_2*x19 + PBz*x31 - PGz[g]*x34;
                        const T x36 = PAx*x32 - PGx[g]*x35;
                        const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                        const T x38 = tpi*x37 + PBz*x6 - PGz[g]*x11;
                        const T x39 = tpi_2*x22 + PBz*x34 - PGz[g]*x38;
                        const T x40 = PAx*x35 - PGx[g]*x39;
                        const T x41 = x32 - x35;
                        block_9[g] = PAx*x36 + tpi*x41 - PGx[g]*x40;
                        const T x42 = PAy*x32 - PGy[g]*x35;
                        const T x43 = PAy*x35 - PGy[g]*x39;
                        block_19[g] = PAx*x42 - PGx[g]*x43;
                        block_39[g] = tpi*x41 + PAy*x42 - PGy[g]*x43;
                        const T x44 = x29 - x31;
                        const T x45 = tpi_3*x44 + PAz*x32 - PGz[g]*x35;
                        const T x46 = x31 - x34;
                        const T x47 = tpi_3*x46 + PAz*x35 - PGz[g]*x39;
                        block_29[g] = PAx*x45 - PGx[g]*x47;
                        block_49[g] = PAy*x45 - PGy[g]*x47;
                        const T x48 = tpi_2*x26 + PAz*x44 - PGz[g]*x46;
                        block_59[g] = tpi*x41 + tpi_3*x48 + PAz*x45 - PGz[g]*x47;
                        const T x49 = PBy*x29 - PGy[g]*x31;
                        const T x50 = PBy*x31 - PGy[g]*x34;
                        const T x51 = PAx*x49 - PGx[g]*x50;
                        const T x52 = PBy*x34 - PGy[g]*x38;
                        const T x53 = PAx*x50 - PGx[g]*x52;
                        const T x54 = x49 - x50;
                        block_8[g] = PAx*x51 + tpi*x54 - PGx[g]*x53;
                        const T x55 = PBx*x29 - PGx[g]*x31;
                        const T x56 = PBx*x31 - PGx[g]*x34;
                        const T x57 = x55 - x56;
                        const T x58 = PAy*x55 - PGy[g]*x56;
                        const T x59 = PBx*x34 - PGx[g]*x38;
                        const T x60 = PAy*x56 - PGy[g]*x59;
                        block_35[g] = tpi*x57 + PAy*x58 - PGy[g]*x60;
                        const T x61 = PAx*x55 + tpi*x44 - PGx[g]*x56;
                        const T x62 = PAx*x56 + tpi*x46 - PGx[g]*x59;
                        block_15[g] = PAy*x61 - PGy[g]*x62;
                        const T x63 = PAx*x44 - PGx[g]*x46;
                        block_5[g] = PAx*x61 + tpi*x57 + tpi*x63 - PGx[g]*x62;
                        const T x64 = tpi_2*x20 + PAz*x55 - PGz[g]*x56;
                        const T x65 = tpi_2*x23 + PAz*x56 - PGz[g]*x59;
                        block_25[g] = PAx*x64 + tpi*x48 - PGx[g]*x65;
                        block_45[g] = PAy*x64 - PGy[g]*x65;
                        const T x66 = tpi*x44 + PAy*x49 - PGy[g]*x50;
                        const T x67 = tpi*x46 + PAy*x50 - PGy[g]*x52;
                        const T x68 = PAy*x44 - PGy[g]*x46;
                        block_38[g] = tpi*x54 + PAy*x66 + tpi*x68 - PGy[g]*x67;
                        block_18[g] = PAx*x66 - PGx[g]*x67;
                        const T x69 = tpi_2*x9 + PAz*x49 - PGz[g]*x50;
                        const T x70 = tpi_2*x14 + PAz*x50 - PGz[g]*x52;
                        block_48[g] = tpi*x48 + PAy*x69 - PGy[g]*x70;
                        block_28[g] = PAx*x69 - PGx[g]*x70;
                        const T x71 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x72 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x73 = x71 - x72;
                        const T x74 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x75 = x72 - x74;
                        const T x76 = x73 - x75;
                        const T x77 = tpi*x76 + PAz*x20 - PGz[g]*x23;
                        block_55[g] = tpi*x57 + tpi_2*x77 + PAz*x64 - PGz[g]*x65;
                        const T x78 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x79 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x80 = x78 - x79;
                        const T x81 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x82 = x79 - x81;
                        const T x83 = x80 - x82;
                        const T x84 = tpi*x83 + PAz*x9 - PGz[g]*x14;
                        block_58[g] = tpi*x54 + tpi_2*x84 + PAz*x69 - PGz[g]*x70;
                        const T x85 = PBx*x80 - PGx[g]*x82;
                        const T x86 = tpi*x85 + PAz*x5 - PGz[g]*x8;
                        const T x87 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x88 = x81 - x87;
                        const T x89 = PBx*x82 - PGx[g]*x88;
                        const T x90 = tpi*x89 + PAz*x8 - PGz[g]*x13;
                        const T x91 = PAz*x85 - PGz[g]*x89;
                        block_54[g] = tpi*x16 + PAz*x86 + tpi*x91 - PGz[g]*x90;
                        block_24[g] = PAx*x86 + tpi*x84 - PGx[g]*x90;
                        block_44[g] = tpi*x77 + PAy*x86 - PGy[g]*x90;
                        const T x92 = tpi*x28 + PBy*x78 - PGy[g]*x79;
                        const T x93 = tpi*x30 + PBy*x79 - PGy[g]*x81;
                        const T x94 = tpi_2*x80 + PBy*x92 - PGy[g]*x93;
                        const T x95 = tpi*x33 + PBy*x81 - PGy[g]*x87;
                        const T x96 = tpi_2*x82 + PBy*x93 - PGy[g]*x95;
                        const T x97 = PAx*x94 - PGx[g]*x96;
                        const T x98 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                        const T x99 = tpi*x37 + PBy*x87 - PGy[g]*x98;
                        const T x100 = tpi_2*x88 + PBy*x95 - PGy[g]*x99;
                        const T x101 = PAx*x96 - PGx[g]*x100;
                        const T x102 = x94 - x96;
                        block_6[g] = PAx*x97 + tpi*x102 - PGx[g]*x101;
                        const T x103 = PBx*x71 + tpi*x28 - PGx[g]*x72;
                        const T x104 = PBx*x72 + tpi*x30 - PGx[g]*x74;
                        const T x105 = PBx*x103 + tpi_2*x73 - PGx[g]*x104;
                        const T x106 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x107 = PBx*x74 + tpi*x33 - PGx[g]*x106;
                        const T x108 = PBx*x104 + tpi_2*x75 - PGx[g]*x107;
                        const T x109 = x105 - x108;
                        const T x110 = PAy*x105 - PGy[g]*x108;
                        const T x111 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                        const T x112 = PBx*x106 + tpi*x37 - PGx[g]*x111;
                        const T x113 = x74 - x106;
                        const T x114 = PBx*x107 + tpi_2*x113 - PGx[g]*x112;
                        const T x115 = PAy*x108 - PGy[g]*x114;
                        block_30[g] = tpi*x109 + PAy*x110 - PGy[g]*x115;
                        const T x116 = PAz*x94 - PGz[g]*x96;
                        const T x117 = PAz*x96 - PGz[g]*x100;
                        block_26[g] = PAx*x116 - PGx[g]*x117;
                        block_56[g] = tpi*x102 + PAz*x116 - PGz[g]*x117;
                        const T x118 = x92 - x93;
                        const T x119 = tpi_3*x118 + PAy*x94 - PGy[g]*x96;
                        const T x120 = x93 - x95;
                        const T x121 = tpi_3*x120 + PAy*x96 - PGy[g]*x100;
                        block_16[g] = PAx*x119 - PGx[g]*x121;
                        block_46[g] = PAz*x119 - PGz[g]*x121;
                        const T x122 = tpi_2*x83 + PAy*x118 - PGy[g]*x120;
                        block_36[g] = tpi*x102 + tpi_3*x122 + PAy*x119 - PGy[g]*x121;
                        const T x123 = PAz*x105 - PGz[g]*x108;
                        const T x124 = PAz*x108 - PGz[g]*x114;
                        block_40[g] = PAy*x123 - PGy[g]*x124;
                        block_50[g] = tpi*x109 + PAz*x123 - PGz[g]*x124;
                        const T x125 = x103 - x104;
                        const T x126 = PAx*x105 + tpi_3*x125 - PGx[g]*x108;
                        const T x127 = x104 - x107;
                        const T x128 = PAx*x108 + tpi_3*x127 - PGx[g]*x114;
                        const T x129 = PAx*x125 + tpi_2*x76 - PGx[g]*x127;
                        block_0[g] = PAx*x126 + tpi*x109 + tpi_3*x129 - PGx[g]*x128;
                        block_10[g] = PAy*x126 - PGy[g]*x128;
                        block_20[g] = PAz*x126 - PGz[g]*x128;
                        const T x130 = PBz*x92 - PGz[g]*x93;
                        const T x131 = PBz*x93 - PGz[g]*x95;
                        const T x132 = tpi_2*x9 + PAy*x130 - PGy[g]*x131;
                        const T x133 = PBz*x95 - PGz[g]*x99;
                        const T x134 = tpi_2*x14 + PAy*x131 - PGy[g]*x133;
                        block_17[g] = PAx*x132 - PGx[g]*x134;
                        const T x135 = x130 - x131;
                        block_37[g] = tpi*x135 + tpi_2*x27 + PAy*x132 - PGy[g]*x134;
                        const T x136 = PAx*x130 - PGx[g]*x131;
                        const T x137 = PAx*x131 - PGx[g]*x133;
                        block_7[g] = PAx*x136 + tpi*x135 - PGx[g]*x137;
                        const T x138 = tpi*x118 + PAz*x130 - PGz[g]*x131;
                        const T x139 = tpi*x120 + PAz*x131 - PGz[g]*x133;
                        block_47[g] = tpi_2*x84 + PAy*x138 - PGy[g]*x139;
                        const T x140 = PAz*x118 - PGz[g]*x120;
                        block_57[g] = tpi*x135 + PAz*x138 + tpi*x140 - PGz[g]*x139;
                        block_27[g] = PAx*x138 - PGx[g]*x139;
                        const T x141 = PBx*x92 - PGx[g]*x93;
                        const T x142 = PBx*x93 - PGx[g]*x95;
                        const T x143 = tpi_2*x85 + PAy*x141 - PGy[g]*x142;
                        const T x144 = PBx*x95 - PGx[g]*x99;
                        const T x145 = tpi_2*x89 + PAy*x142 - PGy[g]*x144;
                        block_13[g] = PAx*x143 + tpi*x122 - PGx[g]*x145;
                        block_43[g] = PAz*x143 - PGz[g]*x145;
                        const T x146 = x141 - x142;
                        const T x147 = tpi*x76 + PAy*x85 - PGy[g]*x89;
                        block_33[g] = tpi*x146 + tpi_2*x147 + PAy*x143 - PGy[g]*x145;
                        const T x148 = PAz*x141 - PGz[g]*x142;
                        const T x149 = PAz*x142 - PGz[g]*x144;
                        block_53[g] = tpi*x146 + PAz*x148 - PGz[g]*x149;
                        const T x150 = PAx*x141 + tpi*x118 - PGx[g]*x142;
                        const T x151 = PAx*x142 + tpi*x120 - PGx[g]*x144;
                        const T x152 = PAx*x118 - PGx[g]*x120;
                        block_3[g] = PAx*x150 + tpi*x146 + tpi*x152 - PGx[g]*x151;
                        block_23[g] = PAz*x150 - PGz[g]*x151;
                        const T x153 = PBy*x103 - PGy[g]*x104;
                        const T x154 = PBy*x104 - PGy[g]*x107;
                        const T x155 = PAx*x153 + tpi_2*x85 - PGx[g]*x154;
                        const T x156 = PBy*x107 - PGy[g]*x112;
                        const T x157 = PAx*x154 + tpi_2*x89 - PGx[g]*x156;
                        const T x158 = x153 - x154;
                        const T x159 = PAx*x85 + tpi*x83 - PGx[g]*x89;
                        block_1[g] = PAx*x155 + tpi*x158 + tpi_2*x159 - PGx[g]*x157;
                        block_21[g] = PAz*x155 - PGz[g]*x157;
                        const T x160 = PBz*x103 - PGz[g]*x104;
                        const T x161 = PBz*x104 - PGz[g]*x107;
                        const T x162 = PAx*x160 + tpi_2*x20 - PGx[g]*x161;
                        const T x163 = PBz*x107 - PGz[g]*x112;
                        const T x164 = PAx*x161 + tpi_2*x23 - PGx[g]*x163;
                        const T x165 = x160 - x161;
                        const T x166 = PAx*x20 + tpi*x26 - PGx[g]*x23;
                        block_2[g] = PAx*x162 + tpi*x165 + tpi_2*x166 - PGx[g]*x164;
                        block_12[g] = PAy*x162 - PGy[g]*x164;
                        const T x167 = PAy*x160 - PGy[g]*x161;
                        const T x168 = PAy*x161 - PGy[g]*x163;
                        block_32[g] = tpi*x165 + PAy*x167 - PGy[g]*x168;
                        const T x169 = tpi*x125 + PAz*x160 - PGz[g]*x161;
                        const T x170 = tpi*x127 + PAz*x161 - PGz[g]*x163;
                        block_22[g] = PAx*x169 + tpi_2*x77 - PGx[g]*x170;
                        const T x171 = PAz*x125 - PGz[g]*x127;
                        block_52[g] = tpi*x165 + PAz*x169 + tpi*x171 - PGz[g]*x170;
                        block_42[g] = PAy*x169 - PGy[g]*x170;
                        const T x172 = PAz*x153 - PGz[g]*x154;
                        const T x173 = PAz*x154 - PGz[g]*x156;
                        block_51[g] = tpi*x158 + PAz*x172 - PGz[g]*x173;
                        const T x174 = tpi*x125 + PAy*x153 - PGy[g]*x154;
                        const T x175 = tpi*x127 + PAy*x154 - PGy[g]*x156;
                        const T x176 = PAy*x125 - PGy[g]*x127;
                        block_31[g] = tpi*x158 + PAy*x174 + tpi*x176 - PGy[g]*x175;
                        block_11[g] = PAx*x174 + tpi_2*x147 - PGx[g]*x175;
                        block_41[g] = PAz*x174 - PGz[g]*x175;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x2 = PBy*x0 - PGy[g]*x1;
                        const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x4 = PBy*x1 - PGy[g]*x3;
                        const T x5 = PBx*x2 - PGx[g]*x4;
                        const T x6 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x7 = PBy*x3 - PGy[g]*x6;
                        const T x8 = PBx*x4 - PGx[g]*x7;
                        const T x9 = x2 - x4;
                        const T x10 = PAx*x5 + tpi*x9 - PGx[g]*x8;
                        const T x11 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                        const T x12 = PBy*x6 - PGy[g]*x11;
                        const T x13 = PBx*x7 - PGx[g]*x12;
                        const T x14 = x4 - x7;
                        const T x15 = PAx*x8 + tpi*x14 - PGx[g]*x13;
                        const T x16 = x5 - x8;
                        const T x17 = PAx*x9 - PGx[g]*x14;
                        block_4[g] += PAx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                        const T x18 = x0 - x1;
                        const T x19 = x1 - x3;
                        const T x20 = PBx*x18 - PGx[g]*x19;
                        const T x21 = tpi*x20 + PAy*x5 - PGy[g]*x8;
                        const T x22 = x3 - x6;
                        const T x23 = PBx*x19 - PGx[g]*x22;
                        const T x24 = tpi*x23 + PAy*x8 - PGy[g]*x13;
                        const T x25 = PAy*x20 - PGy[g]*x23;
                        block_34[g] += tpi*x16 + PAy*x21 + tpi*x25 - PGy[g]*x24;
                        const T x26 = x18 - x19;
                        const T x27 = tpi*x26 + PAy*x9 - PGy[g]*x14;
                        block_14[g] += PAx*x21 + tpi*x27 - PGx[g]*x24;
                        const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x29 = tpi*x28 + PBz*x0 - PGz[g]*x1;
                        const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x31 = tpi*x30 + PBz*x1 - PGz[g]*x3;
                        const T x32 = tpi_2*x18 + PBz*x29 - PGz[g]*x31;
                        const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x34 = tpi*x33 + PBz*x3 - PGz[g]*x6;
                        const T x35 = tpi_2*x19 + PBz*x31 - PGz[g]*x34;
                        const T x36 = PAx*x32 - PGx[g]*x35;
                        const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                        const T x38 = tpi*x37 + PBz*x6 - PGz[g]*x11;
                        const T x39 = tpi_2*x22 + PBz*x34 - PGz[g]*x38;
                        const T x40 = PAx*x35 - PGx[g]*x39;
                        const T x41 = x32 - x35;
                        block_9[g] += PAx*x36 + tpi*x41 - PGx[g]*x40;
                        const T x42 = PAy*x32 - PGy[g]*x35;
                        const T x43 = PAy*x35 - PGy[g]*x39;
                        block_19[g] += PAx*x42 - PGx[g]*x43;
                        block_39[g] += tpi*x41 + PAy*x42 - PGy[g]*x43;
                        const T x44 = x29 - x31;
                        const T x45 = tpi_3*x44 + PAz*x32 - PGz[g]*x35;
                        const T x46 = x31 - x34;
                        const T x47 = tpi_3*x46 + PAz*x35 - PGz[g]*x39;
                        block_29[g] += PAx*x45 - PGx[g]*x47;
                        block_49[g] += PAy*x45 - PGy[g]*x47;
                        const T x48 = tpi_2*x26 + PAz*x44 - PGz[g]*x46;
                        block_59[g] += tpi*x41 + tpi_3*x48 + PAz*x45 - PGz[g]*x47;
                        const T x49 = PBy*x29 - PGy[g]*x31;
                        const T x50 = PBy*x31 - PGy[g]*x34;
                        const T x51 = PAx*x49 - PGx[g]*x50;
                        const T x52 = PBy*x34 - PGy[g]*x38;
                        const T x53 = PAx*x50 - PGx[g]*x52;
                        const T x54 = x49 - x50;
                        block_8[g] += PAx*x51 + tpi*x54 - PGx[g]*x53;
                        const T x55 = PBx*x29 - PGx[g]*x31;
                        const T x56 = PBx*x31 - PGx[g]*x34;
                        const T x57 = x55 - x56;
                        const T x58 = PAy*x55 - PGy[g]*x56;
                        const T x59 = PBx*x34 - PGx[g]*x38;
                        const T x60 = PAy*x56 - PGy[g]*x59;
                        block_35[g] += tpi*x57 + PAy*x58 - PGy[g]*x60;
                        const T x61 = PAx*x55 + tpi*x44 - PGx[g]*x56;
                        const T x62 = PAx*x56 + tpi*x46 - PGx[g]*x59;
                        block_15[g] += PAy*x61 - PGy[g]*x62;
                        const T x63 = PAx*x44 - PGx[g]*x46;
                        block_5[g] += PAx*x61 + tpi*x57 + tpi*x63 - PGx[g]*x62;
                        const T x64 = tpi_2*x20 + PAz*x55 - PGz[g]*x56;
                        const T x65 = tpi_2*x23 + PAz*x56 - PGz[g]*x59;
                        block_25[g] += PAx*x64 + tpi*x48 - PGx[g]*x65;
                        block_45[g] += PAy*x64 - PGy[g]*x65;
                        const T x66 = tpi*x44 + PAy*x49 - PGy[g]*x50;
                        const T x67 = tpi*x46 + PAy*x50 - PGy[g]*x52;
                        const T x68 = PAy*x44 - PGy[g]*x46;
                        block_38[g] += tpi*x54 + PAy*x66 + tpi*x68 - PGy[g]*x67;
                        block_18[g] += PAx*x66 - PGx[g]*x67;
                        const T x69 = tpi_2*x9 + PAz*x49 - PGz[g]*x50;
                        const T x70 = tpi_2*x14 + PAz*x50 - PGz[g]*x52;
                        block_48[g] += tpi*x48 + PAy*x69 - PGy[g]*x70;
                        block_28[g] += PAx*x69 - PGx[g]*x70;
                        const T x71 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x72 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x73 = x71 - x72;
                        const T x74 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x75 = x72 - x74;
                        const T x76 = x73 - x75;
                        const T x77 = tpi*x76 + PAz*x20 - PGz[g]*x23;
                        block_55[g] += tpi*x57 + tpi_2*x77 + PAz*x64 - PGz[g]*x65;
                        const T x78 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x79 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x80 = x78 - x79;
                        const T x81 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x82 = x79 - x81;
                        const T x83 = x80 - x82;
                        const T x84 = tpi*x83 + PAz*x9 - PGz[g]*x14;
                        block_58[g] += tpi*x54 + tpi_2*x84 + PAz*x69 - PGz[g]*x70;
                        const T x85 = PBx*x80 - PGx[g]*x82;
                        const T x86 = tpi*x85 + PAz*x5 - PGz[g]*x8;
                        const T x87 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x88 = x81 - x87;
                        const T x89 = PBx*x82 - PGx[g]*x88;
                        const T x90 = tpi*x89 + PAz*x8 - PGz[g]*x13;
                        const T x91 = PAz*x85 - PGz[g]*x89;
                        block_54[g] += tpi*x16 + PAz*x86 + tpi*x91 - PGz[g]*x90;
                        block_24[g] += PAx*x86 + tpi*x84 - PGx[g]*x90;
                        block_44[g] += tpi*x77 + PAy*x86 - PGy[g]*x90;
                        const T x92 = tpi*x28 + PBy*x78 - PGy[g]*x79;
                        const T x93 = tpi*x30 + PBy*x79 - PGy[g]*x81;
                        const T x94 = tpi_2*x80 + PBy*x92 - PGy[g]*x93;
                        const T x95 = tpi*x33 + PBy*x81 - PGy[g]*x87;
                        const T x96 = tpi_2*x82 + PBy*x93 - PGy[g]*x95;
                        const T x97 = PAx*x94 - PGx[g]*x96;
                        const T x98 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                        const T x99 = tpi*x37 + PBy*x87 - PGy[g]*x98;
                        const T x100 = tpi_2*x88 + PBy*x95 - PGy[g]*x99;
                        const T x101 = PAx*x96 - PGx[g]*x100;
                        const T x102 = x94 - x96;
                        block_6[g] += PAx*x97 + tpi*x102 - PGx[g]*x101;
                        const T x103 = PBx*x71 + tpi*x28 - PGx[g]*x72;
                        const T x104 = PBx*x72 + tpi*x30 - PGx[g]*x74;
                        const T x105 = PBx*x103 + tpi_2*x73 - PGx[g]*x104;
                        const T x106 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x107 = PBx*x74 + tpi*x33 - PGx[g]*x106;
                        const T x108 = PBx*x104 + tpi_2*x75 - PGx[g]*x107;
                        const T x109 = x105 - x108;
                        const T x110 = PAy*x105 - PGy[g]*x108;
                        const T x111 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                        const T x112 = PBx*x106 + tpi*x37 - PGx[g]*x111;
                        const T x113 = x74 - x106;
                        const T x114 = PBx*x107 + tpi_2*x113 - PGx[g]*x112;
                        const T x115 = PAy*x108 - PGy[g]*x114;
                        block_30[g] += tpi*x109 + PAy*x110 - PGy[g]*x115;
                        const T x116 = PAz*x94 - PGz[g]*x96;
                        const T x117 = PAz*x96 - PGz[g]*x100;
                        block_26[g] += PAx*x116 - PGx[g]*x117;
                        block_56[g] += tpi*x102 + PAz*x116 - PGz[g]*x117;
                        const T x118 = x92 - x93;
                        const T x119 = tpi_3*x118 + PAy*x94 - PGy[g]*x96;
                        const T x120 = x93 - x95;
                        const T x121 = tpi_3*x120 + PAy*x96 - PGy[g]*x100;
                        block_16[g] += PAx*x119 - PGx[g]*x121;
                        block_46[g] += PAz*x119 - PGz[g]*x121;
                        const T x122 = tpi_2*x83 + PAy*x118 - PGy[g]*x120;
                        block_36[g] += tpi*x102 + tpi_3*x122 + PAy*x119 - PGy[g]*x121;
                        const T x123 = PAz*x105 - PGz[g]*x108;
                        const T x124 = PAz*x108 - PGz[g]*x114;
                        block_40[g] += PAy*x123 - PGy[g]*x124;
                        block_50[g] += tpi*x109 + PAz*x123 - PGz[g]*x124;
                        const T x125 = x103 - x104;
                        const T x126 = PAx*x105 + tpi_3*x125 - PGx[g]*x108;
                        const T x127 = x104 - x107;
                        const T x128 = PAx*x108 + tpi_3*x127 - PGx[g]*x114;
                        const T x129 = PAx*x125 + tpi_2*x76 - PGx[g]*x127;
                        block_0[g] += PAx*x126 + tpi*x109 + tpi_3*x129 - PGx[g]*x128;
                        block_10[g] += PAy*x126 - PGy[g]*x128;
                        block_20[g] += PAz*x126 - PGz[g]*x128;
                        const T x130 = PBz*x92 - PGz[g]*x93;
                        const T x131 = PBz*x93 - PGz[g]*x95;
                        const T x132 = tpi_2*x9 + PAy*x130 - PGy[g]*x131;
                        const T x133 = PBz*x95 - PGz[g]*x99;
                        const T x134 = tpi_2*x14 + PAy*x131 - PGy[g]*x133;
                        block_17[g] += PAx*x132 - PGx[g]*x134;
                        const T x135 = x130 - x131;
                        block_37[g] += tpi*x135 + tpi_2*x27 + PAy*x132 - PGy[g]*x134;
                        const T x136 = PAx*x130 - PGx[g]*x131;
                        const T x137 = PAx*x131 - PGx[g]*x133;
                        block_7[g] += PAx*x136 + tpi*x135 - PGx[g]*x137;
                        const T x138 = tpi*x118 + PAz*x130 - PGz[g]*x131;
                        const T x139 = tpi*x120 + PAz*x131 - PGz[g]*x133;
                        block_47[g] += tpi_2*x84 + PAy*x138 - PGy[g]*x139;
                        const T x140 = PAz*x118 - PGz[g]*x120;
                        block_57[g] += tpi*x135 + PAz*x138 + tpi*x140 - PGz[g]*x139;
                        block_27[g] += PAx*x138 - PGx[g]*x139;
                        const T x141 = PBx*x92 - PGx[g]*x93;
                        const T x142 = PBx*x93 - PGx[g]*x95;
                        const T x143 = tpi_2*x85 + PAy*x141 - PGy[g]*x142;
                        const T x144 = PBx*x95 - PGx[g]*x99;
                        const T x145 = tpi_2*x89 + PAy*x142 - PGy[g]*x144;
                        block_13[g] += PAx*x143 + tpi*x122 - PGx[g]*x145;
                        block_43[g] += PAz*x143 - PGz[g]*x145;
                        const T x146 = x141 - x142;
                        const T x147 = tpi*x76 + PAy*x85 - PGy[g]*x89;
                        block_33[g] += tpi*x146 + tpi_2*x147 + PAy*x143 - PGy[g]*x145;
                        const T x148 = PAz*x141 - PGz[g]*x142;
                        const T x149 = PAz*x142 - PGz[g]*x144;
                        block_53[g] += tpi*x146 + PAz*x148 - PGz[g]*x149;
                        const T x150 = PAx*x141 + tpi*x118 - PGx[g]*x142;
                        const T x151 = PAx*x142 + tpi*x120 - PGx[g]*x144;
                        const T x152 = PAx*x118 - PGx[g]*x120;
                        block_3[g] += PAx*x150 + tpi*x146 + tpi*x152 - PGx[g]*x151;
                        block_23[g] += PAz*x150 - PGz[g]*x151;
                        const T x153 = PBy*x103 - PGy[g]*x104;
                        const T x154 = PBy*x104 - PGy[g]*x107;
                        const T x155 = PAx*x153 + tpi_2*x85 - PGx[g]*x154;
                        const T x156 = PBy*x107 - PGy[g]*x112;
                        const T x157 = PAx*x154 + tpi_2*x89 - PGx[g]*x156;
                        const T x158 = x153 - x154;
                        const T x159 = PAx*x85 + tpi*x83 - PGx[g]*x89;
                        block_1[g] += PAx*x155 + tpi*x158 + tpi_2*x159 - PGx[g]*x157;
                        block_21[g] += PAz*x155 - PGz[g]*x157;
                        const T x160 = PBz*x103 - PGz[g]*x104;
                        const T x161 = PBz*x104 - PGz[g]*x107;
                        const T x162 = PAx*x160 + tpi_2*x20 - PGx[g]*x161;
                        const T x163 = PBz*x107 - PGz[g]*x112;
                        const T x164 = PAx*x161 + tpi_2*x23 - PGx[g]*x163;
                        const T x165 = x160 - x161;
                        const T x166 = PAx*x20 + tpi*x26 - PGx[g]*x23;
                        block_2[g] += PAx*x162 + tpi*x165 + tpi_2*x166 - PGx[g]*x164;
                        block_12[g] += PAy*x162 - PGy[g]*x164;
                        const T x167 = PAy*x160 - PGy[g]*x161;
                        const T x168 = PAy*x161 - PGy[g]*x163;
                        block_32[g] += tpi*x165 + PAy*x167 - PGy[g]*x168;
                        const T x169 = tpi*x125 + PAz*x160 - PGz[g]*x161;
                        const T x170 = tpi*x127 + PAz*x161 - PGz[g]*x163;
                        block_22[g] += PAx*x169 + tpi_2*x77 - PGx[g]*x170;
                        const T x171 = PAz*x125 - PGz[g]*x127;
                        block_52[g] += tpi*x165 + PAz*x169 + tpi*x171 - PGz[g]*x170;
                        block_42[g] += PAy*x169 - PGy[g]*x170;
                        const T x172 = PAz*x153 - PGz[g]*x154;
                        const T x173 = PAz*x154 - PGz[g]*x156;
                        block_51[g] += tpi*x158 + PAz*x172 - PGz[g]*x173;
                        const T x174 = tpi*x125 + PAy*x153 - PGy[g]*x154;
                        const T x175 = tpi*x127 + PAy*x154 - PGy[g]*x156;
                        const T x176 = PAy*x125 - PGy[g]*x127;
                        block_31[g] += tpi*x158 + PAy*x174 + tpi*x176 - PGy[g]*x175;
                        block_11[g] += PAx*x174 + tpi_2*x147 - PGx[g]*x175;
                        block_41[g] += PAz*x174 - PGz[g]*x175;

                    } // for g
                } // if/else prim
            } // for prim

            // copy the block to the output arrays
            for (int i = 0; i < 60; ++i) {
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

                    const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PBy*x0 - PGy[g]*x1;
                    const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PBy*x1 - PGy[g]*x3;
                    const T x5 = PBx*x2 - PGx[g]*x4;
                    const T x6 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x7 = PBy*x3 - PGy[g]*x6;
                    const T x8 = PBx*x4 - PGx[g]*x7;
                    const T x9 = x2 - x4;
                    const T x10 = PAx*x5 + tpi*x9 - PGx[g]*x8;
                    const T x11 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x12 = PBy*x6 - PGy[g]*x11;
                    const T x13 = PBx*x7 - PGx[g]*x12;
                    const T x14 = x4 - x7;
                    const T x15 = PAx*x8 + tpi*x14 - PGx[g]*x13;
                    const T x16 = x5 - x8;
                    const T x17 = PAx*x9 - PGx[g]*x14;
                    const T x18 = x0 - x1;
                    const T x19 = x1 - x3;
                    const T x20 = PBx*x18 - PGx[g]*x19;
                    const T x21 = tpi*x20 + PAy*x5 - PGy[g]*x8;
                    const T x22 = x3 - x6;
                    const T x23 = PBx*x19 - PGx[g]*x22;
                    const T x24 = tpi*x23 + PAy*x8 - PGy[g]*x13;
                    const T x25 = PAy*x20 - PGy[g]*x23;
                    const T x26 = x18 - x19;
                    const T x27 = tpi*x26 + PAy*x9 - PGy[g]*x14;
                    const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x29 = tpi*x28 + PBz*x0 - PGz[g]*x1;
                    const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x31 = tpi*x30 + PBz*x1 - PGz[g]*x3;
                    const T x32 = tpi_2*x18 + PBz*x29 - PGz[g]*x31;
                    const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x34 = tpi*x33 + PBz*x3 - PGz[g]*x6;
                    const T x35 = tpi_2*x19 + PBz*x31 - PGz[g]*x34;
                    const T x36 = PAx*x32 - PGx[g]*x35;
                    const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x38 = tpi*x37 + PBz*x6 - PGz[g]*x11;
                    const T x39 = tpi_2*x22 + PBz*x34 - PGz[g]*x38;
                    const T x40 = PAx*x35 - PGx[g]*x39;
                    const T x41 = x32 - x35;
                    const T x42 = PAy*x32 - PGy[g]*x35;
                    const T x43 = PAy*x35 - PGy[g]*x39;
                    const T x44 = x29 - x31;
                    const T x45 = tpi_3*x44 + PAz*x32 - PGz[g]*x35;
                    const T x46 = x31 - x34;
                    const T x47 = tpi_3*x46 + PAz*x35 - PGz[g]*x39;
                    const T x48 = tpi_2*x26 + PAz*x44 - PGz[g]*x46;
                    const T x49 = PBy*x29 - PGy[g]*x31;
                    const T x50 = PBy*x31 - PGy[g]*x34;
                    const T x51 = PAx*x49 - PGx[g]*x50;
                    const T x52 = PBy*x34 - PGy[g]*x38;
                    const T x53 = PAx*x50 - PGx[g]*x52;
                    const T x54 = x49 - x50;
                    const T x55 = PBx*x29 - PGx[g]*x31;
                    const T x56 = PBx*x31 - PGx[g]*x34;
                    const T x57 = x55 - x56;
                    const T x58 = PAy*x55 - PGy[g]*x56;
                    const T x59 = PBx*x34 - PGx[g]*x38;
                    const T x60 = PAy*x56 - PGy[g]*x59;
                    const T x61 = PAx*x55 + tpi*x44 - PGx[g]*x56;
                    const T x62 = PAx*x56 + tpi*x46 - PGx[g]*x59;
                    const T x63 = PAx*x44 - PGx[g]*x46;
                    const T x64 = tpi_2*x20 + PAz*x55 - PGz[g]*x56;
                    const T x65 = tpi_2*x23 + PAz*x56 - PGz[g]*x59;
                    const T x66 = tpi*x44 + PAy*x49 - PGy[g]*x50;
                    const T x67 = tpi*x46 + PAy*x50 - PGy[g]*x52;
                    const T x68 = PAy*x44 - PGy[g]*x46;
                    const T x69 = tpi_2*x9 + PAz*x49 - PGz[g]*x50;
                    const T x70 = tpi_2*x14 + PAz*x50 - PGz[g]*x52;
                    const T x71 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x72 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x73 = x71 - x72;
                    const T x74 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x75 = x72 - x74;
                    const T x76 = x73 - x75;
                    const T x77 = tpi*x76 + PAz*x20 - PGz[g]*x23;
                    const T x78 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x79 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x80 = x78 - x79;
                    const T x81 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x82 = x79 - x81;
                    const T x83 = x80 - x82;
                    const T x84 = tpi*x83 + PAz*x9 - PGz[g]*x14;
                    const T x85 = PBx*x80 - PGx[g]*x82;
                    const T x86 = tpi*x85 + PAz*x5 - PGz[g]*x8;
                    const T x87 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x88 = x81 - x87;
                    const T x89 = PBx*x82 - PGx[g]*x88;
                    const T x90 = tpi*x89 + PAz*x8 - PGz[g]*x13;
                    const T x91 = PAz*x85 - PGz[g]*x89;
                    const T x92 = tpi*x28 + PBy*x78 - PGy[g]*x79;
                    const T x93 = tpi*x30 + PBy*x79 - PGy[g]*x81;
                    const T x94 = tpi_2*x80 + PBy*x92 - PGy[g]*x93;
                    const T x95 = tpi*x33 + PBy*x81 - PGy[g]*x87;
                    const T x96 = tpi_2*x82 + PBy*x93 - PGy[g]*x95;
                    const T x97 = PAx*x94 - PGx[g]*x96;
                    const T x98 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x99 = tpi*x37 + PBy*x87 - PGy[g]*x98;
                    const T x100 = tpi_2*x88 + PBy*x95 - PGy[g]*x99;
                    const T x101 = PAx*x96 - PGx[g]*x100;
                    const T x102 = x94 - x96;
                    const T x103 = PBx*x71 + tpi*x28 - PGx[g]*x72;
                    const T x104 = PBx*x72 + tpi*x30 - PGx[g]*x74;
                    const T x105 = PBx*x103 + tpi_2*x73 - PGx[g]*x104;
                    const T x106 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x107 = PBx*x74 + tpi*x33 - PGx[g]*x106;
                    const T x108 = PBx*x104 + tpi_2*x75 - PGx[g]*x107;
                    const T x109 = x105 - x108;
                    const T x110 = PAy*x105 - PGy[g]*x108;
                    const T x111 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x112 = PBx*x106 + tpi*x37 - PGx[g]*x111;
                    const T x113 = x74 - x106;
                    const T x114 = PBx*x107 + tpi_2*x113 - PGx[g]*x112;
                    const T x115 = PAy*x108 - PGy[g]*x114;
                    const T x116 = PAz*x94 - PGz[g]*x96;
                    const T x117 = PAz*x96 - PGz[g]*x100;
                    const T x118 = x92 - x93;
                    const T x119 = tpi_3*x118 + PAy*x94 - PGy[g]*x96;
                    const T x120 = x93 - x95;
                    const T x121 = tpi_3*x120 + PAy*x96 - PGy[g]*x100;
                    const T x122 = tpi_2*x83 + PAy*x118 - PGy[g]*x120;
                    const T x123 = PAz*x105 - PGz[g]*x108;
                    const T x124 = PAz*x108 - PGz[g]*x114;
                    const T x125 = x103 - x104;
                    const T x126 = PAx*x105 + tpi_3*x125 - PGx[g]*x108;
                    const T x127 = x104 - x107;
                    const T x128 = PAx*x108 + tpi_3*x127 - PGx[g]*x114;
                    const T x129 = PAx*x125 + tpi_2*x76 - PGx[g]*x127;
                    const T x130 = PBz*x92 - PGz[g]*x93;
                    const T x131 = PBz*x93 - PGz[g]*x95;
                    const T x132 = tpi_2*x9 + PAy*x130 - PGy[g]*x131;
                    const T x133 = PBz*x95 - PGz[g]*x99;
                    const T x134 = tpi_2*x14 + PAy*x131 - PGy[g]*x133;
                    const T x135 = x130 - x131;
                    const T x136 = PAx*x130 - PGx[g]*x131;
                    const T x137 = PAx*x131 - PGx[g]*x133;
                    const T x138 = tpi*x118 + PAz*x130 - PGz[g]*x131;
                    const T x139 = tpi*x120 + PAz*x131 - PGz[g]*x133;
                    const T x140 = PAz*x118 - PGz[g]*x120;
                    const T x141 = PBx*x92 - PGx[g]*x93;
                    const T x142 = PBx*x93 - PGx[g]*x95;
                    const T x143 = tpi_2*x85 + PAy*x141 - PGy[g]*x142;
                    const T x144 = PBx*x95 - PGx[g]*x99;
                    const T x145 = tpi_2*x89 + PAy*x142 - PGy[g]*x144;
                    const T x146 = x141 - x142;
                    const T x147 = tpi*x76 + PAy*x85 - PGy[g]*x89;
                    const T x148 = PAz*x141 - PGz[g]*x142;
                    const T x149 = PAz*x142 - PGz[g]*x144;
                    const T x150 = PAx*x141 + tpi*x118 - PGx[g]*x142;
                    const T x151 = PAx*x142 + tpi*x120 - PGx[g]*x144;
                    const T x152 = PAx*x118 - PGx[g]*x120;
                    const T x153 = PBy*x103 - PGy[g]*x104;
                    const T x154 = PBy*x104 - PGy[g]*x107;
                    const T x155 = PAx*x153 + tpi_2*x85 - PGx[g]*x154;
                    const T x156 = PBy*x107 - PGy[g]*x112;
                    const T x157 = PAx*x154 + tpi_2*x89 - PGx[g]*x156;
                    const T x158 = x153 - x154;
                    const T x159 = PAx*x85 + tpi*x83 - PGx[g]*x89;
                    const T x160 = PBz*x103 - PGz[g]*x104;
                    const T x161 = PBz*x104 - PGz[g]*x107;
                    const T x162 = PAx*x160 + tpi_2*x20 - PGx[g]*x161;
                    const T x163 = PBz*x107 - PGz[g]*x112;
                    const T x164 = PAx*x161 + tpi_2*x23 - PGx[g]*x163;
                    const T x165 = x160 - x161;
                    const T x166 = PAx*x20 + tpi*x26 - PGx[g]*x23;
                    const T x167 = PAy*x160 - PGy[g]*x161;
                    const T x168 = PAy*x161 - PGy[g]*x163;
                    const T x169 = tpi*x125 + PAz*x160 - PGz[g]*x161;
                    const T x170 = tpi*x127 + PAz*x161 - PGz[g]*x163;
                    const T x171 = PAz*x125 - PGz[g]*x127;
                    const T x172 = PAz*x153 - PGz[g]*x154;
                    const T x173 = PAz*x154 - PGz[g]*x156;
                    const T x174 = tpi*x125 + PAy*x153 - PGy[g]*x154;
                    const T x175 = tpi*x127 + PAy*x154 - PGy[g]*x156;
                    const T x176 = PAy*x125 - PGy[g]*x127;

                    integral_0[g] = PAx*x126 + tpi*x109 + tpi_3*x129 - PGx[g]*x128;
                    integral_1[g] = PAx*x155 + tpi*x158 + tpi_2*x159 - PGx[g]*x157;
                    integral_2[g] = PAx*x162 + tpi*x165 + tpi_2*x166 - PGx[g]*x164;
                    integral_3[g] = PAx*x150 + tpi*x146 + tpi*x152 - PGx[g]*x151;
                    integral_4[g] = PAx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                    integral_5[g] = PAx*x61 + tpi*x57 + tpi*x63 - PGx[g]*x62;
                    integral_6[g] = PAx*x97 + tpi*x102 - PGx[g]*x101;
                    integral_7[g] = PAx*x136 + tpi*x135 - PGx[g]*x137;
                    integral_8[g] = PAx*x51 + tpi*x54 - PGx[g]*x53;
                    integral_9[g] = PAx*x36 + tpi*x41 - PGx[g]*x40;
                    integral_10[g] = PAy*x126 - PGy[g]*x128;
                    integral_11[g] = PAx*x174 + tpi_2*x147 - PGx[g]*x175;
                    integral_12[g] = PAy*x162 - PGy[g]*x164;
                    integral_13[g] = PAx*x143 + tpi*x122 - PGx[g]*x145;
                    integral_14[g] = PAx*x21 + tpi*x27 - PGx[g]*x24;
                    integral_15[g] = PAy*x61 - PGy[g]*x62;
                    integral_16[g] = PAx*x119 - PGx[g]*x121;
                    integral_17[g] = PAx*x132 - PGx[g]*x134;
                    integral_18[g] = PAx*x66 - PGx[g]*x67;
                    integral_19[g] = PAx*x42 - PGx[g]*x43;
                    integral_20[g] = PAz*x126 - PGz[g]*x128;
                    integral_21[g] = PAz*x155 - PGz[g]*x157;
                    integral_22[g] = PAx*x169 + tpi_2*x77 - PGx[g]*x170;
                    integral_23[g] = PAz*x150 - PGz[g]*x151;
                    integral_24[g] = PAx*x86 + tpi*x84 - PGx[g]*x90;
                    integral_25[g] = PAx*x64 + tpi*x48 - PGx[g]*x65;
                    integral_26[g] = PAx*x116 - PGx[g]*x117;
                    integral_27[g] = PAx*x138 - PGx[g]*x139;
                    integral_28[g] = PAx*x69 - PGx[g]*x70;
                    integral_29[g] = PAx*x45 - PGx[g]*x47;
                    integral_30[g] = tpi*x109 + PAy*x110 - PGy[g]*x115;
                    integral_31[g] = tpi*x158 + PAy*x174 + tpi*x176 - PGy[g]*x175;
                    integral_32[g] = tpi*x165 + PAy*x167 - PGy[g]*x168;
                    integral_33[g] = tpi*x146 + tpi_2*x147 + PAy*x143 - PGy[g]*x145;
                    integral_34[g] = tpi*x16 + PAy*x21 + tpi*x25 - PGy[g]*x24;
                    integral_35[g] = tpi*x57 + PAy*x58 - PGy[g]*x60;
                    integral_36[g] = tpi*x102 + tpi_3*x122 + PAy*x119 - PGy[g]*x121;
                    integral_37[g] = tpi*x135 + tpi_2*x27 + PAy*x132 - PGy[g]*x134;
                    integral_38[g] = tpi*x54 + PAy*x66 + tpi*x68 - PGy[g]*x67;
                    integral_39[g] = tpi*x41 + PAy*x42 - PGy[g]*x43;
                    integral_40[g] = PAy*x123 - PGy[g]*x124;
                    integral_41[g] = PAz*x174 - PGz[g]*x175;
                    integral_42[g] = PAy*x169 - PGy[g]*x170;
                    integral_43[g] = PAz*x143 - PGz[g]*x145;
                    integral_44[g] = tpi*x77 + PAy*x86 - PGy[g]*x90;
                    integral_45[g] = PAy*x64 - PGy[g]*x65;
                    integral_46[g] = PAz*x119 - PGz[g]*x121;
                    integral_47[g] = tpi_2*x84 + PAy*x138 - PGy[g]*x139;
                    integral_48[g] = tpi*x48 + PAy*x69 - PGy[g]*x70;
                    integral_49[g] = PAy*x45 - PGy[g]*x47;
                    integral_50[g] = tpi*x109 + PAz*x123 - PGz[g]*x124;
                    integral_51[g] = tpi*x158 + PAz*x172 - PGz[g]*x173;
                    integral_52[g] = tpi*x165 + PAz*x169 + tpi*x171 - PGz[g]*x170;
                    integral_53[g] = tpi*x146 + PAz*x148 - PGz[g]*x149;
                    integral_54[g] = tpi*x16 + PAz*x86 + tpi*x91 - PGz[g]*x90;
                    integral_55[g] = tpi*x57 + tpi_2*x77 + PAz*x64 - PGz[g]*x65;
                    integral_56[g] = tpi*x102 + PAz*x116 - PGz[g]*x117;
                    integral_57[g] = tpi*x135 + PAz*x138 + tpi*x140 - PGz[g]*x139;
                    integral_58[g] = tpi*x54 + tpi_2*x84 + PAz*x69 - PGz[g]*x70;
                    integral_59[g] = tpi*x41 + tpi_3*x48 + PAz*x45 - PGz[g]*x47;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PBy*x0 - PGy[g]*x1;
                    const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PBy*x1 - PGy[g]*x3;
                    const T x5 = PBx*x2 - PGx[g]*x4;
                    const T x6 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x7 = PBy*x3 - PGy[g]*x6;
                    const T x8 = PBx*x4 - PGx[g]*x7;
                    const T x9 = x2 - x4;
                    const T x10 = PAx*x5 + tpi*x9 - PGx[g]*x8;
                    const T x11 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x12 = PBy*x6 - PGy[g]*x11;
                    const T x13 = PBx*x7 - PGx[g]*x12;
                    const T x14 = x4 - x7;
                    const T x15 = PAx*x8 + tpi*x14 - PGx[g]*x13;
                    const T x16 = x5 - x8;
                    const T x17 = PAx*x9 - PGx[g]*x14;
                    const T x18 = x0 - x1;
                    const T x19 = x1 - x3;
                    const T x20 = PBx*x18 - PGx[g]*x19;
                    const T x21 = tpi*x20 + PAy*x5 - PGy[g]*x8;
                    const T x22 = x3 - x6;
                    const T x23 = PBx*x19 - PGx[g]*x22;
                    const T x24 = tpi*x23 + PAy*x8 - PGy[g]*x13;
                    const T x25 = PAy*x20 - PGy[g]*x23;
                    const T x26 = x18 - x19;
                    const T x27 = tpi*x26 + PAy*x9 - PGy[g]*x14;
                    const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x29 = tpi*x28 + PBz*x0 - PGz[g]*x1;
                    const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x31 = tpi*x30 + PBz*x1 - PGz[g]*x3;
                    const T x32 = tpi_2*x18 + PBz*x29 - PGz[g]*x31;
                    const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x34 = tpi*x33 + PBz*x3 - PGz[g]*x6;
                    const T x35 = tpi_2*x19 + PBz*x31 - PGz[g]*x34;
                    const T x36 = PAx*x32 - PGx[g]*x35;
                    const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x38 = tpi*x37 + PBz*x6 - PGz[g]*x11;
                    const T x39 = tpi_2*x22 + PBz*x34 - PGz[g]*x38;
                    const T x40 = PAx*x35 - PGx[g]*x39;
                    const T x41 = x32 - x35;
                    const T x42 = PAy*x32 - PGy[g]*x35;
                    const T x43 = PAy*x35 - PGy[g]*x39;
                    const T x44 = x29 - x31;
                    const T x45 = tpi_3*x44 + PAz*x32 - PGz[g]*x35;
                    const T x46 = x31 - x34;
                    const T x47 = tpi_3*x46 + PAz*x35 - PGz[g]*x39;
                    const T x48 = tpi_2*x26 + PAz*x44 - PGz[g]*x46;
                    const T x49 = PBy*x29 - PGy[g]*x31;
                    const T x50 = PBy*x31 - PGy[g]*x34;
                    const T x51 = PAx*x49 - PGx[g]*x50;
                    const T x52 = PBy*x34 - PGy[g]*x38;
                    const T x53 = PAx*x50 - PGx[g]*x52;
                    const T x54 = x49 - x50;
                    const T x55 = PBx*x29 - PGx[g]*x31;
                    const T x56 = PBx*x31 - PGx[g]*x34;
                    const T x57 = x55 - x56;
                    const T x58 = PAy*x55 - PGy[g]*x56;
                    const T x59 = PBx*x34 - PGx[g]*x38;
                    const T x60 = PAy*x56 - PGy[g]*x59;
                    const T x61 = PAx*x55 + tpi*x44 - PGx[g]*x56;
                    const T x62 = PAx*x56 + tpi*x46 - PGx[g]*x59;
                    const T x63 = PAx*x44 - PGx[g]*x46;
                    const T x64 = tpi_2*x20 + PAz*x55 - PGz[g]*x56;
                    const T x65 = tpi_2*x23 + PAz*x56 - PGz[g]*x59;
                    const T x66 = tpi*x44 + PAy*x49 - PGy[g]*x50;
                    const T x67 = tpi*x46 + PAy*x50 - PGy[g]*x52;
                    const T x68 = PAy*x44 - PGy[g]*x46;
                    const T x69 = tpi_2*x9 + PAz*x49 - PGz[g]*x50;
                    const T x70 = tpi_2*x14 + PAz*x50 - PGz[g]*x52;
                    const T x71 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x72 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x73 = x71 - x72;
                    const T x74 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x75 = x72 - x74;
                    const T x76 = x73 - x75;
                    const T x77 = tpi*x76 + PAz*x20 - PGz[g]*x23;
                    const T x78 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x79 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x80 = x78 - x79;
                    const T x81 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x82 = x79 - x81;
                    const T x83 = x80 - x82;
                    const T x84 = tpi*x83 + PAz*x9 - PGz[g]*x14;
                    const T x85 = PBx*x80 - PGx[g]*x82;
                    const T x86 = tpi*x85 + PAz*x5 - PGz[g]*x8;
                    const T x87 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x88 = x81 - x87;
                    const T x89 = PBx*x82 - PGx[g]*x88;
                    const T x90 = tpi*x89 + PAz*x8 - PGz[g]*x13;
                    const T x91 = PAz*x85 - PGz[g]*x89;
                    const T x92 = tpi*x28 + PBy*x78 - PGy[g]*x79;
                    const T x93 = tpi*x30 + PBy*x79 - PGy[g]*x81;
                    const T x94 = tpi_2*x80 + PBy*x92 - PGy[g]*x93;
                    const T x95 = tpi*x33 + PBy*x81 - PGy[g]*x87;
                    const T x96 = tpi_2*x82 + PBy*x93 - PGy[g]*x95;
                    const T x97 = PAx*x94 - PGx[g]*x96;
                    const T x98 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x99 = tpi*x37 + PBy*x87 - PGy[g]*x98;
                    const T x100 = tpi_2*x88 + PBy*x95 - PGy[g]*x99;
                    const T x101 = PAx*x96 - PGx[g]*x100;
                    const T x102 = x94 - x96;
                    const T x103 = PBx*x71 + tpi*x28 - PGx[g]*x72;
                    const T x104 = PBx*x72 + tpi*x30 - PGx[g]*x74;
                    const T x105 = PBx*x103 + tpi_2*x73 - PGx[g]*x104;
                    const T x106 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x107 = PBx*x74 + tpi*x33 - PGx[g]*x106;
                    const T x108 = PBx*x104 + tpi_2*x75 - PGx[g]*x107;
                    const T x109 = x105 - x108;
                    const T x110 = PAy*x105 - PGy[g]*x108;
                    const T x111 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x112 = PBx*x106 + tpi*x37 - PGx[g]*x111;
                    const T x113 = x74 - x106;
                    const T x114 = PBx*x107 + tpi_2*x113 - PGx[g]*x112;
                    const T x115 = PAy*x108 - PGy[g]*x114;
                    const T x116 = PAz*x94 - PGz[g]*x96;
                    const T x117 = PAz*x96 - PGz[g]*x100;
                    const T x118 = x92 - x93;
                    const T x119 = tpi_3*x118 + PAy*x94 - PGy[g]*x96;
                    const T x120 = x93 - x95;
                    const T x121 = tpi_3*x120 + PAy*x96 - PGy[g]*x100;
                    const T x122 = tpi_2*x83 + PAy*x118 - PGy[g]*x120;
                    const T x123 = PAz*x105 - PGz[g]*x108;
                    const T x124 = PAz*x108 - PGz[g]*x114;
                    const T x125 = x103 - x104;
                    const T x126 = PAx*x105 + tpi_3*x125 - PGx[g]*x108;
                    const T x127 = x104 - x107;
                    const T x128 = PAx*x108 + tpi_3*x127 - PGx[g]*x114;
                    const T x129 = PAx*x125 + tpi_2*x76 - PGx[g]*x127;
                    const T x130 = PBz*x92 - PGz[g]*x93;
                    const T x131 = PBz*x93 - PGz[g]*x95;
                    const T x132 = tpi_2*x9 + PAy*x130 - PGy[g]*x131;
                    const T x133 = PBz*x95 - PGz[g]*x99;
                    const T x134 = tpi_2*x14 + PAy*x131 - PGy[g]*x133;
                    const T x135 = x130 - x131;
                    const T x136 = PAx*x130 - PGx[g]*x131;
                    const T x137 = PAx*x131 - PGx[g]*x133;
                    const T x138 = tpi*x118 + PAz*x130 - PGz[g]*x131;
                    const T x139 = tpi*x120 + PAz*x131 - PGz[g]*x133;
                    const T x140 = PAz*x118 - PGz[g]*x120;
                    const T x141 = PBx*x92 - PGx[g]*x93;
                    const T x142 = PBx*x93 - PGx[g]*x95;
                    const T x143 = tpi_2*x85 + PAy*x141 - PGy[g]*x142;
                    const T x144 = PBx*x95 - PGx[g]*x99;
                    const T x145 = tpi_2*x89 + PAy*x142 - PGy[g]*x144;
                    const T x146 = x141 - x142;
                    const T x147 = tpi*x76 + PAy*x85 - PGy[g]*x89;
                    const T x148 = PAz*x141 - PGz[g]*x142;
                    const T x149 = PAz*x142 - PGz[g]*x144;
                    const T x150 = PAx*x141 + tpi*x118 - PGx[g]*x142;
                    const T x151 = PAx*x142 + tpi*x120 - PGx[g]*x144;
                    const T x152 = PAx*x118 - PGx[g]*x120;
                    const T x153 = PBy*x103 - PGy[g]*x104;
                    const T x154 = PBy*x104 - PGy[g]*x107;
                    const T x155 = PAx*x153 + tpi_2*x85 - PGx[g]*x154;
                    const T x156 = PBy*x107 - PGy[g]*x112;
                    const T x157 = PAx*x154 + tpi_2*x89 - PGx[g]*x156;
                    const T x158 = x153 - x154;
                    const T x159 = PAx*x85 + tpi*x83 - PGx[g]*x89;
                    const T x160 = PBz*x103 - PGz[g]*x104;
                    const T x161 = PBz*x104 - PGz[g]*x107;
                    const T x162 = PAx*x160 + tpi_2*x20 - PGx[g]*x161;
                    const T x163 = PBz*x107 - PGz[g]*x112;
                    const T x164 = PAx*x161 + tpi_2*x23 - PGx[g]*x163;
                    const T x165 = x160 - x161;
                    const T x166 = PAx*x20 + tpi*x26 - PGx[g]*x23;
                    const T x167 = PAy*x160 - PGy[g]*x161;
                    const T x168 = PAy*x161 - PGy[g]*x163;
                    const T x169 = tpi*x125 + PAz*x160 - PGz[g]*x161;
                    const T x170 = tpi*x127 + PAz*x161 - PGz[g]*x163;
                    const T x171 = PAz*x125 - PGz[g]*x127;
                    const T x172 = PAz*x153 - PGz[g]*x154;
                    const T x173 = PAz*x154 - PGz[g]*x156;
                    const T x174 = tpi*x125 + PAy*x153 - PGy[g]*x154;
                    const T x175 = tpi*x127 + PAy*x154 - PGy[g]*x156;
                    const T x176 = PAy*x125 - PGy[g]*x127;

                    integral_0[g] += PAx*x126 + tpi*x109 + tpi_3*x129 - PGx[g]*x128;
                    integral_1[g] += PAx*x155 + tpi*x158 + tpi_2*x159 - PGx[g]*x157;
                    integral_2[g] += PAx*x162 + tpi*x165 + tpi_2*x166 - PGx[g]*x164;
                    integral_3[g] += PAx*x150 + tpi*x146 + tpi*x152 - PGx[g]*x151;
                    integral_4[g] += PAx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                    integral_5[g] += PAx*x61 + tpi*x57 + tpi*x63 - PGx[g]*x62;
                    integral_6[g] += PAx*x97 + tpi*x102 - PGx[g]*x101;
                    integral_7[g] += PAx*x136 + tpi*x135 - PGx[g]*x137;
                    integral_8[g] += PAx*x51 + tpi*x54 - PGx[g]*x53;
                    integral_9[g] += PAx*x36 + tpi*x41 - PGx[g]*x40;
                    integral_10[g] += PAy*x126 - PGy[g]*x128;
                    integral_11[g] += PAx*x174 + tpi_2*x147 - PGx[g]*x175;
                    integral_12[g] += PAy*x162 - PGy[g]*x164;
                    integral_13[g] += PAx*x143 + tpi*x122 - PGx[g]*x145;
                    integral_14[g] += PAx*x21 + tpi*x27 - PGx[g]*x24;
                    integral_15[g] += PAy*x61 - PGy[g]*x62;
                    integral_16[g] += PAx*x119 - PGx[g]*x121;
                    integral_17[g] += PAx*x132 - PGx[g]*x134;
                    integral_18[g] += PAx*x66 - PGx[g]*x67;
                    integral_19[g] += PAx*x42 - PGx[g]*x43;
                    integral_20[g] += PAz*x126 - PGz[g]*x128;
                    integral_21[g] += PAz*x155 - PGz[g]*x157;
                    integral_22[g] += PAx*x169 + tpi_2*x77 - PGx[g]*x170;
                    integral_23[g] += PAz*x150 - PGz[g]*x151;
                    integral_24[g] += PAx*x86 + tpi*x84 - PGx[g]*x90;
                    integral_25[g] += PAx*x64 + tpi*x48 - PGx[g]*x65;
                    integral_26[g] += PAx*x116 - PGx[g]*x117;
                    integral_27[g] += PAx*x138 - PGx[g]*x139;
                    integral_28[g] += PAx*x69 - PGx[g]*x70;
                    integral_29[g] += PAx*x45 - PGx[g]*x47;
                    integral_30[g] += tpi*x109 + PAy*x110 - PGy[g]*x115;
                    integral_31[g] += tpi*x158 + PAy*x174 + tpi*x176 - PGy[g]*x175;
                    integral_32[g] += tpi*x165 + PAy*x167 - PGy[g]*x168;
                    integral_33[g] += tpi*x146 + tpi_2*x147 + PAy*x143 - PGy[g]*x145;
                    integral_34[g] += tpi*x16 + PAy*x21 + tpi*x25 - PGy[g]*x24;
                    integral_35[g] += tpi*x57 + PAy*x58 - PGy[g]*x60;
                    integral_36[g] += tpi*x102 + tpi_3*x122 + PAy*x119 - PGy[g]*x121;
                    integral_37[g] += tpi*x135 + tpi_2*x27 + PAy*x132 - PGy[g]*x134;
                    integral_38[g] += tpi*x54 + PAy*x66 + tpi*x68 - PGy[g]*x67;
                    integral_39[g] += tpi*x41 + PAy*x42 - PGy[g]*x43;
                    integral_40[g] += PAy*x123 - PGy[g]*x124;
                    integral_41[g] += PAz*x174 - PGz[g]*x175;
                    integral_42[g] += PAy*x169 - PGy[g]*x170;
                    integral_43[g] += PAz*x143 - PGz[g]*x145;
                    integral_44[g] += tpi*x77 + PAy*x86 - PGy[g]*x90;
                    integral_45[g] += PAy*x64 - PGy[g]*x65;
                    integral_46[g] += PAz*x119 - PGz[g]*x121;
                    integral_47[g] += tpi_2*x84 + PAy*x138 - PGy[g]*x139;
                    integral_48[g] += tpi*x48 + PAy*x69 - PGy[g]*x70;
                    integral_49[g] += PAy*x45 - PGy[g]*x47;
                    integral_50[g] += tpi*x109 + PAz*x123 - PGz[g]*x124;
                    integral_51[g] += tpi*x158 + PAz*x172 - PGz[g]*x173;
                    integral_52[g] += tpi*x165 + PAz*x169 + tpi*x171 - PGz[g]*x170;
                    integral_53[g] += tpi*x146 + PAz*x148 - PGz[g]*x149;
                    integral_54[g] += tpi*x16 + PAz*x86 + tpi*x91 - PGz[g]*x90;
                    integral_55[g] += tpi*x57 + tpi_2*x77 + PAz*x64 - PGz[g]*x65;
                    integral_56[g] += tpi*x102 + PAz*x116 - PGz[g]*x117;
                    integral_57[g] += tpi*x135 + PAz*x138 + tpi*x140 - PGz[g]*x139;
                    integral_58[g] += tpi*x54 + tpi_2*x84 + PAz*x69 - PGz[g]*x70;
                    integral_59[g] += tpi*x41 + tpi_3*x48 + PAz*x45 - PGz[g]*x47;

                } // for g
            } // if/else prim
        } // for prim

    }
} // VRR_D_F

template void VRR_D_F<float>(
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

template void VRR_D_F<double>(
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

