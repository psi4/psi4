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
        enum { g_block = (sizeof(T) <= 4) ? 32 : 16 };
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
                const T tpi_p2 = tpi*tpi;
                const T tpi_p2_2 = 2*tpi_p2;

                // overwrite integrals for first primitive, accumulate otherwise
                if (prim==0) {
                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x1 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x3 = PBx*x0 + tpi*x2 - PGx[g]*x1;
                        const T x4 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x6 = PBx*x1 + tpi*x5 - PGx[g]*x4;
                        const T x7 = x0 - x1;
                        const T x8 = PBx*x3 + tpi_2*x7 - PGx[g]*x6;
                        const T x9 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x10 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x11 = PBx*x4 + tpi*x10 - PGx[g]*x9;
                        const T x12 = x1 - x4;
                        const T x13 = PBx*x6 + tpi_2*x12 - PGx[g]*x11;
                        const T x14 = x7 - x12;
                        const T x15 = PAz*x8 - PGz[g]*x13;
                        const T x16 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                        const T x17 = scaled_boys_3[g] - scaled_boys_4[g];
                        const T x18 = PBx*x9 + tpi*x17 - PGx[g]*x16;
                        const T x19 = x4 - x9;
                        const T x20 = PBx*x11 + tpi_2*x19 - PGx[g]*x18;
                        const T x21 = PAz*x13 - PGz[g]*x20;
                        block_40[g] = PAy*x15 - PGy[g]*x21;
                        const T x22 = x8 - x13;
                        block_50[g] = tpi*x22 + PAz*x15 - PGz[g]*x21;
                        const T x23 = PAy*x8 - PGy[g]*x13;
                        const T x24 = PAy*x13 - PGy[g]*x20;
                        block_30[g] = tpi*x22 + PAy*x23 - PGy[g]*x24;
                        const T x25 = x3 - x6;
                        const T x26 = tpi*x25;
                        const T x27 = PAx*x8 + 3*x26 - PGx[g]*x13;
                        const T x28 = x6 - x11;
                        const T x29 = tpi*x28;
                        const T x30 = PAx*x13 + 3*x29 - PGx[g]*x20;
                        const T x31 = PAx*x26 + tpi_p2_2*x14 - PGx[g]*x29;
                        block_0[g] = PAx*x27 + tpi*x22 + 3*x31 - PGx[g]*x30;
                        block_10[g] = PAy*x27 - PGy[g]*x30;
                        block_20[g] = PAz*x27 - PGz[g]*x30;
                        const T x32 = PAy*x26 - PGy[g]*x29;
                        const T x33 = PBy*x3 - PGy[g]*x6;
                        const T x34 = PBy*x6 - PGy[g]*x11;
                        const T x35 = x26 + PAy*x33 - PGy[g]*x34;
                        const T x36 = PBy*x11 - PGy[g]*x18;
                        const T x37 = x29 + PAy*x34 - PGy[g]*x36;
                        const T x38 = x33 - x34;
                        block_31[g] = tpi*x38 + PAy*x35 + x32 - PGy[g]*x37;
                        const T x39 = PBz*x3 - PGz[g]*x6;
                        const T x40 = PBz*x6 - PGz[g]*x11;
                        const T x41 = PBz*x11 - PGz[g]*x18;
                        const T x42 = x39 - x40;
                        const T x43 = PAy*x39 - PGy[g]*x40;
                        const T x44 = PAy*x40 - PGy[g]*x41;
                        block_32[g] = tpi*x42 + PAy*x43 - PGy[g]*x44;
                        const T x45 = x26 + PAz*x39 - PGz[g]*x40;
                        const T x46 = x29 + PAz*x40 - PGz[g]*x41;
                        const T x47 = PAz*x26 - PGz[g]*x29;
                        block_52[g] = tpi*x42 + x47 + PAz*x45 - PGz[g]*x46;
                        block_42[g] = PAy*x45 - PGy[g]*x46;
                        const T x48 = PAz*x33 - PGz[g]*x34;
                        const T x49 = PAz*x34 - PGz[g]*x36;
                        block_41[g] = PAy*x48 + x47 - PGy[g]*x49;
                        block_51[g] = tpi*x38 + PAz*x48 - PGz[g]*x49;
                        const T x50 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x51 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x52 = x50 - x51;
                        const T x53 = tpi*x52;
                        const T x54 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x55 = x51 - x54;
                        const T x56 = tpi*x55;
                        const T x57 = x52 - x55;
                        const T x58 = PBx*x53 - PGx[g]*x56;
                        const T x59 = PAx*x33 + 2*x58 - PGx[g]*x34;
                        const T x60 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x61 = x54 - x60;
                        const T x62 = tpi*x61;
                        const T x63 = PBx*x56 - PGx[g]*x62;
                        const T x64 = PAx*x34 + 2*x63 - PGx[g]*x36;
                        const T x65 = PAx*x58 + tpi_p2*x57 - PGx[g]*x63;
                        block_1[g] = PAx*x59 + tpi*x38 + 2*x65 - PGx[g]*x64;
                        block_11[g] = x31 + PAy*x59 - PGy[g]*x64;
                        block_21[g] = PAz*x59 - PGz[g]*x64;
                        const T x66 = tpi*x2 + PBy*x50 - PGy[g]*x51;
                        const T x67 = tpi*x5 + PBy*x51 - PGy[g]*x54;
                        const T x68 = 2*x53 + PBy*x66 - PGy[g]*x67;
                        const T x69 = tpi*x10 + PBy*x54 - PGy[g]*x60;
                        const T x70 = 2*x56 + PBy*x67 - PGy[g]*x69;
                        const T x71 = PAx*x68 - PGx[g]*x70;
                        const T x72 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                        const T x73 = tpi*x17 + PBy*x60 - PGy[g]*x72;
                        const T x74 = 2*x62 + PBy*x69 - PGy[g]*x73;
                        const T x75 = PAx*x70 - PGx[g]*x74;
                        const T x76 = x68 - x70;
                        block_6[g] = PAx*x71 + tpi*x76 - PGx[g]*x75;
                        const T x77 = PAz*x68 - PGz[g]*x70;
                        const T x78 = PAz*x70 - PGz[g]*x74;
                        block_56[g] = tpi*x76 + PAz*x77 - PGz[g]*x78;
                        block_26[g] = PAx*x77 - PGx[g]*x78;
                        const T x79 = x66 - x67;
                        const T x80 = tpi*x79;
                        const T x81 = 3*x80 + PAy*x68 - PGy[g]*x70;
                        const T x82 = x67 - x69;
                        const T x83 = tpi*x82;
                        const T x84 = 3*x83 + PAy*x70 - PGy[g]*x74;
                        const T x85 = tpi_p2_2*x57 + PAy*x80 - PGy[g]*x83;
                        block_36[g] = tpi*x76 + 3*x85 + PAy*x81 - PGy[g]*x84;
                        block_16[g] = PAx*x81 - PGx[g]*x84;
                        block_46[g] = PAz*x81 - PGz[g]*x84;
                        const T x86 = PAx*x80 - PGx[g]*x83;
                        const T x87 = PBx*x66 - PGx[g]*x67;
                        const T x88 = PBx*x67 - PGx[g]*x69;
                        const T x89 = PAx*x87 + x80 - PGx[g]*x88;
                        const T x90 = PBx*x69 - PGx[g]*x73;
                        const T x91 = PAx*x88 + x83 - PGx[g]*x90;
                        const T x92 = x87 - x88;
                        block_3[g] = PAx*x89 + tpi*x92 + x86 - PGx[g]*x91;
                        const T x93 = PBz*x66 - PGz[g]*x67;
                        const T x94 = PBz*x67 - PGz[g]*x69;
                        const T x95 = PBz*x69 - PGz[g]*x73;
                        const T x96 = PAx*x93 - PGx[g]*x94;
                        const T x97 = PAx*x94 - PGx[g]*x95;
                        const T x98 = x93 - x94;
                        block_7[g] = PAx*x96 + tpi*x98 - PGx[g]*x97;
                        const T x99 = x80 + PAz*x93 - PGz[g]*x94;
                        const T x100 = x83 + PAz*x94 - PGz[g]*x95;
                        const T x101 = PAz*x80 - PGz[g]*x83;
                        block_57[g] = tpi*x98 + x101 + PAz*x99 - PGz[g]*x100;
                        block_27[g] = PAx*x99 - PGx[g]*x100;
                        const T x102 = PAz*x87 - PGz[g]*x88;
                        const T x103 = PAz*x88 - PGz[g]*x90;
                        block_23[g] = PAx*x102 + x101 - PGx[g]*x103;
                        block_53[g] = tpi*x92 + PAz*x102 - PGz[g]*x103;
                        const T x104 = 2*x58 + PAy*x87 - PGy[g]*x88;
                        const T x105 = 2*x63 + PAy*x88 - PGy[g]*x90;
                        const T x106 = tpi_p2*x14 + PAy*x58 - PGy[g]*x63;
                        block_33[g] = tpi*x92 + PAy*x104 + 2*x106 - PGy[g]*x105;
                        block_13[g] = PAx*x104 + x85 - PGx[g]*x105;
                        block_43[g] = PAz*x104 - PGz[g]*x105;
                        const T x107 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x108 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x109 = PBy*x107 - PGy[g]*x108;
                        const T x110 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x111 = PBy*x108 - PGy[g]*x110;
                        const T x112 = x109 - x111;
                        const T x113 = tpi*x112;
                        const T x114 = 2*x113 + PAy*x93 - PGy[g]*x94;
                        const T x115 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x116 = PBy*x110 - PGy[g]*x115;
                        const T x117 = x111 - x116;
                        const T x118 = tpi*x117;
                        const T x119 = 2*x118 + PAy*x94 - PGy[g]*x95;
                        block_47[g] = x85 + PAz*x114 - PGz[g]*x119;
                        block_17[g] = PAx*x114 - PGx[g]*x119;
                        const T x120 = x107 - x108;
                        const T x121 = tpi*x2 + PBz*x107 - PGz[g]*x108;
                        const T x122 = x108 - x110;
                        const T x123 = tpi*x5 + PBz*x108 - PGz[g]*x110;
                        const T x124 = x120 - x122;
                        const T x125 = tpi_p2*x124 + PAy*x113 - PGy[g]*x118;
                        block_37[g] = tpi*x98 + 2*x125 + PAy*x114 - PGy[g]*x119;
                        const T x126 = tpi*x120;
                        const T x127 = tpi*x122;
                        const T x128 = PBx*x126 - PGx[g]*x127;
                        const T x129 = PAx*x39 + 2*x128 - PGx[g]*x40;
                        const T x130 = x110 - x115;
                        const T x131 = tpi*x10 + PBz*x110 - PGz[g]*x115;
                        const T x132 = tpi*x130;
                        const T x133 = PBx*x127 - PGx[g]*x132;
                        const T x134 = PAx*x40 + 2*x133 - PGx[g]*x41;
                        block_22[g] = x31 + PAz*x129 - PGz[g]*x134;
                        const T x135 = PAx*x128 + tpi_p2*x124 - PGx[g]*x133;
                        block_2[g] = PAx*x129 + tpi*x42 + 2*x135 - PGx[g]*x134;
                        block_12[g] = PAy*x129 - PGy[g]*x134;
                        const T x136 = 2*x126 + PBz*x121 - PGz[g]*x123;
                        const T x137 = 2*x127 + PBz*x123 - PGz[g]*x131;
                        const T x138 = PAx*x136 - PGx[g]*x137;
                        const T x139 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                        const T x140 = tpi*x17 + PBz*x115 - PGz[g]*x139;
                        const T x141 = PBy*x115 - PGy[g]*x139;
                        const T x142 = 2*x132 + PBz*x131 - PGz[g]*x140;
                        const T x143 = PAx*x137 - PGx[g]*x142;
                        const T x144 = x136 - x137;
                        block_9[g] = PAx*x138 + tpi*x144 - PGx[g]*x143;
                        const T x145 = PBx*x109 - PGx[g]*x111;
                        const T x146 = PBx*x111 - PGx[g]*x116;
                        const T x147 = PBx*x116 - PGx[g]*x141;
                        const T x148 = x128 + PAy*x145 - PGy[g]*x146;
                        const T x149 = x133 + PAy*x146 - PGy[g]*x147;
                        block_14[g] = PAx*x148 + x125 - PGx[g]*x149;
                        const T x150 = PAy*x128 - PGy[g]*x133;
                        const T x151 = x145 - x146;
                        block_34[g] = tpi*x151 + PAy*x148 + x150 - PGy[g]*x149;
                        const T x152 = PAx*x145 + x113 - PGx[g]*x146;
                        const T x153 = PAx*x146 + x118 - PGx[g]*x147;
                        const T x154 = PAx*x113 - PGx[g]*x118;
                        block_4[g] = PAx*x152 + tpi*x151 + x154 - PGx[g]*x153;
                        const T x155 = x58 + PAz*x145 - PGz[g]*x146;
                        const T x156 = x63 + PAz*x146 - PGz[g]*x147;
                        const T x157 = PAz*x58 - PGz[g]*x63;
                        block_54[g] = tpi*x151 + PAz*x155 + x157 - PGz[g]*x156;
                        const T x158 = tpi_p2*x57 + PAz*x113 - PGz[g]*x118;
                        block_24[g] = PAx*x155 + x158 - PGx[g]*x156;
                        const T x159 = tpi_p2*x14 + PAz*x128 - PGz[g]*x133;
                        block_44[g] = PAy*x155 + x159 - PGy[g]*x156;
                        const T x160 = PAy*x136 - PGy[g]*x137;
                        const T x161 = PAy*x137 - PGy[g]*x142;
                        block_39[g] = tpi*x144 + PAy*x160 - PGy[g]*x161;
                        block_19[g] = PAx*x160 - PGx[g]*x161;
                        const T x162 = x121 - x123;
                        const T x163 = tpi*x162;
                        const T x164 = 3*x163 + PAz*x136 - PGz[g]*x137;
                        const T x165 = x123 - x131;
                        const T x166 = tpi*x165;
                        const T x167 = 3*x166 + PAz*x137 - PGz[g]*x142;
                        const T x168 = tpi_p2_2*x124 + PAz*x163 - PGz[g]*x166;
                        block_59[g] = tpi*x144 + 3*x168 + PAz*x164 - PGz[g]*x167;
                        block_29[g] = PAx*x164 - PGx[g]*x167;
                        block_49[g] = PAy*x164 - PGy[g]*x167;
                        const T x169 = PBx*x121 - PGx[g]*x123;
                        const T x170 = PBx*x123 - PGx[g]*x131;
                        const T x171 = 2*x128 + PAz*x169 - PGz[g]*x170;
                        const T x172 = PBx*x131 - PGx[g]*x140;
                        const T x173 = 2*x133 + PAz*x170 - PGz[g]*x172;
                        block_25[g] = PAx*x171 + x168 - PGx[g]*x173;
                        const T x174 = PBy*x121 - PGy[g]*x123;
                        const T x175 = PBy*x123 - PGy[g]*x131;
                        const T x176 = PBy*x131 - PGy[g]*x140;
                        const T x177 = 2*x113 + PAz*x174 - PGz[g]*x175;
                        const T x178 = 2*x118 + PAz*x175 - PGz[g]*x176;
                        block_48[g] = x168 + PAy*x177 - PGy[g]*x178;
                        const T x179 = x169 - x170;
                        block_55[g] = tpi*x179 + PAz*x171 + 2*x159 - PGz[g]*x173;
                        block_45[g] = PAy*x171 - PGy[g]*x173;
                        const T x180 = PAx*x169 + x163 - PGx[g]*x170;
                        const T x181 = PAx*x170 + x166 - PGx[g]*x172;
                        const T x182 = PAx*x163 - PGx[g]*x166;
                        block_5[g] = PAx*x180 + tpi*x179 + x182 - PGx[g]*x181;
                        const T x183 = PAy*x169 - PGy[g]*x170;
                        const T x184 = PAy*x170 - PGy[g]*x172;
                        block_35[g] = tpi*x179 + PAy*x183 - PGy[g]*x184;
                        const T x185 = PAy*x163 - PGy[g]*x166;
                        block_15[g] = PAx*x183 + x185 - PGx[g]*x184;
                        const T x186 = x166 + PAy*x175 - PGy[g]*x176;
                        const T x187 = x163 + PAy*x174 - PGy[g]*x175;
                        const T x188 = x174 - x175;
                        block_38[g] = tpi*x188 + x185 + PAy*x187 - PGy[g]*x186;
                        block_18[g] = PAx*x187 - PGx[g]*x186;
                        block_58[g] = tpi*x188 + 2*x158 + PAz*x177 - PGz[g]*x178;
                        block_28[g] = PAx*x177 - PGx[g]*x178;
                        const T x189 = PAx*x174 - PGx[g]*x175;
                        const T x190 = PAx*x175 - PGx[g]*x176;
                        block_8[g] = PAx*x189 + tpi*x188 - PGx[g]*x190;

                    } // for g
                } else {

                    #pragma omp simd
                    for (int g = 0; g < g_block; ++g) {

                        const T x0 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x1 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x3 = PBx*x0 + tpi*x2 - PGx[g]*x1;
                        const T x4 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x6 = PBx*x1 + tpi*x5 - PGx[g]*x4;
                        const T x7 = x0 - x1;
                        const T x8 = PBx*x3 + tpi_2*x7 - PGx[g]*x6;
                        const T x9 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x10 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x11 = PBx*x4 + tpi*x10 - PGx[g]*x9;
                        const T x12 = x1 - x4;
                        const T x13 = PBx*x6 + tpi_2*x12 - PGx[g]*x11;
                        const T x14 = x7 - x12;
                        const T x15 = PAz*x8 - PGz[g]*x13;
                        const T x16 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                        const T x17 = scaled_boys_3[g] - scaled_boys_4[g];
                        const T x18 = PBx*x9 + tpi*x17 - PGx[g]*x16;
                        const T x19 = x4 - x9;
                        const T x20 = PBx*x11 + tpi_2*x19 - PGx[g]*x18;
                        const T x21 = PAz*x13 - PGz[g]*x20;
                        block_40[g] += PAy*x15 - PGy[g]*x21;
                        const T x22 = x8 - x13;
                        block_50[g] += tpi*x22 + PAz*x15 - PGz[g]*x21;
                        const T x23 = PAy*x8 - PGy[g]*x13;
                        const T x24 = PAy*x13 - PGy[g]*x20;
                        block_30[g] += tpi*x22 + PAy*x23 - PGy[g]*x24;
                        const T x25 = x3 - x6;
                        const T x26 = tpi*x25;
                        const T x27 = PAx*x8 + 3*x26 - PGx[g]*x13;
                        const T x28 = x6 - x11;
                        const T x29 = tpi*x28;
                        const T x30 = PAx*x13 + 3*x29 - PGx[g]*x20;
                        const T x31 = PAx*x26 + tpi_p2_2*x14 - PGx[g]*x29;
                        block_0[g] += PAx*x27 + tpi*x22 + 3*x31 - PGx[g]*x30;
                        block_10[g] += PAy*x27 - PGy[g]*x30;
                        block_20[g] += PAz*x27 - PGz[g]*x30;
                        const T x32 = PAy*x26 - PGy[g]*x29;
                        const T x33 = PBy*x3 - PGy[g]*x6;
                        const T x34 = PBy*x6 - PGy[g]*x11;
                        const T x35 = x26 + PAy*x33 - PGy[g]*x34;
                        const T x36 = PBy*x11 - PGy[g]*x18;
                        const T x37 = x29 + PAy*x34 - PGy[g]*x36;
                        const T x38 = x33 - x34;
                        block_31[g] += tpi*x38 + PAy*x35 + x32 - PGy[g]*x37;
                        const T x39 = PBz*x3 - PGz[g]*x6;
                        const T x40 = PBz*x6 - PGz[g]*x11;
                        const T x41 = PBz*x11 - PGz[g]*x18;
                        const T x42 = x39 - x40;
                        const T x43 = PAy*x39 - PGy[g]*x40;
                        const T x44 = PAy*x40 - PGy[g]*x41;
                        block_32[g] += tpi*x42 + PAy*x43 - PGy[g]*x44;
                        const T x45 = x26 + PAz*x39 - PGz[g]*x40;
                        const T x46 = x29 + PAz*x40 - PGz[g]*x41;
                        const T x47 = PAz*x26 - PGz[g]*x29;
                        block_52[g] += tpi*x42 + x47 + PAz*x45 - PGz[g]*x46;
                        block_42[g] += PAy*x45 - PGy[g]*x46;
                        const T x48 = PAz*x33 - PGz[g]*x34;
                        const T x49 = PAz*x34 - PGz[g]*x36;
                        block_41[g] += PAy*x48 + x47 - PGy[g]*x49;
                        block_51[g] += tpi*x38 + PAz*x48 - PGz[g]*x49;
                        const T x50 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x51 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x52 = x50 - x51;
                        const T x53 = tpi*x52;
                        const T x54 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x55 = x51 - x54;
                        const T x56 = tpi*x55;
                        const T x57 = x52 - x55;
                        const T x58 = PBx*x53 - PGx[g]*x56;
                        const T x59 = PAx*x33 + 2*x58 - PGx[g]*x34;
                        const T x60 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x61 = x54 - x60;
                        const T x62 = tpi*x61;
                        const T x63 = PBx*x56 - PGx[g]*x62;
                        const T x64 = PAx*x34 + 2*x63 - PGx[g]*x36;
                        const T x65 = PAx*x58 + tpi_p2*x57 - PGx[g]*x63;
                        block_1[g] += PAx*x59 + tpi*x38 + 2*x65 - PGx[g]*x64;
                        block_11[g] += x31 + PAy*x59 - PGy[g]*x64;
                        block_21[g] += PAz*x59 - PGz[g]*x64;
                        const T x66 = tpi*x2 + PBy*x50 - PGy[g]*x51;
                        const T x67 = tpi*x5 + PBy*x51 - PGy[g]*x54;
                        const T x68 = 2*x53 + PBy*x66 - PGy[g]*x67;
                        const T x69 = tpi*x10 + PBy*x54 - PGy[g]*x60;
                        const T x70 = 2*x56 + PBy*x67 - PGy[g]*x69;
                        const T x71 = PAx*x68 - PGx[g]*x70;
                        const T x72 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                        const T x73 = tpi*x17 + PBy*x60 - PGy[g]*x72;
                        const T x74 = 2*x62 + PBy*x69 - PGy[g]*x73;
                        const T x75 = PAx*x70 - PGx[g]*x74;
                        const T x76 = x68 - x70;
                        block_6[g] += PAx*x71 + tpi*x76 - PGx[g]*x75;
                        const T x77 = PAz*x68 - PGz[g]*x70;
                        const T x78 = PAz*x70 - PGz[g]*x74;
                        block_56[g] += tpi*x76 + PAz*x77 - PGz[g]*x78;
                        block_26[g] += PAx*x77 - PGx[g]*x78;
                        const T x79 = x66 - x67;
                        const T x80 = tpi*x79;
                        const T x81 = 3*x80 + PAy*x68 - PGy[g]*x70;
                        const T x82 = x67 - x69;
                        const T x83 = tpi*x82;
                        const T x84 = 3*x83 + PAy*x70 - PGy[g]*x74;
                        const T x85 = tpi_p2_2*x57 + PAy*x80 - PGy[g]*x83;
                        block_36[g] += tpi*x76 + 3*x85 + PAy*x81 - PGy[g]*x84;
                        block_16[g] += PAx*x81 - PGx[g]*x84;
                        block_46[g] += PAz*x81 - PGz[g]*x84;
                        const T x86 = PAx*x80 - PGx[g]*x83;
                        const T x87 = PBx*x66 - PGx[g]*x67;
                        const T x88 = PBx*x67 - PGx[g]*x69;
                        const T x89 = PAx*x87 + x80 - PGx[g]*x88;
                        const T x90 = PBx*x69 - PGx[g]*x73;
                        const T x91 = PAx*x88 + x83 - PGx[g]*x90;
                        const T x92 = x87 - x88;
                        block_3[g] += PAx*x89 + tpi*x92 + x86 - PGx[g]*x91;
                        const T x93 = PBz*x66 - PGz[g]*x67;
                        const T x94 = PBz*x67 - PGz[g]*x69;
                        const T x95 = PBz*x69 - PGz[g]*x73;
                        const T x96 = PAx*x93 - PGx[g]*x94;
                        const T x97 = PAx*x94 - PGx[g]*x95;
                        const T x98 = x93 - x94;
                        block_7[g] += PAx*x96 + tpi*x98 - PGx[g]*x97;
                        const T x99 = x80 + PAz*x93 - PGz[g]*x94;
                        const T x100 = x83 + PAz*x94 - PGz[g]*x95;
                        const T x101 = PAz*x80 - PGz[g]*x83;
                        block_57[g] += tpi*x98 + x101 + PAz*x99 - PGz[g]*x100;
                        block_27[g] += PAx*x99 - PGx[g]*x100;
                        const T x102 = PAz*x87 - PGz[g]*x88;
                        const T x103 = PAz*x88 - PGz[g]*x90;
                        block_23[g] += PAx*x102 + x101 - PGx[g]*x103;
                        block_53[g] += tpi*x92 + PAz*x102 - PGz[g]*x103;
                        const T x104 = 2*x58 + PAy*x87 - PGy[g]*x88;
                        const T x105 = 2*x63 + PAy*x88 - PGy[g]*x90;
                        const T x106 = tpi_p2*x14 + PAy*x58 - PGy[g]*x63;
                        block_33[g] += tpi*x92 + PAy*x104 + 2*x106 - PGy[g]*x105;
                        block_13[g] += PAx*x104 + x85 - PGx[g]*x105;
                        block_43[g] += PAz*x104 - PGz[g]*x105;
                        const T x107 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                        const T x108 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                        const T x109 = PBy*x107 - PGy[g]*x108;
                        const T x110 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                        const T x111 = PBy*x108 - PGy[g]*x110;
                        const T x112 = x109 - x111;
                        const T x113 = tpi*x112;
                        const T x114 = 2*x113 + PAy*x93 - PGy[g]*x94;
                        const T x115 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                        const T x116 = PBy*x110 - PGy[g]*x115;
                        const T x117 = x111 - x116;
                        const T x118 = tpi*x117;
                        const T x119 = 2*x118 + PAy*x94 - PGy[g]*x95;
                        block_47[g] += x85 + PAz*x114 - PGz[g]*x119;
                        block_17[g] += PAx*x114 - PGx[g]*x119;
                        const T x120 = x107 - x108;
                        const T x121 = tpi*x2 + PBz*x107 - PGz[g]*x108;
                        const T x122 = x108 - x110;
                        const T x123 = tpi*x5 + PBz*x108 - PGz[g]*x110;
                        const T x124 = x120 - x122;
                        const T x125 = tpi_p2*x124 + PAy*x113 - PGy[g]*x118;
                        block_37[g] += tpi*x98 + 2*x125 + PAy*x114 - PGy[g]*x119;
                        const T x126 = tpi*x120;
                        const T x127 = tpi*x122;
                        const T x128 = PBx*x126 - PGx[g]*x127;
                        const T x129 = PAx*x39 + 2*x128 - PGx[g]*x40;
                        const T x130 = x110 - x115;
                        const T x131 = tpi*x10 + PBz*x110 - PGz[g]*x115;
                        const T x132 = tpi*x130;
                        const T x133 = PBx*x127 - PGx[g]*x132;
                        const T x134 = PAx*x40 + 2*x133 - PGx[g]*x41;
                        block_22[g] += x31 + PAz*x129 - PGz[g]*x134;
                        const T x135 = PAx*x128 + tpi_p2*x124 - PGx[g]*x133;
                        block_2[g] += PAx*x129 + tpi*x42 + 2*x135 - PGx[g]*x134;
                        block_12[g] += PAy*x129 - PGy[g]*x134;
                        const T x136 = 2*x126 + PBz*x121 - PGz[g]*x123;
                        const T x137 = 2*x127 + PBz*x123 - PGz[g]*x131;
                        const T x138 = PAx*x136 - PGx[g]*x137;
                        const T x139 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                        const T x140 = tpi*x17 + PBz*x115 - PGz[g]*x139;
                        const T x141 = PBy*x115 - PGy[g]*x139;
                        const T x142 = 2*x132 + PBz*x131 - PGz[g]*x140;
                        const T x143 = PAx*x137 - PGx[g]*x142;
                        const T x144 = x136 - x137;
                        block_9[g] += PAx*x138 + tpi*x144 - PGx[g]*x143;
                        const T x145 = PBx*x109 - PGx[g]*x111;
                        const T x146 = PBx*x111 - PGx[g]*x116;
                        const T x147 = PBx*x116 - PGx[g]*x141;
                        const T x148 = x128 + PAy*x145 - PGy[g]*x146;
                        const T x149 = x133 + PAy*x146 - PGy[g]*x147;
                        block_14[g] += PAx*x148 + x125 - PGx[g]*x149;
                        const T x150 = PAy*x128 - PGy[g]*x133;
                        const T x151 = x145 - x146;
                        block_34[g] += tpi*x151 + PAy*x148 + x150 - PGy[g]*x149;
                        const T x152 = PAx*x145 + x113 - PGx[g]*x146;
                        const T x153 = PAx*x146 + x118 - PGx[g]*x147;
                        const T x154 = PAx*x113 - PGx[g]*x118;
                        block_4[g] += PAx*x152 + tpi*x151 + x154 - PGx[g]*x153;
                        const T x155 = x58 + PAz*x145 - PGz[g]*x146;
                        const T x156 = x63 + PAz*x146 - PGz[g]*x147;
                        const T x157 = PAz*x58 - PGz[g]*x63;
                        block_54[g] += tpi*x151 + PAz*x155 + x157 - PGz[g]*x156;
                        const T x158 = tpi_p2*x57 + PAz*x113 - PGz[g]*x118;
                        block_24[g] += PAx*x155 + x158 - PGx[g]*x156;
                        const T x159 = tpi_p2*x14 + PAz*x128 - PGz[g]*x133;
                        block_44[g] += PAy*x155 + x159 - PGy[g]*x156;
                        const T x160 = PAy*x136 - PGy[g]*x137;
                        const T x161 = PAy*x137 - PGy[g]*x142;
                        block_39[g] += tpi*x144 + PAy*x160 - PGy[g]*x161;
                        block_19[g] += PAx*x160 - PGx[g]*x161;
                        const T x162 = x121 - x123;
                        const T x163 = tpi*x162;
                        const T x164 = 3*x163 + PAz*x136 - PGz[g]*x137;
                        const T x165 = x123 - x131;
                        const T x166 = tpi*x165;
                        const T x167 = 3*x166 + PAz*x137 - PGz[g]*x142;
                        const T x168 = tpi_p2_2*x124 + PAz*x163 - PGz[g]*x166;
                        block_59[g] += tpi*x144 + 3*x168 + PAz*x164 - PGz[g]*x167;
                        block_29[g] += PAx*x164 - PGx[g]*x167;
                        block_49[g] += PAy*x164 - PGy[g]*x167;
                        const T x169 = PBx*x121 - PGx[g]*x123;
                        const T x170 = PBx*x123 - PGx[g]*x131;
                        const T x171 = 2*x128 + PAz*x169 - PGz[g]*x170;
                        const T x172 = PBx*x131 - PGx[g]*x140;
                        const T x173 = 2*x133 + PAz*x170 - PGz[g]*x172;
                        block_25[g] += PAx*x171 + x168 - PGx[g]*x173;
                        const T x174 = PBy*x121 - PGy[g]*x123;
                        const T x175 = PBy*x123 - PGy[g]*x131;
                        const T x176 = PBy*x131 - PGy[g]*x140;
                        const T x177 = 2*x113 + PAz*x174 - PGz[g]*x175;
                        const T x178 = 2*x118 + PAz*x175 - PGz[g]*x176;
                        block_48[g] += x168 + PAy*x177 - PGy[g]*x178;
                        const T x179 = x169 - x170;
                        block_55[g] += tpi*x179 + PAz*x171 + 2*x159 - PGz[g]*x173;
                        block_45[g] += PAy*x171 - PGy[g]*x173;
                        const T x180 = PAx*x169 + x163 - PGx[g]*x170;
                        const T x181 = PAx*x170 + x166 - PGx[g]*x172;
                        const T x182 = PAx*x163 - PGx[g]*x166;
                        block_5[g] += PAx*x180 + tpi*x179 + x182 - PGx[g]*x181;
                        const T x183 = PAy*x169 - PGy[g]*x170;
                        const T x184 = PAy*x170 - PGy[g]*x172;
                        block_35[g] += tpi*x179 + PAy*x183 - PGy[g]*x184;
                        const T x185 = PAy*x163 - PGy[g]*x166;
                        block_15[g] += PAx*x183 + x185 - PGx[g]*x184;
                        const T x186 = x166 + PAy*x175 - PGy[g]*x176;
                        const T x187 = x163 + PAy*x174 - PGy[g]*x175;
                        const T x188 = x174 - x175;
                        block_38[g] += tpi*x188 + x185 + PAy*x187 - PGy[g]*x186;
                        block_18[g] += PAx*x187 - PGx[g]*x186;
                        block_58[g] += tpi*x188 + 2*x158 + PAz*x177 - PGz[g]*x178;
                        block_28[g] += PAx*x177 - PGx[g]*x178;
                        const T x189 = PAx*x174 - PGx[g]*x175;
                        const T x190 = PAx*x175 - PGx[g]*x176;
                        block_8[g] += PAx*x189 + tpi*x188 - PGx[g]*x190;

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
            const T tpi_p2 = tpi*tpi;
            const T tpi_p2_2 = 2*tpi_p2;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x1 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x3 = PBx*x0 + tpi*x2 - PGx[g]*x1;
                    const T x4 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x6 = PBx*x1 + tpi*x5 - PGx[g]*x4;
                    const T x7 = x0 - x1;
                    const T x8 = PBx*x3 + tpi_2*x7 - PGx[g]*x6;
                    const T x9 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x10 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x11 = PBx*x4 + tpi*x10 - PGx[g]*x9;
                    const T x12 = x1 - x4;
                    const T x13 = PBx*x6 + tpi_2*x12 - PGx[g]*x11;
                    const T x14 = x7 - x12;
                    const T x15 = PAz*x8 - PGz[g]*x13;
                    const T x16 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x17 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x18 = PBx*x9 + tpi*x17 - PGx[g]*x16;
                    const T x19 = x4 - x9;
                    const T x20 = PBx*x11 + tpi_2*x19 - PGx[g]*x18;
                    const T x21 = PAz*x13 - PGz[g]*x20;
                    const T x22 = x8 - x13;
                    const T x23 = PAy*x8 - PGy[g]*x13;
                    const T x24 = PAy*x13 - PGy[g]*x20;
                    const T x25 = x3 - x6;
                    const T x26 = tpi*x25;
                    const T x27 = PAx*x8 + 3*x26 - PGx[g]*x13;
                    const T x28 = x6 - x11;
                    const T x29 = tpi*x28;
                    const T x30 = PAx*x13 + 3*x29 - PGx[g]*x20;
                    const T x31 = PAx*x26 + tpi_p2_2*x14 - PGx[g]*x29;
                    const T x32 = PAy*x26 - PGy[g]*x29;
                    const T x33 = PBy*x3 - PGy[g]*x6;
                    const T x34 = PBy*x6 - PGy[g]*x11;
                    const T x35 = x26 + PAy*x33 - PGy[g]*x34;
                    const T x36 = PBy*x11 - PGy[g]*x18;
                    const T x37 = x29 + PAy*x34 - PGy[g]*x36;
                    const T x38 = x33 - x34;
                    const T x39 = PBz*x3 - PGz[g]*x6;
                    const T x40 = PBz*x6 - PGz[g]*x11;
                    const T x41 = PBz*x11 - PGz[g]*x18;
                    const T x42 = x39 - x40;
                    const T x43 = PAy*x39 - PGy[g]*x40;
                    const T x44 = PAy*x40 - PGy[g]*x41;
                    const T x45 = x26 + PAz*x39 - PGz[g]*x40;
                    const T x46 = x29 + PAz*x40 - PGz[g]*x41;
                    const T x47 = PAz*x26 - PGz[g]*x29;
                    const T x48 = PAz*x33 - PGz[g]*x34;
                    const T x49 = PAz*x34 - PGz[g]*x36;
                    const T x50 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x51 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x52 = x50 - x51;
                    const T x53 = tpi*x52;
                    const T x54 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x55 = x51 - x54;
                    const T x56 = tpi*x55;
                    const T x57 = x52 - x55;
                    const T x58 = PBx*x53 - PGx[g]*x56;
                    const T x59 = PAx*x33 + 2*x58 - PGx[g]*x34;
                    const T x60 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x61 = x54 - x60;
                    const T x62 = tpi*x61;
                    const T x63 = PBx*x56 - PGx[g]*x62;
                    const T x64 = PAx*x34 + 2*x63 - PGx[g]*x36;
                    const T x65 = PAx*x58 + tpi_p2*x57 - PGx[g]*x63;
                    const T x66 = tpi*x2 + PBy*x50 - PGy[g]*x51;
                    const T x67 = tpi*x5 + PBy*x51 - PGy[g]*x54;
                    const T x68 = 2*x53 + PBy*x66 - PGy[g]*x67;
                    const T x69 = tpi*x10 + PBy*x54 - PGy[g]*x60;
                    const T x70 = 2*x56 + PBy*x67 - PGy[g]*x69;
                    const T x71 = PAx*x68 - PGx[g]*x70;
                    const T x72 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x73 = tpi*x17 + PBy*x60 - PGy[g]*x72;
                    const T x74 = 2*x62 + PBy*x69 - PGy[g]*x73;
                    const T x75 = PAx*x70 - PGx[g]*x74;
                    const T x76 = x68 - x70;
                    const T x77 = PAz*x68 - PGz[g]*x70;
                    const T x78 = PAz*x70 - PGz[g]*x74;
                    const T x79 = x66 - x67;
                    const T x80 = tpi*x79;
                    const T x81 = 3*x80 + PAy*x68 - PGy[g]*x70;
                    const T x82 = x67 - x69;
                    const T x83 = tpi*x82;
                    const T x84 = 3*x83 + PAy*x70 - PGy[g]*x74;
                    const T x85 = tpi_p2_2*x57 + PAy*x80 - PGy[g]*x83;
                    const T x86 = PAx*x80 - PGx[g]*x83;
                    const T x87 = PBx*x66 - PGx[g]*x67;
                    const T x88 = PBx*x67 - PGx[g]*x69;
                    const T x89 = PAx*x87 + x80 - PGx[g]*x88;
                    const T x90 = PBx*x69 - PGx[g]*x73;
                    const T x91 = PAx*x88 + x83 - PGx[g]*x90;
                    const T x92 = x87 - x88;
                    const T x93 = PBz*x66 - PGz[g]*x67;
                    const T x94 = PBz*x67 - PGz[g]*x69;
                    const T x95 = PBz*x69 - PGz[g]*x73;
                    const T x96 = PAx*x93 - PGx[g]*x94;
                    const T x97 = PAx*x94 - PGx[g]*x95;
                    const T x98 = x93 - x94;
                    const T x99 = x80 + PAz*x93 - PGz[g]*x94;
                    const T x100 = x83 + PAz*x94 - PGz[g]*x95;
                    const T x101 = PAz*x80 - PGz[g]*x83;
                    const T x102 = PAz*x87 - PGz[g]*x88;
                    const T x103 = PAz*x88 - PGz[g]*x90;
                    const T x104 = 2*x58 + PAy*x87 - PGy[g]*x88;
                    const T x105 = 2*x63 + PAy*x88 - PGy[g]*x90;
                    const T x106 = tpi_p2*x14 + PAy*x58 - PGy[g]*x63;
                    const T x107 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x108 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x109 = PBy*x107 - PGy[g]*x108;
                    const T x110 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x111 = PBy*x108 - PGy[g]*x110;
                    const T x112 = x109 - x111;
                    const T x113 = tpi*x112;
                    const T x114 = 2*x113 + PAy*x93 - PGy[g]*x94;
                    const T x115 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x116 = PBy*x110 - PGy[g]*x115;
                    const T x117 = x111 - x116;
                    const T x118 = tpi*x117;
                    const T x119 = 2*x118 + PAy*x94 - PGy[g]*x95;
                    const T x120 = x107 - x108;
                    const T x121 = tpi*x2 + PBz*x107 - PGz[g]*x108;
                    const T x122 = x108 - x110;
                    const T x123 = tpi*x5 + PBz*x108 - PGz[g]*x110;
                    const T x124 = x120 - x122;
                    const T x125 = tpi_p2*x124 + PAy*x113 - PGy[g]*x118;
                    const T x126 = tpi*x120;
                    const T x127 = tpi*x122;
                    const T x128 = PBx*x126 - PGx[g]*x127;
                    const T x129 = PAx*x39 + 2*x128 - PGx[g]*x40;
                    const T x130 = x110 - x115;
                    const T x131 = tpi*x10 + PBz*x110 - PGz[g]*x115;
                    const T x132 = tpi*x130;
                    const T x133 = PBx*x127 - PGx[g]*x132;
                    const T x134 = PAx*x40 + 2*x133 - PGx[g]*x41;
                    const T x135 = PAx*x128 + tpi_p2*x124 - PGx[g]*x133;
                    const T x136 = 2*x126 + PBz*x121 - PGz[g]*x123;
                    const T x137 = 2*x127 + PBz*x123 - PGz[g]*x131;
                    const T x138 = PAx*x136 - PGx[g]*x137;
                    const T x139 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x140 = tpi*x17 + PBz*x115 - PGz[g]*x139;
                    const T x141 = PBy*x115 - PGy[g]*x139;
                    const T x142 = 2*x132 + PBz*x131 - PGz[g]*x140;
                    const T x143 = PAx*x137 - PGx[g]*x142;
                    const T x144 = x136 - x137;
                    const T x145 = PBx*x109 - PGx[g]*x111;
                    const T x146 = PBx*x111 - PGx[g]*x116;
                    const T x147 = PBx*x116 - PGx[g]*x141;
                    const T x148 = x128 + PAy*x145 - PGy[g]*x146;
                    const T x149 = x133 + PAy*x146 - PGy[g]*x147;
                    const T x150 = PAy*x128 - PGy[g]*x133;
                    const T x151 = x145 - x146;
                    const T x152 = PAx*x145 + x113 - PGx[g]*x146;
                    const T x153 = PAx*x146 + x118 - PGx[g]*x147;
                    const T x154 = PAx*x113 - PGx[g]*x118;
                    const T x155 = x58 + PAz*x145 - PGz[g]*x146;
                    const T x156 = x63 + PAz*x146 - PGz[g]*x147;
                    const T x157 = PAz*x58 - PGz[g]*x63;
                    const T x158 = tpi_p2*x57 + PAz*x113 - PGz[g]*x118;
                    const T x159 = tpi_p2*x14 + PAz*x128 - PGz[g]*x133;
                    const T x160 = PAy*x136 - PGy[g]*x137;
                    const T x161 = PAy*x137 - PGy[g]*x142;
                    const T x162 = x121 - x123;
                    const T x163 = tpi*x162;
                    const T x164 = 3*x163 + PAz*x136 - PGz[g]*x137;
                    const T x165 = x123 - x131;
                    const T x166 = tpi*x165;
                    const T x167 = 3*x166 + PAz*x137 - PGz[g]*x142;
                    const T x168 = tpi_p2_2*x124 + PAz*x163 - PGz[g]*x166;
                    const T x169 = PBx*x121 - PGx[g]*x123;
                    const T x170 = PBx*x123 - PGx[g]*x131;
                    const T x171 = 2*x128 + PAz*x169 - PGz[g]*x170;
                    const T x172 = PBx*x131 - PGx[g]*x140;
                    const T x173 = 2*x133 + PAz*x170 - PGz[g]*x172;
                    const T x174 = PBy*x121 - PGy[g]*x123;
                    const T x175 = PBy*x123 - PGy[g]*x131;
                    const T x176 = PBy*x131 - PGy[g]*x140;
                    const T x177 = 2*x113 + PAz*x174 - PGz[g]*x175;
                    const T x178 = 2*x118 + PAz*x175 - PGz[g]*x176;
                    const T x179 = x169 - x170;
                    const T x180 = PAx*x169 + x163 - PGx[g]*x170;
                    const T x181 = PAx*x170 + x166 - PGx[g]*x172;
                    const T x182 = PAx*x163 - PGx[g]*x166;
                    const T x183 = PAy*x169 - PGy[g]*x170;
                    const T x184 = PAy*x170 - PGy[g]*x172;
                    const T x185 = PAy*x163 - PGy[g]*x166;
                    const T x186 = x166 + PAy*x175 - PGy[g]*x176;
                    const T x187 = x163 + PAy*x174 - PGy[g]*x175;
                    const T x188 = x174 - x175;
                    const T x189 = PAx*x174 - PGx[g]*x175;
                    const T x190 = PAx*x175 - PGx[g]*x176;

                    integral_0[g] = PAx*x27 + tpi*x22 + 3*x31 - PGx[g]*x30;
                    integral_1[g] = PAx*x59 + tpi*x38 + 2*x65 - PGx[g]*x64;
                    integral_2[g] = PAx*x129 + tpi*x42 + 2*x135 - PGx[g]*x134;
                    integral_3[g] = PAx*x89 + tpi*x92 + x86 - PGx[g]*x91;
                    integral_4[g] = PAx*x152 + tpi*x151 + x154 - PGx[g]*x153;
                    integral_5[g] = PAx*x180 + tpi*x179 + x182 - PGx[g]*x181;
                    integral_6[g] = PAx*x71 + tpi*x76 - PGx[g]*x75;
                    integral_7[g] = PAx*x96 + tpi*x98 - PGx[g]*x97;
                    integral_8[g] = PAx*x189 + tpi*x188 - PGx[g]*x190;
                    integral_9[g] = PAx*x138 + tpi*x144 - PGx[g]*x143;
                    integral_10[g] = PAy*x27 - PGy[g]*x30;
                    integral_11[g] = x31 + PAy*x59 - PGy[g]*x64;
                    integral_12[g] = PAy*x129 - PGy[g]*x134;
                    integral_13[g] = PAx*x104 + x85 - PGx[g]*x105;
                    integral_14[g] = PAx*x148 + x125 - PGx[g]*x149;
                    integral_15[g] = PAx*x183 + x185 - PGx[g]*x184;
                    integral_16[g] = PAx*x81 - PGx[g]*x84;
                    integral_17[g] = PAx*x114 - PGx[g]*x119;
                    integral_18[g] = PAx*x187 - PGx[g]*x186;
                    integral_19[g] = PAx*x160 - PGx[g]*x161;
                    integral_20[g] = PAz*x27 - PGz[g]*x30;
                    integral_21[g] = PAz*x59 - PGz[g]*x64;
                    integral_22[g] = x31 + PAz*x129 - PGz[g]*x134;
                    integral_23[g] = PAx*x102 + x101 - PGx[g]*x103;
                    integral_24[g] = PAx*x155 + x158 - PGx[g]*x156;
                    integral_25[g] = PAx*x171 + x168 - PGx[g]*x173;
                    integral_26[g] = PAx*x77 - PGx[g]*x78;
                    integral_27[g] = PAx*x99 - PGx[g]*x100;
                    integral_28[g] = PAx*x177 - PGx[g]*x178;
                    integral_29[g] = PAx*x164 - PGx[g]*x167;
                    integral_30[g] = tpi*x22 + PAy*x23 - PGy[g]*x24;
                    integral_31[g] = tpi*x38 + PAy*x35 + x32 - PGy[g]*x37;
                    integral_32[g] = tpi*x42 + PAy*x43 - PGy[g]*x44;
                    integral_33[g] = tpi*x92 + PAy*x104 + 2*x106 - PGy[g]*x105;
                    integral_34[g] = tpi*x151 + PAy*x148 + x150 - PGy[g]*x149;
                    integral_35[g] = tpi*x179 + PAy*x183 - PGy[g]*x184;
                    integral_36[g] = tpi*x76 + 3*x85 + PAy*x81 - PGy[g]*x84;
                    integral_37[g] = tpi*x98 + 2*x125 + PAy*x114 - PGy[g]*x119;
                    integral_38[g] = tpi*x188 + x185 + PAy*x187 - PGy[g]*x186;
                    integral_39[g] = tpi*x144 + PAy*x160 - PGy[g]*x161;
                    integral_40[g] = PAy*x15 - PGy[g]*x21;
                    integral_41[g] = PAy*x48 + x47 - PGy[g]*x49;
                    integral_42[g] = PAy*x45 - PGy[g]*x46;
                    integral_43[g] = PAz*x104 - PGz[g]*x105;
                    integral_44[g] = PAy*x155 + x159 - PGy[g]*x156;
                    integral_45[g] = PAy*x171 - PGy[g]*x173;
                    integral_46[g] = PAz*x81 - PGz[g]*x84;
                    integral_47[g] = x85 + PAz*x114 - PGz[g]*x119;
                    integral_48[g] = x168 + PAy*x177 - PGy[g]*x178;
                    integral_49[g] = PAy*x164 - PGy[g]*x167;
                    integral_50[g] = tpi*x22 + PAz*x15 - PGz[g]*x21;
                    integral_51[g] = tpi*x38 + PAz*x48 - PGz[g]*x49;
                    integral_52[g] = tpi*x42 + x47 + PAz*x45 - PGz[g]*x46;
                    integral_53[g] = tpi*x92 + PAz*x102 - PGz[g]*x103;
                    integral_54[g] = tpi*x151 + PAz*x155 + x157 - PGz[g]*x156;
                    integral_55[g] = tpi*x179 + PAz*x171 + 2*x159 - PGz[g]*x173;
                    integral_56[g] = tpi*x76 + PAz*x77 - PGz[g]*x78;
                    integral_57[g] = tpi*x98 + x101 + PAz*x99 - PGz[g]*x100;
                    integral_58[g] = tpi*x188 + 2*x158 + PAz*x177 - PGz[g]*x178;
                    integral_59[g] = tpi*x144 + 3*x168 + PAz*x164 - PGz[g]*x167;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < gb_size; ++g) {

                    const T x0 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x1 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x3 = PBx*x0 + tpi*x2 - PGx[g]*x1;
                    const T x4 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x5 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x6 = PBx*x1 + tpi*x5 - PGx[g]*x4;
                    const T x7 = x0 - x1;
                    const T x8 = PBx*x3 + tpi_2*x7 - PGx[g]*x6;
                    const T x9 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x10 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x11 = PBx*x4 + tpi*x10 - PGx[g]*x9;
                    const T x12 = x1 - x4;
                    const T x13 = PBx*x6 + tpi_2*x12 - PGx[g]*x11;
                    const T x14 = x7 - x12;
                    const T x15 = PAz*x8 - PGz[g]*x13;
                    const T x16 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x17 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x18 = PBx*x9 + tpi*x17 - PGx[g]*x16;
                    const T x19 = x4 - x9;
                    const T x20 = PBx*x11 + tpi_2*x19 - PGx[g]*x18;
                    const T x21 = PAz*x13 - PGz[g]*x20;
                    const T x22 = x8 - x13;
                    const T x23 = PAy*x8 - PGy[g]*x13;
                    const T x24 = PAy*x13 - PGy[g]*x20;
                    const T x25 = x3 - x6;
                    const T x26 = tpi*x25;
                    const T x27 = PAx*x8 + 3*x26 - PGx[g]*x13;
                    const T x28 = x6 - x11;
                    const T x29 = tpi*x28;
                    const T x30 = PAx*x13 + 3*x29 - PGx[g]*x20;
                    const T x31 = PAx*x26 + tpi_p2_2*x14 - PGx[g]*x29;
                    const T x32 = PAy*x26 - PGy[g]*x29;
                    const T x33 = PBy*x3 - PGy[g]*x6;
                    const T x34 = PBy*x6 - PGy[g]*x11;
                    const T x35 = x26 + PAy*x33 - PGy[g]*x34;
                    const T x36 = PBy*x11 - PGy[g]*x18;
                    const T x37 = x29 + PAy*x34 - PGy[g]*x36;
                    const T x38 = x33 - x34;
                    const T x39 = PBz*x3 - PGz[g]*x6;
                    const T x40 = PBz*x6 - PGz[g]*x11;
                    const T x41 = PBz*x11 - PGz[g]*x18;
                    const T x42 = x39 - x40;
                    const T x43 = PAy*x39 - PGy[g]*x40;
                    const T x44 = PAy*x40 - PGy[g]*x41;
                    const T x45 = x26 + PAz*x39 - PGz[g]*x40;
                    const T x46 = x29 + PAz*x40 - PGz[g]*x41;
                    const T x47 = PAz*x26 - PGz[g]*x29;
                    const T x48 = PAz*x33 - PGz[g]*x34;
                    const T x49 = PAz*x34 - PGz[g]*x36;
                    const T x50 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x51 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x52 = x50 - x51;
                    const T x53 = tpi*x52;
                    const T x54 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x55 = x51 - x54;
                    const T x56 = tpi*x55;
                    const T x57 = x52 - x55;
                    const T x58 = PBx*x53 - PGx[g]*x56;
                    const T x59 = PAx*x33 + 2*x58 - PGx[g]*x34;
                    const T x60 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x61 = x54 - x60;
                    const T x62 = tpi*x61;
                    const T x63 = PBx*x56 - PGx[g]*x62;
                    const T x64 = PAx*x34 + 2*x63 - PGx[g]*x36;
                    const T x65 = PAx*x58 + tpi_p2*x57 - PGx[g]*x63;
                    const T x66 = tpi*x2 + PBy*x50 - PGy[g]*x51;
                    const T x67 = tpi*x5 + PBy*x51 - PGy[g]*x54;
                    const T x68 = 2*x53 + PBy*x66 - PGy[g]*x67;
                    const T x69 = tpi*x10 + PBy*x54 - PGy[g]*x60;
                    const T x70 = 2*x56 + PBy*x67 - PGy[g]*x69;
                    const T x71 = PAx*x68 - PGx[g]*x70;
                    const T x72 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x73 = tpi*x17 + PBy*x60 - PGy[g]*x72;
                    const T x74 = 2*x62 + PBy*x69 - PGy[g]*x73;
                    const T x75 = PAx*x70 - PGx[g]*x74;
                    const T x76 = x68 - x70;
                    const T x77 = PAz*x68 - PGz[g]*x70;
                    const T x78 = PAz*x70 - PGz[g]*x74;
                    const T x79 = x66 - x67;
                    const T x80 = tpi*x79;
                    const T x81 = 3*x80 + PAy*x68 - PGy[g]*x70;
                    const T x82 = x67 - x69;
                    const T x83 = tpi*x82;
                    const T x84 = 3*x83 + PAy*x70 - PGy[g]*x74;
                    const T x85 = tpi_p2_2*x57 + PAy*x80 - PGy[g]*x83;
                    const T x86 = PAx*x80 - PGx[g]*x83;
                    const T x87 = PBx*x66 - PGx[g]*x67;
                    const T x88 = PBx*x67 - PGx[g]*x69;
                    const T x89 = PAx*x87 + x80 - PGx[g]*x88;
                    const T x90 = PBx*x69 - PGx[g]*x73;
                    const T x91 = PAx*x88 + x83 - PGx[g]*x90;
                    const T x92 = x87 - x88;
                    const T x93 = PBz*x66 - PGz[g]*x67;
                    const T x94 = PBz*x67 - PGz[g]*x69;
                    const T x95 = PBz*x69 - PGz[g]*x73;
                    const T x96 = PAx*x93 - PGx[g]*x94;
                    const T x97 = PAx*x94 - PGx[g]*x95;
                    const T x98 = x93 - x94;
                    const T x99 = x80 + PAz*x93 - PGz[g]*x94;
                    const T x100 = x83 + PAz*x94 - PGz[g]*x95;
                    const T x101 = PAz*x80 - PGz[g]*x83;
                    const T x102 = PAz*x87 - PGz[g]*x88;
                    const T x103 = PAz*x88 - PGz[g]*x90;
                    const T x104 = 2*x58 + PAy*x87 - PGy[g]*x88;
                    const T x105 = 2*x63 + PAy*x88 - PGy[g]*x90;
                    const T x106 = tpi_p2*x14 + PAy*x58 - PGy[g]*x63;
                    const T x107 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x108 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x109 = PBy*x107 - PGy[g]*x108;
                    const T x110 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x111 = PBy*x108 - PGy[g]*x110;
                    const T x112 = x109 - x111;
                    const T x113 = tpi*x112;
                    const T x114 = 2*x113 + PAy*x93 - PGy[g]*x94;
                    const T x115 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x116 = PBy*x110 - PGy[g]*x115;
                    const T x117 = x111 - x116;
                    const T x118 = tpi*x117;
                    const T x119 = 2*x118 + PAy*x94 - PGy[g]*x95;
                    const T x120 = x107 - x108;
                    const T x121 = tpi*x2 + PBz*x107 - PGz[g]*x108;
                    const T x122 = x108 - x110;
                    const T x123 = tpi*x5 + PBz*x108 - PGz[g]*x110;
                    const T x124 = x120 - x122;
                    const T x125 = tpi_p2*x124 + PAy*x113 - PGy[g]*x118;
                    const T x126 = tpi*x120;
                    const T x127 = tpi*x122;
                    const T x128 = PBx*x126 - PGx[g]*x127;
                    const T x129 = PAx*x39 + 2*x128 - PGx[g]*x40;
                    const T x130 = x110 - x115;
                    const T x131 = tpi*x10 + PBz*x110 - PGz[g]*x115;
                    const T x132 = tpi*x130;
                    const T x133 = PBx*x127 - PGx[g]*x132;
                    const T x134 = PAx*x40 + 2*x133 - PGx[g]*x41;
                    const T x135 = PAx*x128 + tpi_p2*x124 - PGx[g]*x133;
                    const T x136 = 2*x126 + PBz*x121 - PGz[g]*x123;
                    const T x137 = 2*x127 + PBz*x123 - PGz[g]*x131;
                    const T x138 = PAx*x136 - PGx[g]*x137;
                    const T x139 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x140 = tpi*x17 + PBz*x115 - PGz[g]*x139;
                    const T x141 = PBy*x115 - PGy[g]*x139;
                    const T x142 = 2*x132 + PBz*x131 - PGz[g]*x140;
                    const T x143 = PAx*x137 - PGx[g]*x142;
                    const T x144 = x136 - x137;
                    const T x145 = PBx*x109 - PGx[g]*x111;
                    const T x146 = PBx*x111 - PGx[g]*x116;
                    const T x147 = PBx*x116 - PGx[g]*x141;
                    const T x148 = x128 + PAy*x145 - PGy[g]*x146;
                    const T x149 = x133 + PAy*x146 - PGy[g]*x147;
                    const T x150 = PAy*x128 - PGy[g]*x133;
                    const T x151 = x145 - x146;
                    const T x152 = PAx*x145 + x113 - PGx[g]*x146;
                    const T x153 = PAx*x146 + x118 - PGx[g]*x147;
                    const T x154 = PAx*x113 - PGx[g]*x118;
                    const T x155 = x58 + PAz*x145 - PGz[g]*x146;
                    const T x156 = x63 + PAz*x146 - PGz[g]*x147;
                    const T x157 = PAz*x58 - PGz[g]*x63;
                    const T x158 = tpi_p2*x57 + PAz*x113 - PGz[g]*x118;
                    const T x159 = tpi_p2*x14 + PAz*x128 - PGz[g]*x133;
                    const T x160 = PAy*x136 - PGy[g]*x137;
                    const T x161 = PAy*x137 - PGy[g]*x142;
                    const T x162 = x121 - x123;
                    const T x163 = tpi*x162;
                    const T x164 = 3*x163 + PAz*x136 - PGz[g]*x137;
                    const T x165 = x123 - x131;
                    const T x166 = tpi*x165;
                    const T x167 = 3*x166 + PAz*x137 - PGz[g]*x142;
                    const T x168 = tpi_p2_2*x124 + PAz*x163 - PGz[g]*x166;
                    const T x169 = PBx*x121 - PGx[g]*x123;
                    const T x170 = PBx*x123 - PGx[g]*x131;
                    const T x171 = 2*x128 + PAz*x169 - PGz[g]*x170;
                    const T x172 = PBx*x131 - PGx[g]*x140;
                    const T x173 = 2*x133 + PAz*x170 - PGz[g]*x172;
                    const T x174 = PBy*x121 - PGy[g]*x123;
                    const T x175 = PBy*x123 - PGy[g]*x131;
                    const T x176 = PBy*x131 - PGy[g]*x140;
                    const T x177 = 2*x113 + PAz*x174 - PGz[g]*x175;
                    const T x178 = 2*x118 + PAz*x175 - PGz[g]*x176;
                    const T x179 = x169 - x170;
                    const T x180 = PAx*x169 + x163 - PGx[g]*x170;
                    const T x181 = PAx*x170 + x166 - PGx[g]*x172;
                    const T x182 = PAx*x163 - PGx[g]*x166;
                    const T x183 = PAy*x169 - PGy[g]*x170;
                    const T x184 = PAy*x170 - PGy[g]*x172;
                    const T x185 = PAy*x163 - PGy[g]*x166;
                    const T x186 = x166 + PAy*x175 - PGy[g]*x176;
                    const T x187 = x163 + PAy*x174 - PGy[g]*x175;
                    const T x188 = x174 - x175;
                    const T x189 = PAx*x174 - PGx[g]*x175;
                    const T x190 = PAx*x175 - PGx[g]*x176;

                    integral_0[g] += PAx*x27 + tpi*x22 + 3*x31 - PGx[g]*x30;
                    integral_1[g] += PAx*x59 + tpi*x38 + 2*x65 - PGx[g]*x64;
                    integral_2[g] += PAx*x129 + tpi*x42 + 2*x135 - PGx[g]*x134;
                    integral_3[g] += PAx*x89 + tpi*x92 + x86 - PGx[g]*x91;
                    integral_4[g] += PAx*x152 + tpi*x151 + x154 - PGx[g]*x153;
                    integral_5[g] += PAx*x180 + tpi*x179 + x182 - PGx[g]*x181;
                    integral_6[g] += PAx*x71 + tpi*x76 - PGx[g]*x75;
                    integral_7[g] += PAx*x96 + tpi*x98 - PGx[g]*x97;
                    integral_8[g] += PAx*x189 + tpi*x188 - PGx[g]*x190;
                    integral_9[g] += PAx*x138 + tpi*x144 - PGx[g]*x143;
                    integral_10[g] += PAy*x27 - PGy[g]*x30;
                    integral_11[g] += x31 + PAy*x59 - PGy[g]*x64;
                    integral_12[g] += PAy*x129 - PGy[g]*x134;
                    integral_13[g] += PAx*x104 + x85 - PGx[g]*x105;
                    integral_14[g] += PAx*x148 + x125 - PGx[g]*x149;
                    integral_15[g] += PAx*x183 + x185 - PGx[g]*x184;
                    integral_16[g] += PAx*x81 - PGx[g]*x84;
                    integral_17[g] += PAx*x114 - PGx[g]*x119;
                    integral_18[g] += PAx*x187 - PGx[g]*x186;
                    integral_19[g] += PAx*x160 - PGx[g]*x161;
                    integral_20[g] += PAz*x27 - PGz[g]*x30;
                    integral_21[g] += PAz*x59 - PGz[g]*x64;
                    integral_22[g] += x31 + PAz*x129 - PGz[g]*x134;
                    integral_23[g] += PAx*x102 + x101 - PGx[g]*x103;
                    integral_24[g] += PAx*x155 + x158 - PGx[g]*x156;
                    integral_25[g] += PAx*x171 + x168 - PGx[g]*x173;
                    integral_26[g] += PAx*x77 - PGx[g]*x78;
                    integral_27[g] += PAx*x99 - PGx[g]*x100;
                    integral_28[g] += PAx*x177 - PGx[g]*x178;
                    integral_29[g] += PAx*x164 - PGx[g]*x167;
                    integral_30[g] += tpi*x22 + PAy*x23 - PGy[g]*x24;
                    integral_31[g] += tpi*x38 + PAy*x35 + x32 - PGy[g]*x37;
                    integral_32[g] += tpi*x42 + PAy*x43 - PGy[g]*x44;
                    integral_33[g] += tpi*x92 + PAy*x104 + 2*x106 - PGy[g]*x105;
                    integral_34[g] += tpi*x151 + PAy*x148 + x150 - PGy[g]*x149;
                    integral_35[g] += tpi*x179 + PAy*x183 - PGy[g]*x184;
                    integral_36[g] += tpi*x76 + 3*x85 + PAy*x81 - PGy[g]*x84;
                    integral_37[g] += tpi*x98 + 2*x125 + PAy*x114 - PGy[g]*x119;
                    integral_38[g] += tpi*x188 + x185 + PAy*x187 - PGy[g]*x186;
                    integral_39[g] += tpi*x144 + PAy*x160 - PGy[g]*x161;
                    integral_40[g] += PAy*x15 - PGy[g]*x21;
                    integral_41[g] += PAy*x48 + x47 - PGy[g]*x49;
                    integral_42[g] += PAy*x45 - PGy[g]*x46;
                    integral_43[g] += PAz*x104 - PGz[g]*x105;
                    integral_44[g] += PAy*x155 + x159 - PGy[g]*x156;
                    integral_45[g] += PAy*x171 - PGy[g]*x173;
                    integral_46[g] += PAz*x81 - PGz[g]*x84;
                    integral_47[g] += x85 + PAz*x114 - PGz[g]*x119;
                    integral_48[g] += x168 + PAy*x177 - PGy[g]*x178;
                    integral_49[g] += PAy*x164 - PGy[g]*x167;
                    integral_50[g] += tpi*x22 + PAz*x15 - PGz[g]*x21;
                    integral_51[g] += tpi*x38 + PAz*x48 - PGz[g]*x49;
                    integral_52[g] += tpi*x42 + x47 + PAz*x45 - PGz[g]*x46;
                    integral_53[g] += tpi*x92 + PAz*x102 - PGz[g]*x103;
                    integral_54[g] += tpi*x151 + PAz*x155 + x157 - PGz[g]*x156;
                    integral_55[g] += tpi*x179 + PAz*x171 + 2*x159 - PGz[g]*x173;
                    integral_56[g] += tpi*x76 + PAz*x77 - PGz[g]*x78;
                    integral_57[g] += tpi*x98 + x101 + PAz*x99 - PGz[g]*x100;
                    integral_58[g] += tpi*x188 + 2*x158 + PAz*x177 - PGz[g]*x178;
                    integral_59[g] += tpi*x144 + 3*x168 + PAz*x164 - PGz[g]*x167;

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

