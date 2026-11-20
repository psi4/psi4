#include "VRR_F_D.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void VRR_F_D(
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
                        const T x10 = PBx*x5 + tpi*x9 - PGx[g]*x8;
                        const T x11 = PAz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                        const T x12 = PAy*x6 - PGy[g]*x11;
                        const T x13 = PAx*x7 - PGx[g]*x12;
                        const T x14 = x4 - x7;
                        const T x15 = PBx*x8 + tpi*x14 - PGx[g]*x13;
                        const T x16 = PBx*x9 - PGx[g]*x14;
                        const T x17 = x5 - x8;
                        block_24[g] = PBx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                        const T x18 = x0 - x1;
                        const T x19 = x1 - x3;
                        const T x20 = PAx*x18 - PGx[g]*x19;
                        const T x21 = tpi*x20 + PBy*x5 - PGy[g]*x8;
                        const T x22 = x3 - x6;
                        const T x23 = PAx*x19 - PGx[g]*x22;
                        const T x24 = tpi*x23 + PBy*x8 - PGy[g]*x13;
                        const T x25 = PBy*x20 - PGy[g]*x23;
                        block_27[g] = tpi*x17 + PBy*x21 + tpi*x25 - PGy[g]*x24;
                        const T x26 = x18 - x19;
                        const T x27 = tpi*x26 + PBy*x9 - PGy[g]*x14;
                        block_25[g] = PBx*x21 + tpi*x27 - PGx[g]*x24;
                        const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x29 = tpi*x28 + PAz*x0 - PGz[g]*x1;
                        const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x31 = tpi*x30 + PAz*x1 - PGz[g]*x3;
                        const T x32 = tpi_2*x18 + PAz*x29 - PGz[g]*x31;
                        const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x34 = tpi*x33 + PAz*x3 - PGz[g]*x6;
                        const T x35 = tpi_2*x19 + PAz*x31 - PGz[g]*x34;
                        const T x36 = PBx*x32 - PGx[g]*x35;
                        const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                        const T x38 = tpi*x37 + PAz*x6 - PGz[g]*x11;
                        const T x39 = tpi_2*x22 + PAz*x34 - PGz[g]*x38;
                        const T x40 = PBx*x35 - PGx[g]*x39;
                        const T x41 = x32 - x35;
                        block_54[g] = PBx*x36 + tpi*x41 - PGx[g]*x40;
                        const T x42 = PBy*x32 - PGy[g]*x35;
                        const T x43 = PBy*x35 - PGy[g]*x39;
                        block_55[g] = PBx*x42 - PGx[g]*x43;
                        block_57[g] = tpi*x41 + PBy*x42 - PGy[g]*x43;
                        const T x44 = x29 - x31;
                        const T x45 = tpi_3*x44 + PBz*x32 - PGz[g]*x35;
                        const T x46 = x31 - x34;
                        const T x47 = tpi_3*x46 + PBz*x35 - PGz[g]*x39;
                        block_56[g] = PBx*x45 - PGx[g]*x47;
                        block_58[g] = PBy*x45 - PGy[g]*x47;
                        const T x48 = tpi_2*x26 + PBz*x44 - PGz[g]*x46;
                        block_59[g] = tpi_3*x48 + tpi*x41 + PBz*x45 - PGz[g]*x47;
                        const T x49 = PBx*x29 - PGx[g]*x31;
                        const T x50 = PBx*x31 - PGx[g]*x34;
                        const T x51 = PBx*x49 + tpi*x44 - PGx[g]*x50;
                        const T x52 = PBx*x34 - PGx[g]*x38;
                        const T x53 = PBx*x50 + tpi*x46 - PGx[g]*x52;
                        const T x54 = x49 - x50;
                        block_30[g] = PAx*x51 + tpi_2*x54 - PGx[g]*x53;
                        block_48[g] = PAy*x51 - PGy[g]*x53;
                        const T x55 = PAx*x29 - PGx[g]*x31;
                        const T x56 = PAx*x31 - PGx[g]*x34;
                        const T x57 = tpi_2*x20 + PBz*x55 - PGz[g]*x56;
                        const T x58 = PAx*x34 - PGx[g]*x38;
                        const T x59 = tpi_2*x23 + PBz*x56 - PGz[g]*x58;
                        block_32[g] = PBx*x57 + tpi*x48 - PGx[g]*x59;
                        block_34[g] = PBy*x57 - PGy[g]*x59;
                        const T x60 = PBy*x55 - PGy[g]*x56;
                        const T x61 = PBy*x56 - PGy[g]*x58;
                        const T x62 = x55 - x56;
                        block_33[g] = PBy*x60 + tpi*x62 - PGy[g]*x61;
                        const T x63 = PBy*x44 - PGy[g]*x46;
                        block_31[g] = PBx*x60 + tpi*x63 - PGx[g]*x61;
                        const T x64 = PAy*x29 - PGy[g]*x31;
                        const T x65 = PAy*x31 - PGy[g]*x34;
                        const T x66 = tpi*x44 + PBy*x64 - PGy[g]*x65;
                        const T x67 = PAy*x34 - PGy[g]*x38;
                        const T x68 = tpi*x46 + PBy*x65 - PGy[g]*x67;
                        block_49[g] = PBx*x66 - PGx[g]*x68;
                        const T x69 = x64 - x65;
                        block_51[g] = tpi*x63 + PBy*x66 + tpi*x69 - PGy[g]*x68;
                        const T x70 = tpi_2*x9 + PBz*x64 - PGz[g]*x65;
                        const T x71 = tpi_2*x14 + PBz*x65 - PGz[g]*x67;
                        block_52[g] = tpi*x48 + PBy*x70 - PGy[g]*x71;
                        block_50[g] = PBx*x70 - PGx[g]*x71;
                        const T x72 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x73 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x74 = x72 - x73;
                        const T x75 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x76 = x73 - x75;
                        const T x77 = x74 - x76;
                        const T x78 = tpi*x77 + PBz*x20 - PGz[g]*x23;
                        block_35[g] = tpi_2*x78 + PBz*x57 + tpi*x62 - PGz[g]*x59;
                        const T x79 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x80 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x81 = x79 - x80;
                        const T x82 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x83 = x80 - x82;
                        const T x84 = x81 - x83;
                        const T x85 = tpi*x84 + PBz*x9 - PGz[g]*x14;
                        block_53[g] = tpi_2*x85 + PBz*x70 + tpi*x69 - PGz[g]*x71;
                        const T x86 = PAx*x81 - PGx[g]*x83;
                        const T x87 = tpi*x86 + PBz*x5 - PGz[g]*x8;
                        const T x88 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x89 = x82 - x88;
                        const T x90 = PAx*x83 - PGx[g]*x89;
                        const T x91 = tpi*x90 + PBz*x8 - PGz[g]*x13;
                        const T x92 = PBz*x86 - PGz[g]*x90;
                        block_29[g] = tpi*x17 + PBz*x87 + tpi*x92 - PGz[g]*x91;
                        block_26[g] = PBx*x87 + tpi*x85 - PGx[g]*x91;
                        block_28[g] = tpi*x78 + PBy*x87 - PGy[g]*x91;
                        const T x93 = tpi*x28 + PAy*x79 - PGy[g]*x80;
                        const T x94 = tpi*x30 + PAy*x80 - PGy[g]*x82;
                        const T x95 = tpi_2*x81 + PAy*x93 - PGy[g]*x94;
                        const T x96 = tpi*x33 + PAy*x82 - PGy[g]*x88;
                        const T x97 = tpi_2*x83 + PAy*x94 - PGy[g]*x96;
                        const T x98 = PBx*x95 - PGx[g]*x97;
                        const T x99 = PAy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                        const T x100 = tpi*x37 + PAy*x88 - PGy[g]*x99;
                        const T x101 = tpi_2*x89 + PAy*x96 - PGy[g]*x100;
                        const T x102 = PBx*x97 - PGx[g]*x101;
                        const T x103 = x95 - x97;
                        block_36[g] = PBx*x98 + tpi*x103 - PGx[g]*x102;
                        const T x104 = PAx*x72 + tpi*x28 - PGx[g]*x73;
                        const T x105 = PAx*x73 + tpi*x30 - PGx[g]*x75;
                        const T x106 = PAx*x104 + tpi_2*x74 - PGx[g]*x105;
                        const T x107 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x108 = PAx*x75 + tpi*x33 - PGx[g]*x107;
                        const T x109 = PAx*x105 + tpi_2*x76 - PGx[g]*x108;
                        const T x110 = x106 - x109;
                        const T x111 = PBy*x106 - PGy[g]*x109;
                        const T x112 = PAx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                        const T x113 = PAx*x107 + tpi*x37 - PGx[g]*x112;
                        const T x114 = x75 - x107;
                        const T x115 = PAx*x108 + tpi_2*x114 - PGx[g]*x113;
                        const T x116 = PBy*x109 - PGy[g]*x115;
                        block_3[g] = tpi*x110 + PBy*x111 - PGy[g]*x116;
                        const T x117 = PBz*x106 - PGz[g]*x109;
                        const T x118 = PBz*x109 - PGz[g]*x115;
                        block_4[g] = PBy*x117 - PGy[g]*x118;
                        block_5[g] = tpi*x110 + PBz*x117 - PGz[g]*x118;
                        const T x119 = x104 - x105;
                        const T x120 = PBx*x106 + tpi_3*x119 - PGx[g]*x109;
                        const T x121 = x105 - x108;
                        const T x122 = PBx*x109 + tpi_3*x121 - PGx[g]*x115;
                        const T x123 = PBx*x119 + tpi_2*x77 - PGx[g]*x121;
                        block_0[g] = PBx*x120 + tpi_3*x123 + tpi*x110 - PGx[g]*x122;
                        block_1[g] = PBy*x120 - PGy[g]*x122;
                        block_2[g] = PBz*x120 - PGz[g]*x122;
                        const T x124 = PBz*x95 - PGz[g]*x97;
                        const T x125 = PBz*x97 - PGz[g]*x101;
                        block_38[g] = PBx*x124 - PGx[g]*x125;
                        block_41[g] = tpi*x103 + PBz*x124 - PGz[g]*x125;
                        const T x126 = x93 - x94;
                        const T x127 = tpi_3*x126 + PBy*x95 - PGy[g]*x97;
                        const T x128 = x94 - x96;
                        const T x129 = tpi_3*x128 + PBy*x97 - PGy[g]*x101;
                        block_37[g] = PBx*x127 - PGx[g]*x129;
                        block_40[g] = PBz*x127 - PGz[g]*x129;
                        const T x130 = tpi_2*x84 + PBy*x126 - PGy[g]*x128;
                        block_39[g] = tpi_3*x130 + tpi*x103 + PBy*x127 - PGy[g]*x129;
                        const T x131 = PBx*x93 - PGx[g]*x94;
                        const T x132 = PBx*x94 - PGx[g]*x96;
                        const T x133 = PBx*x131 + tpi*x126 - PGx[g]*x132;
                        const T x134 = PBx*x96 - PGx[g]*x100;
                        const T x135 = PBx*x132 + tpi*x128 - PGx[g]*x134;
                        const T x136 = x131 - x132;
                        block_18[g] = PAx*x133 + tpi_2*x136 - PGx[g]*x135;
                        block_42[g] = PAz*x133 - PGz[g]*x135;
                        const T x137 = PAz*x93 - PGz[g]*x94;
                        const T x138 = PAz*x94 - PGz[g]*x96;
                        const T x139 = tpi_2*x9 + PBy*x137 - PGy[g]*x138;
                        const T x140 = PAz*x96 - PGz[g]*x100;
                        const T x141 = tpi_2*x14 + PBy*x138 - PGy[g]*x140;
                        block_43[g] = PBx*x139 - PGx[g]*x141;
                        const T x142 = x137 - x138;
                        block_45[g] = tpi_2*x27 + PBy*x139 + tpi*x142 - PGy[g]*x141;
                        const T x143 = PBz*x126 - PGz[g]*x128;
                        const T x144 = tpi*x126 + PBz*x137 - PGz[g]*x138;
                        const T x145 = tpi*x128 + PBz*x138 - PGz[g]*x140;
                        block_47[g] = tpi*x143 + PBz*x144 + tpi*x142 - PGz[g]*x145;
                        block_46[g] = tpi_2*x85 + PBy*x144 - PGy[g]*x145;
                        block_44[g] = PBx*x144 - PGx[g]*x145;
                        const T x146 = PAx*x93 - PGx[g]*x94;
                        const T x147 = PAx*x94 - PGx[g]*x96;
                        const T x148 = tpi_2*x86 + PBy*x146 - PGy[g]*x147;
                        const T x149 = PAx*x96 - PGx[g]*x100;
                        const T x150 = tpi_2*x90 + PBy*x147 - PGy[g]*x149;
                        block_19[g] = PBx*x148 + tpi*x130 - PGx[g]*x150;
                        block_22[g] = PBz*x148 - PGz[g]*x150;
                        const T x151 = tpi*x77 + PBy*x86 - PGy[g]*x90;
                        const T x152 = x146 - x147;
                        block_21[g] = tpi_2*x151 + PBy*x148 + tpi*x152 - PGy[g]*x150;
                        const T x153 = PBz*x146 - PGz[g]*x147;
                        const T x154 = PBz*x147 - PGz[g]*x149;
                        block_20[g] = PBx*x153 + tpi*x143 - PGx[g]*x154;
                        block_23[g] = PBz*x153 + tpi*x152 - PGz[g]*x154;
                        const T x155 = PAy*x104 - PGy[g]*x105;
                        const T x156 = PAy*x105 - PGy[g]*x108;
                        const T x157 = PBx*x155 + tpi_2*x86 - PGx[g]*x156;
                        const T x158 = PAy*x108 - PGy[g]*x113;
                        const T x159 = PBx*x156 + tpi_2*x90 - PGx[g]*x158;
                        const T x160 = PBx*x86 + tpi*x84 - PGx[g]*x90;
                        const T x161 = x155 - x156;
                        block_6[g] = PBx*x157 + tpi_2*x160 + tpi*x161 - PGx[g]*x159;
                        block_8[g] = PBz*x157 - PGz[g]*x159;
                        const T x162 = PAz*x104 - PGz[g]*x105;
                        const T x163 = PAz*x105 - PGz[g]*x108;
                        const T x164 = PBx*x162 + tpi_2*x20 - PGx[g]*x163;
                        const T x165 = PAz*x108 - PGz[g]*x113;
                        const T x166 = PBx*x163 + tpi_2*x23 - PGx[g]*x165;
                        const T x167 = PBx*x20 + tpi*x26 - PGx[g]*x23;
                        const T x168 = x162 - x163;
                        block_12[g] = PBx*x164 + tpi_2*x167 + tpi*x168 - PGx[g]*x166;
                        block_13[g] = PBy*x164 - PGy[g]*x166;
                        const T x169 = PBz*x155 - PGz[g]*x156;
                        const T x170 = PBz*x156 - PGz[g]*x158;
                        block_11[g] = tpi*x161 + PBz*x169 - PGz[g]*x170;
                        const T x171 = tpi*x119 + PBy*x155 - PGy[g]*x156;
                        const T x172 = tpi*x121 + PBy*x156 - PGy[g]*x158;
                        block_7[g] = PBx*x171 + tpi_2*x151 - PGx[g]*x172;
                        const T x173 = PBy*x119 - PGy[g]*x121;
                        block_9[g] = tpi*x161 + PBy*x171 + tpi*x173 - PGy[g]*x172;
                        block_10[g] = PBz*x171 - PGz[g]*x172;
                        const T x174 = PBy*x162 - PGy[g]*x163;
                        const T x175 = PBy*x163 - PGy[g]*x165;
                        block_15[g] = tpi*x168 + PBy*x174 - PGy[g]*x175;
                        const T x176 = tpi*x119 + PBz*x162 - PGz[g]*x163;
                        const T x177 = tpi*x121 + PBz*x163 - PGz[g]*x165;
                        const T x178 = PBz*x119 - PGz[g]*x121;
                        block_17[g] = tpi*x168 + PBz*x176 + tpi*x178 - PGz[g]*x177;
                        block_14[g] = PBx*x176 + tpi_2*x78 - PGx[g]*x177;
                        block_16[g] = PBy*x176 - PGy[g]*x177;

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
                        const T x10 = PBx*x5 + tpi*x9 - PGx[g]*x8;
                        const T x11 = PAz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                        const T x12 = PAy*x6 - PGy[g]*x11;
                        const T x13 = PAx*x7 - PGx[g]*x12;
                        const T x14 = x4 - x7;
                        const T x15 = PBx*x8 + tpi*x14 - PGx[g]*x13;
                        const T x16 = PBx*x9 - PGx[g]*x14;
                        const T x17 = x5 - x8;
                        block_24[g] += PBx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                        const T x18 = x0 - x1;
                        const T x19 = x1 - x3;
                        const T x20 = PAx*x18 - PGx[g]*x19;
                        const T x21 = tpi*x20 + PBy*x5 - PGy[g]*x8;
                        const T x22 = x3 - x6;
                        const T x23 = PAx*x19 - PGx[g]*x22;
                        const T x24 = tpi*x23 + PBy*x8 - PGy[g]*x13;
                        const T x25 = PBy*x20 - PGy[g]*x23;
                        block_27[g] += tpi*x17 + PBy*x21 + tpi*x25 - PGy[g]*x24;
                        const T x26 = x18 - x19;
                        const T x27 = tpi*x26 + PBy*x9 - PGy[g]*x14;
                        block_25[g] += PBx*x21 + tpi*x27 - PGx[g]*x24;
                        const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                        const T x29 = tpi*x28 + PAz*x0 - PGz[g]*x1;
                        const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                        const T x31 = tpi*x30 + PAz*x1 - PGz[g]*x3;
                        const T x32 = tpi_2*x18 + PAz*x29 - PGz[g]*x31;
                        const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                        const T x34 = tpi*x33 + PAz*x3 - PGz[g]*x6;
                        const T x35 = tpi_2*x19 + PAz*x31 - PGz[g]*x34;
                        const T x36 = PBx*x32 - PGx[g]*x35;
                        const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                        const T x38 = tpi*x37 + PAz*x6 - PGz[g]*x11;
                        const T x39 = tpi_2*x22 + PAz*x34 - PGz[g]*x38;
                        const T x40 = PBx*x35 - PGx[g]*x39;
                        const T x41 = x32 - x35;
                        block_54[g] += PBx*x36 + tpi*x41 - PGx[g]*x40;
                        const T x42 = PBy*x32 - PGy[g]*x35;
                        const T x43 = PBy*x35 - PGy[g]*x39;
                        block_55[g] += PBx*x42 - PGx[g]*x43;
                        block_57[g] += tpi*x41 + PBy*x42 - PGy[g]*x43;
                        const T x44 = x29 - x31;
                        const T x45 = tpi_3*x44 + PBz*x32 - PGz[g]*x35;
                        const T x46 = x31 - x34;
                        const T x47 = tpi_3*x46 + PBz*x35 - PGz[g]*x39;
                        block_56[g] += PBx*x45 - PGx[g]*x47;
                        block_58[g] += PBy*x45 - PGy[g]*x47;
                        const T x48 = tpi_2*x26 + PBz*x44 - PGz[g]*x46;
                        block_59[g] += tpi_3*x48 + tpi*x41 + PBz*x45 - PGz[g]*x47;
                        const T x49 = PBx*x29 - PGx[g]*x31;
                        const T x50 = PBx*x31 - PGx[g]*x34;
                        const T x51 = PBx*x49 + tpi*x44 - PGx[g]*x50;
                        const T x52 = PBx*x34 - PGx[g]*x38;
                        const T x53 = PBx*x50 + tpi*x46 - PGx[g]*x52;
                        const T x54 = x49 - x50;
                        block_30[g] += PAx*x51 + tpi_2*x54 - PGx[g]*x53;
                        block_48[g] += PAy*x51 - PGy[g]*x53;
                        const T x55 = PAx*x29 - PGx[g]*x31;
                        const T x56 = PAx*x31 - PGx[g]*x34;
                        const T x57 = tpi_2*x20 + PBz*x55 - PGz[g]*x56;
                        const T x58 = PAx*x34 - PGx[g]*x38;
                        const T x59 = tpi_2*x23 + PBz*x56 - PGz[g]*x58;
                        block_32[g] += PBx*x57 + tpi*x48 - PGx[g]*x59;
                        block_34[g] += PBy*x57 - PGy[g]*x59;
                        const T x60 = PBy*x55 - PGy[g]*x56;
                        const T x61 = PBy*x56 - PGy[g]*x58;
                        const T x62 = x55 - x56;
                        block_33[g] += PBy*x60 + tpi*x62 - PGy[g]*x61;
                        const T x63 = PBy*x44 - PGy[g]*x46;
                        block_31[g] += PBx*x60 + tpi*x63 - PGx[g]*x61;
                        const T x64 = PAy*x29 - PGy[g]*x31;
                        const T x65 = PAy*x31 - PGy[g]*x34;
                        const T x66 = tpi*x44 + PBy*x64 - PGy[g]*x65;
                        const T x67 = PAy*x34 - PGy[g]*x38;
                        const T x68 = tpi*x46 + PBy*x65 - PGy[g]*x67;
                        block_49[g] += PBx*x66 - PGx[g]*x68;
                        const T x69 = x64 - x65;
                        block_51[g] += tpi*x63 + PBy*x66 + tpi*x69 - PGy[g]*x68;
                        const T x70 = tpi_2*x9 + PBz*x64 - PGz[g]*x65;
                        const T x71 = tpi_2*x14 + PBz*x65 - PGz[g]*x67;
                        block_52[g] += tpi*x48 + PBy*x70 - PGy[g]*x71;
                        block_50[g] += PBx*x70 - PGx[g]*x71;
                        const T x72 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                        const T x73 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                        const T x74 = x72 - x73;
                        const T x75 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                        const T x76 = x73 - x75;
                        const T x77 = x74 - x76;
                        const T x78 = tpi*x77 + PBz*x20 - PGz[g]*x23;
                        block_35[g] += tpi_2*x78 + PBz*x57 + tpi*x62 - PGz[g]*x59;
                        const T x79 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                        const T x80 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                        const T x81 = x79 - x80;
                        const T x82 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                        const T x83 = x80 - x82;
                        const T x84 = x81 - x83;
                        const T x85 = tpi*x84 + PBz*x9 - PGz[g]*x14;
                        block_53[g] += tpi_2*x85 + PBz*x70 + tpi*x69 - PGz[g]*x71;
                        const T x86 = PAx*x81 - PGx[g]*x83;
                        const T x87 = tpi*x86 + PBz*x5 - PGz[g]*x8;
                        const T x88 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                        const T x89 = x82 - x88;
                        const T x90 = PAx*x83 - PGx[g]*x89;
                        const T x91 = tpi*x90 + PBz*x8 - PGz[g]*x13;
                        const T x92 = PBz*x86 - PGz[g]*x90;
                        block_29[g] += tpi*x17 + PBz*x87 + tpi*x92 - PGz[g]*x91;
                        block_26[g] += PBx*x87 + tpi*x85 - PGx[g]*x91;
                        block_28[g] += tpi*x78 + PBy*x87 - PGy[g]*x91;
                        const T x93 = tpi*x28 + PAy*x79 - PGy[g]*x80;
                        const T x94 = tpi*x30 + PAy*x80 - PGy[g]*x82;
                        const T x95 = tpi_2*x81 + PAy*x93 - PGy[g]*x94;
                        const T x96 = tpi*x33 + PAy*x82 - PGy[g]*x88;
                        const T x97 = tpi_2*x83 + PAy*x94 - PGy[g]*x96;
                        const T x98 = PBx*x95 - PGx[g]*x97;
                        const T x99 = PAy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                        const T x100 = tpi*x37 + PAy*x88 - PGy[g]*x99;
                        const T x101 = tpi_2*x89 + PAy*x96 - PGy[g]*x100;
                        const T x102 = PBx*x97 - PGx[g]*x101;
                        const T x103 = x95 - x97;
                        block_36[g] += PBx*x98 + tpi*x103 - PGx[g]*x102;
                        const T x104 = PAx*x72 + tpi*x28 - PGx[g]*x73;
                        const T x105 = PAx*x73 + tpi*x30 - PGx[g]*x75;
                        const T x106 = PAx*x104 + tpi_2*x74 - PGx[g]*x105;
                        const T x107 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                        const T x108 = PAx*x75 + tpi*x33 - PGx[g]*x107;
                        const T x109 = PAx*x105 + tpi_2*x76 - PGx[g]*x108;
                        const T x110 = x106 - x109;
                        const T x111 = PBy*x106 - PGy[g]*x109;
                        const T x112 = PAx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                        const T x113 = PAx*x107 + tpi*x37 - PGx[g]*x112;
                        const T x114 = x75 - x107;
                        const T x115 = PAx*x108 + tpi_2*x114 - PGx[g]*x113;
                        const T x116 = PBy*x109 - PGy[g]*x115;
                        block_3[g] += tpi*x110 + PBy*x111 - PGy[g]*x116;
                        const T x117 = PBz*x106 - PGz[g]*x109;
                        const T x118 = PBz*x109 - PGz[g]*x115;
                        block_4[g] += PBy*x117 - PGy[g]*x118;
                        block_5[g] += tpi*x110 + PBz*x117 - PGz[g]*x118;
                        const T x119 = x104 - x105;
                        const T x120 = PBx*x106 + tpi_3*x119 - PGx[g]*x109;
                        const T x121 = x105 - x108;
                        const T x122 = PBx*x109 + tpi_3*x121 - PGx[g]*x115;
                        const T x123 = PBx*x119 + tpi_2*x77 - PGx[g]*x121;
                        block_0[g] += PBx*x120 + tpi_3*x123 + tpi*x110 - PGx[g]*x122;
                        block_1[g] += PBy*x120 - PGy[g]*x122;
                        block_2[g] += PBz*x120 - PGz[g]*x122;
                        const T x124 = PBz*x95 - PGz[g]*x97;
                        const T x125 = PBz*x97 - PGz[g]*x101;
                        block_38[g] += PBx*x124 - PGx[g]*x125;
                        block_41[g] += tpi*x103 + PBz*x124 - PGz[g]*x125;
                        const T x126 = x93 - x94;
                        const T x127 = tpi_3*x126 + PBy*x95 - PGy[g]*x97;
                        const T x128 = x94 - x96;
                        const T x129 = tpi_3*x128 + PBy*x97 - PGy[g]*x101;
                        block_37[g] += PBx*x127 - PGx[g]*x129;
                        block_40[g] += PBz*x127 - PGz[g]*x129;
                        const T x130 = tpi_2*x84 + PBy*x126 - PGy[g]*x128;
                        block_39[g] += tpi_3*x130 + tpi*x103 + PBy*x127 - PGy[g]*x129;
                        const T x131 = PBx*x93 - PGx[g]*x94;
                        const T x132 = PBx*x94 - PGx[g]*x96;
                        const T x133 = PBx*x131 + tpi*x126 - PGx[g]*x132;
                        const T x134 = PBx*x96 - PGx[g]*x100;
                        const T x135 = PBx*x132 + tpi*x128 - PGx[g]*x134;
                        const T x136 = x131 - x132;
                        block_18[g] += PAx*x133 + tpi_2*x136 - PGx[g]*x135;
                        block_42[g] += PAz*x133 - PGz[g]*x135;
                        const T x137 = PAz*x93 - PGz[g]*x94;
                        const T x138 = PAz*x94 - PGz[g]*x96;
                        const T x139 = tpi_2*x9 + PBy*x137 - PGy[g]*x138;
                        const T x140 = PAz*x96 - PGz[g]*x100;
                        const T x141 = tpi_2*x14 + PBy*x138 - PGy[g]*x140;
                        block_43[g] += PBx*x139 - PGx[g]*x141;
                        const T x142 = x137 - x138;
                        block_45[g] += tpi_2*x27 + PBy*x139 + tpi*x142 - PGy[g]*x141;
                        const T x143 = PBz*x126 - PGz[g]*x128;
                        const T x144 = tpi*x126 + PBz*x137 - PGz[g]*x138;
                        const T x145 = tpi*x128 + PBz*x138 - PGz[g]*x140;
                        block_47[g] += tpi*x143 + PBz*x144 + tpi*x142 - PGz[g]*x145;
                        block_46[g] += tpi_2*x85 + PBy*x144 - PGy[g]*x145;
                        block_44[g] += PBx*x144 - PGx[g]*x145;
                        const T x146 = PAx*x93 - PGx[g]*x94;
                        const T x147 = PAx*x94 - PGx[g]*x96;
                        const T x148 = tpi_2*x86 + PBy*x146 - PGy[g]*x147;
                        const T x149 = PAx*x96 - PGx[g]*x100;
                        const T x150 = tpi_2*x90 + PBy*x147 - PGy[g]*x149;
                        block_19[g] += PBx*x148 + tpi*x130 - PGx[g]*x150;
                        block_22[g] += PBz*x148 - PGz[g]*x150;
                        const T x151 = tpi*x77 + PBy*x86 - PGy[g]*x90;
                        const T x152 = x146 - x147;
                        block_21[g] += tpi_2*x151 + PBy*x148 + tpi*x152 - PGy[g]*x150;
                        const T x153 = PBz*x146 - PGz[g]*x147;
                        const T x154 = PBz*x147 - PGz[g]*x149;
                        block_20[g] += PBx*x153 + tpi*x143 - PGx[g]*x154;
                        block_23[g] += PBz*x153 + tpi*x152 - PGz[g]*x154;
                        const T x155 = PAy*x104 - PGy[g]*x105;
                        const T x156 = PAy*x105 - PGy[g]*x108;
                        const T x157 = PBx*x155 + tpi_2*x86 - PGx[g]*x156;
                        const T x158 = PAy*x108 - PGy[g]*x113;
                        const T x159 = PBx*x156 + tpi_2*x90 - PGx[g]*x158;
                        const T x160 = PBx*x86 + tpi*x84 - PGx[g]*x90;
                        const T x161 = x155 - x156;
                        block_6[g] += PBx*x157 + tpi_2*x160 + tpi*x161 - PGx[g]*x159;
                        block_8[g] += PBz*x157 - PGz[g]*x159;
                        const T x162 = PAz*x104 - PGz[g]*x105;
                        const T x163 = PAz*x105 - PGz[g]*x108;
                        const T x164 = PBx*x162 + tpi_2*x20 - PGx[g]*x163;
                        const T x165 = PAz*x108 - PGz[g]*x113;
                        const T x166 = PBx*x163 + tpi_2*x23 - PGx[g]*x165;
                        const T x167 = PBx*x20 + tpi*x26 - PGx[g]*x23;
                        const T x168 = x162 - x163;
                        block_12[g] += PBx*x164 + tpi_2*x167 + tpi*x168 - PGx[g]*x166;
                        block_13[g] += PBy*x164 - PGy[g]*x166;
                        const T x169 = PBz*x155 - PGz[g]*x156;
                        const T x170 = PBz*x156 - PGz[g]*x158;
                        block_11[g] += tpi*x161 + PBz*x169 - PGz[g]*x170;
                        const T x171 = tpi*x119 + PBy*x155 - PGy[g]*x156;
                        const T x172 = tpi*x121 + PBy*x156 - PGy[g]*x158;
                        block_7[g] += PBx*x171 + tpi_2*x151 - PGx[g]*x172;
                        const T x173 = PBy*x119 - PGy[g]*x121;
                        block_9[g] += tpi*x161 + PBy*x171 + tpi*x173 - PGy[g]*x172;
                        block_10[g] += PBz*x171 - PGz[g]*x172;
                        const T x174 = PBy*x162 - PGy[g]*x163;
                        const T x175 = PBy*x163 - PGy[g]*x165;
                        block_15[g] += tpi*x168 + PBy*x174 - PGy[g]*x175;
                        const T x176 = tpi*x119 + PBz*x162 - PGz[g]*x163;
                        const T x177 = tpi*x121 + PBz*x163 - PGz[g]*x165;
                        const T x178 = PBz*x119 - PGz[g]*x121;
                        block_17[g] += tpi*x168 + PBz*x176 + tpi*x178 - PGz[g]*x177;
                        block_14[g] += PBx*x176 + tpi_2*x78 - PGx[g]*x177;
                        block_16[g] += PBy*x176 - PGy[g]*x177;

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
                    const T x10 = PBx*x5 + tpi*x9 - PGx[g]*x8;
                    const T x11 = PAz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x12 = PAy*x6 - PGy[g]*x11;
                    const T x13 = PAx*x7 - PGx[g]*x12;
                    const T x14 = x4 - x7;
                    const T x15 = PBx*x8 + tpi*x14 - PGx[g]*x13;
                    const T x16 = PBx*x9 - PGx[g]*x14;
                    const T x17 = x5 - x8;
                    const T x18 = x0 - x1;
                    const T x19 = x1 - x3;
                    const T x20 = PAx*x18 - PGx[g]*x19;
                    const T x21 = tpi*x20 + PBy*x5 - PGy[g]*x8;
                    const T x22 = x3 - x6;
                    const T x23 = PAx*x19 - PGx[g]*x22;
                    const T x24 = tpi*x23 + PBy*x8 - PGy[g]*x13;
                    const T x25 = PBy*x20 - PGy[g]*x23;
                    const T x26 = x18 - x19;
                    const T x27 = tpi*x26 + PBy*x9 - PGy[g]*x14;
                    const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x29 = tpi*x28 + PAz*x0 - PGz[g]*x1;
                    const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x31 = tpi*x30 + PAz*x1 - PGz[g]*x3;
                    const T x32 = tpi_2*x18 + PAz*x29 - PGz[g]*x31;
                    const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x34 = tpi*x33 + PAz*x3 - PGz[g]*x6;
                    const T x35 = tpi_2*x19 + PAz*x31 - PGz[g]*x34;
                    const T x36 = PBx*x32 - PGx[g]*x35;
                    const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x38 = tpi*x37 + PAz*x6 - PGz[g]*x11;
                    const T x39 = tpi_2*x22 + PAz*x34 - PGz[g]*x38;
                    const T x40 = PBx*x35 - PGx[g]*x39;
                    const T x41 = x32 - x35;
                    const T x42 = PBy*x32 - PGy[g]*x35;
                    const T x43 = PBy*x35 - PGy[g]*x39;
                    const T x44 = x29 - x31;
                    const T x45 = tpi_3*x44 + PBz*x32 - PGz[g]*x35;
                    const T x46 = x31 - x34;
                    const T x47 = tpi_3*x46 + PBz*x35 - PGz[g]*x39;
                    const T x48 = tpi_2*x26 + PBz*x44 - PGz[g]*x46;
                    const T x49 = PBx*x29 - PGx[g]*x31;
                    const T x50 = PBx*x31 - PGx[g]*x34;
                    const T x51 = PBx*x49 + tpi*x44 - PGx[g]*x50;
                    const T x52 = PBx*x34 - PGx[g]*x38;
                    const T x53 = PBx*x50 + tpi*x46 - PGx[g]*x52;
                    const T x54 = x49 - x50;
                    const T x55 = PAx*x29 - PGx[g]*x31;
                    const T x56 = PAx*x31 - PGx[g]*x34;
                    const T x57 = tpi_2*x20 + PBz*x55 - PGz[g]*x56;
                    const T x58 = PAx*x34 - PGx[g]*x38;
                    const T x59 = tpi_2*x23 + PBz*x56 - PGz[g]*x58;
                    const T x60 = PBy*x55 - PGy[g]*x56;
                    const T x61 = PBy*x56 - PGy[g]*x58;
                    const T x62 = x55 - x56;
                    const T x63 = PBy*x44 - PGy[g]*x46;
                    const T x64 = PAy*x29 - PGy[g]*x31;
                    const T x65 = PAy*x31 - PGy[g]*x34;
                    const T x66 = tpi*x44 + PBy*x64 - PGy[g]*x65;
                    const T x67 = PAy*x34 - PGy[g]*x38;
                    const T x68 = tpi*x46 + PBy*x65 - PGy[g]*x67;
                    const T x69 = x64 - x65;
                    const T x70 = tpi_2*x9 + PBz*x64 - PGz[g]*x65;
                    const T x71 = tpi_2*x14 + PBz*x65 - PGz[g]*x67;
                    const T x72 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x73 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x74 = x72 - x73;
                    const T x75 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x76 = x73 - x75;
                    const T x77 = x74 - x76;
                    const T x78 = tpi*x77 + PBz*x20 - PGz[g]*x23;
                    const T x79 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x80 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x81 = x79 - x80;
                    const T x82 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x83 = x80 - x82;
                    const T x84 = x81 - x83;
                    const T x85 = tpi*x84 + PBz*x9 - PGz[g]*x14;
                    const T x86 = PAx*x81 - PGx[g]*x83;
                    const T x87 = tpi*x86 + PBz*x5 - PGz[g]*x8;
                    const T x88 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x89 = x82 - x88;
                    const T x90 = PAx*x83 - PGx[g]*x89;
                    const T x91 = tpi*x90 + PBz*x8 - PGz[g]*x13;
                    const T x92 = PBz*x86 - PGz[g]*x90;
                    const T x93 = tpi*x28 + PAy*x79 - PGy[g]*x80;
                    const T x94 = tpi*x30 + PAy*x80 - PGy[g]*x82;
                    const T x95 = tpi_2*x81 + PAy*x93 - PGy[g]*x94;
                    const T x96 = tpi*x33 + PAy*x82 - PGy[g]*x88;
                    const T x97 = tpi_2*x83 + PAy*x94 - PGy[g]*x96;
                    const T x98 = PBx*x95 - PGx[g]*x97;
                    const T x99 = PAy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x100 = tpi*x37 + PAy*x88 - PGy[g]*x99;
                    const T x101 = tpi_2*x89 + PAy*x96 - PGy[g]*x100;
                    const T x102 = PBx*x97 - PGx[g]*x101;
                    const T x103 = x95 - x97;
                    const T x104 = PAx*x72 + tpi*x28 - PGx[g]*x73;
                    const T x105 = PAx*x73 + tpi*x30 - PGx[g]*x75;
                    const T x106 = PAx*x104 + tpi_2*x74 - PGx[g]*x105;
                    const T x107 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x108 = PAx*x75 + tpi*x33 - PGx[g]*x107;
                    const T x109 = PAx*x105 + tpi_2*x76 - PGx[g]*x108;
                    const T x110 = x106 - x109;
                    const T x111 = PBy*x106 - PGy[g]*x109;
                    const T x112 = PAx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x113 = PAx*x107 + tpi*x37 - PGx[g]*x112;
                    const T x114 = x75 - x107;
                    const T x115 = PAx*x108 + tpi_2*x114 - PGx[g]*x113;
                    const T x116 = PBy*x109 - PGy[g]*x115;
                    const T x117 = PBz*x106 - PGz[g]*x109;
                    const T x118 = PBz*x109 - PGz[g]*x115;
                    const T x119 = x104 - x105;
                    const T x120 = PBx*x106 + tpi_3*x119 - PGx[g]*x109;
                    const T x121 = x105 - x108;
                    const T x122 = PBx*x109 + tpi_3*x121 - PGx[g]*x115;
                    const T x123 = PBx*x119 + tpi_2*x77 - PGx[g]*x121;
                    const T x124 = PBz*x95 - PGz[g]*x97;
                    const T x125 = PBz*x97 - PGz[g]*x101;
                    const T x126 = x93 - x94;
                    const T x127 = tpi_3*x126 + PBy*x95 - PGy[g]*x97;
                    const T x128 = x94 - x96;
                    const T x129 = tpi_3*x128 + PBy*x97 - PGy[g]*x101;
                    const T x130 = tpi_2*x84 + PBy*x126 - PGy[g]*x128;
                    const T x131 = PBx*x93 - PGx[g]*x94;
                    const T x132 = PBx*x94 - PGx[g]*x96;
                    const T x133 = PBx*x131 + tpi*x126 - PGx[g]*x132;
                    const T x134 = PBx*x96 - PGx[g]*x100;
                    const T x135 = PBx*x132 + tpi*x128 - PGx[g]*x134;
                    const T x136 = x131 - x132;
                    const T x137 = PAz*x93 - PGz[g]*x94;
                    const T x138 = PAz*x94 - PGz[g]*x96;
                    const T x139 = tpi_2*x9 + PBy*x137 - PGy[g]*x138;
                    const T x140 = PAz*x96 - PGz[g]*x100;
                    const T x141 = tpi_2*x14 + PBy*x138 - PGy[g]*x140;
                    const T x142 = x137 - x138;
                    const T x143 = PBz*x126 - PGz[g]*x128;
                    const T x144 = tpi*x126 + PBz*x137 - PGz[g]*x138;
                    const T x145 = tpi*x128 + PBz*x138 - PGz[g]*x140;
                    const T x146 = PAx*x93 - PGx[g]*x94;
                    const T x147 = PAx*x94 - PGx[g]*x96;
                    const T x148 = tpi_2*x86 + PBy*x146 - PGy[g]*x147;
                    const T x149 = PAx*x96 - PGx[g]*x100;
                    const T x150 = tpi_2*x90 + PBy*x147 - PGy[g]*x149;
                    const T x151 = tpi*x77 + PBy*x86 - PGy[g]*x90;
                    const T x152 = x146 - x147;
                    const T x153 = PBz*x146 - PGz[g]*x147;
                    const T x154 = PBz*x147 - PGz[g]*x149;
                    const T x155 = PAy*x104 - PGy[g]*x105;
                    const T x156 = PAy*x105 - PGy[g]*x108;
                    const T x157 = PBx*x155 + tpi_2*x86 - PGx[g]*x156;
                    const T x158 = PAy*x108 - PGy[g]*x113;
                    const T x159 = PBx*x156 + tpi_2*x90 - PGx[g]*x158;
                    const T x160 = PBx*x86 + tpi*x84 - PGx[g]*x90;
                    const T x161 = x155 - x156;
                    const T x162 = PAz*x104 - PGz[g]*x105;
                    const T x163 = PAz*x105 - PGz[g]*x108;
                    const T x164 = PBx*x162 + tpi_2*x20 - PGx[g]*x163;
                    const T x165 = PAz*x108 - PGz[g]*x113;
                    const T x166 = PBx*x163 + tpi_2*x23 - PGx[g]*x165;
                    const T x167 = PBx*x20 + tpi*x26 - PGx[g]*x23;
                    const T x168 = x162 - x163;
                    const T x169 = PBz*x155 - PGz[g]*x156;
                    const T x170 = PBz*x156 - PGz[g]*x158;
                    const T x171 = tpi*x119 + PBy*x155 - PGy[g]*x156;
                    const T x172 = tpi*x121 + PBy*x156 - PGy[g]*x158;
                    const T x173 = PBy*x119 - PGy[g]*x121;
                    const T x174 = PBy*x162 - PGy[g]*x163;
                    const T x175 = PBy*x163 - PGy[g]*x165;
                    const T x176 = tpi*x119 + PBz*x162 - PGz[g]*x163;
                    const T x177 = tpi*x121 + PBz*x163 - PGz[g]*x165;
                    const T x178 = PBz*x119 - PGz[g]*x121;

                    integral_0[g] = PBx*x120 + tpi_3*x123 + tpi*x110 - PGx[g]*x122;
                    integral_1[g] = PBy*x120 - PGy[g]*x122;
                    integral_2[g] = PBz*x120 - PGz[g]*x122;
                    integral_3[g] = tpi*x110 + PBy*x111 - PGy[g]*x116;
                    integral_4[g] = PBy*x117 - PGy[g]*x118;
                    integral_5[g] = tpi*x110 + PBz*x117 - PGz[g]*x118;
                    integral_6[g] = PBx*x157 + tpi_2*x160 + tpi*x161 - PGx[g]*x159;
                    integral_7[g] = PBx*x171 + tpi_2*x151 - PGx[g]*x172;
                    integral_8[g] = PBz*x157 - PGz[g]*x159;
                    integral_9[g] = tpi*x161 + PBy*x171 + tpi*x173 - PGy[g]*x172;
                    integral_10[g] = PBz*x171 - PGz[g]*x172;
                    integral_11[g] = tpi*x161 + PBz*x169 - PGz[g]*x170;
                    integral_12[g] = PBx*x164 + tpi_2*x167 + tpi*x168 - PGx[g]*x166;
                    integral_13[g] = PBy*x164 - PGy[g]*x166;
                    integral_14[g] = PBx*x176 + tpi_2*x78 - PGx[g]*x177;
                    integral_15[g] = tpi*x168 + PBy*x174 - PGy[g]*x175;
                    integral_16[g] = PBy*x176 - PGy[g]*x177;
                    integral_17[g] = tpi*x168 + PBz*x176 + tpi*x178 - PGz[g]*x177;
                    integral_18[g] = PAx*x133 + tpi_2*x136 - PGx[g]*x135;
                    integral_19[g] = PBx*x148 + tpi*x130 - PGx[g]*x150;
                    integral_20[g] = PBx*x153 + tpi*x143 - PGx[g]*x154;
                    integral_21[g] = tpi_2*x151 + PBy*x148 + tpi*x152 - PGy[g]*x150;
                    integral_22[g] = PBz*x148 - PGz[g]*x150;
                    integral_23[g] = PBz*x153 + tpi*x152 - PGz[g]*x154;
                    integral_24[g] = PBx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                    integral_25[g] = PBx*x21 + tpi*x27 - PGx[g]*x24;
                    integral_26[g] = PBx*x87 + tpi*x85 - PGx[g]*x91;
                    integral_27[g] = tpi*x17 + PBy*x21 + tpi*x25 - PGy[g]*x24;
                    integral_28[g] = tpi*x78 + PBy*x87 - PGy[g]*x91;
                    integral_29[g] = tpi*x17 + PBz*x87 + tpi*x92 - PGz[g]*x91;
                    integral_30[g] = PAx*x51 + tpi_2*x54 - PGx[g]*x53;
                    integral_31[g] = PBx*x60 + tpi*x63 - PGx[g]*x61;
                    integral_32[g] = PBx*x57 + tpi*x48 - PGx[g]*x59;
                    integral_33[g] = PBy*x60 + tpi*x62 - PGy[g]*x61;
                    integral_34[g] = PBy*x57 - PGy[g]*x59;
                    integral_35[g] = tpi_2*x78 + PBz*x57 + tpi*x62 - PGz[g]*x59;
                    integral_36[g] = PBx*x98 + tpi*x103 - PGx[g]*x102;
                    integral_37[g] = PBx*x127 - PGx[g]*x129;
                    integral_38[g] = PBx*x124 - PGx[g]*x125;
                    integral_39[g] = tpi_3*x130 + tpi*x103 + PBy*x127 - PGy[g]*x129;
                    integral_40[g] = PBz*x127 - PGz[g]*x129;
                    integral_41[g] = tpi*x103 + PBz*x124 - PGz[g]*x125;
                    integral_42[g] = PAz*x133 - PGz[g]*x135;
                    integral_43[g] = PBx*x139 - PGx[g]*x141;
                    integral_44[g] = PBx*x144 - PGx[g]*x145;
                    integral_45[g] = tpi_2*x27 + PBy*x139 + tpi*x142 - PGy[g]*x141;
                    integral_46[g] = tpi_2*x85 + PBy*x144 - PGy[g]*x145;
                    integral_47[g] = tpi*x143 + PBz*x144 + tpi*x142 - PGz[g]*x145;
                    integral_48[g] = PAy*x51 - PGy[g]*x53;
                    integral_49[g] = PBx*x66 - PGx[g]*x68;
                    integral_50[g] = PBx*x70 - PGx[g]*x71;
                    integral_51[g] = tpi*x63 + PBy*x66 + tpi*x69 - PGy[g]*x68;
                    integral_52[g] = tpi*x48 + PBy*x70 - PGy[g]*x71;
                    integral_53[g] = tpi_2*x85 + PBz*x70 + tpi*x69 - PGz[g]*x71;
                    integral_54[g] = PBx*x36 + tpi*x41 - PGx[g]*x40;
                    integral_55[g] = PBx*x42 - PGx[g]*x43;
                    integral_56[g] = PBx*x45 - PGx[g]*x47;
                    integral_57[g] = tpi*x41 + PBy*x42 - PGy[g]*x43;
                    integral_58[g] = PBy*x45 - PGy[g]*x47;
                    integral_59[g] = tpi_3*x48 + tpi*x41 + PBz*x45 - PGz[g]*x47;

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
                    const T x10 = PBx*x5 + tpi*x9 - PGx[g]*x8;
                    const T x11 = PAz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x12 = PAy*x6 - PGy[g]*x11;
                    const T x13 = PAx*x7 - PGx[g]*x12;
                    const T x14 = x4 - x7;
                    const T x15 = PBx*x8 + tpi*x14 - PGx[g]*x13;
                    const T x16 = PBx*x9 - PGx[g]*x14;
                    const T x17 = x5 - x8;
                    const T x18 = x0 - x1;
                    const T x19 = x1 - x3;
                    const T x20 = PAx*x18 - PGx[g]*x19;
                    const T x21 = tpi*x20 + PBy*x5 - PGy[g]*x8;
                    const T x22 = x3 - x6;
                    const T x23 = PAx*x19 - PGx[g]*x22;
                    const T x24 = tpi*x23 + PBy*x8 - PGy[g]*x13;
                    const T x25 = PBy*x20 - PGy[g]*x23;
                    const T x26 = x18 - x19;
                    const T x27 = tpi*x26 + PBy*x9 - PGy[g]*x14;
                    const T x28 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x29 = tpi*x28 + PAz*x0 - PGz[g]*x1;
                    const T x30 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x31 = tpi*x30 + PAz*x1 - PGz[g]*x3;
                    const T x32 = tpi_2*x18 + PAz*x29 - PGz[g]*x31;
                    const T x33 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x34 = tpi*x33 + PAz*x3 - PGz[g]*x6;
                    const T x35 = tpi_2*x19 + PAz*x31 - PGz[g]*x34;
                    const T x36 = PBx*x32 - PGx[g]*x35;
                    const T x37 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x38 = tpi*x37 + PAz*x6 - PGz[g]*x11;
                    const T x39 = tpi_2*x22 + PAz*x34 - PGz[g]*x38;
                    const T x40 = PBx*x35 - PGx[g]*x39;
                    const T x41 = x32 - x35;
                    const T x42 = PBy*x32 - PGy[g]*x35;
                    const T x43 = PBy*x35 - PGy[g]*x39;
                    const T x44 = x29 - x31;
                    const T x45 = tpi_3*x44 + PBz*x32 - PGz[g]*x35;
                    const T x46 = x31 - x34;
                    const T x47 = tpi_3*x46 + PBz*x35 - PGz[g]*x39;
                    const T x48 = tpi_2*x26 + PBz*x44 - PGz[g]*x46;
                    const T x49 = PBx*x29 - PGx[g]*x31;
                    const T x50 = PBx*x31 - PGx[g]*x34;
                    const T x51 = PBx*x49 + tpi*x44 - PGx[g]*x50;
                    const T x52 = PBx*x34 - PGx[g]*x38;
                    const T x53 = PBx*x50 + tpi*x46 - PGx[g]*x52;
                    const T x54 = x49 - x50;
                    const T x55 = PAx*x29 - PGx[g]*x31;
                    const T x56 = PAx*x31 - PGx[g]*x34;
                    const T x57 = tpi_2*x20 + PBz*x55 - PGz[g]*x56;
                    const T x58 = PAx*x34 - PGx[g]*x38;
                    const T x59 = tpi_2*x23 + PBz*x56 - PGz[g]*x58;
                    const T x60 = PBy*x55 - PGy[g]*x56;
                    const T x61 = PBy*x56 - PGy[g]*x58;
                    const T x62 = x55 - x56;
                    const T x63 = PBy*x44 - PGy[g]*x46;
                    const T x64 = PAy*x29 - PGy[g]*x31;
                    const T x65 = PAy*x31 - PGy[g]*x34;
                    const T x66 = tpi*x44 + PBy*x64 - PGy[g]*x65;
                    const T x67 = PAy*x34 - PGy[g]*x38;
                    const T x68 = tpi*x46 + PBy*x65 - PGy[g]*x67;
                    const T x69 = x64 - x65;
                    const T x70 = tpi_2*x9 + PBz*x64 - PGz[g]*x65;
                    const T x71 = tpi_2*x14 + PBz*x65 - PGz[g]*x67;
                    const T x72 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x73 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x74 = x72 - x73;
                    const T x75 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x76 = x73 - x75;
                    const T x77 = x74 - x76;
                    const T x78 = tpi*x77 + PBz*x20 - PGz[g]*x23;
                    const T x79 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x80 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x81 = x79 - x80;
                    const T x82 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x83 = x80 - x82;
                    const T x84 = x81 - x83;
                    const T x85 = tpi*x84 + PBz*x9 - PGz[g]*x14;
                    const T x86 = PAx*x81 - PGx[g]*x83;
                    const T x87 = tpi*x86 + PBz*x5 - PGz[g]*x8;
                    const T x88 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x89 = x82 - x88;
                    const T x90 = PAx*x83 - PGx[g]*x89;
                    const T x91 = tpi*x90 + PBz*x8 - PGz[g]*x13;
                    const T x92 = PBz*x86 - PGz[g]*x90;
                    const T x93 = tpi*x28 + PAy*x79 - PGy[g]*x80;
                    const T x94 = tpi*x30 + PAy*x80 - PGy[g]*x82;
                    const T x95 = tpi_2*x81 + PAy*x93 - PGy[g]*x94;
                    const T x96 = tpi*x33 + PAy*x82 - PGy[g]*x88;
                    const T x97 = tpi_2*x83 + PAy*x94 - PGy[g]*x96;
                    const T x98 = PBx*x95 - PGx[g]*x97;
                    const T x99 = PAy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x100 = tpi*x37 + PAy*x88 - PGy[g]*x99;
                    const T x101 = tpi_2*x89 + PAy*x96 - PGy[g]*x100;
                    const T x102 = PBx*x97 - PGx[g]*x101;
                    const T x103 = x95 - x97;
                    const T x104 = PAx*x72 + tpi*x28 - PGx[g]*x73;
                    const T x105 = PAx*x73 + tpi*x30 - PGx[g]*x75;
                    const T x106 = PAx*x104 + tpi_2*x74 - PGx[g]*x105;
                    const T x107 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x108 = PAx*x75 + tpi*x33 - PGx[g]*x107;
                    const T x109 = PAx*x105 + tpi_2*x76 - PGx[g]*x108;
                    const T x110 = x106 - x109;
                    const T x111 = PBy*x106 - PGy[g]*x109;
                    const T x112 = PAx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x113 = PAx*x107 + tpi*x37 - PGx[g]*x112;
                    const T x114 = x75 - x107;
                    const T x115 = PAx*x108 + tpi_2*x114 - PGx[g]*x113;
                    const T x116 = PBy*x109 - PGy[g]*x115;
                    const T x117 = PBz*x106 - PGz[g]*x109;
                    const T x118 = PBz*x109 - PGz[g]*x115;
                    const T x119 = x104 - x105;
                    const T x120 = PBx*x106 + tpi_3*x119 - PGx[g]*x109;
                    const T x121 = x105 - x108;
                    const T x122 = PBx*x109 + tpi_3*x121 - PGx[g]*x115;
                    const T x123 = PBx*x119 + tpi_2*x77 - PGx[g]*x121;
                    const T x124 = PBz*x95 - PGz[g]*x97;
                    const T x125 = PBz*x97 - PGz[g]*x101;
                    const T x126 = x93 - x94;
                    const T x127 = tpi_3*x126 + PBy*x95 - PGy[g]*x97;
                    const T x128 = x94 - x96;
                    const T x129 = tpi_3*x128 + PBy*x97 - PGy[g]*x101;
                    const T x130 = tpi_2*x84 + PBy*x126 - PGy[g]*x128;
                    const T x131 = PBx*x93 - PGx[g]*x94;
                    const T x132 = PBx*x94 - PGx[g]*x96;
                    const T x133 = PBx*x131 + tpi*x126 - PGx[g]*x132;
                    const T x134 = PBx*x96 - PGx[g]*x100;
                    const T x135 = PBx*x132 + tpi*x128 - PGx[g]*x134;
                    const T x136 = x131 - x132;
                    const T x137 = PAz*x93 - PGz[g]*x94;
                    const T x138 = PAz*x94 - PGz[g]*x96;
                    const T x139 = tpi_2*x9 + PBy*x137 - PGy[g]*x138;
                    const T x140 = PAz*x96 - PGz[g]*x100;
                    const T x141 = tpi_2*x14 + PBy*x138 - PGy[g]*x140;
                    const T x142 = x137 - x138;
                    const T x143 = PBz*x126 - PGz[g]*x128;
                    const T x144 = tpi*x126 + PBz*x137 - PGz[g]*x138;
                    const T x145 = tpi*x128 + PBz*x138 - PGz[g]*x140;
                    const T x146 = PAx*x93 - PGx[g]*x94;
                    const T x147 = PAx*x94 - PGx[g]*x96;
                    const T x148 = tpi_2*x86 + PBy*x146 - PGy[g]*x147;
                    const T x149 = PAx*x96 - PGx[g]*x100;
                    const T x150 = tpi_2*x90 + PBy*x147 - PGy[g]*x149;
                    const T x151 = tpi*x77 + PBy*x86 - PGy[g]*x90;
                    const T x152 = x146 - x147;
                    const T x153 = PBz*x146 - PGz[g]*x147;
                    const T x154 = PBz*x147 - PGz[g]*x149;
                    const T x155 = PAy*x104 - PGy[g]*x105;
                    const T x156 = PAy*x105 - PGy[g]*x108;
                    const T x157 = PBx*x155 + tpi_2*x86 - PGx[g]*x156;
                    const T x158 = PAy*x108 - PGy[g]*x113;
                    const T x159 = PBx*x156 + tpi_2*x90 - PGx[g]*x158;
                    const T x160 = PBx*x86 + tpi*x84 - PGx[g]*x90;
                    const T x161 = x155 - x156;
                    const T x162 = PAz*x104 - PGz[g]*x105;
                    const T x163 = PAz*x105 - PGz[g]*x108;
                    const T x164 = PBx*x162 + tpi_2*x20 - PGx[g]*x163;
                    const T x165 = PAz*x108 - PGz[g]*x113;
                    const T x166 = PBx*x163 + tpi_2*x23 - PGx[g]*x165;
                    const T x167 = PBx*x20 + tpi*x26 - PGx[g]*x23;
                    const T x168 = x162 - x163;
                    const T x169 = PBz*x155 - PGz[g]*x156;
                    const T x170 = PBz*x156 - PGz[g]*x158;
                    const T x171 = tpi*x119 + PBy*x155 - PGy[g]*x156;
                    const T x172 = tpi*x121 + PBy*x156 - PGy[g]*x158;
                    const T x173 = PBy*x119 - PGy[g]*x121;
                    const T x174 = PBy*x162 - PGy[g]*x163;
                    const T x175 = PBy*x163 - PGy[g]*x165;
                    const T x176 = tpi*x119 + PBz*x162 - PGz[g]*x163;
                    const T x177 = tpi*x121 + PBz*x163 - PGz[g]*x165;
                    const T x178 = PBz*x119 - PGz[g]*x121;

                    integral_0[g] += PBx*x120 + tpi_3*x123 + tpi*x110 - PGx[g]*x122;
                    integral_1[g] += PBy*x120 - PGy[g]*x122;
                    integral_2[g] += PBz*x120 - PGz[g]*x122;
                    integral_3[g] += tpi*x110 + PBy*x111 - PGy[g]*x116;
                    integral_4[g] += PBy*x117 - PGy[g]*x118;
                    integral_5[g] += tpi*x110 + PBz*x117 - PGz[g]*x118;
                    integral_6[g] += PBx*x157 + tpi_2*x160 + tpi*x161 - PGx[g]*x159;
                    integral_7[g] += PBx*x171 + tpi_2*x151 - PGx[g]*x172;
                    integral_8[g] += PBz*x157 - PGz[g]*x159;
                    integral_9[g] += tpi*x161 + PBy*x171 + tpi*x173 - PGy[g]*x172;
                    integral_10[g] += PBz*x171 - PGz[g]*x172;
                    integral_11[g] += tpi*x161 + PBz*x169 - PGz[g]*x170;
                    integral_12[g] += PBx*x164 + tpi_2*x167 + tpi*x168 - PGx[g]*x166;
                    integral_13[g] += PBy*x164 - PGy[g]*x166;
                    integral_14[g] += PBx*x176 + tpi_2*x78 - PGx[g]*x177;
                    integral_15[g] += tpi*x168 + PBy*x174 - PGy[g]*x175;
                    integral_16[g] += PBy*x176 - PGy[g]*x177;
                    integral_17[g] += tpi*x168 + PBz*x176 + tpi*x178 - PGz[g]*x177;
                    integral_18[g] += PAx*x133 + tpi_2*x136 - PGx[g]*x135;
                    integral_19[g] += PBx*x148 + tpi*x130 - PGx[g]*x150;
                    integral_20[g] += PBx*x153 + tpi*x143 - PGx[g]*x154;
                    integral_21[g] += tpi_2*x151 + PBy*x148 + tpi*x152 - PGy[g]*x150;
                    integral_22[g] += PBz*x148 - PGz[g]*x150;
                    integral_23[g] += PBz*x153 + tpi*x152 - PGz[g]*x154;
                    integral_24[g] += PBx*x10 + tpi*x16 + tpi*x17 - PGx[g]*x15;
                    integral_25[g] += PBx*x21 + tpi*x27 - PGx[g]*x24;
                    integral_26[g] += PBx*x87 + tpi*x85 - PGx[g]*x91;
                    integral_27[g] += tpi*x17 + PBy*x21 + tpi*x25 - PGy[g]*x24;
                    integral_28[g] += tpi*x78 + PBy*x87 - PGy[g]*x91;
                    integral_29[g] += tpi*x17 + PBz*x87 + tpi*x92 - PGz[g]*x91;
                    integral_30[g] += PAx*x51 + tpi_2*x54 - PGx[g]*x53;
                    integral_31[g] += PBx*x60 + tpi*x63 - PGx[g]*x61;
                    integral_32[g] += PBx*x57 + tpi*x48 - PGx[g]*x59;
                    integral_33[g] += PBy*x60 + tpi*x62 - PGy[g]*x61;
                    integral_34[g] += PBy*x57 - PGy[g]*x59;
                    integral_35[g] += tpi_2*x78 + PBz*x57 + tpi*x62 - PGz[g]*x59;
                    integral_36[g] += PBx*x98 + tpi*x103 - PGx[g]*x102;
                    integral_37[g] += PBx*x127 - PGx[g]*x129;
                    integral_38[g] += PBx*x124 - PGx[g]*x125;
                    integral_39[g] += tpi_3*x130 + tpi*x103 + PBy*x127 - PGy[g]*x129;
                    integral_40[g] += PBz*x127 - PGz[g]*x129;
                    integral_41[g] += tpi*x103 + PBz*x124 - PGz[g]*x125;
                    integral_42[g] += PAz*x133 - PGz[g]*x135;
                    integral_43[g] += PBx*x139 - PGx[g]*x141;
                    integral_44[g] += PBx*x144 - PGx[g]*x145;
                    integral_45[g] += tpi_2*x27 + PBy*x139 + tpi*x142 - PGy[g]*x141;
                    integral_46[g] += tpi_2*x85 + PBy*x144 - PGy[g]*x145;
                    integral_47[g] += tpi*x143 + PBz*x144 + tpi*x142 - PGz[g]*x145;
                    integral_48[g] += PAy*x51 - PGy[g]*x53;
                    integral_49[g] += PBx*x66 - PGx[g]*x68;
                    integral_50[g] += PBx*x70 - PGx[g]*x71;
                    integral_51[g] += tpi*x63 + PBy*x66 + tpi*x69 - PGy[g]*x68;
                    integral_52[g] += tpi*x48 + PBy*x70 - PGy[g]*x71;
                    integral_53[g] += tpi_2*x85 + PBz*x70 + tpi*x69 - PGz[g]*x71;
                    integral_54[g] += PBx*x36 + tpi*x41 - PGx[g]*x40;
                    integral_55[g] += PBx*x42 - PGx[g]*x43;
                    integral_56[g] += PBx*x45 - PGx[g]*x47;
                    integral_57[g] += tpi*x41 + PBy*x42 - PGy[g]*x43;
                    integral_58[g] += PBy*x45 - PGy[g]*x47;
                    integral_59[g] += tpi_3*x48 + tpi*x41 + PBz*x45 - PGz[g]*x47;

                } // for g
            } // if/else prim
        } // for prim

    }
} // VRR_F_D

template void VRR_F_D<float>(
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

template void VRR_F_D<double>(
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

