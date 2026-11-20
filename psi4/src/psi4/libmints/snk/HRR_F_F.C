#include "HRR_F_F.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_F_F(
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
    alignas(64) T base_integral[74*g_block];
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
    T* __restrict base_integral_31 = &base_integral[31*g_block];
    T* __restrict base_integral_32 = &base_integral[32*g_block];
    T* __restrict base_integral_33 = &base_integral[33*g_block];
    T* __restrict base_integral_34 = &base_integral[34*g_block];
    T* __restrict base_integral_35 = &base_integral[35*g_block];
    T* __restrict base_integral_36 = &base_integral[36*g_block];
    T* __restrict base_integral_37 = &base_integral[37*g_block];
    T* __restrict base_integral_38 = &base_integral[38*g_block];
    T* __restrict base_integral_39 = &base_integral[39*g_block];
    T* __restrict base_integral_40 = &base_integral[40*g_block];
    T* __restrict base_integral_41 = &base_integral[41*g_block];
    T* __restrict base_integral_42 = &base_integral[42*g_block];
    T* __restrict base_integral_43 = &base_integral[43*g_block];
    T* __restrict base_integral_44 = &base_integral[44*g_block];
    T* __restrict base_integral_45 = &base_integral[45*g_block];
    T* __restrict base_integral_46 = &base_integral[46*g_block];
    T* __restrict base_integral_47 = &base_integral[47*g_block];
    T* __restrict base_integral_48 = &base_integral[48*g_block];
    T* __restrict base_integral_49 = &base_integral[49*g_block];
    T* __restrict base_integral_50 = &base_integral[50*g_block];
    T* __restrict base_integral_51 = &base_integral[51*g_block];
    T* __restrict base_integral_52 = &base_integral[52*g_block];
    T* __restrict base_integral_53 = &base_integral[53*g_block];
    T* __restrict base_integral_54 = &base_integral[54*g_block];
    T* __restrict base_integral_55 = &base_integral[55*g_block];
    T* __restrict base_integral_56 = &base_integral[56*g_block];
    T* __restrict base_integral_57 = &base_integral[57*g_block];
    T* __restrict base_integral_58 = &base_integral[58*g_block];
    T* __restrict base_integral_59 = &base_integral[59*g_block];
    T* __restrict base_integral_60 = &base_integral[60*g_block];
    T* __restrict base_integral_61 = &base_integral[61*g_block];
    T* __restrict base_integral_62 = &base_integral[62*g_block];
    T* __restrict base_integral_63 = &base_integral[63*g_block];
    T* __restrict base_integral_64 = &base_integral[64*g_block];
    T* __restrict base_integral_65 = &base_integral[65*g_block];
    T* __restrict base_integral_66 = &base_integral[66*g_block];
    T* __restrict base_integral_67 = &base_integral[67*g_block];
    T* __restrict base_integral_68 = &base_integral[68*g_block];
    T* __restrict base_integral_69 = &base_integral[69*g_block];
    T* __restrict base_integral_70 = &base_integral[70*g_block];
    T* __restrict base_integral_71 = &base_integral[71*g_block];
    T* __restrict base_integral_72 = &base_integral[72*g_block];
    T* __restrict base_integral_73 = &base_integral[73*g_block];
    alignas(64) T block[100*g_block];
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
    T* __restrict block_60 = &block[60*g_block];
    T* __restrict block_61 = &block[61*g_block];
    T* __restrict block_62 = &block[62*g_block];
    T* __restrict block_63 = &block[63*g_block];
    T* __restrict block_64 = &block[64*g_block];
    T* __restrict block_65 = &block[65*g_block];
    T* __restrict block_66 = &block[66*g_block];
    T* __restrict block_67 = &block[67*g_block];
    T* __restrict block_68 = &block[68*g_block];
    T* __restrict block_69 = &block[69*g_block];
    T* __restrict block_70 = &block[70*g_block];
    T* __restrict block_71 = &block[71*g_block];
    T* __restrict block_72 = &block[72*g_block];
    T* __restrict block_73 = &block[73*g_block];
    T* __restrict block_74 = &block[74*g_block];
    T* __restrict block_75 = &block[75*g_block];
    T* __restrict block_76 = &block[76*g_block];
    T* __restrict block_77 = &block[77*g_block];
    T* __restrict block_78 = &block[78*g_block];
    T* __restrict block_79 = &block[79*g_block];
    T* __restrict block_80 = &block[80*g_block];
    T* __restrict block_81 = &block[81*g_block];
    T* __restrict block_82 = &block[82*g_block];
    T* __restrict block_83 = &block[83*g_block];
    T* __restrict block_84 = &block[84*g_block];
    T* __restrict block_85 = &block[85*g_block];
    T* __restrict block_86 = &block[86*g_block];
    T* __restrict block_87 = &block[87*g_block];
    T* __restrict block_88 = &block[88*g_block];
    T* __restrict block_89 = &block[89*g_block];
    T* __restrict block_90 = &block[90*g_block];
    T* __restrict block_91 = &block[91*g_block];
    T* __restrict block_92 = &block[92*g_block];
    T* __restrict block_93 = &block[93*g_block];
    T* __restrict block_94 = &block[94*g_block];
    T* __restrict block_95 = &block[95*g_block];
    T* __restrict block_96 = &block[96*g_block];
    T* __restrict block_97 = &block[97*g_block];
    T* __restrict block_98 = &block[98*g_block];
    T* __restrict block_99 = &block[99*g_block];

    // fewer gridpoints than one block: work on a zero padded copy of the inputs
    const bool padded = (gb_size < g_block);
    alignas(64) T pad[10*g_block];
    if (padded) {
        for (int k = 0; k < 10*g_block; ++k) pad[k] = 0;
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
                    pad[3*g_block+k] = v_scaled_boys[prim*7*gb_size+0*gb_size+k];
                    pad[4*g_block+k] = v_scaled_boys[prim*7*gb_size+1*gb_size+k];
                    pad[5*g_block+k] = v_scaled_boys[prim*7*gb_size+2*gb_size+k];
                    pad[6*g_block+k] = v_scaled_boys[prim*7*gb_size+3*gb_size+k];
                    pad[7*g_block+k] = v_scaled_boys[prim*7*gb_size+4*gb_size+k];
                    pad[8*g_block+k] = v_scaled_boys[prim*7*gb_size+5*gb_size+k];
                    pad[9*g_block+k] = v_scaled_boys[prim*7*gb_size+6*gb_size+k];
                }
            }
            const T* __restrict PGx = padded ? &pad[0*g_block] : &v_PGx[prim*gb_size + g_lo];
            const T* __restrict PGy = padded ? &pad[1*g_block] : &v_PGy[prim*gb_size + g_lo];
            const T* __restrict PGz = padded ? &pad[2*g_block] : &v_PGz[prim*gb_size + g_lo];
            const T* __restrict scaled_boys_0 = padded ? &pad[3*g_block] : &v_scaled_boys[prim*7*gb_size+0*gb_size + g_lo];
            const T* __restrict scaled_boys_1 = padded ? &pad[4*g_block] : &v_scaled_boys[prim*7*gb_size+1*gb_size + g_lo];
            const T* __restrict scaled_boys_2 = padded ? &pad[5*g_block] : &v_scaled_boys[prim*7*gb_size+2*gb_size + g_lo];
            const T* __restrict scaled_boys_3 = padded ? &pad[6*g_block] : &v_scaled_boys[prim*7*gb_size+3*gb_size + g_lo];
            const T* __restrict scaled_boys_4 = padded ? &pad[7*g_block] : &v_scaled_boys[prim*7*gb_size+4*gb_size + g_lo];
            const T* __restrict scaled_boys_5 = padded ? &pad[8*g_block] : &v_scaled_boys[prim*7*gb_size+5*gb_size + g_lo];
            const T* __restrict scaled_boys_6 = padded ? &pad[9*g_block] : &v_scaled_boys[prim*7*gb_size+6*gb_size + g_lo];

            const T PAx = v_PAx[prim];
            const T PAy = v_PAy[prim];
            const T PAz = v_PAz[prim];
            const T tpi = 0.5 / v_p[prim];
            const T tpi_2 = 2*tpi;
            const T tpi_3 = 3*tpi;
            const T tpi_4 = 4*tpi;
            const T tpi_5 = 5*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x2 = PAx*x0 - PGx[g]*x1;
                    const T x3 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x4 = PAx*x1 - PGx[g]*x3;
                    base_integral_4[g] = PAz*x2 - PGz[g]*x4;
                    const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x6 = tpi*x5 + PAy*x0 - PGy[g]*x1;
                    const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x8 = tpi*x7 + PAy*x1 - PGy[g]*x3;
                    base_integral_7[g] = PAz*x6 - PGz[g]*x8;
                    const T x9 = x0 - x1;
                    const T x10 = PAy*x6 + tpi_2*x9 - PGy[g]*x8;
                    base_integral_6[g] = x10;
                    const T x11 = PAx*x6 - PGx[g]*x8;
                    base_integral_3[g] = x11;
                    const T x12 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x13 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x14 = PAx*x12 + tpi*x5 - PGx[g]*x13;
                    const T x15 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x16 = PAx*x13 + tpi*x7 - PGx[g]*x15;
                    base_integral_1[g] = PAy*x14 - PGy[g]*x16;
                    const T x17 = x12 - x13;
                    const T x18 = PAx*x14 + tpi_2*x17 - PGx[g]*x16;
                    base_integral_0[g] = x18;
                    const T x19 = PAz*x14 - PGz[g]*x16;
                    base_integral_2[g] = x19;
                    const T x20 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x21 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x22 = tpi*x5 + PAz*x20 - PGz[g]*x21;
                    const T x23 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x24 = tpi*x7 + PAz*x21 - PGz[g]*x23;
                    base_integral_5[g] = PAx*x22 - PGx[g]*x24;
                    const T x25 = x20 - x21;
                    const T x26 = PAz*x22 + tpi_2*x25 - PGz[g]*x24;
                    base_integral_9[g] = x26;
                    const T x27 = PAy*x22 - PGy[g]*x24;
                    base_integral_8[g] = x27;
                    const T x28 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x29 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x30 = PAx*x15 + tpi*x29 - PGx[g]*x28;
                    const T x31 = x13 - x15;
                    const T x32 = PAx*x16 + tpi_2*x31 - PGx[g]*x30;
                    base_integral_11[g] = PAy*x18 - PGy[g]*x32;
                    const T x33 = PAz*x18 - PGz[g]*x32;
                    base_integral_12[g] = x33;
                    const T x34 = x14 - x16;
                    const T x35 = PAx*x18 + tpi_3*x34 - PGx[g]*x32;
                    base_integral_10[g] = x35;
                    const T x36 = PAz*x16 - PGz[g]*x30;
                    base_integral_14[g] = PAy*x19 - PGy[g]*x36;
                    const T x37 = PAz*x19 + tpi*x34 - PGz[g]*x36;
                    base_integral_15[g] = x37;
                    const T x38 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x39 = PAz*x23 + tpi*x29 - PGz[g]*x38;
                    const T x40 = x21 - x23;
                    const T x41 = PAz*x24 + tpi_2*x40 - PGz[g]*x39;
                    base_integral_19[g] = PAx*x26 - PGx[g]*x41;
                    const T x42 = PAy*x26 - PGy[g]*x41;
                    base_integral_23[g] = x42;
                    const T x43 = PAy*x24 - PGy[g]*x39;
                    base_integral_18[g] = PAx*x27 - PGx[g]*x43;
                    const T x44 = x22 - x24;
                    base_integral_22[g] = PAy*x27 + tpi*x44 - PGy[g]*x43;
                    const T x45 = PAz*x26 + tpi_3*x44 - PGz[g]*x41;
                    base_integral_24[g] = x45;
                    const T x46 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x47 = PAy*x3 + tpi*x29 - PGy[g]*x46;
                    const T x48 = x1 - x3;
                    const T x49 = PAy*x8 + tpi_2*x48 - PGy[g]*x47;
                    const T x50 = PAx*x10 - PGx[g]*x49;
                    base_integral_16[g] = x50;
                    const T x51 = PAz*x10 - PGz[g]*x49;
                    base_integral_21[g] = x51;
                    const T x52 = x6 - x8;
                    const T x53 = PAy*x10 + tpi_3*x52 - PGy[g]*x49;
                    base_integral_20[g] = x53;
                    const T x54 = PAx*x8 - PGx[g]*x47;
                    const T x55 = PAx*x11 + tpi*x52 - PGx[g]*x54;
                    base_integral_13[g] = x55;
                    const T x56 = PAz*x11 - PGz[g]*x54;
                    base_integral_17[g] = x56;
                    const T x57 = x11 - x54;
                    const T x58 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x59 = PAy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x60 = PAy*x46 + tpi*x58 - PGy[g]*x59;
                    const T x61 = PAx*x47 - PGx[g]*x60;
                    const T x62 = PAz*x54 - PGz[g]*x61;
                    base_integral_37[g] = PAz*x56 + tpi*x57 - PGz[g]*x62;
                    const T x63 = x8 - x47;
                    const T x64 = PAx*x54 + tpi*x63 - PGx[g]*x61;
                    base_integral_28[g] = PAx*x55 + tpi_2*x57 - PGx[g]*x64;
                    base_integral_32[g] = PAz*x55 - PGz[g]*x64;
                    const T x65 = x3 - x46;
                    const T x66 = PAy*x47 + tpi_2*x65 - PGy[g]*x60;
                    const T x67 = PAy*x49 + tpi_3*x63 - PGy[g]*x66;
                    base_integral_35[g] = PAx*x53 - PGx[g]*x67;
                    const T x68 = PAz*x53 - PGz[g]*x67;
                    base_integral_41[g] = x68;
                    const T x69 = PAz*x49 - PGz[g]*x66;
                    base_integral_36[g] = PAx*x51 - PGx[g]*x69;
                    const T x70 = x10 - x49;
                    const T x71 = PAy*x53 + tpi_4*x70 - PGy[g]*x67;
                    base_integral_40[g] = x71;
                    const T x72 = PAz*x51 + tpi*x70 - PGz[g]*x69;
                    base_integral_42[g] = x72;
                    const T x73 = PAx*x49 - PGx[g]*x66;
                    const T x74 = PAx*x50 + tpi*x70 - PGx[g]*x73;
                    base_integral_31[g] = x74;
                    const T x75 = scaled_boys_4[g] - scaled_boys_5[g];
                    const T x76 = PAy*scaled_boys_5[g] - PGy[g]*scaled_boys_6[g];
                    const T x77 = PAy*x59 + tpi*x75 - PGy[g]*x76;
                    const T x78 = x46 - x59;
                    const T x79 = PAy*x60 + tpi_2*x78 - PGy[g]*x77;
                    const T x80 = PAx*x66 - PGx[g]*x79;
                    const T x81 = x49 - x66;
                    const T x82 = PAx*x73 + tpi*x81 - PGx[g]*x80;
                    const T x83 = x50 - x73;
                    base_integral_52[g] = PAx*x74 + tpi_2*x83 - PGx[g]*x82;
                    const T x84 = x55 - x64;
                    base_integral_56[g] = PAy*x74 + tpi_3*x84 - PGy[g]*x82;
                    base_integral_57[g] = PAz*x74 - PGz[g]*x82;
                    const T x85 = PAz*x66 - PGz[g]*x79;
                    const T x86 = PAz*x69 + tpi*x81 - PGz[g]*x85;
                    const T x87 = x51 - x69;
                    base_integral_70[g] = PAz*x72 + tpi_2*x87 - PGz[g]*x86;
                    const T x88 = x47 - x60;
                    const T x89 = PAy*x66 + tpi_3*x88 - PGy[g]*x79;
                    const T x90 = PAy*x67 + tpi_4*x81 - PGy[g]*x89;
                    base_integral_61[g] = PAx*x71 - PGx[g]*x90;
                    base_integral_63[g] = PAx*x72 - PGx[g]*x86;
                    const T x91 = PAz*x67 - PGz[g]*x89;
                    const T x92 = x53 - x67;
                    base_integral_69[g] = PAz*x68 + tpi*x92 - PGz[g]*x91;
                    base_integral_62[g] = PAx*x68 - PGx[g]*x91;
                    base_integral_67[g] = PAy*x71 + tpi_5*x92 - PGy[g]*x90;
                    base_integral_68[g] = PAz*x71 - PGz[g]*x90;
                    const T x93 = PAz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x94 = PAz*x38 + tpi*x58 - PGz[g]*x93;
                    const T x95 = x23 - x38;
                    const T x96 = PAz*x39 + tpi_2*x95 - PGz[g]*x94;
                    const T x97 = x24 - x39;
                    const T x98 = PAz*x41 + tpi_3*x97 - PGz[g]*x96;
                    base_integral_39[g] = PAx*x45 - PGx[g]*x98;
                    const T x99 = PAy*x45 - PGy[g]*x98;
                    base_integral_44[g] = x99;
                    const T x100 = PAz*scaled_boys_5[g] - PGz[g]*scaled_boys_6[g];
                    const T x101 = PAz*x93 + tpi*x75 - PGz[g]*x100;
                    const T x102 = x38 - x93;
                    const T x103 = PAz*x94 + tpi_2*x102 - PGz[g]*x101;
                    const T x104 = x39 - x94;
                    const T x105 = PAz*x96 + tpi_3*x104 - PGz[g]*x103;
                    const T x106 = PAy*x98 - PGy[g]*x105;
                    base_integral_65[g] = PAx*x99 - PGx[g]*x106;
                    const T x107 = x45 - x98;
                    base_integral_71[g] = PAy*x99 + tpi*x107 - PGy[g]*x106;
                    const T x108 = x26 - x41;
                    const T x109 = PAz*x45 + tpi_4*x108 - PGz[g]*x98;
                    base_integral_45[g] = x109;
                    const T x110 = x41 - x96;
                    const T x111 = PAz*x98 + tpi_4*x110 - PGz[g]*x105;
                    base_integral_73[g] = PAz*x109 + tpi_5*x107 - PGz[g]*x111;
                    const T x112 = PAy*x41 - PGy[g]*x96;
                    const T x113 = PAy*x42 + tpi*x108 - PGy[g]*x112;
                    base_integral_43[g] = x113;
                    const T x114 = PAy*x96 - PGy[g]*x103;
                    const T x115 = PAy*x112 + tpi*x110 - PGy[g]*x114;
                    base_integral_64[g] = PAx*x113 - PGx[g]*x115;
                    base_integral_38[g] = PAx*x42 - PGx[g]*x112;
                    base_integral_66[g] = PAx*x109 - PGx[g]*x111;
                    base_integral_72[g] = PAy*x109 - PGy[g]*x111;
                    const T x116 = PAx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x117 = PAx*x28 + tpi*x58 - PGx[g]*x116;
                    const T x118 = x15 - x28;
                    const T x119 = PAx*x30 + tpi_2*x118 - PGx[g]*x117;
                    const T x120 = PAz*x32 - PGz[g]*x119;
                    const T x121 = PAy*x33 - PGy[g]*x120;
                    base_integral_29[g] = x121;
                    const T x122 = PAx*scaled_boys_5[g] - PGx[g]*scaled_boys_6[g];
                    const T x123 = PAx*x116 + tpi*x75 - PGx[g]*x122;
                    const T x124 = x28 - x116;
                    const T x125 = PAx*x117 + tpi_2*x124 - PGx[g]*x123;
                    const T x126 = PAz*x119 - PGz[g]*x125;
                    const T x127 = PAy*x120 - PGy[g]*x126;
                    const T x128 = x33 - x120;
                    base_integral_53[g] = PAy*x121 + tpi*x128 - PGy[g]*x127;
                    const T x129 = x18 - x32;
                    const T x130 = PAz*x33 + tpi*x129 - PGz[g]*x120;
                    base_integral_30[g] = x130;
                    const T x131 = x32 - x119;
                    const T x132 = PAz*x120 + tpi*x131 - PGz[g]*x126;
                    base_integral_55[g] = PAz*x130 + tpi_2*x128 - PGz[g]*x132;
                    base_integral_54[g] = PAy*x130 - PGy[g]*x132;
                    const T x133 = x16 - x30;
                    const T x134 = PAx*x32 + tpi_3*x133 - PGx[g]*x119;
                    const T x135 = PAx*x35 + tpi_4*x129 - PGx[g]*x134;
                    base_integral_25[g] = x135;
                    const T x136 = PAz*x30 - PGz[g]*x117;
                    const T x137 = PAz*x36 + tpi*x133 - PGz[g]*x136;
                    const T x138 = x19 - x36;
                    const T x139 = PAz*x37 + tpi_2*x138 - PGz[g]*x137;
                    base_integral_34[g] = x139;
                    const T x140 = x30 - x117;
                    const T x141 = PAz*x117 - PGz[g]*x123;
                    const T x142 = PAz*x136 + tpi*x140 - PGz[g]*x141;
                    const T x143 = x36 - x136;
                    const T x144 = PAz*x137 + tpi_2*x143 - PGz[g]*x142;
                    base_integral_59[g] = PAy*x139 - PGy[g]*x144;
                    const T x145 = PAx*x119 + tpi_3*x140 - PGx[g]*x125;
                    const T x146 = PAx*x134 + tpi_4*x131 - PGx[g]*x145;
                    base_integral_47[g] = PAy*x135 - PGy[g]*x146;
                    const T x147 = PAy*x37 - PGy[g]*x137;
                    base_integral_33[g] = x147;
                    const T x148 = x37 - x137;
                    const T x149 = PAy*x137 - PGy[g]*x142;
                    base_integral_58[g] = PAy*x147 + tpi*x148 - PGy[g]*x149;
                    base_integral_51[g] = PAx*x130 + tpi_3*x148 - PGx[g]*x132;
                    base_integral_60[g] = PAz*x139 + tpi_3*x148 - PGz[g]*x144;
                    base_integral_27[g] = PAz*x35 - PGz[g]*x134;
                    base_integral_48[g] = PAz*x135 - PGz[g]*x146;
                    const T x150 = x35 - x134;
                    base_integral_46[g] = PAx*x135 + tpi_5*x150 - PGx[g]*x146;
                    const T x151 = PAy*x35 - PGy[g]*x134;
                    base_integral_26[g] = x151;
                    const T x152 = PAy*x134 - PGy[g]*x145;
                    base_integral_49[g] = PAy*x151 + tpi*x150 - PGy[g]*x152;
                    base_integral_50[g] = PAz*x151 - PGz[g]*x152;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x1 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x2 = PAx*x0 - PGx[g]*x1;
                    const T x3 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x4 = PAx*x1 - PGx[g]*x3;
                    base_integral_4[g] += PAz*x2 - PGz[g]*x4;
                    const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x6 = tpi*x5 + PAy*x0 - PGy[g]*x1;
                    const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x8 = tpi*x7 + PAy*x1 - PGy[g]*x3;
                    base_integral_7[g] += PAz*x6 - PGz[g]*x8;
                    const T x9 = x0 - x1;
                    const T x10 = PAy*x6 + tpi_2*x9 - PGy[g]*x8;
                    base_integral_6[g] += x10;
                    const T x11 = PAx*x6 - PGx[g]*x8;
                    base_integral_3[g] += x11;
                    const T x12 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x13 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x14 = PAx*x12 + tpi*x5 - PGx[g]*x13;
                    const T x15 = PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x16 = PAx*x13 + tpi*x7 - PGx[g]*x15;
                    base_integral_1[g] += PAy*x14 - PGy[g]*x16;
                    const T x17 = x12 - x13;
                    const T x18 = PAx*x14 + tpi_2*x17 - PGx[g]*x16;
                    base_integral_0[g] += x18;
                    const T x19 = PAz*x14 - PGz[g]*x16;
                    base_integral_2[g] += x19;
                    const T x20 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x21 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x22 = tpi*x5 + PAz*x20 - PGz[g]*x21;
                    const T x23 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x24 = tpi*x7 + PAz*x21 - PGz[g]*x23;
                    base_integral_5[g] += PAx*x22 - PGx[g]*x24;
                    const T x25 = x20 - x21;
                    const T x26 = PAz*x22 + tpi_2*x25 - PGz[g]*x24;
                    base_integral_9[g] += x26;
                    const T x27 = PAy*x22 - PGy[g]*x24;
                    base_integral_8[g] += x27;
                    const T x28 = PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x29 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x30 = PAx*x15 + tpi*x29 - PGx[g]*x28;
                    const T x31 = x13 - x15;
                    const T x32 = PAx*x16 + tpi_2*x31 - PGx[g]*x30;
                    base_integral_11[g] += PAy*x18 - PGy[g]*x32;
                    const T x33 = PAz*x18 - PGz[g]*x32;
                    base_integral_12[g] += x33;
                    const T x34 = x14 - x16;
                    const T x35 = PAx*x18 + tpi_3*x34 - PGx[g]*x32;
                    base_integral_10[g] += x35;
                    const T x36 = PAz*x16 - PGz[g]*x30;
                    base_integral_14[g] += PAy*x19 - PGy[g]*x36;
                    const T x37 = PAz*x19 + tpi*x34 - PGz[g]*x36;
                    base_integral_15[g] += x37;
                    const T x38 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x39 = PAz*x23 + tpi*x29 - PGz[g]*x38;
                    const T x40 = x21 - x23;
                    const T x41 = PAz*x24 + tpi_2*x40 - PGz[g]*x39;
                    base_integral_19[g] += PAx*x26 - PGx[g]*x41;
                    const T x42 = PAy*x26 - PGy[g]*x41;
                    base_integral_23[g] += x42;
                    const T x43 = PAy*x24 - PGy[g]*x39;
                    base_integral_18[g] += PAx*x27 - PGx[g]*x43;
                    const T x44 = x22 - x24;
                    base_integral_22[g] += PAy*x27 + tpi*x44 - PGy[g]*x43;
                    const T x45 = PAz*x26 + tpi_3*x44 - PGz[g]*x41;
                    base_integral_24[g] += x45;
                    const T x46 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x47 = PAy*x3 + tpi*x29 - PGy[g]*x46;
                    const T x48 = x1 - x3;
                    const T x49 = PAy*x8 + tpi_2*x48 - PGy[g]*x47;
                    const T x50 = PAx*x10 - PGx[g]*x49;
                    base_integral_16[g] += x50;
                    const T x51 = PAz*x10 - PGz[g]*x49;
                    base_integral_21[g] += x51;
                    const T x52 = x6 - x8;
                    const T x53 = PAy*x10 + tpi_3*x52 - PGy[g]*x49;
                    base_integral_20[g] += x53;
                    const T x54 = PAx*x8 - PGx[g]*x47;
                    const T x55 = PAx*x11 + tpi*x52 - PGx[g]*x54;
                    base_integral_13[g] += x55;
                    const T x56 = PAz*x11 - PGz[g]*x54;
                    base_integral_17[g] += x56;
                    const T x57 = x11 - x54;
                    const T x58 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x59 = PAy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x60 = PAy*x46 + tpi*x58 - PGy[g]*x59;
                    const T x61 = PAx*x47 - PGx[g]*x60;
                    const T x62 = PAz*x54 - PGz[g]*x61;
                    base_integral_37[g] += PAz*x56 + tpi*x57 - PGz[g]*x62;
                    const T x63 = x8 - x47;
                    const T x64 = PAx*x54 + tpi*x63 - PGx[g]*x61;
                    base_integral_28[g] += PAx*x55 + tpi_2*x57 - PGx[g]*x64;
                    base_integral_32[g] += PAz*x55 - PGz[g]*x64;
                    const T x65 = x3 - x46;
                    const T x66 = PAy*x47 + tpi_2*x65 - PGy[g]*x60;
                    const T x67 = PAy*x49 + tpi_3*x63 - PGy[g]*x66;
                    base_integral_35[g] += PAx*x53 - PGx[g]*x67;
                    const T x68 = PAz*x53 - PGz[g]*x67;
                    base_integral_41[g] += x68;
                    const T x69 = PAz*x49 - PGz[g]*x66;
                    base_integral_36[g] += PAx*x51 - PGx[g]*x69;
                    const T x70 = x10 - x49;
                    const T x71 = PAy*x53 + tpi_4*x70 - PGy[g]*x67;
                    base_integral_40[g] += x71;
                    const T x72 = PAz*x51 + tpi*x70 - PGz[g]*x69;
                    base_integral_42[g] += x72;
                    const T x73 = PAx*x49 - PGx[g]*x66;
                    const T x74 = PAx*x50 + tpi*x70 - PGx[g]*x73;
                    base_integral_31[g] += x74;
                    const T x75 = scaled_boys_4[g] - scaled_boys_5[g];
                    const T x76 = PAy*scaled_boys_5[g] - PGy[g]*scaled_boys_6[g];
                    const T x77 = PAy*x59 + tpi*x75 - PGy[g]*x76;
                    const T x78 = x46 - x59;
                    const T x79 = PAy*x60 + tpi_2*x78 - PGy[g]*x77;
                    const T x80 = PAx*x66 - PGx[g]*x79;
                    const T x81 = x49 - x66;
                    const T x82 = PAx*x73 + tpi*x81 - PGx[g]*x80;
                    const T x83 = x50 - x73;
                    base_integral_52[g] += PAx*x74 + tpi_2*x83 - PGx[g]*x82;
                    const T x84 = x55 - x64;
                    base_integral_56[g] += PAy*x74 + tpi_3*x84 - PGy[g]*x82;
                    base_integral_57[g] += PAz*x74 - PGz[g]*x82;
                    const T x85 = PAz*x66 - PGz[g]*x79;
                    const T x86 = PAz*x69 + tpi*x81 - PGz[g]*x85;
                    const T x87 = x51 - x69;
                    base_integral_70[g] += PAz*x72 + tpi_2*x87 - PGz[g]*x86;
                    const T x88 = x47 - x60;
                    const T x89 = PAy*x66 + tpi_3*x88 - PGy[g]*x79;
                    const T x90 = PAy*x67 + tpi_4*x81 - PGy[g]*x89;
                    base_integral_61[g] += PAx*x71 - PGx[g]*x90;
                    base_integral_63[g] += PAx*x72 - PGx[g]*x86;
                    const T x91 = PAz*x67 - PGz[g]*x89;
                    const T x92 = x53 - x67;
                    base_integral_69[g] += PAz*x68 + tpi*x92 - PGz[g]*x91;
                    base_integral_62[g] += PAx*x68 - PGx[g]*x91;
                    base_integral_67[g] += PAy*x71 + tpi_5*x92 - PGy[g]*x90;
                    base_integral_68[g] += PAz*x71 - PGz[g]*x90;
                    const T x93 = PAz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x94 = PAz*x38 + tpi*x58 - PGz[g]*x93;
                    const T x95 = x23 - x38;
                    const T x96 = PAz*x39 + tpi_2*x95 - PGz[g]*x94;
                    const T x97 = x24 - x39;
                    const T x98 = PAz*x41 + tpi_3*x97 - PGz[g]*x96;
                    base_integral_39[g] += PAx*x45 - PGx[g]*x98;
                    const T x99 = PAy*x45 - PGy[g]*x98;
                    base_integral_44[g] += x99;
                    const T x100 = PAz*scaled_boys_5[g] - PGz[g]*scaled_boys_6[g];
                    const T x101 = PAz*x93 + tpi*x75 - PGz[g]*x100;
                    const T x102 = x38 - x93;
                    const T x103 = PAz*x94 + tpi_2*x102 - PGz[g]*x101;
                    const T x104 = x39 - x94;
                    const T x105 = PAz*x96 + tpi_3*x104 - PGz[g]*x103;
                    const T x106 = PAy*x98 - PGy[g]*x105;
                    base_integral_65[g] += PAx*x99 - PGx[g]*x106;
                    const T x107 = x45 - x98;
                    base_integral_71[g] += PAy*x99 + tpi*x107 - PGy[g]*x106;
                    const T x108 = x26 - x41;
                    const T x109 = PAz*x45 + tpi_4*x108 - PGz[g]*x98;
                    base_integral_45[g] += x109;
                    const T x110 = x41 - x96;
                    const T x111 = PAz*x98 + tpi_4*x110 - PGz[g]*x105;
                    base_integral_73[g] += PAz*x109 + tpi_5*x107 - PGz[g]*x111;
                    const T x112 = PAy*x41 - PGy[g]*x96;
                    const T x113 = PAy*x42 + tpi*x108 - PGy[g]*x112;
                    base_integral_43[g] += x113;
                    const T x114 = PAy*x96 - PGy[g]*x103;
                    const T x115 = PAy*x112 + tpi*x110 - PGy[g]*x114;
                    base_integral_64[g] += PAx*x113 - PGx[g]*x115;
                    base_integral_38[g] += PAx*x42 - PGx[g]*x112;
                    base_integral_66[g] += PAx*x109 - PGx[g]*x111;
                    base_integral_72[g] += PAy*x109 - PGy[g]*x111;
                    const T x116 = PAx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x117 = PAx*x28 + tpi*x58 - PGx[g]*x116;
                    const T x118 = x15 - x28;
                    const T x119 = PAx*x30 + tpi_2*x118 - PGx[g]*x117;
                    const T x120 = PAz*x32 - PGz[g]*x119;
                    const T x121 = PAy*x33 - PGy[g]*x120;
                    base_integral_29[g] += x121;
                    const T x122 = PAx*scaled_boys_5[g] - PGx[g]*scaled_boys_6[g];
                    const T x123 = PAx*x116 + tpi*x75 - PGx[g]*x122;
                    const T x124 = x28 - x116;
                    const T x125 = PAx*x117 + tpi_2*x124 - PGx[g]*x123;
                    const T x126 = PAz*x119 - PGz[g]*x125;
                    const T x127 = PAy*x120 - PGy[g]*x126;
                    const T x128 = x33 - x120;
                    base_integral_53[g] += PAy*x121 + tpi*x128 - PGy[g]*x127;
                    const T x129 = x18 - x32;
                    const T x130 = PAz*x33 + tpi*x129 - PGz[g]*x120;
                    base_integral_30[g] += x130;
                    const T x131 = x32 - x119;
                    const T x132 = PAz*x120 + tpi*x131 - PGz[g]*x126;
                    base_integral_55[g] += PAz*x130 + tpi_2*x128 - PGz[g]*x132;
                    base_integral_54[g] += PAy*x130 - PGy[g]*x132;
                    const T x133 = x16 - x30;
                    const T x134 = PAx*x32 + tpi_3*x133 - PGx[g]*x119;
                    const T x135 = PAx*x35 + tpi_4*x129 - PGx[g]*x134;
                    base_integral_25[g] += x135;
                    const T x136 = PAz*x30 - PGz[g]*x117;
                    const T x137 = PAz*x36 + tpi*x133 - PGz[g]*x136;
                    const T x138 = x19 - x36;
                    const T x139 = PAz*x37 + tpi_2*x138 - PGz[g]*x137;
                    base_integral_34[g] += x139;
                    const T x140 = x30 - x117;
                    const T x141 = PAz*x117 - PGz[g]*x123;
                    const T x142 = PAz*x136 + tpi*x140 - PGz[g]*x141;
                    const T x143 = x36 - x136;
                    const T x144 = PAz*x137 + tpi_2*x143 - PGz[g]*x142;
                    base_integral_59[g] += PAy*x139 - PGy[g]*x144;
                    const T x145 = PAx*x119 + tpi_3*x140 - PGx[g]*x125;
                    const T x146 = PAx*x134 + tpi_4*x131 - PGx[g]*x145;
                    base_integral_47[g] += PAy*x135 - PGy[g]*x146;
                    const T x147 = PAy*x37 - PGy[g]*x137;
                    base_integral_33[g] += x147;
                    const T x148 = x37 - x137;
                    const T x149 = PAy*x137 - PGy[g]*x142;
                    base_integral_58[g] += PAy*x147 + tpi*x148 - PGy[g]*x149;
                    base_integral_51[g] += PAx*x130 + tpi_3*x148 - PGx[g]*x132;
                    base_integral_60[g] += PAz*x139 + tpi_3*x148 - PGz[g]*x144;
                    base_integral_27[g] += PAz*x35 - PGz[g]*x134;
                    base_integral_48[g] += PAz*x135 - PGz[g]*x146;
                    const T x150 = x35 - x134;
                    base_integral_46[g] += PAx*x135 + tpi_5*x150 - PGx[g]*x146;
                    const T x151 = PAy*x35 - PGy[g]*x134;
                    base_integral_26[g] += x151;
                    const T x152 = PAy*x134 - PGy[g]*x145;
                    base_integral_49[g] += PAy*x151 + tpi*x150 - PGy[g]*x152;
                    base_integral_50[g] += PAz*x151 - PGz[g]*x152;

                } // for g
            } // if/else prim
        } // for prim

        #pragma omp simd
        for (int g = 0; g < g_block; ++g) {

            const T x0 = base_integral_46[g] + ABx*base_integral_25[g];
            const T x1 = base_integral_25[g] + ABx*base_integral_10[g];
            const T x2 = x0 + ABx*x1;
            const T x3 = base_integral_10[g] + ABx*base_integral_0[g];
            const T x4 = x1 + ABx*x3;
            block_0[g] = x2 + ABx*x4;
            const T x5 = ABz*base_integral_9[g] + base_integral_24[g];
            const T x6 = ABz*base_integral_24[g] + base_integral_45[g];
            const T x7 = ABz*x5 + x6;
            const T x8 = ABz*base_integral_45[g] + base_integral_73[g];
            const T x9 = ABz*x6 + x8;
            block_99[g] = ABz*x7 + x9;
            const T x10 = ABz*base_integral_19[g] + base_integral_39[g];
            const T x11 = ABz*base_integral_39[g] + base_integral_66[g];
            const T x12 = ABz*x10 + x11;
            block_95[g] = x12 + ABx*x7;
            const T x13 = ABz*base_integral_5[g] + base_integral_19[g];
            const T x14 = x10 + ABz*x13;
            block_59[g] = ABz*x14 + x12;
            const T x15 = ABz*base_integral_34[g] + base_integral_60[g];
            const T x16 = x15 + ABx*x10;
            const T x17 = x10 + ABx*x5;
            block_92[g] = x16 + ABx*x17;
            const T x18 = ABz*base_integral_15[g] + base_integral_34[g];
            const T x19 = ABz*x18 + x15;
            block_55[g] = x19 + ABx*x14;
            const T x20 = ABz*base_integral_2[g] + base_integral_15[g];
            const T x21 = ABz*x20 + x18;
            block_29[g] = ABz*x21 + x19;
            const T x22 = base_integral_55[g] + ABx*base_integral_34[g];
            const T x23 = base_integral_34[g] + ABx*base_integral_19[g];
            const T x24 = x22 + ABx*x23;
            const T x25 = base_integral_19[g] + ABx*base_integral_9[g];
            const T x26 = x23 + ABx*x25;
            block_90[g] = x24 + ABx*x26;
            const T x27 = base_integral_30[g] + ABx*base_integral_15[g];
            const T x28 = base_integral_15[g] + ABx*base_integral_5[g];
            const T x29 = x27 + ABx*x28;
            block_52[g] = ABz*x29 + x24;
            const T x30 = base_integral_51[g] + ABx*base_integral_30[g];
            const T x31 = x30 + ABx*x27;
            block_50[g] = x31 + ABx*x29;
            const T x32 = ABz*base_integral_12[g] + base_integral_30[g];
            const T x33 = ABz*base_integral_30[g] + base_integral_55[g];
            const T x34 = ABz*x32 + x33;
            block_25[g] = x34 + ABx*x21;
            const T x35 = ABz*base_integral_0[g] + base_integral_12[g];
            const T x36 = ABz*x35 + x32;
            block_9[g] = ABz*x36 + x34;
            const T x37 = ABz*base_integral_10[g] + base_integral_27[g];
            const T x38 = ABz*base_integral_27[g] + base_integral_51[g];
            const T x39 = ABz*x37 + x38;
            block_5[g] = x39 + ABx*x36;
            const T x40 = base_integral_27[g] + ABx*base_integral_12[g];
            const T x41 = base_integral_12[g] + ABx*base_integral_2[g];
            const T x42 = x40 + ABx*x41;
            block_22[g] = ABz*x42 + x31;
            const T x43 = base_integral_48[g] + ABx*base_integral_27[g];
            const T x44 = x43 + ABx*x40;
            block_2[g] = ABz*x4 + x44;
            block_20[g] = x44 + ABx*x42;
            const T x45 = ABz*base_integral_23[g] + base_integral_44[g];
            const T x46 = ABz*base_integral_44[g] + base_integral_72[g];
            const T x47 = ABz*x45 + x46;
            block_98[g] = x47 + ABy*x7;
            const T x48 = ABz*base_integral_38[g] + base_integral_65[g];
            const T x49 = x48 + ABy*x10;
            const T x50 = x45 + ABy*x5;
            block_94[g] = x49 + ABx*x50;
            const T x51 = ABz*base_integral_18[g] + base_integral_38[g];
            const T x52 = ABz*x51 + x48;
            block_58[g] = x52 + ABy*x14;
            const T x53 = ABz*base_integral_43[g] + base_integral_71[g];
            const T x54 = x53 + ABy*x45;
            block_97[g] = ABy*x50 + x54;
            const T x55 = ABz*base_integral_8[g] + base_integral_23[g];
            const T x56 = ABz*x55 + x45;
            block_89[g] = ABz*x56 + x47;
            block_85[g] = x52 + ABx*x56;
            const T x57 = ABz*base_integral_22[g] + base_integral_43[g];
            const T x58 = ABz*x57 + x53;
            block_88[g] = x58 + ABy*x56;
            const T x59 = base_integral_59[g] + ABy*base_integral_34[g];
            const T x60 = base_integral_38[g] + ABy*base_integral_19[g];
            const T x61 = x59 + ABx*x60;
            const T x62 = base_integral_23[g] + ABy*base_integral_9[g];
            const T x63 = x60 + ABx*x62;
            block_91[g] = x61 + ABx*x63;
            const T x64 = base_integral_64[g] + ABy*base_integral_38[g];
            const T x65 = ABy*x60 + x64;
            const T x66 = base_integral_43[g] + ABy*base_integral_23[g];
            const T x67 = x66 + ABy*x62;
            block_93[g] = x65 + ABx*x67;
            const T x68 = base_integral_70[g] + ABy*base_integral_43[g];
            const T x69 = x68 + ABy*x66;
            block_96[g] = x69 + ABy*x67;
            const T x70 = base_integral_33[g] + ABx*base_integral_18[g];
            const T x71 = x70 + ABy*x28;
            block_54[g] = ABz*x71 + x61;
            const T x72 = ABz*base_integral_33[g] + base_integral_59[g];
            const T x73 = x72 + ABx*x51;
            const T x74 = x51 + ABx*x55;
            block_82[g] = x73 + ABx*x74;
            const T x75 = ABz*base_integral_37[g] + base_integral_64[g];
            const T x76 = x75 + ABx*x57;
            block_84[g] = ABy*x74 + x76;
            const T x77 = base_integral_37[g] + ABy*base_integral_18[g];
            const T x78 = base_integral_18[g] + ABy*base_integral_5[g];
            const T x79 = x77 + ABy*x78;
            block_57[g] = ABz*x79 + x65;
            const T x80 = ABz*base_integral_14[g] + base_integral_33[g];
            const T x81 = ABz*x80 + x72;
            block_28[g] = x81 + ABy*x21;
            const T x82 = ABz*base_integral_4[g] + base_integral_18[g];
            const T x83 = ABz*x82 + x51;
            block_49[g] = ABz*x83 + x52;
            block_45[g] = x81 + ABx*x83;
            const T x84 = base_integral_54[g] + ABx*base_integral_33[g];
            const T x85 = x84 + ABx*x70;
            const T x86 = base_integral_18[g] + ABx*base_integral_8[g];
            const T x87 = x70 + ABx*x86;
            block_80[g] = x85 + ABx*x87;
            block_51[g] = x85 + ABy*x29;
            const T x88 = ABy*base_integral_33[g] + base_integral_58[g];
            const T x89 = x88 + ABx*x77;
            block_53[g] = x89 + ABy*x71;
            const T x90 = base_integral_22[g] + ABy*base_integral_8[g];
            const T x91 = x77 + ABx*x90;
            block_81[g] = x89 + ABx*x91;
            const T x92 = ABy*base_integral_37[g] + base_integral_63[g];
            const T x93 = ABy*x77 + x92;
            block_56[g] = x93 + ABy*x79;
            const T x94 = ABz*base_integral_17[g] + base_integral_37[g];
            const T x95 = ABz*x94 + x75;
            block_48[g] = x95 + ABy*x83;
            const T x96 = ABz*base_integral_29[g] + base_integral_54[g];
            const T x97 = ABy*x32 + x96;
            const T x98 = x80 + ABy*x20;
            block_24[g] = x97 + ABx*x98;
            const T x99 = ABz*base_integral_11[g] + base_integral_29[g];
            const T x100 = ABz*x99 + x96;
            block_8[g] = ABy*x36 + x100;
            const T x101 = ABz*base_integral_26[g] + base_integral_50[g];
            const T x102 = x101 + ABy*x37;
            const T x103 = x99 + ABy*x35;
            block_4[g] = x102 + ABx*x103;
            const T x104 = base_integral_50[g] + ABx*base_integral_29[g];
            const T x105 = base_integral_29[g] + ABx*base_integral_14[g];
            const T x106 = x104 + ABx*x105;
            block_21[g] = x106 + ABy*x42;
            const T x107 = base_integral_14[g] + ABx*base_integral_4[g];
            const T x108 = x105 + ABx*x107;
            block_42[g] = ABz*x108 + x85;
            block_40[g] = x106 + ABx*x108;
            const T x109 = ABz*base_integral_32[g] + base_integral_58[g];
            const T x110 = ABy*x80 + x109;
            block_27[g] = ABy*x98 + x110;
            const T x111 = x94 + ABy*x82;
            block_44[g] = x110 + ABx*x111;
            const T x112 = ABz*base_integral_1[g] + base_integral_14[g];
            const T x113 = ABz*x112 + x80;
            block_15[g] = x100 + ABx*x113;
            block_19[g] = ABz*x113 + x81;
            const T x114 = ABy*base_integral_29[g] + base_integral_53[g];
            const T x115 = ABy*base_integral_12[g] + base_integral_29[g];
            const T x116 = x114 + ABy*x115;
            const T x117 = ABy*base_integral_14[g] + base_integral_32[g];
            const T x118 = base_integral_14[g] + ABy*base_integral_2[g];
            const T x119 = x117 + ABy*x118;
            block_23[g] = x116 + ABx*x119;
            const T x120 = base_integral_53[g] + ABx*base_integral_32[g];
            const T x121 = base_integral_32[g] + ABx*base_integral_17[g];
            const T x122 = x120 + ABx*x121;
            block_41[g] = x122 + ABy*x108;
            const T x123 = base_integral_17[g] + ABx*base_integral_7[g];
            const T x124 = x121 + ABx*x123;
            block_70[g] = x122 + ABx*x124;
            const T x125 = ABz*base_integral_7[g] + base_integral_22[g];
            const T x126 = ABz*x125 + x57;
            block_79[g] = ABz*x126 + x58;
            block_75[g] = x95 + ABx*x126;
            const T x127 = x109 + ABx*x94;
            const T x128 = x94 + ABx*x125;
            block_72[g] = x127 + ABx*x128;
            const T x129 = base_integral_42[g] + ABy*base_integral_22[g];
            const T x130 = x129 + ABy*x90;
            block_83[g] = x93 + ABx*x130;
            block_87[g] = ABz*x130 + x69;
            const T x131 = ABy*base_integral_42[g] + base_integral_69[g];
            const T x132 = x131 + ABy*x129;
            block_86[g] = x132 + ABy*x130;
            const T x133 = ABz*base_integral_21[g] + base_integral_42[g];
            const T x134 = ABz*base_integral_42[g] + base_integral_70[g];
            const T x135 = ABz*x133 + x134;
            block_78[g] = x135 + ABy*x126;
            const T x136 = ABy*base_integral_11[g] + base_integral_28[g];
            const T x137 = ABz*x136 + x114;
            block_7[g] = ABy*x103 + x137;
            const T x138 = ABy*base_integral_26[g] + base_integral_49[g];
            const T x139 = ABy*base_integral_10[g] + base_integral_26[g];
            const T x140 = x138 + ABy*x139;
            const T x141 = ABy*base_integral_0[g] + base_integral_11[g];
            const T x142 = x136 + ABy*x141;
            block_3[g] = x140 + ABx*x142;
            const T x143 = ABy*base_integral_28[g] + base_integral_52[g];
            const T x144 = ABy*x136 + x143;
            block_6[g] = ABy*x142 + x144;
            const T x145 = base_integral_26[g] + ABx*base_integral_11[g];
            const T x146 = base_integral_11[g] + ABx*base_integral_1[g];
            const T x147 = x145 + ABx*x146;
            block_12[g] = ABz*x147 + x106;
            const T x148 = base_integral_47[g] + ABx*base_integral_26[g];
            const T x149 = x148 + ABx*x145;
            block_1[g] = ABy*x4 + x149;
            block_10[g] = x149 + ABx*x147;
            const T x150 = base_integral_49[g] + ABx*base_integral_28[g];
            const T x151 = base_integral_28[g] + ABx*base_integral_13[g];
            const T x152 = x150 + ABx*x151;
            block_11[g] = ABy*x147 + x152;
            const T x153 = ABz*base_integral_13[g] + base_integral_32[g];
            const T x154 = x153 + ABy*x112;
            block_14[g] = x137 + ABx*x154;
            const T x155 = ABz*x153 + x109;
            block_18[g] = ABy*x113 + x155;
            const T x156 = ABy*base_integral_32[g] + base_integral_57[g];
            const T x157 = ABy*x117 + x156;
            block_26[g] = ABy*x119 + x157;
            const T x158 = ABz*base_integral_31[g] + base_integral_57[g];
            const T x159 = ABy*x153 + x158;
            block_17[g] = ABy*x154 + x159;
            const T x160 = ABy*base_integral_13[g] + base_integral_31[g];
            const T x161 = ABy*base_integral_1[g] + base_integral_13[g];
            const T x162 = x160 + ABy*x161;
            block_13[g] = x144 + ABx*x162;
            const T x163 = base_integral_13[g] + ABx*base_integral_3[g];
            const T x164 = x151 + ABx*x163;
            block_30[g] = x152 + ABx*x164;
            block_32[g] = ABz*x164 + x122;
            const T x165 = ABy*base_integral_31[g] + base_integral_56[g];
            const T x166 = ABy*x160 + x165;
            block_16[g] = ABy*x162 + x166;
            const T x167 = base_integral_52[g] + ABx*base_integral_31[g];
            const T x168 = base_integral_31[g] + ABx*base_integral_16[g];
            const T x169 = x167 + ABx*x168;
            block_31[g] = ABy*x164 + x169;
            const T x170 = base_integral_16[g] + ABx*base_integral_6[g];
            const T x171 = x168 + ABx*x170;
            block_60[g] = x169 + ABx*x171;
            const T x172 = ABz*base_integral_3[g] + base_integral_17[g];
            const T x173 = ABz*x172 + x94;
            block_35[g] = x155 + ABx*x173;
            block_39[g] = ABz*x173 + x95;
            const T x174 = ABz*base_integral_36[g] + base_integral_63[g];
            const T x175 = ABy*x94 + x174;
            block_47[g] = ABy*x111 + x175;
            const T x176 = x133 + ABy*x125;
            block_74[g] = x175 + ABx*x176;
            const T x177 = ABy*base_integral_17[g] + base_integral_36[g];
            const T x178 = base_integral_17[g] + ABy*base_integral_4[g];
            const T x179 = x177 + ABy*x178;
            block_43[g] = x157 + ABx*x179;
            const T x180 = x156 + ABx*x177;
            const T x181 = base_integral_21[g] + ABy*base_integral_7[g];
            const T x182 = x177 + ABx*x181;
            block_71[g] = x180 + ABx*x182;
            const T x183 = ABz*base_integral_16[g] + base_integral_36[g];
            const T x184 = x183 + ABy*x172;
            block_34[g] = x159 + ABx*x184;
            const T x185 = ABz*x183 + x174;
            block_38[g] = ABy*x173 + x185;
            const T x186 = ABy*base_integral_36[g] + base_integral_62[g];
            const T x187 = x186 + ABy*x177;
            block_46[g] = ABy*x179 + x187;
            const T x188 = x158 + ABx*x183;
            const T x189 = ABz*base_integral_6[g] + base_integral_21[g];
            const T x190 = x183 + ABx*x189;
            block_62[g] = x188 + ABx*x190;
            const T x191 = ABz*x189 + x133;
            block_65[g] = x185 + ABx*x191;
            block_69[g] = ABz*x191 + x135;
            const T x192 = ABy*base_integral_21[g] + base_integral_41[g];
            const T x193 = x192 + ABy*x181;
            block_73[g] = x187 + ABx*x193;
            block_77[g] = ABz*x193 + x132;
            const T x194 = ABy*base_integral_41[g] + base_integral_68[g];
            const T x195 = x194 + ABy*x192;
            block_76[g] = x195 + ABy*x193;
            const T x196 = ABz*base_integral_20[g] + base_integral_41[g];
            const T x197 = ABz*base_integral_41[g] + base_integral_69[g];
            const T x198 = ABz*x196 + x197;
            block_68[g] = ABy*x191 + x198;
            const T x199 = ABy*base_integral_16[g] + base_integral_35[g];
            const T x200 = ABy*base_integral_3[g] + base_integral_16[g];
            const T x201 = x199 + ABy*x200;
            block_33[g] = x166 + ABx*x201;
            const T x202 = ABz*x199 + x186;
            const T x203 = ABy*x189 + x196;
            block_64[g] = x202 + ABx*x203;
            block_37[g] = ABy*x184 + x202;
            const T x204 = ABy*base_integral_35[g] + base_integral_61[g];
            const T x205 = ABy*x199 + x204;
            block_36[g] = ABy*x201 + x205;
            const T x206 = x165 + ABx*x199;
            const T x207 = ABy*base_integral_6[g] + base_integral_20[g];
            const T x208 = x199 + ABx*x207;
            block_61[g] = x206 + ABx*x208;
            const T x209 = ABy*base_integral_20[g] + base_integral_40[g];
            const T x210 = ABy*x207 + x209;
            const T x211 = ABy*base_integral_40[g] + base_integral_67[g];
            const T x212 = ABy*x209 + x211;
            block_66[g] = ABy*x210 + x212;
            block_63[g] = x205 + ABx*x210;
            block_67[g] = ABz*x210 + x195;

        } // for g

        // copy the block to the output arrays
        for (int i = 0; i < 100; ++i) {
            const T* __restrict src = &block[i*g_block];
            T* __restrict dst = &integrals[i*gb_size + g_lo];
            if (padded) {
                SNK_COPY(dst, src, gb_size);
            } else {
                SNK_COPY(dst, src, g_block);
            }
        }
    } // for g_lo
} // HRR_F_F

template void HRR_F_F<float>(
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

template void HRR_F_F<double>(
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

