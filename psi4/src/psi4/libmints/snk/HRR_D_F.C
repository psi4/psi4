#include "HRR_D_F.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_D_F(
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
    enum { g_block = (sizeof(T) <= 4) ? 16 : 16 };
    alignas(64) T base_integral[46*g_block];
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

            const T PBx = v_PAx[prim] + ABx;
            const T PBy = v_PAy[prim] + ABy;
            const T PBz = v_PAz[prim] + ABz;
            const T tpi = 0.5 / v_p[prim];
            const T tpi_2 = 2*tpi;
            const T tpi_3 = 3*tpi;
            const T tpi_4 = 4*tpi;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PBx*x0 - PGx[g]*x1;
                    const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PBx*x1 - PGx[g]*x3;
                    base_integral_4[g] = PBy*x2 - PGy[g]*x4;
                    const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x6 = tpi*x5 + PBz*x0 - PGz[g]*x1;
                    const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x8 = tpi*x7 + PBz*x1 - PGz[g]*x3;
                    base_integral_5[g] = PBx*x6 - PGx[g]*x8;
                    const T x9 = x0 - x1;
                    const T x10 = PBz*x6 + tpi_2*x9 - PGz[g]*x8;
                    base_integral_9[g] = x10;
                    const T x11 = PBy*x6 - PGy[g]*x8;
                    base_integral_8[g] = x11;
                    const T x12 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x13 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x14 = PBx*x12 + tpi*x5 - PGx[g]*x13;
                    const T x15 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x16 = PBx*x13 + tpi*x7 - PGx[g]*x15;
                    base_integral_1[g] = PBy*x14 - PGy[g]*x16;
                    const T x17 = x12 - x13;
                    const T x18 = PBx*x14 + tpi_2*x17 - PGx[g]*x16;
                    base_integral_0[g] = x18;
                    const T x19 = PBz*x14 - PGz[g]*x16;
                    base_integral_2[g] = x19;
                    const T x20 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x21 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x22 = tpi*x5 + PBy*x20 - PGy[g]*x21;
                    const T x23 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x24 = tpi*x7 + PBy*x21 - PGy[g]*x23;
                    base_integral_7[g] = PBz*x22 - PGz[g]*x24;
                    const T x25 = x20 - x21;
                    const T x26 = PBy*x22 + tpi_2*x25 - PGy[g]*x24;
                    base_integral_6[g] = x26;
                    const T x27 = PBx*x22 - PGx[g]*x24;
                    base_integral_3[g] = x27;
                    const T x28 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x29 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x30 = PBx*x15 + tpi*x29 - PGx[g]*x28;
                    const T x31 = x13 - x15;
                    const T x32 = PBx*x16 + tpi_2*x31 - PGx[g]*x30;
                    base_integral_11[g] = PBy*x18 - PGy[g]*x32;
                    const T x33 = PBz*x18 - PGz[g]*x32;
                    base_integral_12[g] = x33;
                    const T x34 = x14 - x16;
                    const T x35 = PBx*x18 + tpi_3*x34 - PGx[g]*x32;
                    base_integral_10[g] = x35;
                    const T x36 = PBz*x16 - PGz[g]*x30;
                    base_integral_15[g] = PBz*x19 + tpi*x34 - PGz[g]*x36;
                    base_integral_14[g] = PBy*x19 - PGy[g]*x36;
                    const T x37 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x38 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x39 = PBx*x28 + tpi*x38 - PGx[g]*x37;
                    const T x40 = x15 - x28;
                    const T x41 = PBx*x30 + tpi_2*x40 - PGx[g]*x39;
                    const T x42 = x16 - x30;
                    const T x43 = PBx*x32 + tpi_3*x42 - PGx[g]*x41;
                    base_integral_26[g] = PBy*x35 - PGy[g]*x43;
                    const T x44 = x18 - x32;
                    const T x45 = PBz*x32 - PGz[g]*x41;
                    base_integral_30[g] = PBz*x33 + tpi*x44 - PGz[g]*x45;
                    base_integral_29[g] = PBy*x33 - PGy[g]*x45;
                    base_integral_25[g] = PBx*x35 + tpi_4*x44 - PGx[g]*x43;
                    base_integral_27[g] = PBz*x35 - PGz[g]*x43;
                    const T x46 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x47 = PBy*x23 + tpi*x29 - PGy[g]*x46;
                    const T x48 = x21 - x23;
                    const T x49 = PBy*x24 + tpi_2*x48 - PGy[g]*x47;
                    base_integral_21[g] = PBz*x26 - PGz[g]*x49;
                    const T x50 = PBx*x26 - PGx[g]*x49;
                    base_integral_16[g] = x50;
                    const T x51 = PBx*x24 - PGx[g]*x47;
                    base_integral_17[g] = PBz*x27 - PGz[g]*x51;
                    const T x52 = x22 - x24;
                    const T x53 = PBx*x27 + tpi*x52 - PGx[g]*x51;
                    base_integral_13[g] = x53;
                    const T x54 = PBy*x26 + tpi_3*x52 - PGy[g]*x49;
                    base_integral_20[g] = x54;
                    const T x55 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x56 = PBy*x46 + tpi*x38 - PGy[g]*x55;
                    const T x57 = PBx*x47 - PGx[g]*x56;
                    const T x58 = x24 - x47;
                    const T x59 = PBx*x51 + tpi*x58 - PGx[g]*x57;
                    const T x60 = x27 - x51;
                    base_integral_28[g] = PBx*x53 + tpi_2*x60 - PGx[g]*x59;
                    const T x61 = x23 - x46;
                    const T x62 = PBy*x47 + tpi_2*x61 - PGy[g]*x56;
                    const T x63 = PBy*x49 + tpi_3*x58 - PGy[g]*x62;
                    base_integral_35[g] = PBx*x54 - PGx[g]*x63;
                    base_integral_32[g] = PBz*x53 - PGz[g]*x59;
                    const T x64 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x65 = PBz*x3 + tpi*x29 - PGz[g]*x64;
                    const T x66 = x1 - x3;
                    const T x67 = PBz*x8 + tpi_2*x66 - PGz[g]*x65;
                    const T x68 = PBx*x10 - PGx[g]*x67;
                    base_integral_19[g] = x68;
                    const T x69 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x70 = PBz*x64 + tpi*x38 - PGz[g]*x69;
                    const T x71 = x3 - x64;
                    const T x72 = PBz*x65 + tpi_2*x71 - PGz[g]*x70;
                    const T x73 = PBx*x67 - PGx[g]*x72;
                    const T x74 = x10 - x67;
                    base_integral_34[g] = PBx*x68 + tpi*x74 - PGx[g]*x73;
                    const T x75 = PBx*x49 - PGx[g]*x62;
                    const T x76 = x26 - x49;
                    base_integral_31[g] = PBx*x50 + tpi*x76 - PGx[g]*x75;
                    base_integral_36[g] = PBz*x50 - PGz[g]*x75;
                    base_integral_40[g] = PBy*x54 + tpi_4*x76 - PGy[g]*x63;
                    base_integral_41[g] = PBz*x54 - PGz[g]*x63;
                    const T x77 = PBy*x10 - PGy[g]*x67;
                    base_integral_23[g] = x77;
                    const T x78 = x6 - x8;
                    const T x79 = PBz*x10 + tpi_3*x78 - PGz[g]*x67;
                    base_integral_24[g] = x79;
                    const T x80 = PBy*x8 - PGy[g]*x65;
                    const T x81 = PBy*x11 + tpi*x78 - PGy[g]*x80;
                    base_integral_22[g] = x81;
                    const T x82 = PBy*x65 - PGy[g]*x70;
                    const T x83 = x8 - x65;
                    const T x84 = PBy*x80 + tpi*x83 - PGy[g]*x82;
                    base_integral_37[g] = PBx*x81 - PGx[g]*x84;
                    const T x85 = PBx*x11 - PGx[g]*x80;
                    base_integral_18[g] = x85;
                    const T x86 = PBx*x80 - PGx[g]*x82;
                    const T x87 = x11 - x80;
                    base_integral_33[g] = PBx*x85 + tpi*x87 - PGx[g]*x86;
                    base_integral_42[g] = PBy*x81 + tpi_2*x87 - PGy[g]*x84;
                    const T x88 = PBz*x67 + tpi_3*x83 - PGz[g]*x72;
                    base_integral_39[g] = PBx*x79 - PGx[g]*x88;
                    const T x89 = PBy*x67 - PGy[g]*x72;
                    base_integral_38[g] = PBx*x77 - PGx[g]*x89;
                    base_integral_43[g] = PBy*x77 + tpi*x74 - PGy[g]*x89;
                    base_integral_45[g] = PBz*x79 + tpi_4*x74 - PGz[g]*x88;
                    base_integral_44[g] = PBy*x79 - PGy[g]*x88;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                    const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                    const T x2 = PBx*x0 - PGx[g]*x1;
                    const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                    const T x4 = PBx*x1 - PGx[g]*x3;
                    base_integral_4[g] += PBy*x2 - PGy[g]*x4;
                    const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x6 = tpi*x5 + PBz*x0 - PGz[g]*x1;
                    const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x8 = tpi*x7 + PBz*x1 - PGz[g]*x3;
                    base_integral_5[g] += PBx*x6 - PGx[g]*x8;
                    const T x9 = x0 - x1;
                    const T x10 = PBz*x6 + tpi_2*x9 - PGz[g]*x8;
                    base_integral_9[g] += x10;
                    const T x11 = PBy*x6 - PGy[g]*x8;
                    base_integral_8[g] += x11;
                    const T x12 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                    const T x13 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                    const T x14 = PBx*x12 + tpi*x5 - PGx[g]*x13;
                    const T x15 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                    const T x16 = PBx*x13 + tpi*x7 - PGx[g]*x15;
                    base_integral_1[g] += PBy*x14 - PGy[g]*x16;
                    const T x17 = x12 - x13;
                    const T x18 = PBx*x14 + tpi_2*x17 - PGx[g]*x16;
                    base_integral_0[g] += x18;
                    const T x19 = PBz*x14 - PGz[g]*x16;
                    base_integral_2[g] += x19;
                    const T x20 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                    const T x21 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                    const T x22 = tpi*x5 + PBy*x20 - PGy[g]*x21;
                    const T x23 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                    const T x24 = tpi*x7 + PBy*x21 - PGy[g]*x23;
                    base_integral_7[g] += PBz*x22 - PGz[g]*x24;
                    const T x25 = x20 - x21;
                    const T x26 = PBy*x22 + tpi_2*x25 - PGy[g]*x24;
                    base_integral_6[g] += x26;
                    const T x27 = PBx*x22 - PGx[g]*x24;
                    base_integral_3[g] += x27;
                    const T x28 = PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g];
                    const T x29 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x30 = PBx*x15 + tpi*x29 - PGx[g]*x28;
                    const T x31 = x13 - x15;
                    const T x32 = PBx*x16 + tpi_2*x31 - PGx[g]*x30;
                    base_integral_11[g] += PBy*x18 - PGy[g]*x32;
                    const T x33 = PBz*x18 - PGz[g]*x32;
                    base_integral_12[g] += x33;
                    const T x34 = x14 - x16;
                    const T x35 = PBx*x18 + tpi_3*x34 - PGx[g]*x32;
                    base_integral_10[g] += x35;
                    const T x36 = PBz*x16 - PGz[g]*x30;
                    base_integral_15[g] += PBz*x19 + tpi*x34 - PGz[g]*x36;
                    base_integral_14[g] += PBy*x19 - PGy[g]*x36;
                    const T x37 = PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g];
                    const T x38 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x39 = PBx*x28 + tpi*x38 - PGx[g]*x37;
                    const T x40 = x15 - x28;
                    const T x41 = PBx*x30 + tpi_2*x40 - PGx[g]*x39;
                    const T x42 = x16 - x30;
                    const T x43 = PBx*x32 + tpi_3*x42 - PGx[g]*x41;
                    base_integral_26[g] += PBy*x35 - PGy[g]*x43;
                    const T x44 = x18 - x32;
                    const T x45 = PBz*x32 - PGz[g]*x41;
                    base_integral_30[g] += PBz*x33 + tpi*x44 - PGz[g]*x45;
                    base_integral_29[g] += PBy*x33 - PGy[g]*x45;
                    base_integral_25[g] += PBx*x35 + tpi_4*x44 - PGx[g]*x43;
                    base_integral_27[g] += PBz*x35 - PGz[g]*x43;
                    const T x46 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                    const T x47 = PBy*x23 + tpi*x29 - PGy[g]*x46;
                    const T x48 = x21 - x23;
                    const T x49 = PBy*x24 + tpi_2*x48 - PGy[g]*x47;
                    base_integral_21[g] += PBz*x26 - PGz[g]*x49;
                    const T x50 = PBx*x26 - PGx[g]*x49;
                    base_integral_16[g] += x50;
                    const T x51 = PBx*x24 - PGx[g]*x47;
                    base_integral_17[g] += PBz*x27 - PGz[g]*x51;
                    const T x52 = x22 - x24;
                    const T x53 = PBx*x27 + tpi*x52 - PGx[g]*x51;
                    base_integral_13[g] += x53;
                    const T x54 = PBy*x26 + tpi_3*x52 - PGy[g]*x49;
                    base_integral_20[g] += x54;
                    const T x55 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                    const T x56 = PBy*x46 + tpi*x38 - PGy[g]*x55;
                    const T x57 = PBx*x47 - PGx[g]*x56;
                    const T x58 = x24 - x47;
                    const T x59 = PBx*x51 + tpi*x58 - PGx[g]*x57;
                    const T x60 = x27 - x51;
                    base_integral_28[g] += PBx*x53 + tpi_2*x60 - PGx[g]*x59;
                    const T x61 = x23 - x46;
                    const T x62 = PBy*x47 + tpi_2*x61 - PGy[g]*x56;
                    const T x63 = PBy*x49 + tpi_3*x58 - PGy[g]*x62;
                    base_integral_35[g] += PBx*x54 - PGx[g]*x63;
                    base_integral_32[g] += PBz*x53 - PGz[g]*x59;
                    const T x64 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                    const T x65 = PBz*x3 + tpi*x29 - PGz[g]*x64;
                    const T x66 = x1 - x3;
                    const T x67 = PBz*x8 + tpi_2*x66 - PGz[g]*x65;
                    const T x68 = PBx*x10 - PGx[g]*x67;
                    base_integral_19[g] += x68;
                    const T x69 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                    const T x70 = PBz*x64 + tpi*x38 - PGz[g]*x69;
                    const T x71 = x3 - x64;
                    const T x72 = PBz*x65 + tpi_2*x71 - PGz[g]*x70;
                    const T x73 = PBx*x67 - PGx[g]*x72;
                    const T x74 = x10 - x67;
                    base_integral_34[g] += PBx*x68 + tpi*x74 - PGx[g]*x73;
                    const T x75 = PBx*x49 - PGx[g]*x62;
                    const T x76 = x26 - x49;
                    base_integral_31[g] += PBx*x50 + tpi*x76 - PGx[g]*x75;
                    base_integral_36[g] += PBz*x50 - PGz[g]*x75;
                    base_integral_40[g] += PBy*x54 + tpi_4*x76 - PGy[g]*x63;
                    base_integral_41[g] += PBz*x54 - PGz[g]*x63;
                    const T x77 = PBy*x10 - PGy[g]*x67;
                    base_integral_23[g] += x77;
                    const T x78 = x6 - x8;
                    const T x79 = PBz*x10 + tpi_3*x78 - PGz[g]*x67;
                    base_integral_24[g] += x79;
                    const T x80 = PBy*x8 - PGy[g]*x65;
                    const T x81 = PBy*x11 + tpi*x78 - PGy[g]*x80;
                    base_integral_22[g] += x81;
                    const T x82 = PBy*x65 - PGy[g]*x70;
                    const T x83 = x8 - x65;
                    const T x84 = PBy*x80 + tpi*x83 - PGy[g]*x82;
                    base_integral_37[g] += PBx*x81 - PGx[g]*x84;
                    const T x85 = PBx*x11 - PGx[g]*x80;
                    base_integral_18[g] += x85;
                    const T x86 = PBx*x80 - PGx[g]*x82;
                    const T x87 = x11 - x80;
                    base_integral_33[g] += PBx*x85 + tpi*x87 - PGx[g]*x86;
                    base_integral_42[g] += PBy*x81 + tpi_2*x87 - PGy[g]*x84;
                    const T x88 = PBz*x67 + tpi_3*x83 - PGz[g]*x72;
                    base_integral_39[g] += PBx*x79 - PGx[g]*x88;
                    const T x89 = PBy*x67 - PGy[g]*x72;
                    base_integral_38[g] += PBx*x77 - PGx[g]*x89;
                    base_integral_43[g] += PBy*x77 + tpi*x74 - PGy[g]*x89;
                    base_integral_45[g] += PBz*x79 + tpi_4*x74 - PGz[g]*x88;
                    base_integral_44[g] += PBy*x79 - PGy[g]*x88;

                } // for g
            } // if/else prim
        } // for prim

        #pragma omp simd
        for (int g = 0; g < g_block; ++g) {

            const T x0 = base_integral_25[g] - ABx*base_integral_10[g];
            const T x1 = base_integral_10[g] - ABx*base_integral_0[g];
            block_0[g] = x0 - ABx*x1;
            const T x2 = base_integral_26[g] - ABx*base_integral_11[g];
            block_10[g] = x2 - ABy*x1;
            const T x3 = base_integral_11[g] - ABx*base_integral_1[g];
            block_1[g] = x2 - ABx*x3;
            const T x4 = base_integral_27[g] - ABx*base_integral_12[g];
            block_20[g] = x4 - ABz*x1;
            const T x5 = base_integral_12[g] - ABx*base_integral_2[g];
            block_2[g] = x4 - ABx*x5;
            const T x6 = base_integral_29[g] - ABx*base_integral_14[g];
            block_12[g] = x6 - ABy*x5;
            const T x7 = base_integral_14[g] - ABx*base_integral_4[g];
            block_4[g] = x6 - ABx*x7;
            const T x8 = base_integral_29[g] - ABz*base_integral_11[g];
            const T x9 = base_integral_14[g] - ABz*base_integral_1[g];
            block_21[g] = x8 - ABx*x9;
            const T x10 = base_integral_12[g] - ABz*base_integral_0[g];
            block_40[g] = x8 - ABy*x10;
            const T x11 = base_integral_30[g] - ABz*base_integral_12[g];
            block_50[g] = x11 - ABz*x10;
            const T x12 = base_integral_28[g] - ABy*base_integral_11[g];
            const T x13 = base_integral_11[g] - ABy*base_integral_0[g];
            block_30[g] = x12 - ABy*x13;
            const T x14 = base_integral_13[g] - ABy*base_integral_1[g];
            block_11[g] = x12 - ABx*x14;
            const T x15 = base_integral_15[g] - ABz*base_integral_2[g];
            block_22[g] = x11 - ABx*x15;
            const T x16 = base_integral_33[g] - ABz*base_integral_14[g];
            block_51[g] = x16 - ABz*x9;
            block_42[g] = x16 - ABy*x15;
            const T x17 = base_integral_32[g] - ABy*base_integral_14[g];
            const T x18 = base_integral_14[g] - ABy*base_integral_2[g];
            block_32[g] = x17 - ABy*x18;
            block_41[g] = x17 - ABz*x14;
            const T x19 = base_integral_17[g] - ABy*base_integral_4[g];
            block_14[g] = x17 - ABx*x19;
            const T x20 = base_integral_18[g] - ABz*base_integral_4[g];
            block_24[g] = x16 - ABx*x20;
            const T x21 = base_integral_31[g] - ABy*base_integral_13[g];
            block_31[g] = x21 - ABy*x14;
            const T x22 = base_integral_28[g] - ABx*base_integral_13[g];
            const T x23 = base_integral_13[g] - ABx*base_integral_3[g];
            block_3[g] = x22 - ABx*x23;
            const T x24 = base_integral_32[g] - ABx*base_integral_17[g];
            block_23[g] = x24 - ABz*x23;
            const T x25 = base_integral_17[g] - ABx*base_integral_7[g];
            block_7[g] = x24 - ABx*x25;
            const T x26 = base_integral_16[g] - ABy*base_integral_3[g];
            block_13[g] = x21 - ABx*x26;
            const T x27 = base_integral_35[g] - ABy*base_integral_16[g];
            block_33[g] = x27 - ABy*x26;
            const T x28 = base_integral_31[g] - ABx*base_integral_16[g];
            const T x29 = base_integral_16[g] - ABx*base_integral_6[g];
            block_6[g] = x28 - ABx*x29;
            const T x30 = base_integral_20[g] - ABy*base_integral_6[g];
            block_16[g] = x27 - ABx*x30;
            const T x31 = base_integral_40[g] - ABy*base_integral_20[g];
            block_36[g] = x31 - ABy*x30;
            const T x32 = base_integral_36[g] - ABy*base_integral_17[g];
            block_34[g] = x32 - ABy*x19;
            const T x33 = base_integral_36[g] - ABz*base_integral_16[g];
            const T x34 = base_integral_17[g] - ABz*base_integral_3[g];
            block_43[g] = x33 - ABy*x34;
            const T x35 = base_integral_21[g] - ABz*base_integral_6[g];
            block_26[g] = x33 - ABx*x35;
            const T x36 = base_integral_21[g] - ABy*base_integral_7[g];
            block_17[g] = x32 - ABx*x36;
            const T x37 = base_integral_41[g] - ABy*base_integral_21[g];
            block_37[g] = x37 - ABy*x36;
            block_46[g] = x37 - ABz*x30;
            const T x38 = base_integral_37[g] - ABz*base_integral_17[g];
            block_53[g] = x38 - ABz*x34;
            block_44[g] = x38 - ABy*x20;
            const T x39 = base_integral_22[g] - ABz*base_integral_7[g];
            block_27[g] = x38 - ABx*x39;
            const T x40 = base_integral_42[g] - ABz*base_integral_21[g];
            block_56[g] = x40 - ABz*x35;
            block_47[g] = x40 - ABy*x39;
            const T x41 = base_integral_34[g] - ABz*base_integral_15[g];
            block_52[g] = x41 - ABz*x15;
            const T x42 = base_integral_30[g] - ABx*base_integral_15[g];
            const T x43 = base_integral_15[g] - ABx*base_integral_5[g];
            block_5[g] = x42 - ABx*x43;
            const T x44 = base_integral_33[g] - ABx*base_integral_18[g];
            block_15[g] = x44 - ABy*x43;
            const T x45 = base_integral_18[g] - ABx*base_integral_8[g];
            block_8[g] = x44 - ABx*x45;
            const T x46 = base_integral_37[g] - ABy*base_integral_18[g];
            const T x47 = base_integral_18[g] - ABy*base_integral_5[g];
            block_35[g] = x46 - ABy*x47;
            const T x48 = base_integral_22[g] - ABy*base_integral_8[g];
            block_18[g] = x46 - ABx*x48;
            const T x49 = base_integral_42[g] - ABy*base_integral_22[g];
            block_38[g] = x49 - ABy*x48;
            const T x50 = base_integral_19[g] - ABz*base_integral_5[g];
            block_25[g] = x41 - ABx*x50;
            const T x51 = base_integral_38[g] - ABz*base_integral_18[g];
            block_54[g] = x51 - ABz*x20;
            block_45[g] = x51 - ABy*x50;
            const T x52 = base_integral_23[g] - ABz*base_integral_8[g];
            block_28[g] = x51 - ABx*x52;
            const T x53 = base_integral_39[g] - ABz*base_integral_19[g];
            block_55[g] = x53 - ABz*x50;
            const T x54 = base_integral_43[g] - ABz*base_integral_22[g];
            block_57[g] = x54 - ABz*x39;
            block_48[g] = x54 - ABy*x52;
            const T x55 = base_integral_44[g] - ABz*base_integral_23[g];
            block_58[g] = x55 - ABz*x52;
            const T x56 = base_integral_34[g] - ABx*base_integral_19[g];
            const T x57 = base_integral_19[g] - ABx*base_integral_9[g];
            block_9[g] = x56 - ABx*x57;
            const T x58 = base_integral_38[g] - ABy*base_integral_19[g];
            const T x59 = base_integral_23[g] - ABy*base_integral_9[g];
            block_19[g] = x58 - ABx*x59;
            const T x60 = base_integral_43[g] - ABy*base_integral_23[g];
            block_39[g] = x60 - ABy*x59;
            const T x61 = base_integral_24[g] - ABz*base_integral_9[g];
            block_29[g] = x53 - ABx*x61;
            block_49[g] = x55 - ABy*x61;
            const T x62 = base_integral_45[g] - ABz*base_integral_24[g];
            block_59[g] = x62 - ABz*x61;

        } // for g

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
} // HRR_D_F

template void HRR_D_F<float>(
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

template void HRR_D_F<double>(
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

