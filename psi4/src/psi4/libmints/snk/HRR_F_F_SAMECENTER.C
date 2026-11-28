#include "HRR_F_F_SAMECENTER.h"

#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

template <typename T>
void HRR_F_F_SAMECENTER(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        T* base_integrals,         //!< work array for the base [LA+LB|0] integrals of the HRR (not used by all kernels)
        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)
        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively
        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)
        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair
        int gb_size)               //!< number of gridpoints in this grid batch 
{
    if (nprim < 1) return;

    // The gridpoints are processed in blocks.  The integrals of a block are contracted in
    // a buffer that stays in the L1 cache and are then copied to the output arrays.
    enum { g_block = (sizeof(T) <= 4) ? 32 : 32 };
    alignas(64) T base_integral[28*g_block];
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

            const T tpi = 0.5 / v_p[prim];
            const T tpi_p2 = tpi*tpi;
            const T tpi_p3 = tpi*tpi_p2;

            // overwrite integrals for first primitive, accumulate otherwise
            if (prim==0) {
                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x1 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x2 = x0 - x1;
                    const T x3 = tpi_p2*x2;
                    const T x4 = scaled_boys_4[g] - scaled_boys_5[g];
                    const T x5 = tpi*x4;
                    const T x6 = PGz[g]*x5;
                    const T x7 = x3 + PGz[g]*x6;
                    const T x8 = PGz[g]*x3;
                    const T x9 = PGz[g]*x7 + 2*x8;
                    const T x10 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x11 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x12 = x10 - x11;
                    const T x13 = x11 - x0;
                    const T x14 = x12 - x13;
                    const T x15 = tpi_p3*x14;
                    const T x16 = x15 + PGz[g]*x8;
                    const T x17 = PGz[g]*x9 + 3*x16;
                    const T x18 = PGz[g]*scaled_boys_6[g];
                    const T x19 = x5 + PGz[g]*x18;
                    const T x20 = PGz[g]*x19 + 2*x6;
                    const T x21 = PGz[g]*x20 + 3*x7;
                    const T x22 = PGx[g]*x21;
                    base_integral_14[g] = PGx[g]*x22 + x17;
                    const T x23 = 4*x9 + PGz[g]*x21;
                    base_integral_27[g] = 5*x17 + PGz[g]*x23;
                    const T x24 = PGy[g]*x21;
                    base_integral_25[g] = x17 + PGy[g]*x24;
                    base_integral_19[g] = PGx[g]*x24;
                    base_integral_20[g] = PGx[g]*x23;
                    base_integral_26[g] = PGy[g]*x23;
                    const T x25 = PGy[g]*x20;
                    const T x26 = PGy[g]*x25 + x9;
                    base_integral_18[g] = PGx[g]*x26;
                    const T x27 = PGy[g]*x9;
                    base_integral_24[g] = 2*x27 + PGy[g]*x26;
                    const T x28 = PGx[g]*x25;
                    base_integral_13[g] = PGx[g]*x28 + x27;
                    const T x29 = PGx[g]*x5;
                    const T x30 = PGx[g]*x29 + x3;
                    const T x31 = PGx[g]*scaled_boys_6[g];
                    const T x32 = PGx[g]*x31 + x5;
                    const T x33 = PGx[g]*x32 + 2*x29;
                    const T x34 = PGx[g]*x33 + 3*x30;
                    const T x35 = PGz[g]*x34;
                    base_integral_4[g] = PGy[g]*x35;
                    const T x36 = PGx[g]*x3;
                    const T x37 = PGx[g]*x30 + 2*x36;
                    const T x38 = PGx[g]*x36 + x15;
                    const T x39 = PGx[g]*x37 + 3*x38;
                    base_integral_5[g] = x39 + PGz[g]*x35;
                    const T x40 = PGx[g]*x34 + 4*x37;
                    base_integral_0[g] = PGx[g]*x40 + 5*x39;
                    const T x41 = PGz[g]*x33;
                    const T x42 = x37 + PGz[g]*x41;
                    const T x43 = 2*x37 + x42;
                    base_integral_8[g] = PGy[g]*x42;
                    base_integral_9[g] = PGz[g]*x43;
                    base_integral_1[g] = PGy[g]*x40;
                    base_integral_2[g] = PGz[g]*x40;
                    const T x44 = PGy[g]*x5;
                    const T x45 = x3 + PGy[g]*x44;
                    const T x46 = PGy[g]*x3;
                    const T x47 = PGy[g]*x45 + 2*x46;
                    const T x48 = x15 + PGy[g]*x46;
                    const T x49 = PGy[g]*scaled_boys_6[g];
                    const T x50 = x5 + PGy[g]*x49;
                    const T x51 = PGy[g]*x50 + 2*x44;
                    const T x52 = PGz[g]*x51;
                    const T x53 = x47 + PGz[g]*x52;
                    base_integral_17[g] = PGx[g]*x53;
                    const T x54 = PGx[g]*x51;
                    const T x55 = PGx[g]*x54 + x47;
                    base_integral_11[g] = PGz[g]*x55;
                    const T x56 = 2*x47 + x55;
                    base_integral_6[g] = PGx[g]*x56;
                    const T x57 = 3*x48 + PGy[g]*x47;
                    const T x58 = 3*x45 + PGy[g]*x51;
                    const T x59 = 4*x47 + PGy[g]*x58;
                    base_integral_21[g] = PGy[g]*x59 + 5*x57;
                    const T x60 = PGz[g]*x58;
                    base_integral_23[g] = PGz[g]*x60 + x57;
                    base_integral_16[g] = PGx[g]*x60;
                    const T x61 = PGx[g]*x45;
                    const T x62 = PGx[g]*x61 + x48;
                    base_integral_10[g] = 3*x62 + PGy[g]*x55;
                    const T x63 = PGx[g]*x50;
                    const T x64 = PGx[g]*x63 + x45;
                    const T x65 = PGx[g]*x64 + 2*x61;
                    base_integral_3[g] = PGx[g]*x65 + 3*x62;
                    base_integral_7[g] = PGz[g]*x65;
                    const T x66 = PGz[g]*x64;
                    base_integral_12[g] = x62 + PGz[g]*x66;
                    base_integral_15[g] = PGx[g]*x59;
                    base_integral_22[g] = PGz[g]*x59;

                } // for g
            } else {

                #pragma omp simd
                for (int g = 0; g < g_block; ++g) {

                    const T x0 = scaled_boys_2[g] - scaled_boys_3[g];
                    const T x1 = scaled_boys_3[g] - scaled_boys_4[g];
                    const T x2 = x0 - x1;
                    const T x3 = tpi_p2*x2;
                    const T x4 = scaled_boys_4[g] - scaled_boys_5[g];
                    const T x5 = tpi*x4;
                    const T x6 = PGz[g]*x5;
                    const T x7 = x3 + PGz[g]*x6;
                    const T x8 = PGz[g]*x3;
                    const T x9 = PGz[g]*x7 + 2*x8;
                    const T x10 = scaled_boys_0[g] - scaled_boys_1[g];
                    const T x11 = scaled_boys_1[g] - scaled_boys_2[g];
                    const T x12 = x10 - x11;
                    const T x13 = x11 - x0;
                    const T x14 = x12 - x13;
                    const T x15 = tpi_p3*x14;
                    const T x16 = x15 + PGz[g]*x8;
                    const T x17 = PGz[g]*x9 + 3*x16;
                    const T x18 = PGz[g]*scaled_boys_6[g];
                    const T x19 = x5 + PGz[g]*x18;
                    const T x20 = PGz[g]*x19 + 2*x6;
                    const T x21 = PGz[g]*x20 + 3*x7;
                    const T x22 = PGx[g]*x21;
                    base_integral_14[g] += PGx[g]*x22 + x17;
                    const T x23 = 4*x9 + PGz[g]*x21;
                    base_integral_27[g] += 5*x17 + PGz[g]*x23;
                    const T x24 = PGy[g]*x21;
                    base_integral_25[g] += x17 + PGy[g]*x24;
                    base_integral_19[g] += PGx[g]*x24;
                    base_integral_20[g] += PGx[g]*x23;
                    base_integral_26[g] += PGy[g]*x23;
                    const T x25 = PGy[g]*x20;
                    const T x26 = PGy[g]*x25 + x9;
                    base_integral_18[g] += PGx[g]*x26;
                    const T x27 = PGy[g]*x9;
                    base_integral_24[g] += 2*x27 + PGy[g]*x26;
                    const T x28 = PGx[g]*x25;
                    base_integral_13[g] += PGx[g]*x28 + x27;
                    const T x29 = PGx[g]*x5;
                    const T x30 = PGx[g]*x29 + x3;
                    const T x31 = PGx[g]*scaled_boys_6[g];
                    const T x32 = PGx[g]*x31 + x5;
                    const T x33 = PGx[g]*x32 + 2*x29;
                    const T x34 = PGx[g]*x33 + 3*x30;
                    const T x35 = PGz[g]*x34;
                    base_integral_4[g] += PGy[g]*x35;
                    const T x36 = PGx[g]*x3;
                    const T x37 = PGx[g]*x30 + 2*x36;
                    const T x38 = PGx[g]*x36 + x15;
                    const T x39 = PGx[g]*x37 + 3*x38;
                    base_integral_5[g] += x39 + PGz[g]*x35;
                    const T x40 = PGx[g]*x34 + 4*x37;
                    base_integral_0[g] += PGx[g]*x40 + 5*x39;
                    const T x41 = PGz[g]*x33;
                    const T x42 = x37 + PGz[g]*x41;
                    const T x43 = 2*x37 + x42;
                    base_integral_8[g] += PGy[g]*x42;
                    base_integral_9[g] += PGz[g]*x43;
                    base_integral_1[g] += PGy[g]*x40;
                    base_integral_2[g] += PGz[g]*x40;
                    const T x44 = PGy[g]*x5;
                    const T x45 = x3 + PGy[g]*x44;
                    const T x46 = PGy[g]*x3;
                    const T x47 = PGy[g]*x45 + 2*x46;
                    const T x48 = x15 + PGy[g]*x46;
                    const T x49 = PGy[g]*scaled_boys_6[g];
                    const T x50 = x5 + PGy[g]*x49;
                    const T x51 = PGy[g]*x50 + 2*x44;
                    const T x52 = PGz[g]*x51;
                    const T x53 = x47 + PGz[g]*x52;
                    base_integral_17[g] += PGx[g]*x53;
                    const T x54 = PGx[g]*x51;
                    const T x55 = PGx[g]*x54 + x47;
                    base_integral_11[g] += PGz[g]*x55;
                    const T x56 = 2*x47 + x55;
                    base_integral_6[g] += PGx[g]*x56;
                    const T x57 = 3*x48 + PGy[g]*x47;
                    const T x58 = 3*x45 + PGy[g]*x51;
                    const T x59 = 4*x47 + PGy[g]*x58;
                    base_integral_21[g] += PGy[g]*x59 + 5*x57;
                    const T x60 = PGz[g]*x58;
                    base_integral_23[g] += PGz[g]*x60 + x57;
                    base_integral_16[g] += PGx[g]*x60;
                    const T x61 = PGx[g]*x45;
                    const T x62 = PGx[g]*x61 + x48;
                    base_integral_10[g] += 3*x62 + PGy[g]*x55;
                    const T x63 = PGx[g]*x50;
                    const T x64 = PGx[g]*x63 + x45;
                    const T x65 = PGx[g]*x64 + 2*x61;
                    base_integral_3[g] += PGx[g]*x65 + 3*x62;
                    base_integral_7[g] += PGz[g]*x65;
                    const T x66 = PGz[g]*x64;
                    base_integral_12[g] += x62 + PGz[g]*x66;
                    base_integral_15[g] += PGx[g]*x59;
                    base_integral_22[g] += PGz[g]*x59;

                } // for g
            } // if/else prim
        } // for prim

        // copy the block to the output arrays
        static const short row[100] = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 3, 4, 6, 7, 8, 10, 11, 12, 13,
            2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 3, 6, 7, 10, 11, 12, 15, 16, 17, 18,
            4, 7, 8, 11, 12, 13, 16, 17, 18, 19, 5, 8, 9, 12, 13, 14, 17, 18, 19, 20,
            6, 10, 11, 15, 16, 17, 21, 22, 23, 24, 7, 11, 12, 16, 17, 18, 22, 23, 24, 25,
            8, 12, 13, 17, 18, 19, 23, 24, 25, 26, 9, 13, 14, 18, 19, 20, 24, 25, 26, 27,
        };
        for (int i = 0; i < 100; ++i) {
            const T* __restrict src = &base_integral[row[i]*g_block];
            T* __restrict dst = &integrals[i*gb_size + g_lo];
            if (padded) {
                SNK_COPY(dst, src, gb_size);
            } else {
                SNK_COPY(dst, src, g_block);
            }
        }
    } // for g_lo
} // HRR_F_F_SAMECENTER

template void HRR_F_F_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_F_F_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

