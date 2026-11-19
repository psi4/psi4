#include "HRR_P_F_SAMECENTER.h"

template <typename T>
void HRR_P_F_SAMECENTER(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        T* base_integrals,         //!< overwritten, base [LA+LB|0] integrals for HRR
        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)
        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively
        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)
        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair
        int gb_size)               //!< number of gridpoints in this grid batch 
{
    T* __restrict base_integral_0 = &base_integrals[0*gb_size];
    T* __restrict base_integral_1 = &base_integrals[1*gb_size];
    T* __restrict base_integral_2 = &base_integrals[2*gb_size];
    T* __restrict base_integral_3 = &base_integrals[3*gb_size];
    T* __restrict base_integral_4 = &base_integrals[4*gb_size];
    T* __restrict base_integral_5 = &base_integrals[5*gb_size];
    T* __restrict base_integral_6 = &base_integrals[6*gb_size];
    T* __restrict base_integral_7 = &base_integrals[7*gb_size];
    T* __restrict base_integral_8 = &base_integrals[8*gb_size];
    T* __restrict base_integral_9 = &base_integrals[9*gb_size];
    T* __restrict base_integral_10 = &base_integrals[10*gb_size];
    T* __restrict base_integral_11 = &base_integrals[11*gb_size];
    T* __restrict base_integral_12 = &base_integrals[12*gb_size];
    T* __restrict base_integral_13 = &base_integrals[13*gb_size];
    T* __restrict base_integral_14 = &base_integrals[14*gb_size];

    for (int prim = 0; prim < nprim; ++prim) {

        const T tpi = 0.5 / v_p[prim];
        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*5*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*5*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*5*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*5*gb_size+3*gb_size];
        const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*5*gb_size+4*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x2 = -tpi*x1;
                const T x3 = -2*x2;
                const T x4 = scaled_boys_4[g]*x0;
                const T x5 = tpi*x1;
                const T x6 = x4 + x5;
                const T x7 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]) - tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x8 = -x2 + 2*x5;
                const T x9 = PGx[g]*(x4 + x8);
                const T x10 = PGy[g]*PGy[g];
                const T x11 = scaled_boys_4[g]*x10;
                const T x12 = x11 + x5;
                const T x13 = tpi*(scaled_boys_2[g]*x10 - scaled_boys_3[g]*x10 + x7);
                const T x14 = PGy[g]*PGz[g];
                const T x15 = PGz[g]*PGz[g];
                const T x16 = scaled_boys_4[g]*x15;
                const T x17 = x16 + x5;
                const T x18 = tpi*(scaled_boys_2[g]*x15 - scaled_boys_3[g]*x15 + x7);
                const T x19 = 3*x5;
                const T x20 = PGx[g]*PGy[g];
                const T x21 = PGx[g]*PGz[g];
                const T x22 = x16 + x19;

                base_integral_0[g] = 3*tpi*(scaled_boys_2[g]*x0 - scaled_boys_3[g]*x0 + x7) + x0*(x3 + x6);
                base_integral_1[g] = PGy[g]*x9;
                base_integral_2[g] = PGz[g]*x9;
                base_integral_3[g] = x0*x12 + x13;
                base_integral_4[g] = x14*x6;
                base_integral_5[g] = x0*x17 + x18;
                base_integral_6[g] = x20*(x11 + x19);
                base_integral_7[g] = x12*x21;
                base_integral_8[g] = x17*x20;
                base_integral_9[g] = x21*x22;
                base_integral_10[g] = x10*(x12 + x3) + 3*x13;
                base_integral_11[g] = x14*(x11 + x8);
                base_integral_12[g] = x10*x17 + x18;
                base_integral_13[g] = x14*x22;
                base_integral_14[g] = x15*(x17 + x3) + 3*x18;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x2 = -tpi*x1;
                const T x3 = -2*x2;
                const T x4 = scaled_boys_4[g]*x0;
                const T x5 = tpi*x1;
                const T x6 = x4 + x5;
                const T x7 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]) - tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x8 = -x2 + 2*x5;
                const T x9 = PGx[g]*(x4 + x8);
                const T x10 = PGy[g]*PGy[g];
                const T x11 = scaled_boys_4[g]*x10;
                const T x12 = x11 + x5;
                const T x13 = tpi*(scaled_boys_2[g]*x10 - scaled_boys_3[g]*x10 + x7);
                const T x14 = PGy[g]*PGz[g];
                const T x15 = PGz[g]*PGz[g];
                const T x16 = scaled_boys_4[g]*x15;
                const T x17 = x16 + x5;
                const T x18 = tpi*(scaled_boys_2[g]*x15 - scaled_boys_3[g]*x15 + x7);
                const T x19 = 3*x5;
                const T x20 = PGx[g]*PGy[g];
                const T x21 = PGx[g]*PGz[g];
                const T x22 = x16 + x19;

                base_integral_0[g] += 3*tpi*(scaled_boys_2[g]*x0 - scaled_boys_3[g]*x0 + x7) + x0*(x3 + x6);
                base_integral_1[g] += PGy[g]*x9;
                base_integral_2[g] += PGz[g]*x9;
                base_integral_3[g] += x0*x12 + x13;
                base_integral_4[g] += x14*x6;
                base_integral_5[g] += x0*x17 + x18;
                base_integral_6[g] += x20*(x11 + x19);
                base_integral_7[g] += x12*x21;
                base_integral_8[g] += x17*x20;
                base_integral_9[g] += x21*x22;
                base_integral_10[g] += x10*(x12 + x3) + 3*x13;
                base_integral_11[g] += x14*(x11 + x8);
                base_integral_12[g] += x10*x17 + x18;
                base_integral_13[g] += x14*x22;
                base_integral_14[g] += x15*(x17 + x3) + 3*x18;

            } // for g
        } // if/else prim
    } // for prim

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

    #pragma omp simd
    for (int g = 0; g < gb_size; ++g) {

        integral_0[g] = base_integral_0[g];
        integral_1[g] = base_integral_1[g];
        integral_2[g] = base_integral_2[g];
        integral_3[g] = base_integral_3[g];
        integral_4[g] = base_integral_4[g];
        integral_5[g] = base_integral_5[g];
        integral_6[g] = base_integral_6[g];
        integral_7[g] = base_integral_7[g];
        integral_8[g] = base_integral_8[g];
        integral_9[g] = base_integral_9[g];
        integral_10[g] = base_integral_1[g];
        integral_11[g] = base_integral_3[g];
        integral_12[g] = base_integral_4[g];
        integral_13[g] = base_integral_6[g];
        integral_14[g] = base_integral_7[g];
        integral_15[g] = base_integral_8[g];
        integral_16[g] = base_integral_10[g];
        integral_17[g] = base_integral_11[g];
        integral_18[g] = base_integral_12[g];
        integral_19[g] = base_integral_13[g];
        integral_20[g] = base_integral_2[g];
        integral_21[g] = base_integral_4[g];
        integral_22[g] = base_integral_5[g];
        integral_23[g] = base_integral_7[g];
        integral_24[g] = base_integral_8[g];
        integral_25[g] = base_integral_9[g];
        integral_26[g] = base_integral_11[g];
        integral_27[g] = base_integral_12[g];
        integral_28[g] = base_integral_13[g];
        integral_29[g] = base_integral_14[g];

    } // for g
} // HRR_P_F_SAMECENTER

template void HRR_P_F_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_P_F_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

