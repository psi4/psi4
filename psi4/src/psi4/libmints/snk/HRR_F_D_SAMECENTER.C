#include "HRR_F_D_SAMECENTER.h"

template <typename T>
void HRR_F_D_SAMECENTER(
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
    T* __restrict base_integral_15 = &base_integrals[15*gb_size];
    T* __restrict base_integral_16 = &base_integrals[16*gb_size];
    T* __restrict base_integral_17 = &base_integrals[17*gb_size];
    T* __restrict base_integral_18 = &base_integrals[18*gb_size];
    T* __restrict base_integral_19 = &base_integrals[19*gb_size];
    T* __restrict base_integral_20 = &base_integrals[20*gb_size];

    for (int prim = 0; prim < nprim; ++prim) {

        const T tpi = 0.5 / v_p[prim];
        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*6*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*6*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*6*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*6*gb_size+3*gb_size];
        const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*6*gb_size+4*gb_size];
        const T* __restrict scaled_boys_5 = &v_scaled_boys[prim*6*gb_size+5*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_5[g]*x0;
                const T x2 = scaled_boys_3[g] - scaled_boys_4[g];
                const T x3 = tpi*x2;
                const T x4 = -tpi*x2;
                const T x5 = -2*x4;
                const T x6 = scaled_boys_3[g]*x0 - scaled_boys_4[g]*x0;
                const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x8 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x9 = tpi*x7 - tpi*x8;
                const T x10 = x6 + x9;
                const T x11 = 3*tpi;
                const T x12 = -tpi*x7;
                const T x13 = -tpi*x8;
                const T x14 = -2*x12 + 2*x13;
                const T x15 = -x4;
                const T x16 = 2*x3;
                const T x17 = -x12 + x13;
                const T x18 = x0*(-x1 - x15 - x16) - x11*(x17 + x6);
                const T x19 = PGy[g]*PGy[g];
                const T x20 = scaled_boys_5[g]*x19;
                const T x21 = x20 + x3;
                const T x22 = scaled_boys_3[g]*x19 - scaled_boys_4[g]*x19;
                const T x23 = x22 + x9;
                const T x24 = x11*x23;
                const T x25 = 3*x3;
                const T x26 = PGx[g]*PGy[g]*PGz[g];
                const T x27 = PGz[g]*PGz[g];
                const T x28 = scaled_boys_5[g]*x27;
                const T x29 = x28 + x3;
                const T x30 = x0*x29;
                const T x31 = scaled_boys_3[g]*x27 - scaled_boys_4[g]*x27 + x9;
                const T x32 = tpi*x31;
                const T x33 = 3*x32;
                const T x34 = x21 + x5;
                const T x35 = tpi*(x14 + x23);
                const T x36 = x15 + x20;
                const T x37 = tpi*(x17 + x22);
                const T x38 = x29 + x5;
                const T x39 = tpi*(x14 + x31);
                const T x40 = x16 + x36;
                const T x41 = x19*x29;
                const T x42 = x27*x38 + x33;

                base_integral_0[g] = PGx[g]*(-4*tpi*(x10 + x14) + x0*(-x1 - x3 - x5) - x10*x11);
                base_integral_1[g] = PGy[g]*x18;
                base_integral_2[g] = PGz[g]*x18;
                base_integral_3[g] = -PGx[g]*(x0*x21 + x24);
                base_integral_4[g] = -x26*(x1 + x25);
                base_integral_5[g] = -PGx[g]*(x30 + x33);
                base_integral_6[g] = -PGy[g]*(x0*x34 + x35);
                base_integral_7[g] = -PGz[g]*(x0*x36 + x37);
                base_integral_8[g] = -PGy[g]*(x30 + x32);
                base_integral_9[g] = -PGz[g]*(x0*x38 + x39);
                base_integral_10[g] = -PGx[g]*(x19*x34 + x24);
                base_integral_11[g] = -x26*x40;
                base_integral_12[g] = -PGx[g]*(x32 + x41);
                base_integral_13[g] = -x26*(x25 + x28);
                base_integral_14[g] = -PGx[g]*x42;
                base_integral_15[g] = PGy[g]*(-x19*x34 - x24 - 4*x35);
                base_integral_16[g] = PGz[g]*(-x19*x40 - 3*x37);
                base_integral_17[g] = -PGy[g]*(x33 + x41);
                base_integral_18[g] = -PGz[g]*(x19*x38 + x39);
                base_integral_19[g] = -PGy[g]*x42;
                base_integral_20[g] = PGz[g]*(-x27*x38 - x33 - 4*x39);

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_5[g]*x0;
                const T x2 = scaled_boys_3[g] - scaled_boys_4[g];
                const T x3 = tpi*x2;
                const T x4 = -tpi*x2;
                const T x5 = -2*x4;
                const T x6 = scaled_boys_3[g]*x0 - scaled_boys_4[g]*x0;
                const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x8 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x9 = tpi*x7 - tpi*x8;
                const T x10 = x6 + x9;
                const T x11 = 3*tpi;
                const T x12 = -tpi*x7;
                const T x13 = -tpi*x8;
                const T x14 = -2*x12 + 2*x13;
                const T x15 = -x4;
                const T x16 = 2*x3;
                const T x17 = -x12 + x13;
                const T x18 = x0*(-x1 - x15 - x16) - x11*(x17 + x6);
                const T x19 = PGy[g]*PGy[g];
                const T x20 = scaled_boys_5[g]*x19;
                const T x21 = x20 + x3;
                const T x22 = scaled_boys_3[g]*x19 - scaled_boys_4[g]*x19;
                const T x23 = x22 + x9;
                const T x24 = x11*x23;
                const T x25 = 3*x3;
                const T x26 = PGx[g]*PGy[g]*PGz[g];
                const T x27 = PGz[g]*PGz[g];
                const T x28 = scaled_boys_5[g]*x27;
                const T x29 = x28 + x3;
                const T x30 = x0*x29;
                const T x31 = scaled_boys_3[g]*x27 - scaled_boys_4[g]*x27 + x9;
                const T x32 = tpi*x31;
                const T x33 = 3*x32;
                const T x34 = x21 + x5;
                const T x35 = tpi*(x14 + x23);
                const T x36 = x15 + x20;
                const T x37 = tpi*(x17 + x22);
                const T x38 = x29 + x5;
                const T x39 = tpi*(x14 + x31);
                const T x40 = x16 + x36;
                const T x41 = x19*x29;
                const T x42 = x27*x38 + x33;

                base_integral_0[g] += PGx[g]*(-4*tpi*(x10 + x14) + x0*(-x1 - x3 - x5) - x10*x11);
                base_integral_1[g] += PGy[g]*x18;
                base_integral_2[g] += PGz[g]*x18;
                base_integral_3[g] += -PGx[g]*(x0*x21 + x24);
                base_integral_4[g] += -x26*(x1 + x25);
                base_integral_5[g] += -PGx[g]*(x30 + x33);
                base_integral_6[g] += -PGy[g]*(x0*x34 + x35);
                base_integral_7[g] += -PGz[g]*(x0*x36 + x37);
                base_integral_8[g] += -PGy[g]*(x30 + x32);
                base_integral_9[g] += -PGz[g]*(x0*x38 + x39);
                base_integral_10[g] += -PGx[g]*(x19*x34 + x24);
                base_integral_11[g] += -x26*x40;
                base_integral_12[g] += -PGx[g]*(x32 + x41);
                base_integral_13[g] += -x26*(x25 + x28);
                base_integral_14[g] += -PGx[g]*x42;
                base_integral_15[g] += PGy[g]*(-x19*x34 - x24 - 4*x35);
                base_integral_16[g] += PGz[g]*(-x19*x40 - 3*x37);
                base_integral_17[g] += -PGy[g]*(x33 + x41);
                base_integral_18[g] += -PGz[g]*(x19*x38 + x39);
                base_integral_19[g] += -PGy[g]*x42;
                base_integral_20[g] += PGz[g]*(-x27*x38 - x33 - 4*x39);

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

    #pragma omp simd
    for (int g = 0; g < gb_size; ++g) {

        integral_0[g] = base_integral_0[g];
        integral_1[g] = base_integral_1[g];
        integral_2[g] = base_integral_2[g];
        integral_3[g] = base_integral_3[g];
        integral_4[g] = base_integral_4[g];
        integral_5[g] = base_integral_5[g];
        integral_6[g] = base_integral_1[g];
        integral_7[g] = base_integral_3[g];
        integral_8[g] = base_integral_4[g];
        integral_9[g] = base_integral_6[g];
        integral_10[g] = base_integral_7[g];
        integral_11[g] = base_integral_8[g];
        integral_12[g] = base_integral_2[g];
        integral_13[g] = base_integral_4[g];
        integral_14[g] = base_integral_5[g];
        integral_15[g] = base_integral_7[g];
        integral_16[g] = base_integral_8[g];
        integral_17[g] = base_integral_9[g];
        integral_18[g] = base_integral_3[g];
        integral_19[g] = base_integral_6[g];
        integral_20[g] = base_integral_7[g];
        integral_21[g] = base_integral_10[g];
        integral_22[g] = base_integral_11[g];
        integral_23[g] = base_integral_12[g];
        integral_24[g] = base_integral_4[g];
        integral_25[g] = base_integral_7[g];
        integral_26[g] = base_integral_8[g];
        integral_27[g] = base_integral_11[g];
        integral_28[g] = base_integral_12[g];
        integral_29[g] = base_integral_13[g];
        integral_30[g] = base_integral_5[g];
        integral_31[g] = base_integral_8[g];
        integral_32[g] = base_integral_9[g];
        integral_33[g] = base_integral_12[g];
        integral_34[g] = base_integral_13[g];
        integral_35[g] = base_integral_14[g];
        integral_36[g] = base_integral_6[g];
        integral_37[g] = base_integral_10[g];
        integral_38[g] = base_integral_11[g];
        integral_39[g] = base_integral_15[g];
        integral_40[g] = base_integral_16[g];
        integral_41[g] = base_integral_17[g];
        integral_42[g] = base_integral_7[g];
        integral_43[g] = base_integral_11[g];
        integral_44[g] = base_integral_12[g];
        integral_45[g] = base_integral_16[g];
        integral_46[g] = base_integral_17[g];
        integral_47[g] = base_integral_18[g];
        integral_48[g] = base_integral_8[g];
        integral_49[g] = base_integral_12[g];
        integral_50[g] = base_integral_13[g];
        integral_51[g] = base_integral_17[g];
        integral_52[g] = base_integral_18[g];
        integral_53[g] = base_integral_19[g];
        integral_54[g] = base_integral_9[g];
        integral_55[g] = base_integral_13[g];
        integral_56[g] = base_integral_14[g];
        integral_57[g] = base_integral_18[g];
        integral_58[g] = base_integral_19[g];
        integral_59[g] = base_integral_20[g];

    } // for g
} // HRR_F_D_SAMECENTER

template void HRR_F_D_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_F_D_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

