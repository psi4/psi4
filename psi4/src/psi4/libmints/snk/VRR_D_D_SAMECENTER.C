#include "VRR_D_D_SAMECENTER.h"

template <typename T>
void VRR_D_D_SAMECENTER(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)
        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively
        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)
        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair
        int gb_size)               //!< number of gridpoints in this grid batch 
{
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
                const T x10 = PGy[g]*x9;
                const T x11 = PGz[g]*x9;
                const T x12 = PGy[g]*PGy[g];
                const T x13 = scaled_boys_4[g]*x12;
                const T x14 = x13 + x5;
                const T x15 = tpi*(scaled_boys_2[g]*x12 - scaled_boys_3[g]*x12 + x7);
                const T x16 = x0*x14 + x15;
                const T x17 = PGy[g]*PGz[g];
                const T x18 = x17*x6;
                const T x19 = PGz[g]*PGz[g];
                const T x20 = scaled_boys_4[g]*x19;
                const T x21 = x20 + x5;
                const T x22 = tpi*(scaled_boys_2[g]*x19 - scaled_boys_3[g]*x19 + x7);
                const T x23 = x0*x21 + x22;
                const T x24 = 3*x5;
                const T x25 = PGx[g]*PGy[g];
                const T x26 = x25*(x13 + x24);
                const T x27 = PGx[g]*PGz[g];
                const T x28 = x14*x27;
                const T x29 = x21*x25;
                const T x30 = x20 + x24;
                const T x31 = x27*x30;
                const T x32 = x17*(x13 + x8);
                const T x33 = x12*x21 + x22;
                const T x34 = x17*x30;

                integral_0[g] = 3*tpi*(scaled_boys_2[g]*x0 - scaled_boys_3[g]*x0 + x7) + x0*(x3 + x6);
                integral_1[g] = x10;
                integral_2[g] = x11;
                integral_3[g] = x16;
                integral_4[g] = x18;
                integral_5[g] = x23;
                integral_6[g] = x10;
                integral_7[g] = x16;
                integral_8[g] = x18;
                integral_9[g] = x26;
                integral_10[g] = x28;
                integral_11[g] = x29;
                integral_12[g] = x11;
                integral_13[g] = x18;
                integral_14[g] = x23;
                integral_15[g] = x28;
                integral_16[g] = x29;
                integral_17[g] = x31;
                integral_18[g] = x16;
                integral_19[g] = x26;
                integral_20[g] = x28;
                integral_21[g] = x12*(x14 + x3) + 3*x15;
                integral_22[g] = x32;
                integral_23[g] = x33;
                integral_24[g] = x18;
                integral_25[g] = x28;
                integral_26[g] = x29;
                integral_27[g] = x32;
                integral_28[g] = x33;
                integral_29[g] = x34;
                integral_30[g] = x23;
                integral_31[g] = x29;
                integral_32[g] = x31;
                integral_33[g] = x33;
                integral_34[g] = x34;
                integral_35[g] = x19*(x21 + x3) + 3*x22;

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
                const T x10 = PGy[g]*x9;
                const T x11 = PGz[g]*x9;
                const T x12 = PGy[g]*PGy[g];
                const T x13 = scaled_boys_4[g]*x12;
                const T x14 = x13 + x5;
                const T x15 = tpi*(scaled_boys_2[g]*x12 - scaled_boys_3[g]*x12 + x7);
                const T x16 = x0*x14 + x15;
                const T x17 = PGy[g]*PGz[g];
                const T x18 = x17*x6;
                const T x19 = PGz[g]*PGz[g];
                const T x20 = scaled_boys_4[g]*x19;
                const T x21 = x20 + x5;
                const T x22 = tpi*(scaled_boys_2[g]*x19 - scaled_boys_3[g]*x19 + x7);
                const T x23 = x0*x21 + x22;
                const T x24 = 3*x5;
                const T x25 = PGx[g]*PGy[g];
                const T x26 = x25*(x13 + x24);
                const T x27 = PGx[g]*PGz[g];
                const T x28 = x14*x27;
                const T x29 = x21*x25;
                const T x30 = x20 + x24;
                const T x31 = x27*x30;
                const T x32 = x17*(x13 + x8);
                const T x33 = x12*x21 + x22;
                const T x34 = x17*x30;

                integral_0[g] += 3*tpi*(scaled_boys_2[g]*x0 - scaled_boys_3[g]*x0 + x7) + x0*(x3 + x6);
                integral_1[g] += x10;
                integral_2[g] += x11;
                integral_3[g] += x16;
                integral_4[g] += x18;
                integral_5[g] += x23;
                integral_6[g] += x10;
                integral_7[g] += x16;
                integral_8[g] += x18;
                integral_9[g] += x26;
                integral_10[g] += x28;
                integral_11[g] += x29;
                integral_12[g] += x11;
                integral_13[g] += x18;
                integral_14[g] += x23;
                integral_15[g] += x28;
                integral_16[g] += x29;
                integral_17[g] += x31;
                integral_18[g] += x16;
                integral_19[g] += x26;
                integral_20[g] += x28;
                integral_21[g] += x12*(x14 + x3) + 3*x15;
                integral_22[g] += x32;
                integral_23[g] += x33;
                integral_24[g] += x18;
                integral_25[g] += x28;
                integral_26[g] += x29;
                integral_27[g] += x32;
                integral_28[g] += x33;
                integral_29[g] += x34;
                integral_30[g] += x23;
                integral_31[g] += x29;
                integral_32[g] += x31;
                integral_33[g] += x33;
                integral_34[g] += x34;
                integral_35[g] += x19*(x21 + x3) + 3*x22;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_D_D_SAMECENTER

template void VRR_D_D_SAMECENTER<float>(
        float* integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void VRR_D_D_SAMECENTER<double>(
        double* integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

