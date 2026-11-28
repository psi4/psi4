#include "VRR_P_D_SAMECENTER.h"

template <typename T>
void VRR_P_D_SAMECENTER(
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

    for (int prim = 0; prim < nprim; ++prim) {

        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];
        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*4*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*4*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*4*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*4*gb_size+3*gb_size];

        const T tpi = 0.5 / v_p[prim];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGz[g]*scaled_boys_3[g];
                const T x1 = PGx[g]*x0;
                const T x2 = -PGy[g]*x1;
                const T x3 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x4 = tpi*x3;
                const T x5 = x4 + PGz[g]*x0;
                const T x6 = -PGx[g]*x5;
                const T x7 = -PGy[g]*x5;
                const T x8 = 2*x4 + x5;
                const T x9 = PGx[g]*scaled_boys_3[g];
                const T x10 = PGx[g]*x9 + x4;
                const T x11 = -PGy[g]*x10;
                const T x12 = -PGz[g]*x10;
                const T x13 = 2*x4 + x10;
                const T x14 = PGy[g]*scaled_boys_3[g];
                const T x15 = x4 + PGy[g]*x14;
                const T x16 = 2*x4 + x15;
                const T x17 = -PGx[g]*x15;
                const T x18 = -PGz[g]*x15;

                integral_0[g] = -PGx[g]*x13;
                integral_1[g] = x11;
                integral_2[g] = x12;
                integral_3[g] = x17;
                integral_4[g] = x2;
                integral_5[g] = x6;
                integral_6[g] = x11;
                integral_7[g] = x17;
                integral_8[g] = x2;
                integral_9[g] = -PGy[g]*x16;
                integral_10[g] = x18;
                integral_11[g] = x7;
                integral_12[g] = x12;
                integral_13[g] = x2;
                integral_14[g] = x6;
                integral_15[g] = x18;
                integral_16[g] = x7;
                integral_17[g] = -PGz[g]*x8;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGz[g]*scaled_boys_3[g];
                const T x1 = PGx[g]*x0;
                const T x2 = -PGy[g]*x1;
                const T x3 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x4 = tpi*x3;
                const T x5 = x4 + PGz[g]*x0;
                const T x6 = -PGx[g]*x5;
                const T x7 = -PGy[g]*x5;
                const T x8 = 2*x4 + x5;
                const T x9 = PGx[g]*scaled_boys_3[g];
                const T x10 = PGx[g]*x9 + x4;
                const T x11 = -PGy[g]*x10;
                const T x12 = -PGz[g]*x10;
                const T x13 = 2*x4 + x10;
                const T x14 = PGy[g]*scaled_boys_3[g];
                const T x15 = x4 + PGy[g]*x14;
                const T x16 = 2*x4 + x15;
                const T x17 = -PGx[g]*x15;
                const T x18 = -PGz[g]*x15;

                integral_0[g] += -PGx[g]*x13;
                integral_1[g] += x11;
                integral_2[g] += x12;
                integral_3[g] += x17;
                integral_4[g] += x2;
                integral_5[g] += x6;
                integral_6[g] += x11;
                integral_7[g] += x17;
                integral_8[g] += x2;
                integral_9[g] += -PGy[g]*x16;
                integral_10[g] += x18;
                integral_11[g] += x7;
                integral_12[g] += x12;
                integral_13[g] += x2;
                integral_14[g] += x6;
                integral_15[g] += x18;
                integral_16[g] += x7;
                integral_17[g] += -PGz[g]*x8;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_P_D_SAMECENTER

template void VRR_P_D_SAMECENTER<float>(
        float* integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void VRR_P_D_SAMECENTER<double>(
        double* integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

