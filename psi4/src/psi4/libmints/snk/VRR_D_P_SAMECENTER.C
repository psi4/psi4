#include "VRR_D_P_SAMECENTER.h"

template <typename T>
void VRR_D_P_SAMECENTER(
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
        const T tpi_2 = 2*tpi;

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = -PGz[g]*scaled_boys_3[g];
                const T x1 = -PGx[g]*x0;
                const T x2 = -PGy[g]*x1;
                const T x3 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x4 = tpi*x3 - PGz[g]*x0;
                const T x5 = -PGx[g]*x4;
                const T x6 = -PGy[g]*x4;
                const T x7 = -tpi_2*x3 - x4;
                const T x8 = -PGx[g]*scaled_boys_3[g];
                const T x9 = tpi*x3 - PGx[g]*x8;
                const T x10 = -PGy[g]*x9;
                const T x11 = -PGz[g]*x9;
                const T x12 = -tpi_2*x3 - x9;
                const T x13 = -PGy[g]*scaled_boys_3[g];
                const T x14 = tpi*x3 - PGy[g]*x13;
                const T x15 = -tpi_2*x3 - x14;
                const T x16 = -PGx[g]*x14;
                const T x17 = -PGz[g]*x14;

                integral_0[g] = PGx[g]*x12;
                integral_1[g] = x10;
                integral_2[g] = x11;
                integral_3[g] = x10;
                integral_4[g] = x16;
                integral_5[g] = x2;
                integral_6[g] = x11;
                integral_7[g] = x2;
                integral_8[g] = x5;
                integral_9[g] = x16;
                integral_10[g] = PGy[g]*x15;
                integral_11[g] = x17;
                integral_12[g] = x2;
                integral_13[g] = x17;
                integral_14[g] = x6;
                integral_15[g] = x5;
                integral_16[g] = x6;
                integral_17[g] = PGz[g]*x7;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = -PGz[g]*scaled_boys_3[g];
                const T x1 = -PGx[g]*x0;
                const T x2 = -PGy[g]*x1;
                const T x3 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x4 = tpi*x3 - PGz[g]*x0;
                const T x5 = -PGx[g]*x4;
                const T x6 = -PGy[g]*x4;
                const T x7 = -tpi_2*x3 - x4;
                const T x8 = -PGx[g]*scaled_boys_3[g];
                const T x9 = tpi*x3 - PGx[g]*x8;
                const T x10 = -PGy[g]*x9;
                const T x11 = -PGz[g]*x9;
                const T x12 = -tpi_2*x3 - x9;
                const T x13 = -PGy[g]*scaled_boys_3[g];
                const T x14 = tpi*x3 - PGy[g]*x13;
                const T x15 = -tpi_2*x3 - x14;
                const T x16 = -PGx[g]*x14;
                const T x17 = -PGz[g]*x14;

                integral_0[g] += PGx[g]*x12;
                integral_1[g] += x10;
                integral_2[g] += x11;
                integral_3[g] += x10;
                integral_4[g] += x16;
                integral_5[g] += x2;
                integral_6[g] += x11;
                integral_7[g] += x2;
                integral_8[g] += x5;
                integral_9[g] += x16;
                integral_10[g] += PGy[g]*x15;
                integral_11[g] += x17;
                integral_12[g] += x2;
                integral_13[g] += x17;
                integral_14[g] += x6;
                integral_15[g] += x5;
                integral_16[g] += x6;
                integral_17[g] += PGz[g]*x7;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_D_P_SAMECENTER

template void VRR_D_P_SAMECENTER<float>(
        float* integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void VRR_D_P_SAMECENTER<double>(
        double* integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

