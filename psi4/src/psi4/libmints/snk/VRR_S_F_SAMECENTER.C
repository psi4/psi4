#include "VRR_S_F_SAMECENTER.h"

template <typename T>
void VRR_S_F_SAMECENTER(
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
                const T x2 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x3 = tpi*x2;
                const T x4 = x3 + PGz[g]*x0;
                const T x5 = 2*x3 + x4;
                const T x6 = PGx[g]*scaled_boys_3[g];
                const T x7 = PGx[g]*x6 + x3;
                const T x8 = 2*x3 + x7;
                const T x9 = PGy[g]*scaled_boys_3[g];
                const T x10 = x3 + PGy[g]*x9;
                const T x11 = 2*x3 + x10;

                integral_0[g] = -PGx[g]*x8;
                integral_1[g] = -PGy[g]*x7;
                integral_2[g] = -PGz[g]*x7;
                integral_3[g] = -PGx[g]*x10;
                integral_4[g] = -PGy[g]*x1;
                integral_5[g] = -PGx[g]*x4;
                integral_6[g] = -PGy[g]*x11;
                integral_7[g] = -PGz[g]*x10;
                integral_8[g] = -PGy[g]*x4;
                integral_9[g] = -PGz[g]*x5;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGz[g]*scaled_boys_3[g];
                const T x1 = PGx[g]*x0;
                const T x2 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x3 = tpi*x2;
                const T x4 = x3 + PGz[g]*x0;
                const T x5 = 2*x3 + x4;
                const T x6 = PGx[g]*scaled_boys_3[g];
                const T x7 = PGx[g]*x6 + x3;
                const T x8 = 2*x3 + x7;
                const T x9 = PGy[g]*scaled_boys_3[g];
                const T x10 = x3 + PGy[g]*x9;
                const T x11 = 2*x3 + x10;

                integral_0[g] += -PGx[g]*x8;
                integral_1[g] += -PGy[g]*x7;
                integral_2[g] += -PGz[g]*x7;
                integral_3[g] += -PGx[g]*x10;
                integral_4[g] += -PGy[g]*x1;
                integral_5[g] += -PGx[g]*x4;
                integral_6[g] += -PGy[g]*x11;
                integral_7[g] += -PGz[g]*x10;
                integral_8[g] += -PGy[g]*x4;
                integral_9[g] += -PGz[g]*x5;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_S_F_SAMECENTER

template void VRR_S_F_SAMECENTER<float>(
        float* integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void VRR_S_F_SAMECENTER<double>(
        double* integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

