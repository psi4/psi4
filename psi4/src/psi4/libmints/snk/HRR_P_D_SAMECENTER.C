#include "HRR_P_D_SAMECENTER.h"

template <typename T>
void HRR_P_D_SAMECENTER(
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

    for (int prim = 0; prim < nprim; ++prim) {

        const T tpi = 0.5 / v_p[prim];
        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*4*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*4*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*4*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*4*gb_size+3*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g]*scaled_boys_3[g];
                const T x1 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x2 = 3*x1;
                const T x3 = x0 + x1;
                const T x4 = PGy[g]*PGy[g]*scaled_boys_3[g];
                const T x5 = x1 + x4;
                const T x6 = PGz[g]*PGz[g]*scaled_boys_3[g];
                const T x7 = x1 + x6;

                base_integral_0[g] = -PGx[g]*(x0 + x2);
                base_integral_1[g] = -PGy[g]*x3;
                base_integral_2[g] = -PGz[g]*x3;
                base_integral_3[g] = -PGx[g]*x5;
                base_integral_4[g] = -PGx[g]*PGy[g]*PGz[g]*scaled_boys_3[g];
                base_integral_5[g] = -PGx[g]*x7;
                base_integral_6[g] = -PGy[g]*(x2 + x4);
                base_integral_7[g] = -PGz[g]*x5;
                base_integral_8[g] = -PGy[g]*x7;
                base_integral_9[g] = -PGz[g]*(x2 + x6);

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g]*scaled_boys_3[g];
                const T x1 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x2 = 3*x1;
                const T x3 = x0 + x1;
                const T x4 = PGy[g]*PGy[g]*scaled_boys_3[g];
                const T x5 = x1 + x4;
                const T x6 = PGz[g]*PGz[g]*scaled_boys_3[g];
                const T x7 = x1 + x6;

                base_integral_0[g] += -PGx[g]*(x0 + x2);
                base_integral_1[g] += -PGy[g]*x3;
                base_integral_2[g] += -PGz[g]*x3;
                base_integral_3[g] += -PGx[g]*x5;
                base_integral_4[g] += -PGx[g]*PGy[g]*PGz[g]*scaled_boys_3[g];
                base_integral_5[g] += -PGx[g]*x7;
                base_integral_6[g] += -PGy[g]*(x2 + x4);
                base_integral_7[g] += -PGz[g]*x5;
                base_integral_8[g] += -PGy[g]*x7;
                base_integral_9[g] += -PGz[g]*(x2 + x6);

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

    } // for g
} // HRR_P_D_SAMECENTER

template void HRR_P_D_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_P_D_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

