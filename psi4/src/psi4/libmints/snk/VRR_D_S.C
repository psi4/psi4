#include "VRR_D_S.h"

template <typename T>
void VRR_D_S(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
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
    T* __restrict integral_0 = &integrals[0*gb_size];
    T* __restrict integral_1 = &integrals[1*gb_size];
    T* __restrict integral_2 = &integrals[2*gb_size];
    T* __restrict integral_3 = &integrals[3*gb_size];
    T* __restrict integral_4 = &integrals[4*gb_size];
    T* __restrict integral_5 = &integrals[5*gb_size];

    for (int prim = 0; prim < nprim; ++prim) {

        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];
        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*3*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*3*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*3*gb_size+2*gb_size];

        const T PAx = v_PAx[prim];
        const T PAy = v_PAy[prim];
        const T PAz = v_PAz[prim];
        const T tpi = 0.5 / v_p[prim];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x1 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                const T x3 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x4 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                const T x5 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x6 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];

                integral_0[g] = PAx*x3 + tpi*x2 - PGx[g]*x4;
                integral_1[g] = PAx*x0 - PGx[g]*x1;
                integral_2[g] = PAx*x5 - PGx[g]*x6;
                integral_3[g] = tpi*x2 + PAy*x0 - PGy[g]*x1;
                integral_4[g] = PAy*x5 - PGy[g]*x6;
                integral_5[g] = tpi*x2 + PAz*x5 - PGz[g]*x6;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x1 = PAy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                const T x2 = scaled_boys_0[g] - scaled_boys_1[g];
                const T x3 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x4 = PAx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                const T x5 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x6 = PAz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];

                integral_0[g] += PAx*x3 + tpi*x2 - PGx[g]*x4;
                integral_1[g] += PAx*x0 - PGx[g]*x1;
                integral_2[g] += PAx*x5 - PGx[g]*x6;
                integral_3[g] += tpi*x2 + PAy*x0 - PGy[g]*x1;
                integral_4[g] += PAy*x5 - PGy[g]*x6;
                integral_5[g] += tpi*x2 + PAz*x5 - PGz[g]*x6;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_D_S

template void VRR_D_S<float>(
        float* integrals, 
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

template void VRR_D_S<double>(
        double* integrals, 
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

