#include "VRR_P_P.h"

template <typename T>
void VRR_P_P(
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
    T* __restrict integral_6 = &integrals[6*gb_size];
    T* __restrict integral_7 = &integrals[7*gb_size];
    T* __restrict integral_8 = &integrals[8*gb_size];

    for (int prim = 0; prim < nprim; ++prim) {

        const T tpi = 0.5 / v_p[prim];
        const T PAx = v_PAx[prim];
        const T PAy = v_PAy[prim];
        const T PAz = v_PAz[prim];
        const T PBx = PAx + ABx;
        const T PBy = PAy + ABy;
        const T PBz = PAz + ABz;
        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*3*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*3*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*3*gb_size+2*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x1 = -PGy[g]*scaled_boys_1[g];
                const T x2 = PBy*scaled_boys_0[g] + x1;
                const T x3 = -PGy[g]*scaled_boys_2[g];
                const T x4 = PBy*scaled_boys_1[g] + x3;
                const T x5 = -PGz[g]*scaled_boys_1[g];
                const T x6 = PBz*scaled_boys_0[g] + x5;
                const T x7 = -PGz[g]*scaled_boys_2[g];
                const T x8 = PBz*scaled_boys_1[g] + x7;
                const T x9 = PAz*scaled_boys_0[g] + x5;
                const T x10 = PAz*scaled_boys_1[g] + x7;

                integral_0[g] = PAx*(PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g]) - PGx[g]*(PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g]) + x0;
                integral_1[g] = PAx*x2 - PGx[g]*x4;
                integral_2[g] = PAx*x6 - PGx[g]*x8;
                integral_3[g] = PBx*(PAy*scaled_boys_0[g] + x1) - PGx[g]*(PAy*scaled_boys_1[g] + x3);
                integral_4[g] = PAy*x2 - PGy[g]*x4 + x0;
                integral_5[g] = PAy*x6 - PGy[g]*x8;
                integral_6[g] = PBx*x9 - PGx[g]*x10;
                integral_7[g] = PBy*x9 - PGy[g]*x10;
                integral_8[g] = PAz*x6 - PGz[g]*x8 + x0;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x1 = -PGy[g]*scaled_boys_1[g];
                const T x2 = PBy*scaled_boys_0[g] + x1;
                const T x3 = -PGy[g]*scaled_boys_2[g];
                const T x4 = PBy*scaled_boys_1[g] + x3;
                const T x5 = -PGz[g]*scaled_boys_1[g];
                const T x6 = PBz*scaled_boys_0[g] + x5;
                const T x7 = -PGz[g]*scaled_boys_2[g];
                const T x8 = PBz*scaled_boys_1[g] + x7;
                const T x9 = PAz*scaled_boys_0[g] + x5;
                const T x10 = PAz*scaled_boys_1[g] + x7;

                integral_0[g] += PAx*(PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g]) - PGx[g]*(PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g]) + x0;
                integral_1[g] += PAx*x2 - PGx[g]*x4;
                integral_2[g] += PAx*x6 - PGx[g]*x8;
                integral_3[g] += PBx*(PAy*scaled_boys_0[g] + x1) - PGx[g]*(PAy*scaled_boys_1[g] + x3);
                integral_4[g] += PAy*x2 - PGy[g]*x4 + x0;
                integral_5[g] += PAy*x6 - PGy[g]*x8;
                integral_6[g] += PBx*x9 - PGx[g]*x10;
                integral_7[g] += PBy*x9 - PGy[g]*x10;
                integral_8[g] += PAz*x6 - PGz[g]*x8 + x0;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_P_P

template void VRR_P_P<float>(
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

template void VRR_P_P<double>(
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

