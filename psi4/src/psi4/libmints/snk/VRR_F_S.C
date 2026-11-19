#include "VRR_F_S.h"

template <typename T>
void VRR_F_S(
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
    T* __restrict integral_9 = &integrals[9*gb_size];

    for (int prim = 0; prim < nprim; ++prim) {

        const T tpi = 0.5 / v_p[prim];
        const T PAx = v_PAx[prim];
        const T PAy = v_PAy[prim];
        const T PAz = v_PAz[prim];
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

                const T x0 = PGx[g]*scaled_boys_2[g];
                const T x1 = PAx*scaled_boys_1[g];
                const T x2 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x3 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x4 = -x0 + x1;
                const T x5 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x6 = PGy[g]*scaled_boys_2[g];
                const T x7 = PAy*scaled_boys_1[g];
                const T x8 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x9 = tpi*(x6 - x7 + x8);
                const T x10 = -x6 + x7;
                const T x11 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x12 = PGz[g]*scaled_boys_2[g];
                const T x13 = PAz*scaled_boys_1[g];
                const T x14 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x15 = tpi*(x12 - x13 + x14);
                const T x16 = -x12 + x13;
                const T x17 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x18 = PAy*x8 - PGy[g]*x10 + x3;
                const T x19 = PAy*x10 - PGy[g]*x11 + x5;
                const T x20 = PAy*x14 - PGy[g]*x16;
                const T x21 = PAy*x16 - PGy[g]*x17;
                const T x22 = PAz*x14 - PGz[g]*x16 + x3;
                const T x23 = PAz*x16 - PGz[g]*x17 + x5;

                integral_0[g] = PAx*(PAx*x2 - PGx[g]*x4 + x3) - PGx[g]*(PAx*x4 - PGx[g]*(PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g]) + x5) + 2*tpi*(x0 - x1 + x2);
                integral_1[g] = PAx*(PAx*x8 - PGx[g]*x10) - PGx[g]*(PAx*x10 - PGx[g]*x11) + x9;
                integral_2[g] = PAx*(PAx*x14 - PGx[g]*x16) - PGx[g]*(PAx*x16 - PGx[g]*x17) + x15;
                integral_3[g] = PAx*x18 - PGx[g]*x19;
                integral_4[g] = PAx*x20 - PGx[g]*x21;
                integral_5[g] = PAx*x22 - PGx[g]*x23;
                integral_6[g] = PAy*x18 - PGy[g]*x19 + 2*x9;
                integral_7[g] = PAy*x20 - PGy[g]*x21 + x15;
                integral_8[g] = PAy*x22 - PGy[g]*x23;
                integral_9[g] = PAz*x22 - PGz[g]*x23 + 2*x15;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*scaled_boys_2[g];
                const T x1 = PAx*scaled_boys_1[g];
                const T x2 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x3 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x4 = -x0 + x1;
                const T x5 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x6 = PGy[g]*scaled_boys_2[g];
                const T x7 = PAy*scaled_boys_1[g];
                const T x8 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x9 = tpi*(x6 - x7 + x8);
                const T x10 = -x6 + x7;
                const T x11 = PAy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x12 = PGz[g]*scaled_boys_2[g];
                const T x13 = PAz*scaled_boys_1[g];
                const T x14 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x15 = tpi*(x12 - x13 + x14);
                const T x16 = -x12 + x13;
                const T x17 = PAz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x18 = PAy*x8 - PGy[g]*x10 + x3;
                const T x19 = PAy*x10 - PGy[g]*x11 + x5;
                const T x20 = PAy*x14 - PGy[g]*x16;
                const T x21 = PAy*x16 - PGy[g]*x17;
                const T x22 = PAz*x14 - PGz[g]*x16 + x3;
                const T x23 = PAz*x16 - PGz[g]*x17 + x5;

                integral_0[g] += PAx*(PAx*x2 - PGx[g]*x4 + x3) - PGx[g]*(PAx*x4 - PGx[g]*(PAx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g]) + x5) + 2*tpi*(x0 - x1 + x2);
                integral_1[g] += PAx*(PAx*x8 - PGx[g]*x10) - PGx[g]*(PAx*x10 - PGx[g]*x11) + x9;
                integral_2[g] += PAx*(PAx*x14 - PGx[g]*x16) - PGx[g]*(PAx*x16 - PGx[g]*x17) + x15;
                integral_3[g] += PAx*x18 - PGx[g]*x19;
                integral_4[g] += PAx*x20 - PGx[g]*x21;
                integral_5[g] += PAx*x22 - PGx[g]*x23;
                integral_6[g] += PAy*x18 - PGy[g]*x19 + 2*x9;
                integral_7[g] += PAy*x20 - PGy[g]*x21 + x15;
                integral_8[g] += PAy*x22 - PGy[g]*x23;
                integral_9[g] += PAz*x22 - PGz[g]*x23 + 2*x15;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_F_S

template void VRR_F_S<float>(
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

template void VRR_F_S<double>(
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

