#include "VRR_S_F.h"

template <typename T>
void VRR_S_F(
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

        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];
        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*4*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*4*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*4*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*4*gb_size+3*gb_size];

        const T PBx = v_PAx[prim] + ABx;
        const T PBy = v_PAy[prim] + ABy;
        const T PBz = v_PAz[prim] + ABz;
        const T tpi = 0.5 / v_p[prim];
        const T tpi_2 = 2*tpi;

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                const T x2 = PBx*x0 - PGx[g]*x1;
                const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x4 = PBx*x1 - PGx[g]*x3;
                const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                const T x6 = tpi*x5 + PBz*x0 - PGz[g]*x1;
                const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x8 = tpi*x7 + PBz*x1 - PGz[g]*x3;
                const T x9 = x0 - x1;
                const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                const T x12 = PBx*x10 + tpi*x5 - PGx[g]*x11;
                const T x13 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                const T x14 = PBx*x11 + tpi*x7 - PGx[g]*x13;
                const T x15 = x10 - x11;
                const T x16 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x17 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                const T x18 = tpi*x5 + PBy*x16 - PGy[g]*x17;
                const T x19 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x20 = tpi*x7 + PBy*x17 - PGy[g]*x19;
                const T x21 = x16 - x17;

                integral_0[g] = PBx*x12 + tpi_2*x15 - PGx[g]*x14;
                integral_1[g] = PBy*x12 - PGy[g]*x14;
                integral_2[g] = PBz*x12 - PGz[g]*x14;
                integral_3[g] = PBx*x18 - PGx[g]*x20;
                integral_4[g] = PBy*x2 - PGy[g]*x4;
                integral_5[g] = PBx*x6 - PGx[g]*x8;
                integral_6[g] = PBy*x18 + tpi_2*x21 - PGy[g]*x20;
                integral_7[g] = PBz*x18 - PGz[g]*x20;
                integral_8[g] = PBy*x6 - PGy[g]*x8;
                integral_9[g] = PBz*x6 + tpi_2*x9 - PGz[g]*x8;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x1 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                const T x2 = PBx*x0 - PGx[g]*x1;
                const T x3 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x4 = PBx*x1 - PGx[g]*x3;
                const T x5 = scaled_boys_0[g] - scaled_boys_1[g];
                const T x6 = tpi*x5 + PBz*x0 - PGz[g]*x1;
                const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x8 = tpi*x7 + PBz*x1 - PGz[g]*x3;
                const T x9 = x0 - x1;
                const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                const T x12 = PBx*x10 + tpi*x5 - PGx[g]*x11;
                const T x13 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                const T x14 = PBx*x11 + tpi*x7 - PGx[g]*x13;
                const T x15 = x10 - x11;
                const T x16 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x17 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                const T x18 = tpi*x5 + PBy*x16 - PGy[g]*x17;
                const T x19 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x20 = tpi*x7 + PBy*x17 - PGy[g]*x19;
                const T x21 = x16 - x17;

                integral_0[g] += PBx*x12 + tpi_2*x15 - PGx[g]*x14;
                integral_1[g] += PBy*x12 - PGy[g]*x14;
                integral_2[g] += PBz*x12 - PGz[g]*x14;
                integral_3[g] += PBx*x18 - PGx[g]*x20;
                integral_4[g] += PBy*x2 - PGy[g]*x4;
                integral_5[g] += PBx*x6 - PGx[g]*x8;
                integral_6[g] += PBy*x18 + tpi_2*x21 - PGy[g]*x20;
                integral_7[g] += PBz*x18 - PGz[g]*x20;
                integral_8[g] += PBy*x6 - PGy[g]*x8;
                integral_9[g] += PBz*x6 + tpi_2*x9 - PGz[g]*x8;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_S_F

template void VRR_S_F<float>(
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

template void VRR_S_F<double>(
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

