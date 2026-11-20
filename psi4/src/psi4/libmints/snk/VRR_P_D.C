#include "VRR_P_D.h"

template <typename T>
void VRR_P_D(
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

        const T PAx = v_PAx[prim];
        const T PAy = v_PAy[prim];
        const T PAz = v_PAz[prim];
        const T PBx = v_PAx[prim] + ABx;
        const T PBy = v_PAy[prim] + ABy;
        const T PBz = v_PAz[prim] + ABz;
        const T tpi = 0.5 / v_p[prim];
        const T tpi_2 = 2*tpi;

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                const T x1 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x2 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                const T x3 = tpi*x0 + PBy*x1 - PGy[g]*x2;
                const T x4 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x5 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x6 = tpi*x4 + PBy*x2 - PGy[g]*x5;
                const T x7 = PBx*x1 - PGx[g]*x2;
                const T x8 = PBx*x2 - PGx[g]*x5;
                const T x9 = x1 - x2;
                const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                const T x12 = x10 - x11;
                const T x13 = PBx*x10 + tpi*x0 - PGx[g]*x11;
                const T x14 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                const T x15 = PBx*x11 + tpi*x4 - PGx[g]*x14;
                const T x16 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x17 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                const T x18 = tpi*x0 + PBz*x16 - PGz[g]*x17;
                const T x19 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x20 = tpi*x4 + PBz*x17 - PGz[g]*x19;
                const T x21 = x16 - x17;
                const T x22 = PBx*x16 - PGx[g]*x17;
                const T x23 = PBx*x17 - PGx[g]*x19;
                const T x24 = PBy*x16 - PGy[g]*x17;
                const T x25 = PBy*x17 - PGy[g]*x19;

                integral_0[g] = PAx*x13 + tpi_2*x12 - PGx[g]*x15;
                integral_1[g] = PAx*x7 + tpi*x9 - PGx[g]*x8;
                integral_2[g] = PAx*x22 + tpi*x21 - PGx[g]*x23;
                integral_3[g] = PAx*x3 - PGx[g]*x6;
                integral_4[g] = PAx*x24 - PGx[g]*x25;
                integral_5[g] = PAx*x18 - PGx[g]*x20;
                integral_6[g] = PAy*x13 - PGy[g]*x15;
                integral_7[g] = tpi*x12 + PAy*x7 - PGy[g]*x8;
                integral_8[g] = PAy*x22 - PGy[g]*x23;
                integral_9[g] = tpi_2*x9 + PAy*x3 - PGy[g]*x6;
                integral_10[g] = tpi*x21 + PAy*x24 - PGy[g]*x25;
                integral_11[g] = PAy*x18 - PGy[g]*x20;
                integral_12[g] = PAz*x13 - PGz[g]*x15;
                integral_13[g] = PAz*x7 - PGz[g]*x8;
                integral_14[g] = tpi*x12 + PAz*x22 - PGz[g]*x23;
                integral_15[g] = PAz*x3 - PGz[g]*x6;
                integral_16[g] = tpi*x9 + PAz*x24 - PGz[g]*x25;
                integral_17[g] = tpi_2*x21 + PAz*x18 - PGz[g]*x20;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = scaled_boys_0[g] - scaled_boys_1[g];
                const T x1 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x2 = PBy*scaled_boys_1[g] - PGy[g]*scaled_boys_2[g];
                const T x3 = tpi*x0 + PBy*x1 - PGy[g]*x2;
                const T x4 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x5 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x6 = tpi*x4 + PBy*x2 - PGy[g]*x5;
                const T x7 = PBx*x1 - PGx[g]*x2;
                const T x8 = PBx*x2 - PGx[g]*x5;
                const T x9 = x1 - x2;
                const T x10 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x11 = PBx*scaled_boys_1[g] - PGx[g]*scaled_boys_2[g];
                const T x12 = x10 - x11;
                const T x13 = PBx*x10 + tpi*x0 - PGx[g]*x11;
                const T x14 = PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g];
                const T x15 = PBx*x11 + tpi*x4 - PGx[g]*x14;
                const T x16 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x17 = PBz*scaled_boys_1[g] - PGz[g]*scaled_boys_2[g];
                const T x18 = tpi*x0 + PBz*x16 - PGz[g]*x17;
                const T x19 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x20 = tpi*x4 + PBz*x17 - PGz[g]*x19;
                const T x21 = x16 - x17;
                const T x22 = PBx*x16 - PGx[g]*x17;
                const T x23 = PBx*x17 - PGx[g]*x19;
                const T x24 = PBy*x16 - PGy[g]*x17;
                const T x25 = PBy*x17 - PGy[g]*x19;

                integral_0[g] += PAx*x13 + tpi_2*x12 - PGx[g]*x15;
                integral_1[g] += PAx*x7 + tpi*x9 - PGx[g]*x8;
                integral_2[g] += PAx*x22 + tpi*x21 - PGx[g]*x23;
                integral_3[g] += PAx*x3 - PGx[g]*x6;
                integral_4[g] += PAx*x24 - PGx[g]*x25;
                integral_5[g] += PAx*x18 - PGx[g]*x20;
                integral_6[g] += PAy*x13 - PGy[g]*x15;
                integral_7[g] += tpi*x12 + PAy*x7 - PGy[g]*x8;
                integral_8[g] += PAy*x22 - PGy[g]*x23;
                integral_9[g] += tpi_2*x9 + PAy*x3 - PGy[g]*x6;
                integral_10[g] += tpi*x21 + PAy*x24 - PGy[g]*x25;
                integral_11[g] += PAy*x18 - PGy[g]*x20;
                integral_12[g] += PAz*x13 - PGz[g]*x15;
                integral_13[g] += PAz*x7 - PGz[g]*x8;
                integral_14[g] += tpi*x12 + PAz*x22 - PGz[g]*x23;
                integral_15[g] += PAz*x3 - PGz[g]*x6;
                integral_16[g] += tpi*x9 + PAz*x24 - PGz[g]*x25;
                integral_17[g] += tpi_2*x21 + PAz*x18 - PGz[g]*x20;

            } // for g
        } // if/else prim
    } // for prim
} // VRR_P_D

template void VRR_P_D<float>(
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

template void VRR_P_D<double>(
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

