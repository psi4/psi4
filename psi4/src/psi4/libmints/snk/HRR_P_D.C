#include "HRR_P_D.h"

template <typename T>
void HRR_P_D(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        T* base_integrals,         //!< overwritten, base [LA+LB|0] integrals for HRR
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
    T* __restrict base_integral_10 = &base_integrals[10*gb_size];
    T* __restrict base_integral_11 = &base_integrals[11*gb_size];
    T* __restrict base_integral_12 = &base_integrals[12*gb_size];
    T* __restrict base_integral_13 = &base_integrals[13*gb_size];
    T* __restrict base_integral_14 = &base_integrals[14*gb_size];
    T* __restrict base_integral_15 = &base_integrals[15*gb_size];

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

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*4*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*4*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*4*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*4*gb_size+3*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x1 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x2 = PBx*scaled_boys_1[g];
                const T x3 = PGx[g]*scaled_boys_2[g];
                const T x4 = x2 - x3;
                const T x5 = PBx*x1 - PGx[g]*x4 + x0;
                const T x6 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x7 = PBy*scaled_boys_1[g];
                const T x8 = PGy[g]*scaled_boys_2[g];
                const T x9 = x7 - x8;
                const T x10 = PBx*x6 - PGx[g]*x9;
                const T x11 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x12 = PBz*scaled_boys_1[g];
                const T x13 = PGz[g]*scaled_boys_2[g];
                const T x14 = x12 - x13;
                const T x15 = PBx*x11 - PGx[g]*x14;
                const T x16 = PBy*x6 - PGy[g]*x9 + x0;
                const T x17 = PBy*x11 - PGy[g]*x14;
                const T x18 = PBz*x11 - PGz[g]*x14 + x0;
                const T x19 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x20 = tpi*(x6 - x7 + x8);
                const T x21 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x22 = tpi*(x11 - x12 + x13);
                const T x23 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x24 = PBy*x9 - PGy[g]*x21 + x19;
                const T x25 = PBy*x14 - PGy[g]*x23;
                const T x26 = PBz*x14 - PGz[g]*x23 + x19;

                base_integral_0[g] = x5;
                base_integral_1[g] = x10;
                base_integral_2[g] = x15;
                base_integral_3[g] = x16;
                base_integral_4[g] = x17;
                base_integral_5[g] = x18;
                base_integral_6[g] = PBx*x5 - PGx[g]*(PBx*x4 - PGx[g]*(PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g]) + x19) + 2*tpi*(x1 - x2 + x3);
                base_integral_7[g] = PBx*x10 - PGx[g]*(PBx*x9 - PGx[g]*x21) + x20;
                base_integral_8[g] = PBx*x15 - PGx[g]*(PBx*x14 - PGx[g]*x23) + x22;
                base_integral_9[g] = PBx*x16 - PGx[g]*x24;
                base_integral_10[g] = PBx*x17 - PGx[g]*x25;
                base_integral_11[g] = PBx*x18 - PGx[g]*x26;
                base_integral_12[g] = PBy*x16 - PGy[g]*x24 + 2*x20;
                base_integral_13[g] = PBy*x17 - PGy[g]*x25 + x22;
                base_integral_14[g] = PBy*x18 - PGy[g]*x26;
                base_integral_15[g] = PBz*x18 - PGz[g]*x26 + 2*x22;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x1 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x2 = PBx*scaled_boys_1[g];
                const T x3 = PGx[g]*scaled_boys_2[g];
                const T x4 = x2 - x3;
                const T x5 = PBx*x1 - PGx[g]*x4 + x0;
                const T x6 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x7 = PBy*scaled_boys_1[g];
                const T x8 = PGy[g]*scaled_boys_2[g];
                const T x9 = x7 - x8;
                const T x10 = PBx*x6 - PGx[g]*x9;
                const T x11 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x12 = PBz*scaled_boys_1[g];
                const T x13 = PGz[g]*scaled_boys_2[g];
                const T x14 = x12 - x13;
                const T x15 = PBx*x11 - PGx[g]*x14;
                const T x16 = PBy*x6 - PGy[g]*x9 + x0;
                const T x17 = PBy*x11 - PGy[g]*x14;
                const T x18 = PBz*x11 - PGz[g]*x14 + x0;
                const T x19 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x20 = tpi*(x6 - x7 + x8);
                const T x21 = PBy*scaled_boys_2[g] - PGy[g]*scaled_boys_3[g];
                const T x22 = tpi*(x11 - x12 + x13);
                const T x23 = PBz*scaled_boys_2[g] - PGz[g]*scaled_boys_3[g];
                const T x24 = PBy*x9 - PGy[g]*x21 + x19;
                const T x25 = PBy*x14 - PGy[g]*x23;
                const T x26 = PBz*x14 - PGz[g]*x23 + x19;

                base_integral_0[g] += x5;
                base_integral_1[g] += x10;
                base_integral_2[g] += x15;
                base_integral_3[g] += x16;
                base_integral_4[g] += x17;
                base_integral_5[g] += x18;
                base_integral_6[g] += PBx*x5 - PGx[g]*(PBx*x4 - PGx[g]*(PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g]) + x19) + 2*tpi*(x1 - x2 + x3);
                base_integral_7[g] += PBx*x10 - PGx[g]*(PBx*x9 - PGx[g]*x21) + x20;
                base_integral_8[g] += PBx*x15 - PGx[g]*(PBx*x14 - PGx[g]*x23) + x22;
                base_integral_9[g] += PBx*x16 - PGx[g]*x24;
                base_integral_10[g] += PBx*x17 - PGx[g]*x25;
                base_integral_11[g] += PBx*x18 - PGx[g]*x26;
                base_integral_12[g] += PBy*x16 - PGy[g]*x24 + 2*x20;
                base_integral_13[g] += PBy*x17 - PGy[g]*x25 + x22;
                base_integral_14[g] += PBy*x18 - PGy[g]*x26;
                base_integral_15[g] += PBz*x18 - PGz[g]*x26 + 2*x22;

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

        const T x0 = -base_integral_7[g];
        const T x1 = -base_integral_8[g];
        const T x2 = -base_integral_9[g];
        const T x3 = -base_integral_10[g];
        const T x4 = -base_integral_11[g];
        const T x5 = -base_integral_13[g];
        const T x6 = -base_integral_14[g];

        integral_0[g] = -ABx*base_integral_0[g] + base_integral_6[g];
        integral_1[g] = -ABx*base_integral_1[g] - x0;
        integral_2[g] = -ABx*base_integral_2[g] - x1;
        integral_3[g] = -ABx*base_integral_3[g] - x2;
        integral_4[g] = -ABx*base_integral_4[g] - x3;
        integral_5[g] = -ABx*base_integral_5[g] - x4;
        integral_6[g] = -ABy*base_integral_0[g] - x0;
        integral_7[g] = -ABy*base_integral_1[g] - x2;
        integral_8[g] = -ABy*base_integral_2[g] - x3;
        integral_9[g] = -ABy*base_integral_3[g] + base_integral_12[g];
        integral_10[g] = -ABy*base_integral_4[g] - x5;
        integral_11[g] = -ABy*base_integral_5[g] - x6;
        integral_12[g] = -ABz*base_integral_0[g] - x1;
        integral_13[g] = -ABz*base_integral_1[g] - x3;
        integral_14[g] = -ABz*base_integral_2[g] - x4;
        integral_15[g] = -ABz*base_integral_3[g] - x5;
        integral_16[g] = -ABz*base_integral_4[g] - x6;
        integral_17[g] = -ABz*base_integral_5[g] + base_integral_15[g];

    } // for g
} // HRR_P_D

template void HRR_P_D<float>(
        float* integrals, 
        float* base_integrals, 
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

template void HRR_P_D<double>(
        double* integrals, 
        double* base_integrals, 
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

