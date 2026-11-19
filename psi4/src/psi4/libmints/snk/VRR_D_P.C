#include "VRR_D_P.h"

template <typename T>
void VRR_D_P(
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

                const T x0 = PBx*scaled_boys_0[g];
                const T x1 = PBx*scaled_boys_1[g];
                const T x2 = PGx[g]*scaled_boys_1[g];
                const T x3 = PGx[g]*scaled_boys_2[g];
                const T x4 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x5 = x1 - x3;
                const T x6 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x7 = PGy[g]*scaled_boys_2[g];
                const T x8 = PBy*scaled_boys_1[g];
                const T x9 = -x8;
                const T x10 = PBy*scaled_boys_0[g];
                const T x11 = PGy[g]*scaled_boys_1[g];
                const T x12 = -x11;
                const T x13 = x10 + x12;
                const T x14 = -x7;
                const T x15 = x14 + x8;
                const T x16 = -PGy[g]*scaled_boys_3[g];
                const T x17 = PBy*scaled_boys_2[g] + x16;
                const T x18 = PGz[g]*scaled_boys_2[g];
                const T x19 = PBz*scaled_boys_1[g];
                const T x20 = -x19;
                const T x21 = PBz*scaled_boys_0[g];
                const T x22 = PGz[g]*scaled_boys_1[g];
                const T x23 = -x22;
                const T x24 = x21 + x23;
                const T x25 = tpi*(x18 + x20 + x24);
                const T x26 = -x18;
                const T x27 = x19 + x26;
                const T x28 = -PGz[g]*scaled_boys_3[g];
                const T x29 = PBz*scaled_boys_2[g] + x28;
                const T x30 = PAy*scaled_boys_1[g];
                const T x31 = -x30;
                const T x32 = PAy*scaled_boys_0[g];
                const T x33 = x12 + x32;
                const T x34 = x14 + x30;
                const T x35 = PAy*scaled_boys_2[g] + x16;
                const T x36 = PAy*x13 - PGy[g]*x15 + x4;
                const T x37 = PAy*x15 - PGy[g]*x17 + x6;
                const T x38 = PAy*x24 - PGy[g]*x27;
                const T x39 = PAy*x27 - PGy[g]*x29;
                const T x40 = PAz*scaled_boys_1[g];
                const T x41 = -x40;
                const T x42 = PAz*scaled_boys_0[g];
                const T x43 = x23 + x42;
                const T x44 = tpi*(x18 + x41 + x43);
                const T x45 = x26 + x40;
                const T x46 = PAz*scaled_boys_2[g] + x28;
                const T x47 = -PGy[g]*x45;
                const T x48 = PBy*x43 + x47;
                const T x49 = -PGy[g]*x46;
                const T x50 = PBy*x45 + x49;
                const T x51 = PAz*x24 - PGz[g]*x27 + x4;
                const T x52 = PAz*x27 - PGz[g]*x29 + x6;
                const T x53 = PAz*x43 - PGz[g]*x45 + x4;
                const T x54 = PAz*x45 - PGz[g]*x46 + x6;

                integral_0[g] = PAx*(PAx*(x0 - x2) - PGx[g]*x5 + x4) - PGx[g]*(PAx*x5 - PGx[g]*(PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g]) + x6) + tpi*(PAx*scaled_boys_0[g] - PAx*scaled_boys_1[g] + x0 - x1 - 2*x2 + 2*x3);
                integral_1[g] = PAx*(PAx*x13 - PGx[g]*x15) - PGx[g]*(PAx*x15 - PGx[g]*x17) + tpi*(x13 + x7 + x9);
                integral_2[g] = PAx*(PAx*x24 - PGx[g]*x27) - PGx[g]*(PAx*x27 - PGx[g]*x29) + x25;
                integral_3[g] = PAx*(PBx*x33 - PGx[g]*x34) - PGx[g]*(PBx*x34 - PGx[g]*x35) + tpi*(x31 + x33 + x7);
                integral_4[g] = PAx*x36 - PGx[g]*x37;
                integral_5[g] = PAx*x38 - PGx[g]*x39;
                integral_6[g] = PAx*(PBx*x43 - PGx[g]*x45) - PGx[g]*(PBx*x45 - PGx[g]*x46) + x44;
                integral_7[g] = PAx*x48 - PGx[g]*x50;
                integral_8[g] = PAx*x51 - PGx[g]*x52;
                integral_9[g] = PBx*(PAy*x33 - PGy[g]*x34 + x4) - PGx[g]*(PAy*x34 - PGy[g]*x35 + x6);
                integral_10[g] = PAy*x36 - PGy[g]*x37 + tpi*(x10 - 2*x11 + x31 + x32 + 2*x7 + x9);
                integral_11[g] = PAy*x38 - PGy[g]*x39 + x25;
                integral_12[g] = PBx*(PAy*x43 + x47) - PGx[g]*(PAy*x45 + x49);
                integral_13[g] = PAy*x48 - PGy[g]*x50 + x44;
                integral_14[g] = PAy*x51 - PGy[g]*x52;
                integral_15[g] = PBx*x53 - PGx[g]*x54;
                integral_16[g] = PBy*x53 - PGy[g]*x54;
                integral_17[g] = PAz*x51 - PGz[g]*x52 + tpi*(2*x18 + x20 + x21 - 2*x22 + x41 + x42);

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PBx*scaled_boys_0[g];
                const T x1 = PBx*scaled_boys_1[g];
                const T x2 = PGx[g]*scaled_boys_1[g];
                const T x3 = PGx[g]*scaled_boys_2[g];
                const T x4 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x5 = x1 - x3;
                const T x6 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x7 = PGy[g]*scaled_boys_2[g];
                const T x8 = PBy*scaled_boys_1[g];
                const T x9 = -x8;
                const T x10 = PBy*scaled_boys_0[g];
                const T x11 = PGy[g]*scaled_boys_1[g];
                const T x12 = -x11;
                const T x13 = x10 + x12;
                const T x14 = -x7;
                const T x15 = x14 + x8;
                const T x16 = -PGy[g]*scaled_boys_3[g];
                const T x17 = PBy*scaled_boys_2[g] + x16;
                const T x18 = PGz[g]*scaled_boys_2[g];
                const T x19 = PBz*scaled_boys_1[g];
                const T x20 = -x19;
                const T x21 = PBz*scaled_boys_0[g];
                const T x22 = PGz[g]*scaled_boys_1[g];
                const T x23 = -x22;
                const T x24 = x21 + x23;
                const T x25 = tpi*(x18 + x20 + x24);
                const T x26 = -x18;
                const T x27 = x19 + x26;
                const T x28 = -PGz[g]*scaled_boys_3[g];
                const T x29 = PBz*scaled_boys_2[g] + x28;
                const T x30 = PAy*scaled_boys_1[g];
                const T x31 = -x30;
                const T x32 = PAy*scaled_boys_0[g];
                const T x33 = x12 + x32;
                const T x34 = x14 + x30;
                const T x35 = PAy*scaled_boys_2[g] + x16;
                const T x36 = PAy*x13 - PGy[g]*x15 + x4;
                const T x37 = PAy*x15 - PGy[g]*x17 + x6;
                const T x38 = PAy*x24 - PGy[g]*x27;
                const T x39 = PAy*x27 - PGy[g]*x29;
                const T x40 = PAz*scaled_boys_1[g];
                const T x41 = -x40;
                const T x42 = PAz*scaled_boys_0[g];
                const T x43 = x23 + x42;
                const T x44 = tpi*(x18 + x41 + x43);
                const T x45 = x26 + x40;
                const T x46 = PAz*scaled_boys_2[g] + x28;
                const T x47 = -PGy[g]*x45;
                const T x48 = PBy*x43 + x47;
                const T x49 = -PGy[g]*x46;
                const T x50 = PBy*x45 + x49;
                const T x51 = PAz*x24 - PGz[g]*x27 + x4;
                const T x52 = PAz*x27 - PGz[g]*x29 + x6;
                const T x53 = PAz*x43 - PGz[g]*x45 + x4;
                const T x54 = PAz*x45 - PGz[g]*x46 + x6;

                integral_0[g] += PAx*(PAx*(x0 - x2) - PGx[g]*x5 + x4) - PGx[g]*(PAx*x5 - PGx[g]*(PBx*scaled_boys_2[g] - PGx[g]*scaled_boys_3[g]) + x6) + tpi*(PAx*scaled_boys_0[g] - PAx*scaled_boys_1[g] + x0 - x1 - 2*x2 + 2*x3);
                integral_1[g] += PAx*(PAx*x13 - PGx[g]*x15) - PGx[g]*(PAx*x15 - PGx[g]*x17) + tpi*(x13 + x7 + x9);
                integral_2[g] += PAx*(PAx*x24 - PGx[g]*x27) - PGx[g]*(PAx*x27 - PGx[g]*x29) + x25;
                integral_3[g] += PAx*(PBx*x33 - PGx[g]*x34) - PGx[g]*(PBx*x34 - PGx[g]*x35) + tpi*(x31 + x33 + x7);
                integral_4[g] += PAx*x36 - PGx[g]*x37;
                integral_5[g] += PAx*x38 - PGx[g]*x39;
                integral_6[g] += PAx*(PBx*x43 - PGx[g]*x45) - PGx[g]*(PBx*x45 - PGx[g]*x46) + x44;
                integral_7[g] += PAx*x48 - PGx[g]*x50;
                integral_8[g] += PAx*x51 - PGx[g]*x52;
                integral_9[g] += PBx*(PAy*x33 - PGy[g]*x34 + x4) - PGx[g]*(PAy*x34 - PGy[g]*x35 + x6);
                integral_10[g] += PAy*x36 - PGy[g]*x37 + tpi*(x10 - 2*x11 + x31 + x32 + 2*x7 + x9);
                integral_11[g] += PAy*x38 - PGy[g]*x39 + x25;
                integral_12[g] += PBx*(PAy*x43 + x47) - PGx[g]*(PAy*x45 + x49);
                integral_13[g] += PAy*x48 - PGy[g]*x50 + x44;
                integral_14[g] += PAy*x51 - PGy[g]*x52;
                integral_15[g] += PBx*x53 - PGx[g]*x54;
                integral_16[g] += PBy*x53 - PGy[g]*x54;
                integral_17[g] += PAz*x51 - PGz[g]*x52 + tpi*(2*x18 + x20 + x21 - 2*x22 + x41 + x42);

            } // for g
        } // if/else prim
    } // for prim
} // VRR_D_P

template void VRR_D_P<float>(
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

template void VRR_D_P<double>(
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

