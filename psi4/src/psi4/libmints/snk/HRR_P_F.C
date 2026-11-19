#include "HRR_P_F.h"

template <typename T>
void HRR_P_F(
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
    T* __restrict base_integral_16 = &base_integrals[16*gb_size];
    T* __restrict base_integral_17 = &base_integrals[17*gb_size];
    T* __restrict base_integral_18 = &base_integrals[18*gb_size];
    T* __restrict base_integral_19 = &base_integrals[19*gb_size];
    T* __restrict base_integral_20 = &base_integrals[20*gb_size];
    T* __restrict base_integral_21 = &base_integrals[21*gb_size];
    T* __restrict base_integral_22 = &base_integrals[22*gb_size];
    T* __restrict base_integral_23 = &base_integrals[23*gb_size];
    T* __restrict base_integral_24 = &base_integrals[24*gb_size];

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

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*5*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*5*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*5*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*5*gb_size+3*gb_size];
        const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*5*gb_size+4*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*scaled_boys_2[g];
                const T x1 = PBx*scaled_boys_1[g];
                const T x2 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x3 = 2*tpi;
                const T x4 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x5 = -x0 + x1;
                const T x6 = PBx*x2 - PGx[g]*x5 + x4;
                const T x7 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x8 = PBx*x5;
                const T x9 = PBx*scaled_boys_2[g];
                const T x10 = PGx[g]*scaled_boys_3[g];
                const T x11 = -x10 + x9;
                const T x12 = PGx[g]*x11;
                const T x13 = -x12 + x7 + x8;
                const T x14 = PBx*x6 - PGx[g]*x13 + x3*(x0 - x1 + x2);
                const T x15 = PGy[g]*scaled_boys_2[g];
                const T x16 = PBy*scaled_boys_1[g];
                const T x17 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x18 = tpi*(x15 - x16 + x17);
                const T x19 = -x15 + x16;
                const T x20 = PBx*x17 - PGx[g]*x19;
                const T x21 = PBx*x19;
                const T x22 = PBy*scaled_boys_2[g];
                const T x23 = PGy[g]*scaled_boys_3[g];
                const T x24 = x22 - x23;
                const T x25 = PGx[g]*x24;
                const T x26 = x21 - x25;
                const T x27 = PBx*x20 - PGx[g]*x26 + x18;
                const T x28 = PGz[g]*scaled_boys_2[g];
                const T x29 = PBz*scaled_boys_1[g];
                const T x30 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x31 = tpi*(x28 - x29 + x30);
                const T x32 = -x28 + x29;
                const T x33 = PBx*x30 - PGx[g]*x32;
                const T x34 = PBx*x32;
                const T x35 = PBz*scaled_boys_2[g];
                const T x36 = PGz[g]*scaled_boys_3[g];
                const T x37 = x35 - x36;
                const T x38 = PGx[g]*x37;
                const T x39 = x34 - x38;
                const T x40 = PBx*x33 - PGx[g]*x39 + x31;
                const T x41 = PBy*x17 - PGy[g]*x19 + x4;
                const T x42 = PBy*x19;
                const T x43 = PGy[g]*x24;
                const T x44 = x42 - x43 + x7;
                const T x45 = PBx*x41 - PGx[g]*x44;
                const T x46 = PBy*x30 - PGy[g]*x32;
                const T x47 = PBy*x32;
                const T x48 = PGy[g]*x37;
                const T x49 = x47 - x48;
                const T x50 = PBx*x46 - PGx[g]*x49;
                const T x51 = PBz*x30 - PGz[g]*x32 + x4;
                const T x52 = PBz*x32;
                const T x53 = PGz[g]*x37;
                const T x54 = x52 - x53 + x7;
                const T x55 = PBx*x51 - PGx[g]*x54;
                const T x56 = PBy*x41 - PGy[g]*x44 + 2*x18;
                const T x57 = PBy*x46 - PGy[g]*x49 + x31;
                const T x58 = PBy*x51 - PGy[g]*x54;
                const T x59 = PBz*x51 - PGz[g]*x54 + 2*x31;
                const T x60 = -x7;
                const T x61 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x62 = tpi*(x19 - x22 + x23);
                const T x63 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                const T x64 = tpi*(x32 - x35 + x36);
                const T x65 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                const T x66 = tpi*(x41 - x42 + x43 + x60);
                const T x67 = PBy*x24 - PGy[g]*x63 + x61;
                const T x68 = tpi*(x46 - x47 + x48);
                const T x69 = PBy*x37 - PGy[g]*x65;
                const T x70 = tpi*(x51 - x52 + x53 + x60);
                const T x71 = PBz*x37 - PGz[g]*x65 + x61;
                const T x72 = PBy*x44 - PGy[g]*x67 + 2*x62;
                const T x73 = PBy*x49 - PGy[g]*x69 + x64;
                const T x74 = PBy*x54 - PGy[g]*x71;
                const T x75 = PBz*x54 - PGz[g]*x71 + 2*x64;

                base_integral_0[g] = x14;
                base_integral_1[g] = x27;
                base_integral_2[g] = x40;
                base_integral_3[g] = x45;
                base_integral_4[g] = x50;
                base_integral_5[g] = x55;
                base_integral_6[g] = x56;
                base_integral_7[g] = x57;
                base_integral_8[g] = x58;
                base_integral_9[g] = x59;
                base_integral_10[g] = PBx*x14 - PGx[g]*(PBx*x13 - PGx[g]*(PBx*x11 - PGx[g]*(PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x61) + x3*(x10 + x5 - x9)) + 3*tpi*(x12 + x6 + x60 - x8);
                base_integral_11[g] = PBx*x27 - PGx[g]*(PBx*x26 - PGx[g]*(PBx*x24 - PGx[g]*x63) + x62) + x3*(x20 - x21 + x25);
                base_integral_12[g] = PBx*x40 - PGx[g]*(PBx*x39 - PGx[g]*(PBx*x37 - PGx[g]*x65) + x64) + x3*(x33 - x34 + x38);
                base_integral_13[g] = PBx*x45 - PGx[g]*(PBx*x44 - PGx[g]*x67) + x66;
                base_integral_14[g] = PBx*x50 - PGx[g]*(PBx*x49 - PGx[g]*x69) + x68;
                base_integral_15[g] = PBx*x55 - PGx[g]*(PBx*x54 - PGx[g]*x71) + x70;
                base_integral_16[g] = PBx*x56 - PGx[g]*x72;
                base_integral_17[g] = PBx*x57 - PGx[g]*x73;
                base_integral_18[g] = PBx*x58 - PGx[g]*x74;
                base_integral_19[g] = PBx*x59 - PGx[g]*x75;
                base_integral_20[g] = PBy*x56 - PGy[g]*x72 + 3*x66;
                base_integral_21[g] = PBy*x57 - PGy[g]*x73 + 2*x68;
                base_integral_22[g] = PBy*x58 - PGy[g]*x74 + x70;
                base_integral_23[g] = PBy*x59 - PGy[g]*x75;
                base_integral_24[g] = PBz*x59 - PGz[g]*x75 + 3*x70;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*scaled_boys_2[g];
                const T x1 = PBx*scaled_boys_1[g];
                const T x2 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x3 = 2*tpi;
                const T x4 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x5 = -x0 + x1;
                const T x6 = PBx*x2 - PGx[g]*x5 + x4;
                const T x7 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x8 = PBx*x5;
                const T x9 = PBx*scaled_boys_2[g];
                const T x10 = PGx[g]*scaled_boys_3[g];
                const T x11 = -x10 + x9;
                const T x12 = PGx[g]*x11;
                const T x13 = -x12 + x7 + x8;
                const T x14 = PBx*x6 - PGx[g]*x13 + x3*(x0 - x1 + x2);
                const T x15 = PGy[g]*scaled_boys_2[g];
                const T x16 = PBy*scaled_boys_1[g];
                const T x17 = PBy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x18 = tpi*(x15 - x16 + x17);
                const T x19 = -x15 + x16;
                const T x20 = PBx*x17 - PGx[g]*x19;
                const T x21 = PBx*x19;
                const T x22 = PBy*scaled_boys_2[g];
                const T x23 = PGy[g]*scaled_boys_3[g];
                const T x24 = x22 - x23;
                const T x25 = PGx[g]*x24;
                const T x26 = x21 - x25;
                const T x27 = PBx*x20 - PGx[g]*x26 + x18;
                const T x28 = PGz[g]*scaled_boys_2[g];
                const T x29 = PBz*scaled_boys_1[g];
                const T x30 = PBz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x31 = tpi*(x28 - x29 + x30);
                const T x32 = -x28 + x29;
                const T x33 = PBx*x30 - PGx[g]*x32;
                const T x34 = PBx*x32;
                const T x35 = PBz*scaled_boys_2[g];
                const T x36 = PGz[g]*scaled_boys_3[g];
                const T x37 = x35 - x36;
                const T x38 = PGx[g]*x37;
                const T x39 = x34 - x38;
                const T x40 = PBx*x33 - PGx[g]*x39 + x31;
                const T x41 = PBy*x17 - PGy[g]*x19 + x4;
                const T x42 = PBy*x19;
                const T x43 = PGy[g]*x24;
                const T x44 = x42 - x43 + x7;
                const T x45 = PBx*x41 - PGx[g]*x44;
                const T x46 = PBy*x30 - PGy[g]*x32;
                const T x47 = PBy*x32;
                const T x48 = PGy[g]*x37;
                const T x49 = x47 - x48;
                const T x50 = PBx*x46 - PGx[g]*x49;
                const T x51 = PBz*x30 - PGz[g]*x32 + x4;
                const T x52 = PBz*x32;
                const T x53 = PGz[g]*x37;
                const T x54 = x52 - x53 + x7;
                const T x55 = PBx*x51 - PGx[g]*x54;
                const T x56 = PBy*x41 - PGy[g]*x44 + 2*x18;
                const T x57 = PBy*x46 - PGy[g]*x49 + x31;
                const T x58 = PBy*x51 - PGy[g]*x54;
                const T x59 = PBz*x51 - PGz[g]*x54 + 2*x31;
                const T x60 = -x7;
                const T x61 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x62 = tpi*(x19 - x22 + x23);
                const T x63 = PBy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                const T x64 = tpi*(x32 - x35 + x36);
                const T x65 = PBz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                const T x66 = tpi*(x41 - x42 + x43 + x60);
                const T x67 = PBy*x24 - PGy[g]*x63 + x61;
                const T x68 = tpi*(x46 - x47 + x48);
                const T x69 = PBy*x37 - PGy[g]*x65;
                const T x70 = tpi*(x51 - x52 + x53 + x60);
                const T x71 = PBz*x37 - PGz[g]*x65 + x61;
                const T x72 = PBy*x44 - PGy[g]*x67 + 2*x62;
                const T x73 = PBy*x49 - PGy[g]*x69 + x64;
                const T x74 = PBy*x54 - PGy[g]*x71;
                const T x75 = PBz*x54 - PGz[g]*x71 + 2*x64;

                base_integral_0[g] += x14;
                base_integral_1[g] += x27;
                base_integral_2[g] += x40;
                base_integral_3[g] += x45;
                base_integral_4[g] += x50;
                base_integral_5[g] += x55;
                base_integral_6[g] += x56;
                base_integral_7[g] += x57;
                base_integral_8[g] += x58;
                base_integral_9[g] += x59;
                base_integral_10[g] += PBx*x14 - PGx[g]*(PBx*x13 - PGx[g]*(PBx*x11 - PGx[g]*(PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x61) + x3*(x10 + x5 - x9)) + 3*tpi*(x12 + x6 + x60 - x8);
                base_integral_11[g] += PBx*x27 - PGx[g]*(PBx*x26 - PGx[g]*(PBx*x24 - PGx[g]*x63) + x62) + x3*(x20 - x21 + x25);
                base_integral_12[g] += PBx*x40 - PGx[g]*(PBx*x39 - PGx[g]*(PBx*x37 - PGx[g]*x65) + x64) + x3*(x33 - x34 + x38);
                base_integral_13[g] += PBx*x45 - PGx[g]*(PBx*x44 - PGx[g]*x67) + x66;
                base_integral_14[g] += PBx*x50 - PGx[g]*(PBx*x49 - PGx[g]*x69) + x68;
                base_integral_15[g] += PBx*x55 - PGx[g]*(PBx*x54 - PGx[g]*x71) + x70;
                base_integral_16[g] += PBx*x56 - PGx[g]*x72;
                base_integral_17[g] += PBx*x57 - PGx[g]*x73;
                base_integral_18[g] += PBx*x58 - PGx[g]*x74;
                base_integral_19[g] += PBx*x59 - PGx[g]*x75;
                base_integral_20[g] += PBy*x56 - PGy[g]*x72 + 3*x66;
                base_integral_21[g] += PBy*x57 - PGy[g]*x73 + 2*x68;
                base_integral_22[g] += PBy*x58 - PGy[g]*x74 + x70;
                base_integral_23[g] += PBy*x59 - PGy[g]*x75;
                base_integral_24[g] += PBz*x59 - PGz[g]*x75 + 3*x70;

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
    T* __restrict integral_18 = &integrals[18*gb_size];
    T* __restrict integral_19 = &integrals[19*gb_size];
    T* __restrict integral_20 = &integrals[20*gb_size];
    T* __restrict integral_21 = &integrals[21*gb_size];
    T* __restrict integral_22 = &integrals[22*gb_size];
    T* __restrict integral_23 = &integrals[23*gb_size];
    T* __restrict integral_24 = &integrals[24*gb_size];
    T* __restrict integral_25 = &integrals[25*gb_size];
    T* __restrict integral_26 = &integrals[26*gb_size];
    T* __restrict integral_27 = &integrals[27*gb_size];
    T* __restrict integral_28 = &integrals[28*gb_size];
    T* __restrict integral_29 = &integrals[29*gb_size];

    #pragma omp simd
    for (int g = 0; g < gb_size; ++g) {

        const T x0 = -base_integral_11[g];
        const T x1 = -base_integral_12[g];
        const T x2 = -base_integral_13[g];
        const T x3 = -base_integral_14[g];
        const T x4 = -base_integral_15[g];
        const T x5 = -base_integral_16[g];
        const T x6 = -base_integral_17[g];
        const T x7 = -base_integral_18[g];
        const T x8 = -base_integral_19[g];
        const T x9 = -base_integral_21[g];
        const T x10 = -base_integral_22[g];
        const T x11 = -base_integral_23[g];

        integral_0[g] = -ABx*base_integral_0[g] + base_integral_10[g];
        integral_1[g] = -ABx*base_integral_1[g] - x0;
        integral_2[g] = -ABx*base_integral_2[g] - x1;
        integral_3[g] = -ABx*base_integral_3[g] - x2;
        integral_4[g] = -ABx*base_integral_4[g] - x3;
        integral_5[g] = -ABx*base_integral_5[g] - x4;
        integral_6[g] = -ABx*base_integral_6[g] - x5;
        integral_7[g] = -ABx*base_integral_7[g] - x6;
        integral_8[g] = -ABx*base_integral_8[g] - x7;
        integral_9[g] = -ABx*base_integral_9[g] - x8;
        integral_10[g] = -ABy*base_integral_0[g] - x0;
        integral_11[g] = -ABy*base_integral_1[g] - x2;
        integral_12[g] = -ABy*base_integral_2[g] - x3;
        integral_13[g] = -ABy*base_integral_3[g] - x5;
        integral_14[g] = -ABy*base_integral_4[g] - x6;
        integral_15[g] = -ABy*base_integral_5[g] - x7;
        integral_16[g] = -ABy*base_integral_6[g] + base_integral_20[g];
        integral_17[g] = -ABy*base_integral_7[g] - x9;
        integral_18[g] = -ABy*base_integral_8[g] - x10;
        integral_19[g] = -ABy*base_integral_9[g] - x11;
        integral_20[g] = -ABz*base_integral_0[g] - x1;
        integral_21[g] = -ABz*base_integral_1[g] - x3;
        integral_22[g] = -ABz*base_integral_2[g] - x4;
        integral_23[g] = -ABz*base_integral_3[g] - x6;
        integral_24[g] = -ABz*base_integral_4[g] - x7;
        integral_25[g] = -ABz*base_integral_5[g] - x8;
        integral_26[g] = -ABz*base_integral_6[g] - x9;
        integral_27[g] = -ABz*base_integral_7[g] - x10;
        integral_28[g] = -ABz*base_integral_8[g] - x11;
        integral_29[g] = -ABz*base_integral_9[g] + base_integral_24[g];

    } // for g
} // HRR_P_F

template void HRR_P_F<float>(
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

template void HRR_P_F<double>(
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

