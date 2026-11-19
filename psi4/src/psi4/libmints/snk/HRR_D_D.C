#include "HRR_D_D.h"

template <typename T>
void HRR_D_D(
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
    T* __restrict base_integral_25 = &base_integrals[25*gb_size];
    T* __restrict base_integral_26 = &base_integrals[26*gb_size];
    T* __restrict base_integral_27 = &base_integrals[27*gb_size];
    T* __restrict base_integral_28 = &base_integrals[28*gb_size];
    T* __restrict base_integral_29 = &base_integrals[29*gb_size];
    T* __restrict base_integral_30 = &base_integrals[30*gb_size];

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

                const T x0 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x1 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x2 = PAx*scaled_boys_1[g];
                const T x3 = PGx[g]*scaled_boys_2[g];
                const T x4 = x2 - x3;
                const T x5 = PAx*x1 - PGx[g]*x4 + x0;
                const T x6 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x7 = PAy*scaled_boys_1[g];
                const T x8 = PGy[g]*scaled_boys_2[g];
                const T x9 = x7 - x8;
                const T x10 = PAx*x6 - PGx[g]*x9;
                const T x11 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x12 = PAz*scaled_boys_1[g];
                const T x13 = PGz[g]*scaled_boys_2[g];
                const T x14 = x12 - x13;
                const T x15 = PAx*x11 - PGx[g]*x14;
                const T x16 = PAy*x6 - PGy[g]*x9 + x0;
                const T x17 = PAy*x11 - PGy[g]*x14;
                const T x18 = PAz*x11 - PGz[g]*x14 + x0;
                const T x19 = 2*tpi;
                const T x20 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x21 = PAx*x4;
                const T x22 = PAx*scaled_boys_2[g];
                const T x23 = PGx[g]*scaled_boys_3[g];
                const T x24 = x22 - x23;
                const T x25 = PGx[g]*x24;
                const T x26 = x20 + x21 - x25;
                const T x27 = PAx*x5 - PGx[g]*x26 + x19*(x1 - x2 + x3);
                const T x28 = tpi*(x6 - x7 + x8);
                const T x29 = PAx*x9;
                const T x30 = PAy*scaled_boys_2[g];
                const T x31 = PGy[g]*scaled_boys_3[g];
                const T x32 = x30 - x31;
                const T x33 = PGx[g]*x32;
                const T x34 = x29 - x33;
                const T x35 = PAx*x10 - PGx[g]*x34 + x28;
                const T x36 = tpi*(x11 - x12 + x13);
                const T x37 = PAx*x14;
                const T x38 = PAz*scaled_boys_2[g];
                const T x39 = PGz[g]*scaled_boys_3[g];
                const T x40 = x38 - x39;
                const T x41 = PGx[g]*x40;
                const T x42 = x37 - x41;
                const T x43 = PAx*x15 - PGx[g]*x42 + x36;
                const T x44 = PAy*x9;
                const T x45 = PGy[g]*x32;
                const T x46 = x20 + x44 - x45;
                const T x47 = PAx*x16 - PGx[g]*x46;
                const T x48 = PAy*x14;
                const T x49 = PGy[g]*x40;
                const T x50 = x48 - x49;
                const T x51 = PAx*x17 - PGx[g]*x50;
                const T x52 = PAz*x14;
                const T x53 = PGz[g]*x40;
                const T x54 = x20 + x52 - x53;
                const T x55 = PAx*x18 - PGx[g]*x54;
                const T x56 = PAy*x16 - PGy[g]*x46 + 2*x28;
                const T x57 = PAy*x17 - PGy[g]*x50 + x36;
                const T x58 = PAy*x18 - PGy[g]*x54;
                const T x59 = PAz*x18 - PGz[g]*x54 + 2*x36;
                const T x60 = -x20;
                const T x61 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x62 = tpi*(-x30 + x31 + x9);
                const T x63 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                const T x64 = tpi*(x14 - x38 + x39);
                const T x65 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                const T x66 = tpi*(x16 - x44 + x45 + x60);
                const T x67 = PAy*x32 - PGy[g]*x63 + x61;
                const T x68 = tpi*(x17 - x48 + x49);
                const T x69 = PAy*x40 - PGy[g]*x65;
                const T x70 = tpi*(x18 - x52 + x53 + x60);
                const T x71 = PAz*x40 - PGz[g]*x65 + x61;
                const T x72 = PAy*x46 - PGy[g]*x67 + 2*x62;
                const T x73 = PAy*x50 - PGy[g]*x69 + x64;
                const T x74 = PAy*x54 - PGy[g]*x71;
                const T x75 = PAz*x54 - PGz[g]*x71 + 2*x64;

                base_integral_0[g] = x5;
                base_integral_1[g] = x10;
                base_integral_2[g] = x15;
                base_integral_3[g] = x16;
                base_integral_4[g] = x17;
                base_integral_5[g] = x18;
                base_integral_6[g] = x27;
                base_integral_7[g] = x35;
                base_integral_8[g] = x43;
                base_integral_9[g] = x47;
                base_integral_10[g] = x51;
                base_integral_11[g] = x55;
                base_integral_12[g] = x56;
                base_integral_13[g] = x57;
                base_integral_14[g] = x58;
                base_integral_15[g] = x59;
                base_integral_16[g] = PAx*x27 - PGx[g]*(PAx*x26 - PGx[g]*(PAx*x24 - PGx[g]*(PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x61) + x19*(-x22 + x23 + x4)) + 3*tpi*(-x21 + x25 + x5 + x60);
                base_integral_17[g] = PAx*x35 - PGx[g]*(PAx*x34 - PGx[g]*(PAx*x32 - PGx[g]*x63) + x62) + x19*(x10 - x29 + x33);
                base_integral_18[g] = PAx*x43 - PGx[g]*(PAx*x42 - PGx[g]*(PAx*x40 - PGx[g]*x65) + x64) + x19*(x15 - x37 + x41);
                base_integral_19[g] = PAx*x47 - PGx[g]*(PAx*x46 - PGx[g]*x67) + x66;
                base_integral_20[g] = PAx*x51 - PGx[g]*(PAx*x50 - PGx[g]*x69) + x68;
                base_integral_21[g] = PAx*x55 - PGx[g]*(PAx*x54 - PGx[g]*x71) + x70;
                base_integral_22[g] = PAx*x56 - PGx[g]*x72;
                base_integral_23[g] = PAx*x57 - PGx[g]*x73;
                base_integral_24[g] = PAx*x58 - PGx[g]*x74;
                base_integral_25[g] = PAx*x59 - PGx[g]*x75;
                base_integral_26[g] = PAy*x56 - PGy[g]*x72 + 3*x66;
                base_integral_27[g] = PAy*x57 - PGy[g]*x73 + 2*x68;
                base_integral_28[g] = PAy*x58 - PGy[g]*x74 + x70;
                base_integral_29[g] = PAy*x59 - PGy[g]*x75;
                base_integral_30[g] = PAz*x59 - PGz[g]*x75 + 3*x70;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x1 = PAx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x2 = PAx*scaled_boys_1[g];
                const T x3 = PGx[g]*scaled_boys_2[g];
                const T x4 = x2 - x3;
                const T x5 = PAx*x1 - PGx[g]*x4 + x0;
                const T x6 = PAy*scaled_boys_0[g] - PGy[g]*scaled_boys_1[g];
                const T x7 = PAy*scaled_boys_1[g];
                const T x8 = PGy[g]*scaled_boys_2[g];
                const T x9 = x7 - x8;
                const T x10 = PAx*x6 - PGx[g]*x9;
                const T x11 = PAz*scaled_boys_0[g] - PGz[g]*scaled_boys_1[g];
                const T x12 = PAz*scaled_boys_1[g];
                const T x13 = PGz[g]*scaled_boys_2[g];
                const T x14 = x12 - x13;
                const T x15 = PAx*x11 - PGx[g]*x14;
                const T x16 = PAy*x6 - PGy[g]*x9 + x0;
                const T x17 = PAy*x11 - PGy[g]*x14;
                const T x18 = PAz*x11 - PGz[g]*x14 + x0;
                const T x19 = 2*tpi;
                const T x20 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x21 = PAx*x4;
                const T x22 = PAx*scaled_boys_2[g];
                const T x23 = PGx[g]*scaled_boys_3[g];
                const T x24 = x22 - x23;
                const T x25 = PGx[g]*x24;
                const T x26 = x20 + x21 - x25;
                const T x27 = PAx*x5 - PGx[g]*x26 + x19*(x1 - x2 + x3);
                const T x28 = tpi*(x6 - x7 + x8);
                const T x29 = PAx*x9;
                const T x30 = PAy*scaled_boys_2[g];
                const T x31 = PGy[g]*scaled_boys_3[g];
                const T x32 = x30 - x31;
                const T x33 = PGx[g]*x32;
                const T x34 = x29 - x33;
                const T x35 = PAx*x10 - PGx[g]*x34 + x28;
                const T x36 = tpi*(x11 - x12 + x13);
                const T x37 = PAx*x14;
                const T x38 = PAz*scaled_boys_2[g];
                const T x39 = PGz[g]*scaled_boys_3[g];
                const T x40 = x38 - x39;
                const T x41 = PGx[g]*x40;
                const T x42 = x37 - x41;
                const T x43 = PAx*x15 - PGx[g]*x42 + x36;
                const T x44 = PAy*x9;
                const T x45 = PGy[g]*x32;
                const T x46 = x20 + x44 - x45;
                const T x47 = PAx*x16 - PGx[g]*x46;
                const T x48 = PAy*x14;
                const T x49 = PGy[g]*x40;
                const T x50 = x48 - x49;
                const T x51 = PAx*x17 - PGx[g]*x50;
                const T x52 = PAz*x14;
                const T x53 = PGz[g]*x40;
                const T x54 = x20 + x52 - x53;
                const T x55 = PAx*x18 - PGx[g]*x54;
                const T x56 = PAy*x16 - PGy[g]*x46 + 2*x28;
                const T x57 = PAy*x17 - PGy[g]*x50 + x36;
                const T x58 = PAy*x18 - PGy[g]*x54;
                const T x59 = PAz*x18 - PGz[g]*x54 + 2*x36;
                const T x60 = -x20;
                const T x61 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x62 = tpi*(-x30 + x31 + x9);
                const T x63 = PAy*scaled_boys_3[g] - PGy[g]*scaled_boys_4[g];
                const T x64 = tpi*(x14 - x38 + x39);
                const T x65 = PAz*scaled_boys_3[g] - PGz[g]*scaled_boys_4[g];
                const T x66 = tpi*(x16 - x44 + x45 + x60);
                const T x67 = PAy*x32 - PGy[g]*x63 + x61;
                const T x68 = tpi*(x17 - x48 + x49);
                const T x69 = PAy*x40 - PGy[g]*x65;
                const T x70 = tpi*(x18 - x52 + x53 + x60);
                const T x71 = PAz*x40 - PGz[g]*x65 + x61;
                const T x72 = PAy*x46 - PGy[g]*x67 + 2*x62;
                const T x73 = PAy*x50 - PGy[g]*x69 + x64;
                const T x74 = PAy*x54 - PGy[g]*x71;
                const T x75 = PAz*x54 - PGz[g]*x71 + 2*x64;

                base_integral_0[g] += x5;
                base_integral_1[g] += x10;
                base_integral_2[g] += x15;
                base_integral_3[g] += x16;
                base_integral_4[g] += x17;
                base_integral_5[g] += x18;
                base_integral_6[g] += x27;
                base_integral_7[g] += x35;
                base_integral_8[g] += x43;
                base_integral_9[g] += x47;
                base_integral_10[g] += x51;
                base_integral_11[g] += x55;
                base_integral_12[g] += x56;
                base_integral_13[g] += x57;
                base_integral_14[g] += x58;
                base_integral_15[g] += x59;
                base_integral_16[g] += PAx*x27 - PGx[g]*(PAx*x26 - PGx[g]*(PAx*x24 - PGx[g]*(PAx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x61) + x19*(-x22 + x23 + x4)) + 3*tpi*(-x21 + x25 + x5 + x60);
                base_integral_17[g] += PAx*x35 - PGx[g]*(PAx*x34 - PGx[g]*(PAx*x32 - PGx[g]*x63) + x62) + x19*(x10 - x29 + x33);
                base_integral_18[g] += PAx*x43 - PGx[g]*(PAx*x42 - PGx[g]*(PAx*x40 - PGx[g]*x65) + x64) + x19*(x15 - x37 + x41);
                base_integral_19[g] += PAx*x47 - PGx[g]*(PAx*x46 - PGx[g]*x67) + x66;
                base_integral_20[g] += PAx*x51 - PGx[g]*(PAx*x50 - PGx[g]*x69) + x68;
                base_integral_21[g] += PAx*x55 - PGx[g]*(PAx*x54 - PGx[g]*x71) + x70;
                base_integral_22[g] += PAx*x56 - PGx[g]*x72;
                base_integral_23[g] += PAx*x57 - PGx[g]*x73;
                base_integral_24[g] += PAx*x58 - PGx[g]*x74;
                base_integral_25[g] += PAx*x59 - PGx[g]*x75;
                base_integral_26[g] += PAy*x56 - PGy[g]*x72 + 3*x66;
                base_integral_27[g] += PAy*x57 - PGy[g]*x73 + 2*x68;
                base_integral_28[g] += PAy*x58 - PGy[g]*x74 + x70;
                base_integral_29[g] += PAy*x59 - PGy[g]*x75;
                base_integral_30[g] += PAz*x59 - PGz[g]*x75 + 3*x70;

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
    T* __restrict integral_30 = &integrals[30*gb_size];
    T* __restrict integral_31 = &integrals[31*gb_size];
    T* __restrict integral_32 = &integrals[32*gb_size];
    T* __restrict integral_33 = &integrals[33*gb_size];
    T* __restrict integral_34 = &integrals[34*gb_size];
    T* __restrict integral_35 = &integrals[35*gb_size];

    #pragma omp simd
    for (int g = 0; g < gb_size; ++g) {

        const T x0 = ABy*base_integral_0[g] + base_integral_7[g];
        const T x1 = ABz*base_integral_0[g] + base_integral_8[g];
        const T x2 = ABy*base_integral_7[g] + base_integral_19[g];
        const T x3 = ABz*base_integral_7[g] + base_integral_20[g];
        const T x4 = ABz*base_integral_8[g] + base_integral_21[g];
        const T x5 = ABy*base_integral_1[g] + base_integral_9[g];
        const T x6 = ABz*base_integral_1[g] + base_integral_10[g];
        const T x7 = ABy*base_integral_9[g] + base_integral_22[g];
        const T x8 = ABz*base_integral_9[g] + base_integral_23[g];
        const T x9 = ABz*base_integral_10[g] + base_integral_24[g];
        const T x10 = ABy*base_integral_2[g] + base_integral_10[g];
        const T x11 = ABz*base_integral_2[g] + base_integral_11[g];
        const T x12 = ABy*base_integral_10[g] + base_integral_23[g];
        const T x13 = ABz*base_integral_11[g] + base_integral_25[g];
        const T x14 = ABy*base_integral_3[g] + base_integral_12[g];
        const T x15 = ABz*base_integral_3[g] + base_integral_13[g];
        const T x16 = ABz*base_integral_13[g] + base_integral_28[g];
        const T x17 = ABy*base_integral_4[g] + base_integral_13[g];
        const T x18 = ABz*base_integral_4[g] + base_integral_14[g];
        const T x19 = ABz*base_integral_14[g] + base_integral_29[g];
        const T x20 = ABy*base_integral_5[g] + base_integral_14[g];
        const T x21 = ABz*base_integral_5[g] + base_integral_15[g];

        integral_0[g] = ABx*base_integral_6[g] + ABx*(ABx*base_integral_0[g] + base_integral_6[g]) + base_integral_16[g];
        integral_1[g] = ABx*x0 + ABy*base_integral_6[g] + base_integral_17[g];
        integral_2[g] = ABx*x1 + ABz*base_integral_6[g] + base_integral_18[g];
        integral_3[g] = ABy*x0 + x2;
        integral_4[g] = ABy*x1 + x3;
        integral_5[g] = ABz*x1 + x4;
        integral_6[g] = ABx*base_integral_7[g] + ABx*(ABx*base_integral_1[g] + base_integral_7[g]) + base_integral_17[g];
        integral_7[g] = ABx*x5 + x2;
        integral_8[g] = ABx*x6 + x3;
        integral_9[g] = ABy*x5 + x7;
        integral_10[g] = ABy*x6 + x8;
        integral_11[g] = ABz*x6 + x9;
        integral_12[g] = ABx*base_integral_8[g] + ABx*(ABx*base_integral_2[g] + base_integral_8[g]) + base_integral_18[g];
        integral_13[g] = ABx*x10 + ABy*base_integral_8[g] + base_integral_20[g];
        integral_14[g] = ABx*x11 + x4;
        integral_15[g] = ABy*x10 + x12;
        integral_16[g] = ABy*x11 + x9;
        integral_17[g] = ABz*x11 + x13;
        integral_18[g] = ABx*base_integral_9[g] + ABx*(ABx*base_integral_3[g] + base_integral_9[g]) + base_integral_19[g];
        integral_19[g] = ABx*x14 + x7;
        integral_20[g] = ABx*x15 + x8;
        integral_21[g] = ABy*base_integral_12[g] + ABy*x14 + base_integral_26[g];
        integral_22[g] = ABy*x15 + ABz*base_integral_12[g] + base_integral_27[g];
        integral_23[g] = ABz*x15 + x16;
        integral_24[g] = ABx*base_integral_10[g] + ABx*(ABx*base_integral_4[g] + base_integral_10[g]) + base_integral_20[g];
        integral_25[g] = ABx*x17 + x12;
        integral_26[g] = ABx*x18 + x9;
        integral_27[g] = ABy*base_integral_13[g] + ABy*x17 + base_integral_27[g];
        integral_28[g] = ABy*x18 + x16;
        integral_29[g] = ABz*x18 + x19;
        integral_30[g] = ABx*base_integral_11[g] + ABx*(ABx*base_integral_5[g] + base_integral_11[g]) + base_integral_21[g];
        integral_31[g] = ABx*x20 + ABy*base_integral_11[g] + base_integral_24[g];
        integral_32[g] = ABx*x21 + x13;
        integral_33[g] = ABy*base_integral_14[g] + ABy*x20 + base_integral_28[g];
        integral_34[g] = ABy*x21 + x19;
        integral_35[g] = ABz*base_integral_15[g] + ABz*x21 + base_integral_30[g];

    } // for g
} // HRR_D_D

template void HRR_D_D<float>(
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

template void HRR_D_D<double>(
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

