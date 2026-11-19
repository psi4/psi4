#include "HRR_F_F_SAMECENTER.h"

template <typename T>
void HRR_F_F_SAMECENTER(
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

    for (int prim = 0; prim < nprim; ++prim) {

        const T tpi = 0.5 / v_p[prim];
        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*7*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*7*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*7*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*7*gb_size+3*gb_size];
        const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*7*gb_size+4*gb_size];
        const T* __restrict scaled_boys_5 = &v_scaled_boys[prim*7*gb_size+5*gb_size];
        const T* __restrict scaled_boys_6 = &v_scaled_boys[prim*7*gb_size+6*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_6[g]*x0;
                const T x2 = scaled_boys_4[g] - scaled_boys_5[g];
                const T x3 = tpi*x2;
                const T x4 = -tpi*x2;
                const T x5 = -2*x4;
                const T x6 = x1 + x3 + x5;
                const T x7 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x8 = tpi*x7;
                const T x9 = scaled_boys_4[g]*x0;
                const T x10 = x8 + x9;
                const T x11 = scaled_boys_5[g]*x0;
                const T x12 = -x11;
                const T x13 = scaled_boys_3[g] - scaled_boys_4[g];
                const T x14 = tpi*x13;
                const T x15 = -x14;
                const T x16 = x12 + x15;
                const T x17 = x10 + x16;
                const T x18 = -tpi*x13;
                const T x19 = 2*x18;
                const T x20 = -tpi*x7;
                const T x21 = -2*x20;
                const T x22 = x10 + x21;
                const T x23 = 4*tpi;
                const T x24 = -x19;
                const T x25 = scaled_boys_3[g]*x0;
                const T x26 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x27 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]) - x26;
                const T x28 = 3*tpi;
                const T x29 = x26 - x8;
                const T x30 = -x4;
                const T x31 = 2*x3;
                const T x32 = x18 - x20;
                const T x33 = x12 + x32 + x9;
                const T x34 = -2*tpi*x13 + 2*x8;
                const T x35 = PGx[g]*(x0*(x1 + x30 + x31) - x23*(-x33 - x34) + x28*x33);
                const T x36 = PGy[g]*PGy[g];
                const T x37 = scaled_boys_6[g]*x36;
                const T x38 = x3 + x37;
                const T x39 = scaled_boys_4[g]*x36;
                const T x40 = x39 + x8;
                const T x41 = scaled_boys_5[g]*x36;
                const T x42 = -x41;
                const T x43 = x15 + x42;
                const T x44 = x40 + x43;
                const T x45 = tpi*x44;
                const T x46 = 2*tpi;
                const T x47 = x14 + x41;
                const T x48 = scaled_boys_3[g]*x36;
                const T x49 = tpi*(scaled_boys_2[g]*x36 + x27 - x48);
                const T x50 = tpi*(x29 - x39 + x48);
                const T x51 = PGy[g]*PGz[g];
                const T x52 = PGz[g]*PGz[g];
                const T x53 = scaled_boys_6[g]*x52 + x3;
                const T x54 = x0*x53;
                const T x55 = scaled_boys_4[g]*x52;
                const T x56 = x55 + x8;
                const T x57 = scaled_boys_5[g]*x52;
                const T x58 = x15 - x57;
                const T x59 = x56 + x58;
                const T x60 = tpi*x59;
                const T x61 = -tpi*x59;
                const T x62 = -2*x61;
                const T x63 = x14 + x57;
                const T x64 = scaled_boys_3[g]*x52;
                const T x65 = tpi*(scaled_boys_2[g]*x52 + x27 - x64);
                const T x66 = tpi*(x29 - x55 + x64);
                const T x67 = x65 - x66;
                const T x68 = x38 + x5;
                const T x69 = x21 + x40;
                const T x70 = x19 + x43 + x69;
                const T x71 = tpi*x70;
                const T x72 = -x70;
                const T x73 = PGx[g]*PGy[g];
                const T x74 = x30 + x37;
                const T x75 = x32 + x39 + x42;
                const T x76 = -x75;
                const T x77 = PGx[g]*PGz[g];
                const T x78 = x5 + x53;
                const T x79 = x0*x78;
                const T x80 = x21 + x56;
                const T x81 = x19 + x58 + x80;
                const T x82 = tpi*x81;
                const T x83 = -tpi*x81;
                const T x84 = -2*x83;
                const T x85 = -x36*x68;
                const T x86 = 3*x45;
                const T x87 = x85 - x86;
                const T x88 = tpi*(-x36*x69 - x36*(-x24 - x47) - 3*x49 + 3*x50);
                const T x89 = x31 + x74;
                const T x90 = -x89;
                const T x91 = tpi*(-x34 - x75);
                const T x92 = x36*x53;
                const T x93 = x60 + x92;
                const T x94 = tpi*(x36*x56 - x36*x63 + x67);
                const T x95 = -x52*x78;
                const T x96 = 3*x60;
                const T x97 = x95 - x96;
                const T x98 = tpi*(-x52*x80 - x52*(-x24 - x63) - 3*x65 + 3*x66);
                const T x99 = x36*x78 + x82;
                const T x100 = 4*x82 - x95 + x96;

                base_integral_0[g] = -5*tpi*(-x0*x22 - x0*(-x11 - x14 - x24) + x28*(x25 + x29 - x9) - x28*(scaled_boys_2[g]*x0 - x25 + x27)) + x0*(3*tpi*x17 + x0*x6 - x23*(-x16 - x19 - x22));
                base_integral_1[g] = PGy[g]*x35;
                base_integral_2[g] = PGz[g]*x35;
                base_integral_3[g] = x0*(x0*x38 + x44*x46 + x45) + x28*(x0*x40 - x0*x47 + x49 - x50);
                base_integral_4[g] = x51*(x0*x6 + x17*x28);
                base_integral_5[g] = x0*(x54 + x60 + x62) + x28*(x0*x56 - x0*x63 + x67);
                base_integral_6[g] = x73*(x0*x68 - x46*x72 + x71);
                base_integral_7[g] = x77*(tpi*x75 + x0*x74 - x46*x76);
                base_integral_8[g] = x73*(x54 + 2*x60 - x61);
                base_integral_9[g] = x77*(x79 + x82 + x84);
                base_integral_10[g] = -x0*x87 - x88;
                base_integral_11[g] = -x51*(x0*x90 + x91);
                base_integral_12[g] = x0*x93 + x94;
                base_integral_13[g] = x51*(x79 - x83);
                base_integral_14[g] = -x0*x97 - x98;
                base_integral_15[g] = x73*(4*x71 - x85 + x86);
                base_integral_16[g] = x77*(3*tpi*x75 - x36*x90);
                base_integral_17[g] = x73*(x92 + x96);
                base_integral_18[g] = x77*x99;
                base_integral_19[g] = x73*(x52*x78 + x96);
                base_integral_20[g] = x100*x77;
                base_integral_21[g] = x36*(-x23*x72 - x87) - 5*x88;
                base_integral_22[g] = x51*(-x28*x76 + x36*x89 - 4*x91);
                base_integral_23[g] = x36*(x62 + x93) + 3*x94;
                base_integral_24[g] = x51*(x84 + x99);
                base_integral_25[g] = -x36*x97 - x98;
                base_integral_26[g] = x100*x51;
                base_integral_27[g] = x52*(-4*x83 - x97) - 5*x98;

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_6[g]*x0;
                const T x2 = scaled_boys_4[g] - scaled_boys_5[g];
                const T x3 = tpi*x2;
                const T x4 = -tpi*x2;
                const T x5 = -2*x4;
                const T x6 = x1 + x3 + x5;
                const T x7 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x8 = tpi*x7;
                const T x9 = scaled_boys_4[g]*x0;
                const T x10 = x8 + x9;
                const T x11 = scaled_boys_5[g]*x0;
                const T x12 = -x11;
                const T x13 = scaled_boys_3[g] - scaled_boys_4[g];
                const T x14 = tpi*x13;
                const T x15 = -x14;
                const T x16 = x12 + x15;
                const T x17 = x10 + x16;
                const T x18 = -tpi*x13;
                const T x19 = 2*x18;
                const T x20 = -tpi*x7;
                const T x21 = -2*x20;
                const T x22 = x10 + x21;
                const T x23 = 4*tpi;
                const T x24 = -x19;
                const T x25 = scaled_boys_3[g]*x0;
                const T x26 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x27 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]) - x26;
                const T x28 = 3*tpi;
                const T x29 = x26 - x8;
                const T x30 = -x4;
                const T x31 = 2*x3;
                const T x32 = x18 - x20;
                const T x33 = x12 + x32 + x9;
                const T x34 = -2*tpi*x13 + 2*x8;
                const T x35 = PGx[g]*(x0*(x1 + x30 + x31) - x23*(-x33 - x34) + x28*x33);
                const T x36 = PGy[g]*PGy[g];
                const T x37 = scaled_boys_6[g]*x36;
                const T x38 = x3 + x37;
                const T x39 = scaled_boys_4[g]*x36;
                const T x40 = x39 + x8;
                const T x41 = scaled_boys_5[g]*x36;
                const T x42 = -x41;
                const T x43 = x15 + x42;
                const T x44 = x40 + x43;
                const T x45 = tpi*x44;
                const T x46 = 2*tpi;
                const T x47 = x14 + x41;
                const T x48 = scaled_boys_3[g]*x36;
                const T x49 = tpi*(scaled_boys_2[g]*x36 + x27 - x48);
                const T x50 = tpi*(x29 - x39 + x48);
                const T x51 = PGy[g]*PGz[g];
                const T x52 = PGz[g]*PGz[g];
                const T x53 = scaled_boys_6[g]*x52 + x3;
                const T x54 = x0*x53;
                const T x55 = scaled_boys_4[g]*x52;
                const T x56 = x55 + x8;
                const T x57 = scaled_boys_5[g]*x52;
                const T x58 = x15 - x57;
                const T x59 = x56 + x58;
                const T x60 = tpi*x59;
                const T x61 = -tpi*x59;
                const T x62 = -2*x61;
                const T x63 = x14 + x57;
                const T x64 = scaled_boys_3[g]*x52;
                const T x65 = tpi*(scaled_boys_2[g]*x52 + x27 - x64);
                const T x66 = tpi*(x29 - x55 + x64);
                const T x67 = x65 - x66;
                const T x68 = x38 + x5;
                const T x69 = x21 + x40;
                const T x70 = x19 + x43 + x69;
                const T x71 = tpi*x70;
                const T x72 = -x70;
                const T x73 = PGx[g]*PGy[g];
                const T x74 = x30 + x37;
                const T x75 = x32 + x39 + x42;
                const T x76 = -x75;
                const T x77 = PGx[g]*PGz[g];
                const T x78 = x5 + x53;
                const T x79 = x0*x78;
                const T x80 = x21 + x56;
                const T x81 = x19 + x58 + x80;
                const T x82 = tpi*x81;
                const T x83 = -tpi*x81;
                const T x84 = -2*x83;
                const T x85 = -x36*x68;
                const T x86 = 3*x45;
                const T x87 = x85 - x86;
                const T x88 = tpi*(-x36*x69 - x36*(-x24 - x47) - 3*x49 + 3*x50);
                const T x89 = x31 + x74;
                const T x90 = -x89;
                const T x91 = tpi*(-x34 - x75);
                const T x92 = x36*x53;
                const T x93 = x60 + x92;
                const T x94 = tpi*(x36*x56 - x36*x63 + x67);
                const T x95 = -x52*x78;
                const T x96 = 3*x60;
                const T x97 = x95 - x96;
                const T x98 = tpi*(-x52*x80 - x52*(-x24 - x63) - 3*x65 + 3*x66);
                const T x99 = x36*x78 + x82;
                const T x100 = 4*x82 - x95 + x96;

                base_integral_0[g] += -5*tpi*(-x0*x22 - x0*(-x11 - x14 - x24) + x28*(x25 + x29 - x9) - x28*(scaled_boys_2[g]*x0 - x25 + x27)) + x0*(3*tpi*x17 + x0*x6 - x23*(-x16 - x19 - x22));
                base_integral_1[g] += PGy[g]*x35;
                base_integral_2[g] += PGz[g]*x35;
                base_integral_3[g] += x0*(x0*x38 + x44*x46 + x45) + x28*(x0*x40 - x0*x47 + x49 - x50);
                base_integral_4[g] += x51*(x0*x6 + x17*x28);
                base_integral_5[g] += x0*(x54 + x60 + x62) + x28*(x0*x56 - x0*x63 + x67);
                base_integral_6[g] += x73*(x0*x68 - x46*x72 + x71);
                base_integral_7[g] += x77*(tpi*x75 + x0*x74 - x46*x76);
                base_integral_8[g] += x73*(x54 + 2*x60 - x61);
                base_integral_9[g] += x77*(x79 + x82 + x84);
                base_integral_10[g] += -x0*x87 - x88;
                base_integral_11[g] += -x51*(x0*x90 + x91);
                base_integral_12[g] += x0*x93 + x94;
                base_integral_13[g] += x51*(x79 - x83);
                base_integral_14[g] += -x0*x97 - x98;
                base_integral_15[g] += x73*(4*x71 - x85 + x86);
                base_integral_16[g] += x77*(3*tpi*x75 - x36*x90);
                base_integral_17[g] += x73*(x92 + x96);
                base_integral_18[g] += x77*x99;
                base_integral_19[g] += x73*(x52*x78 + x96);
                base_integral_20[g] += x100*x77;
                base_integral_21[g] += x36*(-x23*x72 - x87) - 5*x88;
                base_integral_22[g] += x51*(-x28*x76 + x36*x89 - 4*x91);
                base_integral_23[g] += x36*(x62 + x93) + 3*x94;
                base_integral_24[g] += x51*(x84 + x99);
                base_integral_25[g] += -x36*x97 - x98;
                base_integral_26[g] += x100*x51;
                base_integral_27[g] += x52*(-4*x83 - x97) - 5*x98;

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
    T* __restrict integral_36 = &integrals[36*gb_size];
    T* __restrict integral_37 = &integrals[37*gb_size];
    T* __restrict integral_38 = &integrals[38*gb_size];
    T* __restrict integral_39 = &integrals[39*gb_size];
    T* __restrict integral_40 = &integrals[40*gb_size];
    T* __restrict integral_41 = &integrals[41*gb_size];
    T* __restrict integral_42 = &integrals[42*gb_size];
    T* __restrict integral_43 = &integrals[43*gb_size];
    T* __restrict integral_44 = &integrals[44*gb_size];
    T* __restrict integral_45 = &integrals[45*gb_size];
    T* __restrict integral_46 = &integrals[46*gb_size];
    T* __restrict integral_47 = &integrals[47*gb_size];
    T* __restrict integral_48 = &integrals[48*gb_size];
    T* __restrict integral_49 = &integrals[49*gb_size];
    T* __restrict integral_50 = &integrals[50*gb_size];
    T* __restrict integral_51 = &integrals[51*gb_size];
    T* __restrict integral_52 = &integrals[52*gb_size];
    T* __restrict integral_53 = &integrals[53*gb_size];
    T* __restrict integral_54 = &integrals[54*gb_size];
    T* __restrict integral_55 = &integrals[55*gb_size];
    T* __restrict integral_56 = &integrals[56*gb_size];
    T* __restrict integral_57 = &integrals[57*gb_size];
    T* __restrict integral_58 = &integrals[58*gb_size];
    T* __restrict integral_59 = &integrals[59*gb_size];
    T* __restrict integral_60 = &integrals[60*gb_size];
    T* __restrict integral_61 = &integrals[61*gb_size];
    T* __restrict integral_62 = &integrals[62*gb_size];
    T* __restrict integral_63 = &integrals[63*gb_size];
    T* __restrict integral_64 = &integrals[64*gb_size];
    T* __restrict integral_65 = &integrals[65*gb_size];
    T* __restrict integral_66 = &integrals[66*gb_size];
    T* __restrict integral_67 = &integrals[67*gb_size];
    T* __restrict integral_68 = &integrals[68*gb_size];
    T* __restrict integral_69 = &integrals[69*gb_size];
    T* __restrict integral_70 = &integrals[70*gb_size];
    T* __restrict integral_71 = &integrals[71*gb_size];
    T* __restrict integral_72 = &integrals[72*gb_size];
    T* __restrict integral_73 = &integrals[73*gb_size];
    T* __restrict integral_74 = &integrals[74*gb_size];
    T* __restrict integral_75 = &integrals[75*gb_size];
    T* __restrict integral_76 = &integrals[76*gb_size];
    T* __restrict integral_77 = &integrals[77*gb_size];
    T* __restrict integral_78 = &integrals[78*gb_size];
    T* __restrict integral_79 = &integrals[79*gb_size];
    T* __restrict integral_80 = &integrals[80*gb_size];
    T* __restrict integral_81 = &integrals[81*gb_size];
    T* __restrict integral_82 = &integrals[82*gb_size];
    T* __restrict integral_83 = &integrals[83*gb_size];
    T* __restrict integral_84 = &integrals[84*gb_size];
    T* __restrict integral_85 = &integrals[85*gb_size];
    T* __restrict integral_86 = &integrals[86*gb_size];
    T* __restrict integral_87 = &integrals[87*gb_size];
    T* __restrict integral_88 = &integrals[88*gb_size];
    T* __restrict integral_89 = &integrals[89*gb_size];
    T* __restrict integral_90 = &integrals[90*gb_size];
    T* __restrict integral_91 = &integrals[91*gb_size];
    T* __restrict integral_92 = &integrals[92*gb_size];
    T* __restrict integral_93 = &integrals[93*gb_size];
    T* __restrict integral_94 = &integrals[94*gb_size];
    T* __restrict integral_95 = &integrals[95*gb_size];
    T* __restrict integral_96 = &integrals[96*gb_size];
    T* __restrict integral_97 = &integrals[97*gb_size];
    T* __restrict integral_98 = &integrals[98*gb_size];
    T* __restrict integral_99 = &integrals[99*gb_size];

    #pragma omp simd
    for (int g = 0; g < gb_size; ++g) {

        integral_0[g] = base_integral_0[g];
        integral_1[g] = base_integral_1[g];
        integral_2[g] = base_integral_2[g];
        integral_3[g] = base_integral_3[g];
        integral_4[g] = base_integral_4[g];
        integral_5[g] = base_integral_5[g];
        integral_6[g] = base_integral_6[g];
        integral_7[g] = base_integral_7[g];
        integral_8[g] = base_integral_8[g];
        integral_9[g] = base_integral_9[g];
        integral_10[g] = base_integral_1[g];
        integral_11[g] = base_integral_3[g];
        integral_12[g] = base_integral_4[g];
        integral_13[g] = base_integral_6[g];
        integral_14[g] = base_integral_7[g];
        integral_15[g] = base_integral_8[g];
        integral_16[g] = base_integral_10[g];
        integral_17[g] = base_integral_11[g];
        integral_18[g] = base_integral_12[g];
        integral_19[g] = base_integral_13[g];
        integral_20[g] = base_integral_2[g];
        integral_21[g] = base_integral_4[g];
        integral_22[g] = base_integral_5[g];
        integral_23[g] = base_integral_7[g];
        integral_24[g] = base_integral_8[g];
        integral_25[g] = base_integral_9[g];
        integral_26[g] = base_integral_11[g];
        integral_27[g] = base_integral_12[g];
        integral_28[g] = base_integral_13[g];
        integral_29[g] = base_integral_14[g];
        integral_30[g] = base_integral_3[g];
        integral_31[g] = base_integral_6[g];
        integral_32[g] = base_integral_7[g];
        integral_33[g] = base_integral_10[g];
        integral_34[g] = base_integral_11[g];
        integral_35[g] = base_integral_12[g];
        integral_36[g] = base_integral_15[g];
        integral_37[g] = base_integral_16[g];
        integral_38[g] = base_integral_17[g];
        integral_39[g] = base_integral_18[g];
        integral_40[g] = base_integral_4[g];
        integral_41[g] = base_integral_7[g];
        integral_42[g] = base_integral_8[g];
        integral_43[g] = base_integral_11[g];
        integral_44[g] = base_integral_12[g];
        integral_45[g] = base_integral_13[g];
        integral_46[g] = base_integral_16[g];
        integral_47[g] = base_integral_17[g];
        integral_48[g] = base_integral_18[g];
        integral_49[g] = base_integral_19[g];
        integral_50[g] = base_integral_5[g];
        integral_51[g] = base_integral_8[g];
        integral_52[g] = base_integral_9[g];
        integral_53[g] = base_integral_12[g];
        integral_54[g] = base_integral_13[g];
        integral_55[g] = base_integral_14[g];
        integral_56[g] = base_integral_17[g];
        integral_57[g] = base_integral_18[g];
        integral_58[g] = base_integral_19[g];
        integral_59[g] = base_integral_20[g];
        integral_60[g] = base_integral_6[g];
        integral_61[g] = base_integral_10[g];
        integral_62[g] = base_integral_11[g];
        integral_63[g] = base_integral_15[g];
        integral_64[g] = base_integral_16[g];
        integral_65[g] = base_integral_17[g];
        integral_66[g] = base_integral_21[g];
        integral_67[g] = base_integral_22[g];
        integral_68[g] = base_integral_23[g];
        integral_69[g] = base_integral_24[g];
        integral_70[g] = base_integral_7[g];
        integral_71[g] = base_integral_11[g];
        integral_72[g] = base_integral_12[g];
        integral_73[g] = base_integral_16[g];
        integral_74[g] = base_integral_17[g];
        integral_75[g] = base_integral_18[g];
        integral_76[g] = base_integral_22[g];
        integral_77[g] = base_integral_23[g];
        integral_78[g] = base_integral_24[g];
        integral_79[g] = base_integral_25[g];
        integral_80[g] = base_integral_8[g];
        integral_81[g] = base_integral_12[g];
        integral_82[g] = base_integral_13[g];
        integral_83[g] = base_integral_17[g];
        integral_84[g] = base_integral_18[g];
        integral_85[g] = base_integral_19[g];
        integral_86[g] = base_integral_23[g];
        integral_87[g] = base_integral_24[g];
        integral_88[g] = base_integral_25[g];
        integral_89[g] = base_integral_26[g];
        integral_90[g] = base_integral_9[g];
        integral_91[g] = base_integral_13[g];
        integral_92[g] = base_integral_14[g];
        integral_93[g] = base_integral_18[g];
        integral_94[g] = base_integral_19[g];
        integral_95[g] = base_integral_20[g];
        integral_96[g] = base_integral_24[g];
        integral_97[g] = base_integral_25[g];
        integral_98[g] = base_integral_26[g];
        integral_99[g] = base_integral_27[g];

    } // for g
} // HRR_F_F_SAMECENTER

template void HRR_F_F_SAMECENTER<float>(
        float* integrals, 
        float* base_integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void HRR_F_F_SAMECENTER<double>(
        double* integrals, 
        double* base_integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

