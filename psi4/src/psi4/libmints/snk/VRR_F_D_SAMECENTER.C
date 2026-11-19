#include "VRR_F_D_SAMECENTER.h"

template <typename T>
void VRR_F_D_SAMECENTER(
        T* integrals,              //!< overwritten, cartesian integrals for this grid batch
        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)
        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively
        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)
        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)
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

    for (int prim = 0; prim < nprim; ++prim) {

        const T tpi = 0.5 / v_p[prim];
        const T* __restrict PGx = &v_PGx[prim*gb_size];
        const T* __restrict PGy = &v_PGy[prim*gb_size];
        const T* __restrict PGz = &v_PGz[prim*gb_size];

        const T* __restrict scaled_boys_0 = &v_scaled_boys[prim*6*gb_size+0*gb_size];
        const T* __restrict scaled_boys_1 = &v_scaled_boys[prim*6*gb_size+1*gb_size];
        const T* __restrict scaled_boys_2 = &v_scaled_boys[prim*6*gb_size+2*gb_size];
        const T* __restrict scaled_boys_3 = &v_scaled_boys[prim*6*gb_size+3*gb_size];
        const T* __restrict scaled_boys_4 = &v_scaled_boys[prim*6*gb_size+4*gb_size];
        const T* __restrict scaled_boys_5 = &v_scaled_boys[prim*6*gb_size+5*gb_size];

        // overwrite integrals for first primitive, accumulate otherwise
        if (prim==0) {
            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_5[g]*x0;
                const T x2 = scaled_boys_3[g] - scaled_boys_4[g];
                const T x3 = tpi*x2;
                const T x4 = -tpi*x2;
                const T x5 = -2*x4;
                const T x6 = scaled_boys_3[g]*x0 - scaled_boys_4[g]*x0;
                const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x8 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x9 = tpi*x7 - tpi*x8;
                const T x10 = x6 + x9;
                const T x11 = 3*tpi;
                const T x12 = -tpi*x7;
                const T x13 = -tpi*x8;
                const T x14 = -2*x12 + 2*x13;
                const T x15 = -x4;
                const T x16 = 2*x3;
                const T x17 = -x12 + x13;
                const T x18 = x0*(-x1 - x15 - x16) - x11*(x17 + x6);
                const T x19 = PGy[g]*x18;
                const T x20 = PGz[g]*x18;
                const T x21 = PGy[g]*PGy[g];
                const T x22 = scaled_boys_5[g]*x21;
                const T x23 = x22 + x3;
                const T x24 = scaled_boys_3[g]*x21 - scaled_boys_4[g]*x21;
                const T x25 = x24 + x9;
                const T x26 = x11*x25;
                const T x27 = -PGx[g]*(x0*x23 + x26);
                const T x28 = 3*x3;
                const T x29 = PGx[g]*PGy[g]*PGz[g];
                const T x30 = -x29*(x1 + x28);
                const T x31 = PGz[g]*PGz[g];
                const T x32 = scaled_boys_5[g]*x31;
                const T x33 = x3 + x32;
                const T x34 = x0*x33;
                const T x35 = scaled_boys_3[g]*x31 - scaled_boys_4[g]*x31 + x9;
                const T x36 = tpi*x35;
                const T x37 = 3*x36;
                const T x38 = -PGx[g]*(x34 + x37);
                const T x39 = x23 + x5;
                const T x40 = tpi*(x14 + x25);
                const T x41 = -PGy[g]*(x0*x39 + x40);
                const T x42 = x15 + x22;
                const T x43 = tpi*(x17 + x24);
                const T x44 = -PGz[g]*(x0*x42 + x43);
                const T x45 = -PGy[g]*(x34 + x36);
                const T x46 = x33 + x5;
                const T x47 = tpi*(x14 + x35);
                const T x48 = -PGz[g]*(x0*x46 + x47);
                const T x49 = -PGx[g]*(x21*x39 + x26);
                const T x50 = x16 + x42;
                const T x51 = -x29*x50;
                const T x52 = x21*x33;
                const T x53 = -PGx[g]*(x36 + x52);
                const T x54 = -x29*(x28 + x32);
                const T x55 = x31*x46 + x37;
                const T x56 = -PGx[g]*x55;
                const T x57 = PGz[g]*(-x21*x50 - 3*x43);
                const T x58 = -PGy[g]*(x37 + x52);
                const T x59 = -PGz[g]*(x21*x46 + x47);
                const T x60 = -PGy[g]*x55;

                integral_0[g] = PGx[g]*(-4*tpi*(x10 + x14) + x0*(-x1 - x3 - x5) - x10*x11);
                integral_1[g] = x19;
                integral_2[g] = x20;
                integral_3[g] = x27;
                integral_4[g] = x30;
                integral_5[g] = x38;
                integral_6[g] = x19;
                integral_7[g] = x27;
                integral_8[g] = x30;
                integral_9[g] = x41;
                integral_10[g] = x44;
                integral_11[g] = x45;
                integral_12[g] = x20;
                integral_13[g] = x30;
                integral_14[g] = x38;
                integral_15[g] = x44;
                integral_16[g] = x45;
                integral_17[g] = x48;
                integral_18[g] = x27;
                integral_19[g] = x41;
                integral_20[g] = x44;
                integral_21[g] = x49;
                integral_22[g] = x51;
                integral_23[g] = x53;
                integral_24[g] = x30;
                integral_25[g] = x44;
                integral_26[g] = x45;
                integral_27[g] = x51;
                integral_28[g] = x53;
                integral_29[g] = x54;
                integral_30[g] = x38;
                integral_31[g] = x45;
                integral_32[g] = x48;
                integral_33[g] = x53;
                integral_34[g] = x54;
                integral_35[g] = x56;
                integral_36[g] = x41;
                integral_37[g] = x49;
                integral_38[g] = x51;
                integral_39[g] = PGy[g]*(-x21*x39 - x26 - 4*x40);
                integral_40[g] = x57;
                integral_41[g] = x58;
                integral_42[g] = x44;
                integral_43[g] = x51;
                integral_44[g] = x53;
                integral_45[g] = x57;
                integral_46[g] = x58;
                integral_47[g] = x59;
                integral_48[g] = x45;
                integral_49[g] = x53;
                integral_50[g] = x54;
                integral_51[g] = x58;
                integral_52[g] = x59;
                integral_53[g] = x60;
                integral_54[g] = x48;
                integral_55[g] = x54;
                integral_56[g] = x56;
                integral_57[g] = x59;
                integral_58[g] = x60;
                integral_59[g] = PGz[g]*(-x31*x46 - x37 - 4*x47);

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*PGx[g];
                const T x1 = scaled_boys_5[g]*x0;
                const T x2 = scaled_boys_3[g] - scaled_boys_4[g];
                const T x3 = tpi*x2;
                const T x4 = -tpi*x2;
                const T x5 = -2*x4;
                const T x6 = scaled_boys_3[g]*x0 - scaled_boys_4[g]*x0;
                const T x7 = scaled_boys_1[g] - scaled_boys_2[g];
                const T x8 = scaled_boys_2[g] - scaled_boys_3[g];
                const T x9 = tpi*x7 - tpi*x8;
                const T x10 = x6 + x9;
                const T x11 = 3*tpi;
                const T x12 = -tpi*x7;
                const T x13 = -tpi*x8;
                const T x14 = -2*x12 + 2*x13;
                const T x15 = -x4;
                const T x16 = 2*x3;
                const T x17 = -x12 + x13;
                const T x18 = x0*(-x1 - x15 - x16) - x11*(x17 + x6);
                const T x19 = PGy[g]*x18;
                const T x20 = PGz[g]*x18;
                const T x21 = PGy[g]*PGy[g];
                const T x22 = scaled_boys_5[g]*x21;
                const T x23 = x22 + x3;
                const T x24 = scaled_boys_3[g]*x21 - scaled_boys_4[g]*x21;
                const T x25 = x24 + x9;
                const T x26 = x11*x25;
                const T x27 = -PGx[g]*(x0*x23 + x26);
                const T x28 = 3*x3;
                const T x29 = PGx[g]*PGy[g]*PGz[g];
                const T x30 = -x29*(x1 + x28);
                const T x31 = PGz[g]*PGz[g];
                const T x32 = scaled_boys_5[g]*x31;
                const T x33 = x3 + x32;
                const T x34 = x0*x33;
                const T x35 = scaled_boys_3[g]*x31 - scaled_boys_4[g]*x31 + x9;
                const T x36 = tpi*x35;
                const T x37 = 3*x36;
                const T x38 = -PGx[g]*(x34 + x37);
                const T x39 = x23 + x5;
                const T x40 = tpi*(x14 + x25);
                const T x41 = -PGy[g]*(x0*x39 + x40);
                const T x42 = x15 + x22;
                const T x43 = tpi*(x17 + x24);
                const T x44 = -PGz[g]*(x0*x42 + x43);
                const T x45 = -PGy[g]*(x34 + x36);
                const T x46 = x33 + x5;
                const T x47 = tpi*(x14 + x35);
                const T x48 = -PGz[g]*(x0*x46 + x47);
                const T x49 = -PGx[g]*(x21*x39 + x26);
                const T x50 = x16 + x42;
                const T x51 = -x29*x50;
                const T x52 = x21*x33;
                const T x53 = -PGx[g]*(x36 + x52);
                const T x54 = -x29*(x28 + x32);
                const T x55 = x31*x46 + x37;
                const T x56 = -PGx[g]*x55;
                const T x57 = PGz[g]*(-x21*x50 - 3*x43);
                const T x58 = -PGy[g]*(x37 + x52);
                const T x59 = -PGz[g]*(x21*x46 + x47);
                const T x60 = -PGy[g]*x55;

                integral_0[g] += PGx[g]*(-4*tpi*(x10 + x14) + x0*(-x1 - x3 - x5) - x10*x11);
                integral_1[g] += x19;
                integral_2[g] += x20;
                integral_3[g] += x27;
                integral_4[g] += x30;
                integral_5[g] += x38;
                integral_6[g] += x19;
                integral_7[g] += x27;
                integral_8[g] += x30;
                integral_9[g] += x41;
                integral_10[g] += x44;
                integral_11[g] += x45;
                integral_12[g] += x20;
                integral_13[g] += x30;
                integral_14[g] += x38;
                integral_15[g] += x44;
                integral_16[g] += x45;
                integral_17[g] += x48;
                integral_18[g] += x27;
                integral_19[g] += x41;
                integral_20[g] += x44;
                integral_21[g] += x49;
                integral_22[g] += x51;
                integral_23[g] += x53;
                integral_24[g] += x30;
                integral_25[g] += x44;
                integral_26[g] += x45;
                integral_27[g] += x51;
                integral_28[g] += x53;
                integral_29[g] += x54;
                integral_30[g] += x38;
                integral_31[g] += x45;
                integral_32[g] += x48;
                integral_33[g] += x53;
                integral_34[g] += x54;
                integral_35[g] += x56;
                integral_36[g] += x41;
                integral_37[g] += x49;
                integral_38[g] += x51;
                integral_39[g] += PGy[g]*(-x21*x39 - x26 - 4*x40);
                integral_40[g] += x57;
                integral_41[g] += x58;
                integral_42[g] += x44;
                integral_43[g] += x51;
                integral_44[g] += x53;
                integral_45[g] += x57;
                integral_46[g] += x58;
                integral_47[g] += x59;
                integral_48[g] += x45;
                integral_49[g] += x53;
                integral_50[g] += x54;
                integral_51[g] += x58;
                integral_52[g] += x59;
                integral_53[g] += x60;
                integral_54[g] += x48;
                integral_55[g] += x54;
                integral_56[g] += x56;
                integral_57[g] += x59;
                integral_58[g] += x60;
                integral_59[g] += PGz[g]*(-x31*x46 - x37 - 4*x47);

            } // for g
        } // if/else prim
    } // for prim
} // VRR_F_D_SAMECENTER

template void VRR_F_D_SAMECENTER<float>(
        float* integrals, 
        const float* v_scaled_boys, 
        const float* v_p, 
        const float* v_PGx, 
        const float* v_PGy, 
        const float* v_PGz, 
        int nprim, 
        int gb_size);

template void VRR_F_D_SAMECENTER<double>(
        double* integrals, 
        const double* v_scaled_boys, 
        const double* v_p, 
        const double* v_PGx, 
        const double* v_PGy, 
        const double* v_PGz, 
        int nprim, 
        int gb_size);

