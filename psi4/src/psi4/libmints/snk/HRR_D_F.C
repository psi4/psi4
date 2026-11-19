#include "HRR_D_F.h"

template <typename T>
void HRR_D_F(
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
    T* __restrict base_integral_31 = &base_integrals[31*gb_size];
    T* __restrict base_integral_32 = &base_integrals[32*gb_size];
    T* __restrict base_integral_33 = &base_integrals[33*gb_size];
    T* __restrict base_integral_34 = &base_integrals[34*gb_size];
    T* __restrict base_integral_35 = &base_integrals[35*gb_size];
    T* __restrict base_integral_36 = &base_integrals[36*gb_size];
    T* __restrict base_integral_37 = &base_integrals[37*gb_size];
    T* __restrict base_integral_38 = &base_integrals[38*gb_size];
    T* __restrict base_integral_39 = &base_integrals[39*gb_size];
    T* __restrict base_integral_40 = &base_integrals[40*gb_size];
    T* __restrict base_integral_41 = &base_integrals[41*gb_size];
    T* __restrict base_integral_42 = &base_integrals[42*gb_size];
    T* __restrict base_integral_43 = &base_integrals[43*gb_size];
    T* __restrict base_integral_44 = &base_integrals[44*gb_size];
    T* __restrict base_integral_45 = &base_integrals[45*gb_size];

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
                const T x61 = 3*tpi;
                const T x62 = x3*(x10 + x5 - x9);
                const T x63 = PBx*x13;
                const T x64 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x65 = PBx*x11;
                const T x66 = PBx*scaled_boys_3[g];
                const T x67 = PGx[g]*scaled_boys_4[g];
                const T x68 = x66 - x67;
                const T x69 = PGx[g]*x68;
                const T x70 = x64 + x65 - x69;
                const T x71 = PGx[g]*x70;
                const T x72 = x62 + x63 - x71;
                const T x73 = PBx*x14 - PGx[g]*x72 + x61*(x12 + x6 + x60 - x8);
                const T x74 = tpi*(x19 - x22 + x23);
                const T x75 = PBx*x26;
                const T x76 = PBx*x24;
                const T x77 = PBy*scaled_boys_3[g];
                const T x78 = PGy[g]*scaled_boys_4[g];
                const T x79 = x77 - x78;
                const T x80 = PGx[g]*x79;
                const T x81 = x76 - x80;
                const T x82 = PGx[g]*x81;
                const T x83 = x74 + x75 - x82;
                const T x84 = PBx*x27 - PGx[g]*x83 + x3*(x20 - x21 + x25);
                const T x85 = tpi*(x32 - x35 + x36);
                const T x86 = PBx*x39;
                const T x87 = PBx*x37;
                const T x88 = PBz*scaled_boys_3[g];
                const T x89 = PGz[g]*scaled_boys_4[g];
                const T x90 = x88 - x89;
                const T x91 = PGx[g]*x90;
                const T x92 = x87 - x91;
                const T x93 = PGx[g]*x92;
                const T x94 = x85 + x86 - x93;
                const T x95 = PBx*x40 - PGx[g]*x94 + x3*(x33 - x34 + x38);
                const T x96 = tpi*(x41 - x42 + x43 + x60);
                const T x97 = PBx*x44;
                const T x98 = PBy*x24;
                const T x99 = PGy[g]*x79;
                const T x100 = x64 + x98 - x99;
                const T x101 = PGx[g]*x100;
                const T x102 = -x101 + x97;
                const T x103 = PBx*x45 - PGx[g]*x102 + x96;
                const T x104 = tpi*(x46 - x47 + x48);
                const T x105 = PBx*x49;
                const T x106 = PBy*x37;
                const T x107 = PGy[g]*x90;
                const T x108 = x106 - x107;
                const T x109 = PGx[g]*x108;
                const T x110 = x105 - x109;
                const T x111 = PBx*x50 - PGx[g]*x110 + x104;
                const T x112 = tpi*(x51 - x52 + x53 + x60);
                const T x113 = PBx*x54;
                const T x114 = PBz*x37;
                const T x115 = PGz[g]*x90;
                const T x116 = x114 - x115 + x64;
                const T x117 = PGx[g]*x116;
                const T x118 = x113 - x117;
                const T x119 = PBx*x55 - PGx[g]*x118 + x112;
                const T x120 = 2*x74;
                const T x121 = PBy*x44;
                const T x122 = PGy[g]*x100;
                const T x123 = x120 + x121 - x122;
                const T x124 = PBx*x56 - PGx[g]*x123;
                const T x125 = PBy*x49;
                const T x126 = PGy[g]*x108;
                const T x127 = x125 - x126 + x85;
                const T x128 = PBx*x57 - PGx[g]*x127;
                const T x129 = PBy*x54;
                const T x130 = PGy[g]*x116;
                const T x131 = x129 - x130;
                const T x132 = PBx*x58 - PGx[g]*x131;
                const T x133 = 2*x85;
                const T x134 = PBz*x54;
                const T x135 = PGz[g]*x116;
                const T x136 = x133 + x134 - x135;
                const T x137 = PBx*x59 - PGx[g]*x136;
                const T x138 = PBy*x56 - PGy[g]*x123 + 3*x96;
                const T x139 = PBy*x57 - PGy[g]*x127 + 2*x104;
                const T x140 = PBy*x58 - PGy[g]*x131 + x112;
                const T x141 = PBy*x59 - PGy[g]*x136;
                const T x142 = PBz*x59 - PGz[g]*x136 + 3*x112;
                const T x143 = -x64;
                const T x144 = tpi*(scaled_boys_3[g] - scaled_boys_4[g]);
                const T x145 = tpi*(x24 - x77 + x78);
                const T x146 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                const T x147 = -x85;
                const T x148 = tpi*(x37 - x88 + x89);
                const T x149 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                const T x150 = tpi*(x143 + x44 - x98 + x99);
                const T x151 = PBy*x79 - PGy[g]*x146 + x144;
                const T x152 = tpi*(-x106 + x107 + x49);
                const T x153 = PBy*x90 - PGy[g]*x149;
                const T x154 = tpi*(-x114 + x115 + x143 + x54);
                const T x155 = PBz*x90 - PGz[g]*x149 + x144;
                const T x156 = tpi*(-x120 - x121 + x122 + x56);
                const T x157 = PBy*x100 - PGy[g]*x151 + 2*x145;
                const T x158 = tpi*(-x125 + x126 + x147 + x57);
                const T x159 = PBy*x108 - PGy[g]*x153 + x148;
                const T x160 = tpi*(-x129 + x130 + x58);
                const T x161 = PBy*x116 - PGy[g]*x155;
                const T x162 = tpi*(-x133 - x134 + x135 + x59);
                const T x163 = PBz*x116 - PGz[g]*x155 + 2*x148;
                const T x164 = PBy*x123 - PGy[g]*x157 + 3*x150;
                const T x165 = PBy*x127 - PGy[g]*x159 + 2*x152;
                const T x166 = PBy*x131 - PGy[g]*x161 + x154;
                const T x167 = PBy*x136 - PGy[g]*x163;
                const T x168 = PBz*x136 - PGz[g]*x163 + 3*x154;

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
                base_integral_10[g] = x73;
                base_integral_11[g] = x84;
                base_integral_12[g] = x95;
                base_integral_13[g] = x103;
                base_integral_14[g] = x111;
                base_integral_15[g] = x119;
                base_integral_16[g] = x124;
                base_integral_17[g] = x128;
                base_integral_18[g] = x132;
                base_integral_19[g] = x137;
                base_integral_20[g] = x138;
                base_integral_21[g] = x139;
                base_integral_22[g] = x140;
                base_integral_23[g] = x141;
                base_integral_24[g] = x142;
                base_integral_25[g] = PBx*x73 - PGx[g]*(PBx*x72 - PGx[g]*(PBx*x70 - PGx[g]*(PBx*x68 - PGx[g]*(PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g]) + x144) + x3*(x11 - x66 + x67)) + x61*(x13 + x143 - x65 + x69)) + 4*tpi*(x14 - x62 - x63 + x71);
                base_integral_26[g] = PBx*x84 - PGx[g]*(PBx*x83 - PGx[g]*(PBx*x81 - PGx[g]*(PBx*x79 - PGx[g]*x146) + x145) + x3*(x26 - x76 + x80)) + x61*(x27 - x74 - x75 + x82);
                base_integral_27[g] = PBx*x95 - PGx[g]*(PBx*x94 - PGx[g]*(PBx*x92 - PGx[g]*(PBx*x90 - PGx[g]*x149) + x148) + x3*(x39 - x87 + x91)) + x61*(x147 + x40 - x86 + x93);
                base_integral_28[g] = PBx*x103 - PGx[g]*(PBx*x102 - PGx[g]*(PBx*x100 - PGx[g]*x151) + x150) + x3*(x101 + x45 - x97);
                base_integral_29[g] = PBx*x111 - PGx[g]*(PBx*x110 - PGx[g]*(PBx*x108 - PGx[g]*x153) + x152) + x3*(-x105 + x109 + x50);
                base_integral_30[g] = PBx*x119 - PGx[g]*(PBx*x118 - PGx[g]*(PBx*x116 - PGx[g]*x155) + x154) + x3*(-x113 + x117 + x55);
                base_integral_31[g] = PBx*x124 - PGx[g]*(PBx*x123 - PGx[g]*x157) + x156;
                base_integral_32[g] = PBx*x128 - PGx[g]*(PBx*x127 - PGx[g]*x159) + x158;
                base_integral_33[g] = PBx*x132 - PGx[g]*(PBx*x131 - PGx[g]*x161) + x160;
                base_integral_34[g] = PBx*x137 - PGx[g]*(PBx*x136 - PGx[g]*x163) + x162;
                base_integral_35[g] = PBx*x138 - PGx[g]*x164;
                base_integral_36[g] = PBx*x139 - PGx[g]*x165;
                base_integral_37[g] = PBx*x140 - PGx[g]*x166;
                base_integral_38[g] = PBx*x141 - PGx[g]*x167;
                base_integral_39[g] = PBx*x142 - PGx[g]*x168;
                base_integral_40[g] = PBy*x138 - PGy[g]*x164 + 4*x156;
                base_integral_41[g] = PBy*x139 - PGy[g]*x165 + 3*x158;
                base_integral_42[g] = PBy*x140 - PGy[g]*x166 + 2*x160;
                base_integral_43[g] = PBy*x141 - PGy[g]*x167 + x162;
                base_integral_44[g] = PBy*x142 - PGy[g]*x168;
                base_integral_45[g] = PBz*x142 - PGz[g]*x168 + 4*x162;

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
                const T x61 = 3*tpi;
                const T x62 = x3*(x10 + x5 - x9);
                const T x63 = PBx*x13;
                const T x64 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x65 = PBx*x11;
                const T x66 = PBx*scaled_boys_3[g];
                const T x67 = PGx[g]*scaled_boys_4[g];
                const T x68 = x66 - x67;
                const T x69 = PGx[g]*x68;
                const T x70 = x64 + x65 - x69;
                const T x71 = PGx[g]*x70;
                const T x72 = x62 + x63 - x71;
                const T x73 = PBx*x14 - PGx[g]*x72 + x61*(x12 + x6 + x60 - x8);
                const T x74 = tpi*(x19 - x22 + x23);
                const T x75 = PBx*x26;
                const T x76 = PBx*x24;
                const T x77 = PBy*scaled_boys_3[g];
                const T x78 = PGy[g]*scaled_boys_4[g];
                const T x79 = x77 - x78;
                const T x80 = PGx[g]*x79;
                const T x81 = x76 - x80;
                const T x82 = PGx[g]*x81;
                const T x83 = x74 + x75 - x82;
                const T x84 = PBx*x27 - PGx[g]*x83 + x3*(x20 - x21 + x25);
                const T x85 = tpi*(x32 - x35 + x36);
                const T x86 = PBx*x39;
                const T x87 = PBx*x37;
                const T x88 = PBz*scaled_boys_3[g];
                const T x89 = PGz[g]*scaled_boys_4[g];
                const T x90 = x88 - x89;
                const T x91 = PGx[g]*x90;
                const T x92 = x87 - x91;
                const T x93 = PGx[g]*x92;
                const T x94 = x85 + x86 - x93;
                const T x95 = PBx*x40 - PGx[g]*x94 + x3*(x33 - x34 + x38);
                const T x96 = tpi*(x41 - x42 + x43 + x60);
                const T x97 = PBx*x44;
                const T x98 = PBy*x24;
                const T x99 = PGy[g]*x79;
                const T x100 = x64 + x98 - x99;
                const T x101 = PGx[g]*x100;
                const T x102 = -x101 + x97;
                const T x103 = PBx*x45 - PGx[g]*x102 + x96;
                const T x104 = tpi*(x46 - x47 + x48);
                const T x105 = PBx*x49;
                const T x106 = PBy*x37;
                const T x107 = PGy[g]*x90;
                const T x108 = x106 - x107;
                const T x109 = PGx[g]*x108;
                const T x110 = x105 - x109;
                const T x111 = PBx*x50 - PGx[g]*x110 + x104;
                const T x112 = tpi*(x51 - x52 + x53 + x60);
                const T x113 = PBx*x54;
                const T x114 = PBz*x37;
                const T x115 = PGz[g]*x90;
                const T x116 = x114 - x115 + x64;
                const T x117 = PGx[g]*x116;
                const T x118 = x113 - x117;
                const T x119 = PBx*x55 - PGx[g]*x118 + x112;
                const T x120 = 2*x74;
                const T x121 = PBy*x44;
                const T x122 = PGy[g]*x100;
                const T x123 = x120 + x121 - x122;
                const T x124 = PBx*x56 - PGx[g]*x123;
                const T x125 = PBy*x49;
                const T x126 = PGy[g]*x108;
                const T x127 = x125 - x126 + x85;
                const T x128 = PBx*x57 - PGx[g]*x127;
                const T x129 = PBy*x54;
                const T x130 = PGy[g]*x116;
                const T x131 = x129 - x130;
                const T x132 = PBx*x58 - PGx[g]*x131;
                const T x133 = 2*x85;
                const T x134 = PBz*x54;
                const T x135 = PGz[g]*x116;
                const T x136 = x133 + x134 - x135;
                const T x137 = PBx*x59 - PGx[g]*x136;
                const T x138 = PBy*x56 - PGy[g]*x123 + 3*x96;
                const T x139 = PBy*x57 - PGy[g]*x127 + 2*x104;
                const T x140 = PBy*x58 - PGy[g]*x131 + x112;
                const T x141 = PBy*x59 - PGy[g]*x136;
                const T x142 = PBz*x59 - PGz[g]*x136 + 3*x112;
                const T x143 = -x64;
                const T x144 = tpi*(scaled_boys_3[g] - scaled_boys_4[g]);
                const T x145 = tpi*(x24 - x77 + x78);
                const T x146 = PBy*scaled_boys_4[g] - PGy[g]*scaled_boys_5[g];
                const T x147 = -x85;
                const T x148 = tpi*(x37 - x88 + x89);
                const T x149 = PBz*scaled_boys_4[g] - PGz[g]*scaled_boys_5[g];
                const T x150 = tpi*(x143 + x44 - x98 + x99);
                const T x151 = PBy*x79 - PGy[g]*x146 + x144;
                const T x152 = tpi*(-x106 + x107 + x49);
                const T x153 = PBy*x90 - PGy[g]*x149;
                const T x154 = tpi*(-x114 + x115 + x143 + x54);
                const T x155 = PBz*x90 - PGz[g]*x149 + x144;
                const T x156 = tpi*(-x120 - x121 + x122 + x56);
                const T x157 = PBy*x100 - PGy[g]*x151 + 2*x145;
                const T x158 = tpi*(-x125 + x126 + x147 + x57);
                const T x159 = PBy*x108 - PGy[g]*x153 + x148;
                const T x160 = tpi*(-x129 + x130 + x58);
                const T x161 = PBy*x116 - PGy[g]*x155;
                const T x162 = tpi*(-x133 - x134 + x135 + x59);
                const T x163 = PBz*x116 - PGz[g]*x155 + 2*x148;
                const T x164 = PBy*x123 - PGy[g]*x157 + 3*x150;
                const T x165 = PBy*x127 - PGy[g]*x159 + 2*x152;
                const T x166 = PBy*x131 - PGy[g]*x161 + x154;
                const T x167 = PBy*x136 - PGy[g]*x163;
                const T x168 = PBz*x136 - PGz[g]*x163 + 3*x154;

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
                base_integral_10[g] += x73;
                base_integral_11[g] += x84;
                base_integral_12[g] += x95;
                base_integral_13[g] += x103;
                base_integral_14[g] += x111;
                base_integral_15[g] += x119;
                base_integral_16[g] += x124;
                base_integral_17[g] += x128;
                base_integral_18[g] += x132;
                base_integral_19[g] += x137;
                base_integral_20[g] += x138;
                base_integral_21[g] += x139;
                base_integral_22[g] += x140;
                base_integral_23[g] += x141;
                base_integral_24[g] += x142;
                base_integral_25[g] += PBx*x73 - PGx[g]*(PBx*x72 - PGx[g]*(PBx*x70 - PGx[g]*(PBx*x68 - PGx[g]*(PBx*scaled_boys_4[g] - PGx[g]*scaled_boys_5[g]) + x144) + x3*(x11 - x66 + x67)) + x61*(x13 + x143 - x65 + x69)) + 4*tpi*(x14 - x62 - x63 + x71);
                base_integral_26[g] += PBx*x84 - PGx[g]*(PBx*x83 - PGx[g]*(PBx*x81 - PGx[g]*(PBx*x79 - PGx[g]*x146) + x145) + x3*(x26 - x76 + x80)) + x61*(x27 - x74 - x75 + x82);
                base_integral_27[g] += PBx*x95 - PGx[g]*(PBx*x94 - PGx[g]*(PBx*x92 - PGx[g]*(PBx*x90 - PGx[g]*x149) + x148) + x3*(x39 - x87 + x91)) + x61*(x147 + x40 - x86 + x93);
                base_integral_28[g] += PBx*x103 - PGx[g]*(PBx*x102 - PGx[g]*(PBx*x100 - PGx[g]*x151) + x150) + x3*(x101 + x45 - x97);
                base_integral_29[g] += PBx*x111 - PGx[g]*(PBx*x110 - PGx[g]*(PBx*x108 - PGx[g]*x153) + x152) + x3*(-x105 + x109 + x50);
                base_integral_30[g] += PBx*x119 - PGx[g]*(PBx*x118 - PGx[g]*(PBx*x116 - PGx[g]*x155) + x154) + x3*(-x113 + x117 + x55);
                base_integral_31[g] += PBx*x124 - PGx[g]*(PBx*x123 - PGx[g]*x157) + x156;
                base_integral_32[g] += PBx*x128 - PGx[g]*(PBx*x127 - PGx[g]*x159) + x158;
                base_integral_33[g] += PBx*x132 - PGx[g]*(PBx*x131 - PGx[g]*x161) + x160;
                base_integral_34[g] += PBx*x137 - PGx[g]*(PBx*x136 - PGx[g]*x163) + x162;
                base_integral_35[g] += PBx*x138 - PGx[g]*x164;
                base_integral_36[g] += PBx*x139 - PGx[g]*x165;
                base_integral_37[g] += PBx*x140 - PGx[g]*x166;
                base_integral_38[g] += PBx*x141 - PGx[g]*x167;
                base_integral_39[g] += PBx*x142 - PGx[g]*x168;
                base_integral_40[g] += PBy*x138 - PGy[g]*x164 + 4*x156;
                base_integral_41[g] += PBy*x139 - PGy[g]*x165 + 3*x158;
                base_integral_42[g] += PBy*x140 - PGy[g]*x166 + 2*x160;
                base_integral_43[g] += PBy*x141 - PGy[g]*x167 + x162;
                base_integral_44[g] += PBy*x142 - PGy[g]*x168;
                base_integral_45[g] += PBz*x142 - PGz[g]*x168 + 4*x162;

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
        const T x9 = ABy*base_integral_0[g] + x0;
        const T x10 = ABy*base_integral_1[g] + x2;
        const T x11 = -ABy*base_integral_11[g] + base_integral_28[g];
        const T x12 = ABy*base_integral_2[g] + x3;
        const T x13 = ABy*base_integral_3[g] + x5;
        const T x14 = -ABy*base_integral_13[g] + base_integral_31[g];
        const T x15 = ABy*base_integral_4[g] + x6;
        const T x16 = -ABy*base_integral_14[g] + base_integral_32[g];
        const T x17 = ABy*base_integral_5[g] + x7;
        const T x18 = ABy*base_integral_6[g] - base_integral_20[g];
        const T x19 = -ABy*base_integral_16[g] + base_integral_35[g];
        const T x20 = -base_integral_21[g];
        const T x21 = ABy*base_integral_7[g] + x20;
        const T x22 = -ABy*base_integral_17[g] + base_integral_36[g];
        const T x23 = -base_integral_22[g];
        const T x24 = ABy*base_integral_8[g] + x23;
        const T x25 = -ABy*base_integral_18[g] + base_integral_37[g];
        const T x26 = -base_integral_23[g];
        const T x27 = ABy*base_integral_9[g] + x26;
        const T x28 = ABz*base_integral_0[g] + x1;
        const T x29 = ABz*base_integral_1[g] + x3;
        const T x30 = -ABz*base_integral_11[g] + base_integral_29[g];
        const T x31 = ABz*base_integral_2[g] + x4;
        const T x32 = -ABz*base_integral_12[g] + base_integral_30[g];
        const T x33 = ABz*base_integral_3[g] + x6;
        const T x34 = -ABz*base_integral_13[g] + base_integral_32[g];
        const T x35 = ABz*base_integral_4[g] + x7;
        const T x36 = -ABz*base_integral_14[g] + base_integral_33[g];
        const T x37 = ABz*base_integral_5[g] + x8;
        const T x38 = -ABz*base_integral_15[g] + base_integral_34[g];
        const T x39 = ABz*base_integral_6[g] + x20;
        const T x40 = -ABz*base_integral_16[g] + base_integral_36[g];
        const T x41 = ABz*base_integral_7[g] + x23;
        const T x42 = -ABz*base_integral_17[g] + base_integral_37[g];
        const T x43 = ABz*base_integral_8[g] + x26;
        const T x44 = -ABz*base_integral_18[g] + base_integral_38[g];
        const T x45 = ABz*base_integral_9[g] - base_integral_24[g];
        const T x46 = -ABz*base_integral_19[g] + base_integral_39[g];
        const T x47 = -ABz*base_integral_21[g] + base_integral_42[g];
        const T x48 = -ABz*base_integral_22[g] + base_integral_43[g];
        const T x49 = -ABz*base_integral_23[g] + base_integral_44[g];

        integral_0[g] = -ABx*base_integral_10[g] + ABx*(ABx*base_integral_0[g] - base_integral_10[g]) + base_integral_25[g];
        integral_1[g] = -ABx*base_integral_11[g] + ABx*(ABx*base_integral_1[g] + x0) + base_integral_26[g];
        integral_2[g] = -ABx*base_integral_12[g] + ABx*(ABx*base_integral_2[g] + x1) + base_integral_27[g];
        integral_3[g] = -ABx*base_integral_13[g] + ABx*(ABx*base_integral_3[g] + x2) + base_integral_28[g];
        integral_4[g] = -ABx*base_integral_14[g] + ABx*(ABx*base_integral_4[g] + x3) + base_integral_29[g];
        integral_5[g] = -ABx*base_integral_15[g] + ABx*(ABx*base_integral_5[g] + x4) + base_integral_30[g];
        integral_6[g] = -ABx*base_integral_16[g] + ABx*(ABx*base_integral_6[g] + x5) + base_integral_31[g];
        integral_7[g] = -ABx*base_integral_17[g] + ABx*(ABx*base_integral_7[g] + x6) + base_integral_32[g];
        integral_8[g] = -ABx*base_integral_18[g] + ABx*(ABx*base_integral_8[g] + x7) + base_integral_33[g];
        integral_9[g] = -ABx*base_integral_19[g] + ABx*(ABx*base_integral_9[g] + x8) + base_integral_34[g];
        integral_10[g] = ABx*x9 - ABy*base_integral_10[g] + base_integral_26[g];
        integral_11[g] = ABx*x10 + x11;
        integral_12[g] = ABx*x12 - ABy*base_integral_12[g] + base_integral_29[g];
        integral_13[g] = ABx*x13 + x14;
        integral_14[g] = ABx*x15 + x16;
        integral_15[g] = ABx*x17 - ABy*base_integral_15[g] + base_integral_33[g];
        integral_16[g] = ABx*x18 + x19;
        integral_17[g] = ABx*x21 + x22;
        integral_18[g] = ABx*x24 + x25;
        integral_19[g] = ABx*x27 - ABy*base_integral_19[g] + base_integral_38[g];
        integral_20[g] = ABx*x28 - ABz*base_integral_10[g] + base_integral_27[g];
        integral_21[g] = ABx*x29 + x30;
        integral_22[g] = ABx*x31 + x32;
        integral_23[g] = ABx*x33 + x34;
        integral_24[g] = ABx*x35 + x36;
        integral_25[g] = ABx*x37 + x38;
        integral_26[g] = ABx*x39 + x40;
        integral_27[g] = ABx*x41 + x42;
        integral_28[g] = ABx*x43 + x44;
        integral_29[g] = ABx*x45 + x46;
        integral_30[g] = ABy*x9 + x11;
        integral_31[g] = ABy*x10 + x14;
        integral_32[g] = ABy*x12 + x16;
        integral_33[g] = ABy*x13 + x19;
        integral_34[g] = ABy*x15 + x22;
        integral_35[g] = ABy*x17 + x25;
        integral_36[g] = -ABy*base_integral_20[g] + ABy*x18 + base_integral_40[g];
        integral_37[g] = -ABy*base_integral_21[g] + ABy*x21 + base_integral_41[g];
        integral_38[g] = -ABy*base_integral_22[g] + ABy*x24 + base_integral_42[g];
        integral_39[g] = -ABy*base_integral_23[g] + ABy*x27 + base_integral_43[g];
        integral_40[g] = ABy*x28 + x30;
        integral_41[g] = ABy*x29 + x34;
        integral_42[g] = ABy*x31 + x36;
        integral_43[g] = ABy*x33 + x40;
        integral_44[g] = ABy*x35 + x42;
        integral_45[g] = ABy*x37 + x44;
        integral_46[g] = ABy*x39 - ABz*base_integral_20[g] + base_integral_41[g];
        integral_47[g] = ABy*x41 + x47;
        integral_48[g] = ABy*x43 + x48;
        integral_49[g] = ABy*x45 + x49;
        integral_50[g] = ABz*x28 + x32;
        integral_51[g] = ABz*x29 + x36;
        integral_52[g] = ABz*x31 + x38;
        integral_53[g] = ABz*x33 + x42;
        integral_54[g] = ABz*x35 + x44;
        integral_55[g] = ABz*x37 + x46;
        integral_56[g] = ABz*x39 + x47;
        integral_57[g] = ABz*x41 + x48;
        integral_58[g] = ABz*x43 + x49;
        integral_59[g] = -ABz*base_integral_24[g] + ABz*x45 + base_integral_45[g];

    } // for g
} // HRR_D_F

template void HRR_D_F<float>(
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

template void HRR_D_F<double>(
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

