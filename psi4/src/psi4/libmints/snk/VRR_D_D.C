#include "VRR_D_D.h"

template <typename T>
void VRR_D_D(
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
                const T x5 = PBx*x2;
                const T x6 = -x0 + x1;
                const T x7 = PGx[g]*x6;
                const T x8 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x9 = PBx*x6;
                const T x10 = PBx*scaled_boys_2[g];
                const T x11 = PGx[g]*scaled_boys_3[g];
                const T x12 = x10 - x11;
                const T x13 = PGx[g]*x12;
                const T x14 = -x13 + x8 + x9;
                const T x15 = 2*PAx;
                const T x16 = 3*x4 - 3*x8;
                const T x17 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x18 = PBy*scaled_boys_0[g];
                const T x19 = PGy[g]*scaled_boys_1[g];
                const T x20 = -x19;
                const T x21 = x18 + x20;
                const T x22 = PBx*x21;
                const T x23 = PBy*scaled_boys_1[g];
                const T x24 = PGy[g]*scaled_boys_2[g];
                const T x25 = -x24;
                const T x26 = x23 + x25;
                const T x27 = PBx*x26;
                const T x28 = PGx[g]*x26;
                const T x29 = PBy*scaled_boys_2[g];
                const T x30 = PGy[g]*scaled_boys_3[g];
                const T x31 = -x30;
                const T x32 = x29 + x31;
                const T x33 = PGx[g]*x32;
                const T x34 = -x23;
                const T x35 = tpi*(x21 + x24 + x34);
                const T x36 = x27 - x33;
                const T x37 = -x29;
                const T x38 = tpi*(x26 + x30 + x37);
                const T x39 = -PGy[g]*scaled_boys_4[g];
                const T x40 = PBy*scaled_boys_3[g] + x39;
                const T x41 = PBz*scaled_boys_0[g];
                const T x42 = PGz[g]*scaled_boys_1[g];
                const T x43 = -x42;
                const T x44 = x41 + x43;
                const T x45 = PBx*x44;
                const T x46 = PBz*scaled_boys_1[g];
                const T x47 = PGz[g]*scaled_boys_2[g];
                const T x48 = -x47;
                const T x49 = x46 + x48;
                const T x50 = PBx*x49;
                const T x51 = PGx[g]*x49;
                const T x52 = PBz*scaled_boys_2[g];
                const T x53 = PGz[g]*scaled_boys_3[g];
                const T x54 = -x53;
                const T x55 = x52 + x54;
                const T x56 = PGx[g]*x55;
                const T x57 = -x46;
                const T x58 = tpi*(x44 + x47 + x57);
                const T x59 = x50 - x56;
                const T x60 = -x52;
                const T x61 = tpi*(x49 + x53 + x60);
                const T x62 = -PGz[g]*scaled_boys_4[g];
                const T x63 = PBz*scaled_boys_3[g] + x62;
                const T x64 = PBy*x26;
                const T x65 = -x64;
                const T x66 = PBy*x21;
                const T x67 = PGy[g]*x26;
                const T x68 = x4 - x67;
                const T x69 = x66 + x68;
                const T x70 = PGy[g]*x32;
                const T x71 = -x8;
                const T x72 = x70 + x71;
                const T x73 = -x70 + x8;
                const T x74 = x64 + x73;
                const T x75 = -PGy[g]*x40 + x17;
                const T x76 = PBy*x32 + x75;
                const T x77 = PGy[g]*x55;
                const T x78 = PBy*x49;
                const T x79 = -x78;
                const T x80 = PBy*x44;
                const T x81 = PGy[g]*x49;
                const T x82 = -x81;
                const T x83 = x80 + x82;
                const T x84 = -x77;
                const T x85 = x78 + x84;
                const T x86 = -PGy[g]*x63;
                const T x87 = PBy*x55 + x86;
                const T x88 = PBz*x49;
                const T x89 = -x88;
                const T x90 = PBz*x44;
                const T x91 = PGz[g]*x49;
                const T x92 = x4 - x91;
                const T x93 = x90 + x92;
                const T x94 = PGz[g]*x55;
                const T x95 = x71 + x94;
                const T x96 = tpi*(x89 + x93 + x95);
                const T x97 = x8 - x94;
                const T x98 = x88 + x97;
                const T x99 = -PGz[g]*x63 + x17;
                const T x100 = PBz*x55 + x99;
                const T x101 = PAy*scaled_boys_2[g];
                const T x102 = x101 + x31;
                const T x103 = PGx[g]*x102;
                const T x104 = PAy*scaled_boys_1[g];
                const T x105 = x104 + x25;
                const T x106 = PBx*x105;
                const T x107 = PAy*scaled_boys_0[g];
                const T x108 = x107 + x20;
                const T x109 = PBx*x108 - PGx[g]*x105;
                const T x110 = -x104;
                const T x111 = -x103 + x106;
                const T x112 = -x101;
                const T x113 = PAy*scaled_boys_3[g] + x39;
                const T x114 = PAy*x26;
                const T x115 = PAy*x21;
                const T x116 = x115 + x68;
                const T x117 = x114 + x73;
                const T x118 = PAy*x32 + x75;
                const T x119 = PAy*x49;
                const T x120 = -x119;
                const T x121 = PAy*x44;
                const T x122 = x121 + x82;
                const T x123 = x119 + x84;
                const T x124 = PAy*x55 + x86;
                const T x125 = PAy*x69 - PGy[g]*x74 + 2*x35;
                const T x126 = PAy*x74 - PGy[g]*x76 + 2*x38;
                const T x127 = PAy*x83 - PGy[g]*x85 + x58;
                const T x128 = PAy*x85 - PGy[g]*x87 + x61;
                const T x129 = PAy*x93 - PGy[g]*x98;
                const T x130 = PAy*x98 - PGy[g]*x100;
                const T x131 = PAz*scaled_boys_2[g];
                const T x132 = x131 + x54;
                const T x133 = PGx[g]*x132;
                const T x134 = PAz*scaled_boys_1[g];
                const T x135 = x134 + x48;
                const T x136 = PBx*x135;
                const T x137 = PAz*scaled_boys_0[g];
                const T x138 = x137 + x43;
                const T x139 = PBx*x138 - PGx[g]*x135;
                const T x140 = -x134;
                const T x141 = tpi*(x138 + x140 + x47);
                const T x142 = -x133 + x136;
                const T x143 = -x131;
                const T x144 = tpi*(x135 + x143 + x53);
                const T x145 = PAz*scaled_boys_3[g] + x62;
                const T x146 = PGy[g]*x132;
                const T x147 = PBy*x135;
                const T x148 = -PGy[g]*x135;
                const T x149 = PBy*x138 + x148;
                const T x150 = tpi*(x146 - x147 + x149);
                const T x151 = -x146;
                const T x152 = x147 + x151;
                const T x153 = -PGy[g]*x145;
                const T x154 = PBy*x132 + x153;
                const T x155 = PAz*x49;
                const T x156 = PAz*x44;
                const T x157 = x156 + x92;
                const T x158 = tpi*(-x155 + x157 + x95);
                const T x159 = x155 + x97;
                const T x160 = PAz*x55 + x99;
                const T x161 = -PGy[g]*x152 + x141;
                const T x162 = PBy*x149 + x161;
                const T x163 = -PGy[g]*x154 + x144;
                const T x164 = PBy*x152 + x163;
                const T x165 = -PGy[g]*x159;
                const T x166 = PBy*x157 + x165;
                const T x167 = -PGy[g]*x160;
                const T x168 = PBy*x159 + x167;
                const T x169 = PAz*x93 - PGz[g]*x98 + 2*x58;
                const T x170 = PAz*x98 - PGz[g]*x100 + 2*x61;
                const T x171 = PGy[g]*x102;
                const T x172 = PAy*x105;
                const T x173 = PAy*x108 - PGy[g]*x105 + x4;
                const T x174 = -x171 + x172 + x8;
                const T x175 = 2*x24;
                const T x176 = PAy*x135;
                const T x177 = PAy*x138 + x148;
                const T x178 = x151 + x176;
                const T x179 = PGz[g]*x132;
                const T x180 = PAz*x135;
                const T x181 = PAz*x138 - PGz[g]*x135 + x4;
                const T x182 = tpi*(x179 - x180 + x181 + x71);
                const T x183 = -x179 + x180 + x8;
                const T x184 = PAz*x132 - PGz[g]*x145 + x17;
                const T x185 = PBy*x181 - PGy[g]*x183;
                const T x186 = PBy*x183 - PGy[g]*x184;
                const T x187 = 2*x47;
                const T x188 = PAz*x157 - PGz[g]*x159 + tpi*(x137 + x140 + x187 + x41 - 2*x42 + x57);
                const T x189 = PAz*x159 - PGz[g]*x160 + tpi*(x134 + x143 - x187 + x46 + 2*x53 + x60);

                integral_0[g] = PAx*(PAx*(x4 + x5 - x7) - PGx[g]*x14 + x3*(x0 - x1 + x2)) - PGx[g]*(PAx*x14 - PGx[g]*(PBx*x12 - PGx[g]*(PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x17) + x3*(-x10 + x11 + x6)) + tpi*(3*x13 + x15*x2 - x15*x6 + x16 + x5 - 3*x7 - x9);
                integral_1[g] = PAx*(PAx*(x22 - x28) - PGx[g]*x36 + x35) - PGx[g]*(PAx*x36 - PGx[g]*(PBx*x32 - PGx[g]*x40) + x38) + tpi*(PAx*x21 - PAx*x26 + x22 - x27 - 2*x28 + 2*x33);
                integral_2[g] = PAx*(PAx*(x45 - x51) - PGx[g]*x59 + x58) - PGx[g]*(PAx*x59 - PGx[g]*(PBx*x55 - PGx[g]*x63) + x61) + tpi*(PAx*x44 - PAx*x49 + x45 - x50 - 2*x51 + 2*x56);
                integral_3[g] = PAx*(PAx*x69 - PGx[g]*x74) - PGx[g]*(PAx*x74 - PGx[g]*x76) + tpi*(x65 + x69 + x72);
                integral_4[g] = PAx*(PAx*x83 - PGx[g]*x85) - PGx[g]*(PAx*x85 - PGx[g]*x87) + tpi*(x77 + x79 + x83);
                integral_5[g] = PAx*(PAx*x93 - PGx[g]*x98) - PGx[g]*(PAx*x98 - PGx[g]*x100) + x96;
                integral_6[g] = PAx*(PBx*x109 - PGx[g]*x111 + tpi*(x108 + x110 + x24)) - PGx[g]*(PBx*x111 - PGx[g]*(PBx*x102 - PGx[g]*x113) + tpi*(x105 + x112 + x30)) + x3*(x103 - x106 + x109);
                integral_7[g] = PAx*(PBx*x116 - PGx[g]*x117) - PGx[g]*(PBx*x117 - PGx[g]*x118) + tpi*(-x114 + x116 + x72);
                integral_8[g] = PAx*(PBx*x122 - PGx[g]*x123) - PGx[g]*(PBx*x123 - PGx[g]*x124) + tpi*(x120 + x122 + x77);
                integral_9[g] = PAx*x125 - PGx[g]*x126;
                integral_10[g] = PAx*x127 - PGx[g]*x128;
                integral_11[g] = PAx*x129 - PGx[g]*x130;
                integral_12[g] = PAx*(PBx*x139 - PGx[g]*x142 + x141) - PGx[g]*(PBx*x142 - PGx[g]*(PBx*x132 - PGx[g]*x145) + x144) + x3*(x133 - x136 + x139);
                integral_13[g] = PAx*(PBx*x149 - PGx[g]*x152) - PGx[g]*(PBx*x152 - PGx[g]*x154) + x150;
                integral_14[g] = PAx*(PBx*x157 - PGx[g]*x159) - PGx[g]*(PBx*x159 - PGx[g]*x160) + x158;
                integral_15[g] = PAx*x162 - PGx[g]*x164;
                integral_16[g] = PAx*x166 - PGx[g]*x168;
                integral_17[g] = PAx*x169 - PGx[g]*x170;
                integral_18[g] = PBx*(PBx*x173 - PGx[g]*x174) - PGx[g]*(PBx*x174 - PGx[g]*(PAy*x102 - PGy[g]*x113 + x17)) + tpi*(x171 - x172 + x173 + x71);
                integral_19[g] = PBx*(PAy*x116 - PGy[g]*x117 + tpi*(x107 + x110 + x175 + x18 - 2*x19 + x34)) - PGx[g]*(PAy*x117 - PGy[g]*x118 + tpi*(x104 + x112 - x175 + x23 + 2*x30 + x37));
                integral_20[g] = PBx*(PAy*x122 - PGy[g]*x123 + x58) - PGx[g]*(PAy*x123 - PGy[g]*x124 + x61);
                integral_21[g] = PAy*x125 - PGy[g]*x126 + tpi*(-2*x114 + 2*x115 + x16 + x65 + x66 - 3*x67 + 3*x70);
                integral_22[g] = PAy*x127 - PGy[g]*x128 + tpi*(x120 + x121 + 2*x77 + x79 + x80 - 2*x81);
                integral_23[g] = PAy*x129 - PGy[g]*x130 + x96;
                integral_24[g] = PBx*(PBx*x177 - PGx[g]*x178) - PGx[g]*(PBx*x178 - PGx[g]*(PAy*x132 + x153)) + tpi*(x146 - x176 + x177);
                integral_25[g] = PBx*(PAy*x149 + x161) - PGx[g]*(PAy*x152 + x163);
                integral_26[g] = PBx*(PAy*x157 + x165) - PGx[g]*(PAy*x159 + x167);
                integral_27[g] = PAy*x162 - PGy[g]*x164 + 2*x150;
                integral_28[g] = PAy*x166 - PGy[g]*x168 + x158;
                integral_29[g] = PAy*x169 - PGy[g]*x170;
                integral_30[g] = PBx*(PBx*x181 - PGx[g]*x183) - PGx[g]*(PBx*x183 - PGx[g]*x184) + x182;
                integral_31[g] = PBx*x185 - PGx[g]*x186;
                integral_32[g] = PBx*x188 - PGx[g]*x189;
                integral_33[g] = PBy*x185 - PGy[g]*x186 + x182;
                integral_34[g] = PBy*x188 - PGy[g]*x189;
                integral_35[g] = PAz*x169 - PGz[g]*x170 + tpi*(-2*x155 + 2*x156 + x16 + x89 + x90 - 3*x91 + 3*x94);

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PGx[g]*scaled_boys_2[g];
                const T x1 = PBx*scaled_boys_1[g];
                const T x2 = PBx*scaled_boys_0[g] - PGx[g]*scaled_boys_1[g];
                const T x3 = 2*tpi;
                const T x4 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x5 = PBx*x2;
                const T x6 = -x0 + x1;
                const T x7 = PGx[g]*x6;
                const T x8 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x9 = PBx*x6;
                const T x10 = PBx*scaled_boys_2[g];
                const T x11 = PGx[g]*scaled_boys_3[g];
                const T x12 = x10 - x11;
                const T x13 = PGx[g]*x12;
                const T x14 = -x13 + x8 + x9;
                const T x15 = 2*PAx;
                const T x16 = 3*x4 - 3*x8;
                const T x17 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x18 = PBy*scaled_boys_0[g];
                const T x19 = PGy[g]*scaled_boys_1[g];
                const T x20 = -x19;
                const T x21 = x18 + x20;
                const T x22 = PBx*x21;
                const T x23 = PBy*scaled_boys_1[g];
                const T x24 = PGy[g]*scaled_boys_2[g];
                const T x25 = -x24;
                const T x26 = x23 + x25;
                const T x27 = PBx*x26;
                const T x28 = PGx[g]*x26;
                const T x29 = PBy*scaled_boys_2[g];
                const T x30 = PGy[g]*scaled_boys_3[g];
                const T x31 = -x30;
                const T x32 = x29 + x31;
                const T x33 = PGx[g]*x32;
                const T x34 = -x23;
                const T x35 = tpi*(x21 + x24 + x34);
                const T x36 = x27 - x33;
                const T x37 = -x29;
                const T x38 = tpi*(x26 + x30 + x37);
                const T x39 = -PGy[g]*scaled_boys_4[g];
                const T x40 = PBy*scaled_boys_3[g] + x39;
                const T x41 = PBz*scaled_boys_0[g];
                const T x42 = PGz[g]*scaled_boys_1[g];
                const T x43 = -x42;
                const T x44 = x41 + x43;
                const T x45 = PBx*x44;
                const T x46 = PBz*scaled_boys_1[g];
                const T x47 = PGz[g]*scaled_boys_2[g];
                const T x48 = -x47;
                const T x49 = x46 + x48;
                const T x50 = PBx*x49;
                const T x51 = PGx[g]*x49;
                const T x52 = PBz*scaled_boys_2[g];
                const T x53 = PGz[g]*scaled_boys_3[g];
                const T x54 = -x53;
                const T x55 = x52 + x54;
                const T x56 = PGx[g]*x55;
                const T x57 = -x46;
                const T x58 = tpi*(x44 + x47 + x57);
                const T x59 = x50 - x56;
                const T x60 = -x52;
                const T x61 = tpi*(x49 + x53 + x60);
                const T x62 = -PGz[g]*scaled_boys_4[g];
                const T x63 = PBz*scaled_boys_3[g] + x62;
                const T x64 = PBy*x26;
                const T x65 = -x64;
                const T x66 = PBy*x21;
                const T x67 = PGy[g]*x26;
                const T x68 = x4 - x67;
                const T x69 = x66 + x68;
                const T x70 = PGy[g]*x32;
                const T x71 = -x8;
                const T x72 = x70 + x71;
                const T x73 = -x70 + x8;
                const T x74 = x64 + x73;
                const T x75 = -PGy[g]*x40 + x17;
                const T x76 = PBy*x32 + x75;
                const T x77 = PGy[g]*x55;
                const T x78 = PBy*x49;
                const T x79 = -x78;
                const T x80 = PBy*x44;
                const T x81 = PGy[g]*x49;
                const T x82 = -x81;
                const T x83 = x80 + x82;
                const T x84 = -x77;
                const T x85 = x78 + x84;
                const T x86 = -PGy[g]*x63;
                const T x87 = PBy*x55 + x86;
                const T x88 = PBz*x49;
                const T x89 = -x88;
                const T x90 = PBz*x44;
                const T x91 = PGz[g]*x49;
                const T x92 = x4 - x91;
                const T x93 = x90 + x92;
                const T x94 = PGz[g]*x55;
                const T x95 = x71 + x94;
                const T x96 = tpi*(x89 + x93 + x95);
                const T x97 = x8 - x94;
                const T x98 = x88 + x97;
                const T x99 = -PGz[g]*x63 + x17;
                const T x100 = PBz*x55 + x99;
                const T x101 = PAy*scaled_boys_2[g];
                const T x102 = x101 + x31;
                const T x103 = PGx[g]*x102;
                const T x104 = PAy*scaled_boys_1[g];
                const T x105 = x104 + x25;
                const T x106 = PBx*x105;
                const T x107 = PAy*scaled_boys_0[g];
                const T x108 = x107 + x20;
                const T x109 = PBx*x108 - PGx[g]*x105;
                const T x110 = -x104;
                const T x111 = -x103 + x106;
                const T x112 = -x101;
                const T x113 = PAy*scaled_boys_3[g] + x39;
                const T x114 = PAy*x26;
                const T x115 = PAy*x21;
                const T x116 = x115 + x68;
                const T x117 = x114 + x73;
                const T x118 = PAy*x32 + x75;
                const T x119 = PAy*x49;
                const T x120 = -x119;
                const T x121 = PAy*x44;
                const T x122 = x121 + x82;
                const T x123 = x119 + x84;
                const T x124 = PAy*x55 + x86;
                const T x125 = PAy*x69 - PGy[g]*x74 + 2*x35;
                const T x126 = PAy*x74 - PGy[g]*x76 + 2*x38;
                const T x127 = PAy*x83 - PGy[g]*x85 + x58;
                const T x128 = PAy*x85 - PGy[g]*x87 + x61;
                const T x129 = PAy*x93 - PGy[g]*x98;
                const T x130 = PAy*x98 - PGy[g]*x100;
                const T x131 = PAz*scaled_boys_2[g];
                const T x132 = x131 + x54;
                const T x133 = PGx[g]*x132;
                const T x134 = PAz*scaled_boys_1[g];
                const T x135 = x134 + x48;
                const T x136 = PBx*x135;
                const T x137 = PAz*scaled_boys_0[g];
                const T x138 = x137 + x43;
                const T x139 = PBx*x138 - PGx[g]*x135;
                const T x140 = -x134;
                const T x141 = tpi*(x138 + x140 + x47);
                const T x142 = -x133 + x136;
                const T x143 = -x131;
                const T x144 = tpi*(x135 + x143 + x53);
                const T x145 = PAz*scaled_boys_3[g] + x62;
                const T x146 = PGy[g]*x132;
                const T x147 = PBy*x135;
                const T x148 = -PGy[g]*x135;
                const T x149 = PBy*x138 + x148;
                const T x150 = tpi*(x146 - x147 + x149);
                const T x151 = -x146;
                const T x152 = x147 + x151;
                const T x153 = -PGy[g]*x145;
                const T x154 = PBy*x132 + x153;
                const T x155 = PAz*x49;
                const T x156 = PAz*x44;
                const T x157 = x156 + x92;
                const T x158 = tpi*(-x155 + x157 + x95);
                const T x159 = x155 + x97;
                const T x160 = PAz*x55 + x99;
                const T x161 = -PGy[g]*x152 + x141;
                const T x162 = PBy*x149 + x161;
                const T x163 = -PGy[g]*x154 + x144;
                const T x164 = PBy*x152 + x163;
                const T x165 = -PGy[g]*x159;
                const T x166 = PBy*x157 + x165;
                const T x167 = -PGy[g]*x160;
                const T x168 = PBy*x159 + x167;
                const T x169 = PAz*x93 - PGz[g]*x98 + 2*x58;
                const T x170 = PAz*x98 - PGz[g]*x100 + 2*x61;
                const T x171 = PGy[g]*x102;
                const T x172 = PAy*x105;
                const T x173 = PAy*x108 - PGy[g]*x105 + x4;
                const T x174 = -x171 + x172 + x8;
                const T x175 = 2*x24;
                const T x176 = PAy*x135;
                const T x177 = PAy*x138 + x148;
                const T x178 = x151 + x176;
                const T x179 = PGz[g]*x132;
                const T x180 = PAz*x135;
                const T x181 = PAz*x138 - PGz[g]*x135 + x4;
                const T x182 = tpi*(x179 - x180 + x181 + x71);
                const T x183 = -x179 + x180 + x8;
                const T x184 = PAz*x132 - PGz[g]*x145 + x17;
                const T x185 = PBy*x181 - PGy[g]*x183;
                const T x186 = PBy*x183 - PGy[g]*x184;
                const T x187 = 2*x47;
                const T x188 = PAz*x157 - PGz[g]*x159 + tpi*(x137 + x140 + x187 + x41 - 2*x42 + x57);
                const T x189 = PAz*x159 - PGz[g]*x160 + tpi*(x134 + x143 - x187 + x46 + 2*x53 + x60);

                integral_0[g] += PAx*(PAx*(x4 + x5 - x7) - PGx[g]*x14 + x3*(x0 - x1 + x2)) - PGx[g]*(PAx*x14 - PGx[g]*(PBx*x12 - PGx[g]*(PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x17) + x3*(-x10 + x11 + x6)) + tpi*(3*x13 + x15*x2 - x15*x6 + x16 + x5 - 3*x7 - x9);
                integral_1[g] += PAx*(PAx*(x22 - x28) - PGx[g]*x36 + x35) - PGx[g]*(PAx*x36 - PGx[g]*(PBx*x32 - PGx[g]*x40) + x38) + tpi*(PAx*x21 - PAx*x26 + x22 - x27 - 2*x28 + 2*x33);
                integral_2[g] += PAx*(PAx*(x45 - x51) - PGx[g]*x59 + x58) - PGx[g]*(PAx*x59 - PGx[g]*(PBx*x55 - PGx[g]*x63) + x61) + tpi*(PAx*x44 - PAx*x49 + x45 - x50 - 2*x51 + 2*x56);
                integral_3[g] += PAx*(PAx*x69 - PGx[g]*x74) - PGx[g]*(PAx*x74 - PGx[g]*x76) + tpi*(x65 + x69 + x72);
                integral_4[g] += PAx*(PAx*x83 - PGx[g]*x85) - PGx[g]*(PAx*x85 - PGx[g]*x87) + tpi*(x77 + x79 + x83);
                integral_5[g] += PAx*(PAx*x93 - PGx[g]*x98) - PGx[g]*(PAx*x98 - PGx[g]*x100) + x96;
                integral_6[g] += PAx*(PBx*x109 - PGx[g]*x111 + tpi*(x108 + x110 + x24)) - PGx[g]*(PBx*x111 - PGx[g]*(PBx*x102 - PGx[g]*x113) + tpi*(x105 + x112 + x30)) + x3*(x103 - x106 + x109);
                integral_7[g] += PAx*(PBx*x116 - PGx[g]*x117) - PGx[g]*(PBx*x117 - PGx[g]*x118) + tpi*(-x114 + x116 + x72);
                integral_8[g] += PAx*(PBx*x122 - PGx[g]*x123) - PGx[g]*(PBx*x123 - PGx[g]*x124) + tpi*(x120 + x122 + x77);
                integral_9[g] += PAx*x125 - PGx[g]*x126;
                integral_10[g] += PAx*x127 - PGx[g]*x128;
                integral_11[g] += PAx*x129 - PGx[g]*x130;
                integral_12[g] += PAx*(PBx*x139 - PGx[g]*x142 + x141) - PGx[g]*(PBx*x142 - PGx[g]*(PBx*x132 - PGx[g]*x145) + x144) + x3*(x133 - x136 + x139);
                integral_13[g] += PAx*(PBx*x149 - PGx[g]*x152) - PGx[g]*(PBx*x152 - PGx[g]*x154) + x150;
                integral_14[g] += PAx*(PBx*x157 - PGx[g]*x159) - PGx[g]*(PBx*x159 - PGx[g]*x160) + x158;
                integral_15[g] += PAx*x162 - PGx[g]*x164;
                integral_16[g] += PAx*x166 - PGx[g]*x168;
                integral_17[g] += PAx*x169 - PGx[g]*x170;
                integral_18[g] += PBx*(PBx*x173 - PGx[g]*x174) - PGx[g]*(PBx*x174 - PGx[g]*(PAy*x102 - PGy[g]*x113 + x17)) + tpi*(x171 - x172 + x173 + x71);
                integral_19[g] += PBx*(PAy*x116 - PGy[g]*x117 + tpi*(x107 + x110 + x175 + x18 - 2*x19 + x34)) - PGx[g]*(PAy*x117 - PGy[g]*x118 + tpi*(x104 + x112 - x175 + x23 + 2*x30 + x37));
                integral_20[g] += PBx*(PAy*x122 - PGy[g]*x123 + x58) - PGx[g]*(PAy*x123 - PGy[g]*x124 + x61);
                integral_21[g] += PAy*x125 - PGy[g]*x126 + tpi*(-2*x114 + 2*x115 + x16 + x65 + x66 - 3*x67 + 3*x70);
                integral_22[g] += PAy*x127 - PGy[g]*x128 + tpi*(x120 + x121 + 2*x77 + x79 + x80 - 2*x81);
                integral_23[g] += PAy*x129 - PGy[g]*x130 + x96;
                integral_24[g] += PBx*(PBx*x177 - PGx[g]*x178) - PGx[g]*(PBx*x178 - PGx[g]*(PAy*x132 + x153)) + tpi*(x146 - x176 + x177);
                integral_25[g] += PBx*(PAy*x149 + x161) - PGx[g]*(PAy*x152 + x163);
                integral_26[g] += PBx*(PAy*x157 + x165) - PGx[g]*(PAy*x159 + x167);
                integral_27[g] += PAy*x162 - PGy[g]*x164 + 2*x150;
                integral_28[g] += PAy*x166 - PGy[g]*x168 + x158;
                integral_29[g] += PAy*x169 - PGy[g]*x170;
                integral_30[g] += PBx*(PBx*x181 - PGx[g]*x183) - PGx[g]*(PBx*x183 - PGx[g]*x184) + x182;
                integral_31[g] += PBx*x185 - PGx[g]*x186;
                integral_32[g] += PBx*x188 - PGx[g]*x189;
                integral_33[g] += PBy*x185 - PGy[g]*x186 + x182;
                integral_34[g] += PBy*x188 - PGy[g]*x189;
                integral_35[g] += PAz*x169 - PGz[g]*x170 + tpi*(-2*x155 + 2*x156 + x16 + x89 + x90 - 3*x91 + 3*x94);

            } // for g
        } // if/else prim
    } // for prim
} // VRR_D_D

template void VRR_D_D<float>(
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

template void VRR_D_D<double>(
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

