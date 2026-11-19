#include "VRR_F_P.h"

template <typename T>
void VRR_F_P(
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

                const T x0 = PAx*scaled_boys_0[g];
                const T x1 = PBx*scaled_boys_0[g];
                const T x2 = PAx*scaled_boys_1[g];
                const T x3 = PBx*scaled_boys_1[g];
                const T x4 = PGx[g]*scaled_boys_1[g];
                const T x5 = PGx[g]*scaled_boys_2[g];
                const T x6 = 2*x5;
                const T x7 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x8 = -x4;
                const T x9 = PAx*(x1 + x8);
                const T x10 = -x5;
                const T x11 = x10 + x3;
                const T x12 = PGx[g]*x11;
                const T x13 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x14 = PAx*x11;
                const T x15 = PBx*scaled_boys_2[g];
                const T x16 = PGx[g]*scaled_boys_3[g];
                const T x17 = -x16;
                const T x18 = x15 + x17;
                const T x19 = PGx[g]*x18;
                const T x20 = x13 + x14 - x19;
                const T x21 = PAx*scaled_boys_2[g];
                const T x22 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x23 = x10 + x2;
                const T x24 = -3*x13 + 3*x7;
                const T x25 = PBy*scaled_boys_2[g];
                const T x26 = PGy[g]*scaled_boys_3[g];
                const T x27 = -x26;
                const T x28 = x25 + x27;
                const T x29 = PGx[g]*x28;
                const T x30 = PBy*scaled_boys_1[g];
                const T x31 = PGy[g]*scaled_boys_2[g];
                const T x32 = -x31;
                const T x33 = x30 + x32;
                const T x34 = PAx*x33;
                const T x35 = PBy*scaled_boys_0[g];
                const T x36 = PGy[g]*scaled_boys_1[g];
                const T x37 = -x36;
                const T x38 = x35 + x37;
                const T x39 = PAx*x38 - PGx[g]*x33;
                const T x40 = 2*tpi;
                const T x41 = -x30;
                const T x42 = -x29 + x34;
                const T x43 = -x25;
                const T x44 = -PGy[g]*scaled_boys_4[g];
                const T x45 = PBy*scaled_boys_3[g] + x44;
                const T x46 = PBz*scaled_boys_2[g];
                const T x47 = PGz[g]*scaled_boys_3[g];
                const T x48 = -x47;
                const T x49 = x46 + x48;
                const T x50 = PGx[g]*x49;
                const T x51 = PBz*scaled_boys_1[g];
                const T x52 = PGz[g]*scaled_boys_2[g];
                const T x53 = -x52;
                const T x54 = x51 + x53;
                const T x55 = PAx*x54;
                const T x56 = PBz*scaled_boys_0[g];
                const T x57 = PGz[g]*scaled_boys_1[g];
                const T x58 = -x57;
                const T x59 = x56 + x58;
                const T x60 = PAx*x59 - PGx[g]*x54;
                const T x61 = -x51;
                const T x62 = tpi*(x52 + x59 + x61);
                const T x63 = -x50 + x55;
                const T x64 = -x46;
                const T x65 = tpi*(x47 + x54 + x64);
                const T x66 = -PGz[g]*scaled_boys_4[g];
                const T x67 = PBz*scaled_boys_3[g] + x66;
                const T x68 = PAy*scaled_boys_0[g];
                const T x69 = x37 + x68;
                const T x70 = PBx*x69;
                const T x71 = PAy*scaled_boys_1[g];
                const T x72 = x32 + x71;
                const T x73 = PBx*x72;
                const T x74 = PGx[g]*x72;
                const T x75 = PAy*scaled_boys_2[g];
                const T x76 = x27 + x75;
                const T x77 = PGx[g]*x76;
                const T x78 = -x71;
                const T x79 = tpi*(x31 + x69 + x78);
                const T x80 = x73 - x77;
                const T x81 = -x75;
                const T x82 = tpi*(x26 + x72 + x81);
                const T x83 = PAy*scaled_boys_3[g] + x44;
                const T x84 = PGy[g]*x28;
                const T x85 = PAy*x33;
                const T x86 = -x13;
                const T x87 = PAy*x38;
                const T x88 = PGy[g]*x33;
                const T x89 = x7 + x87 - x88;
                const T x90 = x13 - x84 + x85;
                const T x91 = PAy*x28 - PGy[g]*x45 + x22;
                const T x92 = PGy[g]*x49;
                const T x93 = PAy*x54;
                const T x94 = PAy*x59 - PGy[g]*x54;
                const T x95 = tpi*(x92 - x93 + x94);
                const T x96 = -x92 + x93;
                const T x97 = PAy*x49 - PGy[g]*x67;
                const T x98 = PAz*scaled_boys_0[g];
                const T x99 = x58 + x98;
                const T x100 = PBx*x99;
                const T x101 = PAz*scaled_boys_1[g];
                const T x102 = x101 + x53;
                const T x103 = PBx*x102;
                const T x104 = PGx[g]*x102;
                const T x105 = PAz*scaled_boys_2[g];
                const T x106 = x105 + x48;
                const T x107 = PGx[g]*x106;
                const T x108 = -x101;
                const T x109 = tpi*(x108 + x52 + x99);
                const T x110 = x103 - x107;
                const T x111 = -x105;
                const T x112 = tpi*(x102 + x111 + x47);
                const T x113 = PAz*scaled_boys_3[g] + x66;
                const T x114 = PGy[g]*x106;
                const T x115 = PBy*x102;
                const T x116 = -x115;
                const T x117 = PBy*x99;
                const T x118 = PGy[g]*x102;
                const T x119 = -x118;
                const T x120 = x117 + x119;
                const T x121 = -x114;
                const T x122 = x115 + x121;
                const T x123 = -PGy[g]*x113;
                const T x124 = PBy*x106 + x123;
                const T x125 = PGz[g]*x49;
                const T x126 = PAz*x54;
                const T x127 = PAz*x59;
                const T x128 = PGz[g]*x54;
                const T x129 = x127 - x128 + x7;
                const T x130 = tpi*(x125 - x126 + x129 + x86);
                const T x131 = -x125 + x126 + x13;
                const T x132 = PAz*x49 - PGz[g]*x67 + x22;
                const T x133 = PAy*x69 - PGy[g]*x72;
                const T x134 = x133 + x7;
                const T x135 = PGy[g]*x76;
                const T x136 = PAy*x72;
                const T x137 = x135 - x136;
                const T x138 = x13 - x135 + x136;
                const T x139 = PAy*x76 - PGy[g]*x83 + x22;
                const T x140 = 2*x31;
                const T x141 = PAy*x89 - PGy[g]*x90 + tpi*(x140 + x35 - 2*x36 + x41 + x68 + x78);
                const T x142 = PAy*x90 - PGy[g]*x91 + tpi*(-x140 + 2*x26 + x30 + x43 + x71 + x81);
                const T x143 = PAy*x94 - PGy[g]*x96 + x62;
                const T x144 = PAy*x96 - PGy[g]*x97 + x65;
                const T x145 = PAy*x102;
                const T x146 = -x145;
                const T x147 = PAy*x99;
                const T x148 = x119 + x147;
                const T x149 = x121 + x145;
                const T x150 = PAy*x106 + x123;
                const T x151 = PAy*x120 - PGy[g]*x122 + x109;
                const T x152 = PAy*x122 - PGy[g]*x124 + x112;
                const T x153 = PAy*x129 - PGy[g]*x131;
                const T x154 = PAy*x131 - PGy[g]*x132;
                const T x155 = PAz*x99 - PGz[g]*x102;
                const T x156 = x155 + x7;
                const T x157 = PGz[g]*x106;
                const T x158 = PAz*x102;
                const T x159 = x157 - x158;
                const T x160 = tpi*(x156 + x159 + x86);
                const T x161 = x13 - x157 + x158;
                const T x162 = PAz*x106 - PGz[g]*x113 + x22;
                const T x163 = -PGy[g]*x161;
                const T x164 = PBy*x156 + x163;
                const T x165 = -PGy[g]*x162;
                const T x166 = PBy*x161 + x165;
                const T x167 = 2*x52;
                const T x168 = PAz*x129 - PGz[g]*x131 + tpi*(x108 + x167 + x56 - 2*x57 + x61 + x98);
                const T x169 = PAz*x131 - PGz[g]*x132 + tpi*(x101 + x111 - x167 + 2*x47 + x51 + x64);
                const T x170 = PAz*x156 - PGz[g]*x161 + 2*x109;
                const T x171 = PAz*x161 - PGz[g]*x162 + 2*x112;

                integral_0[g] = PAx*(PAx*(-x12 + x7 + x9) - PGx[g]*x20 + tpi*(x0 + x1 - x2 - x3 - 2*x4 + x6)) - PGx[g]*(PAx*x20 - PGx[g]*(PAx*x18 - PGx[g]*(PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x22) + tpi*(-x15 + 2*x16 + x2 - x21 + x3 - x6)) + tpi*(-PAx*x23 + PAx*(x0 + x8) - PGx[g]*x23 + PGx[g]*(x17 + x21) - 2*x12 - 2*x14 + 2*x19 + x24 + 2*x9);
                integral_1[g] = PAx*(PAx*x39 - PGx[g]*x42 + tpi*(x31 + x38 + x41)) - PGx[g]*(PAx*x42 - PGx[g]*(PAx*x28 - PGx[g]*x45) + tpi*(x26 + x33 + x43)) + x40*(x29 - x34 + x39);
                integral_2[g] = PAx*(PAx*x60 - PGx[g]*x63 + x62) - PGx[g]*(PAx*x63 - PGx[g]*(PAx*x49 - PGx[g]*x67) + x65) + x40*(x50 - x55 + x60);
                integral_3[g] = PAx*(PAx*(x70 - x74) - PGx[g]*x80 + x79) - PGx[g]*(PAx*x80 - PGx[g]*(PBx*x76 - PGx[g]*x83) + x82) + tpi*(PAx*x69 - PAx*x72 + x70 - x73 - 2*x74 + 2*x77);
                integral_4[g] = PAx*(PAx*x89 - PGx[g]*x90) - PGx[g]*(PAx*x90 - PGx[g]*x91) + tpi*(x84 - x85 + x86 + x89);
                integral_5[g] = PAx*(PAx*x94 - PGx[g]*x96) - PGx[g]*(PAx*x96 - PGx[g]*x97) + x95;
                integral_6[g] = PAx*(PAx*(x100 - x104) - PGx[g]*x110 + x109) - PGx[g]*(PAx*x110 - PGx[g]*(PBx*x106 - PGx[g]*x113) + x112) + tpi*(-PAx*x102 + PAx*x99 + x100 - x103 - 2*x104 + 2*x107);
                integral_7[g] = PAx*(PAx*x120 - PGx[g]*x122) - PGx[g]*(PAx*x122 - PGx[g]*x124) + tpi*(x114 + x116 + x120);
                integral_8[g] = PAx*(PAx*x129 - PGx[g]*x131) - PGx[g]*(PAx*x131 - PGx[g]*x132) + x130;
                integral_9[g] = PAx*(PBx*x134 - PGx[g]*x138) - PGx[g]*(PBx*x138 - PGx[g]*x139) + tpi*(x134 + x137 + x86);
                integral_10[g] = PAx*x141 - PGx[g]*x142;
                integral_11[g] = PAx*x143 - PGx[g]*x144;
                integral_12[g] = PAx*(PBx*x148 - PGx[g]*x149) - PGx[g]*(PBx*x149 - PGx[g]*x150) + tpi*(x114 + x146 + x148);
                integral_13[g] = PAx*x151 - PGx[g]*x152;
                integral_14[g] = PAx*x153 - PGx[g]*x154;
                integral_15[g] = PAx*(PBx*x156 - PGx[g]*x161) - PGx[g]*(PBx*x161 - PGx[g]*x162) + x160;
                integral_16[g] = PAx*x164 - PGx[g]*x166;
                integral_17[g] = PAx*x168 - PGx[g]*x169;
                integral_18[g] = PBx*(PAy*x134 - PGy[g]*x138 + 2*x79) - PGx[g]*(PAy*x138 - PGy[g]*x139 + 2*x82);
                integral_19[g] = PAy*x141 - PGy[g]*x142 + tpi*(x133 + x137 + x24 + 2*x84 - 2*x85 + 2*x87 - 2*x88);
                integral_20[g] = PAy*x143 - PGy[g]*x144 + 2*x95;
                integral_21[g] = PBx*(PAy*x148 - PGy[g]*x149 + x109) - PGx[g]*(PAy*x149 - PGy[g]*x150 + x112);
                integral_22[g] = PAy*x151 - PGy[g]*x152 + tpi*(2*x114 + x116 + x117 - 2*x118 + x146 + x147);
                integral_23[g] = PAy*x153 - PGy[g]*x154 + x130;
                integral_24[g] = PBx*(PAy*x156 + x163) - PGx[g]*(PAy*x161 + x165);
                integral_25[g] = PAy*x164 - PGy[g]*x166 + x160;
                integral_26[g] = PAy*x168 - PGy[g]*x169;
                integral_27[g] = PBx*x170 - PGx[g]*x171;
                integral_28[g] = PBy*x170 - PGy[g]*x171;
                integral_29[g] = PAz*x168 - PGz[g]*x169 + tpi*(2*x125 - 2*x126 + 2*x127 - 2*x128 + x155 + x159 + x24);

            } // for g
        } else {

            #pragma omp simd
            for (int g = 0; g < gb_size; ++g) {

                const T x0 = PAx*scaled_boys_0[g];
                const T x1 = PBx*scaled_boys_0[g];
                const T x2 = PAx*scaled_boys_1[g];
                const T x3 = PBx*scaled_boys_1[g];
                const T x4 = PGx[g]*scaled_boys_1[g];
                const T x5 = PGx[g]*scaled_boys_2[g];
                const T x6 = 2*x5;
                const T x7 = tpi*(scaled_boys_0[g] - scaled_boys_1[g]);
                const T x8 = -x4;
                const T x9 = PAx*(x1 + x8);
                const T x10 = -x5;
                const T x11 = x10 + x3;
                const T x12 = PGx[g]*x11;
                const T x13 = tpi*(scaled_boys_1[g] - scaled_boys_2[g]);
                const T x14 = PAx*x11;
                const T x15 = PBx*scaled_boys_2[g];
                const T x16 = PGx[g]*scaled_boys_3[g];
                const T x17 = -x16;
                const T x18 = x15 + x17;
                const T x19 = PGx[g]*x18;
                const T x20 = x13 + x14 - x19;
                const T x21 = PAx*scaled_boys_2[g];
                const T x22 = tpi*(scaled_boys_2[g] - scaled_boys_3[g]);
                const T x23 = x10 + x2;
                const T x24 = -3*x13 + 3*x7;
                const T x25 = PBy*scaled_boys_2[g];
                const T x26 = PGy[g]*scaled_boys_3[g];
                const T x27 = -x26;
                const T x28 = x25 + x27;
                const T x29 = PGx[g]*x28;
                const T x30 = PBy*scaled_boys_1[g];
                const T x31 = PGy[g]*scaled_boys_2[g];
                const T x32 = -x31;
                const T x33 = x30 + x32;
                const T x34 = PAx*x33;
                const T x35 = PBy*scaled_boys_0[g];
                const T x36 = PGy[g]*scaled_boys_1[g];
                const T x37 = -x36;
                const T x38 = x35 + x37;
                const T x39 = PAx*x38 - PGx[g]*x33;
                const T x40 = 2*tpi;
                const T x41 = -x30;
                const T x42 = -x29 + x34;
                const T x43 = -x25;
                const T x44 = -PGy[g]*scaled_boys_4[g];
                const T x45 = PBy*scaled_boys_3[g] + x44;
                const T x46 = PBz*scaled_boys_2[g];
                const T x47 = PGz[g]*scaled_boys_3[g];
                const T x48 = -x47;
                const T x49 = x46 + x48;
                const T x50 = PGx[g]*x49;
                const T x51 = PBz*scaled_boys_1[g];
                const T x52 = PGz[g]*scaled_boys_2[g];
                const T x53 = -x52;
                const T x54 = x51 + x53;
                const T x55 = PAx*x54;
                const T x56 = PBz*scaled_boys_0[g];
                const T x57 = PGz[g]*scaled_boys_1[g];
                const T x58 = -x57;
                const T x59 = x56 + x58;
                const T x60 = PAx*x59 - PGx[g]*x54;
                const T x61 = -x51;
                const T x62 = tpi*(x52 + x59 + x61);
                const T x63 = -x50 + x55;
                const T x64 = -x46;
                const T x65 = tpi*(x47 + x54 + x64);
                const T x66 = -PGz[g]*scaled_boys_4[g];
                const T x67 = PBz*scaled_boys_3[g] + x66;
                const T x68 = PAy*scaled_boys_0[g];
                const T x69 = x37 + x68;
                const T x70 = PBx*x69;
                const T x71 = PAy*scaled_boys_1[g];
                const T x72 = x32 + x71;
                const T x73 = PBx*x72;
                const T x74 = PGx[g]*x72;
                const T x75 = PAy*scaled_boys_2[g];
                const T x76 = x27 + x75;
                const T x77 = PGx[g]*x76;
                const T x78 = -x71;
                const T x79 = tpi*(x31 + x69 + x78);
                const T x80 = x73 - x77;
                const T x81 = -x75;
                const T x82 = tpi*(x26 + x72 + x81);
                const T x83 = PAy*scaled_boys_3[g] + x44;
                const T x84 = PGy[g]*x28;
                const T x85 = PAy*x33;
                const T x86 = -x13;
                const T x87 = PAy*x38;
                const T x88 = PGy[g]*x33;
                const T x89 = x7 + x87 - x88;
                const T x90 = x13 - x84 + x85;
                const T x91 = PAy*x28 - PGy[g]*x45 + x22;
                const T x92 = PGy[g]*x49;
                const T x93 = PAy*x54;
                const T x94 = PAy*x59 - PGy[g]*x54;
                const T x95 = tpi*(x92 - x93 + x94);
                const T x96 = -x92 + x93;
                const T x97 = PAy*x49 - PGy[g]*x67;
                const T x98 = PAz*scaled_boys_0[g];
                const T x99 = x58 + x98;
                const T x100 = PBx*x99;
                const T x101 = PAz*scaled_boys_1[g];
                const T x102 = x101 + x53;
                const T x103 = PBx*x102;
                const T x104 = PGx[g]*x102;
                const T x105 = PAz*scaled_boys_2[g];
                const T x106 = x105 + x48;
                const T x107 = PGx[g]*x106;
                const T x108 = -x101;
                const T x109 = tpi*(x108 + x52 + x99);
                const T x110 = x103 - x107;
                const T x111 = -x105;
                const T x112 = tpi*(x102 + x111 + x47);
                const T x113 = PAz*scaled_boys_3[g] + x66;
                const T x114 = PGy[g]*x106;
                const T x115 = PBy*x102;
                const T x116 = -x115;
                const T x117 = PBy*x99;
                const T x118 = PGy[g]*x102;
                const T x119 = -x118;
                const T x120 = x117 + x119;
                const T x121 = -x114;
                const T x122 = x115 + x121;
                const T x123 = -PGy[g]*x113;
                const T x124 = PBy*x106 + x123;
                const T x125 = PGz[g]*x49;
                const T x126 = PAz*x54;
                const T x127 = PAz*x59;
                const T x128 = PGz[g]*x54;
                const T x129 = x127 - x128 + x7;
                const T x130 = tpi*(x125 - x126 + x129 + x86);
                const T x131 = -x125 + x126 + x13;
                const T x132 = PAz*x49 - PGz[g]*x67 + x22;
                const T x133 = PAy*x69 - PGy[g]*x72;
                const T x134 = x133 + x7;
                const T x135 = PGy[g]*x76;
                const T x136 = PAy*x72;
                const T x137 = x135 - x136;
                const T x138 = x13 - x135 + x136;
                const T x139 = PAy*x76 - PGy[g]*x83 + x22;
                const T x140 = 2*x31;
                const T x141 = PAy*x89 - PGy[g]*x90 + tpi*(x140 + x35 - 2*x36 + x41 + x68 + x78);
                const T x142 = PAy*x90 - PGy[g]*x91 + tpi*(-x140 + 2*x26 + x30 + x43 + x71 + x81);
                const T x143 = PAy*x94 - PGy[g]*x96 + x62;
                const T x144 = PAy*x96 - PGy[g]*x97 + x65;
                const T x145 = PAy*x102;
                const T x146 = -x145;
                const T x147 = PAy*x99;
                const T x148 = x119 + x147;
                const T x149 = x121 + x145;
                const T x150 = PAy*x106 + x123;
                const T x151 = PAy*x120 - PGy[g]*x122 + x109;
                const T x152 = PAy*x122 - PGy[g]*x124 + x112;
                const T x153 = PAy*x129 - PGy[g]*x131;
                const T x154 = PAy*x131 - PGy[g]*x132;
                const T x155 = PAz*x99 - PGz[g]*x102;
                const T x156 = x155 + x7;
                const T x157 = PGz[g]*x106;
                const T x158 = PAz*x102;
                const T x159 = x157 - x158;
                const T x160 = tpi*(x156 + x159 + x86);
                const T x161 = x13 - x157 + x158;
                const T x162 = PAz*x106 - PGz[g]*x113 + x22;
                const T x163 = -PGy[g]*x161;
                const T x164 = PBy*x156 + x163;
                const T x165 = -PGy[g]*x162;
                const T x166 = PBy*x161 + x165;
                const T x167 = 2*x52;
                const T x168 = PAz*x129 - PGz[g]*x131 + tpi*(x108 + x167 + x56 - 2*x57 + x61 + x98);
                const T x169 = PAz*x131 - PGz[g]*x132 + tpi*(x101 + x111 - x167 + 2*x47 + x51 + x64);
                const T x170 = PAz*x156 - PGz[g]*x161 + 2*x109;
                const T x171 = PAz*x161 - PGz[g]*x162 + 2*x112;

                integral_0[g] += PAx*(PAx*(-x12 + x7 + x9) - PGx[g]*x20 + tpi*(x0 + x1 - x2 - x3 - 2*x4 + x6)) - PGx[g]*(PAx*x20 - PGx[g]*(PAx*x18 - PGx[g]*(PBx*scaled_boys_3[g] - PGx[g]*scaled_boys_4[g]) + x22) + tpi*(-x15 + 2*x16 + x2 - x21 + x3 - x6)) + tpi*(-PAx*x23 + PAx*(x0 + x8) - PGx[g]*x23 + PGx[g]*(x17 + x21) - 2*x12 - 2*x14 + 2*x19 + x24 + 2*x9);
                integral_1[g] += PAx*(PAx*x39 - PGx[g]*x42 + tpi*(x31 + x38 + x41)) - PGx[g]*(PAx*x42 - PGx[g]*(PAx*x28 - PGx[g]*x45) + tpi*(x26 + x33 + x43)) + x40*(x29 - x34 + x39);
                integral_2[g] += PAx*(PAx*x60 - PGx[g]*x63 + x62) - PGx[g]*(PAx*x63 - PGx[g]*(PAx*x49 - PGx[g]*x67) + x65) + x40*(x50 - x55 + x60);
                integral_3[g] += PAx*(PAx*(x70 - x74) - PGx[g]*x80 + x79) - PGx[g]*(PAx*x80 - PGx[g]*(PBx*x76 - PGx[g]*x83) + x82) + tpi*(PAx*x69 - PAx*x72 + x70 - x73 - 2*x74 + 2*x77);
                integral_4[g] += PAx*(PAx*x89 - PGx[g]*x90) - PGx[g]*(PAx*x90 - PGx[g]*x91) + tpi*(x84 - x85 + x86 + x89);
                integral_5[g] += PAx*(PAx*x94 - PGx[g]*x96) - PGx[g]*(PAx*x96 - PGx[g]*x97) + x95;
                integral_6[g] += PAx*(PAx*(x100 - x104) - PGx[g]*x110 + x109) - PGx[g]*(PAx*x110 - PGx[g]*(PBx*x106 - PGx[g]*x113) + x112) + tpi*(-PAx*x102 + PAx*x99 + x100 - x103 - 2*x104 + 2*x107);
                integral_7[g] += PAx*(PAx*x120 - PGx[g]*x122) - PGx[g]*(PAx*x122 - PGx[g]*x124) + tpi*(x114 + x116 + x120);
                integral_8[g] += PAx*(PAx*x129 - PGx[g]*x131) - PGx[g]*(PAx*x131 - PGx[g]*x132) + x130;
                integral_9[g] += PAx*(PBx*x134 - PGx[g]*x138) - PGx[g]*(PBx*x138 - PGx[g]*x139) + tpi*(x134 + x137 + x86);
                integral_10[g] += PAx*x141 - PGx[g]*x142;
                integral_11[g] += PAx*x143 - PGx[g]*x144;
                integral_12[g] += PAx*(PBx*x148 - PGx[g]*x149) - PGx[g]*(PBx*x149 - PGx[g]*x150) + tpi*(x114 + x146 + x148);
                integral_13[g] += PAx*x151 - PGx[g]*x152;
                integral_14[g] += PAx*x153 - PGx[g]*x154;
                integral_15[g] += PAx*(PBx*x156 - PGx[g]*x161) - PGx[g]*(PBx*x161 - PGx[g]*x162) + x160;
                integral_16[g] += PAx*x164 - PGx[g]*x166;
                integral_17[g] += PAx*x168 - PGx[g]*x169;
                integral_18[g] += PBx*(PAy*x134 - PGy[g]*x138 + 2*x79) - PGx[g]*(PAy*x138 - PGy[g]*x139 + 2*x82);
                integral_19[g] += PAy*x141 - PGy[g]*x142 + tpi*(x133 + x137 + x24 + 2*x84 - 2*x85 + 2*x87 - 2*x88);
                integral_20[g] += PAy*x143 - PGy[g]*x144 + 2*x95;
                integral_21[g] += PBx*(PAy*x148 - PGy[g]*x149 + x109) - PGx[g]*(PAy*x149 - PGy[g]*x150 + x112);
                integral_22[g] += PAy*x151 - PGy[g]*x152 + tpi*(2*x114 + x116 + x117 - 2*x118 + x146 + x147);
                integral_23[g] += PAy*x153 - PGy[g]*x154 + x130;
                integral_24[g] += PBx*(PAy*x156 + x163) - PGx[g]*(PAy*x161 + x165);
                integral_25[g] += PAy*x164 - PGy[g]*x166 + x160;
                integral_26[g] += PAy*x168 - PGy[g]*x169;
                integral_27[g] += PBx*x170 - PGx[g]*x171;
                integral_28[g] += PBy*x170 - PGy[g]*x171;
                integral_29[g] += PAz*x168 - PGz[g]*x169 + tpi*(2*x125 - 2*x126 + 2*x127 - 2*x128 + x155 + x159 + x24);

            } // for g
        } // if/else prim
    } // for prim
} // VRR_F_P

template void VRR_F_P<float>(
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

template void VRR_F_P<double>(
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

