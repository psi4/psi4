#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|fp) integrals */

void d12hrr_order_00fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][4][11] = int_stack + 0;
 Libderiv->deriv_classes[0][4][10] = int_stack + 15;
 Libderiv->deriv_classes[0][4][9] = int_stack + 30;
 Libderiv->deriv_classes[0][4][8] = int_stack + 45;
 Libderiv->deriv_classes[0][4][7] = int_stack + 60;
 Libderiv->dvrr_classes[0][3] = int_stack + 75;
 Libderiv->deriv_classes[0][4][6] = int_stack + 85;
 Libderiv->deriv_classes[0][4][2] = int_stack + 100;
 Libderiv->deriv_classes[0][4][1] = int_stack + 115;
 Libderiv->deriv_classes[0][4][0] = int_stack + 130;
 Libderiv->deriv2_classes[0][3][143] = int_stack + 145;
 Libderiv->deriv2_classes[0][4][143] = int_stack + 155;
 Libderiv->deriv2_classes[0][3][131] = int_stack + 170;
 Libderiv->deriv2_classes[0][4][131] = int_stack + 180;
 Libderiv->deriv2_classes[0][3][130] = int_stack + 195;
 Libderiv->deriv2_classes[0][4][130] = int_stack + 205;
 Libderiv->deriv2_classes[0][3][119] = int_stack + 220;
 Libderiv->deriv2_classes[0][4][119] = int_stack + 230;
 Libderiv->deriv2_classes[0][3][118] = int_stack + 245;
 Libderiv->deriv2_classes[0][4][118] = int_stack + 255;
 Libderiv->deriv2_classes[0][3][117] = int_stack + 270;
 Libderiv->deriv2_classes[0][4][117] = int_stack + 280;
 Libderiv->deriv2_classes[0][3][107] = int_stack + 295;
 Libderiv->deriv2_classes[0][4][107] = int_stack + 305;
 Libderiv->deriv2_classes[0][3][106] = int_stack + 320;
 Libderiv->deriv2_classes[0][4][106] = int_stack + 330;
 Libderiv->deriv2_classes[0][3][105] = int_stack + 345;
 Libderiv->deriv2_classes[0][4][105] = int_stack + 355;
 Libderiv->deriv2_classes[0][3][104] = int_stack + 370;
 Libderiv->deriv2_classes[0][4][104] = int_stack + 380;
 Libderiv->deriv2_classes[0][3][95] = int_stack + 395;
 Libderiv->deriv2_classes[0][4][95] = int_stack + 405;
 Libderiv->deriv2_classes[0][3][94] = int_stack + 420;
 Libderiv->deriv2_classes[0][4][94] = int_stack + 430;
 Libderiv->deriv2_classes[0][3][93] = int_stack + 445;
 Libderiv->deriv2_classes[0][4][93] = int_stack + 455;
 Libderiv->deriv2_classes[0][3][92] = int_stack + 470;
 Libderiv->deriv2_classes[0][4][92] = int_stack + 480;
 Libderiv->deriv2_classes[0][3][91] = int_stack + 495;
 Libderiv->deriv2_classes[0][4][91] = int_stack + 505;
 Libderiv->deriv_classes[0][3][11] = int_stack + 520;
 Libderiv->deriv2_classes[0][3][83] = int_stack + 530;
 Libderiv->deriv2_classes[0][4][83] = int_stack + 540;
 Libderiv->deriv_classes[0][3][10] = int_stack + 555;
 Libderiv->deriv2_classes[0][3][82] = int_stack + 565;
 Libderiv->deriv2_classes[0][4][82] = int_stack + 575;
 Libderiv->deriv_classes[0][3][9] = int_stack + 590;
 Libderiv->deriv2_classes[0][3][81] = int_stack + 600;
 Libderiv->deriv2_classes[0][4][81] = int_stack + 610;
 Libderiv->deriv_classes[0][3][8] = int_stack + 625;
 Libderiv->deriv2_classes[0][3][80] = int_stack + 635;
 Libderiv->deriv2_classes[0][4][80] = int_stack + 645;
 Libderiv->deriv_classes[0][3][7] = int_stack + 660;
 Libderiv->deriv2_classes[0][3][79] = int_stack + 670;
 Libderiv->deriv2_classes[0][4][79] = int_stack + 680;
 Libderiv->deriv_classes[0][3][6] = int_stack + 695;
 Libderiv->deriv2_classes[0][3][78] = int_stack + 705;
 Libderiv->deriv2_classes[0][4][78] = int_stack + 715;
 Libderiv->deriv2_classes[0][3][35] = int_stack + 730;
 Libderiv->deriv2_classes[0][4][35] = int_stack + 740;
 Libderiv->deriv2_classes[0][3][34] = int_stack + 755;
 Libderiv->deriv2_classes[0][4][34] = int_stack + 765;
 Libderiv->deriv2_classes[0][3][33] = int_stack + 780;
 Libderiv->deriv2_classes[0][4][33] = int_stack + 790;
 Libderiv->deriv2_classes[0][3][32] = int_stack + 805;
 Libderiv->deriv2_classes[0][4][32] = int_stack + 815;
 Libderiv->deriv2_classes[0][3][31] = int_stack + 830;
 Libderiv->deriv2_classes[0][4][31] = int_stack + 840;
 Libderiv->deriv_classes[0][3][2] = int_stack + 855;
 Libderiv->deriv2_classes[0][3][30] = int_stack + 865;
 Libderiv->deriv2_classes[0][4][30] = int_stack + 875;
 Libderiv->deriv2_classes[0][3][26] = int_stack + 890;
 Libderiv->deriv2_classes[0][4][26] = int_stack + 900;
 Libderiv->deriv2_classes[0][3][23] = int_stack + 915;
 Libderiv->deriv2_classes[0][4][23] = int_stack + 925;
 Libderiv->deriv2_classes[0][3][22] = int_stack + 940;
 Libderiv->deriv2_classes[0][4][22] = int_stack + 950;
 Libderiv->deriv2_classes[0][3][21] = int_stack + 965;
 Libderiv->deriv2_classes[0][4][21] = int_stack + 975;
 Libderiv->deriv2_classes[0][3][20] = int_stack + 990;
 Libderiv->deriv2_classes[0][4][20] = int_stack + 1000;
 Libderiv->deriv2_classes[0][3][19] = int_stack + 1015;
 Libderiv->deriv2_classes[0][4][19] = int_stack + 1025;
 Libderiv->deriv_classes[0][3][1] = int_stack + 1040;
 Libderiv->deriv2_classes[0][3][18] = int_stack + 1050;
 Libderiv->deriv2_classes[0][4][18] = int_stack + 1060;
 Libderiv->deriv2_classes[0][3][14] = int_stack + 1075;
 Libderiv->deriv2_classes[0][4][14] = int_stack + 1085;
 Libderiv->deriv2_classes[0][3][13] = int_stack + 1100;
 Libderiv->deriv2_classes[0][4][13] = int_stack + 1110;
 Libderiv->deriv2_classes[0][3][11] = int_stack + 1125;
 Libderiv->deriv2_classes[0][4][11] = int_stack + 1135;
 Libderiv->deriv2_classes[0][3][10] = int_stack + 1150;
 Libderiv->deriv2_classes[0][4][10] = int_stack + 1160;
 Libderiv->deriv2_classes[0][3][9] = int_stack + 1175;
 Libderiv->deriv2_classes[0][4][9] = int_stack + 1185;
 Libderiv->deriv2_classes[0][3][8] = int_stack + 1200;
 Libderiv->deriv2_classes[0][4][8] = int_stack + 1210;
 Libderiv->deriv2_classes[0][3][7] = int_stack + 1225;
 Libderiv->deriv2_classes[0][4][7] = int_stack + 1235;
 Libderiv->deriv_classes[0][3][0] = int_stack + 1250;
 Libderiv->deriv2_classes[0][3][6] = int_stack + 1260;
 Libderiv->deriv2_classes[0][4][6] = int_stack + 1270;
 Libderiv->deriv2_classes[0][3][2] = int_stack + 1285;
 Libderiv->deriv2_classes[0][4][2] = int_stack + 1295;
 Libderiv->deriv2_classes[0][3][1] = int_stack + 1310;
 Libderiv->deriv2_classes[0][4][1] = int_stack + 1320;
 Libderiv->deriv2_classes[0][3][0] = int_stack + 1335;
 Libderiv->deriv2_classes[0][4][0] = int_stack + 1345;
 memset(int_stack,0,10880);

 Libderiv->dvrr_stack = int_stack + 1660;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1360,int_stack+0,int_stack+520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75,1);
     Libderiv->ABCD[11] = int_stack + 1360;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1390,int_stack+15,int_stack+555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 1390;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+30,int_stack+590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1420,int_stack+45,int_stack+625, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 1420;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30,int_stack+60,int_stack+660, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 30;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1450,int_stack+85,int_stack+695, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 1450;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+60,int_stack+100,int_stack+855,1);
     Libderiv->ABCD[2] = int_stack + 60;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1480,int_stack+115,int_stack+1040,1);
     Libderiv->ABCD[1] = int_stack + 1480;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+130,int_stack+1250,1);
     Libderiv->ABCD[0] = int_stack + 90;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1510,int_stack+155,int_stack+145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+520,1);
     Libderiv->ABCD[155] = int_stack + 1510;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+120,int_stack+180,int_stack+170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+520, 1.0, int_stack+555,1);
     Libderiv->ABCD[143] = int_stack + 120;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+150,int_stack+205,int_stack+195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+555, 0.0, zero_stack,1);
     Libderiv->ABCD[142] = int_stack + 150;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+230,int_stack+220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+520, 0.0, zero_stack, 1.0, int_stack+590,1);
     Libderiv->ABCD[131] = int_stack + 180;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+210,int_stack+255,int_stack+245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+555, 1.0, int_stack+590, 0.0, zero_stack,1);
     Libderiv->ABCD[130] = int_stack + 210;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+240,int_stack+280,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+590, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[129] = int_stack + 240;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1540,int_stack+305,int_stack+295, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+625,1);
     Libderiv->ABCD[119] = int_stack + 1540;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+330,int_stack+320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+555, 0.0, zero_stack, 1.0, int_stack+625, 0.0, zero_stack,1);
     Libderiv->ABCD[118] = int_stack + 270;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+355,int_stack+345, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+590, 1.0, int_stack+625, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[117] = int_stack + 300;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+330,int_stack+380,int_stack+370, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[116] = int_stack + 330;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+405,int_stack+395, 0.0, zero_stack, 1.0, int_stack+520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+660,1);
     Libderiv->ABCD[107] = int_stack + 360;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+390,int_stack+430,int_stack+420, 0.0, zero_stack, 1.0, int_stack+555, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+660, 0.0, zero_stack,1);
     Libderiv->ABCD[106] = int_stack + 390;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1570,int_stack+455,int_stack+445, 0.0, zero_stack, 1.0, int_stack+590, 0.0, zero_stack, 1.0, int_stack+660, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[105] = int_stack + 1570;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+420,int_stack+480,int_stack+470, 0.0, zero_stack, 1.0, int_stack+625, 1.0, int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[104] = int_stack + 420;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+505,int_stack+495, 0.0, zero_stack, 2.0, int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[103] = int_stack + 450;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+480,int_stack+540,int_stack+530, 1.0, int_stack+520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+695,1);
     Libderiv->ABCD[95] = int_stack + 480;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+510,int_stack+575,int_stack+565, 1.0, int_stack+555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+695, 0.0, zero_stack,1);
     Libderiv->ABCD[94] = int_stack + 510;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+610,int_stack+600, 1.0, int_stack+590, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+695, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[93] = int_stack + 540;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+570,int_stack+645,int_stack+635, 1.0, int_stack+625, 0.0, zero_stack, 1.0, int_stack+695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[92] = int_stack + 570;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+680,int_stack+670, 1.0, int_stack+660, 1.0, int_stack+695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[91] = int_stack + 600;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+715,int_stack+705, 2.0, int_stack+695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[90] = int_stack + 630;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+660,int_stack+740,int_stack+730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+855,1);
     Libderiv->ABCD[47] = int_stack + 660;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+690,int_stack+765,int_stack+755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+855, 0.0, zero_stack,1);
     Libderiv->ABCD[46] = int_stack + 690;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+790,int_stack+780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+855, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[45] = int_stack + 720;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+750,int_stack+815,int_stack+805, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[44] = int_stack + 750;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+780,int_stack+840,int_stack+830, 0.0, zero_stack, 1.0, int_stack+855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[43] = int_stack + 780;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+810,int_stack+875,int_stack+865, 1.0, int_stack+855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[42] = int_stack + 810;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+840,int_stack+900,int_stack+890,1);
     Libderiv->ABCD[38] = int_stack + 840;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+870,int_stack+925,int_stack+915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1040,1);
     Libderiv->ABCD[35] = int_stack + 870;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+950,int_stack+940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1040, 0.0, zero_stack,1);
     Libderiv->ABCD[34] = int_stack + 900;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+930,int_stack+975,int_stack+965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1040, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[33] = int_stack + 930;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+960,int_stack+1000,int_stack+990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[32] = int_stack + 960;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1600,int_stack+1025,int_stack+1015, 0.0, zero_stack, 1.0, int_stack+1040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[31] = int_stack + 1600;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+990,int_stack+1060,int_stack+1050, 1.0, int_stack+1040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[30] = int_stack + 990;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1020,int_stack+1085,int_stack+1075,1);
     Libderiv->ABCD[26] = int_stack + 1020;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1050,int_stack+1110,int_stack+1100,1);
     Libderiv->ABCD[25] = int_stack + 1050;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1135,int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250,1);
     Libderiv->ABCD[23] = int_stack + 1080;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1110,int_stack+1160,int_stack+1150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack,1);
     Libderiv->ABCD[22] = int_stack + 1110;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1140,int_stack+1185,int_stack+1175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[21] = int_stack + 1140;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1170,int_stack+1210,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[20] = int_stack + 1170;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1630,int_stack+1235,int_stack+1225, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[19] = int_stack + 1630;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+1270,int_stack+1260, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[18] = int_stack + 1200;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1230,int_stack+1295,int_stack+1285,1);
     Libderiv->ABCD[14] = int_stack + 1230;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1260,int_stack+1320,int_stack+1310,1);
     Libderiv->ABCD[13] = int_stack + 1260;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1290,int_stack+1345,int_stack+1335,1);
     Libderiv->ABCD[12] = int_stack + 1290;

}
