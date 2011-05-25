#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00pp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|pp) integrals */

void d12hrr_order_00pp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][2][11] = int_stack + 0;
 Libderiv->deriv_classes[0][2][10] = int_stack + 6;
 Libderiv->deriv_classes[0][2][9] = int_stack + 12;
 Libderiv->deriv_classes[0][2][8] = int_stack + 18;
 Libderiv->deriv_classes[0][2][7] = int_stack + 24;
 Libderiv->dvrr_classes[0][1] = int_stack + 30;
 Libderiv->deriv_classes[0][2][6] = int_stack + 33;
 Libderiv->deriv_classes[0][2][2] = int_stack + 39;
 Libderiv->deriv_classes[0][2][1] = int_stack + 45;
 Libderiv->deriv_classes[0][2][0] = int_stack + 51;
 Libderiv->deriv2_classes[0][1][143] = int_stack + 57;
 Libderiv->deriv2_classes[0][2][143] = int_stack + 60;
 Libderiv->deriv2_classes[0][1][131] = int_stack + 66;
 Libderiv->deriv2_classes[0][2][131] = int_stack + 69;
 Libderiv->deriv2_classes[0][1][130] = int_stack + 75;
 Libderiv->deriv2_classes[0][2][130] = int_stack + 78;
 Libderiv->deriv2_classes[0][1][119] = int_stack + 84;
 Libderiv->deriv2_classes[0][2][119] = int_stack + 87;
 Libderiv->deriv2_classes[0][1][118] = int_stack + 93;
 Libderiv->deriv2_classes[0][2][118] = int_stack + 96;
 Libderiv->deriv2_classes[0][1][117] = int_stack + 102;
 Libderiv->deriv2_classes[0][2][117] = int_stack + 105;
 Libderiv->deriv2_classes[0][1][107] = int_stack + 111;
 Libderiv->deriv2_classes[0][2][107] = int_stack + 114;
 Libderiv->deriv2_classes[0][1][106] = int_stack + 120;
 Libderiv->deriv2_classes[0][2][106] = int_stack + 123;
 Libderiv->deriv2_classes[0][1][105] = int_stack + 129;
 Libderiv->deriv2_classes[0][2][105] = int_stack + 132;
 Libderiv->deriv2_classes[0][1][104] = int_stack + 138;
 Libderiv->deriv2_classes[0][2][104] = int_stack + 141;
 Libderiv->deriv2_classes[0][1][95] = int_stack + 147;
 Libderiv->deriv2_classes[0][2][95] = int_stack + 150;
 Libderiv->deriv2_classes[0][1][94] = int_stack + 156;
 Libderiv->deriv2_classes[0][2][94] = int_stack + 159;
 Libderiv->deriv2_classes[0][1][93] = int_stack + 165;
 Libderiv->deriv2_classes[0][2][93] = int_stack + 168;
 Libderiv->deriv2_classes[0][1][92] = int_stack + 174;
 Libderiv->deriv2_classes[0][2][92] = int_stack + 177;
 Libderiv->deriv2_classes[0][1][91] = int_stack + 183;
 Libderiv->deriv2_classes[0][2][91] = int_stack + 186;
 Libderiv->deriv_classes[0][1][11] = int_stack + 192;
 Libderiv->deriv2_classes[0][1][83] = int_stack + 195;
 Libderiv->deriv2_classes[0][2][83] = int_stack + 198;
 Libderiv->deriv_classes[0][1][10] = int_stack + 204;
 Libderiv->deriv2_classes[0][1][82] = int_stack + 207;
 Libderiv->deriv2_classes[0][2][82] = int_stack + 210;
 Libderiv->deriv_classes[0][1][9] = int_stack + 216;
 Libderiv->deriv2_classes[0][1][81] = int_stack + 219;
 Libderiv->deriv2_classes[0][2][81] = int_stack + 222;
 Libderiv->deriv_classes[0][1][8] = int_stack + 228;
 Libderiv->deriv2_classes[0][1][80] = int_stack + 231;
 Libderiv->deriv2_classes[0][2][80] = int_stack + 234;
 Libderiv->deriv_classes[0][1][7] = int_stack + 240;
 Libderiv->deriv2_classes[0][1][79] = int_stack + 243;
 Libderiv->deriv2_classes[0][2][79] = int_stack + 246;
 Libderiv->deriv_classes[0][1][6] = int_stack + 252;
 Libderiv->deriv2_classes[0][1][78] = int_stack + 255;
 Libderiv->deriv2_classes[0][2][78] = int_stack + 258;
 Libderiv->deriv2_classes[0][1][35] = int_stack + 264;
 Libderiv->deriv2_classes[0][2][35] = int_stack + 267;
 Libderiv->deriv2_classes[0][1][34] = int_stack + 273;
 Libderiv->deriv2_classes[0][2][34] = int_stack + 276;
 Libderiv->deriv2_classes[0][1][33] = int_stack + 282;
 Libderiv->deriv2_classes[0][2][33] = int_stack + 285;
 Libderiv->deriv2_classes[0][1][32] = int_stack + 291;
 Libderiv->deriv2_classes[0][2][32] = int_stack + 294;
 Libderiv->deriv2_classes[0][1][31] = int_stack + 300;
 Libderiv->deriv2_classes[0][2][31] = int_stack + 303;
 Libderiv->deriv_classes[0][1][2] = int_stack + 309;
 Libderiv->deriv2_classes[0][1][30] = int_stack + 312;
 Libderiv->deriv2_classes[0][2][30] = int_stack + 315;
 Libderiv->deriv2_classes[0][1][26] = int_stack + 321;
 Libderiv->deriv2_classes[0][2][26] = int_stack + 324;
 Libderiv->deriv2_classes[0][1][23] = int_stack + 330;
 Libderiv->deriv2_classes[0][2][23] = int_stack + 333;
 Libderiv->deriv2_classes[0][1][22] = int_stack + 339;
 Libderiv->deriv2_classes[0][2][22] = int_stack + 342;
 Libderiv->deriv2_classes[0][1][21] = int_stack + 348;
 Libderiv->deriv2_classes[0][2][21] = int_stack + 351;
 Libderiv->deriv2_classes[0][1][20] = int_stack + 357;
 Libderiv->deriv2_classes[0][2][20] = int_stack + 360;
 Libderiv->deriv2_classes[0][1][19] = int_stack + 366;
 Libderiv->deriv2_classes[0][2][19] = int_stack + 369;
 Libderiv->deriv_classes[0][1][1] = int_stack + 375;
 Libderiv->deriv2_classes[0][1][18] = int_stack + 378;
 Libderiv->deriv2_classes[0][2][18] = int_stack + 381;
 Libderiv->deriv2_classes[0][1][14] = int_stack + 387;
 Libderiv->deriv2_classes[0][2][14] = int_stack + 390;
 Libderiv->deriv2_classes[0][1][13] = int_stack + 396;
 Libderiv->deriv2_classes[0][2][13] = int_stack + 399;
 Libderiv->deriv2_classes[0][1][11] = int_stack + 405;
 Libderiv->deriv2_classes[0][2][11] = int_stack + 408;
 Libderiv->deriv2_classes[0][1][10] = int_stack + 414;
 Libderiv->deriv2_classes[0][2][10] = int_stack + 417;
 Libderiv->deriv2_classes[0][1][9] = int_stack + 423;
 Libderiv->deriv2_classes[0][2][9] = int_stack + 426;
 Libderiv->deriv2_classes[0][1][8] = int_stack + 432;
 Libderiv->deriv2_classes[0][2][8] = int_stack + 435;
 Libderiv->deriv2_classes[0][1][7] = int_stack + 441;
 Libderiv->deriv2_classes[0][2][7] = int_stack + 444;
 Libderiv->deriv_classes[0][1][0] = int_stack + 450;
 Libderiv->deriv2_classes[0][1][6] = int_stack + 453;
 Libderiv->deriv2_classes[0][2][6] = int_stack + 456;
 Libderiv->deriv2_classes[0][1][2] = int_stack + 462;
 Libderiv->deriv2_classes[0][2][2] = int_stack + 465;
 Libderiv->deriv2_classes[0][1][1] = int_stack + 471;
 Libderiv->deriv2_classes[0][2][1] = int_stack + 474;
 Libderiv->deriv2_classes[0][1][0] = int_stack + 480;
 Libderiv->deriv2_classes[0][2][0] = int_stack + 483;
 memset(int_stack,0,3912);

 Libderiv->dvrr_stack = int_stack + 525;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00pp(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+489,int_stack+0,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30,1);
     Libderiv->ABCD[11] = int_stack + 489;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+498,int_stack+6,int_stack+204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 498;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+0,int_stack+12,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+9,int_stack+18,int_stack+228, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 9;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+507,int_stack+24,int_stack+240, 0.0, zero_stack, 1.0, int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 507;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+18,int_stack+33,int_stack+252, 1.0, int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 18;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+27,int_stack+39,int_stack+309,1);
     Libderiv->ABCD[2] = int_stack + 27;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+36,int_stack+45,int_stack+375,1);
     Libderiv->ABCD[1] = int_stack + 36;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+516,int_stack+51,int_stack+450,1);
     Libderiv->ABCD[0] = int_stack + 516;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+45,int_stack+60,int_stack+57, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+192,1);
     Libderiv->ABCD[155] = int_stack + 45;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+54,int_stack+69,int_stack+66, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192, 1.0, int_stack+204,1);
     Libderiv->ABCD[143] = int_stack + 54;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+63,int_stack+78,int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+204, 0.0, zero_stack,1);
     Libderiv->ABCD[142] = int_stack + 63;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+72,int_stack+87,int_stack+84, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192, 0.0, zero_stack, 1.0, int_stack+216,1);
     Libderiv->ABCD[131] = int_stack + 72;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+81,int_stack+96,int_stack+93, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+204, 1.0, int_stack+216, 0.0, zero_stack,1);
     Libderiv->ABCD[130] = int_stack + 81;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+90,int_stack+105,int_stack+102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[129] = int_stack + 90;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+99,int_stack+114,int_stack+111, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+228,1);
     Libderiv->ABCD[119] = int_stack + 99;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+108,int_stack+123,int_stack+120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+204, 0.0, zero_stack, 1.0, int_stack+228, 0.0, zero_stack,1);
     Libderiv->ABCD[118] = int_stack + 108;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+117,int_stack+132,int_stack+129, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 1.0, int_stack+228, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[117] = int_stack + 117;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+126,int_stack+141,int_stack+138, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[116] = int_stack + 126;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+135,int_stack+150,int_stack+147, 0.0, zero_stack, 1.0, int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240,1);
     Libderiv->ABCD[107] = int_stack + 135;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+144,int_stack+159,int_stack+156, 0.0, zero_stack, 1.0, int_stack+204, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240, 0.0, zero_stack,1);
     Libderiv->ABCD[106] = int_stack + 144;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+153,int_stack+168,int_stack+165, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 1.0, int_stack+240, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[105] = int_stack + 153;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+162,int_stack+177,int_stack+174, 0.0, zero_stack, 1.0, int_stack+228, 1.0, int_stack+240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[104] = int_stack + 162;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+171,int_stack+186,int_stack+183, 0.0, zero_stack, 2.0, int_stack+240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[103] = int_stack + 171;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+180,int_stack+198,int_stack+195, 1.0, int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+252,1);
     Libderiv->ABCD[95] = int_stack + 180;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+189,int_stack+210,int_stack+207, 1.0, int_stack+204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+252, 0.0, zero_stack,1);
     Libderiv->ABCD[94] = int_stack + 189;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+198,int_stack+222,int_stack+219, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+252, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[93] = int_stack + 198;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+207,int_stack+234,int_stack+231, 1.0, int_stack+228, 0.0, zero_stack, 1.0, int_stack+252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[92] = int_stack + 207;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+216,int_stack+246,int_stack+243, 1.0, int_stack+240, 1.0, int_stack+252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[91] = int_stack + 216;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+225,int_stack+258,int_stack+255, 2.0, int_stack+252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[90] = int_stack + 225;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+234,int_stack+267,int_stack+264, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309,1);
     Libderiv->ABCD[47] = int_stack + 234;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+243,int_stack+276,int_stack+273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309, 0.0, zero_stack,1);
     Libderiv->ABCD[46] = int_stack + 243;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+252,int_stack+285,int_stack+282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[45] = int_stack + 252;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+261,int_stack+294,int_stack+291, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[44] = int_stack + 261;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+270,int_stack+303,int_stack+300, 0.0, zero_stack, 1.0, int_stack+309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[43] = int_stack + 270;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+279,int_stack+315,int_stack+312, 1.0, int_stack+309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[42] = int_stack + 279;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+288,int_stack+324,int_stack+321,1);
     Libderiv->ABCD[38] = int_stack + 288;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+297,int_stack+333,int_stack+330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375,1);
     Libderiv->ABCD[35] = int_stack + 297;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+306,int_stack+342,int_stack+339, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack,1);
     Libderiv->ABCD[34] = int_stack + 306;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+315,int_stack+351,int_stack+348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[33] = int_stack + 315;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+324,int_stack+360,int_stack+357, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[32] = int_stack + 324;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+333,int_stack+369,int_stack+366, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[31] = int_stack + 333;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+342,int_stack+381,int_stack+378, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[30] = int_stack + 342;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+351,int_stack+390,int_stack+387,1);
     Libderiv->ABCD[26] = int_stack + 351;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+360,int_stack+399,int_stack+396,1);
     Libderiv->ABCD[25] = int_stack + 360;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+369,int_stack+408,int_stack+405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,1);
     Libderiv->ABCD[23] = int_stack + 369;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+378,int_stack+417,int_stack+414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,1);
     Libderiv->ABCD[22] = int_stack + 378;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+387,int_stack+426,int_stack+423, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[21] = int_stack + 387;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+396,int_stack+435,int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[20] = int_stack + 396;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+405,int_stack+444,int_stack+441, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[19] = int_stack + 405;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+414,int_stack+456,int_stack+453, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[18] = int_stack + 414;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+423,int_stack+465,int_stack+462,1);
     Libderiv->ABCD[14] = int_stack + 423;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+432,int_stack+474,int_stack+471,1);
     Libderiv->ABCD[13] = int_stack + 432;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+441,int_stack+483,int_stack+480,1);
     Libderiv->ABCD[12] = int_stack + 441;

}
