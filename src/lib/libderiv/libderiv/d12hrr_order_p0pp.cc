#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0pp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|pp) integrals */

void d12hrr_order_p0pp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[1][2][10] = int_stack + 18;
 Libderiv->deriv_classes[1][2][9] = int_stack + 36;
 Libderiv->deriv_classes[1][2][8] = int_stack + 54;
 Libderiv->deriv_classes[1][2][7] = int_stack + 72;
 Libderiv->dvrr_classes[1][1] = int_stack + 90;
 Libderiv->deriv_classes[1][2][6] = int_stack + 99;
 Libderiv->deriv_classes[1][2][2] = int_stack + 117;
 Libderiv->deriv_classes[1][2][1] = int_stack + 135;
 Libderiv->deriv_classes[1][2][0] = int_stack + 153;
 Libderiv->deriv2_classes[1][1][143] = int_stack + 171;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 180;
 Libderiv->deriv2_classes[1][1][131] = int_stack + 198;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 207;
 Libderiv->deriv2_classes[1][1][130] = int_stack + 225;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 234;
 Libderiv->deriv2_classes[1][1][119] = int_stack + 252;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 261;
 Libderiv->deriv2_classes[1][1][118] = int_stack + 279;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 288;
 Libderiv->deriv2_classes[1][1][117] = int_stack + 306;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 315;
 Libderiv->deriv2_classes[1][1][107] = int_stack + 333;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 342;
 Libderiv->deriv2_classes[1][1][106] = int_stack + 360;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 369;
 Libderiv->deriv2_classes[1][1][105] = int_stack + 387;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 396;
 Libderiv->deriv2_classes[1][1][104] = int_stack + 414;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 423;
 Libderiv->deriv2_classes[1][1][95] = int_stack + 441;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 450;
 Libderiv->deriv2_classes[1][1][94] = int_stack + 468;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 477;
 Libderiv->deriv2_classes[1][1][93] = int_stack + 495;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 504;
 Libderiv->deriv2_classes[1][1][92] = int_stack + 522;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 531;
 Libderiv->deriv2_classes[1][1][91] = int_stack + 549;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 558;
 Libderiv->deriv_classes[1][1][11] = int_stack + 576;
 Libderiv->deriv2_classes[1][1][83] = int_stack + 585;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 594;
 Libderiv->deriv_classes[1][1][10] = int_stack + 612;
 Libderiv->deriv2_classes[1][1][82] = int_stack + 621;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 630;
 Libderiv->deriv_classes[1][1][9] = int_stack + 648;
 Libderiv->deriv2_classes[1][1][81] = int_stack + 657;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 666;
 Libderiv->deriv_classes[1][1][8] = int_stack + 684;
 Libderiv->deriv2_classes[1][1][80] = int_stack + 693;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 702;
 Libderiv->deriv_classes[1][1][7] = int_stack + 720;
 Libderiv->deriv2_classes[1][1][79] = int_stack + 729;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 738;
 Libderiv->deriv_classes[1][1][6] = int_stack + 756;
 Libderiv->deriv2_classes[1][1][78] = int_stack + 765;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 774;
 Libderiv->deriv2_classes[1][1][35] = int_stack + 792;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 801;
 Libderiv->deriv2_classes[1][1][34] = int_stack + 819;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 828;
 Libderiv->deriv2_classes[1][1][33] = int_stack + 846;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 855;
 Libderiv->deriv2_classes[1][1][32] = int_stack + 873;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 882;
 Libderiv->deriv2_classes[1][1][31] = int_stack + 900;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 909;
 Libderiv->deriv_classes[1][1][2] = int_stack + 927;
 Libderiv->deriv2_classes[1][1][30] = int_stack + 936;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 945;
 Libderiv->deriv2_classes[1][1][26] = int_stack + 963;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 972;
 Libderiv->deriv2_classes[1][1][23] = int_stack + 990;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 999;
 Libderiv->deriv2_classes[1][1][22] = int_stack + 1017;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 1026;
 Libderiv->deriv2_classes[1][1][21] = int_stack + 1044;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 1053;
 Libderiv->deriv2_classes[1][1][20] = int_stack + 1071;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 1080;
 Libderiv->deriv2_classes[1][1][19] = int_stack + 1098;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 1107;
 Libderiv->deriv_classes[1][1][1] = int_stack + 1125;
 Libderiv->deriv2_classes[1][1][18] = int_stack + 1134;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 1143;
 Libderiv->deriv2_classes[1][1][14] = int_stack + 1161;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 1170;
 Libderiv->deriv2_classes[1][1][13] = int_stack + 1188;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 1197;
 Libderiv->deriv2_classes[1][1][11] = int_stack + 1215;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 1224;
 Libderiv->deriv2_classes[1][1][10] = int_stack + 1242;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 1251;
 Libderiv->deriv2_classes[1][1][9] = int_stack + 1269;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 1278;
 Libderiv->deriv2_classes[1][1][8] = int_stack + 1296;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 1305;
 Libderiv->deriv2_classes[1][1][7] = int_stack + 1323;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 1332;
 Libderiv->deriv_classes[1][1][0] = int_stack + 1350;
 Libderiv->deriv2_classes[1][1][6] = int_stack + 1359;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 1368;
 Libderiv->deriv2_classes[1][1][2] = int_stack + 1386;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 1395;
 Libderiv->deriv2_classes[1][1][1] = int_stack + 1413;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 1422;
 Libderiv->deriv2_classes[1][1][0] = int_stack + 1440;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 1449;
 memset(int_stack,0,11736);

 Libderiv->dvrr_stack = int_stack + 1575;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0pp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1467,int_stack+0,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90,3);
     Libderiv->ABCD[11] = int_stack + 1467;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1494,int_stack+18,int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 1494;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+0,int_stack+36,int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+27,int_stack+54,int_stack+684, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 27;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1521,int_stack+72,int_stack+720, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 1521;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+54,int_stack+99,int_stack+756, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 54;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+81,int_stack+117,int_stack+927,3);
     Libderiv->ABCD[2] = int_stack + 81;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+108,int_stack+135,int_stack+1125,3);
     Libderiv->ABCD[1] = int_stack + 108;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+1548,int_stack+153,int_stack+1350,3);
     Libderiv->ABCD[0] = int_stack + 1548;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+135,int_stack+180,int_stack+171, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+576,3);
     Libderiv->ABCD[155] = int_stack + 135;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+162,int_stack+207,int_stack+198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+576, 1.0, int_stack+612,3);
     Libderiv->ABCD[143] = int_stack + 162;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+189,int_stack+234,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+612, 0.0, zero_stack,3);
     Libderiv->ABCD[142] = int_stack + 189;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+216,int_stack+261,int_stack+252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+576, 0.0, zero_stack, 1.0, int_stack+648,3);
     Libderiv->ABCD[131] = int_stack + 216;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+243,int_stack+288,int_stack+279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+612, 1.0, int_stack+648, 0.0, zero_stack,3);
     Libderiv->ABCD[130] = int_stack + 243;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+270,int_stack+315,int_stack+306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[129] = int_stack + 270;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+297,int_stack+342,int_stack+333, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+684,3);
     Libderiv->ABCD[119] = int_stack + 297;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+324,int_stack+369,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+612, 0.0, zero_stack, 1.0, int_stack+684, 0.0, zero_stack,3);
     Libderiv->ABCD[118] = int_stack + 324;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+351,int_stack+396,int_stack+387, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648, 1.0, int_stack+684, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[117] = int_stack + 351;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+378,int_stack+423,int_stack+414, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[116] = int_stack + 378;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+405,int_stack+450,int_stack+441, 0.0, zero_stack, 1.0, int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720,3);
     Libderiv->ABCD[107] = int_stack + 405;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+432,int_stack+477,int_stack+468, 0.0, zero_stack, 1.0, int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack,3);
     Libderiv->ABCD[106] = int_stack + 432;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+459,int_stack+504,int_stack+495, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[105] = int_stack + 459;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+486,int_stack+531,int_stack+522, 0.0, zero_stack, 1.0, int_stack+684, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[104] = int_stack + 486;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+513,int_stack+558,int_stack+549, 0.0, zero_stack, 2.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[103] = int_stack + 513;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+540,int_stack+594,int_stack+585, 1.0, int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+756,3);
     Libderiv->ABCD[95] = int_stack + 540;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+567,int_stack+630,int_stack+621, 1.0, int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+756, 0.0, zero_stack,3);
     Libderiv->ABCD[94] = int_stack + 567;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+594,int_stack+666,int_stack+657, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+756, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[93] = int_stack + 594;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+621,int_stack+702,int_stack+693, 1.0, int_stack+684, 0.0, zero_stack, 1.0, int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[92] = int_stack + 621;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+648,int_stack+738,int_stack+729, 1.0, int_stack+720, 1.0, int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[91] = int_stack + 648;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+675,int_stack+774,int_stack+765, 2.0, int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[90] = int_stack + 675;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+702,int_stack+801,int_stack+792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+927,3);
     Libderiv->ABCD[47] = int_stack + 702;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+729,int_stack+828,int_stack+819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+927, 0.0, zero_stack,3);
     Libderiv->ABCD[46] = int_stack + 729;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+756,int_stack+855,int_stack+846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+927, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[45] = int_stack + 756;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+783,int_stack+882,int_stack+873, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[44] = int_stack + 783;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+810,int_stack+909,int_stack+900, 0.0, zero_stack, 1.0, int_stack+927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[43] = int_stack + 810;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+837,int_stack+945,int_stack+936, 1.0, int_stack+927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[42] = int_stack + 837;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+864,int_stack+972,int_stack+963,3);
     Libderiv->ABCD[38] = int_stack + 864;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+891,int_stack+999,int_stack+990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125,3);
     Libderiv->ABCD[35] = int_stack + 891;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+918,int_stack+1026,int_stack+1017, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack,3);
     Libderiv->ABCD[34] = int_stack + 918;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+945,int_stack+1053,int_stack+1044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[33] = int_stack + 945;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+972,int_stack+1080,int_stack+1071, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[32] = int_stack + 972;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+999,int_stack+1107,int_stack+1098, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[31] = int_stack + 999;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1026,int_stack+1143,int_stack+1134, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[30] = int_stack + 1026;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+1053,int_stack+1170,int_stack+1161,3);
     Libderiv->ABCD[26] = int_stack + 1053;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+1080,int_stack+1197,int_stack+1188,3);
     Libderiv->ABCD[25] = int_stack + 1080;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1107,int_stack+1224,int_stack+1215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350,3);
     Libderiv->ABCD[23] = int_stack + 1107;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1134,int_stack+1251,int_stack+1242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack,3);
     Libderiv->ABCD[22] = int_stack + 1134;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1161,int_stack+1278,int_stack+1269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[21] = int_stack + 1161;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1188,int_stack+1305,int_stack+1296, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[20] = int_stack + 1188;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1215,int_stack+1332,int_stack+1323, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[19] = int_stack + 1215;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1242,int_stack+1368,int_stack+1359, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[18] = int_stack + 1242;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+1269,int_stack+1395,int_stack+1386,3);
     Libderiv->ABCD[14] = int_stack + 1269;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+1296,int_stack+1422,int_stack+1413,3);
     Libderiv->ABCD[13] = int_stack + 1296;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+1323,int_stack+1449,int_stack+1440,3);
     Libderiv->ABCD[12] = int_stack + 1323;

}
