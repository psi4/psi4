#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|fd) integrals */

void d1hrr_order_ddfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[2][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][3][11] = int_stack + 276;
 Libderiv->deriv_classes[3][4][11] = int_stack + 376;
 Libderiv->deriv_classes[3][5][11] = int_stack + 526;
 Libderiv->deriv_classes[4][3][11] = int_stack + 736;
 Libderiv->deriv_classes[4][4][11] = int_stack + 886;
 Libderiv->deriv_classes[4][5][11] = int_stack + 1111;
 Libderiv->deriv_classes[2][3][10] = int_stack + 1426;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1486;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1576;
 Libderiv->deriv_classes[3][3][10] = int_stack + 1702;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1802;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1952;
 Libderiv->deriv_classes[4][3][10] = int_stack + 2162;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2312;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2537;
 Libderiv->deriv_classes[2][3][9] = int_stack + 2852;
 Libderiv->deriv_classes[2][4][9] = int_stack + 2912;
 Libderiv->deriv_classes[2][5][9] = int_stack + 3002;
 Libderiv->deriv_classes[3][3][9] = int_stack + 3128;
 Libderiv->deriv_classes[3][4][9] = int_stack + 3228;
 Libderiv->deriv_classes[3][5][9] = int_stack + 3378;
 Libderiv->deriv_classes[4][3][9] = int_stack + 3588;
 Libderiv->deriv_classes[4][4][9] = int_stack + 3738;
 Libderiv->deriv_classes[4][5][9] = int_stack + 3963;
 Libderiv->deriv_classes[2][3][8] = int_stack + 4278;
 Libderiv->deriv_classes[2][4][8] = int_stack + 4338;
 Libderiv->deriv_classes[2][5][8] = int_stack + 4428;
 Libderiv->deriv_classes[3][3][8] = int_stack + 4554;
 Libderiv->deriv_classes[3][4][8] = int_stack + 4654;
 Libderiv->deriv_classes[3][5][8] = int_stack + 4804;
 Libderiv->deriv_classes[4][3][8] = int_stack + 5014;
 Libderiv->deriv_classes[4][4][8] = int_stack + 5164;
 Libderiv->deriv_classes[4][5][8] = int_stack + 5389;
 Libderiv->deriv_classes[2][3][7] = int_stack + 5704;
 Libderiv->deriv_classes[2][4][7] = int_stack + 5764;
 Libderiv->deriv_classes[2][5][7] = int_stack + 5854;
 Libderiv->deriv_classes[3][3][7] = int_stack + 5980;
 Libderiv->deriv_classes[3][4][7] = int_stack + 6080;
 Libderiv->deriv_classes[3][5][7] = int_stack + 6230;
 Libderiv->deriv_classes[4][3][7] = int_stack + 6440;
 Libderiv->deriv_classes[4][4][7] = int_stack + 6590;
 Libderiv->deriv_classes[4][5][7] = int_stack + 6815;
 Libderiv->deriv_classes[2][3][6] = int_stack + 7130;
 Libderiv->deriv_classes[2][4][6] = int_stack + 7190;
 Libderiv->deriv_classes[2][5][6] = int_stack + 7280;
 Libderiv->deriv_classes[3][3][6] = int_stack + 7406;
 Libderiv->deriv_classes[3][4][6] = int_stack + 7506;
 Libderiv->deriv_classes[3][5][6] = int_stack + 7656;
 Libderiv->dvrr_classes[4][3] = int_stack + 7866;
 Libderiv->deriv_classes[4][3][6] = int_stack + 8016;
 Libderiv->dvrr_classes[4][4] = int_stack + 8166;
 Libderiv->deriv_classes[4][4][6] = int_stack + 8391;
 Libderiv->deriv_classes[4][5][6] = int_stack + 8616;
 Libderiv->deriv_classes[2][3][2] = int_stack + 8931;
 Libderiv->deriv_classes[2][4][2] = int_stack + 8991;
 Libderiv->deriv_classes[2][5][2] = int_stack + 9081;
 Libderiv->deriv_classes[3][3][2] = int_stack + 9207;
 Libderiv->deriv_classes[3][4][2] = int_stack + 9307;
 Libderiv->deriv_classes[3][5][2] = int_stack + 9457;
 Libderiv->deriv_classes[4][3][2] = int_stack + 9667;
 Libderiv->deriv_classes[4][4][2] = int_stack + 9817;
 Libderiv->deriv_classes[4][5][2] = int_stack + 10042;
 Libderiv->deriv_classes[2][3][1] = int_stack + 10357;
 Libderiv->deriv_classes[2][4][1] = int_stack + 10417;
 Libderiv->deriv_classes[2][5][1] = int_stack + 10507;
 Libderiv->deriv_classes[3][3][1] = int_stack + 10633;
 Libderiv->deriv_classes[3][4][1] = int_stack + 10733;
 Libderiv->deriv_classes[3][5][1] = int_stack + 10883;
 Libderiv->deriv_classes[4][3][1] = int_stack + 11093;
 Libderiv->deriv_classes[4][4][1] = int_stack + 11243;
 Libderiv->deriv_classes[4][5][1] = int_stack + 11468;
 Libderiv->dvrr_classes[2][3] = int_stack + 11783;
 Libderiv->dvrr_classes[2][4] = int_stack + 11843;
 Libderiv->dvrr_classes[2][5] = int_stack + 11933;
 Libderiv->deriv_classes[2][3][0] = int_stack + 12059;
 Libderiv->deriv_classes[2][4][0] = int_stack + 12119;
 Libderiv->deriv_classes[2][5][0] = int_stack + 12209;
 Libderiv->dvrr_classes[3][3] = int_stack + 12335;
 Libderiv->dvrr_classes[3][4] = int_stack + 12435;
 Libderiv->dvrr_classes[3][5] = int_stack + 12585;
 Libderiv->deriv_classes[3][3][0] = int_stack + 12795;
 Libderiv->deriv_classes[3][4][0] = int_stack + 12895;
 Libderiv->deriv_classes[3][5][0] = int_stack + 13045;
 Libderiv->deriv_classes[4][3][0] = int_stack + 13255;
 Libderiv->deriv_classes[4][4][0] = int_stack + 13405;
 Libderiv->deriv_classes[4][5][0] = int_stack + 13630;
 memset(int_stack,0,111560);

 Libderiv->dvrr_stack = int_stack + 35050;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13945,int_stack+11843,int_stack+11783,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14125,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11783,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14305,int_stack+150,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11843,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14575,int_stack+14305,int_stack+14125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13945,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14125,int_stack+12435,int_stack+12335,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14935,int_stack+376,int_stack+276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12335,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15235,int_stack+526,int_stack+376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12435,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+15235,int_stack+14935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14125,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14935,int_stack+0,int_stack+14575,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16015,int_stack+8166,int_stack+7866,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16465,int_stack+886,int_stack+736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7866,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16915,int_stack+1111,int_stack+886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8166,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17590,int_stack+16915,int_stack+16465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16015,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+18490,int_stack+17590,int_stack+0,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1486,int_stack+1426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11783, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180,int_stack+1576,int_stack+1486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11843, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13945, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1802,int_stack+1702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12335, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+810,int_stack+1952,int_stack+1802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12435, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1260,int_stack+810,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14125, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+16465,int_stack+1260,int_stack+450,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2312,int_stack+2162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7866, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+2537,int_stack+2312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8166, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1860,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16015, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+20290,int_stack+1860,int_stack+1260,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2912,int_stack+2852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11783, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180,int_stack+3002,int_stack+2912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11843, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13945, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+3228,int_stack+3128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12335, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+810,int_stack+3378,int_stack+3228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12435, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1260,int_stack+810,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14125, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1860,int_stack+1260,int_stack+450,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+3738,int_stack+3588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7866, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+3963,int_stack+3738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8166, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2940,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16015, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+22090,int_stack+2940,int_stack+1260,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2940,int_stack+4338,int_stack+4278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11783, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3120,int_stack+4428,int_stack+4338, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11843, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3390,int_stack+3120,int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2940,int_stack+4654,int_stack+4554, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3750,int_stack+4804,int_stack+4654, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4200,int_stack+3750,int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+4200,int_stack+3390,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2940,int_stack+5164,int_stack+5014, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7866, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3390,int_stack+5389,int_stack+5164, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4800,int_stack+3390,int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+23890,int_stack+4800,int_stack+4200,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2940,int_stack+5764,int_stack+5704, 0.0, zero_stack, 1.0, int_stack+11783, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3120,int_stack+5854,int_stack+5764, 0.0, zero_stack, 1.0, int_stack+11843, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3390,int_stack+3120,int_stack+2940, 0.0, zero_stack, 1.0, int_stack+13945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2940,int_stack+6080,int_stack+5980, 0.0, zero_stack, 1.0, int_stack+12335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3750,int_stack+6230,int_stack+6080, 0.0, zero_stack, 1.0, int_stack+12435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4200,int_stack+3750,int_stack+2940, 0.0, zero_stack, 1.0, int_stack+14125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4800,int_stack+4200,int_stack+3390,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2940,int_stack+6590,int_stack+6440, 0.0, zero_stack, 1.0, int_stack+7866, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3390,int_stack+6815,int_stack+6590, 0.0, zero_stack, 1.0, int_stack+8166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5880,int_stack+3390,int_stack+2940, 0.0, zero_stack, 1.0, int_stack+16015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+25690,int_stack+5880,int_stack+4200,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5880,int_stack+7190,int_stack+7130, 1.0, int_stack+11783, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6060,int_stack+7280,int_stack+7190, 1.0, int_stack+11843, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6330,int_stack+6060,int_stack+5880, 1.0, int_stack+13945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5880,int_stack+7506,int_stack+7406, 1.0, int_stack+12335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6690,int_stack+7656,int_stack+7506, 1.0, int_stack+12435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7140,int_stack+6690,int_stack+5880, 1.0, int_stack+14125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2940,int_stack+7140,int_stack+6330,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5880,int_stack+8391,int_stack+8016, 1.0, int_stack+7866, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6330,int_stack+8616,int_stack+8391, 1.0, int_stack+8166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7740,int_stack+6330,int_stack+5880, 1.0, int_stack+16015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+27490,int_stack+7740,int_stack+7140,60);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16015,int_stack+11933,int_stack+11843,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5880,int_stack+16015,int_stack+13945,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16015,int_stack+12585,int_stack+12435,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+6240,int_stack+16015,int_stack+14125,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6840,int_stack+6240,int_stack+5880,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16015,int_stack+8991,int_stack+8931,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16195,int_stack+9081,int_stack+8991,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+13945,int_stack+16195,int_stack+16015,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16015,int_stack+9307,int_stack+9207,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14305,int_stack+9457,int_stack+9307,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7920,int_stack+14305,int_stack+16015,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+8520,int_stack+7920,int_stack+13945, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16015,int_stack+9817,int_stack+9667,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13945,int_stack+10042,int_stack+9817,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17545,int_stack+13945,int_stack+16015,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+29290,int_stack+17545,int_stack+7920, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7920,int_stack+10417,int_stack+10357,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+10507,int_stack+10417,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17545,int_stack+8100,int_stack+7920,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7920,int_stack+10733,int_stack+10633,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16015,int_stack+10883,int_stack+10733,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+13945,int_stack+16015,int_stack+7920,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9600,int_stack+13945,int_stack+17545, 0.0, zero_stack, 1.0, int_stack+5880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16015,int_stack+11243,int_stack+11093,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17545,int_stack+11468,int_stack+11243,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10680,int_stack+17545,int_stack+16015,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+31090,int_stack+10680,int_stack+13945, 0.0, zero_stack, 1.0, int_stack+6240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13945,int_stack+12119,int_stack+12059,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14125,int_stack+12209,int_stack+12119,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14395,int_stack+14125,int_stack+13945,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13945,int_stack+12895,int_stack+12795,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16015,int_stack+13045,int_stack+12895,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7920,int_stack+16015,int_stack+13945,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10680,int_stack+7920,int_stack+14395, 1.0, int_stack+5880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16015,int_stack+13405,int_stack+13255,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13945,int_stack+13630,int_stack+13405,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11760,int_stack+13945,int_stack+16015,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+12660,int_stack+11760,int_stack+7920, 1.0, int_stack+6240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+32890,int_stack+18490,int_stack+14935,60);
     Libderiv->ABCD[11] = int_stack + 32890;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+17545,int_stack+20290,int_stack+16465,60);
     Libderiv->ABCD[10] = int_stack + 17545;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+14460,int_stack+22090,int_stack+1860,60);
     Libderiv->ABCD[9] = int_stack + 14460;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+19705,int_stack+23890,int_stack+0,60);
     Libderiv->ABCD[8] = int_stack + 19705;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+25690,int_stack+4800,60);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+21865,int_stack+27490,int_stack+2940,60);
     Libderiv->ABCD[6] = int_stack + 21865;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+2160,int_stack+29290,int_stack+8520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 2160;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+4320,int_stack+31090,int_stack+9600, 0.0, zero_stack, 1.0, int_stack+6840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 4320;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+7920,int_stack+12660,int_stack+10680, 1.0, int_stack+6840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 7920;

}
