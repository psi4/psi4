#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpdp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|dp) integrals */

void d1hrr_order_dpdp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][11] = int_stack + 36;
 Libderiv->deriv_classes[3][2][11] = int_stack + 96;
 Libderiv->deriv_classes[3][3][11] = int_stack + 156;
 Libderiv->deriv_classes[2][2][10] = int_stack + 256;
 Libderiv->deriv_classes[2][3][10] = int_stack + 292;
 Libderiv->deriv_classes[3][2][10] = int_stack + 352;
 Libderiv->deriv_classes[3][3][10] = int_stack + 412;
 Libderiv->deriv_classes[2][2][9] = int_stack + 512;
 Libderiv->deriv_classes[2][3][9] = int_stack + 548;
 Libderiv->deriv_classes[3][2][9] = int_stack + 608;
 Libderiv->deriv_classes[3][3][9] = int_stack + 668;
 Libderiv->deriv_classes[2][2][8] = int_stack + 768;
 Libderiv->deriv_classes[2][3][8] = int_stack + 804;
 Libderiv->deriv_classes[3][2][8] = int_stack + 864;
 Libderiv->deriv_classes[3][3][8] = int_stack + 924;
 Libderiv->deriv_classes[2][2][7] = int_stack + 1024;
 Libderiv->deriv_classes[2][3][7] = int_stack + 1060;
 Libderiv->deriv_classes[3][2][7] = int_stack + 1120;
 Libderiv->deriv_classes[3][3][7] = int_stack + 1180;
 Libderiv->deriv_classes[2][2][6] = int_stack + 1280;
 Libderiv->deriv_classes[2][3][6] = int_stack + 1316;
 Libderiv->dvrr_classes[3][2] = int_stack + 1376;
 Libderiv->deriv_classes[3][2][6] = int_stack + 1436;
 Libderiv->deriv_classes[3][3][6] = int_stack + 1496;
 Libderiv->deriv_classes[2][2][2] = int_stack + 1596;
 Libderiv->deriv_classes[2][3][2] = int_stack + 1632;
 Libderiv->deriv_classes[3][2][2] = int_stack + 1692;
 Libderiv->deriv_classes[3][3][2] = int_stack + 1752;
 Libderiv->deriv_classes[2][2][1] = int_stack + 1852;
 Libderiv->deriv_classes[2][3][1] = int_stack + 1888;
 Libderiv->deriv_classes[3][2][1] = int_stack + 1948;
 Libderiv->deriv_classes[3][3][1] = int_stack + 2008;
 Libderiv->dvrr_classes[2][2] = int_stack + 2108;
 Libderiv->dvrr_classes[2][3] = int_stack + 2144;
 Libderiv->deriv_classes[2][2][0] = int_stack + 2204;
 Libderiv->deriv_classes[2][3][0] = int_stack + 2240;
 Libderiv->deriv_classes[3][2][0] = int_stack + 2300;
 Libderiv->deriv_classes[3][3][0] = int_stack + 2360;
 memset(int_stack,0,19680);

 Libderiv->dvrr_stack = int_stack + 3900;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpdp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2460,int_stack+36,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2108,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2568,int_stack+156,int_stack+96, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1376,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+292,int_stack+256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2108, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+412,int_stack+352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1376, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+288,int_stack+548,int_stack+512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2108, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+396,int_stack+668,int_stack+608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1376, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+576,int_stack+804,int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+684,int_stack+924,int_stack+864, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+864,int_stack+1060,int_stack+1024, 0.0, zero_stack, 1.0, int_stack+2108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2748,int_stack+1180,int_stack+1120, 0.0, zero_stack, 1.0, int_stack+1376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+972,int_stack+1316,int_stack+1280, 1.0, int_stack+2108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1496,int_stack+1436, 1.0, int_stack+1376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1260,int_stack+2144,int_stack+2108,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1368,int_stack+1632,int_stack+1596,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1476,int_stack+1752,int_stack+1692,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1656,int_stack+1888,int_stack+1852,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1764,int_stack+2008,int_stack+1948,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1944,int_stack+2240,int_stack+2204,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2052,int_stack+2360,int_stack+2300,10);
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2928,int_stack+2568,int_stack+2460,18);
     Libderiv->ABCD[11] = int_stack + 2928;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2232,int_stack+108,int_stack+0,18);
     Libderiv->ABCD[10] = int_stack + 2232;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3252,int_stack+396,int_stack+288,18);
     Libderiv->ABCD[9] = int_stack + 3252;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+684,int_stack+576,18);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+324,int_stack+2748,int_stack+864,18);
     Libderiv->ABCD[7] = int_stack + 324;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+648,int_stack+1080,int_stack+972,18);
     Libderiv->ABCD[6] = int_stack + 648;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+2556,int_stack+1476,int_stack+1368, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[2] = int_stack + 2556;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+3576,int_stack+1764,int_stack+1656, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[1] = int_stack + 3576;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1368,int_stack+2052,int_stack+1944, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[0] = int_stack + 1368;

}
