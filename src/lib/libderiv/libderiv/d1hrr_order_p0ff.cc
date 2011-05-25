#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|ff) integrals */

void d1hrr_order_p0ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][11] = int_stack + 30;
 Libderiv->deriv_classes[1][5][11] = int_stack + 75;
 Libderiv->deriv_classes[1][6][11] = int_stack + 138;
 Libderiv->deriv_classes[1][3][10] = int_stack + 222;
 Libderiv->deriv_classes[1][4][10] = int_stack + 252;
 Libderiv->deriv_classes[1][5][10] = int_stack + 297;
 Libderiv->deriv_classes[1][6][10] = int_stack + 360;
 Libderiv->deriv_classes[1][3][9] = int_stack + 444;
 Libderiv->deriv_classes[1][4][9] = int_stack + 474;
 Libderiv->deriv_classes[1][5][9] = int_stack + 519;
 Libderiv->deriv_classes[1][6][9] = int_stack + 582;
 Libderiv->deriv_classes[1][3][8] = int_stack + 666;
 Libderiv->deriv_classes[1][4][8] = int_stack + 696;
 Libderiv->deriv_classes[1][5][8] = int_stack + 741;
 Libderiv->deriv_classes[1][6][8] = int_stack + 804;
 Libderiv->deriv_classes[1][3][7] = int_stack + 888;
 Libderiv->deriv_classes[1][4][7] = int_stack + 918;
 Libderiv->deriv_classes[1][5][7] = int_stack + 963;
 Libderiv->deriv_classes[1][6][7] = int_stack + 1026;
 Libderiv->dvrr_classes[1][3] = int_stack + 1110;
 Libderiv->deriv_classes[1][3][6] = int_stack + 1140;
 Libderiv->dvrr_classes[1][4] = int_stack + 1170;
 Libderiv->deriv_classes[1][4][6] = int_stack + 1215;
 Libderiv->dvrr_classes[1][5] = int_stack + 1260;
 Libderiv->deriv_classes[1][5][6] = int_stack + 1323;
 Libderiv->deriv_classes[1][6][6] = int_stack + 1386;
 Libderiv->deriv_classes[1][3][2] = int_stack + 1470;
 Libderiv->deriv_classes[1][4][2] = int_stack + 1500;
 Libderiv->deriv_classes[1][5][2] = int_stack + 1545;
 Libderiv->deriv_classes[1][6][2] = int_stack + 1608;
 Libderiv->deriv_classes[1][3][1] = int_stack + 1692;
 Libderiv->deriv_classes[1][4][1] = int_stack + 1722;
 Libderiv->deriv_classes[1][5][1] = int_stack + 1767;
 Libderiv->deriv_classes[1][6][1] = int_stack + 1830;
 Libderiv->deriv_classes[1][3][0] = int_stack + 1914;
 Libderiv->deriv_classes[1][4][0] = int_stack + 1944;
 Libderiv->deriv_classes[1][5][0] = int_stack + 1989;
 Libderiv->deriv_classes[1][6][0] = int_stack + 2052;
 memset(int_stack,0,17088);

 Libderiv->dvrr_stack = int_stack + 5910;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2136,int_stack+1170,int_stack+1110,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2226,int_stack+1260,int_stack+1170,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2361,int_stack+2226,int_stack+2136,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2541,int_stack+30,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1110,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2631,int_stack+75,int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2766,int_stack+2631,int_stack+2541, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2136,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2946,int_stack+138,int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3135,int_stack+2946,int_stack+2631, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2226,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2946,int_stack+252,int_stack+222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1110, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2541,int_stack+297,int_stack+252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+2541,int_stack+2946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2136, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2946,int_stack+360,int_stack+297, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3405,int_stack+2946,int_stack+2541, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2226, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2541,int_stack+474,int_stack+444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1110, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2631,int_stack+519,int_stack+474, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2946,int_stack+2631,int_stack+2541, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2136, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+180,int_stack+582,int_stack+519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+369,int_stack+180,int_stack+2631, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2226, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+696,int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2541,int_stack+741,int_stack+696, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3675,int_stack+2541,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+180,int_stack+804,int_stack+741, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3855,int_stack+180,int_stack+2541, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2541,int_stack+918,int_stack+888, 0.0, zero_stack, 1.0, int_stack+1110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2631,int_stack+963,int_stack+918, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+180,int_stack+2631,int_stack+2541, 0.0, zero_stack, 1.0, int_stack+2136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+639,int_stack+1026,int_stack+963, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+828,int_stack+639,int_stack+2631, 0.0, zero_stack, 1.0, int_stack+2226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+639,int_stack+1215,int_stack+1140, 1.0, int_stack+1110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2541,int_stack+1323,int_stack+1215, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4125,int_stack+2541,int_stack+639, 1.0, int_stack+2136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+639,int_stack+1386,int_stack+1323, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1098,int_stack+639,int_stack+2541, 1.0, int_stack+2226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2541,int_stack+1500,int_stack+1470,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2631,int_stack+1545,int_stack+1500,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+639,int_stack+2631,int_stack+2541,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2136,int_stack+1608,int_stack+1545,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1368,int_stack+2136,int_stack+2631,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2136,int_stack+1722,int_stack+1692,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2226,int_stack+1767,int_stack+1722,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2541,int_stack+2226,int_stack+2136,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4305,int_stack+1830,int_stack+1767,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1638,int_stack+4305,int_stack+2226,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4305,int_stack+1944,int_stack+1914,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4395,int_stack+1989,int_stack+1944,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4530,int_stack+4395,int_stack+4305,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4710,int_stack+2052,int_stack+1989,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1908,int_stack+4710,int_stack+4395,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4710,int_stack+3135,int_stack+2766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2361,3);
     Libderiv->ABCD[11] = int_stack + 4710;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5010,int_stack+3405,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2361, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 5010;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5310,int_stack+369,int_stack+2946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2361, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 5310;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5610,int_stack+3855,int_stack+3675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2361, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 5610;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2721,int_stack+828,int_stack+180, 0.0, zero_stack, 1.0, int_stack+2361, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 2721;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+1098,int_stack+4125, 1.0, int_stack+2361, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+300,int_stack+1368,int_stack+639,3);
     Libderiv->ABCD[2] = int_stack + 300;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+600,int_stack+1638,int_stack+2541,3);
     Libderiv->ABCD[1] = int_stack + 600;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+900,int_stack+1908,int_stack+4530,3);
     Libderiv->ABCD[0] = int_stack + 900;

}
