#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|ff) integrals */

void d1hrr_order_g0ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][3][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][11] = int_stack + 150;
 Libderiv->deriv_classes[4][5][11] = int_stack + 375;
 Libderiv->deriv_classes[4][6][11] = int_stack + 690;
 Libderiv->deriv_classes[4][3][10] = int_stack + 1110;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1260;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1485;
 Libderiv->deriv_classes[4][6][10] = int_stack + 1800;
 Libderiv->deriv_classes[4][3][9] = int_stack + 2220;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2370;
 Libderiv->deriv_classes[4][5][9] = int_stack + 2595;
 Libderiv->deriv_classes[4][6][9] = int_stack + 2910;
 Libderiv->deriv_classes[4][3][8] = int_stack + 3330;
 Libderiv->deriv_classes[4][4][8] = int_stack + 3480;
 Libderiv->deriv_classes[4][5][8] = int_stack + 3705;
 Libderiv->deriv_classes[4][6][8] = int_stack + 4020;
 Libderiv->deriv_classes[4][3][7] = int_stack + 4440;
 Libderiv->deriv_classes[4][4][7] = int_stack + 4590;
 Libderiv->deriv_classes[4][5][7] = int_stack + 4815;
 Libderiv->deriv_classes[4][6][7] = int_stack + 5130;
 Libderiv->dvrr_classes[4][3] = int_stack + 5550;
 Libderiv->deriv_classes[4][3][6] = int_stack + 5700;
 Libderiv->dvrr_classes[4][4] = int_stack + 5850;
 Libderiv->deriv_classes[4][4][6] = int_stack + 6075;
 Libderiv->dvrr_classes[4][5] = int_stack + 6300;
 Libderiv->deriv_classes[4][5][6] = int_stack + 6615;
 Libderiv->deriv_classes[4][6][6] = int_stack + 6930;
 Libderiv->deriv_classes[4][3][2] = int_stack + 7350;
 Libderiv->deriv_classes[4][4][2] = int_stack + 7500;
 Libderiv->deriv_classes[4][5][2] = int_stack + 7725;
 Libderiv->deriv_classes[4][6][2] = int_stack + 8040;
 Libderiv->deriv_classes[4][3][1] = int_stack + 8460;
 Libderiv->deriv_classes[4][4][1] = int_stack + 8610;
 Libderiv->deriv_classes[4][5][1] = int_stack + 8835;
 Libderiv->deriv_classes[4][6][1] = int_stack + 9150;
 Libderiv->deriv_classes[4][3][0] = int_stack + 9570;
 Libderiv->deriv_classes[4][4][0] = int_stack + 9720;
 Libderiv->deriv_classes[4][5][0] = int_stack + 9945;
 Libderiv->deriv_classes[4][6][0] = int_stack + 10260;
 memset(int_stack,0,85440);

 Libderiv->dvrr_stack = int_stack + 26445;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10680,int_stack+5850,int_stack+5550,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11130,int_stack+6300,int_stack+5850,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11805,int_stack+11130,int_stack+10680,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12705,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5550,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13155,int_stack+375,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5850,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13830,int_stack+13155,int_stack+12705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10680,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14730,int_stack+690,int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6300,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15675,int_stack+14730,int_stack+13155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11130,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14730,int_stack+1260,int_stack+1110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5550, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12705,int_stack+1485,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5850, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+12705,int_stack+14730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10680, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14730,int_stack+1800,int_stack+1485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6300, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17025,int_stack+14730,int_stack+12705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11130, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12705,int_stack+2370,int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5550, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13155,int_stack+2595,int_stack+2370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5850, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14730,int_stack+13155,int_stack+12705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10680, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+900,int_stack+2910,int_stack+2595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6300, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1845,int_stack+900,int_stack+13155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11130, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+3480,int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12705,int_stack+3705,int_stack+3480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18375,int_stack+12705,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+900,int_stack+4020,int_stack+3705, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19275,int_stack+900,int_stack+12705, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12705,int_stack+4590,int_stack+4440, 0.0, zero_stack, 1.0, int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13155,int_stack+4815,int_stack+4590, 0.0, zero_stack, 1.0, int_stack+5850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+900,int_stack+13155,int_stack+12705, 0.0, zero_stack, 1.0, int_stack+10680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3195,int_stack+5130,int_stack+4815, 0.0, zero_stack, 1.0, int_stack+6300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4140,int_stack+3195,int_stack+13155, 0.0, zero_stack, 1.0, int_stack+11130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3195,int_stack+6075,int_stack+5700, 1.0, int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12705,int_stack+6615,int_stack+6075, 1.0, int_stack+5850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20625,int_stack+12705,int_stack+3195, 1.0, int_stack+10680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3195,int_stack+6930,int_stack+6615, 1.0, int_stack+6300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5490,int_stack+3195,int_stack+12705, 1.0, int_stack+11130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12705,int_stack+7500,int_stack+7350,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13155,int_stack+7725,int_stack+7500,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3195,int_stack+13155,int_stack+12705,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10680,int_stack+8040,int_stack+7725,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6840,int_stack+10680,int_stack+13155,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10680,int_stack+8610,int_stack+8460,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11130,int_stack+8835,int_stack+8610,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12705,int_stack+11130,int_stack+10680,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+21525,int_stack+9150,int_stack+8835,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8190,int_stack+21525,int_stack+11130,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+21525,int_stack+9720,int_stack+9570,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+21975,int_stack+9945,int_stack+9720,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10680,int_stack+21975,int_stack+21525,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+22650,int_stack+10260,int_stack+9945,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23595,int_stack+22650,int_stack+21975,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+21525,int_stack+15675,int_stack+13830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11805,15);
     Libderiv->ABCD[11] = int_stack + 21525;
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+24945,int_stack+17025,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11805, 0.0, zero_stack,15);
     Libderiv->ABCD[10] = int_stack + 24945;
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+15630,int_stack+1845,int_stack+14730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11805, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[9] = int_stack + 15630;
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13605,int_stack+19275,int_stack+18375, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[8] = int_stack + 13605;
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+17130,int_stack+4140,int_stack+900, 0.0, zero_stack, 1.0, int_stack+11805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[7] = int_stack + 17130;
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+5490,int_stack+20625, 1.0, int_stack+11805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1500,int_stack+6840,int_stack+3195,15);
     Libderiv->ABCD[2] = int_stack + 1500;
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+3000,int_stack+8190,int_stack+12705,15);
     Libderiv->ABCD[1] = int_stack + 3000;
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+4500,int_stack+23595,int_stack+10680,15);
     Libderiv->ABCD[0] = int_stack + 4500;

}
