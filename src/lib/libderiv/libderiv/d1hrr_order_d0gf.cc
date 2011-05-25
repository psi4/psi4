#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0gf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|gf) integrals */

void d1hrr_order_d0gf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][6][11] = int_stack + 216;
 Libderiv->deriv_classes[2][7][11] = int_stack + 384;
 Libderiv->deriv_classes[2][4][10] = int_stack + 600;
 Libderiv->deriv_classes[2][5][10] = int_stack + 690;
 Libderiv->deriv_classes[2][6][10] = int_stack + 816;
 Libderiv->deriv_classes[2][7][10] = int_stack + 984;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1200;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1290;
 Libderiv->deriv_classes[2][6][9] = int_stack + 1416;
 Libderiv->deriv_classes[2][7][9] = int_stack + 1584;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1800;
 Libderiv->deriv_classes[2][5][8] = int_stack + 1890;
 Libderiv->deriv_classes[2][6][8] = int_stack + 2016;
 Libderiv->deriv_classes[2][7][8] = int_stack + 2184;
 Libderiv->deriv_classes[2][4][7] = int_stack + 2400;
 Libderiv->deriv_classes[2][5][7] = int_stack + 2490;
 Libderiv->deriv_classes[2][6][7] = int_stack + 2616;
 Libderiv->deriv_classes[2][7][7] = int_stack + 2784;
 Libderiv->dvrr_classes[2][4] = int_stack + 3000;
 Libderiv->deriv_classes[2][4][6] = int_stack + 3090;
 Libderiv->dvrr_classes[2][5] = int_stack + 3180;
 Libderiv->deriv_classes[2][5][6] = int_stack + 3306;
 Libderiv->dvrr_classes[2][6] = int_stack + 3432;
 Libderiv->deriv_classes[2][6][6] = int_stack + 3600;
 Libderiv->deriv_classes[2][7][6] = int_stack + 3768;
 Libderiv->deriv_classes[2][4][2] = int_stack + 3984;
 Libderiv->deriv_classes[2][5][2] = int_stack + 4074;
 Libderiv->deriv_classes[2][6][2] = int_stack + 4200;
 Libderiv->deriv_classes[2][7][2] = int_stack + 4368;
 Libderiv->deriv_classes[2][4][1] = int_stack + 4584;
 Libderiv->deriv_classes[2][5][1] = int_stack + 4674;
 Libderiv->deriv_classes[2][6][1] = int_stack + 4800;
 Libderiv->deriv_classes[2][7][1] = int_stack + 4968;
 Libderiv->deriv_classes[2][4][0] = int_stack + 5184;
 Libderiv->deriv_classes[2][5][0] = int_stack + 5274;
 Libderiv->deriv_classes[2][6][0] = int_stack + 5400;
 Libderiv->deriv_classes[2][7][0] = int_stack + 5568;
 memset(int_stack,0,46272);

 Libderiv->dvrr_stack = int_stack + 14964;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0gf(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5784,int_stack+3180,int_stack+3000,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6054,int_stack+3432,int_stack+3180,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6432,int_stack+6054,int_stack+5784,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6972,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7242,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3180,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7620,int_stack+7242,int_stack+6972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5784,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+384,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3432,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8664,int_stack+8160,int_stack+7242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6054,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8160,int_stack+690,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6972,int_stack+816,int_stack+690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3180, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+6972,int_stack+8160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5784, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+984,int_stack+816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3432, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+9420,int_stack+8160,int_stack+6972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6054, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6972,int_stack+1290,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7242,int_stack+1416,int_stack+1290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3180, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+540,int_stack+7242,int_stack+6972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5784, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+1584,int_stack+1416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3432, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10176,int_stack+8160,int_stack+7242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6054, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8160,int_stack+1890,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6972,int_stack+2016,int_stack+1890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+6972,int_stack+8160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+2184,int_stack+2016, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1620,int_stack+8160,int_stack+6972, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6972,int_stack+2490,int_stack+2400, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7242,int_stack+2616,int_stack+2490, 0.0, zero_stack, 1.0, int_stack+3180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10932,int_stack+7242,int_stack+6972, 0.0, zero_stack, 1.0, int_stack+5784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+2784,int_stack+2616, 0.0, zero_stack, 1.0, int_stack+3432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+11472,int_stack+8160,int_stack+7242, 0.0, zero_stack, 1.0, int_stack+6054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8160,int_stack+3306,int_stack+3090, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6972,int_stack+3600,int_stack+3306, 1.0, int_stack+3180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2376,int_stack+6972,int_stack+8160, 1.0, int_stack+5784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+3768,int_stack+3600, 1.0, int_stack+3432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2916,int_stack+8160,int_stack+6972, 1.0, int_stack+6054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6972,int_stack+4074,int_stack+3984,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7242,int_stack+4200,int_stack+4074,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5784,int_stack+7242,int_stack+6972,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+4368,int_stack+4200,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+3672,int_stack+8160,int_stack+7242,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8160,int_stack+4674,int_stack+4584,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6972,int_stack+4800,int_stack+4674,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12228,int_stack+6972,int_stack+8160,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+4968,int_stack+4800,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+4428,int_stack+8160,int_stack+6972,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6972,int_stack+5274,int_stack+5184,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7242,int_stack+5400,int_stack+5274,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12768,int_stack+7242,int_stack+6972,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8160,int_stack+5568,int_stack+5400,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+13308,int_stack+8160,int_stack+7242,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+14064,int_stack+8664,int_stack+7620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6432,6);
     Libderiv->ABCD[11] = int_stack + 14064;
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6972,int_stack+9420,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6432, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 6972;
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+7872,int_stack+10176,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6432, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 7872;
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+1620,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+900,int_stack+11472,int_stack+10932, 0.0, zero_stack, 1.0, int_stack+6432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 900;
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8772,int_stack+2916,int_stack+2376, 1.0, int_stack+6432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 8772;
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+1800,int_stack+3672,int_stack+5784,6);
     Libderiv->ABCD[2] = int_stack + 1800;
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+2700,int_stack+4428,int_stack+12228,6);
     Libderiv->ABCD[1] = int_stack + 2700;
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+3600,int_stack+13308,int_stack+12768,6);
     Libderiv->ABCD[0] = int_stack + 3600;

}
