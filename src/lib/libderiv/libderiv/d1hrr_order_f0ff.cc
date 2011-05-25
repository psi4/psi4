#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|ff) integrals */

void d1hrr_order_f0ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 100;
 Libderiv->deriv_classes[3][5][11] = int_stack + 250;
 Libderiv->deriv_classes[3][6][11] = int_stack + 460;
 Libderiv->deriv_classes[3][3][10] = int_stack + 740;
 Libderiv->deriv_classes[3][4][10] = int_stack + 840;
 Libderiv->deriv_classes[3][5][10] = int_stack + 990;
 Libderiv->deriv_classes[3][6][10] = int_stack + 1200;
 Libderiv->deriv_classes[3][3][9] = int_stack + 1480;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1580;
 Libderiv->deriv_classes[3][5][9] = int_stack + 1730;
 Libderiv->deriv_classes[3][6][9] = int_stack + 1940;
 Libderiv->deriv_classes[3][3][8] = int_stack + 2220;
 Libderiv->deriv_classes[3][4][8] = int_stack + 2320;
 Libderiv->deriv_classes[3][5][8] = int_stack + 2470;
 Libderiv->deriv_classes[3][6][8] = int_stack + 2680;
 Libderiv->deriv_classes[3][3][7] = int_stack + 2960;
 Libderiv->deriv_classes[3][4][7] = int_stack + 3060;
 Libderiv->deriv_classes[3][5][7] = int_stack + 3210;
 Libderiv->deriv_classes[3][6][7] = int_stack + 3420;
 Libderiv->dvrr_classes[3][3] = int_stack + 3700;
 Libderiv->deriv_classes[3][3][6] = int_stack + 3800;
 Libderiv->dvrr_classes[3][4] = int_stack + 3900;
 Libderiv->deriv_classes[3][4][6] = int_stack + 4050;
 Libderiv->dvrr_classes[3][5] = int_stack + 4200;
 Libderiv->deriv_classes[3][5][6] = int_stack + 4410;
 Libderiv->deriv_classes[3][6][6] = int_stack + 4620;
 Libderiv->deriv_classes[3][3][2] = int_stack + 4900;
 Libderiv->deriv_classes[3][4][2] = int_stack + 5000;
 Libderiv->deriv_classes[3][5][2] = int_stack + 5150;
 Libderiv->deriv_classes[3][6][2] = int_stack + 5360;
 Libderiv->deriv_classes[3][3][1] = int_stack + 5640;
 Libderiv->deriv_classes[3][4][1] = int_stack + 5740;
 Libderiv->deriv_classes[3][5][1] = int_stack + 5890;
 Libderiv->deriv_classes[3][6][1] = int_stack + 6100;
 Libderiv->deriv_classes[3][3][0] = int_stack + 6380;
 Libderiv->deriv_classes[3][4][0] = int_stack + 6480;
 Libderiv->deriv_classes[3][5][0] = int_stack + 6630;
 Libderiv->deriv_classes[3][6][0] = int_stack + 6840;
 memset(int_stack,0,56960);

 Libderiv->dvrr_stack = int_stack + 17780;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7120,int_stack+3900,int_stack+3700,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7420,int_stack+4200,int_stack+3900,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7870,int_stack+7420,int_stack+7120,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8470,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8770,int_stack+250,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3900,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9220,int_stack+8770,int_stack+8470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7120,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9820,int_stack+460,int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4200,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10450,int_stack+9820,int_stack+8770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7420,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9820,int_stack+840,int_stack+740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8470,int_stack+990,int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3900, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+8470,int_stack+9820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7120, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9820,int_stack+1200,int_stack+990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4200, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11350,int_stack+9820,int_stack+8470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7420, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8470,int_stack+1580,int_stack+1480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8770,int_stack+1730,int_stack+1580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3900, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9820,int_stack+8770,int_stack+8470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7120, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+600,int_stack+1940,int_stack+1730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1230,int_stack+600,int_stack+8770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7420, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+2320,int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8470,int_stack+2470,int_stack+2320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12250,int_stack+8470,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+600,int_stack+2680,int_stack+2470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12850,int_stack+600,int_stack+8470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8470,int_stack+3060,int_stack+2960, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8770,int_stack+3210,int_stack+3060, 0.0, zero_stack, 1.0, int_stack+3900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+8770,int_stack+8470, 0.0, zero_stack, 1.0, int_stack+7120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2130,int_stack+3420,int_stack+3210, 0.0, zero_stack, 1.0, int_stack+4200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2760,int_stack+2130,int_stack+8770, 0.0, zero_stack, 1.0, int_stack+7420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2130,int_stack+4050,int_stack+3800, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8470,int_stack+4410,int_stack+4050, 1.0, int_stack+3900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13750,int_stack+8470,int_stack+2130, 1.0, int_stack+7120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2130,int_stack+4620,int_stack+4410, 1.0, int_stack+4200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3660,int_stack+2130,int_stack+8470, 1.0, int_stack+7420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8470,int_stack+5000,int_stack+4900,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8770,int_stack+5150,int_stack+5000,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2130,int_stack+8770,int_stack+8470,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7120,int_stack+5360,int_stack+5150,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4560,int_stack+7120,int_stack+8770,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7120,int_stack+5740,int_stack+5640,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7420,int_stack+5890,int_stack+5740,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8470,int_stack+7420,int_stack+7120,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+14350,int_stack+6100,int_stack+5890,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5460,int_stack+14350,int_stack+7420,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14350,int_stack+6480,int_stack+6380,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7120,int_stack+6630,int_stack+6480,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14650,int_stack+7120,int_stack+14350,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15250,int_stack+6840,int_stack+6630,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15880,int_stack+15250,int_stack+7120,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6360,int_stack+10450,int_stack+9220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7870,10);
     Libderiv->ABCD[11] = int_stack + 6360;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16780,int_stack+11350,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7870, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 16780;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10420,int_stack+1230,int_stack+9820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7870, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 10420;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9070,int_stack+12850,int_stack+12250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 9070;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11420,int_stack+2760,int_stack+600, 0.0, zero_stack, 1.0, int_stack+7870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 11420;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+3660,int_stack+13750, 1.0, int_stack+7870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1000,int_stack+4560,int_stack+2130,10);
     Libderiv->ABCD[2] = int_stack + 1000;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+2000,int_stack+5460,int_stack+8470,10);
     Libderiv->ABCD[1] = int_stack + 2000;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+3000,int_stack+15880,int_stack+14650,10);
     Libderiv->ABCD[0] = int_stack + 3000;

}
