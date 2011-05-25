#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gpgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gp|gd) integrals */

void d1hrr_order_gpgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[5][4][11] = int_stack + 960;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1275;
 Libderiv->deriv_classes[5][6][11] = int_stack + 1716;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2304;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2529;
 Libderiv->deriv_classes[4][6][10] = int_stack + 2844;
 Libderiv->deriv_classes[5][4][10] = int_stack + 3264;
 Libderiv->deriv_classes[5][5][10] = int_stack + 3579;
 Libderiv->deriv_classes[5][6][10] = int_stack + 4020;
 Libderiv->deriv_classes[4][4][9] = int_stack + 4608;
 Libderiv->deriv_classes[4][5][9] = int_stack + 4833;
 Libderiv->deriv_classes[4][6][9] = int_stack + 5148;
 Libderiv->deriv_classes[5][4][9] = int_stack + 5568;
 Libderiv->deriv_classes[5][5][9] = int_stack + 5883;
 Libderiv->deriv_classes[5][6][9] = int_stack + 6324;
 Libderiv->deriv_classes[4][4][8] = int_stack + 6912;
 Libderiv->deriv_classes[4][5][8] = int_stack + 7137;
 Libderiv->deriv_classes[4][6][8] = int_stack + 7452;
 Libderiv->deriv_classes[5][4][8] = int_stack + 7872;
 Libderiv->deriv_classes[5][5][8] = int_stack + 8187;
 Libderiv->deriv_classes[5][6][8] = int_stack + 8628;
 Libderiv->deriv_classes[4][4][7] = int_stack + 9216;
 Libderiv->deriv_classes[4][5][7] = int_stack + 9441;
 Libderiv->deriv_classes[4][6][7] = int_stack + 9756;
 Libderiv->deriv_classes[5][4][7] = int_stack + 10176;
 Libderiv->deriv_classes[5][5][7] = int_stack + 10491;
 Libderiv->deriv_classes[5][6][7] = int_stack + 10932;
 Libderiv->deriv_classes[4][4][6] = int_stack + 11520;
 Libderiv->deriv_classes[4][5][6] = int_stack + 11745;
 Libderiv->deriv_classes[4][6][6] = int_stack + 12060;
 Libderiv->dvrr_classes[5][4] = int_stack + 12480;
 Libderiv->deriv_classes[5][4][6] = int_stack + 12795;
 Libderiv->dvrr_classes[5][5] = int_stack + 13110;
 Libderiv->deriv_classes[5][5][6] = int_stack + 13551;
 Libderiv->deriv_classes[5][6][6] = int_stack + 13992;
 Libderiv->deriv_classes[4][4][2] = int_stack + 14580;
 Libderiv->deriv_classes[4][5][2] = int_stack + 14805;
 Libderiv->deriv_classes[4][6][2] = int_stack + 15120;
 Libderiv->deriv_classes[5][4][2] = int_stack + 15540;
 Libderiv->deriv_classes[5][5][2] = int_stack + 15855;
 Libderiv->deriv_classes[5][6][2] = int_stack + 16296;
 Libderiv->deriv_classes[4][4][1] = int_stack + 16884;
 Libderiv->deriv_classes[4][5][1] = int_stack + 17109;
 Libderiv->deriv_classes[4][6][1] = int_stack + 17424;
 Libderiv->deriv_classes[5][4][1] = int_stack + 17844;
 Libderiv->deriv_classes[5][5][1] = int_stack + 18159;
 Libderiv->deriv_classes[5][6][1] = int_stack + 18600;
 Libderiv->dvrr_classes[4][4] = int_stack + 19188;
 Libderiv->dvrr_classes[4][5] = int_stack + 19413;
 Libderiv->dvrr_classes[4][6] = int_stack + 19728;
 Libderiv->deriv_classes[4][4][0] = int_stack + 20148;
 Libderiv->deriv_classes[4][5][0] = int_stack + 20373;
 Libderiv->deriv_classes[4][6][0] = int_stack + 20688;
 Libderiv->deriv_classes[5][4][0] = int_stack + 21108;
 Libderiv->deriv_classes[5][5][0] = int_stack + 21423;
 Libderiv->deriv_classes[5][6][0] = int_stack + 21864;
 memset(int_stack,0,179616);

 Libderiv->dvrr_stack = int_stack + 46590;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gpgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22452,int_stack+19413,int_stack+19188,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23127,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19188,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23802,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19413,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24747,int_stack+23802,int_stack+23127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22452,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23127,int_stack+13110,int_stack+12480,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1275,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12480,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+1716,int_stack+1275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13110,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27420,int_stack+26097,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23127,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24072,int_stack+2529,int_stack+2304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19188, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+2844,int_stack+2529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19413, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+945,int_stack+0,int_stack+24072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22452, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3579,int_stack+3264, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12480, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+4020,int_stack+3579, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13110, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2295,int_stack+26097,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23127, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24072,int_stack+4833,int_stack+4608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19188, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+5148,int_stack+4833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19413, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4185,int_stack+0,int_stack+24072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22452, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+5883,int_stack+5568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12480, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+6324,int_stack+5883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13110, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29310,int_stack+26097,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23127, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24072,int_stack+7137,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+7452,int_stack+7137, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19413, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5535,int_stack+0,int_stack+24072, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+8187,int_stack+7872, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+8628,int_stack+8187, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6885,int_stack+26097,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24072,int_stack+9441,int_stack+9216, 0.0, zero_stack, 1.0, int_stack+19188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+9756,int_stack+9441, 0.0, zero_stack, 1.0, int_stack+19413, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8775,int_stack+0,int_stack+24072, 0.0, zero_stack, 1.0, int_stack+22452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+10491,int_stack+10176, 0.0, zero_stack, 1.0, int_stack+12480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+10932,int_stack+10491, 0.0, zero_stack, 1.0, int_stack+13110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31200,int_stack+26097,int_stack+0, 0.0, zero_stack, 1.0, int_stack+23127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24072,int_stack+11745,int_stack+11520, 1.0, int_stack+19188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+12060,int_stack+11745, 1.0, int_stack+19413, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10125,int_stack+0,int_stack+24072, 1.0, int_stack+22452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+13551,int_stack+12795, 1.0, int_stack+12480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+13992,int_stack+13551, 1.0, int_stack+13110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11475,int_stack+26097,int_stack+0, 1.0, int_stack+23127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+19728,int_stack+19413,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23127,int_stack+0,int_stack+22452,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22452,int_stack+14805,int_stack+14580,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+15120,int_stack+14805,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+13365,int_stack+0,int_stack+22452,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+15855,int_stack+15540,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+16296,int_stack+15855,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+14715,int_stack+26097,int_stack+0,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22452,int_stack+17109,int_stack+16884,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+17424,int_stack+17109,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+33090,int_stack+0,int_stack+22452,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+18159,int_stack+17844,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+18600,int_stack+18159,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16605,int_stack+26097,int_stack+0,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22452,int_stack+20373,int_stack+20148,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+20688,int_stack+20373,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18495,int_stack+0,int_stack+22452,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+21423,int_stack+21108,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+26097,int_stack+21864,int_stack+21423,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+19845,int_stack+26097,int_stack+0,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+34440,int_stack+27420,int_stack+24747,90);
     Libderiv->ABCD[11] = int_stack + 34440;
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+24477,int_stack+2295,int_stack+945,90);
     Libderiv->ABCD[10] = int_stack + 24477;
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+29310,int_stack+4185,90);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+38490,int_stack+6885,int_stack+5535,90);
     Libderiv->ABCD[8] = int_stack + 38490;
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+4050,int_stack+31200,int_stack+8775,90);
     Libderiv->ABCD[7] = int_stack + 4050;
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+28527,int_stack+11475,int_stack+10125,90);
     Libderiv->ABCD[6] = int_stack + 28527;
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+8100,int_stack+14715,int_stack+13365, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 8100;
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+12150,int_stack+16605,int_stack+33090, 0.0, zero_stack, 1.0, int_stack+23127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 12150;
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+42540,int_stack+19845,int_stack+18495, 1.0, int_stack+23127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 42540;

}
