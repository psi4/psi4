#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fdgp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|gp) integrals */

void d1hrr_order_fdgp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[4][4][11] = int_stack + 360;
 Libderiv->deriv_classes[4][5][11] = int_stack + 585;
 Libderiv->deriv_classes[5][4][11] = int_stack + 900;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1215;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1656;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1806;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2016;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2241;
 Libderiv->deriv_classes[5][4][10] = int_stack + 2556;
 Libderiv->deriv_classes[5][5][10] = int_stack + 2871;
 Libderiv->deriv_classes[3][4][9] = int_stack + 3312;
 Libderiv->deriv_classes[3][5][9] = int_stack + 3462;
 Libderiv->deriv_classes[4][4][9] = int_stack + 3672;
 Libderiv->deriv_classes[4][5][9] = int_stack + 3897;
 Libderiv->deriv_classes[5][4][9] = int_stack + 4212;
 Libderiv->deriv_classes[5][5][9] = int_stack + 4527;
 Libderiv->deriv_classes[3][4][8] = int_stack + 4968;
 Libderiv->deriv_classes[3][5][8] = int_stack + 5118;
 Libderiv->deriv_classes[4][4][8] = int_stack + 5328;
 Libderiv->deriv_classes[4][5][8] = int_stack + 5553;
 Libderiv->deriv_classes[5][4][8] = int_stack + 5868;
 Libderiv->deriv_classes[5][5][8] = int_stack + 6183;
 Libderiv->deriv_classes[3][4][7] = int_stack + 6624;
 Libderiv->deriv_classes[3][5][7] = int_stack + 6774;
 Libderiv->deriv_classes[4][4][7] = int_stack + 6984;
 Libderiv->deriv_classes[4][5][7] = int_stack + 7209;
 Libderiv->deriv_classes[5][4][7] = int_stack + 7524;
 Libderiv->deriv_classes[5][5][7] = int_stack + 7839;
 Libderiv->deriv_classes[3][4][6] = int_stack + 8280;
 Libderiv->deriv_classes[3][5][6] = int_stack + 8430;
 Libderiv->deriv_classes[4][4][6] = int_stack + 8640;
 Libderiv->deriv_classes[4][5][6] = int_stack + 8865;
 Libderiv->dvrr_classes[5][4] = int_stack + 9180;
 Libderiv->deriv_classes[5][4][6] = int_stack + 9495;
 Libderiv->deriv_classes[5][5][6] = int_stack + 9810;
 Libderiv->deriv_classes[3][4][2] = int_stack + 10251;
 Libderiv->deriv_classes[3][5][2] = int_stack + 10401;
 Libderiv->deriv_classes[4][4][2] = int_stack + 10611;
 Libderiv->deriv_classes[4][5][2] = int_stack + 10836;
 Libderiv->deriv_classes[5][4][2] = int_stack + 11151;
 Libderiv->deriv_classes[5][5][2] = int_stack + 11466;
 Libderiv->deriv_classes[3][4][1] = int_stack + 11907;
 Libderiv->deriv_classes[3][5][1] = int_stack + 12057;
 Libderiv->deriv_classes[4][4][1] = int_stack + 12267;
 Libderiv->deriv_classes[4][5][1] = int_stack + 12492;
 Libderiv->deriv_classes[5][4][1] = int_stack + 12807;
 Libderiv->deriv_classes[5][5][1] = int_stack + 13122;
 Libderiv->dvrr_classes[3][4] = int_stack + 13563;
 Libderiv->dvrr_classes[3][5] = int_stack + 13713;
 Libderiv->deriv_classes[3][4][0] = int_stack + 13923;
 Libderiv->deriv_classes[3][5][0] = int_stack + 14073;
 Libderiv->dvrr_classes[4][4] = int_stack + 14283;
 Libderiv->dvrr_classes[4][5] = int_stack + 14508;
 Libderiv->deriv_classes[4][4][0] = int_stack + 14823;
 Libderiv->deriv_classes[4][5][0] = int_stack + 15048;
 Libderiv->deriv_classes[5][4][0] = int_stack + 15363;
 Libderiv->deriv_classes[5][5][0] = int_stack + 15678;
 memset(int_stack,0,128952);

 Libderiv->dvrr_stack = int_stack + 37764;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fdgp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16119,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13563,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+585,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14283,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+17244,int_stack+16569,int_stack+16119,45);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+1215,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180,21);
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+19539,int_stack+18594,int_stack+16569,45);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+1806,int_stack+1656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13563, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16119,int_stack+2241,int_stack+2016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14283, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+16119,int_stack+18594,45);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+2871,int_stack+2556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack,21);
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+21564,int_stack+18594,int_stack+16119,45);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16119,int_stack+3462,int_stack+3312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13563, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+3897,int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14283, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1350,int_stack+16569,int_stack+16119,45);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+4527,int_stack+4212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+2700,int_stack+18594,int_stack+16569,45);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+5118,int_stack+4968, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13563, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16119,int_stack+5553,int_stack+5328, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+23589,int_stack+16119,int_stack+18594,45);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+6183,int_stack+5868, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+24939,int_stack+18594,int_stack+16119,45);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16119,int_stack+6774,int_stack+6624, 0.0, zero_stack, 1.0, int_stack+13563, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+7209,int_stack+6984, 0.0, zero_stack, 1.0, int_stack+14283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4725,int_stack+16569,int_stack+16119,45);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+7839,int_stack+7524, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+6075,int_stack+18594,int_stack+16569,45);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+8430,int_stack+8280, 1.0, int_stack+13563, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16119,int_stack+8865,int_stack+8640, 1.0, int_stack+14283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+26964,int_stack+16119,int_stack+18594,45);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+9810,int_stack+9495, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+8100,int_stack+18594,int_stack+16119,45);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16119,int_stack+13713,int_stack+13563,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+14508,int_stack+14283,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28314,int_stack+16569,int_stack+16119,45);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14283,int_stack+10401,int_stack+10251,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+10836,int_stack+10611,15);
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+29664,int_stack+18594,int_stack+14283, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16119, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10125,int_stack+11466,int_stack+11151,21);
 /*--- compute (gp|gp) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+31014,int_stack+10125,int_stack+18594, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+12057,int_stack+11907,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10125,int_stack+12492,int_stack+12267,15);
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+10800,int_stack+10125,int_stack+18594, 0.0, zero_stack, 1.0, int_stack+16119, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+13122,int_stack+12807,21);
 /*--- compute (gp|gp) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+33039,int_stack+18594,int_stack+10125, 0.0, zero_stack, 1.0, int_stack+16569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10125,int_stack+14073,int_stack+13923,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+15048,int_stack+14823,15);
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+12150,int_stack+18594,int_stack+10125, 1.0, int_stack+16119, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13500,int_stack+15678,int_stack+15363,21);
 /*--- compute (gp|gp) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+14445,int_stack+13500,int_stack+18594, 1.0, int_stack+16569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (fd|gp) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+35064,int_stack+19539,int_stack+17244,45);
     Libderiv->ABCD[11] = int_stack + 35064;
 /*--- compute (fd|gp) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+16470,int_stack+21564,int_stack+0,45);
     Libderiv->ABCD[10] = int_stack + 16470;
 /*--- compute (fd|gp) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+19170,int_stack+2700,int_stack+1350,45);
     Libderiv->ABCD[9] = int_stack + 19170;
 /*--- compute (fd|gp) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+0,int_stack+24939,int_stack+23589,45);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (fd|gp) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+21870,int_stack+6075,int_stack+4725,45);
     Libderiv->ABCD[7] = int_stack + 21870;
 /*--- compute (fd|gp) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+2700,int_stack+8100,int_stack+26964,45);
     Libderiv->ABCD[6] = int_stack + 2700;
 /*--- compute (fd|gp) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+5400,int_stack+31014,int_stack+29664, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[2] = int_stack + 5400;
 /*--- compute (fd|gp) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+8100,int_stack+33039,int_stack+10800, 0.0, zero_stack, 1.0, int_stack+28314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[1] = int_stack + 8100;
 /*--- compute (fd|gp) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+29664,int_stack+14445,int_stack+12150, 1.0, int_stack+28314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[0] = int_stack + 29664;

}
