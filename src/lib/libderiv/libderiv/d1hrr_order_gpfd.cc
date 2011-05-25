#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gpfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gp|fd) integrals */

void d1hrr_order_gpfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][3][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][11] = int_stack + 150;
 Libderiv->deriv_classes[4][5][11] = int_stack + 375;
 Libderiv->deriv_classes[5][3][11] = int_stack + 690;
 Libderiv->deriv_classes[5][4][11] = int_stack + 900;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1215;
 Libderiv->deriv_classes[4][3][10] = int_stack + 1656;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1806;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2031;
 Libderiv->deriv_classes[5][3][10] = int_stack + 2346;
 Libderiv->deriv_classes[5][4][10] = int_stack + 2556;
 Libderiv->deriv_classes[5][5][10] = int_stack + 2871;
 Libderiv->deriv_classes[4][3][9] = int_stack + 3312;
 Libderiv->deriv_classes[4][4][9] = int_stack + 3462;
 Libderiv->deriv_classes[4][5][9] = int_stack + 3687;
 Libderiv->deriv_classes[5][3][9] = int_stack + 4002;
 Libderiv->deriv_classes[5][4][9] = int_stack + 4212;
 Libderiv->deriv_classes[5][5][9] = int_stack + 4527;
 Libderiv->deriv_classes[4][3][8] = int_stack + 4968;
 Libderiv->deriv_classes[4][4][8] = int_stack + 5118;
 Libderiv->deriv_classes[4][5][8] = int_stack + 5343;
 Libderiv->deriv_classes[5][3][8] = int_stack + 5658;
 Libderiv->deriv_classes[5][4][8] = int_stack + 5868;
 Libderiv->deriv_classes[5][5][8] = int_stack + 6183;
 Libderiv->deriv_classes[4][3][7] = int_stack + 6624;
 Libderiv->deriv_classes[4][4][7] = int_stack + 6774;
 Libderiv->deriv_classes[4][5][7] = int_stack + 6999;
 Libderiv->deriv_classes[5][3][7] = int_stack + 7314;
 Libderiv->deriv_classes[5][4][7] = int_stack + 7524;
 Libderiv->deriv_classes[5][5][7] = int_stack + 7839;
 Libderiv->deriv_classes[4][3][6] = int_stack + 8280;
 Libderiv->deriv_classes[4][4][6] = int_stack + 8430;
 Libderiv->deriv_classes[4][5][6] = int_stack + 8655;
 Libderiv->dvrr_classes[5][3] = int_stack + 8970;
 Libderiv->deriv_classes[5][3][6] = int_stack + 9180;
 Libderiv->dvrr_classes[5][4] = int_stack + 9390;
 Libderiv->deriv_classes[5][4][6] = int_stack + 9705;
 Libderiv->deriv_classes[5][5][6] = int_stack + 10020;
 Libderiv->deriv_classes[4][3][2] = int_stack + 10461;
 Libderiv->deriv_classes[4][4][2] = int_stack + 10611;
 Libderiv->deriv_classes[4][5][2] = int_stack + 10836;
 Libderiv->deriv_classes[5][3][2] = int_stack + 11151;
 Libderiv->deriv_classes[5][4][2] = int_stack + 11361;
 Libderiv->deriv_classes[5][5][2] = int_stack + 11676;
 Libderiv->deriv_classes[4][3][1] = int_stack + 12117;
 Libderiv->deriv_classes[4][4][1] = int_stack + 12267;
 Libderiv->deriv_classes[4][5][1] = int_stack + 12492;
 Libderiv->deriv_classes[5][3][1] = int_stack + 12807;
 Libderiv->deriv_classes[5][4][1] = int_stack + 13017;
 Libderiv->deriv_classes[5][5][1] = int_stack + 13332;
 Libderiv->dvrr_classes[4][3] = int_stack + 13773;
 Libderiv->dvrr_classes[4][4] = int_stack + 13923;
 Libderiv->dvrr_classes[4][5] = int_stack + 14148;
 Libderiv->deriv_classes[4][3][0] = int_stack + 14463;
 Libderiv->deriv_classes[4][4][0] = int_stack + 14613;
 Libderiv->deriv_classes[4][5][0] = int_stack + 14838;
 Libderiv->deriv_classes[5][3][0] = int_stack + 15153;
 Libderiv->deriv_classes[5][4][0] = int_stack + 15363;
 Libderiv->deriv_classes[5][5][0] = int_stack + 15678;
 memset(int_stack,0,128952);

 Libderiv->dvrr_stack = int_stack + 31059;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gpfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16119,int_stack+13923,int_stack+13773,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16569,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13773,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17019,int_stack+375,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13923,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17694,int_stack+17019,int_stack+16569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16119,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16569,int_stack+9390,int_stack+8970,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+900,int_stack+690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8970,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+1215,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9390,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19539,int_stack+18594,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16569,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1806,int_stack+1656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13773, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+2031,int_stack+1806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13923, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16119, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2556,int_stack+2346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8970, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18594,int_stack+2871,int_stack+2556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9390, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2025,int_stack+18594,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16569, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+3462,int_stack+3312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13773, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+3687,int_stack+3462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13923, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18594,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16119, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+4212,int_stack+4002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8970, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20799,int_stack+4527,int_stack+4212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9390, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3285,int_stack+20799,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16569, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+5118,int_stack+4968, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+5343,int_stack+5118, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20799,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16119, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+5868,int_stack+5658, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4545,int_stack+6183,int_stack+5868, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21699,int_stack+4545,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+6774,int_stack+6624, 0.0, zero_stack, 1.0, int_stack+13773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+6999,int_stack+6774, 0.0, zero_stack, 1.0, int_stack+13923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4545,int_stack+450,int_stack+0, 0.0, zero_stack, 1.0, int_stack+16119, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+7524,int_stack+7314, 0.0, zero_stack, 1.0, int_stack+8970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5445,int_stack+7839,int_stack+7524, 0.0, zero_stack, 1.0, int_stack+9390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6390,int_stack+5445,int_stack+0, 0.0, zero_stack, 1.0, int_stack+16569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+8430,int_stack+8280, 1.0, int_stack+13773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+8655,int_stack+8430, 1.0, int_stack+13923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5445,int_stack+450,int_stack+0, 1.0, int_stack+16119, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9705,int_stack+9180, 1.0, int_stack+8970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7650,int_stack+10020,int_stack+9705, 1.0, int_stack+9390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8595,int_stack+7650,int_stack+0, 1.0, int_stack+16569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+14148,int_stack+13923,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+16569,int_stack+16119,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16119,int_stack+10611,int_stack+10461,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+10836,int_stack+10611,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7650,int_stack+16569,int_stack+16119,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16119,int_stack+11361,int_stack+11151,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16749,int_stack+11676,int_stack+11361,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+9855,int_stack+16749,int_stack+16119,21);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16119,int_stack+12267,int_stack+12117,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+12492,int_stack+12267,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11115,int_stack+16569,int_stack+16119,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16119,int_stack+13017,int_stack+12807,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16749,int_stack+13332,int_stack+13017,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12015,int_stack+16749,int_stack+16119,21);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16119,int_stack+14613,int_stack+14463,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16569,int_stack+14838,int_stack+14613,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+13275,int_stack+16569,int_stack+16119,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16119,int_stack+15363,int_stack+15153,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16749,int_stack+15678,int_stack+15363,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14175,int_stack+16749,int_stack+16119,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+22959,int_stack+19539,int_stack+17694,60);
     Libderiv->ABCD[11] = int_stack + 22959;
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+15435,int_stack+2025,int_stack+1125,60);
     Libderiv->ABCD[10] = int_stack + 15435;
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+25659,int_stack+3285,int_stack+18594,60);
     Libderiv->ABCD[9] = int_stack + 25659;
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+900,int_stack+21699,int_stack+20799,60);
     Libderiv->ABCD[8] = int_stack + 900;
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+18135,int_stack+6390,int_stack+4545,60);
     Libderiv->ABCD[7] = int_stack + 18135;
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+28359,int_stack+8595,int_stack+5445,60);
     Libderiv->ABCD[6] = int_stack + 28359;
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+3600,int_stack+9855,int_stack+7650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 3600;
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+6300,int_stack+12015,int_stack+11115, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 6300;
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+9000,int_stack+14175,int_stack+13275, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 9000;

}
