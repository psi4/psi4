#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpf0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|f0) integrals */

void d1hrr_order_dpf0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][3][11] = int_stack + 60;
 Libderiv->deriv_classes[2][3][10] = int_stack + 160;
 Libderiv->deriv_classes[3][3][10] = int_stack + 220;
 Libderiv->deriv_classes[2][3][9] = int_stack + 320;
 Libderiv->deriv_classes[3][3][9] = int_stack + 380;
 Libderiv->deriv_classes[2][3][8] = int_stack + 480;
 Libderiv->deriv_classes[3][3][8] = int_stack + 540;
 Libderiv->deriv_classes[2][3][7] = int_stack + 640;
 Libderiv->deriv_classes[3][3][7] = int_stack + 700;
 Libderiv->deriv_classes[2][3][6] = int_stack + 800;
 Libderiv->deriv_classes[3][3][6] = int_stack + 860;
 Libderiv->deriv_classes[2][3][2] = int_stack + 960;
 Libderiv->deriv_classes[3][3][2] = int_stack + 1020;
 Libderiv->deriv_classes[2][3][1] = int_stack + 1120;
 Libderiv->deriv_classes[3][3][1] = int_stack + 1180;
 Libderiv->dvrr_classes[2][3] = int_stack + 1280;
 Libderiv->deriv_classes[2][3][0] = int_stack + 1340;
 Libderiv->deriv_classes[3][3][0] = int_stack + 1400;
 memset(int_stack,0,12000);

 Libderiv->dvrr_stack = int_stack + 1860;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpf0(Libderiv, Data);
   Data++;
 }

 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1500,int_stack+60,int_stack+0,10);
     Libderiv->ABCD[11] = int_stack + 1500;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1680,int_stack+220,int_stack+160,10);
     Libderiv->ABCD[10] = int_stack + 1680;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+380,int_stack+320,10);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+180,int_stack+540,int_stack+480,10);
     Libderiv->ABCD[8] = int_stack + 180;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+360,int_stack+700,int_stack+640,10);
     Libderiv->ABCD[7] = int_stack + 360;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+540,int_stack+860,int_stack+800,10);
     Libderiv->ABCD[6] = int_stack + 540;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+720,int_stack+1020,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[2] = int_stack + 720;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+900,int_stack+1180,int_stack+1120, 0.0, zero_stack, 1.0, int_stack+1280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[1] = int_stack + 900;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1080,int_stack+1400,int_stack+1340, 1.0, int_stack+1280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[0] = int_stack + 1080;

}
