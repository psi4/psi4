#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppg0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|g0) integrals */

void d1hrr_order_ppg0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 45;
 Libderiv->deriv_classes[1][4][10] = int_stack + 135;
 Libderiv->deriv_classes[2][4][10] = int_stack + 180;
 Libderiv->deriv_classes[1][4][9] = int_stack + 270;
 Libderiv->deriv_classes[2][4][9] = int_stack + 315;
 Libderiv->deriv_classes[1][4][8] = int_stack + 405;
 Libderiv->deriv_classes[2][4][8] = int_stack + 450;
 Libderiv->deriv_classes[1][4][7] = int_stack + 540;
 Libderiv->deriv_classes[2][4][7] = int_stack + 585;
 Libderiv->deriv_classes[1][4][6] = int_stack + 675;
 Libderiv->deriv_classes[2][4][6] = int_stack + 720;
 Libderiv->deriv_classes[1][4][2] = int_stack + 810;
 Libderiv->deriv_classes[2][4][2] = int_stack + 855;
 Libderiv->deriv_classes[1][4][1] = int_stack + 945;
 Libderiv->deriv_classes[2][4][1] = int_stack + 990;
 Libderiv->dvrr_classes[1][4] = int_stack + 1080;
 Libderiv->deriv_classes[1][4][0] = int_stack + 1125;
 Libderiv->deriv_classes[2][4][0] = int_stack + 1170;
 memset(int_stack,0,10080);

 Libderiv->dvrr_stack = int_stack + 1395;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppg0(Libderiv, Data);
   Data++;
 }

 /*--- compute (pp|g0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1260,int_stack+45,int_stack+0,15);
     Libderiv->ABCD[11] = int_stack + 1260;
 /*--- compute (pp|g0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+180,int_stack+135,15);
     Libderiv->ABCD[10] = int_stack + 0;
 /*--- compute (pp|g0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+135,int_stack+315,int_stack+270,15);
     Libderiv->ABCD[9] = int_stack + 135;
 /*--- compute (pp|g0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+270,int_stack+450,int_stack+405,15);
     Libderiv->ABCD[8] = int_stack + 270;
 /*--- compute (pp|g0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+405,int_stack+585,int_stack+540,15);
     Libderiv->ABCD[7] = int_stack + 405;
 /*--- compute (pp|g0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+540,int_stack+720,int_stack+675,15);
     Libderiv->ABCD[6] = int_stack + 540;
 /*--- compute (pp|g0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+675,int_stack+855,int_stack+810, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[2] = int_stack + 675;
 /*--- compute (pp|g0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+810,int_stack+990,int_stack+945, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[1] = int_stack + 810;
 /*--- compute (pp|g0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+945,int_stack+1170,int_stack+1125, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[0] = int_stack + 945;

}
