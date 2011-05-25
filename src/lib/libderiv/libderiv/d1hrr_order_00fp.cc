#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|fp) integrals */

void d1hrr_order_00fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][3][11] = int_stack + 0;
 Libderiv->deriv_classes[0][4][11] = int_stack + 10;
 Libderiv->deriv_classes[0][3][10] = int_stack + 25;
 Libderiv->deriv_classes[0][4][10] = int_stack + 35;
 Libderiv->deriv_classes[0][3][9] = int_stack + 50;
 Libderiv->deriv_classes[0][4][9] = int_stack + 60;
 Libderiv->deriv_classes[0][3][8] = int_stack + 75;
 Libderiv->deriv_classes[0][4][8] = int_stack + 85;
 Libderiv->deriv_classes[0][3][7] = int_stack + 100;
 Libderiv->deriv_classes[0][4][7] = int_stack + 110;
 Libderiv->dvrr_classes[0][3] = int_stack + 125;
 Libderiv->deriv_classes[0][3][6] = int_stack + 135;
 Libderiv->deriv_classes[0][4][6] = int_stack + 145;
 Libderiv->deriv_classes[0][3][2] = int_stack + 160;
 Libderiv->deriv_classes[0][4][2] = int_stack + 170;
 Libderiv->deriv_classes[0][3][1] = int_stack + 185;
 Libderiv->deriv_classes[0][4][1] = int_stack + 195;
 Libderiv->deriv_classes[0][3][0] = int_stack + 210;
 Libderiv->deriv_classes[0][4][0] = int_stack + 220;
 memset(int_stack,0,1880);

 Libderiv->dvrr_stack = int_stack + 295;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+235,int_stack+10,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125,1);
     Libderiv->ABCD[11] = int_stack + 235;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+265,int_stack+35,int_stack+25, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 265;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+60,int_stack+50, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30,int_stack+85,int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 30;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60,int_stack+110,int_stack+100, 0.0, zero_stack, 1.0, int_stack+125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 60;
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+145,int_stack+135, 1.0, int_stack+125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 90;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+120,int_stack+170,int_stack+160,1);
     Libderiv->ABCD[2] = int_stack + 120;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+150,int_stack+195,int_stack+185,1);
     Libderiv->ABCD[1] = int_stack + 150;
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+220,int_stack+210,1);
     Libderiv->ABCD[0] = int_stack + 180;

}
