#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|fd) integrals */

void d1hrr_order_00fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][3][11] = int_stack + 0;
 Libderiv->deriv_classes[0][4][11] = int_stack + 10;
 Libderiv->deriv_classes[0][5][11] = int_stack + 25;
 Libderiv->deriv_classes[0][3][10] = int_stack + 46;
 Libderiv->deriv_classes[0][4][10] = int_stack + 56;
 Libderiv->deriv_classes[0][5][10] = int_stack + 71;
 Libderiv->deriv_classes[0][3][9] = int_stack + 92;
 Libderiv->deriv_classes[0][4][9] = int_stack + 102;
 Libderiv->deriv_classes[0][5][9] = int_stack + 117;
 Libderiv->deriv_classes[0][3][8] = int_stack + 138;
 Libderiv->deriv_classes[0][4][8] = int_stack + 148;
 Libderiv->deriv_classes[0][5][8] = int_stack + 163;
 Libderiv->deriv_classes[0][3][7] = int_stack + 184;
 Libderiv->deriv_classes[0][4][7] = int_stack + 194;
 Libderiv->deriv_classes[0][5][7] = int_stack + 209;
 Libderiv->dvrr_classes[0][3] = int_stack + 230;
 Libderiv->deriv_classes[0][3][6] = int_stack + 240;
 Libderiv->dvrr_classes[0][4] = int_stack + 250;
 Libderiv->deriv_classes[0][4][6] = int_stack + 265;
 Libderiv->deriv_classes[0][5][6] = int_stack + 280;
 Libderiv->deriv_classes[0][3][2] = int_stack + 301;
 Libderiv->deriv_classes[0][4][2] = int_stack + 311;
 Libderiv->deriv_classes[0][5][2] = int_stack + 326;
 Libderiv->deriv_classes[0][3][1] = int_stack + 347;
 Libderiv->deriv_classes[0][4][1] = int_stack + 357;
 Libderiv->deriv_classes[0][5][1] = int_stack + 372;
 Libderiv->deriv_classes[0][3][0] = int_stack + 393;
 Libderiv->deriv_classes[0][4][0] = int_stack + 403;
 Libderiv->deriv_classes[0][5][0] = int_stack + 418;
 memset(int_stack,0,3512);

 Libderiv->dvrr_stack = int_stack + 769;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+439,int_stack+250,int_stack+230,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+469,int_stack+10,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+499,int_stack+25,int_stack+10, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+250,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+56,int_stack+46, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+544,int_stack+71,int_stack+56, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+250, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30,int_stack+102,int_stack+92, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+589,int_stack+117,int_stack+102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+250, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60,int_stack+148,int_stack+138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+163,int_stack+148, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135,int_stack+194,int_stack+184, 0.0, zero_stack, 1.0, int_stack+230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+634,int_stack+209,int_stack+194, 0.0, zero_stack, 1.0, int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+165,int_stack+265,int_stack+240, 1.0, int_stack+230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195,int_stack+280,int_stack+265, 1.0, int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+240,int_stack+311,int_stack+301,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+679,int_stack+326,int_stack+311,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+357,int_stack+347,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+300,int_stack+372,int_stack+357,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+345,int_stack+403,int_stack+393,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+724,int_stack+418,int_stack+403,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+375,int_stack+499,int_stack+469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+439,1);
     Libderiv->ABCD[11] = int_stack + 375;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+469,int_stack+544,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+439, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 469;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+529,int_stack+589,int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+439, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 529;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+90,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+60,int_stack+634,int_stack+135, 0.0, zero_stack, 1.0, int_stack+439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 60;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+589,int_stack+195,int_stack+165, 1.0, int_stack+439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 589;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+120,int_stack+679,int_stack+240,1);
     Libderiv->ABCD[2] = int_stack + 120;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+180,int_stack+300,int_stack+270,1);
     Libderiv->ABCD[1] = int_stack + 180;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+240,int_stack+724,int_stack+345,1);
     Libderiv->ABCD[0] = int_stack + 240;

}
