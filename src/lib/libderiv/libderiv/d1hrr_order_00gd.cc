#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00gd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|gd) integrals */

void d1hrr_order_00gd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][4][11] = int_stack + 0;
 Libderiv->deriv_classes[0][5][11] = int_stack + 15;
 Libderiv->deriv_classes[0][6][11] = int_stack + 36;
 Libderiv->deriv_classes[0][4][10] = int_stack + 64;
 Libderiv->deriv_classes[0][5][10] = int_stack + 79;
 Libderiv->deriv_classes[0][6][10] = int_stack + 100;
 Libderiv->deriv_classes[0][4][9] = int_stack + 128;
 Libderiv->deriv_classes[0][5][9] = int_stack + 143;
 Libderiv->deriv_classes[0][6][9] = int_stack + 164;
 Libderiv->deriv_classes[0][4][8] = int_stack + 192;
 Libderiv->deriv_classes[0][5][8] = int_stack + 207;
 Libderiv->deriv_classes[0][6][8] = int_stack + 228;
 Libderiv->deriv_classes[0][4][7] = int_stack + 256;
 Libderiv->deriv_classes[0][5][7] = int_stack + 271;
 Libderiv->deriv_classes[0][6][7] = int_stack + 292;
 Libderiv->dvrr_classes[0][4] = int_stack + 320;
 Libderiv->deriv_classes[0][4][6] = int_stack + 335;
 Libderiv->dvrr_classes[0][5] = int_stack + 350;
 Libderiv->deriv_classes[0][5][6] = int_stack + 371;
 Libderiv->deriv_classes[0][6][6] = int_stack + 392;
 Libderiv->deriv_classes[0][4][2] = int_stack + 420;
 Libderiv->deriv_classes[0][5][2] = int_stack + 435;
 Libderiv->deriv_classes[0][6][2] = int_stack + 456;
 Libderiv->deriv_classes[0][4][1] = int_stack + 484;
 Libderiv->deriv_classes[0][5][1] = int_stack + 499;
 Libderiv->deriv_classes[0][6][1] = int_stack + 520;
 Libderiv->deriv_classes[0][4][0] = int_stack + 548;
 Libderiv->deriv_classes[0][5][0] = int_stack + 563;
 Libderiv->deriv_classes[0][6][0] = int_stack + 584;
 memset(int_stack,0,4896);

 Libderiv->dvrr_stack = int_stack + 1260;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00gd(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+612,int_stack+350,int_stack+320,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+657,int_stack+15,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+320,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+702,int_stack+36,int_stack+15, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+350,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+79,int_stack+64, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+320, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+765,int_stack+100,int_stack+79, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+350, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45,int_stack+143,int_stack+128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+320, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+828,int_stack+164,int_stack+143, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+350, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+207,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+135,int_stack+228,int_stack+207, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+198,int_stack+271,int_stack+256, 0.0, zero_stack, 1.0, int_stack+320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+891,int_stack+292,int_stack+271, 0.0, zero_stack, 1.0, int_stack+350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+243,int_stack+371,int_stack+335, 1.0, int_stack+320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+954,int_stack+392,int_stack+371, 1.0, int_stack+350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+288,int_stack+435,int_stack+420,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+333,int_stack+456,int_stack+435,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+396,int_stack+499,int_stack+484,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1017,int_stack+520,int_stack+499,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+441,int_stack+563,int_stack+548,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+486,int_stack+584,int_stack+563,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+702,int_stack+657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+612,1);
     Libderiv->ABCD[11] = int_stack + 1080;
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+657,int_stack+765,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+612, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 657;
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1170,int_stack+828,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+612, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 1170;
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+135,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+90,int_stack+891,int_stack+198, 0.0, zero_stack, 1.0, int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 90;
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+747,int_stack+954,int_stack+243, 1.0, int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 747;
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+180,int_stack+333,int_stack+288,1);
     Libderiv->ABCD[2] = int_stack + 180;
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+270,int_stack+1017,int_stack+396,1);
     Libderiv->ABCD[1] = int_stack + 270;
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+837,int_stack+486,int_stack+441,1);
     Libderiv->ABCD[0] = int_stack + 837;

}
