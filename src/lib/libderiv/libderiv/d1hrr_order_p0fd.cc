#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|fd) integrals */

void d1hrr_order_p0fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][11] = int_stack + 30;
 Libderiv->deriv_classes[1][5][11] = int_stack + 75;
 Libderiv->deriv_classes[1][3][10] = int_stack + 138;
 Libderiv->deriv_classes[1][4][10] = int_stack + 168;
 Libderiv->deriv_classes[1][5][10] = int_stack + 213;
 Libderiv->deriv_classes[1][3][9] = int_stack + 276;
 Libderiv->deriv_classes[1][4][9] = int_stack + 306;
 Libderiv->deriv_classes[1][5][9] = int_stack + 351;
 Libderiv->deriv_classes[1][3][8] = int_stack + 414;
 Libderiv->deriv_classes[1][4][8] = int_stack + 444;
 Libderiv->deriv_classes[1][5][8] = int_stack + 489;
 Libderiv->deriv_classes[1][3][7] = int_stack + 552;
 Libderiv->deriv_classes[1][4][7] = int_stack + 582;
 Libderiv->deriv_classes[1][5][7] = int_stack + 627;
 Libderiv->dvrr_classes[1][3] = int_stack + 690;
 Libderiv->deriv_classes[1][3][6] = int_stack + 720;
 Libderiv->dvrr_classes[1][4] = int_stack + 750;
 Libderiv->deriv_classes[1][4][6] = int_stack + 795;
 Libderiv->deriv_classes[1][5][6] = int_stack + 840;
 Libderiv->deriv_classes[1][3][2] = int_stack + 903;
 Libderiv->deriv_classes[1][4][2] = int_stack + 933;
 Libderiv->deriv_classes[1][5][2] = int_stack + 978;
 Libderiv->deriv_classes[1][3][1] = int_stack + 1041;
 Libderiv->deriv_classes[1][4][1] = int_stack + 1071;
 Libderiv->deriv_classes[1][5][1] = int_stack + 1116;
 Libderiv->deriv_classes[1][3][0] = int_stack + 1179;
 Libderiv->deriv_classes[1][4][0] = int_stack + 1209;
 Libderiv->deriv_classes[1][5][0] = int_stack + 1254;
 memset(int_stack,0,10536);

 Libderiv->dvrr_stack = int_stack + 2307;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1317,int_stack+750,int_stack+690,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1407,int_stack+30,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+690,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1497,int_stack+75,int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+168,int_stack+138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+690, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1632,int_stack+213,int_stack+168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+306,int_stack+276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+690, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1767,int_stack+351,int_stack+306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+444,int_stack+414, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+489,int_stack+444, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+405,int_stack+582,int_stack+552, 0.0, zero_stack, 1.0, int_stack+690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1902,int_stack+627,int_stack+582, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+495,int_stack+795,int_stack+720, 1.0, int_stack+690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+585,int_stack+840,int_stack+795, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+933,int_stack+903,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2037,int_stack+978,int_stack+933,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+810,int_stack+1071,int_stack+1041,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+1116,int_stack+1071,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1035,int_stack+1209,int_stack+1179,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2172,int_stack+1254,int_stack+1209,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+1497,int_stack+1407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1317,3);
     Libderiv->ABCD[11] = int_stack + 1125;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1407,int_stack+1632,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1317, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 1407;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1587,int_stack+1767,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1317, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 1587;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+270,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+180,int_stack+1902,int_stack+405, 0.0, zero_stack, 1.0, int_stack+1317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 180;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1767,int_stack+585,int_stack+495, 1.0, int_stack+1317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 1767;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+2037,int_stack+720,3);
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+540,int_stack+900,int_stack+810,3);
     Libderiv->ABCD[1] = int_stack + 540;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+720,int_stack+2172,int_stack+1035,3);
     Libderiv->ABCD[0] = int_stack + 720;

}
