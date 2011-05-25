#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|fd) integrals */

void d1hrr_order_fpfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 100;
 Libderiv->deriv_classes[3][5][11] = int_stack + 250;
 Libderiv->deriv_classes[4][3][11] = int_stack + 460;
 Libderiv->deriv_classes[4][4][11] = int_stack + 610;
 Libderiv->deriv_classes[4][5][11] = int_stack + 835;
 Libderiv->deriv_classes[3][3][10] = int_stack + 1150;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1250;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1400;
 Libderiv->deriv_classes[4][3][10] = int_stack + 1610;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1760;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1985;
 Libderiv->deriv_classes[3][3][9] = int_stack + 2300;
 Libderiv->deriv_classes[3][4][9] = int_stack + 2400;
 Libderiv->deriv_classes[3][5][9] = int_stack + 2550;
 Libderiv->deriv_classes[4][3][9] = int_stack + 2760;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2910;
 Libderiv->deriv_classes[4][5][9] = int_stack + 3135;
 Libderiv->deriv_classes[3][3][8] = int_stack + 3450;
 Libderiv->deriv_classes[3][4][8] = int_stack + 3550;
 Libderiv->deriv_classes[3][5][8] = int_stack + 3700;
 Libderiv->deriv_classes[4][3][8] = int_stack + 3910;
 Libderiv->deriv_classes[4][4][8] = int_stack + 4060;
 Libderiv->deriv_classes[4][5][8] = int_stack + 4285;
 Libderiv->deriv_classes[3][3][7] = int_stack + 4600;
 Libderiv->deriv_classes[3][4][7] = int_stack + 4700;
 Libderiv->deriv_classes[3][5][7] = int_stack + 4850;
 Libderiv->deriv_classes[4][3][7] = int_stack + 5060;
 Libderiv->deriv_classes[4][4][7] = int_stack + 5210;
 Libderiv->deriv_classes[4][5][7] = int_stack + 5435;
 Libderiv->deriv_classes[3][3][6] = int_stack + 5750;
 Libderiv->deriv_classes[3][4][6] = int_stack + 5850;
 Libderiv->deriv_classes[3][5][6] = int_stack + 6000;
 Libderiv->dvrr_classes[4][3] = int_stack + 6210;
 Libderiv->deriv_classes[4][3][6] = int_stack + 6360;
 Libderiv->dvrr_classes[4][4] = int_stack + 6510;
 Libderiv->deriv_classes[4][4][6] = int_stack + 6735;
 Libderiv->deriv_classes[4][5][6] = int_stack + 6960;
 Libderiv->deriv_classes[3][3][2] = int_stack + 7275;
 Libderiv->deriv_classes[3][4][2] = int_stack + 7375;
 Libderiv->deriv_classes[3][5][2] = int_stack + 7525;
 Libderiv->deriv_classes[4][3][2] = int_stack + 7735;
 Libderiv->deriv_classes[4][4][2] = int_stack + 7885;
 Libderiv->deriv_classes[4][5][2] = int_stack + 8110;
 Libderiv->deriv_classes[3][3][1] = int_stack + 8425;
 Libderiv->deriv_classes[3][4][1] = int_stack + 8525;
 Libderiv->deriv_classes[3][5][1] = int_stack + 8675;
 Libderiv->deriv_classes[4][3][1] = int_stack + 8885;
 Libderiv->deriv_classes[4][4][1] = int_stack + 9035;
 Libderiv->deriv_classes[4][5][1] = int_stack + 9260;
 Libderiv->dvrr_classes[3][3] = int_stack + 9575;
 Libderiv->dvrr_classes[3][4] = int_stack + 9675;
 Libderiv->dvrr_classes[3][5] = int_stack + 9825;
 Libderiv->deriv_classes[3][3][0] = int_stack + 10035;
 Libderiv->deriv_classes[3][4][0] = int_stack + 10135;
 Libderiv->deriv_classes[3][5][0] = int_stack + 10285;
 Libderiv->deriv_classes[4][3][0] = int_stack + 10495;
 Libderiv->deriv_classes[4][4][0] = int_stack + 10645;
 Libderiv->deriv_classes[4][5][0] = int_stack + 10870;
 memset(int_stack,0,89480);

 Libderiv->dvrr_stack = int_stack + 21310;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11185,int_stack+9675,int_stack+9575,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11485,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9575,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11785,int_stack+250,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9675,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12235,int_stack+11785,int_stack+11485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11185,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11485,int_stack+6510,int_stack+6210,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+610,int_stack+460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6210,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12835,int_stack+835,int_stack+610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6510,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13510,int_stack+12835,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11485,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11935,int_stack+1250,int_stack+1150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9575, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1400,int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9675, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+0,int_stack+11935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11185, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1760,int_stack+1610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6210, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12835,int_stack+1985,int_stack+1760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6510, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1050,int_stack+12835,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11485, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11935,int_stack+2400,int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9575, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+2550,int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9675, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12835,int_stack+0,int_stack+11935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11185, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2910,int_stack+2760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6210, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1950,int_stack+3135,int_stack+2910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6510, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14410,int_stack+1950,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11485, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11935,int_stack+3550,int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3700,int_stack+3550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1950,int_stack+0,int_stack+11935, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+4060,int_stack+3910, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2550,int_stack+4285,int_stack+4060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3225,int_stack+2550,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11935,int_stack+4700,int_stack+4600, 0.0, zero_stack, 1.0, int_stack+9575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+4850,int_stack+4700, 0.0, zero_stack, 1.0, int_stack+9675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2550,int_stack+0,int_stack+11935, 0.0, zero_stack, 1.0, int_stack+11185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+5210,int_stack+5060, 0.0, zero_stack, 1.0, int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4125,int_stack+5435,int_stack+5210, 0.0, zero_stack, 1.0, int_stack+6510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4800,int_stack+4125,int_stack+0, 0.0, zero_stack, 1.0, int_stack+11485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11935,int_stack+5850,int_stack+5750, 1.0, int_stack+9575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+6000,int_stack+5850, 1.0, int_stack+9675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4125,int_stack+0,int_stack+11935, 1.0, int_stack+11185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+6735,int_stack+6360, 1.0, int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5700,int_stack+6960,int_stack+6735, 1.0, int_stack+6510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6375,int_stack+5700,int_stack+0, 1.0, int_stack+11485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+9825,int_stack+9675,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11485,int_stack+0,int_stack+11185,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11185,int_stack+7375,int_stack+7275,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+7525,int_stack+7375,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5700,int_stack+0,int_stack+11185,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+7885,int_stack+7735,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15310,int_stack+8110,int_stack+7885,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7275,int_stack+15310,int_stack+0,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11185,int_stack+8525,int_stack+8425,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+8675,int_stack+8525,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+15310,int_stack+0,int_stack+11185,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9035,int_stack+8885,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8175,int_stack+9260,int_stack+9035,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8850,int_stack+8175,int_stack+0,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11185,int_stack+10135,int_stack+10035,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+10285,int_stack+10135,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8175,int_stack+0,int_stack+11185,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+10645,int_stack+10495,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9750,int_stack+10870,int_stack+10645,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10425,int_stack+9750,int_stack+0,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15910,int_stack+13510,int_stack+12235,60);
     Libderiv->ABCD[11] = int_stack + 15910;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+17710,int_stack+1050,int_stack+450,60);
     Libderiv->ABCD[10] = int_stack + 17710;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+14410,int_stack+12835,60);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12085,int_stack+3225,int_stack+1950,60);
     Libderiv->ABCD[8] = int_stack + 12085;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+19510,int_stack+4800,int_stack+2550,60);
     Libderiv->ABCD[7] = int_stack + 19510;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1800,int_stack+6375,int_stack+4125,60);
     Libderiv->ABCD[6] = int_stack + 1800;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+3600,int_stack+7275,int_stack+5700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 3600;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+5400,int_stack+8850,int_stack+15310, 0.0, zero_stack, 1.0, int_stack+11485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 5400;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+13885,int_stack+10425,int_stack+8175, 1.0, int_stack+11485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 13885;

}
