#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00gf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|gf) integrals */

void d1hrr_order_00gf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][4][11] = int_stack + 0;
 Libderiv->deriv_classes[0][5][11] = int_stack + 15;
 Libderiv->deriv_classes[0][6][11] = int_stack + 36;
 Libderiv->deriv_classes[0][7][11] = int_stack + 64;
 Libderiv->deriv_classes[0][4][10] = int_stack + 100;
 Libderiv->deriv_classes[0][5][10] = int_stack + 115;
 Libderiv->deriv_classes[0][6][10] = int_stack + 136;
 Libderiv->deriv_classes[0][7][10] = int_stack + 164;
 Libderiv->deriv_classes[0][4][9] = int_stack + 200;
 Libderiv->deriv_classes[0][5][9] = int_stack + 215;
 Libderiv->deriv_classes[0][6][9] = int_stack + 236;
 Libderiv->deriv_classes[0][7][9] = int_stack + 264;
 Libderiv->deriv_classes[0][4][8] = int_stack + 300;
 Libderiv->deriv_classes[0][5][8] = int_stack + 315;
 Libderiv->deriv_classes[0][6][8] = int_stack + 336;
 Libderiv->deriv_classes[0][7][8] = int_stack + 364;
 Libderiv->deriv_classes[0][4][7] = int_stack + 400;
 Libderiv->deriv_classes[0][5][7] = int_stack + 415;
 Libderiv->deriv_classes[0][6][7] = int_stack + 436;
 Libderiv->deriv_classes[0][7][7] = int_stack + 464;
 Libderiv->dvrr_classes[0][4] = int_stack + 500;
 Libderiv->deriv_classes[0][4][6] = int_stack + 515;
 Libderiv->dvrr_classes[0][5] = int_stack + 530;
 Libderiv->deriv_classes[0][5][6] = int_stack + 551;
 Libderiv->dvrr_classes[0][6] = int_stack + 572;
 Libderiv->deriv_classes[0][6][6] = int_stack + 600;
 Libderiv->deriv_classes[0][7][6] = int_stack + 628;
 Libderiv->deriv_classes[0][4][2] = int_stack + 664;
 Libderiv->deriv_classes[0][5][2] = int_stack + 679;
 Libderiv->deriv_classes[0][6][2] = int_stack + 700;
 Libderiv->deriv_classes[0][7][2] = int_stack + 728;
 Libderiv->deriv_classes[0][4][1] = int_stack + 764;
 Libderiv->deriv_classes[0][5][1] = int_stack + 779;
 Libderiv->deriv_classes[0][6][1] = int_stack + 800;
 Libderiv->deriv_classes[0][7][1] = int_stack + 828;
 Libderiv->deriv_classes[0][4][0] = int_stack + 864;
 Libderiv->deriv_classes[0][5][0] = int_stack + 879;
 Libderiv->deriv_classes[0][6][0] = int_stack + 900;
 Libderiv->deriv_classes[0][7][0] = int_stack + 928;
 memset(int_stack,0,7712);

 Libderiv->dvrr_stack = int_stack + 2494;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00gf(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+964,int_stack+530,int_stack+500,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1009,int_stack+572,int_stack+530,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1072,int_stack+1009,int_stack+964,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1162,int_stack+15,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1207,int_stack+36,int_stack+15, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+530,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1270,int_stack+1207,int_stack+1162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+964,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+64,int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+572,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1444,int_stack+1360,int_stack+1207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1009,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1360,int_stack+115,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1162,int_stack+136,int_stack+115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+530, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+1162,int_stack+1360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+964, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+164,int_stack+136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+572, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1570,int_stack+1360,int_stack+1162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1009, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1162,int_stack+215,int_stack+200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1207,int_stack+236,int_stack+215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+530, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+90,int_stack+1207,int_stack+1162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+964, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+264,int_stack+236, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+572, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1696,int_stack+1360,int_stack+1207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1009, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1360,int_stack+315,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1162,int_stack+336,int_stack+315, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+180,int_stack+1162,int_stack+1360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+364,int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+572, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+270,int_stack+1360,int_stack+1162, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1009, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1162,int_stack+415,int_stack+400, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1207,int_stack+436,int_stack+415, 0.0, zero_stack, 1.0, int_stack+530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1822,int_stack+1207,int_stack+1162, 0.0, zero_stack, 1.0, int_stack+964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+464,int_stack+436, 0.0, zero_stack, 1.0, int_stack+572, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1912,int_stack+1360,int_stack+1207, 0.0, zero_stack, 1.0, int_stack+1009, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1360,int_stack+551,int_stack+515, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1162,int_stack+600,int_stack+551, 1.0, int_stack+530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+396,int_stack+1162,int_stack+1360, 1.0, int_stack+964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+628,int_stack+600, 1.0, int_stack+572, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+486,int_stack+1360,int_stack+1162, 1.0, int_stack+1009, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1162,int_stack+679,int_stack+664,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1207,int_stack+700,int_stack+679,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+964,int_stack+1207,int_stack+1162,1);
 /*--- compute (00|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+728,int_stack+700,1);
 /*--- compute (00|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+612,int_stack+1360,int_stack+1207,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1360,int_stack+779,int_stack+764,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1162,int_stack+800,int_stack+779,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2038,int_stack+1162,int_stack+1360,1);
 /*--- compute (00|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+828,int_stack+800,1);
 /*--- compute (00|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+738,int_stack+1360,int_stack+1162,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1162,int_stack+879,int_stack+864,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1207,int_stack+900,int_stack+879,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2128,int_stack+1207,int_stack+1162,1);
 /*--- compute (00|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+1360,int_stack+928,int_stack+900,1);
 /*--- compute (00|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+2218,int_stack+1360,int_stack+1207,1);
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2344,int_stack+1444,int_stack+1270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1072,1);
     Libderiv->ABCD[11] = int_stack + 2344;
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1162,int_stack+1570,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1072, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 1162;
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1312,int_stack+1696,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1072, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 1312;
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+270,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+150,int_stack+1912,int_stack+1822, 0.0, zero_stack, 1.0, int_stack+1072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 150;
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1462,int_stack+486,int_stack+396, 1.0, int_stack+1072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 1462;
 /*--- compute (00|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+300,int_stack+612,int_stack+964,1);
     Libderiv->ABCD[2] = int_stack + 300;
 /*--- compute (00|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+450,int_stack+738,int_stack+2038,1);
     Libderiv->ABCD[1] = int_stack + 450;
 /*--- compute (00|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+600,int_stack+2218,int_stack+2128,1);
     Libderiv->ABCD[0] = int_stack + 600;

}
