#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0gf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|gf) integrals */

void d1hrr_order_p0gf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][5][11] = int_stack + 45;
 Libderiv->deriv_classes[1][6][11] = int_stack + 108;
 Libderiv->deriv_classes[1][7][11] = int_stack + 192;
 Libderiv->deriv_classes[1][4][10] = int_stack + 300;
 Libderiv->deriv_classes[1][5][10] = int_stack + 345;
 Libderiv->deriv_classes[1][6][10] = int_stack + 408;
 Libderiv->deriv_classes[1][7][10] = int_stack + 492;
 Libderiv->deriv_classes[1][4][9] = int_stack + 600;
 Libderiv->deriv_classes[1][5][9] = int_stack + 645;
 Libderiv->deriv_classes[1][6][9] = int_stack + 708;
 Libderiv->deriv_classes[1][7][9] = int_stack + 792;
 Libderiv->deriv_classes[1][4][8] = int_stack + 900;
 Libderiv->deriv_classes[1][5][8] = int_stack + 945;
 Libderiv->deriv_classes[1][6][8] = int_stack + 1008;
 Libderiv->deriv_classes[1][7][8] = int_stack + 1092;
 Libderiv->deriv_classes[1][4][7] = int_stack + 1200;
 Libderiv->deriv_classes[1][5][7] = int_stack + 1245;
 Libderiv->deriv_classes[1][6][7] = int_stack + 1308;
 Libderiv->deriv_classes[1][7][7] = int_stack + 1392;
 Libderiv->dvrr_classes[1][4] = int_stack + 1500;
 Libderiv->deriv_classes[1][4][6] = int_stack + 1545;
 Libderiv->dvrr_classes[1][5] = int_stack + 1590;
 Libderiv->deriv_classes[1][5][6] = int_stack + 1653;
 Libderiv->dvrr_classes[1][6] = int_stack + 1716;
 Libderiv->deriv_classes[1][6][6] = int_stack + 1800;
 Libderiv->deriv_classes[1][7][6] = int_stack + 1884;
 Libderiv->deriv_classes[1][4][2] = int_stack + 1992;
 Libderiv->deriv_classes[1][5][2] = int_stack + 2037;
 Libderiv->deriv_classes[1][6][2] = int_stack + 2100;
 Libderiv->deriv_classes[1][7][2] = int_stack + 2184;
 Libderiv->deriv_classes[1][4][1] = int_stack + 2292;
 Libderiv->deriv_classes[1][5][1] = int_stack + 2337;
 Libderiv->deriv_classes[1][6][1] = int_stack + 2400;
 Libderiv->deriv_classes[1][7][1] = int_stack + 2484;
 Libderiv->deriv_classes[1][4][0] = int_stack + 2592;
 Libderiv->deriv_classes[1][5][0] = int_stack + 2637;
 Libderiv->deriv_classes[1][6][0] = int_stack + 2700;
 Libderiv->deriv_classes[1][7][0] = int_stack + 2784;
 memset(int_stack,0,23136);

 Libderiv->dvrr_stack = int_stack + 7482;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0gf(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2892,int_stack+1590,int_stack+1500,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3027,int_stack+1716,int_stack+1590,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3216,int_stack+3027,int_stack+2892,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3486,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3621,int_stack+108,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1590,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3810,int_stack+3621,int_stack+3486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2892,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+192,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1716,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4332,int_stack+4080,int_stack+3621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3027,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4080,int_stack+345,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3486,int_stack+408,int_stack+345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1590, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+3486,int_stack+4080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2892, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+492,int_stack+408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1716, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4710,int_stack+4080,int_stack+3486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3027, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3486,int_stack+645,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3621,int_stack+708,int_stack+645, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1590, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+270,int_stack+3621,int_stack+3486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2892, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+792,int_stack+708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1716, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+5088,int_stack+4080,int_stack+3621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3027, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4080,int_stack+945,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3486,int_stack+1008,int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+540,int_stack+3486,int_stack+4080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+1092,int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+810,int_stack+4080,int_stack+3486, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3027, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3486,int_stack+1245,int_stack+1200, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3621,int_stack+1308,int_stack+1245, 0.0, zero_stack, 1.0, int_stack+1590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5466,int_stack+3621,int_stack+3486, 0.0, zero_stack, 1.0, int_stack+2892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+1392,int_stack+1308, 0.0, zero_stack, 1.0, int_stack+1716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+5736,int_stack+4080,int_stack+3621, 0.0, zero_stack, 1.0, int_stack+3027, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4080,int_stack+1653,int_stack+1545, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3486,int_stack+1800,int_stack+1653, 1.0, int_stack+1590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1188,int_stack+3486,int_stack+4080, 1.0, int_stack+2892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+1884,int_stack+1800, 1.0, int_stack+1716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1458,int_stack+4080,int_stack+3486, 1.0, int_stack+3027, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3486,int_stack+2037,int_stack+1992,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3621,int_stack+2100,int_stack+2037,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2892,int_stack+3621,int_stack+3486,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+2184,int_stack+2100,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+1836,int_stack+4080,int_stack+3621,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4080,int_stack+2337,int_stack+2292,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3486,int_stack+2400,int_stack+2337,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6114,int_stack+3486,int_stack+4080,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+2484,int_stack+2400,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+2214,int_stack+4080,int_stack+3486,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3486,int_stack+2637,int_stack+2592,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3621,int_stack+2700,int_stack+2637,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6384,int_stack+3621,int_stack+3486,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+4080,int_stack+2784,int_stack+2700,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+6654,int_stack+4080,int_stack+3621,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+7032,int_stack+4332,int_stack+3810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3216,3);
     Libderiv->ABCD[11] = int_stack + 7032;
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3486,int_stack+4710,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3216, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 3486;
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3936,int_stack+5088,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3216, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 3936;
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+810,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+450,int_stack+5736,int_stack+5466, 0.0, zero_stack, 1.0, int_stack+3216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 450;
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4386,int_stack+1458,int_stack+1188, 1.0, int_stack+3216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 4386;
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+900,int_stack+1836,int_stack+2892,3);
     Libderiv->ABCD[2] = int_stack + 900;
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+1350,int_stack+2214,int_stack+6114,3);
     Libderiv->ABCD[1] = int_stack + 1350;
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+1800,int_stack+6654,int_stack+6384,3);
     Libderiv->ABCD[0] = int_stack + 1800;

}
