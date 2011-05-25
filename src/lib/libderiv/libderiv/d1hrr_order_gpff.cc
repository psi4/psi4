#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gpff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gp|ff) integrals */

void d1hrr_order_gpff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][3][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][11] = int_stack + 150;
 Libderiv->deriv_classes[4][5][11] = int_stack + 375;
 Libderiv->deriv_classes[4][6][11] = int_stack + 690;
 Libderiv->deriv_classes[5][3][11] = int_stack + 1110;
 Libderiv->deriv_classes[5][4][11] = int_stack + 1320;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1635;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2076;
 Libderiv->deriv_classes[4][3][10] = int_stack + 2664;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2814;
 Libderiv->deriv_classes[4][5][10] = int_stack + 3039;
 Libderiv->deriv_classes[4][6][10] = int_stack + 3354;
 Libderiv->deriv_classes[5][3][10] = int_stack + 3774;
 Libderiv->deriv_classes[5][4][10] = int_stack + 3984;
 Libderiv->deriv_classes[5][5][10] = int_stack + 4299;
 Libderiv->deriv_classes[5][6][10] = int_stack + 4740;
 Libderiv->deriv_classes[4][3][9] = int_stack + 5328;
 Libderiv->deriv_classes[4][4][9] = int_stack + 5478;
 Libderiv->deriv_classes[4][5][9] = int_stack + 5703;
 Libderiv->deriv_classes[4][6][9] = int_stack + 6018;
 Libderiv->deriv_classes[5][3][9] = int_stack + 6438;
 Libderiv->deriv_classes[5][4][9] = int_stack + 6648;
 Libderiv->deriv_classes[5][5][9] = int_stack + 6963;
 Libderiv->deriv_classes[5][6][9] = int_stack + 7404;
 Libderiv->deriv_classes[4][3][8] = int_stack + 7992;
 Libderiv->deriv_classes[4][4][8] = int_stack + 8142;
 Libderiv->deriv_classes[4][5][8] = int_stack + 8367;
 Libderiv->deriv_classes[4][6][8] = int_stack + 8682;
 Libderiv->deriv_classes[5][3][8] = int_stack + 9102;
 Libderiv->deriv_classes[5][4][8] = int_stack + 9312;
 Libderiv->deriv_classes[5][5][8] = int_stack + 9627;
 Libderiv->deriv_classes[5][6][8] = int_stack + 10068;
 Libderiv->deriv_classes[4][3][7] = int_stack + 10656;
 Libderiv->deriv_classes[4][4][7] = int_stack + 10806;
 Libderiv->deriv_classes[4][5][7] = int_stack + 11031;
 Libderiv->deriv_classes[4][6][7] = int_stack + 11346;
 Libderiv->deriv_classes[5][3][7] = int_stack + 11766;
 Libderiv->deriv_classes[5][4][7] = int_stack + 11976;
 Libderiv->deriv_classes[5][5][7] = int_stack + 12291;
 Libderiv->deriv_classes[5][6][7] = int_stack + 12732;
 Libderiv->deriv_classes[4][3][6] = int_stack + 13320;
 Libderiv->deriv_classes[4][4][6] = int_stack + 13470;
 Libderiv->deriv_classes[4][5][6] = int_stack + 13695;
 Libderiv->deriv_classes[4][6][6] = int_stack + 14010;
 Libderiv->dvrr_classes[5][3] = int_stack + 14430;
 Libderiv->deriv_classes[5][3][6] = int_stack + 14640;
 Libderiv->dvrr_classes[5][4] = int_stack + 14850;
 Libderiv->deriv_classes[5][4][6] = int_stack + 15165;
 Libderiv->dvrr_classes[5][5] = int_stack + 15480;
 Libderiv->deriv_classes[5][5][6] = int_stack + 15921;
 Libderiv->deriv_classes[5][6][6] = int_stack + 16362;
 Libderiv->deriv_classes[4][3][2] = int_stack + 16950;
 Libderiv->deriv_classes[4][4][2] = int_stack + 17100;
 Libderiv->deriv_classes[4][5][2] = int_stack + 17325;
 Libderiv->deriv_classes[4][6][2] = int_stack + 17640;
 Libderiv->deriv_classes[5][3][2] = int_stack + 18060;
 Libderiv->deriv_classes[5][4][2] = int_stack + 18270;
 Libderiv->deriv_classes[5][5][2] = int_stack + 18585;
 Libderiv->deriv_classes[5][6][2] = int_stack + 19026;
 Libderiv->deriv_classes[4][3][1] = int_stack + 19614;
 Libderiv->deriv_classes[4][4][1] = int_stack + 19764;
 Libderiv->deriv_classes[4][5][1] = int_stack + 19989;
 Libderiv->deriv_classes[4][6][1] = int_stack + 20304;
 Libderiv->deriv_classes[5][3][1] = int_stack + 20724;
 Libderiv->deriv_classes[5][4][1] = int_stack + 20934;
 Libderiv->deriv_classes[5][5][1] = int_stack + 21249;
 Libderiv->deriv_classes[5][6][1] = int_stack + 21690;
 Libderiv->dvrr_classes[4][3] = int_stack + 22278;
 Libderiv->dvrr_classes[4][4] = int_stack + 22428;
 Libderiv->dvrr_classes[4][5] = int_stack + 22653;
 Libderiv->dvrr_classes[4][6] = int_stack + 22968;
 Libderiv->deriv_classes[4][3][0] = int_stack + 23388;
 Libderiv->deriv_classes[4][4][0] = int_stack + 23538;
 Libderiv->deriv_classes[4][5][0] = int_stack + 23763;
 Libderiv->deriv_classes[4][6][0] = int_stack + 24078;
 Libderiv->deriv_classes[5][3][0] = int_stack + 24498;
 Libderiv->deriv_classes[5][4][0] = int_stack + 24708;
 Libderiv->deriv_classes[5][5][0] = int_stack + 25023;
 Libderiv->deriv_classes[5][6][0] = int_stack + 25464;
 memset(int_stack,0,208416);

 Libderiv->dvrr_stack = int_stack + 56622;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gpff(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+26052,int_stack+22428,int_stack+22278,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+26502,int_stack+22653,int_stack+22428,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+27177,int_stack+26502,int_stack+26052,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28077,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22278,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+28527,int_stack+375,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22428,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29202,int_stack+28527,int_stack+28077, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26052,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30102,int_stack+690,int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22653,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31047,int_stack+30102,int_stack+28527, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26502,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+32397,int_stack+31047,int_stack+29202, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27177,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+28077,int_stack+14850,int_stack+14430,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28707,int_stack+15480,int_stack+14850,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+29652,int_stack+28707,int_stack+28077,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30912,int_stack+1320,int_stack+1110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14430,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1635,int_stack+1320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14850,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33897,int_stack+0,int_stack+30912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28077,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30912,int_stack+2076,int_stack+1635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15480,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+35157,int_stack+30912,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28707,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+35157,int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29652,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33897,int_stack+2814,int_stack+2664, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22278, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34347,int_stack+3039,int_stack+2814, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22428, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35022,int_stack+34347,int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26052, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35922,int_stack+3354,int_stack+3039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22653, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2100,int_stack+35922,int_stack+34347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26502, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+35922,int_stack+2100,int_stack+35022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27177, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+3984,int_stack+3774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14430, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2730,int_stack+4299,int_stack+3984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14850, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33897,int_stack+2730,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28077, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30912,int_stack+4740,int_stack+4299, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15480, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+37422,int_stack+30912,int_stack+2730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28707, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2100,int_stack+37422,int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29652, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33897,int_stack+5478,int_stack+5328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22278, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34347,int_stack+5703,int_stack+5478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22428, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35022,int_stack+34347,int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26052, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+37422,int_stack+6018,int_stack+5703, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22653, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+38367,int_stack+37422,int_stack+34347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26502, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4200,int_stack+38367,int_stack+35022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27177, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37422,int_stack+6648,int_stack+6438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14430, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+38052,int_stack+6963,int_stack+6648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14850, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33897,int_stack+38052,int_stack+37422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28077, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30912,int_stack+7404,int_stack+6963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15480, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5700,int_stack+30912,int_stack+38052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28707, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+37422,int_stack+5700,int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29652, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33897,int_stack+8142,int_stack+7992, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34347,int_stack+8367,int_stack+8142, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35022,int_stack+34347,int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5700,int_stack+8682,int_stack+8367, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22653, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6645,int_stack+5700,int_stack+34347, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+39522,int_stack+6645,int_stack+35022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5700,int_stack+9312,int_stack+9102, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6330,int_stack+9627,int_stack+9312, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7275,int_stack+6330,int_stack+5700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28077, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+33897,int_stack+10068,int_stack+9627, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8535,int_stack+33897,int_stack+6330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28707, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+41022,int_stack+8535,int_stack+7275, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29652, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33897,int_stack+10806,int_stack+10656, 0.0, zero_stack, 1.0, int_stack+22278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34347,int_stack+11031,int_stack+10806, 0.0, zero_stack, 1.0, int_stack+22428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35022,int_stack+34347,int_stack+33897, 0.0, zero_stack, 1.0, int_stack+26052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5700,int_stack+11346,int_stack+11031, 0.0, zero_stack, 1.0, int_stack+22653, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6645,int_stack+5700,int_stack+34347, 0.0, zero_stack, 1.0, int_stack+26502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7995,int_stack+6645,int_stack+35022, 0.0, zero_stack, 1.0, int_stack+27177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5700,int_stack+11976,int_stack+11766, 0.0, zero_stack, 1.0, int_stack+14430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6330,int_stack+12291,int_stack+11976, 0.0, zero_stack, 1.0, int_stack+14850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33897,int_stack+6330,int_stack+5700, 0.0, zero_stack, 1.0, int_stack+28077, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9495,int_stack+12732,int_stack+12291, 0.0, zero_stack, 1.0, int_stack+15480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10818,int_stack+9495,int_stack+6330, 0.0, zero_stack, 1.0, int_stack+28707, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5700,int_stack+10818,int_stack+33897, 0.0, zero_stack, 1.0, int_stack+29652, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33897,int_stack+13470,int_stack+13320, 1.0, int_stack+22278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34347,int_stack+13695,int_stack+13470, 1.0, int_stack+22428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35022,int_stack+34347,int_stack+33897, 1.0, int_stack+26052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9495,int_stack+14010,int_stack+13695, 1.0, int_stack+22653, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10440,int_stack+9495,int_stack+34347, 1.0, int_stack+26502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11790,int_stack+10440,int_stack+35022, 1.0, int_stack+27177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9495,int_stack+15165,int_stack+14640, 1.0, int_stack+14430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10125,int_stack+15921,int_stack+15165, 1.0, int_stack+14850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33897,int_stack+10125,int_stack+9495, 1.0, int_stack+28077, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13290,int_stack+16362,int_stack+15921, 1.0, int_stack+15480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14613,int_stack+13290,int_stack+10125, 1.0, int_stack+28707, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9495,int_stack+14613,int_stack+33897, 1.0, int_stack+29652, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+33897,int_stack+22968,int_stack+22653,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+13290,int_stack+33897,int_stack+26502,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+33897,int_stack+13290,int_stack+27177,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13290,int_stack+17100,int_stack+16950,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13740,int_stack+17325,int_stack+17100,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14415,int_stack+13740,int_stack+13290,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15315,int_stack+17640,int_stack+17325,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16260,int_stack+15315,int_stack+13740,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+26052,int_stack+16260,int_stack+14415,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13290,int_stack+18270,int_stack+18060,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13920,int_stack+18585,int_stack+18270,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14865,int_stack+13920,int_stack+13290,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16125,int_stack+19026,int_stack+18585,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17448,int_stack+16125,int_stack+13920,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+27552,int_stack+17448,int_stack+14865,21);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13290,int_stack+19764,int_stack+19614,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13740,int_stack+19989,int_stack+19764,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14415,int_stack+13740,int_stack+13290,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15315,int_stack+20304,int_stack+19989,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16260,int_stack+15315,int_stack+13740,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+17610,int_stack+16260,int_stack+14415,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13290,int_stack+20934,int_stack+20724,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13920,int_stack+21249,int_stack+20934,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14865,int_stack+13920,int_stack+13290,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16125,int_stack+21690,int_stack+21249,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+19110,int_stack+16125,int_stack+13920,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+21000,int_stack+19110,int_stack+14865,21);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19110,int_stack+23538,int_stack+23388,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19560,int_stack+23763,int_stack+23538,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+13290,int_stack+19560,int_stack+19110,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+14190,int_stack+24078,int_stack+23763,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15135,int_stack+14190,int_stack+19560,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+19110,int_stack+15135,int_stack+13290,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13290,int_stack+24708,int_stack+24498,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13920,int_stack+25023,int_stack+24708,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14865,int_stack+13920,int_stack+13290,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16125,int_stack+25464,int_stack+25023,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23100,int_stack+16125,int_stack+13920,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+29652,int_stack+23100,int_stack+14865,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+43122,int_stack+0,int_stack+32397,100);
     Libderiv->ABCD[11] = int_stack + 43122;
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+47622,int_stack+2100,int_stack+35922,100);
     Libderiv->ABCD[10] = int_stack + 47622;
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+52122,int_stack+37422,int_stack+4200,100);
     Libderiv->ABCD[9] = int_stack + 52122;
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+41022,int_stack+39522,100);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+35397,int_stack+5700,int_stack+7995,100);
     Libderiv->ABCD[7] = int_stack + 35397;
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+4500,int_stack+9495,int_stack+11790,100);
     Libderiv->ABCD[6] = int_stack + 4500;
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+9000,int_stack+27552,int_stack+26052, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 9000;
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+23100,int_stack+21000,int_stack+17610, 0.0, zero_stack, 1.0, int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 23100;
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+13500,int_stack+29652,int_stack+19110, 1.0, int_stack+33897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 13500;

}
