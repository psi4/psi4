#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|ff) integrals */

void d1hrr_order_00ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][3][11] = int_stack + 0;
 Libderiv->deriv_classes[0][4][11] = int_stack + 10;
 Libderiv->deriv_classes[0][5][11] = int_stack + 25;
 Libderiv->deriv_classes[0][6][11] = int_stack + 46;
 Libderiv->deriv_classes[0][3][10] = int_stack + 74;
 Libderiv->deriv_classes[0][4][10] = int_stack + 84;
 Libderiv->deriv_classes[0][5][10] = int_stack + 99;
 Libderiv->deriv_classes[0][6][10] = int_stack + 120;
 Libderiv->deriv_classes[0][3][9] = int_stack + 148;
 Libderiv->deriv_classes[0][4][9] = int_stack + 158;
 Libderiv->deriv_classes[0][5][9] = int_stack + 173;
 Libderiv->deriv_classes[0][6][9] = int_stack + 194;
 Libderiv->deriv_classes[0][3][8] = int_stack + 222;
 Libderiv->deriv_classes[0][4][8] = int_stack + 232;
 Libderiv->deriv_classes[0][5][8] = int_stack + 247;
 Libderiv->deriv_classes[0][6][8] = int_stack + 268;
 Libderiv->deriv_classes[0][3][7] = int_stack + 296;
 Libderiv->deriv_classes[0][4][7] = int_stack + 306;
 Libderiv->deriv_classes[0][5][7] = int_stack + 321;
 Libderiv->deriv_classes[0][6][7] = int_stack + 342;
 Libderiv->dvrr_classes[0][3] = int_stack + 370;
 Libderiv->deriv_classes[0][3][6] = int_stack + 380;
 Libderiv->dvrr_classes[0][4] = int_stack + 390;
 Libderiv->deriv_classes[0][4][6] = int_stack + 405;
 Libderiv->dvrr_classes[0][5] = int_stack + 420;
 Libderiv->deriv_classes[0][5][6] = int_stack + 441;
 Libderiv->deriv_classes[0][6][6] = int_stack + 462;
 Libderiv->deriv_classes[0][3][2] = int_stack + 490;
 Libderiv->deriv_classes[0][4][2] = int_stack + 500;
 Libderiv->deriv_classes[0][5][2] = int_stack + 515;
 Libderiv->deriv_classes[0][6][2] = int_stack + 536;
 Libderiv->deriv_classes[0][3][1] = int_stack + 564;
 Libderiv->deriv_classes[0][4][1] = int_stack + 574;
 Libderiv->deriv_classes[0][5][1] = int_stack + 589;
 Libderiv->deriv_classes[0][6][1] = int_stack + 610;
 Libderiv->deriv_classes[0][3][0] = int_stack + 638;
 Libderiv->deriv_classes[0][4][0] = int_stack + 648;
 Libderiv->deriv_classes[0][5][0] = int_stack + 663;
 Libderiv->deriv_classes[0][6][0] = int_stack + 684;
 memset(int_stack,0,5696);

 Libderiv->dvrr_stack = int_stack + 1923;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+712,int_stack+390,int_stack+370,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+742,int_stack+420,int_stack+390,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+787,int_stack+742,int_stack+712,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+847,int_stack+10,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+370,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+877,int_stack+25,int_stack+10, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+390,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+922,int_stack+877,int_stack+847, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+712,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+982,int_stack+46,int_stack+25, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1045,int_stack+982,int_stack+877, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+742,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+982,int_stack+84,int_stack+74, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+370, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+847,int_stack+99,int_stack+84, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+390, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+847,int_stack+982, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+712, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+982,int_stack+120,int_stack+99, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1135,int_stack+982,int_stack+847, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+742, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+847,int_stack+158,int_stack+148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+370, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+877,int_stack+173,int_stack+158, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+390, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+982,int_stack+877,int_stack+847, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+712, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60,int_stack+194,int_stack+173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+123,int_stack+60,int_stack+877, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+742, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60,int_stack+232,int_stack+222, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+847,int_stack+247,int_stack+232, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1225,int_stack+847,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60,int_stack+268,int_stack+247, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1285,int_stack+60,int_stack+847, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+847,int_stack+306,int_stack+296, 0.0, zero_stack, 1.0, int_stack+370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+877,int_stack+321,int_stack+306, 0.0, zero_stack, 1.0, int_stack+390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+60,int_stack+877,int_stack+847, 0.0, zero_stack, 1.0, int_stack+712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+213,int_stack+342,int_stack+321, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+276,int_stack+213,int_stack+877, 0.0, zero_stack, 1.0, int_stack+742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+213,int_stack+405,int_stack+380, 1.0, int_stack+370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+847,int_stack+441,int_stack+405, 1.0, int_stack+390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1375,int_stack+847,int_stack+213, 1.0, int_stack+712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+213,int_stack+462,int_stack+441, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+366,int_stack+213,int_stack+847, 1.0, int_stack+742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+847,int_stack+500,int_stack+490,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+877,int_stack+515,int_stack+500,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+213,int_stack+877,int_stack+847,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+712,int_stack+536,int_stack+515,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+456,int_stack+712,int_stack+877,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+712,int_stack+574,int_stack+564,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+742,int_stack+589,int_stack+574,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+847,int_stack+742,int_stack+712,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1435,int_stack+610,int_stack+589,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+546,int_stack+1435,int_stack+742,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1435,int_stack+648,int_stack+638,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1465,int_stack+663,int_stack+648,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1510,int_stack+1465,int_stack+1435,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1570,int_stack+684,int_stack+663,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1633,int_stack+1570,int_stack+1465,1);
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1723,int_stack+1045,int_stack+922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+787,1);
     Libderiv->ABCD[11] = int_stack + 1723;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1823,int_stack+1135,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+787, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 1823;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1042,int_stack+123,int_stack+982, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+787, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 1042;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+907,int_stack+1285,int_stack+1225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 907;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1142,int_stack+276,int_stack+60, 0.0, zero_stack, 1.0, int_stack+787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 1142;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+366,int_stack+1375, 1.0, int_stack+787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+100,int_stack+456,int_stack+213,1);
     Libderiv->ABCD[2] = int_stack + 100;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+200,int_stack+546,int_stack+847,1);
     Libderiv->ABCD[1] = int_stack + 200;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+300,int_stack+1633,int_stack+1510,1);
     Libderiv->ABCD[0] = int_stack + 300;

}
