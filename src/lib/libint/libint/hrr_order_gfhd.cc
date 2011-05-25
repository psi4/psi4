#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gfhd(Libint_t*, prim_data*);

  /* Computes quartets of (gf|hd) integrals */

REALTYPE *hrr_order_gfhd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[5][5] = int_stack + 1275;
 Libint->vrr_classes[5][6] = int_stack + 1716;
 Libint->vrr_classes[5][7] = int_stack + 2304;
 Libint->vrr_classes[6][5] = int_stack + 3060;
 Libint->vrr_classes[6][6] = int_stack + 3648;
 Libint->vrr_classes[6][7] = int_stack + 4432;
 Libint->vrr_classes[7][5] = int_stack + 5440;
 Libint->vrr_classes[7][6] = int_stack + 6196;
 Libint->vrr_classes[7][7] = int_stack + 7204;
 memset(int_stack,0,8500*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 8500;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gfhd(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8500,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9445,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10705,int_stack+9445,int_stack+8500,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8500,int_stack+1716,int_stack+1275,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12595,int_stack+2304,int_stack+1716,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+12595,int_stack+8500,21);
 /*--- compute (gp|hd) ---*/
 hrr1_build_gp(Libint->AB,int_stack+12595,int_stack+0,int_stack+10705,126);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8500,int_stack+3648,int_stack+3060,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+18265,int_stack+4432,int_stack+3648,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+20617,int_stack+18265,int_stack+8500,28);
 /*--- compute (hp|hd) ---*/
 hrr1_build_hp(Libint->AB,int_stack+24145,int_stack+20617,int_stack+0,126);
 /*--- compute (gd|hd) ---*/
 hrr1_build_gd(Libint->AB,int_stack+32083,int_stack+24145,int_stack+12595,126);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+6196,int_stack+5440,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+7204,int_stack+6196,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (ip|hd) ---*/
 hrr1_build_ip(Libint->AB,int_stack+9828,int_stack+5292,int_stack+20617,126);
 /*--- compute (hd|hd) ---*/
 hrr1_build_hd(Libint->AB,int_stack+43423,int_stack+9828,int_stack+24145,126);
 /*--- compute (gf|hd) ---*/
 hrr1_build_gf(Libint->AB,int_stack+0,int_stack+43423,int_stack+32083,126);
 return int_stack+0;}
