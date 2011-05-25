#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gdhd(Libint_t*, prim_data*);

  /* Computes quartets of (gd|hd) integrals */

REALTYPE *hrr_order_gdhd(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,5440*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 5440;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gdhd(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5440,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6385,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7645,int_stack+6385,int_stack+5440,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5440,int_stack+1716,int_stack+1275,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9535,int_stack+2304,int_stack+1716,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+9535,int_stack+5440,21);
 /*--- compute (gp|hd) ---*/
 hrr1_build_gp(Libint->AB,int_stack+9535,int_stack+0,int_stack+7645,126);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5440,int_stack+3648,int_stack+3060,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+15205,int_stack+4432,int_stack+3648,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+17557,int_stack+15205,int_stack+5440,28);
 /*--- compute (hp|hd) ---*/
 hrr1_build_hp(Libint->AB,int_stack+21085,int_stack+17557,int_stack+0,126);
 /*--- compute (gd|hd) ---*/
 hrr1_build_gd(Libint->AB,int_stack+29023,int_stack+21085,int_stack+9535,126);
 return int_stack+29023;}
