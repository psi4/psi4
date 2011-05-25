#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ffhd(Libint_t*, prim_data*);

  /* Computes quartets of (ff|hd) integrals */

REALTYPE *hrr_order_ffhd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[3][7] = int_stack + 490;
 Libint->vrr_classes[4][5] = int_stack + 850;
 Libint->vrr_classes[4][6] = int_stack + 1165;
 Libint->vrr_classes[4][7] = int_stack + 1585;
 Libint->vrr_classes[5][5] = int_stack + 2125;
 Libint->vrr_classes[5][6] = int_stack + 2566;
 Libint->vrr_classes[5][7] = int_stack + 3154;
 Libint->vrr_classes[6][5] = int_stack + 3910;
 Libint->vrr_classes[6][6] = int_stack + 4498;
 Libint->vrr_classes[6][7] = int_stack + 5282;
 memset(int_stack,0,6290*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 6290;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ffhd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6290,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6920,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7760,int_stack+6920,int_stack+6290,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6290,int_stack+1165,int_stack+850,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9020,int_stack+1585,int_stack+1165,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+9020,int_stack+6290,15);
 /*--- compute (fp|hd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+9020,int_stack+0,int_stack+7760,126);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6290,int_stack+2566,int_stack+2125,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12800,int_stack+3154,int_stack+2566,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+14564,int_stack+12800,int_stack+6290,21);
 /*--- compute (gp|hd) ---*/
 hrr1_build_gp(Libint->AB,int_stack+17210,int_stack+14564,int_stack+0,126);
 /*--- compute (fd|hd) ---*/
 hrr1_build_fd(Libint->AB,int_stack+22880,int_stack+17210,int_stack+9020,126);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+4498,int_stack+3910,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+5282,int_stack+4498,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4116,int_stack+1764,int_stack+0,28);
 /*--- compute (hp|hd) ---*/
 hrr1_build_hp(Libint->AB,int_stack+30440,int_stack+4116,int_stack+14564,126);
 /*--- compute (gd|hd) ---*/
 hrr1_build_gd(Libint->AB,int_stack+0,int_stack+30440,int_stack+17210,126);
 /*--- compute (ff|hd) ---*/
 hrr1_build_ff(Libint->AB,int_stack+30440,int_stack+0,int_stack+22880,126);
 return int_stack+30440;}
