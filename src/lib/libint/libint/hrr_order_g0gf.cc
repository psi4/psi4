#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_g0gf(Libint_t*, prim_data*);

  /* Computes quartets of (g0|gf) integrals */

REALTYPE *hrr_order_g0gf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][4] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 225;
 Libint->vrr_classes[4][6] = int_stack + 540;
 Libint->vrr_classes[4][7] = int_stack + 960;
 memset(int_stack,0,1500*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1500;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_g0gf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1500,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2175,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3120,int_stack+2175,int_stack+1500,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4470,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+4470,int_stack+2175,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+4470,int_stack+0,int_stack+3120,15);
 return int_stack+4470;}
