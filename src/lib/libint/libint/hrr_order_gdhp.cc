#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gdhp(Libint_t*, prim_data*);

  /* Computes quartets of (gd|hp) integrals */

REALTYPE *hrr_order_gdhp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[5][5] = int_stack + 735;
 Libint->vrr_classes[5][6] = int_stack + 1176;
 Libint->vrr_classes[6][5] = int_stack + 1764;
 Libint->vrr_classes[6][6] = int_stack + 2352;
 memset(int_stack,0,3136*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3136;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gdhp(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3136,int_stack+315,int_stack+0,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4081,int_stack+1176,int_stack+735,21);
 /*--- compute (gp|hp) ---*/
 hrr1_build_gp(Libint->AB,int_stack+5404,int_stack+4081,int_stack+3136,63);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+2352,int_stack+1764,28);
 /*--- compute (hp|hp) ---*/
 hrr1_build_hp(Libint->AB,int_stack+8239,int_stack+0,int_stack+4081,63);
 /*--- compute (gd|hp) ---*/
 hrr1_build_gd(Libint->AB,int_stack+12208,int_stack+8239,int_stack+5404,63);
 return int_stack+12208;}
