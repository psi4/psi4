#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_pphd(Libint_t*, prim_data*);

  /* Computes quartets of (pp|hd) integrals */

REALTYPE *hrr_order_pphd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 Libint->vrr_classes[1][7] = int_stack + 147;
 Libint->vrr_classes[2][5] = int_stack + 255;
 Libint->vrr_classes[2][6] = int_stack + 381;
 Libint->vrr_classes[2][7] = int_stack + 549;
 memset(int_stack,0,765*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 765;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_pphd(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+765,int_stack+63,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+954,int_stack+147,int_stack+63,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+1206,int_stack+954,int_stack+765,3);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+765,int_stack+381,int_stack+255,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1584,int_stack+549,int_stack+381,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+1584,int_stack+765,6);
 /*--- compute (pp|hd) ---*/
 hrr1_build_pp(Libint->AB,int_stack+1584,int_stack+0,int_stack+1206,126);
 return int_stack+1584;}
