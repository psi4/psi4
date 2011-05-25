#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_g0hd(Libint_t*, prim_data*);

  /* Computes quartets of (g0|hd) integrals */

REALTYPE *hrr_order_g0hd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 memset(int_stack,0,1275*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1275;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_g0hd(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1275,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2220,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3480,int_stack+2220,int_stack+1275,15);
 return int_stack+3480;}
