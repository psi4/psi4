#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gpgp(Libint_t*, prim_data*);

  /* Computes quartets of (gp|gp) integrals */

REALTYPE *hrr_order_gpgp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][4] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 225;
 Libint->vrr_classes[5][4] = int_stack + 540;
 Libint->vrr_classes[5][5] = int_stack + 855;
 memset(int_stack,0,1296*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1296;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gpgp(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1296,int_stack+225,int_stack+0,15);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1971,int_stack+855,int_stack+540,21);
 /*--- compute (gp|gp) ---*/
 hrr1_build_gp(Libint->AB,int_stack+2916,int_stack+1971,int_stack+1296,45);
 return int_stack+2916;}
