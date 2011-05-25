#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fpgp(Libint_t*, prim_data*);

  /* Computes quartets of (fp|gp) integrals */

REALTYPE *hrr_order_fpgp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[4][4] = int_stack + 360;
 Libint->vrr_classes[4][5] = int_stack + 585;
 memset(int_stack,0,900*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 900;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fpgp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+900,int_stack+150,int_stack+0,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1350,int_stack+585,int_stack+360,15);
 /*--- compute (fp|gp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+2025,int_stack+1350,int_stack+900,45);
 return int_stack+2025;}
