#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fpfp(Libint_t*, prim_data*);

  /* Computes quartets of (fp|fp) integrals */

REALTYPE *hrr_order_fpfp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][3] = int_stack + 0;
 Libint->vrr_classes[3][4] = int_stack + 100;
 Libint->vrr_classes[4][3] = int_stack + 250;
 Libint->vrr_classes[4][4] = int_stack + 400;
 memset(int_stack,0,625*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 625;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fpfp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+625,int_stack+100,int_stack+0,10);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+925,int_stack+400,int_stack+250,15);
 /*--- compute (fp|fp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+1375,int_stack+925,int_stack+625,30);
 return int_stack+1375;}
