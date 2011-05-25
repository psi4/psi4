#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0gd(Libint_t*, prim_data*);

  /* Computes quartets of (f0|gd) integrals */

REALTYPE *hrr_order_f0gd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[3][6] = int_stack + 360;
 memset(int_stack,0,640*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 640;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0gd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+640,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1090,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1720,int_stack+1090,int_stack+640,10);
 return int_stack+1720;}
