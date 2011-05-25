#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gpfd(Libint_t*, prim_data*);

  /* Computes quartets of (gp|fd) integrals */

REALTYPE *hrr_order_gpfd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][3] = int_stack + 0;
 Libint->vrr_classes[4][4] = int_stack + 150;
 Libint->vrr_classes[4][5] = int_stack + 375;
 Libint->vrr_classes[5][3] = int_stack + 690;
 Libint->vrr_classes[5][4] = int_stack + 900;
 Libint->vrr_classes[5][5] = int_stack + 1215;
 memset(int_stack,0,1656*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1656;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gpfd(Libint, Data);
   Data++;
 }
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1656,int_stack+150,int_stack+0,15);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2106,int_stack+375,int_stack+150,15);
 /*--- compute (g0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+2781,int_stack+2106,int_stack+1656,15);
 /*--- compute (h0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1656,int_stack+900,int_stack+690,21);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3681,int_stack+1215,int_stack+900,21);
 /*--- compute (h0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+0,int_stack+3681,int_stack+1656,21);
 /*--- compute (gp|fd) ---*/
 hrr1_build_gp(Libint->AB,int_stack+3681,int_stack+0,int_stack+2781,60);
 return int_stack+3681;}
