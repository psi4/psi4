#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdgp(Libint_t*, prim_data*);

  /* Computes quartets of (fd|gp) integrals */

REALTYPE *hrr_order_fdgp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[4][4] = int_stack + 360;
 Libint->vrr_classes[4][5] = int_stack + 585;
 Libint->vrr_classes[5][4] = int_stack + 900;
 Libint->vrr_classes[5][5] = int_stack + 1215;
 memset(int_stack,0,1656*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1656;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdgp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1656,int_stack+150,int_stack+0,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2106,int_stack+585,int_stack+360,15);
 /*--- compute (fp|gp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+2781,int_stack+2106,int_stack+1656,45);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4131,int_stack+1215,int_stack+900,21);
 /*--- compute (gp|gp) ---*/
 hrr1_build_gp(Libint->AB,int_stack+0,int_stack+4131,int_stack+2106,45);
 /*--- compute (fd|gp) ---*/
 hrr1_build_fd(Libint->AB,int_stack+4131,int_stack+0,int_stack+2781,45);
 return int_stack+4131;}
