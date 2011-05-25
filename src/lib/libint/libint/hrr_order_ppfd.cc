#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppfd(Libint_t*, prim_data*);

  /* Computes quartets of (pp|fd) integrals */

REALTYPE *hrr_order_ppfd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][3] = int_stack + 0;
 Libint->vrr_classes[1][4] = int_stack + 30;
 Libint->vrr_classes[1][5] = int_stack + 75;
 Libint->vrr_classes[2][3] = int_stack + 138;
 Libint->vrr_classes[2][4] = int_stack + 198;
 Libint->vrr_classes[2][5] = int_stack + 288;
 memset(int_stack,0,414*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 414;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppfd(Libint, Data);
   Data++;
 }
 /*--- compute (p0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+414,int_stack+30,int_stack+0,3);
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+504,int_stack+75,int_stack+30,3);
 /*--- compute (p0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+639,int_stack+504,int_stack+414,3);
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+414,int_stack+198,int_stack+138,6);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+819,int_stack+288,int_stack+198,6);
 /*--- compute (d0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+0,int_stack+819,int_stack+414,6);
 /*--- compute (pp|fd) ---*/
 hrr1_build_pp(Libint->AB,int_stack+819,int_stack+0,int_stack+639,60);
 return int_stack+819;}
