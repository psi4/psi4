#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fpff(Libint_t*, prim_data*);

  /* Computes quartets of (fp|ff) integrals */

REALTYPE *hrr_order_fpff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][3] = int_stack + 0;
 Libint->vrr_classes[3][4] = int_stack + 100;
 Libint->vrr_classes[3][5] = int_stack + 250;
 Libint->vrr_classes[3][6] = int_stack + 460;
 Libint->vrr_classes[4][3] = int_stack + 740;
 Libint->vrr_classes[4][4] = int_stack + 890;
 Libint->vrr_classes[4][5] = int_stack + 1115;
 Libint->vrr_classes[4][6] = int_stack + 1430;
 memset(int_stack,0,1850*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1850;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fpff(Libint, Data);
   Data++;
 }
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1850,int_stack+100,int_stack+0,10);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2150,int_stack+250,int_stack+100,10);
 /*--- compute (f0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+2600,int_stack+2150,int_stack+1850,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3200,int_stack+460,int_stack+250,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3830,int_stack+3200,int_stack+2150,10);
 /*--- compute (f0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+4730,int_stack+3830,int_stack+2600,10);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1850,int_stack+890,int_stack+740,15);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2300,int_stack+1115,int_stack+890,15);
 /*--- compute (g0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+2975,int_stack+2300,int_stack+1850,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+1430,int_stack+1115,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+945,int_stack+0,int_stack+2300,15);
 /*--- compute (g0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+5730,int_stack+945,int_stack+2975,15);
 /*--- compute (fp|ff) ---*/
 hrr1_build_fp(Libint->AB,int_stack+0,int_stack+5730,int_stack+4730,100);
 return int_stack+0;}
