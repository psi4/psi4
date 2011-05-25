#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hphd(Libint_t*, prim_data*);

  /* Computes quartets of (hp|hd) integrals */

REALTYPE *hrr_order_hphd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 Libint->vrr_classes[5][7] = int_stack + 1029;
 Libint->vrr_classes[6][5] = int_stack + 1785;
 Libint->vrr_classes[6][6] = int_stack + 2373;
 Libint->vrr_classes[6][7] = int_stack + 3157;
 memset(int_stack,0,4165*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4165;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hphd(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4165,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+5488,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7252,int_stack+5488,int_stack+4165,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4165,int_stack+2373,int_stack+1785,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+3157,int_stack+2373,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9898,int_stack+0,int_stack+4165,28);
 /*--- compute (hp|hd) ---*/
 hrr1_build_hp(Libint->AB,int_stack+13426,int_stack+9898,int_stack+7252,126);
 return int_stack+13426;}
