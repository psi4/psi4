#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hdhd(Libint_t*, prim_data*);

  /* Computes quartets of (hd|hd) integrals */

REALTYPE *hrr_order_hdhd(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[7][5] = int_stack + 4165;
 Libint->vrr_classes[7][6] = int_stack + 4921;
 Libint->vrr_classes[7][7] = int_stack + 5929;
 memset(int_stack,0,7225*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 7225;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hdhd(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7225,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+8548,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10312,int_stack+8548,int_stack+7225,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7225,int_stack+2373,int_stack+1785,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+3157,int_stack+2373,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+12958,int_stack+0,int_stack+7225,28);
 /*--- compute (hp|hd) ---*/
 hrr1_build_hp(Libint->AB,int_stack+16486,int_stack+12958,int_stack+10312,126);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7225,int_stack+4921,int_stack+4165,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9493,int_stack+5929,int_stack+4921,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+9493,int_stack+7225,36);
 /*--- compute (ip|hd) ---*/
 hrr1_build_ip(Libint->AB,int_stack+24424,int_stack+0,int_stack+12958,126);
 /*--- compute (hd|hd) ---*/
 hrr1_build_hd(Libint->AB,int_stack+0,int_stack+24424,int_stack+16486,126);
 return int_stack+0;}
