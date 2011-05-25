#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdhd(Libint_t*, prim_data*);

  /* Computes quartets of (fd|hd) integrals */

REALTYPE *hrr_order_fdhd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[3][7] = int_stack + 490;
 Libint->vrr_classes[4][5] = int_stack + 850;
 Libint->vrr_classes[4][6] = int_stack + 1165;
 Libint->vrr_classes[4][7] = int_stack + 1585;
 Libint->vrr_classes[5][5] = int_stack + 2125;
 Libint->vrr_classes[5][6] = int_stack + 2566;
 Libint->vrr_classes[5][7] = int_stack + 3154;
 memset(int_stack,0,3910*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3910;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdhd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3910,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4540,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5380,int_stack+4540,int_stack+3910,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3910,int_stack+1165,int_stack+850,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6640,int_stack+1585,int_stack+1165,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+6640,int_stack+3910,15);
 /*--- compute (fp|hd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+6640,int_stack+0,int_stack+5380,126);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3910,int_stack+2566,int_stack+2125,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10420,int_stack+3154,int_stack+2566,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+12184,int_stack+10420,int_stack+3910,21);
 /*--- compute (gp|hd) ---*/
 hrr1_build_gp(Libint->AB,int_stack+14830,int_stack+12184,int_stack+0,126);
 /*--- compute (fd|hd) ---*/
 hrr1_build_fd(Libint->AB,int_stack+20500,int_stack+14830,int_stack+6640,126);
 return int_stack+20500;}
