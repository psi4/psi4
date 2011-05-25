#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gdhh(Libint_t*, prim_data*);

  /* Computes quartets of (gd|hh) integrals */

REALTYPE *hrr_order_gdhh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[4][9] = int_stack + 1950;
 Libint->vrr_classes[4][10] = int_stack + 2775;
 Libint->vrr_classes[5][5] = int_stack + 3765;
 Libint->vrr_classes[5][6] = int_stack + 4206;
 Libint->vrr_classes[5][7] = int_stack + 4794;
 Libint->vrr_classes[5][8] = int_stack + 5550;
 Libint->vrr_classes[5][9] = int_stack + 6495;
 Libint->vrr_classes[5][10] = int_stack + 7650;
 Libint->vrr_classes[6][5] = int_stack + 9036;
 Libint->vrr_classes[6][6] = int_stack + 9624;
 Libint->vrr_classes[6][7] = int_stack + 10408;
 Libint->vrr_classes[6][8] = int_stack + 11416;
 Libint->vrr_classes[6][9] = int_stack + 12676;
 Libint->vrr_classes[6][10] = int_stack + 14216;
 memset(int_stack,0,16064*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 16064;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gdhh(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+16064,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+17009,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+18269,int_stack+17009,int_stack+16064,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+20159,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+21779,int_stack+20159,int_stack+17009,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+24299,int_stack+21779,int_stack+18269,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+16064,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+27449,int_stack+16064,int_stack+20159,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+30689,int_stack+27449,int_stack+21779,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+18089,int_stack+30689,int_stack+24299,15);
 /*--- compute (g0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+22814,int_stack+2775,int_stack+1950,15);
 /*--- compute (g0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+34889,int_stack+22814,int_stack+16064,15);
 /*--- compute (g0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+38939,int_stack+34889,int_stack+27449,15);
 /*--- compute (g0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+22814,int_stack+38939,int_stack+30689,15);
 /*--- compute (g0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+29114,int_stack+22814,int_stack+18089,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+16064,int_stack+4206,int_stack+3765,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+17387,int_stack+4794,int_stack+4206,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+19151,int_stack+17387,int_stack+16064,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+21797,int_stack+5550,int_stack+4794,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24065,int_stack+21797,int_stack+17387,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+35729,int_stack+24065,int_stack+19151,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+16064,int_stack+6495,int_stack+5550,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+40139,int_stack+16064,int_stack+21797,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+40139,int_stack+24065,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+18899,int_stack+0,int_stack+35729,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+35729,int_stack+7650,int_stack+6495,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+44675,int_stack+35729,int_stack+16064,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+50345,int_stack+44675,int_stack+40139,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+35729,int_stack+50345,int_stack+0,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+44549,int_stack+35729,int_stack+18899,21);
 /*--- compute (gp|hh) ---*/
 hrr1_build_gp(Libint->AB,int_stack+53810,int_stack+44549,int_stack+29114,441);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+9624,int_stack+9036,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+10408,int_stack+9624,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4116,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+16064,int_stack+11416,int_stack+10408,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+19088,int_stack+16064,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+23792,int_stack+19088,int_stack+4116,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+12676,int_stack+11416,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+3780,int_stack+0,int_stack+16064,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+29672,int_stack+3780,int_stack+19088,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+73655,int_stack+29672,int_stack+23792,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+16064,int_stack+14216,int_stack+12676,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+20684,int_stack+16064,int_stack+0,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+9828,int_stack+20684,int_stack+3780,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+82475,int_stack+9828,int_stack+29672,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+82475,int_stack+73655,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+12348,int_stack+0,int_stack+44549,441);
 /*--- compute (gd|hh) ---*/
 hrr1_build_gd(Libint->AB,int_stack+73655,int_stack+12348,int_stack+53810,441);
 return int_stack+73655;}
