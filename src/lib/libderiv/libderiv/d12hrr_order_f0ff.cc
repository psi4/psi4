#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_f0ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|ff) integrals */

void d12hrr_order_f0ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][6][11] = int_stack + 0;
 Libderiv->deriv_classes[3][6][10] = int_stack + 280;
 Libderiv->deriv_classes[3][6][9] = int_stack + 560;
 Libderiv->deriv_classes[3][6][8] = int_stack + 840;
 Libderiv->deriv_classes[3][6][7] = int_stack + 1120;
 Libderiv->dvrr_classes[3][5] = int_stack + 1400;
 Libderiv->deriv_classes[3][6][6] = int_stack + 1610;
 Libderiv->deriv_classes[3][6][2] = int_stack + 1890;
 Libderiv->deriv_classes[3][6][1] = int_stack + 2170;
 Libderiv->deriv_classes[3][6][0] = int_stack + 2450;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 2730;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 2830;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 2980;
 Libderiv->deriv2_classes[3][6][143] = int_stack + 3190;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 3470;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 3570;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 3720;
 Libderiv->deriv2_classes[3][6][131] = int_stack + 3930;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 4210;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 4310;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 4460;
 Libderiv->deriv2_classes[3][6][130] = int_stack + 4670;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 4950;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 5050;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 5200;
 Libderiv->deriv2_classes[3][6][119] = int_stack + 5410;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 5690;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 5790;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 5940;
 Libderiv->deriv2_classes[3][6][118] = int_stack + 6150;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 6430;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 6530;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 6680;
 Libderiv->deriv2_classes[3][6][117] = int_stack + 6890;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 7170;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 7270;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 7420;
 Libderiv->deriv2_classes[3][6][107] = int_stack + 7630;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 7910;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 8010;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 8160;
 Libderiv->deriv2_classes[3][6][106] = int_stack + 8370;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 8650;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 8750;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 8900;
 Libderiv->deriv2_classes[3][6][105] = int_stack + 9110;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 9390;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 9490;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 9640;
 Libderiv->deriv2_classes[3][6][104] = int_stack + 9850;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 10130;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 10230;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 10380;
 Libderiv->deriv2_classes[3][6][95] = int_stack + 10590;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 10870;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 10970;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 11120;
 Libderiv->deriv2_classes[3][6][94] = int_stack + 11330;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 11610;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 11710;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 11860;
 Libderiv->deriv2_classes[3][6][93] = int_stack + 12070;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 12350;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 12450;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 12600;
 Libderiv->deriv2_classes[3][6][92] = int_stack + 12810;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 13090;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 13190;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 13340;
 Libderiv->deriv2_classes[3][6][91] = int_stack + 13550;
 Libderiv->deriv_classes[3][3][11] = int_stack + 13830;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 13930;
 Libderiv->deriv_classes[3][4][11] = int_stack + 14030;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 14180;
 Libderiv->deriv_classes[3][5][11] = int_stack + 14330;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 14540;
 Libderiv->deriv2_classes[3][6][83] = int_stack + 14750;
 Libderiv->deriv_classes[3][3][10] = int_stack + 15030;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 15130;
 Libderiv->deriv_classes[3][4][10] = int_stack + 15230;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 15380;
 Libderiv->deriv_classes[3][5][10] = int_stack + 15530;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 15740;
 Libderiv->deriv2_classes[3][6][82] = int_stack + 15950;
 Libderiv->deriv_classes[3][3][9] = int_stack + 16230;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 16330;
 Libderiv->deriv_classes[3][4][9] = int_stack + 16430;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 16580;
 Libderiv->deriv_classes[3][5][9] = int_stack + 16730;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 16940;
 Libderiv->deriv2_classes[3][6][81] = int_stack + 17150;
 Libderiv->deriv_classes[3][3][8] = int_stack + 17430;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 17530;
 Libderiv->deriv_classes[3][4][8] = int_stack + 17630;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 17780;
 Libderiv->deriv_classes[3][5][8] = int_stack + 17930;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 18140;
 Libderiv->deriv2_classes[3][6][80] = int_stack + 18350;
 Libderiv->deriv_classes[3][3][7] = int_stack + 18630;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 18730;
 Libderiv->deriv_classes[3][4][7] = int_stack + 18830;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 18980;
 Libderiv->deriv_classes[3][5][7] = int_stack + 19130;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 19340;
 Libderiv->deriv2_classes[3][6][79] = int_stack + 19550;
 Libderiv->dvrr_classes[3][3] = int_stack + 19830;
 Libderiv->deriv_classes[3][3][6] = int_stack + 19930;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 20030;
 Libderiv->dvrr_classes[3][4] = int_stack + 20130;
 Libderiv->deriv_classes[3][4][6] = int_stack + 20280;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 20430;
 Libderiv->deriv_classes[3][5][6] = int_stack + 20580;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 20790;
 Libderiv->deriv2_classes[3][6][78] = int_stack + 21000;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 21280;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 21380;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 21530;
 Libderiv->deriv2_classes[3][6][35] = int_stack + 21740;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 22020;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 22120;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 22270;
 Libderiv->deriv2_classes[3][6][34] = int_stack + 22480;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 22760;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 22860;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 23010;
 Libderiv->deriv2_classes[3][6][33] = int_stack + 23220;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 23500;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 23600;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 23750;
 Libderiv->deriv2_classes[3][6][32] = int_stack + 23960;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 24240;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 24340;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 24490;
 Libderiv->deriv2_classes[3][6][31] = int_stack + 24700;
 Libderiv->deriv_classes[3][3][2] = int_stack + 24980;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 25080;
 Libderiv->deriv_classes[3][4][2] = int_stack + 25180;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 25330;
 Libderiv->deriv_classes[3][5][2] = int_stack + 25480;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 25690;
 Libderiv->deriv2_classes[3][6][30] = int_stack + 25900;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 26180;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 26280;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 26430;
 Libderiv->deriv2_classes[3][6][26] = int_stack + 26640;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 26920;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 27020;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 27170;
 Libderiv->deriv2_classes[3][6][23] = int_stack + 27380;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 27660;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 27760;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 27910;
 Libderiv->deriv2_classes[3][6][22] = int_stack + 28120;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 28400;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 28500;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 28650;
 Libderiv->deriv2_classes[3][6][21] = int_stack + 28860;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 29140;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 29240;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 29390;
 Libderiv->deriv2_classes[3][6][20] = int_stack + 29600;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 29880;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 29980;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 30130;
 Libderiv->deriv2_classes[3][6][19] = int_stack + 30340;
 Libderiv->deriv_classes[3][3][1] = int_stack + 30620;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 30720;
 Libderiv->deriv_classes[3][4][1] = int_stack + 30820;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 30970;
 Libderiv->deriv_classes[3][5][1] = int_stack + 31120;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 31330;
 Libderiv->deriv2_classes[3][6][18] = int_stack + 31540;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 31820;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 31920;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 32070;
 Libderiv->deriv2_classes[3][6][14] = int_stack + 32280;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 32560;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 32660;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 32810;
 Libderiv->deriv2_classes[3][6][13] = int_stack + 33020;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 33300;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 33400;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 33550;
 Libderiv->deriv2_classes[3][6][11] = int_stack + 33760;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 34040;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 34140;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 34290;
 Libderiv->deriv2_classes[3][6][10] = int_stack + 34500;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 34780;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 34880;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 35030;
 Libderiv->deriv2_classes[3][6][9] = int_stack + 35240;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 35520;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 35620;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 35770;
 Libderiv->deriv2_classes[3][6][8] = int_stack + 35980;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 36260;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 36360;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 36510;
 Libderiv->deriv2_classes[3][6][7] = int_stack + 36720;
 Libderiv->deriv_classes[3][3][0] = int_stack + 37000;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 37100;
 Libderiv->deriv_classes[3][4][0] = int_stack + 37200;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 37350;
 Libderiv->deriv_classes[3][5][0] = int_stack + 37500;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 37710;
 Libderiv->deriv2_classes[3][6][6] = int_stack + 37920;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 38200;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 38300;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 38450;
 Libderiv->deriv2_classes[3][6][2] = int_stack + 38660;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 38940;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 39040;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 39190;
 Libderiv->deriv2_classes[3][6][1] = int_stack + 39400;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 39680;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 39780;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 39930;
 Libderiv->deriv2_classes[3][6][0] = int_stack + 40140;
 memset(int_stack,0,323360);

 Libderiv->dvrr_stack = int_stack + 93470;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_f0ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40420,int_stack+20130,int_stack+19830,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+40720,int_stack+1400,int_stack+20130,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+41170,int_stack+40720,int_stack+40420,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41770,int_stack+14030,int_stack+13830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19830,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42070,int_stack+14330,int_stack+14030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42520,int_stack+42070,int_stack+41770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43120,int_stack+0,int_stack+14330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+43750,int_stack+43120,int_stack+42070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40720,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43120,int_stack+15230,int_stack+15030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19830, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44650,int_stack+15530,int_stack+15230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45100,int_stack+44650,int_stack+43120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45700,int_stack+280,int_stack+15530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+46330,int_stack+45700,int_stack+44650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40720, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45700,int_stack+16430,int_stack+16230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19830, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+16730,int_stack+16430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+47230,int_stack+0,int_stack+45700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47830,int_stack+560,int_stack+16730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48460,int_stack+47830,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40720, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47830,int_stack+17630,int_stack+17430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+17930,int_stack+17630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+49810,int_stack+49360,int_stack+47830, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+50410,int_stack+840,int_stack+17930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+51040,int_stack+50410,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+50410,int_stack+18830,int_stack+18630, 0.0, zero_stack, 1.0, int_stack+19830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+19130,int_stack+18830, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+51940,int_stack+450,int_stack+50410, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52540,int_stack+1120,int_stack+19130, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+53170,int_stack+52540,int_stack+450, 0.0, zero_stack, 1.0, int_stack+40720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+52540,int_stack+20280,int_stack+19930, 1.0, int_stack+19830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+20580,int_stack+20280, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54070,int_stack+900,int_stack+52540, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+54670,int_stack+1610,int_stack+20580, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55300,int_stack+54670,int_stack+900, 1.0, int_stack+40720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+54670,int_stack+25180,int_stack+24980,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+40420,int_stack+25480,int_stack+25180,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56200,int_stack+40420,int_stack+54670,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+56800,int_stack+1890,int_stack+25480,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+57430,int_stack+56800,int_stack+40420,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40870,int_stack+30820,int_stack+30620,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+56800,int_stack+31120,int_stack+30820,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1350,int_stack+56800,int_stack+40870,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+58330,int_stack+2170,int_stack+31120,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+58960,int_stack+58330,int_stack+56800,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+58330,int_stack+37200,int_stack+37000,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1950,int_stack+37500,int_stack+37200,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+59860,int_stack+1950,int_stack+58330,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+2450,int_stack+37500,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+61090,int_stack+60460,int_stack+1950,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+2830,int_stack+2730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13830,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+2980,int_stack+2830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14030,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+62440,int_stack+61990,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41770,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+3190,int_stack+2980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14330,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2400,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42070,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+3570,int_stack+3470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13830, 1.0, int_stack+15030,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+3720,int_stack+3570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14030, 1.0, int_stack+15230,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+63040,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41770, 1.0, int_stack+43120,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+63640,int_stack+3930,int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14330, 1.0, int_stack+15530,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3300,int_stack+63640,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42070, 1.0, int_stack+44650,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+4310,int_stack+4210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15030, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+4460,int_stack+4310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15230, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+63640,int_stack+61990,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43120, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+4670,int_stack+4460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15530, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+64240,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+44650, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+5050,int_stack+4950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13830, 0.0, zero_stack, 1.0, int_stack+16230,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+5200,int_stack+5050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14030, 0.0, zero_stack, 1.0, int_stack+16430,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65140,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41770, 0.0, zero_stack, 1.0, int_stack+45700,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4200,int_stack+5410,int_stack+5200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14330, 0.0, zero_stack, 1.0, int_stack+16730,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+65740,int_stack+4200,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42070, 0.0, zero_stack, 1.0, int_stack+0,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+5790,int_stack+5690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15030, 1.0, int_stack+16230, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+5940,int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15230, 1.0, int_stack+16430, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4200,int_stack+61990,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43120, 1.0, int_stack+45700, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+6150,int_stack+5940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15530, 1.0, int_stack+16730, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4800,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44650, 1.0, int_stack+0, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+6530,int_stack+6430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16230, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+6680,int_stack+6530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16430, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5700,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45700, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66640,int_stack+6890,int_stack+6680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16730, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67270,int_stack+66640,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+7270,int_stack+7170, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13830, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17430,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+7420,int_stack+7270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17630,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+66640,int_stack+61990,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41770, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47830,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+7630,int_stack+7420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17930,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6300,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49360,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+8010,int_stack+7910, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15030, 0.0, zero_stack, 1.0, int_stack+17430, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+8160,int_stack+8010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15230, 0.0, zero_stack, 1.0, int_stack+17630, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7200,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43120, 0.0, zero_stack, 1.0, int_stack+47830, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+68170,int_stack+8370,int_stack+8160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15530, 0.0, zero_stack, 1.0, int_stack+17930, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+68800,int_stack+68170,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44650, 0.0, zero_stack, 1.0, int_stack+49360, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+8750,int_stack+8650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16230, 1.0, int_stack+17430, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+8900,int_stack+8750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16430, 1.0, int_stack+17630, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+68170,int_stack+61990,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45700, 1.0, int_stack+47830, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+9110,int_stack+8900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16730, 1.0, int_stack+17930, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7800,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+49360, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+9490,int_stack+9390, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+9640,int_stack+9490, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8700,int_stack+60460,int_stack+61990, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+47830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69700,int_stack+9850,int_stack+9640, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+70330,int_stack+69700,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+10230,int_stack+10130, 0.0, zero_stack, 1.0, int_stack+13830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18630,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+10380,int_stack+10230, 0.0, zero_stack, 1.0, int_stack+14030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18830,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+69700,int_stack+61990,int_stack+60460, 0.0, zero_stack, 1.0, int_stack+41770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50410,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+10590,int_stack+10380, 0.0, zero_stack, 1.0, int_stack+14330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19130,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9300,int_stack+60460,int_stack+61990, 0.0, zero_stack, 1.0, int_stack+42070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+10970,int_stack+10870, 0.0, zero_stack, 1.0, int_stack+15030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18630, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+11120,int_stack+10970, 0.0, zero_stack, 1.0, int_stack+15230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18830, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10200,int_stack+60460,int_stack+61990, 0.0, zero_stack, 1.0, int_stack+43120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50410, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+71230,int_stack+11330,int_stack+11120, 0.0, zero_stack, 1.0, int_stack+15530, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19130, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+71860,int_stack+71230,int_stack+60460, 0.0, zero_stack, 1.0, int_stack+44650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+11710,int_stack+11610, 0.0, zero_stack, 1.0, int_stack+16230, 0.0, zero_stack, 1.0, int_stack+18630, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+11860,int_stack+11710, 0.0, zero_stack, 1.0, int_stack+16430, 0.0, zero_stack, 1.0, int_stack+18830, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+71230,int_stack+61990,int_stack+60460, 0.0, zero_stack, 1.0, int_stack+45700, 0.0, zero_stack, 1.0, int_stack+50410, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+12070,int_stack+11860, 0.0, zero_stack, 1.0, int_stack+16730, 0.0, zero_stack, 1.0, int_stack+19130, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10800,int_stack+60460,int_stack+61990, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+12450,int_stack+12350, 0.0, zero_stack, 1.0, int_stack+17430, 1.0, int_stack+18630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+12600,int_stack+12450, 0.0, zero_stack, 1.0, int_stack+17630, 1.0, int_stack+18830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11700,int_stack+60460,int_stack+61990, 0.0, zero_stack, 1.0, int_stack+47830, 1.0, int_stack+50410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+72760,int_stack+12810,int_stack+12600, 0.0, zero_stack, 1.0, int_stack+17930, 1.0, int_stack+19130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+73390,int_stack+72760,int_stack+60460, 0.0, zero_stack, 1.0, int_stack+49360, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+13190,int_stack+13090, 0.0, zero_stack, 2.0, int_stack+18630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+13340,int_stack+13190, 0.0, zero_stack, 2.0, int_stack+18830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+72760,int_stack+61990,int_stack+60460, 0.0, zero_stack, 2.0, int_stack+50410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+13550,int_stack+13340, 0.0, zero_stack, 2.0, int_stack+19130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12300,int_stack+60460,int_stack+61990, 0.0, zero_stack, 2.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+14180,int_stack+13930, 1.0, int_stack+13830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19930,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60460,int_stack+14540,int_stack+14180, 1.0, int_stack+14030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20280,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13200,int_stack+60460,int_stack+61990, 1.0, int_stack+41770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52540,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+74290,int_stack+14750,int_stack+14540, 1.0, int_stack+14330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20580,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13800,int_stack+74290,int_stack+60460, 1.0, int_stack+42070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60460,int_stack+15380,int_stack+15130, 1.0, int_stack+15030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19930, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+15740,int_stack+15380, 1.0, int_stack+15230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20280, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+74290,int_stack+61990,int_stack+60460, 1.0, int_stack+43120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52540, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43120,int_stack+15950,int_stack+15740, 1.0, int_stack+15530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20580, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14700,int_stack+43120,int_stack+61990, 1.0, int_stack+44650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44650,int_stack+16580,int_stack+16330, 1.0, int_stack+16230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+16940,int_stack+16580, 1.0, int_stack+16430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20280, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+43120,int_stack+61990,int_stack+44650, 1.0, int_stack+45700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52540, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45700,int_stack+17150,int_stack+16940, 1.0, int_stack+16730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20580, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15600,int_stack+45700,int_stack+61990, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+17780,int_stack+17530, 1.0, int_stack+17430, 0.0, zero_stack, 1.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+18140,int_stack+17780, 1.0, int_stack+17630, 0.0, zero_stack, 1.0, int_stack+20280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45700,int_stack+61990,int_stack+0, 1.0, int_stack+47830, 0.0, zero_stack, 1.0, int_stack+52540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47830,int_stack+18350,int_stack+18140, 1.0, int_stack+17930, 0.0, zero_stack, 1.0, int_stack+20580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16500,int_stack+47830,int_stack+61990, 1.0, int_stack+49360, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+18980,int_stack+18730, 1.0, int_stack+18630, 1.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61990,int_stack+19340,int_stack+18980, 1.0, int_stack+18830, 1.0, int_stack+20280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+47830,int_stack+61990,int_stack+49360, 1.0, int_stack+50410, 1.0, int_stack+52540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+50410,int_stack+19550,int_stack+19340, 1.0, int_stack+19130, 1.0, int_stack+20580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17400,int_stack+50410,int_stack+61990, 1.0, int_stack+450, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61990,int_stack+20430,int_stack+20030, 2.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+20790,int_stack+20430, 2.0, int_stack+20280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+50410,int_stack+49360,int_stack+61990, 2.0, int_stack+52540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52540,int_stack+21000,int_stack+20790, 2.0, int_stack+20580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+52540,int_stack+49360, 2.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+21380,int_stack+21280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24980,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+21530,int_stack+21380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25180,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+52540,int_stack+49360,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54670,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60460,int_stack+21740,int_stack+21530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25480,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18300,int_stack+60460,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+22120,int_stack+22020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24980, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+22270,int_stack+22120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25180, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+60460,int_stack+900,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54670, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41770,int_stack+22480,int_stack+22270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25480, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19200,int_stack+41770,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+22860,int_stack+22760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24980, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+23010,int_stack+22860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25180, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41770,int_stack+49360,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54670, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20100,int_stack+23220,int_stack+23010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25480, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20730,int_stack+20100,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+23600,int_stack+23500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+23750,int_stack+23600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20100,int_stack+900,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21630,int_stack+23960,int_stack+23750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22260,int_stack+21630,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+24340,int_stack+24240, 0.0, zero_stack, 1.0, int_stack+24980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+24490,int_stack+24340, 0.0, zero_stack, 1.0, int_stack+25180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21630,int_stack+49360,int_stack+900, 0.0, zero_stack, 1.0, int_stack+54670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23160,int_stack+24700,int_stack+24490, 0.0, zero_stack, 1.0, int_stack+25480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23790,int_stack+23160,int_stack+49360, 0.0, zero_stack, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+25330,int_stack+25080, 1.0, int_stack+24980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+25690,int_stack+25330, 1.0, int_stack+25180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23160,int_stack+900,int_stack+49360, 1.0, int_stack+54670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+54670,int_stack+25900,int_stack+25690, 1.0, int_stack+25480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24690,int_stack+54670,int_stack+900, 1.0, int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40420,int_stack+26280,int_stack+26180,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+26430,int_stack+26280,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+54670,int_stack+900,int_stack+40420,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+25590,int_stack+26640,int_stack+26430,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+74890,int_stack+25590,int_stack+900,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+27020,int_stack+26920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30620,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40420,int_stack+27170,int_stack+27020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30820,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25590,int_stack+40420,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40870,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26190,int_stack+27380,int_stack+27170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31120,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+75790,int_stack+26190,int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56800,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40420,int_stack+27760,int_stack+27660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30620, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+27910,int_stack+27760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30820, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26190,int_stack+900,int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40870, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26790,int_stack+28120,int_stack+27910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31120, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27420,int_stack+26790,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56800, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+28500,int_stack+28400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30620, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40420,int_stack+28650,int_stack+28500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30820, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26790,int_stack+40420,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40870, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+76690,int_stack+28860,int_stack+28650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31120, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+77320,int_stack+76690,int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56800, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40420,int_stack+29240,int_stack+29140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+29390,int_stack+29240, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+76690,int_stack+900,int_stack+40420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+28320,int_stack+29600,int_stack+29390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+28950,int_stack+28320,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+29980,int_stack+29880, 0.0, zero_stack, 1.0, int_stack+30620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40420,int_stack+30130,int_stack+29980, 0.0, zero_stack, 1.0, int_stack+30820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28320,int_stack+40420,int_stack+900, 0.0, zero_stack, 1.0, int_stack+40870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+78220,int_stack+30340,int_stack+30130, 0.0, zero_stack, 1.0, int_stack+31120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+78850,int_stack+78220,int_stack+40420, 0.0, zero_stack, 1.0, int_stack+56800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40420,int_stack+30970,int_stack+30720, 1.0, int_stack+30620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+31330,int_stack+30970, 1.0, int_stack+30820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+78220,int_stack+900,int_stack+40420, 1.0, int_stack+40870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+40420,int_stack+31540,int_stack+31330, 1.0, int_stack+31120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29850,int_stack+40420,int_stack+900, 1.0, int_stack+56800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+56800,int_stack+31920,int_stack+31820,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+32070,int_stack+31920,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+40420,int_stack+900,int_stack+56800,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+56800,int_stack+32280,int_stack+32070,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+30750,int_stack+56800,int_stack+900,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+32660,int_stack+32560,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+32810,int_stack+32660,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56800,int_stack+49360,int_stack+900,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+31650,int_stack+33020,int_stack+32810,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+32280,int_stack+31650,int_stack+49360,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+33400,int_stack+33300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37000,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+33550,int_stack+33400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37200,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31650,int_stack+900,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58330,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+79750,int_stack+33760,int_stack+33550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37500,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+80380,int_stack+79750,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+34140,int_stack+34040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37000, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+34290,int_stack+34140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37200, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+79750,int_stack+49360,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58330, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+33180,int_stack+34500,int_stack+34290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37500, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33810,int_stack+33180,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+34880,int_stack+34780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37000, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+35030,int_stack+34880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33180,int_stack+900,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81280,int_stack+35240,int_stack+35030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81910,int_stack+81280,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+35620,int_stack+35520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+35770,int_stack+35620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81280,int_stack+49360,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+34710,int_stack+35980,int_stack+35770, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+35340,int_stack+34710,int_stack+49360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+36360,int_stack+36260, 0.0, zero_stack, 1.0, int_stack+37000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+36510,int_stack+36360, 0.0, zero_stack, 1.0, int_stack+37200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+34710,int_stack+900,int_stack+49360, 0.0, zero_stack, 1.0, int_stack+58330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+82810,int_stack+36720,int_stack+36510, 0.0, zero_stack, 1.0, int_stack+37500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83440,int_stack+82810,int_stack+900, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+37350,int_stack+37100, 1.0, int_stack+37000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+37710,int_stack+37350, 1.0, int_stack+37200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+82810,int_stack+49360,int_stack+900, 1.0, int_stack+58330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+58330,int_stack+37920,int_stack+37710, 1.0, int_stack+37500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+36240,int_stack+58330,int_stack+49360, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+38300,int_stack+38200,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+38450,int_stack+38300,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+58330,int_stack+49360,int_stack+1950,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+37140,int_stack+38660,int_stack+38450,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+37770,int_stack+37140,int_stack+49360,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+49360,int_stack+39040,int_stack+38940,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1950,int_stack+39190,int_stack+39040,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+37140,int_stack+1950,int_stack+49360,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+84340,int_stack+39400,int_stack+39190,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+38670,int_stack+84340,int_stack+1950,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+39780,int_stack+39680,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+49360,int_stack+39930,int_stack+39780,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+84340,int_stack+49360,int_stack+1950,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+84940,int_stack+40140,int_stack+39930,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+85570,int_stack+84940,int_stack+49360,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+86470,int_stack+43750,int_stack+42520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41170,10);
     Libderiv->ABCD[11] = int_stack + 86470;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+43720,int_stack+46330,int_stack+45100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41170, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 43720;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+87470,int_stack+48460,int_stack+47230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41170, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 87470;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+48430,int_stack+51040,int_stack+49810, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 48430;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+88470,int_stack+53170,int_stack+51940, 0.0, zero_stack, 1.0, int_stack+41170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 88470;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+89470,int_stack+55300,int_stack+54070, 1.0, int_stack+41170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 89470;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+90470,int_stack+57430,int_stack+56200,10);
     Libderiv->ABCD[2] = int_stack + 90470;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+91470,int_stack+58960,int_stack+1350,10);
     Libderiv->ABCD[1] = int_stack + 91470;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+92470,int_stack+61090,int_stack+59860,10);
     Libderiv->ABCD[0] = int_stack + 92470;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+61060,int_stack+2400,int_stack+62440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42520,10);
     Libderiv->ABCD[155] = int_stack + 61060;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1950,int_stack+3300,int_stack+63040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42520, 1.0, int_stack+45100,10);
     Libderiv->ABCD[143] = int_stack + 1950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2950,int_stack+64240,int_stack+63640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45100, 0.0, zero_stack,10);
     Libderiv->ABCD[142] = int_stack + 2950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+62060,int_stack+65740,int_stack+65140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42520, 0.0, zero_stack, 1.0, int_stack+47230,10);
     Libderiv->ABCD[131] = int_stack + 62060;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+63060,int_stack+4800,int_stack+4200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45100, 1.0, int_stack+47230, 0.0, zero_stack,10);
     Libderiv->ABCD[130] = int_stack + 63060;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3950,int_stack+67270,int_stack+5700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+47230, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[129] = int_stack + 3950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4950,int_stack+6300,int_stack+66640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49810,10);
     Libderiv->ABCD[119] = int_stack + 4950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5950,int_stack+68800,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45100, 0.0, zero_stack, 1.0, int_stack+49810, 0.0, zero_stack,10);
     Libderiv->ABCD[118] = int_stack + 5950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+64060,int_stack+7800,int_stack+68170, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47230, 1.0, int_stack+49810, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[117] = int_stack + 64060;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6950,int_stack+70330,int_stack+8700, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[116] = int_stack + 6950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7950,int_stack+9300,int_stack+69700, 0.0, zero_stack, 1.0, int_stack+42520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51940,10);
     Libderiv->ABCD[107] = int_stack + 7950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8950,int_stack+71860,int_stack+10200, 0.0, zero_stack, 1.0, int_stack+45100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51940, 0.0, zero_stack,10);
     Libderiv->ABCD[106] = int_stack + 8950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+65060,int_stack+10800,int_stack+71230, 0.0, zero_stack, 1.0, int_stack+47230, 0.0, zero_stack, 1.0, int_stack+51940, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[105] = int_stack + 65060;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9950,int_stack+73390,int_stack+11700, 0.0, zero_stack, 1.0, int_stack+49810, 1.0, int_stack+51940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[104] = int_stack + 9950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10950,int_stack+12300,int_stack+72760, 0.0, zero_stack, 2.0, int_stack+51940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[103] = int_stack + 10950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11950,int_stack+13800,int_stack+13200, 1.0, int_stack+42520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54070,10);
     Libderiv->ABCD[95] = int_stack + 11950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12950,int_stack+14700,int_stack+74290, 1.0, int_stack+45100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54070, 0.0, zero_stack,10);
     Libderiv->ABCD[94] = int_stack + 12950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13950,int_stack+15600,int_stack+43120, 1.0, int_stack+47230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54070, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[93] = int_stack + 13950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14950,int_stack+16500,int_stack+45700, 1.0, int_stack+49810, 0.0, zero_stack, 1.0, int_stack+54070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[92] = int_stack + 14950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+15950,int_stack+17400,int_stack+47830, 1.0, int_stack+51940, 1.0, int_stack+54070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[91] = int_stack + 15950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16950,int_stack+0,int_stack+50410, 2.0, int_stack+54070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[90] = int_stack + 16950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+18300,int_stack+52540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56200,10);
     Libderiv->ABCD[47] = int_stack + 0;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+17950,int_stack+19200,int_stack+60460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56200, 0.0, zero_stack,10);
     Libderiv->ABCD[46] = int_stack + 17950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18950,int_stack+20730,int_stack+41770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56200, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[45] = int_stack + 18950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+41020,int_stack+22260,int_stack+20100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[44] = int_stack + 41020;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+19950,int_stack+23790,int_stack+21630, 0.0, zero_stack, 1.0, int_stack+56200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[43] = int_stack + 19950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+20950,int_stack+24690,int_stack+23160, 1.0, int_stack+56200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[42] = int_stack + 20950;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+21950,int_stack+74890,int_stack+54670,10);
     Libderiv->ABCD[38] = int_stack + 21950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+22950,int_stack+75790,int_stack+25590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350,10);
     Libderiv->ABCD[35] = int_stack + 22950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23950,int_stack+27420,int_stack+26190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack,10);
     Libderiv->ABCD[34] = int_stack + 23950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+24950,int_stack+77320,int_stack+26790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[33] = int_stack + 24950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+25950,int_stack+28950,int_stack+76690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[32] = int_stack + 25950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26950,int_stack+78850,int_stack+28320, 0.0, zero_stack, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[31] = int_stack + 26950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27950,int_stack+29850,int_stack+78220, 1.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[30] = int_stack + 27950;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+28950,int_stack+30750,int_stack+40420,10);
     Libderiv->ABCD[26] = int_stack + 28950;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+29950,int_stack+32280,int_stack+56800,10);
     Libderiv->ABCD[25] = int_stack + 29950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+42020,int_stack+80380,int_stack+31650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59860,10);
     Libderiv->ABCD[23] = int_stack + 42020;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+30950,int_stack+33810,int_stack+79750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59860, 0.0, zero_stack,10);
     Libderiv->ABCD[22] = int_stack + 30950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+31950,int_stack+81910,int_stack+33180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59860, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[21] = int_stack + 31950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+32950,int_stack+35340,int_stack+81280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[20] = int_stack + 32950;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+39570,int_stack+83440,int_stack+34710, 0.0, zero_stack, 1.0, int_stack+59860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[19] = int_stack + 39570;
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+33950,int_stack+36240,int_stack+82810, 1.0, int_stack+59860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[18] = int_stack + 33950;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+34950,int_stack+37770,int_stack+58330,10);
     Libderiv->ABCD[14] = int_stack + 34950;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+35950,int_stack+38670,int_stack+37140,10);
     Libderiv->ABCD[13] = int_stack + 35950;
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+36950,int_stack+85570,int_stack+84340,10);
     Libderiv->ABCD[12] = int_stack + 36950;

}
