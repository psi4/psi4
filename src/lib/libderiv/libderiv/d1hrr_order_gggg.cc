#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gggg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gg|gg) integrals */

void d1hrr_order_gggg(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[4][7][11] = int_stack + 960;
 Libderiv->deriv_classes[4][8][11] = int_stack + 1500;
 Libderiv->deriv_classes[5][4][11] = int_stack + 2175;
 Libderiv->deriv_classes[5][5][11] = int_stack + 2490;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2931;
 Libderiv->deriv_classes[5][7][11] = int_stack + 3519;
 Libderiv->deriv_classes[5][8][11] = int_stack + 4275;
 Libderiv->deriv_classes[6][4][11] = int_stack + 5220;
 Libderiv->deriv_classes[6][5][11] = int_stack + 5640;
 Libderiv->deriv_classes[6][6][11] = int_stack + 6228;
 Libderiv->deriv_classes[6][7][11] = int_stack + 7012;
 Libderiv->deriv_classes[6][8][11] = int_stack + 8020;
 Libderiv->deriv_classes[7][4][11] = int_stack + 9280;
 Libderiv->deriv_classes[7][5][11] = int_stack + 9820;
 Libderiv->deriv_classes[7][6][11] = int_stack + 10576;
 Libderiv->deriv_classes[7][7][11] = int_stack + 11584;
 Libderiv->deriv_classes[7][8][11] = int_stack + 12880;
 Libderiv->deriv_classes[8][4][11] = int_stack + 14500;
 Libderiv->deriv_classes[8][5][11] = int_stack + 15175;
 Libderiv->deriv_classes[8][6][11] = int_stack + 16120;
 Libderiv->deriv_classes[8][7][11] = int_stack + 17380;
 Libderiv->deriv_classes[8][8][11] = int_stack + 19000;
 Libderiv->deriv_classes[4][4][10] = int_stack + 21025;
 Libderiv->deriv_classes[4][5][10] = int_stack + 21250;
 Libderiv->deriv_classes[4][6][10] = int_stack + 21565;
 Libderiv->deriv_classes[4][7][10] = int_stack + 21985;
 Libderiv->deriv_classes[4][8][10] = int_stack + 22525;
 Libderiv->deriv_classes[5][4][10] = int_stack + 23200;
 Libderiv->deriv_classes[5][5][10] = int_stack + 23515;
 Libderiv->deriv_classes[5][6][10] = int_stack + 23956;
 Libderiv->deriv_classes[5][7][10] = int_stack + 24544;
 Libderiv->deriv_classes[5][8][10] = int_stack + 25300;
 Libderiv->deriv_classes[6][4][10] = int_stack + 26245;
 Libderiv->deriv_classes[6][5][10] = int_stack + 26665;
 Libderiv->deriv_classes[6][6][10] = int_stack + 27253;
 Libderiv->deriv_classes[6][7][10] = int_stack + 28037;
 Libderiv->deriv_classes[6][8][10] = int_stack + 29045;
 Libderiv->deriv_classes[7][4][10] = int_stack + 30305;
 Libderiv->deriv_classes[7][5][10] = int_stack + 30845;
 Libderiv->deriv_classes[7][6][10] = int_stack + 31601;
 Libderiv->deriv_classes[7][7][10] = int_stack + 32609;
 Libderiv->deriv_classes[7][8][10] = int_stack + 33905;
 Libderiv->deriv_classes[8][4][10] = int_stack + 35525;
 Libderiv->deriv_classes[8][5][10] = int_stack + 36200;
 Libderiv->deriv_classes[8][6][10] = int_stack + 37145;
 Libderiv->deriv_classes[8][7][10] = int_stack + 38405;
 Libderiv->deriv_classes[8][8][10] = int_stack + 40025;
 Libderiv->deriv_classes[4][4][9] = int_stack + 42050;
 Libderiv->deriv_classes[4][5][9] = int_stack + 42275;
 Libderiv->deriv_classes[4][6][9] = int_stack + 42590;
 Libderiv->deriv_classes[4][7][9] = int_stack + 43010;
 Libderiv->deriv_classes[4][8][9] = int_stack + 43550;
 Libderiv->deriv_classes[5][4][9] = int_stack + 44225;
 Libderiv->deriv_classes[5][5][9] = int_stack + 44540;
 Libderiv->deriv_classes[5][6][9] = int_stack + 44981;
 Libderiv->deriv_classes[5][7][9] = int_stack + 45569;
 Libderiv->deriv_classes[5][8][9] = int_stack + 46325;
 Libderiv->deriv_classes[6][4][9] = int_stack + 47270;
 Libderiv->deriv_classes[6][5][9] = int_stack + 47690;
 Libderiv->deriv_classes[6][6][9] = int_stack + 48278;
 Libderiv->deriv_classes[6][7][9] = int_stack + 49062;
 Libderiv->deriv_classes[6][8][9] = int_stack + 50070;
 Libderiv->deriv_classes[7][4][9] = int_stack + 51330;
 Libderiv->deriv_classes[7][5][9] = int_stack + 51870;
 Libderiv->deriv_classes[7][6][9] = int_stack + 52626;
 Libderiv->deriv_classes[7][7][9] = int_stack + 53634;
 Libderiv->deriv_classes[7][8][9] = int_stack + 54930;
 Libderiv->deriv_classes[8][4][9] = int_stack + 56550;
 Libderiv->deriv_classes[8][5][9] = int_stack + 57225;
 Libderiv->deriv_classes[8][6][9] = int_stack + 58170;
 Libderiv->deriv_classes[8][7][9] = int_stack + 59430;
 Libderiv->deriv_classes[8][8][9] = int_stack + 61050;
 Libderiv->deriv_classes[4][4][8] = int_stack + 63075;
 Libderiv->deriv_classes[4][5][8] = int_stack + 63300;
 Libderiv->deriv_classes[4][6][8] = int_stack + 63615;
 Libderiv->deriv_classes[4][7][8] = int_stack + 64035;
 Libderiv->deriv_classes[4][8][8] = int_stack + 64575;
 Libderiv->deriv_classes[5][4][8] = int_stack + 65250;
 Libderiv->deriv_classes[5][5][8] = int_stack + 65565;
 Libderiv->deriv_classes[5][6][8] = int_stack + 66006;
 Libderiv->deriv_classes[5][7][8] = int_stack + 66594;
 Libderiv->deriv_classes[5][8][8] = int_stack + 67350;
 Libderiv->deriv_classes[6][4][8] = int_stack + 68295;
 Libderiv->deriv_classes[6][5][8] = int_stack + 68715;
 Libderiv->deriv_classes[6][6][8] = int_stack + 69303;
 Libderiv->deriv_classes[6][7][8] = int_stack + 70087;
 Libderiv->deriv_classes[6][8][8] = int_stack + 71095;
 Libderiv->deriv_classes[7][4][8] = int_stack + 72355;
 Libderiv->deriv_classes[7][5][8] = int_stack + 72895;
 Libderiv->deriv_classes[7][6][8] = int_stack + 73651;
 Libderiv->deriv_classes[7][7][8] = int_stack + 74659;
 Libderiv->deriv_classes[7][8][8] = int_stack + 75955;
 Libderiv->deriv_classes[8][4][8] = int_stack + 77575;
 Libderiv->deriv_classes[8][5][8] = int_stack + 78250;
 Libderiv->deriv_classes[8][6][8] = int_stack + 79195;
 Libderiv->deriv_classes[8][7][8] = int_stack + 80455;
 Libderiv->deriv_classes[8][8][8] = int_stack + 82075;
 Libderiv->deriv_classes[4][4][7] = int_stack + 84100;
 Libderiv->deriv_classes[4][5][7] = int_stack + 84325;
 Libderiv->deriv_classes[4][6][7] = int_stack + 84640;
 Libderiv->deriv_classes[4][7][7] = int_stack + 85060;
 Libderiv->deriv_classes[4][8][7] = int_stack + 85600;
 Libderiv->deriv_classes[5][4][7] = int_stack + 86275;
 Libderiv->deriv_classes[5][5][7] = int_stack + 86590;
 Libderiv->deriv_classes[5][6][7] = int_stack + 87031;
 Libderiv->deriv_classes[5][7][7] = int_stack + 87619;
 Libderiv->deriv_classes[5][8][7] = int_stack + 88375;
 Libderiv->deriv_classes[6][4][7] = int_stack + 89320;
 Libderiv->deriv_classes[6][5][7] = int_stack + 89740;
 Libderiv->deriv_classes[6][6][7] = int_stack + 90328;
 Libderiv->deriv_classes[6][7][7] = int_stack + 91112;
 Libderiv->deriv_classes[6][8][7] = int_stack + 92120;
 Libderiv->deriv_classes[7][4][7] = int_stack + 93380;
 Libderiv->deriv_classes[7][5][7] = int_stack + 93920;
 Libderiv->deriv_classes[7][6][7] = int_stack + 94676;
 Libderiv->deriv_classes[7][7][7] = int_stack + 95684;
 Libderiv->deriv_classes[7][8][7] = int_stack + 96980;
 Libderiv->deriv_classes[8][4][7] = int_stack + 98600;
 Libderiv->deriv_classes[8][5][7] = int_stack + 99275;
 Libderiv->deriv_classes[8][6][7] = int_stack + 100220;
 Libderiv->deriv_classes[8][7][7] = int_stack + 101480;
 Libderiv->deriv_classes[8][8][7] = int_stack + 103100;
 Libderiv->deriv_classes[4][4][6] = int_stack + 105125;
 Libderiv->deriv_classes[4][5][6] = int_stack + 105350;
 Libderiv->deriv_classes[4][6][6] = int_stack + 105665;
 Libderiv->deriv_classes[4][7][6] = int_stack + 106085;
 Libderiv->deriv_classes[4][8][6] = int_stack + 106625;
 Libderiv->deriv_classes[5][4][6] = int_stack + 107300;
 Libderiv->deriv_classes[5][5][6] = int_stack + 107615;
 Libderiv->deriv_classes[5][6][6] = int_stack + 108056;
 Libderiv->deriv_classes[5][7][6] = int_stack + 108644;
 Libderiv->deriv_classes[5][8][6] = int_stack + 109400;
 Libderiv->deriv_classes[6][4][6] = int_stack + 110345;
 Libderiv->deriv_classes[6][5][6] = int_stack + 110765;
 Libderiv->deriv_classes[6][6][6] = int_stack + 111353;
 Libderiv->deriv_classes[6][7][6] = int_stack + 112137;
 Libderiv->deriv_classes[6][8][6] = int_stack + 113145;
 Libderiv->deriv_classes[7][4][6] = int_stack + 114405;
 Libderiv->deriv_classes[7][5][6] = int_stack + 114945;
 Libderiv->deriv_classes[7][6][6] = int_stack + 115701;
 Libderiv->deriv_classes[7][7][6] = int_stack + 116709;
 Libderiv->deriv_classes[7][8][6] = int_stack + 118005;
 Libderiv->dvrr_classes[8][4] = int_stack + 119625;
 Libderiv->deriv_classes[8][4][6] = int_stack + 120300;
 Libderiv->dvrr_classes[8][5] = int_stack + 120975;
 Libderiv->deriv_classes[8][5][6] = int_stack + 121920;
 Libderiv->dvrr_classes[8][6] = int_stack + 122865;
 Libderiv->deriv_classes[8][6][6] = int_stack + 124125;
 Libderiv->dvrr_classes[8][7] = int_stack + 125385;
 Libderiv->deriv_classes[8][7][6] = int_stack + 127005;
 Libderiv->deriv_classes[8][8][6] = int_stack + 128625;
 Libderiv->deriv_classes[4][4][2] = int_stack + 130650;
 Libderiv->deriv_classes[4][5][2] = int_stack + 130875;
 Libderiv->deriv_classes[4][6][2] = int_stack + 131190;
 Libderiv->deriv_classes[4][7][2] = int_stack + 131610;
 Libderiv->deriv_classes[4][8][2] = int_stack + 132150;
 Libderiv->deriv_classes[5][4][2] = int_stack + 132825;
 Libderiv->deriv_classes[5][5][2] = int_stack + 133140;
 Libderiv->deriv_classes[5][6][2] = int_stack + 133581;
 Libderiv->deriv_classes[5][7][2] = int_stack + 134169;
 Libderiv->deriv_classes[5][8][2] = int_stack + 134925;
 Libderiv->deriv_classes[6][4][2] = int_stack + 135870;
 Libderiv->deriv_classes[6][5][2] = int_stack + 136290;
 Libderiv->deriv_classes[6][6][2] = int_stack + 136878;
 Libderiv->deriv_classes[6][7][2] = int_stack + 137662;
 Libderiv->deriv_classes[6][8][2] = int_stack + 138670;
 Libderiv->deriv_classes[7][4][2] = int_stack + 139930;
 Libderiv->deriv_classes[7][5][2] = int_stack + 140470;
 Libderiv->deriv_classes[7][6][2] = int_stack + 141226;
 Libderiv->deriv_classes[7][7][2] = int_stack + 142234;
 Libderiv->deriv_classes[7][8][2] = int_stack + 143530;
 Libderiv->deriv_classes[8][4][2] = int_stack + 145150;
 Libderiv->deriv_classes[8][5][2] = int_stack + 145825;
 Libderiv->deriv_classes[8][6][2] = int_stack + 146770;
 Libderiv->deriv_classes[8][7][2] = int_stack + 148030;
 Libderiv->deriv_classes[8][8][2] = int_stack + 149650;
 Libderiv->deriv_classes[4][4][1] = int_stack + 151675;
 Libderiv->deriv_classes[4][5][1] = int_stack + 151900;
 Libderiv->deriv_classes[4][6][1] = int_stack + 152215;
 Libderiv->deriv_classes[4][7][1] = int_stack + 152635;
 Libderiv->deriv_classes[4][8][1] = int_stack + 153175;
 Libderiv->deriv_classes[5][4][1] = int_stack + 153850;
 Libderiv->deriv_classes[5][5][1] = int_stack + 154165;
 Libderiv->deriv_classes[5][6][1] = int_stack + 154606;
 Libderiv->deriv_classes[5][7][1] = int_stack + 155194;
 Libderiv->deriv_classes[5][8][1] = int_stack + 155950;
 Libderiv->deriv_classes[6][4][1] = int_stack + 156895;
 Libderiv->deriv_classes[6][5][1] = int_stack + 157315;
 Libderiv->deriv_classes[6][6][1] = int_stack + 157903;
 Libderiv->deriv_classes[6][7][1] = int_stack + 158687;
 Libderiv->deriv_classes[6][8][1] = int_stack + 159695;
 Libderiv->deriv_classes[7][4][1] = int_stack + 160955;
 Libderiv->deriv_classes[7][5][1] = int_stack + 161495;
 Libderiv->deriv_classes[7][6][1] = int_stack + 162251;
 Libderiv->deriv_classes[7][7][1] = int_stack + 163259;
 Libderiv->deriv_classes[7][8][1] = int_stack + 164555;
 Libderiv->deriv_classes[8][4][1] = int_stack + 166175;
 Libderiv->deriv_classes[8][5][1] = int_stack + 166850;
 Libderiv->deriv_classes[8][6][1] = int_stack + 167795;
 Libderiv->deriv_classes[8][7][1] = int_stack + 169055;
 Libderiv->deriv_classes[8][8][1] = int_stack + 170675;
 Libderiv->dvrr_classes[4][4] = int_stack + 172700;
 Libderiv->dvrr_classes[4][5] = int_stack + 172925;
 Libderiv->dvrr_classes[4][6] = int_stack + 173240;
 Libderiv->dvrr_classes[4][7] = int_stack + 173660;
 Libderiv->dvrr_classes[4][8] = int_stack + 174200;
 Libderiv->deriv_classes[4][4][0] = int_stack + 174875;
 Libderiv->deriv_classes[4][5][0] = int_stack + 175100;
 Libderiv->deriv_classes[4][6][0] = int_stack + 175415;
 Libderiv->deriv_classes[4][7][0] = int_stack + 175835;
 Libderiv->deriv_classes[4][8][0] = int_stack + 176375;
 Libderiv->dvrr_classes[5][4] = int_stack + 177050;
 Libderiv->dvrr_classes[5][5] = int_stack + 177365;
 Libderiv->dvrr_classes[5][6] = int_stack + 177806;
 Libderiv->dvrr_classes[5][7] = int_stack + 178394;
 Libderiv->dvrr_classes[5][8] = int_stack + 179150;
 Libderiv->deriv_classes[5][4][0] = int_stack + 180095;
 Libderiv->deriv_classes[5][5][0] = int_stack + 180410;
 Libderiv->deriv_classes[5][6][0] = int_stack + 180851;
 Libderiv->deriv_classes[5][7][0] = int_stack + 181439;
 Libderiv->deriv_classes[5][8][0] = int_stack + 182195;
 Libderiv->dvrr_classes[6][4] = int_stack + 183140;
 Libderiv->dvrr_classes[6][5] = int_stack + 183560;
 Libderiv->dvrr_classes[6][6] = int_stack + 184148;
 Libderiv->dvrr_classes[6][7] = int_stack + 184932;
 Libderiv->dvrr_classes[6][8] = int_stack + 185940;
 Libderiv->deriv_classes[6][4][0] = int_stack + 187200;
 Libderiv->deriv_classes[6][5][0] = int_stack + 187620;
 Libderiv->deriv_classes[6][6][0] = int_stack + 188208;
 Libderiv->deriv_classes[6][7][0] = int_stack + 188992;
 Libderiv->deriv_classes[6][8][0] = int_stack + 190000;
 Libderiv->dvrr_classes[7][4] = int_stack + 191260;
 Libderiv->dvrr_classes[7][5] = int_stack + 191800;
 Libderiv->dvrr_classes[7][6] = int_stack + 192556;
 Libderiv->dvrr_classes[7][7] = int_stack + 193564;
 Libderiv->dvrr_classes[7][8] = int_stack + 194860;
 Libderiv->deriv_classes[7][4][0] = int_stack + 196480;
 Libderiv->deriv_classes[7][5][0] = int_stack + 197020;
 Libderiv->deriv_classes[7][6][0] = int_stack + 197776;
 Libderiv->deriv_classes[7][7][0] = int_stack + 198784;
 Libderiv->deriv_classes[7][8][0] = int_stack + 200080;
 Libderiv->deriv_classes[8][4][0] = int_stack + 201700;
 Libderiv->deriv_classes[8][5][0] = int_stack + 202375;
 Libderiv->deriv_classes[8][6][0] = int_stack + 203320;
 Libderiv->deriv_classes[8][7][0] = int_stack + 204580;
 Libderiv->deriv_classes[8][8][0] = int_stack + 206200;
 memset(int_stack,0,1665800);

 Libderiv->dvrr_stack = int_stack + 1022680;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gggg(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+208225,int_stack+172925,int_stack+172700,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+208900,int_stack+173240,int_stack+172925,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+209845,int_stack+208900,int_stack+208225,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+211195,int_stack+173660,int_stack+173240,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+212455,int_stack+211195,int_stack+208900,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+214345,int_stack+212455,int_stack+209845,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+216595,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172700,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+217270,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172925,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+218215,int_stack+217270,int_stack+216595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208225,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+219565,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173240,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+220825,int_stack+219565,int_stack+217270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208900,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+222715,int_stack+220825,int_stack+218215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+209845,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+216595,int_stack+1500,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173660,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+224965,int_stack+216595,int_stack+219565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211195,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+216595,int_stack+224965,int_stack+220825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212455,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+224965,int_stack+216595,int_stack+222715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+214345,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+216595,int_stack+177365,int_stack+177050,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+217540,int_stack+177806,int_stack+177365,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+218863,int_stack+217540,int_stack+216595,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+220753,int_stack+178394,int_stack+177806,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+228340,int_stack+220753,int_stack+217540,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+230986,int_stack+228340,int_stack+218863,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+222517,int_stack+2490,int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177050,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+223462,int_stack+2931,int_stack+2490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177365,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+223462,int_stack+222517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216595,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+234136,int_stack+3519,int_stack+2931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177806,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+235900,int_stack+234136,int_stack+223462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217540,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+238546,int_stack+235900,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218863,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+4275,int_stack+3519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178394,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+241696,int_stack+0,int_stack+234136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220753,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+241696,int_stack+235900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+228340,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+241696,int_stack+0,int_stack+238546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230986,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+246421,int_stack+241696,int_stack+224965,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+183560,int_stack+183140,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1260,int_stack+184148,int_stack+183560,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+234136,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+236656,int_stack+184932,int_stack+184148,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+222517,int_stack+236656,int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+256546,int_stack+222517,int_stack+234136,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+5640,int_stack+5220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183140,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239008,int_stack+6228,int_stack+5640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183560,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+260746,int_stack+239008,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3024,int_stack+7012,int_stack+6228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184148,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+263266,int_stack+3024,int_stack+239008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+266794,int_stack+263266,int_stack+260746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234136,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+270994,int_stack+8020,int_stack+7012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184932,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+274018,int_stack+270994,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+236656,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3024,int_stack+274018,int_stack+263266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+222517,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+270994,int_stack+3024,int_stack+266794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+256546,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+277294,int_stack+270994,int_stack+241696,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+291469,int_stack+277294,int_stack+246421,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+191800,int_stack+191260,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4644,int_stack+192556,int_stack+191800,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+260746,int_stack+4644,int_stack+3024,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+263986,int_stack+193564,int_stack+192556,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+239008,int_stack+263986,int_stack+4644,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+243544,int_stack+239008,int_stack+260746,36);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6912,int_stack+9820,int_stack+9280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191260,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+267010,int_stack+10576,int_stack+9820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191800,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+248944,int_stack+267010,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+6912,int_stack+11584,int_stack+10576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192556,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+311719,int_stack+6912,int_stack+267010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+316255,int_stack+311719,int_stack+248944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260746,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+248944,int_stack+12880,int_stack+11584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193564,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+321655,int_stack+248944,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+263986,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+6912,int_stack+321655,int_stack+311719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+239008,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+321655,int_stack+6912,int_stack+316255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243544,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+329755,int_stack+321655,int_stack+270994,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+348655,int_stack+329755,int_stack+277294,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+377005,int_stack+348655,int_stack+291469,225);
 /*--- compute (l0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6912,int_stack+120975,int_stack+119625,45);
 /*--- compute (l0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8937,int_stack+122865,int_stack+120975,45);
 /*--- compute (l0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+248944,int_stack+8937,int_stack+6912,45);
 /*--- compute (l0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+267010,int_stack+125385,int_stack+122865,45);
 /*--- compute (l0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+270790,int_stack+267010,int_stack+8937,45);
 /*--- compute (l0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+276460,int_stack+270790,int_stack+248944,45);
 /*--- compute (l0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11772,int_stack+15175,int_stack+14500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119625,45);
 /*--- compute (l0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+252994,int_stack+16120,int_stack+15175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120975,45);
 /*--- compute (l0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+283210,int_stack+252994,int_stack+11772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6912,45);
 /*--- compute (l0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+11772,int_stack+17380,int_stack+16120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122865,45);
 /*--- compute (l0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+287260,int_stack+11772,int_stack+252994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8937,45);
 /*--- compute (l0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+292930,int_stack+287260,int_stack+283210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+248944,45);
 /*--- compute (l0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+299680,int_stack+19000,int_stack+17380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125385,45);
 /*--- compute (l0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+304540,int_stack+299680,int_stack+11772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+267010,45);
 /*--- compute (l0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312100,int_stack+304540,int_stack+287260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270790,45);
 /*--- compute (l0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+299680,int_stack+312100,int_stack+292930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276460,45);
 /*--- compute (kp|gg) ---*/
   hrr1_build_kp(Libderiv->AB,int_stack+410755,int_stack+299680,int_stack+321655,225);
 /*--- compute (id|gg) ---*/
   hrr1_build_id(Libderiv->AB,int_stack+283210,int_stack+410755,int_stack+329755,225);
 /*--- compute (hf|gg) ---*/
   hrr1_build_hf(Libderiv->AB,int_stack+410755,int_stack+283210,int_stack+348655,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+283210,int_stack+21250,int_stack+21025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172700, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+283885,int_stack+21565,int_stack+21250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172925, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+284830,int_stack+283885,int_stack+283210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208225, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+286180,int_stack+21985,int_stack+21565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173240, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+287440,int_stack+286180,int_stack+283885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208900, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+289330,int_stack+287440,int_stack+284830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+209845, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+283210,int_stack+22525,int_stack+21985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173660, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+291580,int_stack+283210,int_stack+286180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211195, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+283210,int_stack+291580,int_stack+287440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212455, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+291580,int_stack+283210,int_stack+289330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+214345, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+283210,int_stack+23515,int_stack+23200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177050, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+284155,int_stack+23956,int_stack+23515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177365, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+285478,int_stack+284155,int_stack+283210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216595, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+287368,int_stack+24544,int_stack+23956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177806, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+294955,int_stack+287368,int_stack+284155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217540, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+297601,int_stack+294955,int_stack+285478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218863, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+283210,int_stack+25300,int_stack+24544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178394, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+300751,int_stack+283210,int_stack+287368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220753, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+283210,int_stack+300751,int_stack+294955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+228340, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+300751,int_stack+283210,int_stack+297601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230986, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+305476,int_stack+300751,int_stack+291580,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+283210,int_stack+26665,int_stack+26245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183140, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+284470,int_stack+27253,int_stack+26665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183560, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+286234,int_stack+284470,int_stack+283210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+288754,int_stack+28037,int_stack+27253, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184148, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+291106,int_stack+288754,int_stack+284470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+294634,int_stack+291106,int_stack+286234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234136, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+283210,int_stack+29045,int_stack+28037, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184932, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+315601,int_stack+283210,int_stack+288754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+236656, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+283210,int_stack+315601,int_stack+291106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+222517, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+315601,int_stack+283210,int_stack+294634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+256546, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+283210,int_stack+315601,int_stack+300751,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+321901,int_stack+283210,int_stack+305476,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+297385,int_stack+30845,int_stack+30305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191260, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+299005,int_stack+31601,int_stack+30845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191800, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+301273,int_stack+299005,int_stack+297385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+304513,int_stack+32609,int_stack+31601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192556, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+307537,int_stack+304513,int_stack+299005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+342151,int_stack+307537,int_stack+301273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260746, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+297385,int_stack+33905,int_stack+32609, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193564, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+347551,int_stack+297385,int_stack+304513, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+263986, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+297385,int_stack+347551,int_stack+307537, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+239008, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+347551,int_stack+297385,int_stack+342151, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243544, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+355651,int_stack+347551,int_stack+315601,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+458005,int_stack+355651,int_stack+283210,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+283210,int_stack+458005,int_stack+321901,225);
 /*--- compute (l0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+316960,int_stack+36200,int_stack+35525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119625, 0.0, zero_stack,45);
 /*--- compute (l0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+318985,int_stack+37145,int_stack+36200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120975, 0.0, zero_stack,45);
 /*--- compute (l0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+321820,int_stack+318985,int_stack+316960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6912, 0.0, zero_stack,45);
 /*--- compute (l0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+325870,int_stack+38405,int_stack+37145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122865, 0.0, zero_stack,45);
 /*--- compute (l0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+329650,int_stack+325870,int_stack+318985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8937, 0.0, zero_stack,45);
 /*--- compute (l0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+335320,int_stack+329650,int_stack+321820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+248944, 0.0, zero_stack,45);
 /*--- compute (l0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+316960,int_stack+40025,int_stack+38405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125385, 0.0, zero_stack,45);
 /*--- compute (l0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+11772,int_stack+316960,int_stack+325870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+267010, 0.0, zero_stack,45);
 /*--- compute (l0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+316960,int_stack+11772,int_stack+329650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270790, 0.0, zero_stack,45);
 /*--- compute (l0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+11772,int_stack+316960,int_stack+335320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276460, 0.0, zero_stack,45);
 /*--- compute (kp|gg) ---*/
   hrr1_build_kp(Libderiv->AB,int_stack+316960,int_stack+11772,int_stack+347551,225);
 /*--- compute (id|gg) ---*/
   hrr1_build_id(Libderiv->AB,int_stack+486355,int_stack+316960,int_stack+355651,225);
 /*--- compute (hf|gg) ---*/
   hrr1_build_hf(Libderiv->AB,int_stack+316960,int_stack+486355,int_stack+458005,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+458005,int_stack+42275,int_stack+42050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172700, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+458680,int_stack+42590,int_stack+42275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172925, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+459625,int_stack+458680,int_stack+458005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208225, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+460975,int_stack+43010,int_stack+42590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173240, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+462235,int_stack+460975,int_stack+458680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208900, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+464125,int_stack+462235,int_stack+459625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+209845, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+458005,int_stack+43550,int_stack+43010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173660, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+466375,int_stack+458005,int_stack+460975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211195, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+458005,int_stack+466375,int_stack+462235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212455, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+466375,int_stack+458005,int_stack+464125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+214345, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+458005,int_stack+44540,int_stack+44225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177050, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+458950,int_stack+44981,int_stack+44540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177365, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+460273,int_stack+458950,int_stack+458005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216595, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+462163,int_stack+45569,int_stack+44981, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177806, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+469750,int_stack+462163,int_stack+458950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217540, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+472396,int_stack+469750,int_stack+460273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218863, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+458005,int_stack+46325,int_stack+45569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178394, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+475546,int_stack+458005,int_stack+462163, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220753, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+458005,int_stack+475546,int_stack+469750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+228340, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+475546,int_stack+458005,int_stack+472396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230986, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+480271,int_stack+475546,int_stack+466375,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+458005,int_stack+47690,int_stack+47270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183140, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+459265,int_stack+48278,int_stack+47690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183560, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+461029,int_stack+459265,int_stack+458005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+463549,int_stack+49062,int_stack+48278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184148, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+465901,int_stack+463549,int_stack+459265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+469429,int_stack+465901,int_stack+461029, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234136, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+458005,int_stack+50070,int_stack+49062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184932, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+490396,int_stack+458005,int_stack+463549, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+236656, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+458005,int_stack+490396,int_stack+465901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+222517, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+490396,int_stack+458005,int_stack+469429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+256546, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+458005,int_stack+490396,int_stack+475546,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+496696,int_stack+458005,int_stack+480271,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+472180,int_stack+51870,int_stack+51330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191260, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+473800,int_stack+52626,int_stack+51870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191800, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+476068,int_stack+473800,int_stack+472180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+479308,int_stack+53634,int_stack+52626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192556, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+482332,int_stack+479308,int_stack+473800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+516946,int_stack+482332,int_stack+476068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260746, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+472180,int_stack+54930,int_stack+53634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193564, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+364210,int_stack+472180,int_stack+479308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+263986, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+472180,int_stack+364210,int_stack+482332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+239008, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+364210,int_stack+472180,int_stack+516946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243544, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+11772,int_stack+364210,int_stack+490396,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+516946,int_stack+11772,int_stack+458005,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+458005,int_stack+516946,int_stack+496696,225);
 /*--- compute (l0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+491755,int_stack+57225,int_stack+56550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119625, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+493780,int_stack+58170,int_stack+57225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120975, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+496615,int_stack+493780,int_stack+491755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6912, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+500665,int_stack+59430,int_stack+58170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122865, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+504445,int_stack+500665,int_stack+493780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8937, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+510115,int_stack+504445,int_stack+496615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+248944, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+491755,int_stack+61050,int_stack+59430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125385, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+30672,int_stack+491755,int_stack+500665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+267010, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+491755,int_stack+30672,int_stack+504445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270790, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+30672,int_stack+491755,int_stack+510115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276460, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (kp|gg) ---*/
   hrr1_build_kp(Libderiv->AB,int_stack+491755,int_stack+30672,int_stack+364210,225);
 /*--- compute (id|gg) ---*/
   hrr1_build_id(Libderiv->AB,int_stack+545296,int_stack+491755,int_stack+11772,225);
 /*--- compute (hf|gg) ---*/
   hrr1_build_hf(Libderiv->AB,int_stack+11772,int_stack+545296,int_stack+516946,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59022,int_stack+63300,int_stack+63075, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+59697,int_stack+63615,int_stack+63300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+60642,int_stack+59697,int_stack+59022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+61992,int_stack+64035,int_stack+63615, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+491755,int_stack+61992,int_stack+59697, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+493645,int_stack+491755,int_stack+60642, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+209845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+59022,int_stack+64575,int_stack+64035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+173660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+495895,int_stack+59022,int_stack+61992, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+59022,int_stack+495895,int_stack+491755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+495895,int_stack+59022,int_stack+493645, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+214345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59022,int_stack+65565,int_stack+65250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+59967,int_stack+66006,int_stack+65565, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+61290,int_stack+59967,int_stack+59022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+63180,int_stack+66594,int_stack+66006, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+491755,int_stack+63180,int_stack+59967, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+499270,int_stack+491755,int_stack+61290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+59022,int_stack+67350,int_stack+66594, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+502420,int_stack+59022,int_stack+63180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220753, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+59022,int_stack+502420,int_stack+491755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+228340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+502420,int_stack+59022,int_stack+499270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+230986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+507145,int_stack+502420,int_stack+495895,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59022,int_stack+68715,int_stack+68295, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60282,int_stack+69303,int_stack+68715, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+62046,int_stack+60282,int_stack+59022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+64566,int_stack+70087,int_stack+69303, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+491755,int_stack+64566,int_stack+60282, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+495283,int_stack+491755,int_stack+62046, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+59022,int_stack+71095,int_stack+70087, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+184932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+66918,int_stack+59022,int_stack+64566, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+236656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+59022,int_stack+66918,int_stack+491755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+222517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+64902,int_stack+59022,int_stack+495283, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+256546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+517270,int_stack+64902,int_stack+502420,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+531445,int_stack+517270,int_stack+507145,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59022,int_stack+72895,int_stack+72355, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60642,int_stack+73651,int_stack+72895, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+491755,int_stack+60642,int_stack+59022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+494995,int_stack+74659,int_stack+73651, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+498019,int_stack+494995,int_stack+60642, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+59022,int_stack+498019,int_stack+491755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+502555,int_stack+75955,int_stack+74659, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+506443,int_stack+502555,int_stack+494995, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+263986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+551695,int_stack+506443,int_stack+498019, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+239008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+491755,int_stack+551695,int_stack+59022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+551695,int_stack+491755,int_stack+64902,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+570595,int_stack+551695,int_stack+517270,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+598945,int_stack+570595,int_stack+531445,225);
 /*--- compute (l0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59022,int_stack+78250,int_stack+77575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+61047,int_stack+79195,int_stack+78250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+63882,int_stack+61047,int_stack+59022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+67932,int_stack+80455,int_stack+79195, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+71712,int_stack+67932,int_stack+61047, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8937, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+499855,int_stack+71712,int_stack+63882, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+248944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+59022,int_stack+82075,int_stack+80455, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+506605,int_stack+59022,int_stack+67932, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+267010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+59022,int_stack+506605,int_stack+71712, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+506605,int_stack+59022,int_stack+499855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (kp|gg) ---*/
   hrr1_build_kp(Libderiv->AB,int_stack+59022,int_stack+506605,int_stack+491755,225);
 /*--- compute (id|gg) ---*/
   hrr1_build_id(Libderiv->AB,int_stack+491755,int_stack+59022,int_stack+551695,225);
 /*--- compute (hf|gg) ---*/
   hrr1_build_hf(Libderiv->AB,int_stack+632695,int_stack+491755,int_stack+570595,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+491755,int_stack+84325,int_stack+84100, 0.0, zero_stack, 1.0, int_stack+172700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+492430,int_stack+84640,int_stack+84325, 0.0, zero_stack, 1.0, int_stack+172925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+493375,int_stack+492430,int_stack+491755, 0.0, zero_stack, 1.0, int_stack+208225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+494725,int_stack+85060,int_stack+84640, 0.0, zero_stack, 1.0, int_stack+173240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+495985,int_stack+494725,int_stack+492430, 0.0, zero_stack, 1.0, int_stack+208900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+497875,int_stack+495985,int_stack+493375, 0.0, zero_stack, 1.0, int_stack+209845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+491755,int_stack+85600,int_stack+85060, 0.0, zero_stack, 1.0, int_stack+173660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+500125,int_stack+491755,int_stack+494725, 0.0, zero_stack, 1.0, int_stack+211195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+491755,int_stack+500125,int_stack+495985, 0.0, zero_stack, 1.0, int_stack+212455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+500125,int_stack+491755,int_stack+497875, 0.0, zero_stack, 1.0, int_stack+214345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+491755,int_stack+86590,int_stack+86275, 0.0, zero_stack, 1.0, int_stack+177050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+492700,int_stack+87031,int_stack+86590, 0.0, zero_stack, 1.0, int_stack+177365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+494023,int_stack+492700,int_stack+491755, 0.0, zero_stack, 1.0, int_stack+216595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+495913,int_stack+87619,int_stack+87031, 0.0, zero_stack, 1.0, int_stack+177806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+503500,int_stack+495913,int_stack+492700, 0.0, zero_stack, 1.0, int_stack+217540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+506146,int_stack+503500,int_stack+494023, 0.0, zero_stack, 1.0, int_stack+218863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+491755,int_stack+88375,int_stack+87619, 0.0, zero_stack, 1.0, int_stack+178394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+509296,int_stack+491755,int_stack+495913, 0.0, zero_stack, 1.0, int_stack+220753, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+491755,int_stack+509296,int_stack+503500, 0.0, zero_stack, 1.0, int_stack+228340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+509296,int_stack+491755,int_stack+506146, 0.0, zero_stack, 1.0, int_stack+230986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+514021,int_stack+509296,int_stack+500125,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+491755,int_stack+89740,int_stack+89320, 0.0, zero_stack, 1.0, int_stack+183140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+493015,int_stack+90328,int_stack+89740, 0.0, zero_stack, 1.0, int_stack+183560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+494779,int_stack+493015,int_stack+491755, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+497299,int_stack+91112,int_stack+90328, 0.0, zero_stack, 1.0, int_stack+184148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+499651,int_stack+497299,int_stack+493015, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+503179,int_stack+499651,int_stack+494779, 0.0, zero_stack, 1.0, int_stack+234136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+491755,int_stack+92120,int_stack+91112, 0.0, zero_stack, 1.0, int_stack+184932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+524146,int_stack+491755,int_stack+497299, 0.0, zero_stack, 1.0, int_stack+236656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+491755,int_stack+524146,int_stack+499651, 0.0, zero_stack, 1.0, int_stack+222517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+524146,int_stack+491755,int_stack+503179, 0.0, zero_stack, 1.0, int_stack+256546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+491755,int_stack+524146,int_stack+509296,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+530446,int_stack+491755,int_stack+514021,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+505930,int_stack+93920,int_stack+93380, 0.0, zero_stack, 1.0, int_stack+191260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+507550,int_stack+94676,int_stack+93920, 0.0, zero_stack, 1.0, int_stack+191800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+509818,int_stack+507550,int_stack+505930, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+513058,int_stack+95684,int_stack+94676, 0.0, zero_stack, 1.0, int_stack+192556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+516082,int_stack+513058,int_stack+507550, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+550696,int_stack+516082,int_stack+509818, 0.0, zero_stack, 1.0, int_stack+260746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+505930,int_stack+96980,int_stack+95684, 0.0, zero_stack, 1.0, int_stack+193564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+556096,int_stack+505930,int_stack+513058, 0.0, zero_stack, 1.0, int_stack+263986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+505930,int_stack+556096,int_stack+516082, 0.0, zero_stack, 1.0, int_stack+239008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+556096,int_stack+505930,int_stack+550696, 0.0, zero_stack, 1.0, int_stack+243544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+564196,int_stack+556096,int_stack+524146,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+59022,int_stack+564196,int_stack+491755,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+491755,int_stack+59022,int_stack+530446,225);
 /*--- compute (l0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+525505,int_stack+99275,int_stack+98600, 0.0, zero_stack, 1.0, int_stack+119625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+527530,int_stack+100220,int_stack+99275, 0.0, zero_stack, 1.0, int_stack+120975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+530365,int_stack+527530,int_stack+525505, 0.0, zero_stack, 1.0, int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+534415,int_stack+101480,int_stack+100220, 0.0, zero_stack, 1.0, int_stack+122865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+538195,int_stack+534415,int_stack+527530, 0.0, zero_stack, 1.0, int_stack+8937, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+543865,int_stack+538195,int_stack+530365, 0.0, zero_stack, 1.0, int_stack+248944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+525505,int_stack+103100,int_stack+101480, 0.0, zero_stack, 1.0, int_stack+125385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+583096,int_stack+525505,int_stack+534415, 0.0, zero_stack, 1.0, int_stack+267010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+525505,int_stack+583096,int_stack+538195, 0.0, zero_stack, 1.0, int_stack+270790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+583096,int_stack+525505,int_stack+543865, 0.0, zero_stack, 1.0, int_stack+276460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (kp|gg) ---*/
   hrr1_build_kp(Libderiv->AB,int_stack+525505,int_stack+583096,int_stack+556096,225);
 /*--- compute (id|gg) ---*/
   hrr1_build_id(Libderiv->AB,int_stack+679945,int_stack+525505,int_stack+564196,225);
 /*--- compute (hf|gg) ---*/
   hrr1_build_hf(Libderiv->AB,int_stack+525505,int_stack+679945,int_stack+59022,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59022,int_stack+105350,int_stack+105125, 1.0, int_stack+172700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+59697,int_stack+105665,int_stack+105350, 1.0, int_stack+172925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+60642,int_stack+59697,int_stack+59022, 1.0, int_stack+208225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+61992,int_stack+106085,int_stack+105665, 1.0, int_stack+173240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+63252,int_stack+61992,int_stack+59697, 1.0, int_stack+208900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+65142,int_stack+63252,int_stack+60642, 1.0, int_stack+209845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+208225,int_stack+106625,int_stack+106085, 1.0, int_stack+173660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+59022,int_stack+208225,int_stack+61992, 1.0, int_stack+211195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+67392,int_stack+59022,int_stack+63252, 1.0, int_stack+212455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+59022,int_stack+67392,int_stack+65142, 1.0, int_stack+214345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+62397,int_stack+107615,int_stack+107300, 1.0, int_stack+177050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+63342,int_stack+108056,int_stack+107615, 1.0, int_stack+177365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+64665,int_stack+63342,int_stack+62397, 1.0, int_stack+216595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+66555,int_stack+108644,int_stack+108056, 1.0, int_stack+177806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+68319,int_stack+66555,int_stack+63342, 1.0, int_stack+217540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+70965,int_stack+68319,int_stack+64665, 1.0, int_stack+218863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+216595,int_stack+109400,int_stack+108644, 1.0, int_stack+178394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+62397,int_stack+216595,int_stack+66555, 1.0, int_stack+220753, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+74115,int_stack+62397,int_stack+68319, 1.0, int_stack+228340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+62397,int_stack+74115,int_stack+70965, 1.0, int_stack+230986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+67122,int_stack+62397,int_stack+59022,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59022,int_stack+110765,int_stack+110345, 1.0, int_stack+183140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60282,int_stack+111353,int_stack+110765, 1.0, int_stack+183560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+77247,int_stack+60282,int_stack+59022, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+79767,int_stack+112137,int_stack+111353, 1.0, int_stack+184148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+82119,int_stack+79767,int_stack+60282, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+85647,int_stack+82119,int_stack+77247, 1.0, int_stack+234136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+113145,int_stack+112137, 1.0, int_stack+184932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+89847,int_stack+0,int_stack+79767, 1.0, int_stack+236656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+94551,int_stack+89847,int_stack+82119, 1.0, int_stack+222517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+77247,int_stack+94551,int_stack+85647, 1.0, int_stack+256546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+83547,int_stack+77247,int_stack+62397,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+679945,int_stack+83547,int_stack+67122,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+114945,int_stack+114405, 1.0, int_stack+191260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234136,int_stack+115701,int_stack+114945, 1.0, int_stack+191800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+97722,int_stack+234136,int_stack+0, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+116709,int_stack+115701, 1.0, int_stack+192556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+100962,int_stack+0,int_stack+234136, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+105498,int_stack+100962,int_stack+97722, 1.0, int_stack+260746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+3024,int_stack+118005,int_stack+116709, 1.0, int_stack+193564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+110898,int_stack+3024,int_stack+0, 1.0, int_stack+263986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+59022,int_stack+110898,int_stack+100962, 1.0, int_stack+239008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+110898,int_stack+59022,int_stack+105498, 1.0, int_stack+243544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+572755,int_stack+110898,int_stack+77247,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+700195,int_stack+572755,int_stack+83547,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+59022,int_stack+700195,int_stack+679945,225);
 /*--- compute (l0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+679945,int_stack+121920,int_stack+120300, 1.0, int_stack+119625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+681970,int_stack+124125,int_stack+121920, 1.0, int_stack+120975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+684805,int_stack+681970,int_stack+679945, 1.0, int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+688855,int_stack+127005,int_stack+124125, 1.0, int_stack+122865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+692635,int_stack+688855,int_stack+681970, 1.0, int_stack+8937, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+92772,int_stack+692635,int_stack+684805, 1.0, int_stack+248944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+248944,int_stack+128625,int_stack+127005, 1.0, int_stack+125385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+679945,int_stack+248944,int_stack+688855, 1.0, int_stack+267010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+99522,int_stack+679945,int_stack+692635, 1.0, int_stack+270790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (l0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+679945,int_stack+99522,int_stack+92772, 1.0, int_stack+276460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (kp|gg) ---*/
   hrr1_build_kp(Libderiv->AB,int_stack+728545,int_stack+679945,int_stack+110898,225);
 /*--- compute (id|gg) ---*/
   hrr1_build_id(Libderiv->AB,int_stack+92772,int_stack+728545,int_stack+572755,225);
 /*--- compute (hf|gg) ---*/
   hrr1_build_hf(Libderiv->AB,int_stack+728545,int_stack+92772,int_stack+700195,225);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+92772,int_stack+174200,int_stack+173660,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+234136,int_stack+92772,int_stack+211195,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+92772,int_stack+234136,int_stack+212455,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+95922,int_stack+92772,int_stack+214345,15);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+92772,int_stack+179150,int_stack+178394,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+99297,int_stack+92772,int_stack+220753,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+102825,int_stack+99297,int_stack+228340,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+107235,int_stack+102825,int_stack+230986,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+111960,int_stack+107235,int_stack+95922,225);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+99297,int_stack+185940,int_stack+184932,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+102321,int_stack+99297,int_stack+236656,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+122085,int_stack+102321,int_stack+222517,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+99297,int_stack+122085,int_stack+256546,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+572755,int_stack+99297,int_stack+107235,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+679945,int_stack+572755,int_stack+111960,225);
 /*--- compute (k0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+122085,int_stack+194860,int_stack+193564,36);
 /*--- compute (k0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+586930,int_stack+122085,int_stack+263986,36);
 /*--- compute (k0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+122085,int_stack+586930,int_stack+239008,36);
 /*--- compute (k0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+586930,int_stack+122085,int_stack+243544,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+700195,int_stack+586930,int_stack+99297,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+208225,int_stack+700195,int_stack+572755,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+236575,int_stack+208225,int_stack+679945,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+122085,int_stack+130875,int_stack+130650,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+122760,int_stack+131190,int_stack+130875,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+123705,int_stack+122760,int_stack+122085,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+125055,int_stack+131610,int_stack+131190,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+126315,int_stack+125055,int_stack+122760,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+128205,int_stack+126315,int_stack+123705,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+122085,int_stack+132150,int_stack+131610,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+595030,int_stack+122085,int_stack+125055,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+92772,int_stack+595030,int_stack+126315,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+595030,int_stack+92772,int_stack+128205,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+92772,int_stack+133140,int_stack+132825,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+93717,int_stack+133581,int_stack+133140,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+122085,int_stack+93717,int_stack+92772,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+123975,int_stack+134169,int_stack+133581,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+125739,int_stack+123975,int_stack+93717,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+92772,int_stack+125739,int_stack+122085,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+128385,int_stack+134925,int_stack+134169,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+130653,int_stack+128385,int_stack+123975,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+719095,int_stack+130653,int_stack+125739,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+122085,int_stack+719095,int_stack+92772,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+122085,int_stack+595030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+595030,int_stack+136290,int_stack+135870,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+596290,int_stack+136878,int_stack+136290,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+92772,int_stack+596290,int_stack+595030,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+719095,int_stack+137662,int_stack+136878,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+721447,int_stack+719095,int_stack+596290,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+126810,int_stack+721447,int_stack+92772,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+92772,int_stack+138670,int_stack+137662,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+131010,int_stack+92772,int_stack+719095,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+270325,int_stack+131010,int_stack+721447,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+131010,int_stack+270325,int_stack+126810,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+775795,int_stack+131010,int_stack+122085, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+789970,int_stack+775795,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+140470,int_stack+139930,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1620,int_stack+141226,int_stack+140470,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3888,int_stack+1620,int_stack+0,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+7128,int_stack+142234,int_stack+141226,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+122085,int_stack+7128,int_stack+1620,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+270325,int_stack+122085,int_stack+3888,36);
 /*--- compute (k0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+143530,int_stack+142234,36);
 /*--- compute (k0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+275725,int_stack+0,int_stack+7128,36);
 /*--- compute (k0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+275725,int_stack+122085,36);
 /*--- compute (k0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+122085,int_stack+0,int_stack+270325,36);
 /*--- compute (ip|gg) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+810220,int_stack+122085,int_stack+131010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99297, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hd|gg) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+829120,int_stack+810220,int_stack+775795, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+572755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gf|gg) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+857470,int_stack+829120,int_stack+789970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+679945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (l0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+775795,int_stack+145825,int_stack+145150,45);
 /*--- compute (l0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+777820,int_stack+146770,int_stack+145825,45);
 /*--- compute (l0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+780655,int_stack+777820,int_stack+775795,45);
 /*--- compute (l0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+784705,int_stack+148030,int_stack+146770,45);
 /*--- compute (l0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+788485,int_stack+784705,int_stack+777820,45);
 /*--- compute (l0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+794155,int_stack+788485,int_stack+780655,45);
 /*--- compute (l0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+775795,int_stack+149650,int_stack+148030,45);
 /*--- compute (l0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+800905,int_stack+775795,int_stack+784705,45);
 /*--- compute (l0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+719095,int_stack+800905,int_stack+788485,45);
 /*--- compute (l0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+775795,int_stack+719095,int_stack+794155,45);
 /*--- compute (kp|gg) ---*/
   d1hrr1_build_kp(Libderiv->AB,int_stack+785920,int_stack+775795,int_stack+122085, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+586930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (id|gg) ---*/
   d1hrr1_build_id(Libderiv->AB,int_stack+891220,int_stack+785920,int_stack+810220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+700195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hf|gg) ---*/
   d1hrr1_build_hf(Libderiv->AB,int_stack+775795,int_stack+891220,int_stack+829120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+891220,int_stack+151900,int_stack+151675,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+891895,int_stack+152215,int_stack+151900,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+892840,int_stack+891895,int_stack+891220,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+894190,int_stack+152635,int_stack+152215,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+895450,int_stack+894190,int_stack+891895,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+897340,int_stack+895450,int_stack+892840,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+891220,int_stack+153175,int_stack+152635,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+899590,int_stack+891220,int_stack+894190,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+92772,int_stack+899590,int_stack+895450,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+899590,int_stack+92772,int_stack+897340,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+92772,int_stack+154165,int_stack+153850,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+93717,int_stack+154606,int_stack+154165,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+902965,int_stack+93717,int_stack+92772,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+904855,int_stack+155194,int_stack+154606,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+906619,int_stack+904855,int_stack+93717,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+92772,int_stack+906619,int_stack+902965,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+909265,int_stack+155950,int_stack+155194,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+911533,int_stack+909265,int_stack+904855,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+915061,int_stack+911533,int_stack+906619,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+902965,int_stack+915061,int_stack+92772,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+907690,int_stack+902965,int_stack+899590, 0.0, zero_stack, 1.0, int_stack+95922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+92772,int_stack+157315,int_stack+156895,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+94032,int_stack+157903,int_stack+157315,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+917815,int_stack+94032,int_stack+92772,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+920335,int_stack+158687,int_stack+157903,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+922687,int_stack+920335,int_stack+94032,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+891220,int_stack+922687,int_stack+917815,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+92772,int_stack+159695,int_stack+158687,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+895420,int_stack+92772,int_stack+920335,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+122085,int_stack+895420,int_stack+922687,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+895420,int_stack+122085,int_stack+891220,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+122085,int_stack+895420,int_stack+902965, 0.0, zero_stack, 1.0, int_stack+107235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+136260,int_stack+122085,int_stack+907690, 0.0, zero_stack, 1.0, int_stack+111960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+891220,int_stack+161495,int_stack+160955,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+892840,int_stack+162251,int_stack+161495,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+156510,int_stack+892840,int_stack+891220,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+901720,int_stack+163259,int_stack+162251,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+904744,int_stack+901720,int_stack+892840,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+909280,int_stack+904744,int_stack+156510,36);
 /*--- compute (k0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+156510,int_stack+164555,int_stack+163259,36);
 /*--- compute (k0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+914680,int_stack+156510,int_stack+901720,36);
 /*--- compute (k0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+156510,int_stack+914680,int_stack+904744,36);
 /*--- compute (k0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+914680,int_stack+156510,int_stack+909280,36);
 /*--- compute (ip|gg) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+823045,int_stack+914680,int_stack+895420, 0.0, zero_stack, 1.0, int_stack+99297, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hd|gg) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+922780,int_stack+823045,int_stack+122085, 0.0, zero_stack, 1.0, int_stack+572755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gf|gg) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+951130,int_stack+922780,int_stack+136260, 0.0, zero_stack, 1.0, int_stack+679945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (l0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+122085,int_stack+166850,int_stack+166175,45);
 /*--- compute (l0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+124110,int_stack+167795,int_stack+166850,45);
 /*--- compute (l0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+126945,int_stack+124110,int_stack+122085,45);
 /*--- compute (l0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+130995,int_stack+169055,int_stack+167795,45);
 /*--- compute (l0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+134775,int_stack+130995,int_stack+124110,45);
 /*--- compute (l0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+140445,int_stack+134775,int_stack+126945,45);
 /*--- compute (l0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+122085,int_stack+170675,int_stack+169055,45);
 /*--- compute (l0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+147195,int_stack+122085,int_stack+130995,45);
 /*--- compute (l0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+719095,int_stack+147195,int_stack+134775,45);
 /*--- compute (l0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+147195,int_stack+719095,int_stack+140445,45);
 /*--- compute (kp|gg) ---*/
   d1hrr1_build_kp(Libderiv->AB,int_stack+122085,int_stack+147195,int_stack+914680, 0.0, zero_stack, 1.0, int_stack+586930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (id|gg) ---*/
   d1hrr1_build_id(Libderiv->AB,int_stack+984880,int_stack+122085,int_stack+823045, 0.0, zero_stack, 1.0, int_stack+700195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hf|gg) ---*/
   d1hrr1_build_hf(Libderiv->AB,int_stack+122085,int_stack+984880,int_stack+922780, 0.0, zero_stack, 1.0, int_stack+208225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+984880,int_stack+175100,int_stack+174875,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+985555,int_stack+175415,int_stack+175100,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+986500,int_stack+985555,int_stack+984880,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+987850,int_stack+175835,int_stack+175415,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+989110,int_stack+987850,int_stack+985555,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+991000,int_stack+989110,int_stack+986500,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+984880,int_stack+176375,int_stack+175835,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+993250,int_stack+984880,int_stack+987850,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+92772,int_stack+993250,int_stack+989110,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+993250,int_stack+92772,int_stack+991000,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+92772,int_stack+180410,int_stack+180095,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+93717,int_stack+180851,int_stack+180410,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+996625,int_stack+93717,int_stack+92772,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+998515,int_stack+181439,int_stack+180851,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+1000279,int_stack+998515,int_stack+93717,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+92772,int_stack+1000279,int_stack+996625,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+1002925,int_stack+182195,int_stack+181439,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+1005193,int_stack+1002925,int_stack+998515,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+1008721,int_stack+1005193,int_stack+1000279,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+996625,int_stack+1008721,int_stack+92772,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+1001350,int_stack+996625,int_stack+993250, 1.0, int_stack+95922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+92772,int_stack+187620,int_stack+187200,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+94032,int_stack+188208,int_stack+187620,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+95796,int_stack+94032,int_stack+92772,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+1011475,int_stack+188992,int_stack+188208,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+1013827,int_stack+1011475,int_stack+94032,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+1017355,int_stack+1013827,int_stack+95796,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+92772,int_stack+190000,int_stack+188992,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+984880,int_stack+92772,int_stack+1011475,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+92772,int_stack+984880,int_stack+1013827,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+984880,int_stack+92772,int_stack+1017355,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+823045,int_stack+984880,int_stack+996625, 1.0, int_stack+107235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+837220,int_stack+823045,int_stack+1001350, 1.0, int_stack+111960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+92772,int_stack+197020,int_stack+196480,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+94392,int_stack+197776,int_stack+197020,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+991180,int_stack+94392,int_stack+92772,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+994420,int_stack+198784,int_stack+197776,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+997444,int_stack+994420,int_stack+94392,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+92772,int_stack+997444,int_stack+991180,36);
 /*--- compute (k0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+1001980,int_stack+200080,int_stack+198784,36);
 /*--- compute (k0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+1005868,int_stack+1001980,int_stack+994420,36);
 /*--- compute (k0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+1011916,int_stack+1005868,int_stack+997444,36);
 /*--- compute (k0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+991180,int_stack+1011916,int_stack+92772,36);
 /*--- compute (ip|gg) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+999280,int_stack+991180,int_stack+984880, 1.0, int_stack+99297, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hd|gg) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+92772,int_stack+999280,int_stack+823045, 1.0, int_stack+572755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gf|gg) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+891220,int_stack+92772,int_stack+837220, 1.0, int_stack+679945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (l0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+679945,int_stack+202375,int_stack+201700,45);
 /*--- compute (l0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+681970,int_stack+203320,int_stack+202375,45);
 /*--- compute (l0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+684805,int_stack+681970,int_stack+679945,45);
 /*--- compute (l0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+688855,int_stack+204580,int_stack+203320,45);
 /*--- compute (l0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+692635,int_stack+688855,int_stack+681970,45);
 /*--- compute (l0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+572755,int_stack+692635,int_stack+684805,45);
 /*--- compute (l0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+679945,int_stack+206200,int_stack+204580,45);
 /*--- compute (l0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+823045,int_stack+679945,int_stack+688855,45);
 /*--- compute (l0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+719095,int_stack+823045,int_stack+692635,45);
 /*--- compute (l0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+823045,int_stack+719095,int_stack+572755,45);
 /*--- compute (kp|gg) ---*/
   d1hrr1_build_kp(Libderiv->AB,int_stack+833170,int_stack+823045,int_stack+991180, 1.0, int_stack+586930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (id|gg) ---*/
   d1hrr1_build_id(Libderiv->AB,int_stack+169335,int_stack+833170,int_stack+999280, 1.0, int_stack+700195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hf|gg) ---*/
   d1hrr1_build_hf(Libderiv->AB,int_stack+679945,int_stack+169335,int_stack+92772, 1.0, int_stack+208225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gg|gg) ---*/
   hrr1_build_gg(Libderiv->AB,int_stack+169335,int_stack+410755,int_stack+377005,225);
     Libderiv->ABCD[11] = int_stack + 169335;
 /*--- compute (gg|gg) ---*/
   hrr1_build_gg(Libderiv->AB,int_stack+364210,int_stack+316960,int_stack+283210,225);
     Libderiv->ABCD[10] = int_stack + 364210;
 /*--- compute (gg|gg) ---*/
   hrr1_build_gg(Libderiv->AB,int_stack+270325,int_stack+11772,int_stack+458005,225);
     Libderiv->ABCD[9] = int_stack + 270325;
 /*--- compute (gg|gg) ---*/
   hrr1_build_gg(Libderiv->AB,int_stack+0,int_stack+632695,int_stack+598945,225);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (gg|gg) ---*/
   hrr1_build_gg(Libderiv->AB,int_stack+572755,int_stack+525505,int_stack+491755,225);
     Libderiv->ABCD[7] = int_stack + 572755;
 /*--- compute (gg|gg) ---*/
   hrr1_build_gg(Libderiv->AB,int_stack+623380,int_stack+728545,int_stack+59022,225);
     Libderiv->ABCD[6] = int_stack + 623380;
 /*--- compute (gg|gg) ---*/
   d1hrr1_build_gg(Libderiv->AB,int_stack+50625,int_stack+775795,int_stack+857470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+236575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 50625;
 /*--- compute (gg|gg) ---*/
   d1hrr1_build_gg(Libderiv->AB,int_stack+727195,int_stack+122085,int_stack+951130, 0.0, zero_stack, 1.0, int_stack+236575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 727195;
 /*--- compute (gg|gg) ---*/
   d1hrr1_build_gg(Libderiv->AB,int_stack+777820,int_stack+679945,int_stack+891220, 1.0, int_stack+236575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 777820;

}
