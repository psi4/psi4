#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_fdff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|ff) integrals */

void d12hrr_order_fdff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[5][6][11] = int_stack + 0;
 Libderiv->deriv_classes[5][6][10] = int_stack + 588;
 Libderiv->deriv_classes[5][6][9] = int_stack + 1176;
 Libderiv->deriv_classes[5][6][8] = int_stack + 1764;
 Libderiv->deriv_classes[5][6][7] = int_stack + 2352;
 Libderiv->dvrr_classes[5][5] = int_stack + 2940;
 Libderiv->deriv_classes[5][6][6] = int_stack + 3381;
 Libderiv->deriv_classes[5][6][2] = int_stack + 3969;
 Libderiv->deriv_classes[5][6][1] = int_stack + 4557;
 Libderiv->dvrr_classes[4][6] = int_stack + 5145;
 Libderiv->deriv_classes[5][6][0] = int_stack + 5565;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 6153;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 6253;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 6403;
 Libderiv->deriv2_classes[3][6][143] = int_stack + 6613;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 6893;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 7043;
 Libderiv->deriv2_classes[4][5][143] = int_stack + 7268;
 Libderiv->deriv2_classes[4][6][143] = int_stack + 7583;
 Libderiv->deriv2_classes[5][3][143] = int_stack + 8003;
 Libderiv->deriv2_classes[5][4][143] = int_stack + 8213;
 Libderiv->deriv2_classes[5][5][143] = int_stack + 8528;
 Libderiv->deriv2_classes[5][6][143] = int_stack + 8969;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 9557;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 9657;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 9807;
 Libderiv->deriv2_classes[3][6][131] = int_stack + 10017;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 10297;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 10447;
 Libderiv->deriv2_classes[4][5][131] = int_stack + 10672;
 Libderiv->deriv2_classes[4][6][131] = int_stack + 10987;
 Libderiv->deriv2_classes[5][3][131] = int_stack + 11407;
 Libderiv->deriv2_classes[5][4][131] = int_stack + 11617;
 Libderiv->deriv2_classes[5][5][131] = int_stack + 11932;
 Libderiv->deriv2_classes[5][6][131] = int_stack + 12373;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 12961;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 13061;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 13211;
 Libderiv->deriv2_classes[3][6][130] = int_stack + 13421;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 13701;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 13851;
 Libderiv->deriv2_classes[4][5][130] = int_stack + 14076;
 Libderiv->deriv2_classes[4][6][130] = int_stack + 14391;
 Libderiv->deriv2_classes[5][3][130] = int_stack + 14811;
 Libderiv->deriv2_classes[5][4][130] = int_stack + 15021;
 Libderiv->deriv2_classes[5][5][130] = int_stack + 15336;
 Libderiv->deriv2_classes[5][6][130] = int_stack + 15777;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 16365;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 16465;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 16615;
 Libderiv->deriv2_classes[3][6][119] = int_stack + 16825;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 17105;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 17255;
 Libderiv->deriv2_classes[4][5][119] = int_stack + 17480;
 Libderiv->deriv2_classes[4][6][119] = int_stack + 17795;
 Libderiv->deriv2_classes[5][3][119] = int_stack + 18215;
 Libderiv->deriv2_classes[5][4][119] = int_stack + 18425;
 Libderiv->deriv2_classes[5][5][119] = int_stack + 18740;
 Libderiv->deriv2_classes[5][6][119] = int_stack + 19181;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 19769;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 19869;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 20019;
 Libderiv->deriv2_classes[3][6][118] = int_stack + 20229;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 20509;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 20659;
 Libderiv->deriv2_classes[4][5][118] = int_stack + 20884;
 Libderiv->deriv2_classes[4][6][118] = int_stack + 21199;
 Libderiv->deriv2_classes[5][3][118] = int_stack + 21619;
 Libderiv->deriv2_classes[5][4][118] = int_stack + 21829;
 Libderiv->deriv2_classes[5][5][118] = int_stack + 22144;
 Libderiv->deriv2_classes[5][6][118] = int_stack + 22585;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 23173;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 23273;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 23423;
 Libderiv->deriv2_classes[3][6][117] = int_stack + 23633;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 23913;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 24063;
 Libderiv->deriv2_classes[4][5][117] = int_stack + 24288;
 Libderiv->deriv2_classes[4][6][117] = int_stack + 24603;
 Libderiv->deriv2_classes[5][3][117] = int_stack + 25023;
 Libderiv->deriv2_classes[5][4][117] = int_stack + 25233;
 Libderiv->deriv2_classes[5][5][117] = int_stack + 25548;
 Libderiv->deriv2_classes[5][6][117] = int_stack + 25989;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 26577;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 26677;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 26827;
 Libderiv->deriv2_classes[3][6][107] = int_stack + 27037;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 27317;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 27467;
 Libderiv->deriv2_classes[4][5][107] = int_stack + 27692;
 Libderiv->deriv2_classes[4][6][107] = int_stack + 28007;
 Libderiv->deriv2_classes[5][3][107] = int_stack + 28427;
 Libderiv->deriv2_classes[5][4][107] = int_stack + 28637;
 Libderiv->deriv2_classes[5][5][107] = int_stack + 28952;
 Libderiv->deriv2_classes[5][6][107] = int_stack + 29393;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 29981;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 30081;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 30231;
 Libderiv->deriv2_classes[3][6][106] = int_stack + 30441;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 30721;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 30871;
 Libderiv->deriv2_classes[4][5][106] = int_stack + 31096;
 Libderiv->deriv2_classes[4][6][106] = int_stack + 31411;
 Libderiv->deriv2_classes[5][3][106] = int_stack + 31831;
 Libderiv->deriv2_classes[5][4][106] = int_stack + 32041;
 Libderiv->deriv2_classes[5][5][106] = int_stack + 32356;
 Libderiv->deriv2_classes[5][6][106] = int_stack + 32797;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 33385;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 33485;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 33635;
 Libderiv->deriv2_classes[3][6][105] = int_stack + 33845;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 34125;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 34275;
 Libderiv->deriv2_classes[4][5][105] = int_stack + 34500;
 Libderiv->deriv2_classes[4][6][105] = int_stack + 34815;
 Libderiv->deriv2_classes[5][3][105] = int_stack + 35235;
 Libderiv->deriv2_classes[5][4][105] = int_stack + 35445;
 Libderiv->deriv2_classes[5][5][105] = int_stack + 35760;
 Libderiv->deriv2_classes[5][6][105] = int_stack + 36201;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 36789;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 36889;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 37039;
 Libderiv->deriv2_classes[3][6][104] = int_stack + 37249;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 37529;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 37679;
 Libderiv->deriv2_classes[4][5][104] = int_stack + 37904;
 Libderiv->deriv2_classes[4][6][104] = int_stack + 38219;
 Libderiv->deriv2_classes[5][3][104] = int_stack + 38639;
 Libderiv->deriv2_classes[5][4][104] = int_stack + 38849;
 Libderiv->deriv2_classes[5][5][104] = int_stack + 39164;
 Libderiv->deriv2_classes[5][6][104] = int_stack + 39605;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 40193;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 40293;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 40443;
 Libderiv->deriv2_classes[3][6][95] = int_stack + 40653;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 40933;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 41083;
 Libderiv->deriv2_classes[4][5][95] = int_stack + 41308;
 Libderiv->deriv2_classes[4][6][95] = int_stack + 41623;
 Libderiv->deriv2_classes[5][3][95] = int_stack + 42043;
 Libderiv->deriv2_classes[5][4][95] = int_stack + 42253;
 Libderiv->deriv2_classes[5][5][95] = int_stack + 42568;
 Libderiv->deriv2_classes[5][6][95] = int_stack + 43009;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 43597;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 43697;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 43847;
 Libderiv->deriv2_classes[3][6][94] = int_stack + 44057;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 44337;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 44487;
 Libderiv->deriv2_classes[4][5][94] = int_stack + 44712;
 Libderiv->deriv2_classes[4][6][94] = int_stack + 45027;
 Libderiv->deriv2_classes[5][3][94] = int_stack + 45447;
 Libderiv->deriv2_classes[5][4][94] = int_stack + 45657;
 Libderiv->deriv2_classes[5][5][94] = int_stack + 45972;
 Libderiv->deriv2_classes[5][6][94] = int_stack + 46413;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 47001;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 47101;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 47251;
 Libderiv->deriv2_classes[3][6][93] = int_stack + 47461;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 47741;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 47891;
 Libderiv->deriv2_classes[4][5][93] = int_stack + 48116;
 Libderiv->deriv2_classes[4][6][93] = int_stack + 48431;
 Libderiv->deriv2_classes[5][3][93] = int_stack + 48851;
 Libderiv->deriv2_classes[5][4][93] = int_stack + 49061;
 Libderiv->deriv2_classes[5][5][93] = int_stack + 49376;
 Libderiv->deriv2_classes[5][6][93] = int_stack + 49817;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 50405;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 50505;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 50655;
 Libderiv->deriv2_classes[3][6][92] = int_stack + 50865;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 51145;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 51295;
 Libderiv->deriv2_classes[4][5][92] = int_stack + 51520;
 Libderiv->deriv2_classes[4][6][92] = int_stack + 51835;
 Libderiv->deriv2_classes[5][3][92] = int_stack + 52255;
 Libderiv->deriv2_classes[5][4][92] = int_stack + 52465;
 Libderiv->deriv2_classes[5][5][92] = int_stack + 52780;
 Libderiv->deriv2_classes[5][6][92] = int_stack + 53221;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 53809;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 53909;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 54059;
 Libderiv->deriv2_classes[3][6][91] = int_stack + 54269;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 54549;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 54699;
 Libderiv->deriv2_classes[4][5][91] = int_stack + 54924;
 Libderiv->deriv2_classes[4][6][91] = int_stack + 55239;
 Libderiv->deriv2_classes[5][3][91] = int_stack + 55659;
 Libderiv->deriv2_classes[5][4][91] = int_stack + 55869;
 Libderiv->deriv2_classes[5][5][91] = int_stack + 56184;
 Libderiv->deriv2_classes[5][6][91] = int_stack + 56625;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 57213;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 57313;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 57463;
 Libderiv->deriv2_classes[3][6][83] = int_stack + 57673;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 57953;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 58103;
 Libderiv->deriv2_classes[4][5][83] = int_stack + 58328;
 Libderiv->deriv2_classes[4][6][83] = int_stack + 58643;
 Libderiv->deriv_classes[5][3][11] = int_stack + 59063;
 Libderiv->deriv2_classes[5][3][83] = int_stack + 59273;
 Libderiv->deriv_classes[5][4][11] = int_stack + 59483;
 Libderiv->deriv2_classes[5][4][83] = int_stack + 59798;
 Libderiv->deriv_classes[5][5][11] = int_stack + 60113;
 Libderiv->deriv2_classes[5][5][83] = int_stack + 60554;
 Libderiv->deriv2_classes[5][6][83] = int_stack + 60995;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 61583;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 61683;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 61833;
 Libderiv->deriv2_classes[3][6][82] = int_stack + 62043;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 62323;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 62473;
 Libderiv->deriv2_classes[4][5][82] = int_stack + 62698;
 Libderiv->deriv2_classes[4][6][82] = int_stack + 63013;
 Libderiv->deriv_classes[5][3][10] = int_stack + 63433;
 Libderiv->deriv2_classes[5][3][82] = int_stack + 63643;
 Libderiv->deriv_classes[5][4][10] = int_stack + 63853;
 Libderiv->deriv2_classes[5][4][82] = int_stack + 64168;
 Libderiv->deriv_classes[5][5][10] = int_stack + 64483;
 Libderiv->deriv2_classes[5][5][82] = int_stack + 64924;
 Libderiv->deriv2_classes[5][6][82] = int_stack + 65365;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 65953;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 66053;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 66203;
 Libderiv->deriv2_classes[3][6][81] = int_stack + 66413;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 66693;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 66843;
 Libderiv->deriv2_classes[4][5][81] = int_stack + 67068;
 Libderiv->deriv2_classes[4][6][81] = int_stack + 67383;
 Libderiv->deriv_classes[5][3][9] = int_stack + 67803;
 Libderiv->deriv2_classes[5][3][81] = int_stack + 68013;
 Libderiv->deriv_classes[5][4][9] = int_stack + 68223;
 Libderiv->deriv2_classes[5][4][81] = int_stack + 68538;
 Libderiv->deriv_classes[5][5][9] = int_stack + 68853;
 Libderiv->deriv2_classes[5][5][81] = int_stack + 69294;
 Libderiv->deriv2_classes[5][6][81] = int_stack + 69735;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 70323;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 70423;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 70573;
 Libderiv->deriv2_classes[3][6][80] = int_stack + 70783;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 71063;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 71213;
 Libderiv->deriv2_classes[4][5][80] = int_stack + 71438;
 Libderiv->deriv2_classes[4][6][80] = int_stack + 71753;
 Libderiv->deriv_classes[5][3][8] = int_stack + 72173;
 Libderiv->deriv2_classes[5][3][80] = int_stack + 72383;
 Libderiv->deriv_classes[5][4][8] = int_stack + 72593;
 Libderiv->deriv2_classes[5][4][80] = int_stack + 72908;
 Libderiv->deriv_classes[5][5][8] = int_stack + 73223;
 Libderiv->deriv2_classes[5][5][80] = int_stack + 73664;
 Libderiv->deriv2_classes[5][6][80] = int_stack + 74105;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 74693;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 74793;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 74943;
 Libderiv->deriv2_classes[3][6][79] = int_stack + 75153;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 75433;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 75583;
 Libderiv->deriv2_classes[4][5][79] = int_stack + 75808;
 Libderiv->deriv2_classes[4][6][79] = int_stack + 76123;
 Libderiv->deriv_classes[5][3][7] = int_stack + 76543;
 Libderiv->deriv2_classes[5][3][79] = int_stack + 76753;
 Libderiv->deriv_classes[5][4][7] = int_stack + 76963;
 Libderiv->deriv2_classes[5][4][79] = int_stack + 77278;
 Libderiv->deriv_classes[5][5][7] = int_stack + 77593;
 Libderiv->deriv2_classes[5][5][79] = int_stack + 78034;
 Libderiv->deriv2_classes[5][6][79] = int_stack + 78475;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 79063;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 79163;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 79313;
 Libderiv->deriv2_classes[3][6][78] = int_stack + 79523;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 79803;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 79953;
 Libderiv->deriv2_classes[4][5][78] = int_stack + 80178;
 Libderiv->deriv2_classes[4][6][78] = int_stack + 80493;
 Libderiv->dvrr_classes[5][3] = int_stack + 80913;
 Libderiv->deriv_classes[5][3][6] = int_stack + 81123;
 Libderiv->deriv2_classes[5][3][78] = int_stack + 81333;
 Libderiv->dvrr_classes[5][4] = int_stack + 81543;
 Libderiv->deriv_classes[5][4][6] = int_stack + 81858;
 Libderiv->deriv2_classes[5][4][78] = int_stack + 82173;
 Libderiv->deriv_classes[5][5][6] = int_stack + 82488;
 Libderiv->deriv2_classes[5][5][78] = int_stack + 82929;
 Libderiv->deriv2_classes[5][6][78] = int_stack + 83370;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 83958;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 84058;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 84208;
 Libderiv->deriv2_classes[3][6][35] = int_stack + 84418;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 84698;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 84848;
 Libderiv->deriv2_classes[4][5][35] = int_stack + 85073;
 Libderiv->deriv2_classes[4][6][35] = int_stack + 85388;
 Libderiv->deriv2_classes[5][3][35] = int_stack + 85808;
 Libderiv->deriv2_classes[5][4][35] = int_stack + 86018;
 Libderiv->deriv2_classes[5][5][35] = int_stack + 86333;
 Libderiv->deriv2_classes[5][6][35] = int_stack + 86774;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 87362;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 87462;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 87612;
 Libderiv->deriv2_classes[3][6][34] = int_stack + 87822;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 88102;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 88252;
 Libderiv->deriv2_classes[4][5][34] = int_stack + 88477;
 Libderiv->deriv2_classes[4][6][34] = int_stack + 88792;
 Libderiv->deriv2_classes[5][3][34] = int_stack + 89212;
 Libderiv->deriv2_classes[5][4][34] = int_stack + 89422;
 Libderiv->deriv2_classes[5][5][34] = int_stack + 89737;
 Libderiv->deriv2_classes[5][6][34] = int_stack + 90178;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 90766;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 90866;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 91016;
 Libderiv->deriv2_classes[3][6][33] = int_stack + 91226;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 91506;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 91656;
 Libderiv->deriv2_classes[4][5][33] = int_stack + 91881;
 Libderiv->deriv2_classes[4][6][33] = int_stack + 92196;
 Libderiv->deriv2_classes[5][3][33] = int_stack + 92616;
 Libderiv->deriv2_classes[5][4][33] = int_stack + 92826;
 Libderiv->deriv2_classes[5][5][33] = int_stack + 93141;
 Libderiv->deriv2_classes[5][6][33] = int_stack + 93582;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 94170;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 94270;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 94420;
 Libderiv->deriv2_classes[3][6][32] = int_stack + 94630;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 94910;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 95060;
 Libderiv->deriv2_classes[4][5][32] = int_stack + 95285;
 Libderiv->deriv2_classes[4][6][32] = int_stack + 95600;
 Libderiv->deriv2_classes[5][3][32] = int_stack + 96020;
 Libderiv->deriv2_classes[5][4][32] = int_stack + 96230;
 Libderiv->deriv2_classes[5][5][32] = int_stack + 96545;
 Libderiv->deriv2_classes[5][6][32] = int_stack + 96986;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 97574;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 97674;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 97824;
 Libderiv->deriv2_classes[3][6][31] = int_stack + 98034;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 98314;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 98464;
 Libderiv->deriv2_classes[4][5][31] = int_stack + 98689;
 Libderiv->deriv2_classes[4][6][31] = int_stack + 99004;
 Libderiv->deriv2_classes[5][3][31] = int_stack + 99424;
 Libderiv->deriv2_classes[5][4][31] = int_stack + 99634;
 Libderiv->deriv2_classes[5][5][31] = int_stack + 99949;
 Libderiv->deriv2_classes[5][6][31] = int_stack + 100390;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 100978;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 101078;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 101228;
 Libderiv->deriv2_classes[3][6][30] = int_stack + 101438;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 101718;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 101868;
 Libderiv->deriv2_classes[4][5][30] = int_stack + 102093;
 Libderiv->deriv2_classes[4][6][30] = int_stack + 102408;
 Libderiv->deriv_classes[5][3][2] = int_stack + 102828;
 Libderiv->deriv2_classes[5][3][30] = int_stack + 103038;
 Libderiv->deriv_classes[5][4][2] = int_stack + 103248;
 Libderiv->deriv2_classes[5][4][30] = int_stack + 103563;
 Libderiv->deriv_classes[5][5][2] = int_stack + 103878;
 Libderiv->deriv2_classes[5][5][30] = int_stack + 104319;
 Libderiv->deriv2_classes[5][6][30] = int_stack + 104760;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 105348;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 105448;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 105598;
 Libderiv->deriv2_classes[3][6][26] = int_stack + 105808;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 106088;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 106238;
 Libderiv->deriv2_classes[4][5][26] = int_stack + 106463;
 Libderiv->deriv2_classes[4][6][26] = int_stack + 106778;
 Libderiv->deriv2_classes[5][3][26] = int_stack + 107198;
 Libderiv->deriv2_classes[5][4][26] = int_stack + 107408;
 Libderiv->deriv2_classes[5][5][26] = int_stack + 107723;
 Libderiv->deriv2_classes[5][6][26] = int_stack + 108164;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 108752;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 108852;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 109002;
 Libderiv->deriv2_classes[3][6][23] = int_stack + 109212;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 109492;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 109642;
 Libderiv->deriv2_classes[4][5][23] = int_stack + 109867;
 Libderiv->deriv2_classes[4][6][23] = int_stack + 110182;
 Libderiv->deriv2_classes[5][3][23] = int_stack + 110602;
 Libderiv->deriv2_classes[5][4][23] = int_stack + 110812;
 Libderiv->deriv2_classes[5][5][23] = int_stack + 111127;
 Libderiv->deriv2_classes[5][6][23] = int_stack + 111568;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 112156;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 112256;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 112406;
 Libderiv->deriv2_classes[3][6][22] = int_stack + 112616;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 112896;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 113046;
 Libderiv->deriv2_classes[4][5][22] = int_stack + 113271;
 Libderiv->deriv2_classes[4][6][22] = int_stack + 113586;
 Libderiv->deriv2_classes[5][3][22] = int_stack + 114006;
 Libderiv->deriv2_classes[5][4][22] = int_stack + 114216;
 Libderiv->deriv2_classes[5][5][22] = int_stack + 114531;
 Libderiv->deriv2_classes[5][6][22] = int_stack + 114972;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 115560;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 115660;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 115810;
 Libderiv->deriv2_classes[3][6][21] = int_stack + 116020;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 116300;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 116450;
 Libderiv->deriv2_classes[4][5][21] = int_stack + 116675;
 Libderiv->deriv2_classes[4][6][21] = int_stack + 116990;
 Libderiv->deriv2_classes[5][3][21] = int_stack + 117410;
 Libderiv->deriv2_classes[5][4][21] = int_stack + 117620;
 Libderiv->deriv2_classes[5][5][21] = int_stack + 117935;
 Libderiv->deriv2_classes[5][6][21] = int_stack + 118376;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 118964;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 119064;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 119214;
 Libderiv->deriv2_classes[3][6][20] = int_stack + 119424;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 119704;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 119854;
 Libderiv->deriv2_classes[4][5][20] = int_stack + 120079;
 Libderiv->deriv2_classes[4][6][20] = int_stack + 120394;
 Libderiv->deriv2_classes[5][3][20] = int_stack + 120814;
 Libderiv->deriv2_classes[5][4][20] = int_stack + 121024;
 Libderiv->deriv2_classes[5][5][20] = int_stack + 121339;
 Libderiv->deriv2_classes[5][6][20] = int_stack + 121780;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 122368;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 122468;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 122618;
 Libderiv->deriv2_classes[3][6][19] = int_stack + 122828;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 123108;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 123258;
 Libderiv->deriv2_classes[4][5][19] = int_stack + 123483;
 Libderiv->deriv2_classes[4][6][19] = int_stack + 123798;
 Libderiv->deriv2_classes[5][3][19] = int_stack + 124218;
 Libderiv->deriv2_classes[5][4][19] = int_stack + 124428;
 Libderiv->deriv2_classes[5][5][19] = int_stack + 124743;
 Libderiv->deriv2_classes[5][6][19] = int_stack + 125184;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 125772;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 125872;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 126022;
 Libderiv->deriv2_classes[3][6][18] = int_stack + 126232;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 126512;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 126662;
 Libderiv->deriv2_classes[4][5][18] = int_stack + 126887;
 Libderiv->deriv2_classes[4][6][18] = int_stack + 127202;
 Libderiv->deriv_classes[5][3][1] = int_stack + 127622;
 Libderiv->deriv2_classes[5][3][18] = int_stack + 127832;
 Libderiv->deriv_classes[5][4][1] = int_stack + 128042;
 Libderiv->deriv2_classes[5][4][18] = int_stack + 128357;
 Libderiv->deriv_classes[5][5][1] = int_stack + 128672;
 Libderiv->deriv2_classes[5][5][18] = int_stack + 129113;
 Libderiv->deriv2_classes[5][6][18] = int_stack + 129554;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 130142;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 130242;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 130392;
 Libderiv->deriv2_classes[3][6][14] = int_stack + 130602;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 130882;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 131032;
 Libderiv->deriv2_classes[4][5][14] = int_stack + 131257;
 Libderiv->deriv2_classes[4][6][14] = int_stack + 131572;
 Libderiv->deriv2_classes[5][3][14] = int_stack + 131992;
 Libderiv->deriv2_classes[5][4][14] = int_stack + 132202;
 Libderiv->deriv2_classes[5][5][14] = int_stack + 132517;
 Libderiv->deriv2_classes[5][6][14] = int_stack + 132958;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 133546;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 133646;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 133796;
 Libderiv->deriv2_classes[3][6][13] = int_stack + 134006;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 134286;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 134436;
 Libderiv->deriv2_classes[4][5][13] = int_stack + 134661;
 Libderiv->deriv2_classes[4][6][13] = int_stack + 134976;
 Libderiv->deriv2_classes[5][3][13] = int_stack + 135396;
 Libderiv->deriv2_classes[5][4][13] = int_stack + 135606;
 Libderiv->deriv2_classes[5][5][13] = int_stack + 135921;
 Libderiv->deriv2_classes[5][6][13] = int_stack + 136362;
 Libderiv->deriv_classes[3][3][11] = int_stack + 136950;
 Libderiv->deriv_classes[3][4][11] = int_stack + 137050;
 Libderiv->deriv_classes[3][5][11] = int_stack + 137200;
 Libderiv->deriv_classes[3][6][11] = int_stack + 137410;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 137690;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 137790;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 137940;
 Libderiv->deriv2_classes[3][6][11] = int_stack + 138150;
 Libderiv->deriv_classes[4][3][11] = int_stack + 138430;
 Libderiv->deriv_classes[4][4][11] = int_stack + 138580;
 Libderiv->deriv_classes[4][5][11] = int_stack + 138805;
 Libderiv->deriv_classes[4][6][11] = int_stack + 139120;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 139540;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 139690;
 Libderiv->deriv2_classes[4][5][11] = int_stack + 139915;
 Libderiv->deriv2_classes[4][6][11] = int_stack + 140230;
 Libderiv->deriv2_classes[5][3][11] = int_stack + 140650;
 Libderiv->deriv2_classes[5][4][11] = int_stack + 140860;
 Libderiv->deriv2_classes[5][5][11] = int_stack + 141175;
 Libderiv->deriv2_classes[5][6][11] = int_stack + 141616;
 Libderiv->deriv_classes[3][3][10] = int_stack + 142204;
 Libderiv->deriv_classes[3][4][10] = int_stack + 142304;
 Libderiv->deriv_classes[3][5][10] = int_stack + 142454;
 Libderiv->deriv_classes[3][6][10] = int_stack + 142664;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 142944;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 143044;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 143194;
 Libderiv->deriv2_classes[3][6][10] = int_stack + 143404;
 Libderiv->deriv_classes[4][3][10] = int_stack + 143684;
 Libderiv->deriv_classes[4][4][10] = int_stack + 143834;
 Libderiv->deriv_classes[4][5][10] = int_stack + 144059;
 Libderiv->deriv_classes[4][6][10] = int_stack + 144374;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 144794;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 144944;
 Libderiv->deriv2_classes[4][5][10] = int_stack + 145169;
 Libderiv->deriv2_classes[4][6][10] = int_stack + 145484;
 Libderiv->deriv2_classes[5][3][10] = int_stack + 145904;
 Libderiv->deriv2_classes[5][4][10] = int_stack + 146114;
 Libderiv->deriv2_classes[5][5][10] = int_stack + 146429;
 Libderiv->deriv2_classes[5][6][10] = int_stack + 146870;
 Libderiv->deriv_classes[3][3][9] = int_stack + 147458;
 Libderiv->deriv_classes[3][4][9] = int_stack + 147558;
 Libderiv->deriv_classes[3][5][9] = int_stack + 147708;
 Libderiv->deriv_classes[3][6][9] = int_stack + 147918;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 148198;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 148298;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 148448;
 Libderiv->deriv2_classes[3][6][9] = int_stack + 148658;
 Libderiv->deriv_classes[4][3][9] = int_stack + 148938;
 Libderiv->deriv_classes[4][4][9] = int_stack + 149088;
 Libderiv->deriv_classes[4][5][9] = int_stack + 149313;
 Libderiv->deriv_classes[4][6][9] = int_stack + 149628;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 150048;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 150198;
 Libderiv->deriv2_classes[4][5][9] = int_stack + 150423;
 Libderiv->deriv2_classes[4][6][9] = int_stack + 150738;
 Libderiv->deriv2_classes[5][3][9] = int_stack + 151158;
 Libderiv->deriv2_classes[5][4][9] = int_stack + 151368;
 Libderiv->deriv2_classes[5][5][9] = int_stack + 151683;
 Libderiv->deriv2_classes[5][6][9] = int_stack + 152124;
 Libderiv->deriv_classes[3][3][8] = int_stack + 152712;
 Libderiv->deriv_classes[3][4][8] = int_stack + 152812;
 Libderiv->deriv_classes[3][5][8] = int_stack + 152962;
 Libderiv->deriv_classes[3][6][8] = int_stack + 153172;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 153452;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 153552;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 153702;
 Libderiv->deriv2_classes[3][6][8] = int_stack + 153912;
 Libderiv->deriv_classes[4][3][8] = int_stack + 154192;
 Libderiv->deriv_classes[4][4][8] = int_stack + 154342;
 Libderiv->deriv_classes[4][5][8] = int_stack + 154567;
 Libderiv->deriv_classes[4][6][8] = int_stack + 154882;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 155302;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 155452;
 Libderiv->deriv2_classes[4][5][8] = int_stack + 155677;
 Libderiv->deriv2_classes[4][6][8] = int_stack + 155992;
 Libderiv->deriv2_classes[5][3][8] = int_stack + 156412;
 Libderiv->deriv2_classes[5][4][8] = int_stack + 156622;
 Libderiv->deriv2_classes[5][5][8] = int_stack + 156937;
 Libderiv->deriv2_classes[5][6][8] = int_stack + 157378;
 Libderiv->deriv_classes[3][3][7] = int_stack + 157966;
 Libderiv->deriv_classes[3][4][7] = int_stack + 158066;
 Libderiv->deriv_classes[3][5][7] = int_stack + 158216;
 Libderiv->deriv_classes[3][6][7] = int_stack + 158426;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 158706;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 158806;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 158956;
 Libderiv->deriv2_classes[3][6][7] = int_stack + 159166;
 Libderiv->deriv_classes[4][3][7] = int_stack + 159446;
 Libderiv->deriv_classes[4][4][7] = int_stack + 159596;
 Libderiv->deriv_classes[4][5][7] = int_stack + 159821;
 Libderiv->deriv_classes[4][6][7] = int_stack + 160136;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 160556;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 160706;
 Libderiv->deriv2_classes[4][5][7] = int_stack + 160931;
 Libderiv->deriv2_classes[4][6][7] = int_stack + 161246;
 Libderiv->deriv2_classes[5][3][7] = int_stack + 161666;
 Libderiv->deriv2_classes[5][4][7] = int_stack + 161876;
 Libderiv->deriv2_classes[5][5][7] = int_stack + 162191;
 Libderiv->deriv2_classes[5][6][7] = int_stack + 162632;
 Libderiv->deriv_classes[3][3][6] = int_stack + 163220;
 Libderiv->deriv_classes[3][4][6] = int_stack + 163320;
 Libderiv->deriv_classes[3][5][6] = int_stack + 163470;
 Libderiv->deriv_classes[3][6][6] = int_stack + 163680;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 163960;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 164060;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 164210;
 Libderiv->deriv2_classes[3][6][6] = int_stack + 164420;
 Libderiv->dvrr_classes[4][3] = int_stack + 164700;
 Libderiv->deriv_classes[4][3][6] = int_stack + 164850;
 Libderiv->dvrr_classes[4][4] = int_stack + 165000;
 Libderiv->deriv_classes[4][4][6] = int_stack + 165225;
 Libderiv->dvrr_classes[4][5] = int_stack + 165450;
 Libderiv->deriv_classes[4][5][6] = int_stack + 165765;
 Libderiv->deriv_classes[4][6][6] = int_stack + 166080;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 166500;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 166650;
 Libderiv->deriv2_classes[4][5][6] = int_stack + 166875;
 Libderiv->deriv2_classes[4][6][6] = int_stack + 167190;
 Libderiv->deriv_classes[5][3][0] = int_stack + 167610;
 Libderiv->deriv2_classes[5][3][6] = int_stack + 167820;
 Libderiv->deriv_classes[5][4][0] = int_stack + 168030;
 Libderiv->deriv2_classes[5][4][6] = int_stack + 168345;
 Libderiv->deriv_classes[5][5][0] = int_stack + 168660;
 Libderiv->deriv2_classes[5][5][6] = int_stack + 169101;
 Libderiv->deriv2_classes[5][6][6] = int_stack + 169542;
 Libderiv->deriv_classes[3][3][2] = int_stack + 170130;
 Libderiv->deriv_classes[3][4][2] = int_stack + 170230;
 Libderiv->deriv_classes[3][5][2] = int_stack + 170380;
 Libderiv->deriv_classes[3][6][2] = int_stack + 170590;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 170870;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 170970;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 171120;
 Libderiv->deriv2_classes[3][6][2] = int_stack + 171330;
 Libderiv->deriv_classes[4][3][2] = int_stack + 171610;
 Libderiv->deriv_classes[4][4][2] = int_stack + 171760;
 Libderiv->deriv_classes[4][5][2] = int_stack + 171985;
 Libderiv->deriv_classes[4][6][2] = int_stack + 172300;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 172720;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 172870;
 Libderiv->deriv2_classes[4][5][2] = int_stack + 173095;
 Libderiv->deriv2_classes[4][6][2] = int_stack + 173410;
 Libderiv->deriv2_classes[5][3][2] = int_stack + 173830;
 Libderiv->deriv2_classes[5][4][2] = int_stack + 174040;
 Libderiv->deriv2_classes[5][5][2] = int_stack + 174355;
 Libderiv->deriv2_classes[5][6][2] = int_stack + 174796;
 Libderiv->deriv_classes[3][3][1] = int_stack + 175384;
 Libderiv->deriv_classes[3][4][1] = int_stack + 175484;
 Libderiv->deriv_classes[3][5][1] = int_stack + 175634;
 Libderiv->deriv_classes[3][6][1] = int_stack + 175844;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 176124;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 176224;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 176374;
 Libderiv->deriv2_classes[3][6][1] = int_stack + 176584;
 Libderiv->deriv_classes[4][3][1] = int_stack + 176864;
 Libderiv->deriv_classes[4][4][1] = int_stack + 177014;
 Libderiv->deriv_classes[4][5][1] = int_stack + 177239;
 Libderiv->deriv_classes[4][6][1] = int_stack + 177554;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 177974;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 178124;
 Libderiv->deriv2_classes[4][5][1] = int_stack + 178349;
 Libderiv->deriv2_classes[4][6][1] = int_stack + 178664;
 Libderiv->deriv2_classes[5][3][1] = int_stack + 179084;
 Libderiv->deriv2_classes[5][4][1] = int_stack + 179294;
 Libderiv->deriv2_classes[5][5][1] = int_stack + 179609;
 Libderiv->deriv2_classes[5][6][1] = int_stack + 180050;
 Libderiv->dvrr_classes[3][3] = int_stack + 180638;
 Libderiv->dvrr_classes[3][4] = int_stack + 180738;
 Libderiv->dvrr_classes[3][5] = int_stack + 180888;
 Libderiv->dvrr_classes[3][6] = int_stack + 181098;
 Libderiv->deriv_classes[3][3][0] = int_stack + 181378;
 Libderiv->deriv_classes[3][4][0] = int_stack + 181478;
 Libderiv->deriv_classes[3][5][0] = int_stack + 181628;
 Libderiv->deriv_classes[3][6][0] = int_stack + 181838;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 182118;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 182218;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 182368;
 Libderiv->deriv2_classes[3][6][0] = int_stack + 182578;
 Libderiv->deriv_classes[4][3][0] = int_stack + 182858;
 Libderiv->deriv_classes[4][4][0] = int_stack + 183008;
 Libderiv->deriv_classes[4][5][0] = int_stack + 183233;
 Libderiv->deriv_classes[4][6][0] = int_stack + 183548;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 183968;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 184118;
 Libderiv->deriv2_classes[4][5][0] = int_stack + 184343;
 Libderiv->deriv2_classes[4][6][0] = int_stack + 184658;
 Libderiv->deriv2_classes[5][3][0] = int_stack + 185078;
 Libderiv->deriv2_classes[5][4][0] = int_stack + 185288;
 Libderiv->deriv2_classes[5][5][0] = int_stack + 185603;
 Libderiv->deriv2_classes[5][6][0] = int_stack + 186044;
 memset(int_stack,0,1493056);

 Libderiv->dvrr_stack = int_stack + 477319;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_fdff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+186632,int_stack+180738,int_stack+180638,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+186932,int_stack+180888,int_stack+180738,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+187382,int_stack+186932,int_stack+186632,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+187982,int_stack+137050,int_stack+136950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180638,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+188282,int_stack+137200,int_stack+137050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180738,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+188732,int_stack+188282,int_stack+187982, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+189332,int_stack+137410,int_stack+137200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180888,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+189962,int_stack+189332,int_stack+188282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186932,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+190862,int_stack+189962,int_stack+188732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187382,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+189332,int_stack+165000,int_stack+164700,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+189782,int_stack+165450,int_stack+165000,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+191862,int_stack+189782,int_stack+189332,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+192762,int_stack+138580,int_stack+138430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+164700,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+193212,int_stack+138805,int_stack+138580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165000,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+193887,int_stack+193212,int_stack+192762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+194787,int_stack+139120,int_stack+138805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165450,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+195732,int_stack+194787,int_stack+193212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189782,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+197082,int_stack+195732,int_stack+193887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+198582,int_stack+197082,int_stack+190862,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+194787,int_stack+81543,int_stack+80913,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+195417,int_stack+2940,int_stack+81543,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+201582,int_stack+195417,int_stack+194787,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+196362,int_stack+59483,int_stack+59063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80913,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+202842,int_stack+60113,int_stack+59483, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81543,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+203787,int_stack+202842,int_stack+196362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+205047,int_stack+0,int_stack+60113, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2940,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+206370,int_stack+205047,int_stack+202842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195417,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+208260,int_stack+206370,int_stack+203787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+210360,int_stack+208260,int_stack+197082,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+205047,int_stack+142304,int_stack+142204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180638, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+205347,int_stack+142454,int_stack+142304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180738, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+205797,int_stack+205347,int_stack+205047, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+206397,int_stack+142664,int_stack+142454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180888, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+207027,int_stack+206397,int_stack+205347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186932, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+207927,int_stack+207027,int_stack+205797, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187382, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+206397,int_stack+143834,int_stack+143684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+164700, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+206847,int_stack+144059,int_stack+143834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165000, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+208927,int_stack+206847,int_stack+206397, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+214860,int_stack+144374,int_stack+144059, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165450, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+215805,int_stack+214860,int_stack+206847, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189782, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+217155,int_stack+215805,int_stack+208927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+218655,int_stack+217155,int_stack+207927,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+214860,int_stack+63853,int_stack+63433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80913, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+215490,int_stack+64483,int_stack+63853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81543, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+221655,int_stack+215490,int_stack+214860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+222915,int_stack+588,int_stack+64483, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2940, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+224238,int_stack+222915,int_stack+215490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195417, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+226128,int_stack+224238,int_stack+221655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+228228,int_stack+226128,int_stack+217155,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+222915,int_stack+147558,int_stack+147458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180638, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+223215,int_stack+147708,int_stack+147558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180738, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+223665,int_stack+223215,int_stack+222915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+224265,int_stack+147918,int_stack+147708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180888, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+224895,int_stack+224265,int_stack+223215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186932, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+225795,int_stack+224895,int_stack+223665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187382, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+224265,int_stack+149088,int_stack+148938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+164700, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+224715,int_stack+149313,int_stack+149088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165000, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+226795,int_stack+224715,int_stack+224265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+149628,int_stack+149313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165450, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+232728,int_stack+0,int_stack+224715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189782, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234078,int_stack+232728,int_stack+226795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+235578,int_stack+234078,int_stack+225795,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+232728,int_stack+68223,int_stack+67803, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80913, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+68853,int_stack+68223, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81543, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+238578,int_stack+0,int_stack+232728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239838,int_stack+1176,int_stack+68853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2940, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+241161,int_stack+239838,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195417, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+243051,int_stack+241161,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+245151,int_stack+243051,int_stack+234078,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239838,int_stack+152812,int_stack+152712, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+240138,int_stack+152962,int_stack+152812, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+240588,int_stack+240138,int_stack+239838, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+241188,int_stack+153172,int_stack+152962, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+241818,int_stack+241188,int_stack+240138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+242718,int_stack+241818,int_stack+240588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187382, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+241188,int_stack+154342,int_stack+154192, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+164700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+241638,int_stack+154567,int_stack+154342, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+243718,int_stack+241638,int_stack+241188, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+249651,int_stack+154882,int_stack+154567, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+250596,int_stack+249651,int_stack+241638, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+251946,int_stack+250596,int_stack+243718, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+253446,int_stack+251946,int_stack+242718,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+249651,int_stack+72593,int_stack+72173, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80913, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+250281,int_stack+73223,int_stack+72593, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+256446,int_stack+250281,int_stack+249651, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+257706,int_stack+1764,int_stack+73223, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+259029,int_stack+257706,int_stack+250281, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195417, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+260919,int_stack+259029,int_stack+256446, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+263019,int_stack+260919,int_stack+251946,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+257706,int_stack+158066,int_stack+157966, 0.0, zero_stack, 1.0, int_stack+180638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+258006,int_stack+158216,int_stack+158066, 0.0, zero_stack, 1.0, int_stack+180738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+258456,int_stack+258006,int_stack+257706, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+259056,int_stack+158426,int_stack+158216, 0.0, zero_stack, 1.0, int_stack+180888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+259686,int_stack+259056,int_stack+258006, 0.0, zero_stack, 1.0, int_stack+186932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+260586,int_stack+259686,int_stack+258456, 0.0, zero_stack, 1.0, int_stack+187382, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+259056,int_stack+159596,int_stack+159446, 0.0, zero_stack, 1.0, int_stack+164700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+259506,int_stack+159821,int_stack+159596, 0.0, zero_stack, 1.0, int_stack+165000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+261586,int_stack+259506,int_stack+259056, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+945,int_stack+160136,int_stack+159821, 0.0, zero_stack, 1.0, int_stack+165450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+945,int_stack+259506, 0.0, zero_stack, 1.0, int_stack+189782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+268869,int_stack+267519,int_stack+261586, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+270369,int_stack+268869,int_stack+260586,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+267519,int_stack+76963,int_stack+76543, 0.0, zero_stack, 1.0, int_stack+80913, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+945,int_stack+77593,int_stack+76963, 0.0, zero_stack, 1.0, int_stack+81543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+273369,int_stack+945,int_stack+267519, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+274629,int_stack+2352,int_stack+77593, 0.0, zero_stack, 1.0, int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+275952,int_stack+274629,int_stack+945, 0.0, zero_stack, 1.0, int_stack+195417, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+277842,int_stack+275952,int_stack+273369, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+279942,int_stack+277842,int_stack+268869,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+274629,int_stack+163320,int_stack+163220, 1.0, int_stack+180638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+274929,int_stack+163470,int_stack+163320, 1.0, int_stack+180738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+275379,int_stack+274929,int_stack+274629, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+275979,int_stack+163680,int_stack+163470, 1.0, int_stack+180888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+276609,int_stack+275979,int_stack+274929, 1.0, int_stack+186932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+277509,int_stack+276609,int_stack+275379, 1.0, int_stack+187382, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+275979,int_stack+165225,int_stack+164850, 1.0, int_stack+164700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+276429,int_stack+165765,int_stack+165225, 1.0, int_stack+165000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+278509,int_stack+276429,int_stack+275979, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1890,int_stack+166080,int_stack+165765, 1.0, int_stack+165450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+284442,int_stack+1890,int_stack+276429, 1.0, int_stack+189782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+285792,int_stack+284442,int_stack+278509, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+287292,int_stack+285792,int_stack+277509,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+284442,int_stack+81858,int_stack+81123, 1.0, int_stack+80913, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1890,int_stack+82488,int_stack+81858, 1.0, int_stack+81543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+290292,int_stack+1890,int_stack+284442, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+291552,int_stack+3381,int_stack+82488, 1.0, int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+292875,int_stack+291552,int_stack+1890, 1.0, int_stack+195417, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+294765,int_stack+292875,int_stack+290292, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+296865,int_stack+294765,int_stack+285792,100);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+201582,int_stack+181098,int_stack+180888,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+291552,int_stack+201582,int_stack+186932,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+201582,int_stack+291552,int_stack+187382,10);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+291552,int_stack+5145,int_stack+165450,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+186632,int_stack+291552,int_stack+189782,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+186632,int_stack+191862,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+293052,int_stack+291552,int_stack+201582,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+170230,int_stack+170130,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+192162,int_stack+170380,int_stack+170230,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+186632,int_stack+192162,int_stack+191862,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+187232,int_stack+170590,int_stack+170380,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+194787,int_stack+187232,int_stack+192162,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+2835,int_stack+194787,int_stack+186632,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+194787,int_stack+171760,int_stack+171610,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+195237,int_stack+171985,int_stack+171760,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+189332,int_stack+195237,int_stack+194787,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+301365,int_stack+172300,int_stack+171985,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+302310,int_stack+301365,int_stack+195237,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+303660,int_stack+302310,int_stack+189332,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+305160,int_stack+303660,int_stack+2835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+190232,int_stack+103248,int_stack+102828,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+301365,int_stack+103878,int_stack+103248,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+302310,int_stack+301365,int_stack+190232,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+308160,int_stack+3969,int_stack+103878,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+309483,int_stack+308160,int_stack+301365,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+311373,int_stack+309483,int_stack+302310,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+313473,int_stack+311373,int_stack+303660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+308160,int_stack+175484,int_stack+175384,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+195912,int_stack+175634,int_stack+175484,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+308460,int_stack+195912,int_stack+308160,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+309060,int_stack+175844,int_stack+175634,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+309690,int_stack+309060,int_stack+195912,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+310590,int_stack+309690,int_stack+308460,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+309060,int_stack+177014,int_stack+176864,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+309510,int_stack+177239,int_stack+177014,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+311590,int_stack+309510,int_stack+309060,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+312490,int_stack+177554,int_stack+177239,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+317973,int_stack+312490,int_stack+309510,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+319323,int_stack+317973,int_stack+311590,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+320823,int_stack+319323,int_stack+310590, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+317973,int_stack+128042,int_stack+127622,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+312490,int_stack+128672,int_stack+128042,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+312490,int_stack+317973,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+325083,int_stack+4557,int_stack+128672,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+326406,int_stack+325083,int_stack+312490,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+328296,int_stack+326406,int_stack+323823,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+330396,int_stack+328296,int_stack+319323, 0.0, zero_stack, 1.0, int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+325083,int_stack+181478,int_stack+181378,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+325383,int_stack+181628,int_stack+181478,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+325833,int_stack+325383,int_stack+325083,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+326433,int_stack+181838,int_stack+181628,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+327063,int_stack+326433,int_stack+325383,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+327963,int_stack+327063,int_stack+325833,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+326433,int_stack+183008,int_stack+182858,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+326883,int_stack+183233,int_stack+183008,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+328963,int_stack+326883,int_stack+326433,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3835,int_stack+183548,int_stack+183233,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+334896,int_stack+3835,int_stack+326883,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+3835,int_stack+334896,int_stack+328963,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+334896,int_stack+3835,int_stack+327963, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+201582,int_stack+168030,int_stack+167610,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+337896,int_stack+168660,int_stack+168030,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+338841,int_stack+337896,int_stack+201582,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+340101,int_stack+5565,int_stack+168660,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+341424,int_stack+340101,int_stack+337896,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+343314,int_stack+341424,int_stack+338841,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+345414,int_stack+343314,int_stack+3835, 1.0, int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+6253,int_stack+6153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+136950,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+6403,int_stack+6253, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+137050,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+187982,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+6613,int_stack+6403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+137200,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+340101,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+188282,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+341001,int_stack+340101,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+188732,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+340101,int_stack+7043,int_stack+6893, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+138430,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+7268,int_stack+7043, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+138580,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+342001,int_stack+291552,int_stack+340101, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+192762,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+342901,int_stack+7583,int_stack+7268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+138805,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+343846,int_stack+342901,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+193212,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+343846,int_stack+342001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+193887,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+342001,int_stack+291552,int_stack+341001,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+8213,int_stack+8003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+59063,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+340101,int_stack+8528,int_stack+8213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+59483,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5335,int_stack+340101,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+196362,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6595,int_stack+8969,int_stack+8528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+60113,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+349914,int_stack+6595,int_stack+340101, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+202842,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6595,int_stack+349914,int_stack+5335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+203787,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+349914,int_stack+6595,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+9657,int_stack+9557, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136950, 1.0, int_stack+142204,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+9807,int_stack+9657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137050, 1.0, int_stack+142304,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187982, 1.0, int_stack+205047,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+10017,int_stack+9807, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137200, 1.0, int_stack+142454,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5335,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+188282, 1.0, int_stack+205347,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6235,int_stack+5335,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+188732, 1.0, int_stack+205797,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5335,int_stack+10447,int_stack+10297, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138430, 1.0, int_stack+143684,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+10672,int_stack+10447, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138580, 1.0, int_stack+143834,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7235,int_stack+291552,int_stack+5335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192762, 1.0, int_stack+206397,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8135,int_stack+10987,int_stack+10672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138805, 1.0, int_stack+144059,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9080,int_stack+8135,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193212, 1.0, int_stack+206847,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+9080,int_stack+7235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193887, 1.0, int_stack+208927,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7235,int_stack+291552,int_stack+6235,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+11617,int_stack+11407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59063, 1.0, int_stack+63433,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10235,int_stack+11932,int_stack+11617, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59483, 1.0, int_stack+63853,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5335,int_stack+10235,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+196362, 1.0, int_stack+214860,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+340101,int_stack+12373,int_stack+11932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60113, 1.0, int_stack+64483,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+354414,int_stack+340101,int_stack+10235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202842, 1.0, int_stack+215490,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10235,int_stack+354414,int_stack+5335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203787, 1.0, int_stack+221655,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+354414,int_stack+10235,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+13061,int_stack+12961, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+142204, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+13211,int_stack+13061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+142304, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+205047, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+13421,int_stack+13211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+142454, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10235,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+205347, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11135,int_stack+10235,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+205797, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10235,int_stack+13851,int_stack+13701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+143684, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+14076,int_stack+13851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+143834, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12135,int_stack+291552,int_stack+10235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+206397, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13035,int_stack+14391,int_stack+14076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+144059, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5335,int_stack+13035,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+206847, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+5335,int_stack+12135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+208927, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+358914,int_stack+291552,int_stack+11135,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+15021,int_stack+14811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+63433, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5335,int_stack+15336,int_stack+15021, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+63853, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10235,int_stack+5335,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+214860, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11495,int_stack+15777,int_stack+15336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+64483, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12818,int_stack+11495,int_stack+5335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+215490, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+361914,int_stack+12818,int_stack+10235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+221655, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+10235,int_stack+361914,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+16465,int_stack+16365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136950, 0.0, zero_stack, 1.0, int_stack+147458,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+16615,int_stack+16465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137050, 0.0, zero_stack, 1.0, int_stack+147558,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187982, 0.0, zero_stack, 1.0, int_stack+222915,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+16825,int_stack+16615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137200, 0.0, zero_stack, 1.0, int_stack+147708,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+361914,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+188282, 0.0, zero_stack, 1.0, int_stack+223215,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+362814,int_stack+361914,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+188732, 0.0, zero_stack, 1.0, int_stack+223665,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+361914,int_stack+17255,int_stack+17105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138430, 0.0, zero_stack, 1.0, int_stack+148938,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+17480,int_stack+17255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138580, 0.0, zero_stack, 1.0, int_stack+149088,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14735,int_stack+291552,int_stack+361914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192762, 0.0, zero_stack, 1.0, int_stack+224265,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15635,int_stack+17795,int_stack+17480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138805, 0.0, zero_stack, 1.0, int_stack+149313,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16580,int_stack+15635,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193212, 0.0, zero_stack, 1.0, int_stack+224715,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+16580,int_stack+14735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193887, 0.0, zero_stack, 1.0, int_stack+226795,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+14735,int_stack+291552,int_stack+362814,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+18425,int_stack+18215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59063, 0.0, zero_stack, 1.0, int_stack+67803,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+361914,int_stack+18740,int_stack+18425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59483, 0.0, zero_stack, 1.0, int_stack+68223,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+362859,int_stack+361914,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+196362, 0.0, zero_stack, 1.0, int_stack+232728,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5335,int_stack+19181,int_stack+18740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60113, 0.0, zero_stack, 1.0, int_stack+68853,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17735,int_stack+5335,int_stack+361914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202842, 0.0, zero_stack, 1.0, int_stack+0,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+364119,int_stack+17735,int_stack+362859, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203787, 0.0, zero_stack, 1.0, int_stack+238578,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+366219,int_stack+364119,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+19869,int_stack+19769, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+142204, 1.0, int_stack+147458, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+20019,int_stack+19869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+142304, 1.0, int_stack+147558, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+205047, 1.0, int_stack+222915, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+20229,int_stack+20019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+142454, 1.0, int_stack+147708, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17735,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+205347, 1.0, int_stack+223215, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18635,int_stack+17735,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+205797, 1.0, int_stack+223665, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17735,int_stack+20659,int_stack+20509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143684, 1.0, int_stack+148938, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+20884,int_stack+20659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143834, 1.0, int_stack+149088, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19635,int_stack+291552,int_stack+17735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+206397, 1.0, int_stack+224265, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+361914,int_stack+21199,int_stack+20884, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144059, 1.0, int_stack+149313, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+362859,int_stack+361914,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+206847, 1.0, int_stack+224715, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+362859,int_stack+19635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208927, 1.0, int_stack+226795, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+361914,int_stack+291552,int_stack+18635,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+21829,int_stack+21619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63433, 1.0, int_stack+67803, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+364914,int_stack+22144,int_stack+21829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63853, 1.0, int_stack+68223, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17735,int_stack+364914,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+214860, 1.0, int_stack+232728, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18995,int_stack+22585,int_stack+22144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64483, 1.0, int_stack+68853, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20318,int_stack+18995,int_stack+364914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+215490, 1.0, int_stack+0, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+370719,int_stack+20318,int_stack+17735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+221655, 1.0, int_stack+238578, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+17735,int_stack+370719,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+23273,int_stack+23173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+147458, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+23423,int_stack+23273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+147558, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+222915, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+23633,int_stack+23423, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+147708, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+370719,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+223215, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+371619,int_stack+370719,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+223665, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+370719,int_stack+24063,int_stack+23913, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+148938, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+24288,int_stack+24063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+149088, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22235,int_stack+291552,int_stack+370719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+224265, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23135,int_stack+24603,int_stack+24288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+149313, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5335,int_stack+23135,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+224715, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+5335,int_stack+22235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+226795, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+372619,int_stack+291552,int_stack+371619,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+25233,int_stack+25023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+67803, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22235,int_stack+25548,int_stack+25233, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+68223, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23180,int_stack+22235,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+232728, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5335,int_stack+25989,int_stack+25548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+68853, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24440,int_stack+5335,int_stack+22235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+24440,int_stack+23180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+238578, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+377719,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+26677,int_stack+26577, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152712,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+26827,int_stack+26677, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152812,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187982, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+239838,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+27037,int_stack+26827, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152962,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+188282, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240138,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+188732, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240588,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+27467,int_stack+27317, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+154192,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+27692,int_stack+27467, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138580, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+154342,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22235,int_stack+291552,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+241188,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23135,int_stack+28007,int_stack+27692, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138805, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+154567,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24080,int_stack+23135,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+241638,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+24080,int_stack+22235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193887, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243718,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+22235,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+28637,int_stack+28427, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59063, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72173,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25235,int_stack+28952,int_stack+28637, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59483, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72593,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26180,int_stack+25235,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+196362, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+249651,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27440,int_stack+29393,int_stack+28952, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60113, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73223,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+27440,int_stack+25235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202842, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+250281,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27440,int_stack+375619,int_stack+26180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203787, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+256446,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+382219,int_stack+27440,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+30081,int_stack+29981, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+142204, 0.0, zero_stack, 1.0, int_stack+152712, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+30231,int_stack+30081, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+142304, 0.0, zero_stack, 1.0, int_stack+152812, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+205047, 0.0, zero_stack, 1.0, int_stack+239838, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+30441,int_stack+30231, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+142454, 0.0, zero_stack, 1.0, int_stack+152962, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+205347, 0.0, zero_stack, 1.0, int_stack+240138, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+205797, 0.0, zero_stack, 1.0, int_stack+240588, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+30871,int_stack+30721, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143684, 0.0, zero_stack, 1.0, int_stack+154192, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+31096,int_stack+30871, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143834, 0.0, zero_stack, 1.0, int_stack+154342, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25235,int_stack+291552,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+206397, 0.0, zero_stack, 1.0, int_stack+241188, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26135,int_stack+31411,int_stack+31096, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144059, 0.0, zero_stack, 1.0, int_stack+154567, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27080,int_stack+26135,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+206847, 0.0, zero_stack, 1.0, int_stack+241638, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+27080,int_stack+25235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+208927, 0.0, zero_stack, 1.0, int_stack+243718, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+25235,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+32041,int_stack+31831, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63433, 0.0, zero_stack, 1.0, int_stack+72173, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+28235,int_stack+32356,int_stack+32041, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63853, 0.0, zero_stack, 1.0, int_stack+72593, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29180,int_stack+28235,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+214860, 0.0, zero_stack, 1.0, int_stack+249651, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30440,int_stack+32797,int_stack+32356, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64483, 0.0, zero_stack, 1.0, int_stack+73223, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+30440,int_stack+28235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+215490, 0.0, zero_stack, 1.0, int_stack+250281, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+30440,int_stack+375619,int_stack+29180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+221655, 0.0, zero_stack, 1.0, int_stack+256446, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+386719,int_stack+30440,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+33485,int_stack+33385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+147458, 1.0, int_stack+152712, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+33635,int_stack+33485, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+147558, 1.0, int_stack+152812, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+222915, 1.0, int_stack+239838, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+33845,int_stack+33635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+147708, 1.0, int_stack+152962, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+223215, 1.0, int_stack+240138, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+223665, 1.0, int_stack+240588, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+34275,int_stack+34125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148938, 1.0, int_stack+154192, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+34500,int_stack+34275, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149088, 1.0, int_stack+154342, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28235,int_stack+291552,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+224265, 1.0, int_stack+241188, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29135,int_stack+34815,int_stack+34500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149313, 1.0, int_stack+154567, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+30080,int_stack+29135,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+224715, 1.0, int_stack+241638, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+30080,int_stack+28235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+226795, 1.0, int_stack+243718, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28235,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+35445,int_stack+35235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67803, 1.0, int_stack+72173, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31235,int_stack+35760,int_stack+35445, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68223, 1.0, int_stack+72593, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32180,int_stack+31235,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+232728, 1.0, int_stack+249651, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+33440,int_stack+36201,int_stack+35760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68853, 1.0, int_stack+73223, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+34763,int_stack+33440,int_stack+31235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+250281, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+34763,int_stack+32180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+238578, 1.0, int_stack+256446, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+31235,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+36889,int_stack+36789, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+152712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+37039,int_stack+36889, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+152812, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+239838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+37249,int_stack+37039, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+152962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+240138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+240588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+37679,int_stack+37529, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+154192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+37904,int_stack+37679, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+154342, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35735,int_stack+291552,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+241188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36635,int_stack+38219,int_stack+37904, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+154567, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5335,int_stack+36635,int_stack+291552, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+241638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+5335,int_stack+35735, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+243718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+391219,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+38849,int_stack+38639, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+72173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35735,int_stack+39164,int_stack+38849, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+72593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+36680,int_stack+35735,int_stack+202212, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+249651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5335,int_stack+39605,int_stack+39164, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+73223, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+37940,int_stack+5335,int_stack+35735, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+250281, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+37940,int_stack+36680, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+256446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+394219,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+40293,int_stack+40193, 0.0, zero_stack, 1.0, int_stack+136950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+157966,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+40443,int_stack+40293, 0.0, zero_stack, 1.0, int_stack+137050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+158066,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+187982, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+257706,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+40653,int_stack+40443, 0.0, zero_stack, 1.0, int_stack+137200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+158216,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 1.0, int_stack+188282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+258006,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 1.0, int_stack+188732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+258456,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+41083,int_stack+40933, 0.0, zero_stack, 1.0, int_stack+138430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159446,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+41308,int_stack+41083, 0.0, zero_stack, 1.0, int_stack+138580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159596,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35735,int_stack+291552,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+192762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+259056,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36635,int_stack+41623,int_stack+41308, 0.0, zero_stack, 1.0, int_stack+138805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159821,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+37580,int_stack+36635,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+193212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+259506,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+37580,int_stack+35735, 0.0, zero_stack, 1.0, int_stack+193887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+261586,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+35735,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+42253,int_stack+42043, 0.0, zero_stack, 1.0, int_stack+59063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76543,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+38735,int_stack+42568,int_stack+42253, 0.0, zero_stack, 1.0, int_stack+59483, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76963,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+39680,int_stack+38735,int_stack+202212, 0.0, zero_stack, 1.0, int_stack+196362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+267519,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+40940,int_stack+43009,int_stack+42568, 0.0, zero_stack, 1.0, int_stack+60113, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77593,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+40940,int_stack+38735, 0.0, zero_stack, 1.0, int_stack+202842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+40940,int_stack+375619,int_stack+39680, 0.0, zero_stack, 1.0, int_stack+203787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273369,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+398719,int_stack+40940,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+43697,int_stack+43597, 0.0, zero_stack, 1.0, int_stack+142204, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+157966, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+43847,int_stack+43697, 0.0, zero_stack, 1.0, int_stack+142304, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+158066, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+205047, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+257706, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+44057,int_stack+43847, 0.0, zero_stack, 1.0, int_stack+142454, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+158216, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 1.0, int_stack+205347, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+258006, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 1.0, int_stack+205797, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+258456, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+44487,int_stack+44337, 0.0, zero_stack, 1.0, int_stack+143684, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159446, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+44712,int_stack+44487, 0.0, zero_stack, 1.0, int_stack+143834, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159596, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+38735,int_stack+291552,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+206397, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+259056, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+39635,int_stack+45027,int_stack+44712, 0.0, zero_stack, 1.0, int_stack+144059, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159821, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+40580,int_stack+39635,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+206847, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+259506, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+40580,int_stack+38735, 0.0, zero_stack, 1.0, int_stack+208927, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+261586, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+38735,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+45657,int_stack+45447, 0.0, zero_stack, 1.0, int_stack+63433, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76543, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41735,int_stack+45972,int_stack+45657, 0.0, zero_stack, 1.0, int_stack+63853, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76963, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42680,int_stack+41735,int_stack+202212, 0.0, zero_stack, 1.0, int_stack+214860, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+267519, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43940,int_stack+46413,int_stack+45972, 0.0, zero_stack, 1.0, int_stack+64483, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77593, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+43940,int_stack+41735, 0.0, zero_stack, 1.0, int_stack+215490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+43940,int_stack+375619,int_stack+42680, 0.0, zero_stack, 1.0, int_stack+221655, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273369, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+403219,int_stack+43940,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+47101,int_stack+47001, 0.0, zero_stack, 1.0, int_stack+147458, 0.0, zero_stack, 1.0, int_stack+157966, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+47251,int_stack+47101, 0.0, zero_stack, 1.0, int_stack+147558, 0.0, zero_stack, 1.0, int_stack+158066, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+222915, 0.0, zero_stack, 1.0, int_stack+257706, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+47461,int_stack+47251, 0.0, zero_stack, 1.0, int_stack+147708, 0.0, zero_stack, 1.0, int_stack+158216, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 1.0, int_stack+223215, 0.0, zero_stack, 1.0, int_stack+258006, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 1.0, int_stack+223665, 0.0, zero_stack, 1.0, int_stack+258456, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+47891,int_stack+47741, 0.0, zero_stack, 1.0, int_stack+148938, 0.0, zero_stack, 1.0, int_stack+159446, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+48116,int_stack+47891, 0.0, zero_stack, 1.0, int_stack+149088, 0.0, zero_stack, 1.0, int_stack+159596, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41735,int_stack+291552,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+224265, 0.0, zero_stack, 1.0, int_stack+259056, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42635,int_stack+48431,int_stack+48116, 0.0, zero_stack, 1.0, int_stack+149313, 0.0, zero_stack, 1.0, int_stack+159821, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+43580,int_stack+42635,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+224715, 0.0, zero_stack, 1.0, int_stack+259506, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+43580,int_stack+41735, 0.0, zero_stack, 1.0, int_stack+226795, 0.0, zero_stack, 1.0, int_stack+261586, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+41735,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+49061,int_stack+48851, 0.0, zero_stack, 1.0, int_stack+67803, 0.0, zero_stack, 1.0, int_stack+76543, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44735,int_stack+49376,int_stack+49061, 0.0, zero_stack, 1.0, int_stack+68223, 0.0, zero_stack, 1.0, int_stack+76963, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45680,int_stack+44735,int_stack+202212, 0.0, zero_stack, 1.0, int_stack+232728, 0.0, zero_stack, 1.0, int_stack+267519, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+46940,int_stack+49817,int_stack+49376, 0.0, zero_stack, 1.0, int_stack+68853, 0.0, zero_stack, 1.0, int_stack+77593, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48263,int_stack+46940,int_stack+44735, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+48263,int_stack+45680, 0.0, zero_stack, 1.0, int_stack+238578, 0.0, zero_stack, 1.0, int_stack+273369, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+44735,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+50505,int_stack+50405, 0.0, zero_stack, 1.0, int_stack+152712, 1.0, int_stack+157966, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+50655,int_stack+50505, 0.0, zero_stack, 1.0, int_stack+152812, 1.0, int_stack+158066, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+239838, 1.0, int_stack+257706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+50865,int_stack+50655, 0.0, zero_stack, 1.0, int_stack+152962, 1.0, int_stack+158216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 0.0, zero_stack, 1.0, int_stack+240138, 1.0, int_stack+258006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 0.0, zero_stack, 1.0, int_stack+240588, 1.0, int_stack+258456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+51295,int_stack+51145, 0.0, zero_stack, 1.0, int_stack+154192, 1.0, int_stack+159446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+51520,int_stack+51295, 0.0, zero_stack, 1.0, int_stack+154342, 1.0, int_stack+159596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+49235,int_stack+291552,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+241188, 1.0, int_stack+259056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+50135,int_stack+51835,int_stack+51520, 0.0, zero_stack, 1.0, int_stack+154567, 1.0, int_stack+159821, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5335,int_stack+50135,int_stack+291552, 0.0, zero_stack, 1.0, int_stack+241638, 1.0, int_stack+259506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+5335,int_stack+49235, 0.0, zero_stack, 1.0, int_stack+243718, 1.0, int_stack+261586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+49235,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+52465,int_stack+52255, 0.0, zero_stack, 1.0, int_stack+72173, 1.0, int_stack+76543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5335,int_stack+52780,int_stack+52465, 0.0, zero_stack, 1.0, int_stack+72593, 1.0, int_stack+76963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+375619,int_stack+5335,int_stack+202212, 0.0, zero_stack, 1.0, int_stack+249651, 1.0, int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+370719,int_stack+53221,int_stack+52780, 0.0, zero_stack, 1.0, int_stack+73223, 1.0, int_stack+77593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+340101,int_stack+370719,int_stack+5335, 0.0, zero_stack, 1.0, int_stack+250281, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+340101,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+256446, 1.0, int_stack+273369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+409819,int_stack+407719,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+53909,int_stack+53809, 0.0, zero_stack, 2.0, int_stack+157966, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+54059,int_stack+53909, 0.0, zero_stack, 2.0, int_stack+158066, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 0.0, zero_stack, 2.0, int_stack+257706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+54269,int_stack+54059, 0.0, zero_stack, 2.0, int_stack+158216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+202212,int_stack+291852, 0.0, zero_stack, 2.0, int_stack+258006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+408619,int_stack+407719,int_stack+292302, 0.0, zero_stack, 2.0, int_stack+258456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+54699,int_stack+54549, 0.0, zero_stack, 2.0, int_stack+159446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+54924,int_stack+54699, 0.0, zero_stack, 2.0, int_stack+159596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+375619,int_stack+291552,int_stack+407719, 0.0, zero_stack, 2.0, int_stack+259056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376519,int_stack+55239,int_stack+54924, 0.0, zero_stack, 2.0, int_stack+159821, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+340101,int_stack+376519,int_stack+291552, 0.0, zero_stack, 2.0, int_stack+259506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+340101,int_stack+375619, 0.0, zero_stack, 2.0, int_stack+261586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+52235,int_stack+291552,int_stack+408619,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+55869,int_stack+55659, 0.0, zero_stack, 2.0, int_stack+76543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375619,int_stack+56184,int_stack+55869, 0.0, zero_stack, 2.0, int_stack+76963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+340101,int_stack+375619,int_stack+202212, 0.0, zero_stack, 2.0, int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+56625,int_stack+56184, 0.0, zero_stack, 2.0, int_stack+77593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5335,int_stack+407719,int_stack+375619, 0.0, zero_stack, 2.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+5335,int_stack+340101, 0.0, zero_stack, 2.0, int_stack+273369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+414319,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+57313,int_stack+57213, 1.0, int_stack+136950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163220,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+57463,int_stack+57313, 1.0, int_stack+137050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163320,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 1.0, int_stack+187982, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+274629,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202212,int_stack+57673,int_stack+57463, 1.0, int_stack+137200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163470,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+202212,int_stack+291852, 1.0, int_stack+188282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+274929,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+376519,int_stack+375619,int_stack+292302, 1.0, int_stack+188732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275379,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+58103,int_stack+57953, 1.0, int_stack+138430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+164850,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291552,int_stack+58328,int_stack+58103, 1.0, int_stack+138580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165225,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+340101,int_stack+291552,int_stack+375619, 1.0, int_stack+192762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275979,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+341001,int_stack+58643,int_stack+58328, 1.0, int_stack+138805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165765,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5335,int_stack+341001,int_stack+291552, 1.0, int_stack+193212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276429,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+5335,int_stack+340101, 1.0, int_stack+193887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+278509,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+55235,int_stack+291552,int_stack+376519,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+202212,int_stack+59798,int_stack+59273, 1.0, int_stack+59063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81123,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+340101,int_stack+60554,int_stack+59798, 1.0, int_stack+59483, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81858,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5335,int_stack+340101,int_stack+202212, 1.0, int_stack+196362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284442,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+60995,int_stack+60554, 1.0, int_stack+60113, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82488,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+340101, 1.0, int_stack+202842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1890,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+5335, 1.0, int_stack+203787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290292,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+418819,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+61683,int_stack+61583, 1.0, int_stack+142204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163220, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+61833,int_stack+61683, 1.0, int_stack+142304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163320, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 1.0, int_stack+205047, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+274629, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+62043,int_stack+61833, 1.0, int_stack+142454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163470, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+376249,int_stack+375619,int_stack+291852, 1.0, int_stack+205347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+274929, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5335,int_stack+376249,int_stack+292302, 1.0, int_stack+205797, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275379, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+62473,int_stack+62323, 1.0, int_stack+143684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+164850, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+62698,int_stack+62473, 1.0, int_stack+143834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165225, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6335,int_stack+376069,int_stack+375619, 1.0, int_stack+206397, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275979, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376744,int_stack+63013,int_stack+62698, 1.0, int_stack+144059, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165765, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+291552,int_stack+376744,int_stack+376069, 1.0, int_stack+206847, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276429, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+291552,int_stack+6335, 1.0, int_stack+208927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+278509, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+202212,int_stack+375619,int_stack+5335,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5335,int_stack+64168,int_stack+63643, 1.0, int_stack+63433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81123, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5965,int_stack+64924,int_stack+64168, 1.0, int_stack+63853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81858, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+208927,int_stack+5965,int_stack+5335, 1.0, int_stack+214860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284442, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+291552,int_stack+65365,int_stack+64924, 1.0, int_stack+64483, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82488, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+291552,int_stack+5965, 1.0, int_stack+215490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1890, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+187232,int_stack+407719,int_stack+208927, 1.0, int_stack+221655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290292, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+58235,int_stack+187232,int_stack+375619,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+66053,int_stack+65953, 1.0, int_stack+147458, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163220, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+66203,int_stack+66053, 1.0, int_stack+147558, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163320, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 1.0, int_stack+222915, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+274629, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+66413,int_stack+66203, 1.0, int_stack+147708, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163470, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+187232,int_stack+376969,int_stack+375919, 1.0, int_stack+223215, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+274929, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+188132,int_stack+187232,int_stack+376369, 1.0, int_stack+223665, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275379, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+187232,int_stack+66843,int_stack+66693, 1.0, int_stack+148938, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+164850, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375619,int_stack+67068,int_stack+66843, 1.0, int_stack+149088, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165225, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376294,int_stack+375619,int_stack+187232, 1.0, int_stack+224265, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275979, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+221655,int_stack+67383,int_stack+67068, 1.0, int_stack+149313, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+165765, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+222600,int_stack+221655,int_stack+375619, 1.0, int_stack+224715, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276429, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+222600,int_stack+376294, 1.0, int_stack+226795, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+278509, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+221655,int_stack+291552,int_stack+188132,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+226795,int_stack+68538,int_stack+68013, 1.0, int_stack+67803, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81123, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375619,int_stack+69294,int_stack+68538, 1.0, int_stack+68223, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81858, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+187232,int_stack+375619,int_stack+226795, 1.0, int_stack+232728, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284442, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+232728,int_stack+69735,int_stack+69294, 1.0, int_stack+68853, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82488, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+232728,int_stack+375619, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1890, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+187232, 1.0, int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290292, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+62735,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+70423,int_stack+70323, 1.0, int_stack+152712, 0.0, zero_stack, 1.0, int_stack+163220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+70573,int_stack+70423, 1.0, int_stack+152812, 0.0, zero_stack, 1.0, int_stack+163320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 1.0, int_stack+239838, 0.0, zero_stack, 1.0, int_stack+274629, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+70783,int_stack+70573, 1.0, int_stack+152962, 0.0, zero_stack, 1.0, int_stack+163470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+376249,int_stack+375619,int_stack+291852, 1.0, int_stack+240138, 0.0, zero_stack, 1.0, int_stack+274929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+376249,int_stack+292302, 1.0, int_stack+240588, 0.0, zero_stack, 1.0, int_stack+275379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+71213,int_stack+71063, 1.0, int_stack+154192, 0.0, zero_stack, 1.0, int_stack+164850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+71438,int_stack+71213, 1.0, int_stack+154342, 0.0, zero_stack, 1.0, int_stack+165225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619, 1.0, int_stack+241188, 0.0, zero_stack, 1.0, int_stack+275979, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+71753,int_stack+71438, 1.0, int_stack+154567, 0.0, zero_stack, 1.0, int_stack+165765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+232728,int_stack+0,int_stack+376069, 1.0, int_stack+241638, 0.0, zero_stack, 1.0, int_stack+276429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+232728,int_stack+376744, 1.0, int_stack+243718, 0.0, zero_stack, 1.0, int_stack+278509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+239578,int_stack+291552,int_stack+238578,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+72908,int_stack+72383, 1.0, int_stack+72173, 0.0, zero_stack, 1.0, int_stack+81123, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+73664,int_stack+72908, 1.0, int_stack+72593, 0.0, zero_stack, 1.0, int_stack+81858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+243718,int_stack+0,int_stack+238578, 1.0, int_stack+249651, 0.0, zero_stack, 1.0, int_stack+284442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+232728,int_stack+74105,int_stack+73664, 1.0, int_stack+73223, 0.0, zero_stack, 1.0, int_stack+82488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+232728,int_stack+0, 1.0, int_stack+250281, 0.0, zero_stack, 1.0, int_stack+1890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+187232,int_stack+375619,int_stack+243718, 1.0, int_stack+256446, 0.0, zero_stack, 1.0, int_stack+290292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+67235,int_stack+187232,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+74793,int_stack+74693, 1.0, int_stack+157966, 1.0, int_stack+163220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+74943,int_stack+74793, 1.0, int_stack+158066, 1.0, int_stack+163320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 1.0, int_stack+257706, 1.0, int_stack+274629, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+187232,int_stack+75153,int_stack+74943, 1.0, int_stack+158216, 1.0, int_stack+163470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+187862,int_stack+187232,int_stack+291852, 1.0, int_stack+258006, 1.0, int_stack+274929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+187862,int_stack+292302, 1.0, int_stack+258456, 1.0, int_stack+275379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+187232,int_stack+75583,int_stack+75433, 1.0, int_stack+159446, 1.0, int_stack+164850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+187682,int_stack+75808,int_stack+75583, 1.0, int_stack+159596, 1.0, int_stack+165225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+188357,int_stack+187682,int_stack+187232, 1.0, int_stack+259056, 1.0, int_stack+275979, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+76123,int_stack+75808, 1.0, int_stack+159821, 1.0, int_stack+165765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+232728,int_stack+0,int_stack+187682, 1.0, int_stack+259506, 1.0, int_stack+276429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+232728,int_stack+188357, 1.0, int_stack+261586, 1.0, int_stack+278509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+256446,int_stack+291552,int_stack+238578,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+77278,int_stack+76753, 1.0, int_stack+76543, 1.0, int_stack+81123, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+78034,int_stack+77278, 1.0, int_stack+76963, 1.0, int_stack+81858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+261586,int_stack+0,int_stack+238578, 1.0, int_stack+267519, 1.0, int_stack+284442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+267519,int_stack+78475,int_stack+78034, 1.0, int_stack+77593, 1.0, int_stack+82488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+187232,int_stack+267519,int_stack+0, 1.0, int_stack+945, 1.0, int_stack+1890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+187232,int_stack+261586, 1.0, int_stack+273369, 1.0, int_stack+290292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+71735,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+291552,int_stack+79163,int_stack+79063, 2.0, int_stack+163220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+291852,int_stack+79313,int_stack+79163, 2.0, int_stack+163320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+292302,int_stack+291852,int_stack+291552, 2.0, int_stack+274629, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+79523,int_stack+79313, 2.0, int_stack+163470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+376249,int_stack+375619,int_stack+291852, 2.0, int_stack+274929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+376249,int_stack+292302, 2.0, int_stack+275379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+79953,int_stack+79803, 2.0, int_stack+164850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+80178,int_stack+79953, 2.0, int_stack+165225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619, 2.0, int_stack+275979, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+291552,int_stack+80493,int_stack+80178, 2.0, int_stack+165765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+291552,int_stack+376069, 2.0, int_stack+276429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+291552,int_stack+267519,int_stack+376744, 2.0, int_stack+278509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+273369,int_stack+291552,int_stack+238578,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+82173,int_stack+81333, 2.0, int_stack+81123, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+278509,int_stack+82929,int_stack+82173, 2.0, int_stack+81858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+267519,int_stack+278509,int_stack+238578, 2.0, int_stack+284442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+284442,int_stack+83370,int_stack+82929, 2.0, int_stack+82488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+284442,int_stack+278509, 2.0, int_stack+1890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+0,int_stack+267519, 2.0, int_stack+290292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+76235,int_stack+375619,int_stack+291552,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+84058,int_stack+83958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170130,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+84208,int_stack+84058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170230,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+84418,int_stack+84208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170380,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+290292,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192162,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+290292,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+290292,int_stack+84848,int_stack+84698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171610,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+290742,int_stack+85073,int_stack+84848, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171760,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+291417,int_stack+290742,int_stack+290292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+85388,int_stack+85073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171985,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+375619,int_stack+290742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195237,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+267519,int_stack+291417, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+80735,int_stack+375619,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+86018,int_stack+85808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102828,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+86333,int_stack+86018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103248,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+290292,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190232,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+291552,int_stack+86774,int_stack+86333, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103878,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+291552,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301365,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+187232,int_stack+0,int_stack+290292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302310,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+423319,int_stack+187232,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+197082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+87462,int_stack+87362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170130, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+87612,int_stack+87462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170230, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+87822,int_stack+87612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170380, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+187232,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192162, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+187232,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+187232,int_stack+88252,int_stack+88102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171610, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+187682,int_stack+88477,int_stack+88252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171760, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+188357,int_stack+187682,int_stack+187232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+88792,int_stack+88477, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171985, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+375619,int_stack+187682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195237, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+267519,int_stack+188357, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+83735,int_stack+375619,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+207927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+89422,int_stack+89212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102828, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+89737,int_stack+89422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103248, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+187232,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190232, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+290292,int_stack+90178,int_stack+89737, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103878, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+290292,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301365, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+0,int_stack+187232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302310, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+427819,int_stack+407719,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+90866,int_stack+90766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170130, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+91016,int_stack+90866, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170230, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+91226,int_stack+91016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170380, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192162, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+407719,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+91656,int_stack+91506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171610, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+408169,int_stack+91881,int_stack+91656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171760, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+408844,int_stack+408169,int_stack+407719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+92196,int_stack+91881, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171985, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+375619,int_stack+408169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195237, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+267519,int_stack+408844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+86735,int_stack+375619,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+92826,int_stack+92616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102828, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+93141,int_stack+92826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103248, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+407719,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190232, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+187232,int_stack+93582,int_stack+93141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103878, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+187232,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301365, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+187232,int_stack+0,int_stack+407719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302310, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+432319,int_stack+187232,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+94270,int_stack+94170, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+94420,int_stack+94270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+94630,int_stack+94420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+187232,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+192162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+187232,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+187232,int_stack+95060,int_stack+94910, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+187682,int_stack+95285,int_stack+95060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+188357,int_stack+187682,int_stack+187232, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+95600,int_stack+95285, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+375619,int_stack+187682, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+267519,int_stack+188357, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+89735,int_stack+375619,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+242718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+96230,int_stack+96020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+96545,int_stack+96230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+187232,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+96986,int_stack+96545, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+407719,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+0,int_stack+187232, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+92735,int_stack+407719,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+251946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+97674,int_stack+97574, 0.0, zero_stack, 1.0, int_stack+170130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+97824,int_stack+97674, 0.0, zero_stack, 1.0, int_stack+170230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+98034,int_stack+97824, 0.0, zero_stack, 1.0, int_stack+170380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+376969,int_stack+375919, 0.0, zero_stack, 1.0, int_stack+192162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+407719,int_stack+376369, 0.0, zero_stack, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+98464,int_stack+98314, 0.0, zero_stack, 1.0, int_stack+171610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+408169,int_stack+98689,int_stack+98464, 0.0, zero_stack, 1.0, int_stack+171760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+408844,int_stack+408169,int_stack+407719, 0.0, zero_stack, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+99004,int_stack+98689, 0.0, zero_stack, 1.0, int_stack+171985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+375619,int_stack+408169, 0.0, zero_stack, 1.0, int_stack+195237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+267519,int_stack+408844, 0.0, zero_stack, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+436819,int_stack+375619,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+99634,int_stack+99424, 0.0, zero_stack, 1.0, int_stack+102828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+99949,int_stack+99634, 0.0, zero_stack, 1.0, int_stack+103248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+407719,int_stack+267519,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+190232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+187232,int_stack+100390,int_stack+99949, 0.0, zero_stack, 1.0, int_stack+103878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+187232,int_stack+267519, 0.0, zero_stack, 1.0, int_stack+301365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+187232,int_stack+0,int_stack+407719, 0.0, zero_stack, 1.0, int_stack+302310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+439819,int_stack+187232,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+268869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+101078,int_stack+100978, 1.0, int_stack+170130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+101228,int_stack+101078, 1.0, int_stack+170230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 1.0, int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+101438,int_stack+101228, 1.0, int_stack+170380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+187232,int_stack+376969,int_stack+375919, 1.0, int_stack+192162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+187232,int_stack+376369, 1.0, int_stack+186632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+186632,int_stack+101868,int_stack+101718, 1.0, int_stack+171610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+187082,int_stack+102093,int_stack+101868, 1.0, int_stack+171760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+187757,int_stack+187082,int_stack+186632, 1.0, int_stack+194787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+191862,int_stack+102408,int_stack+102093, 1.0, int_stack+171985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+191862,int_stack+187082, 1.0, int_stack+195237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+191862,int_stack+267519,int_stack+187757, 1.0, int_stack+189332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+186632,int_stack+191862,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+277509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+103563,int_stack+103038, 1.0, int_stack+102828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+104319,int_stack+103563, 1.0, int_stack+103248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+193362,int_stack+267519,int_stack+238578, 1.0, int_stack+190232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+104760,int_stack+104319, 1.0, int_stack+103878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+267519, 1.0, int_stack+301365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+193362, 1.0, int_stack+302310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+97235,int_stack+375619,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+285792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+105448,int_stack+105348,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+192162,int_stack+105598,int_stack+105448,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+192612,int_stack+192162,int_stack+191862,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+193212,int_stack+105808,int_stack+105598,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+193842,int_stack+193212,int_stack+192162,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+193842,int_stack+192612,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+106238,int_stack+106088,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+192312,int_stack+106463,int_stack+106238,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+192987,int_stack+192312,int_stack+191862,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+193887,int_stack+106778,int_stack+106463,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+193887,int_stack+192312,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+193887,int_stack+267519,int_stack+192987,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+101735,int_stack+193887,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+107408,int_stack+107198,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+107723,int_stack+107408,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+191862,int_stack+267519,int_stack+238578,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+108164,int_stack+107723,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+267519,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+191862,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+444319,int_stack+375619,int_stack+193887, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+303660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+108852,int_stack+108752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175384,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+109002,int_stack+108852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175484,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308160,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+109212,int_stack+109002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175634,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+191862,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195912,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+191862,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308460,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+109642,int_stack+109492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+176864,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+192312,int_stack+109867,int_stack+109642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177014,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+192987,int_stack+192312,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309060,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+193887,int_stack+110182,int_stack+109867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177239,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+193887,int_stack+192312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309510,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+193887,int_stack+267519,int_stack+192987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+311590,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+104735,int_stack+193887,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+190862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+110812,int_stack+110602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127622,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+111127,int_stack+110812, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128042,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+191862,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+317973,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+111568,int_stack+111127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128672,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+312490,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+323823,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+448819,int_stack+375619,int_stack+193887, 0.0, zero_stack, 1.0, int_stack+197082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+112256,int_stack+112156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175384, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+112406,int_stack+112256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175484, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308160, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+112616,int_stack+112406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175634, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+191862,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195912, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+191862,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308460, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+113046,int_stack+112896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+176864, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+192312,int_stack+113271,int_stack+113046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177014, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+192987,int_stack+192312,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309060, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+193887,int_stack+113586,int_stack+113271, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177239, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+193887,int_stack+192312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309510, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+193887,int_stack+267519,int_stack+192987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+311590, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+107735,int_stack+193887,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+207927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+114216,int_stack+114006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127622, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+114531,int_stack+114216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128042, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+191862,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+317973, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+114972,int_stack+114531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128672, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+312490, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+323823, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+110735,int_stack+375619,int_stack+193887, 0.0, zero_stack, 1.0, int_stack+217155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+115660,int_stack+115560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175384, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+115810,int_stack+115660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175484, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308160, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+116020,int_stack+115810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175634, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+191862,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195912, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+191862,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308460, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+116450,int_stack+116300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+176864, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+192312,int_stack+116675,int_stack+116450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177014, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+192987,int_stack+192312,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309060, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+193887,int_stack+116990,int_stack+116675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177239, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+193887,int_stack+192312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309510, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+193887,int_stack+267519,int_stack+192987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+311590, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+453319,int_stack+193887,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+225795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+117620,int_stack+117410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127622, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+117935,int_stack+117620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128042, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+191862,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+317973, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+118376,int_stack+117935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128672, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+312490, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+323823, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+456319,int_stack+375619,int_stack+193887, 0.0, zero_stack, 1.0, int_stack+234078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+119064,int_stack+118964, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+119214,int_stack+119064, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+119424,int_stack+119214, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+175634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+191862,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+191862,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+119854,int_stack+119704, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+176864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+192312,int_stack+120079,int_stack+119854, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+192987,int_stack+192312,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+193887,int_stack+120394,int_stack+120079, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+193887,int_stack+192312, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+193887,int_stack+267519,int_stack+192987, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+311590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+115235,int_stack+193887,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+242718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+121024,int_stack+120814, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+121339,int_stack+121024, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+191862,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+317973, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+121780,int_stack+121339, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+312490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+191862, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+460819,int_stack+375619,int_stack+193887, 0.0, zero_stack, 1.0, int_stack+251946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+122468,int_stack+122368, 0.0, zero_stack, 1.0, int_stack+175384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+122618,int_stack+122468, 0.0, zero_stack, 1.0, int_stack+175484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+308160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+122828,int_stack+122618, 0.0, zero_stack, 1.0, int_stack+175634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+191862,int_stack+376969,int_stack+375919, 0.0, zero_stack, 1.0, int_stack+195912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+191862,int_stack+376369, 0.0, zero_stack, 1.0, int_stack+308460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+123258,int_stack+123108, 0.0, zero_stack, 1.0, int_stack+176864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+192312,int_stack+123483,int_stack+123258, 0.0, zero_stack, 1.0, int_stack+177014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+192987,int_stack+192312,int_stack+191862, 0.0, zero_stack, 1.0, int_stack+309060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+193887,int_stack+123798,int_stack+123483, 0.0, zero_stack, 1.0, int_stack+177239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+193887,int_stack+192312, 0.0, zero_stack, 1.0, int_stack+309510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+193887,int_stack+267519,int_stack+192987, 0.0, zero_stack, 1.0, int_stack+311590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+118235,int_stack+193887,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+260586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+124428,int_stack+124218, 0.0, zero_stack, 1.0, int_stack+127622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+124743,int_stack+124428, 0.0, zero_stack, 1.0, int_stack+128042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+191862,int_stack+267519,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+317973, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+125184,int_stack+124743, 0.0, zero_stack, 1.0, int_stack+128672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+267519, 0.0, zero_stack, 1.0, int_stack+312490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+191862, 0.0, zero_stack, 1.0, int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+121235,int_stack+375619,int_stack+193887, 0.0, zero_stack, 1.0, int_stack+268869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+125872,int_stack+125772, 1.0, int_stack+175384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+126022,int_stack+125872, 1.0, int_stack+175484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 1.0, int_stack+308160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+126232,int_stack+126022, 1.0, int_stack+175634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+191862,int_stack+376969,int_stack+375919, 1.0, int_stack+195912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+191862,int_stack+376369, 1.0, int_stack+308460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+191862,int_stack+126662,int_stack+126512, 1.0, int_stack+176864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+192312,int_stack+126887,int_stack+126662, 1.0, int_stack+177014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+308160,int_stack+192312,int_stack+191862, 1.0, int_stack+309060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+192987,int_stack+127202,int_stack+126887, 1.0, int_stack+177239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+192987,int_stack+192312, 1.0, int_stack+309510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+309060,int_stack+267519,int_stack+308160, 1.0, int_stack+311590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+191862,int_stack+309060,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+277509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+128357,int_stack+127832, 1.0, int_stack+127622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+129113,int_stack+128357, 1.0, int_stack+128042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+194862,int_stack+267519,int_stack+238578, 1.0, int_stack+317973, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+317973,int_stack+129554,int_stack+129113, 1.0, int_stack+128672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+317973,int_stack+267519, 1.0, int_stack+312490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+194862, 1.0, int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+465319,int_stack+407719,int_stack+309060, 0.0, zero_stack, 1.0, int_stack+285792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+130242,int_stack+130142,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+408019,int_stack+130392,int_stack+130242,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+408469,int_stack+408019,int_stack+407719,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+409069,int_stack+130602,int_stack+130392,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+323823,int_stack+409069,int_stack+408019,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+323823,int_stack+408469,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+323823,int_stack+131032,int_stack+130882,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+324273,int_stack+131257,int_stack+131032,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+407719,int_stack+324273,int_stack+323823,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+408619,int_stack+131572,int_stack+131257,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+408619,int_stack+324273,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+194862,int_stack+267519,int_stack+407719,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+125735,int_stack+194862,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+2835, 1.0, int_stack+310590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+132202,int_stack+131992,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+407719,int_stack+132517,int_stack+132202,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+407719,int_stack+238578,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+267519,int_stack+132958,int_stack+132517,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+267519,int_stack+407719,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+323823,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+128735,int_stack+407719,int_stack+194862, 0.0, zero_stack, 1.0, int_stack+303660, 1.0, int_stack+319323, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+194862,int_stack+133646,int_stack+133546,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+195162,int_stack+133796,int_stack+133646,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+195612,int_stack+195162,int_stack+194862,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+196212,int_stack+134006,int_stack+133796,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+196212,int_stack+195162,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+407719,int_stack+195612,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+134436,int_stack+134286,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+408169,int_stack+134661,int_stack+134436,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+408844,int_stack+408169,int_stack+407719,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+194862,int_stack+134976,int_stack+134661,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+194862,int_stack+408169,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+194862,int_stack+267519,int_stack+408844,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+469819,int_stack+194862,int_stack+238578, 0.0, zero_stack, 2.0, int_stack+310590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+135606,int_stack+135396,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+135921,int_stack+135606,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+267519,int_stack+238578,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+136362,int_stack+135921,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+267519,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+323823,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+472819,int_stack+407719,int_stack+194862, 0.0, zero_stack, 2.0, int_stack+319323, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+194862,int_stack+137790,int_stack+137690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181378,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195162,int_stack+137940,int_stack+137790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181478,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+195612,int_stack+195162,int_stack+194862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325083,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+196212,int_stack+138150,int_stack+137940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181628,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+196212,int_stack+195162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325383,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+407719,int_stack+195612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325833,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+139690,int_stack+139540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182858,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+408169,int_stack+139915,int_stack+139690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183008,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+408844,int_stack+408169,int_stack+407719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326433,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+194862,int_stack+140230,int_stack+139915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183233,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+194862,int_stack+408169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326883,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+194862,int_stack+267519,int_stack+408844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+328963,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+133235,int_stack+194862,int_stack+238578, 1.0, int_stack+190862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+140860,int_stack+140650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167610,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+141175,int_stack+140860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168030,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+141616,int_stack+141175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168660,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+337896,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+338841,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+136235,int_stack+407719,int_stack+194862, 1.0, int_stack+197082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+194862,int_stack+143044,int_stack+142944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181378, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195162,int_stack+143194,int_stack+143044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181478, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+195612,int_stack+195162,int_stack+194862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325083, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+196212,int_stack+143404,int_stack+143194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181628, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+196842,int_stack+196212,int_stack+195162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325383, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+196842,int_stack+195612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325833, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+194862,int_stack+144944,int_stack+144794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182858, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195312,int_stack+145169,int_stack+144944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183008, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+195987,int_stack+195312,int_stack+194862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326433, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+196887,int_stack+145484,int_stack+145169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183233, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+196887,int_stack+195312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326883, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+196887,int_stack+267519,int_stack+195987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+328963, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+140735,int_stack+196887,int_stack+238578, 1.0, int_stack+207927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+146114,int_stack+145904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167610, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+146429,int_stack+146114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168030, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+194862,int_stack+146870,int_stack+146429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168660, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+407719,int_stack+194862,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+337896, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+338841, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+205212,int_stack+375619,int_stack+196887, 1.0, int_stack+217155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+148298,int_stack+148198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181378, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375919,int_stack+148448,int_stack+148298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181478, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376369,int_stack+375919,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325083, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376969,int_stack+148658,int_stack+148448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181628, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+323823,int_stack+376969,int_stack+375919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325383, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+323823,int_stack+376369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325833, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+323823,int_stack+150198,int_stack+150048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182858, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+324273,int_stack+150423,int_stack+150198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183008, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+375619,int_stack+324273,int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326433, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+376519,int_stack+150738,int_stack+150423, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183233, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+376519,int_stack+324273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326883, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+267519,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+328963, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+194862,int_stack+407719,int_stack+238578, 1.0, int_stack+225795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+151368,int_stack+151158, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167610, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+375619,int_stack+151683,int_stack+151368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168030, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+375619,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+267519,int_stack+152124,int_stack+151683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168660, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+308160,int_stack+267519,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+337896, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+308160,int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+338841, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+143735,int_stack+375619,int_stack+407719, 1.0, int_stack+234078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+153552,int_stack+153452, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+408019,int_stack+153702,int_stack+153552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+408469,int_stack+408019,int_stack+407719, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325083, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+409069,int_stack+153912,int_stack+153702, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+409069,int_stack+408019, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+375619,int_stack+408469, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+325833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+155452,int_stack+155302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+155677,int_stack+155452, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+155992,int_stack+155677, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+183233, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+407719,int_stack+376069, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+326883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+267519,int_stack+376744, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+328963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+224655,int_stack+407719,int_stack+238578, 1.0, int_stack+242718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+156622,int_stack+156412, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+156937,int_stack+156622, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+267519,int_stack+238578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+157378,int_stack+156937, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+308160,int_stack+375619,int_stack+267519, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+337896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+308160,int_stack+323823, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+338841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+148235,int_stack+375619,int_stack+407719, 1.0, int_stack+251946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+158806,int_stack+158706, 0.0, zero_stack, 1.0, int_stack+181378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+408019,int_stack+158956,int_stack+158806, 0.0, zero_stack, 1.0, int_stack+181478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+408469,int_stack+408019,int_stack+407719, 0.0, zero_stack, 1.0, int_stack+325083, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+409069,int_stack+159166,int_stack+158956, 0.0, zero_stack, 1.0, int_stack+181628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+409069,int_stack+408019, 0.0, zero_stack, 1.0, int_stack+325383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+375619,int_stack+408469, 0.0, zero_stack, 1.0, int_stack+325833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+160706,int_stack+160556, 0.0, zero_stack, 1.0, int_stack+182858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+160931,int_stack+160706, 0.0, zero_stack, 1.0, int_stack+183008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619, 0.0, zero_stack, 1.0, int_stack+326433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+161246,int_stack+160931, 0.0, zero_stack, 1.0, int_stack+183233, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+267519,int_stack+407719,int_stack+376069, 0.0, zero_stack, 1.0, int_stack+326883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+267519,int_stack+376744, 0.0, zero_stack, 1.0, int_stack+328963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+249651,int_stack+407719,int_stack+238578, 1.0, int_stack+260586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+161876,int_stack+161666, 0.0, zero_stack, 1.0, int_stack+167610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+267519,int_stack+162191,int_stack+161876, 0.0, zero_stack, 1.0, int_stack+168030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+323823,int_stack+267519,int_stack+238578, 0.0, zero_stack, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+162632,int_stack+162191, 0.0, zero_stack, 1.0, int_stack+168660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+308160,int_stack+375619,int_stack+267519, 0.0, zero_stack, 1.0, int_stack+337896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+308160,int_stack+323823, 0.0, zero_stack, 1.0, int_stack+338841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+152735,int_stack+375619,int_stack+407719, 1.0, int_stack+268869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+164060,int_stack+163960, 1.0, int_stack+181378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+408019,int_stack+164210,int_stack+164060, 1.0, int_stack+181478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+408469,int_stack+408019,int_stack+407719, 1.0, int_stack+325083, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+409069,int_stack+164420,int_stack+164210, 1.0, int_stack+181628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+409069,int_stack+408019, 1.0, int_stack+325383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+375619,int_stack+408469, 1.0, int_stack+325833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+166650,int_stack+166500, 1.0, int_stack+182858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+166875,int_stack+166650, 1.0, int_stack+183008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619, 1.0, int_stack+326433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+167190,int_stack+166875, 1.0, int_stack+183233, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+317973,int_stack+407719,int_stack+376069, 1.0, int_stack+326883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+317973,int_stack+376744, 1.0, int_stack+328963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+323823,int_stack+407719,int_stack+238578, 1.0, int_stack+277509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+238578,int_stack+168345,int_stack+167820, 1.0, int_stack+167610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+328963,int_stack+169101,int_stack+168345, 1.0, int_stack+168030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+317973,int_stack+328963,int_stack+238578, 1.0, int_stack+201582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+169542,int_stack+169101, 1.0, int_stack+168660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+308160,int_stack+375619,int_stack+328963, 1.0, int_stack+337896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+308160,int_stack+317973, 1.0, int_stack+338841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+157235,int_stack+375619,int_stack+407719, 1.0, int_stack+285792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+170970,int_stack+170870,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+408019,int_stack+171120,int_stack+170970,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+408469,int_stack+408019,int_stack+407719,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+201582,int_stack+171330,int_stack+171120,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+201582,int_stack+408019,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+375619,int_stack+408469,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+172870,int_stack+172720,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+173095,int_stack+172870,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+173410,int_stack+173095,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+317973,int_stack+407719,int_stack+376069,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+317973,int_stack+376744,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+337896,int_stack+407719,int_stack+238578, 1.0, int_stack+2835, 0.0, zero_stack, 1.0, int_stack+327963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+201582,int_stack+174040,int_stack+173830,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+238578,int_stack+174355,int_stack+174040,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+317973,int_stack+238578,int_stack+201582,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+174796,int_stack+174355,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+308160,int_stack+375619,int_stack+238578,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+308160,int_stack+317973,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+161735,int_stack+375619,int_stack+407719, 1.0, int_stack+303660, 0.0, zero_stack, 1.0, int_stack+3835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+176224,int_stack+176124,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+408019,int_stack+176374,int_stack+176224,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+408469,int_stack+408019,int_stack+407719,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+201582,int_stack+176584,int_stack+176374,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+201582,int_stack+408019,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+375619,int_stack+408469,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+178124,int_stack+177974,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+178349,int_stack+178124,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+178664,int_stack+178349,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+317973,int_stack+407719,int_stack+376069,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+317973,int_stack+376744,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+301365,int_stack+407719,int_stack+238578, 1.0, int_stack+310590, 1.0, int_stack+327963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+201582,int_stack+179294,int_stack+179084,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+238578,int_stack+179609,int_stack+179294,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+317973,int_stack+238578,int_stack+201582,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+375619,int_stack+180050,int_stack+179609,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+308160,int_stack+375619,int_stack+238578,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+375619,int_stack+308160,int_stack+317973,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+308160,int_stack+375619,int_stack+407719, 1.0, int_stack+319323, 1.0, int_stack+3835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+407719,int_stack+182218,int_stack+182118,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+408019,int_stack+182368,int_stack+182218,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+408469,int_stack+408019,int_stack+407719,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+201582,int_stack+182578,int_stack+182368,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+201582,int_stack+408019,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+238578,int_stack+375619,int_stack+408469,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+375619,int_stack+184118,int_stack+183968,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+376069,int_stack+184343,int_stack+184118,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+376744,int_stack+376069,int_stack+375619,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+184658,int_stack+184343,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+317973,int_stack+407719,int_stack+376069,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+319323,int_stack+317973,int_stack+376744,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+319323,int_stack+238578, 2.0, int_stack+327963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+201582,int_stack+185288,int_stack+185078,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+238578,int_stack+185603,int_stack+185288,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+317973,int_stack+238578,int_stack+201582,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+407719,int_stack+186044,int_stack+185603,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+375619,int_stack+407719,int_stack+238578,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+407719,int_stack+375619,int_stack+317973,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+166235,int_stack+407719,int_stack+319323, 2.0, int_stack+3835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+170735,int_stack+210360,int_stack+198582,100);
     Libderiv->ABCD[11] = int_stack + 170735;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+209712,int_stack+228228,int_stack+218655,100);
     Libderiv->ABCD[10] = int_stack + 209712;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+227655,int_stack+245151,int_stack+235578,100);
     Libderiv->ABCD[9] = int_stack + 227655;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+242578,int_stack+263019,int_stack+253446,100);
     Libderiv->ABCD[8] = int_stack + 242578;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+259446,int_stack+279942,int_stack+270369,100);
     Libderiv->ABCD[7] = int_stack + 259446;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+276369,int_stack+296865,int_stack+287292,100);
     Libderiv->ABCD[6] = int_stack + 276369;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+176735,int_stack+313473,int_stack+305160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 176735;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+312660,int_stack+330396,int_stack+320823, 0.0, zero_stack, 1.0, int_stack+293052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 312660;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+326823,int_stack+345414,int_stack+334896, 1.0, int_stack+293052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 326823;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+290292,int_stack+349914,int_stack+342001,100);
     Libderiv->ABCD[155] = int_stack + 290292;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+340896,int_stack+354414,int_stack+7235,100);
     Libderiv->ABCD[143] = int_stack + 340896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+346896,int_stack+10235,int_stack+358914,100);
     Libderiv->ABCD[142] = int_stack + 346896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+352896,int_stack+366219,int_stack+14735,100);
     Libderiv->ABCD[131] = int_stack + 352896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+3000,int_stack+17735,int_stack+361914,100);
     Libderiv->ABCD[130] = int_stack + 3000;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+358896,int_stack+377719,int_stack+372619,100);
     Libderiv->ABCD[129] = int_stack + 358896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+364896,int_stack+382219,int_stack+22235,100);
     Libderiv->ABCD[119] = int_stack + 364896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+370896,int_stack+386719,int_stack+25235,100);
     Libderiv->ABCD[118] = int_stack + 370896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+376896,int_stack+31235,int_stack+28235,100);
     Libderiv->ABCD[117] = int_stack + 376896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+382896,int_stack+394219,int_stack+391219,100);
     Libderiv->ABCD[116] = int_stack + 382896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+388896,int_stack+398719,int_stack+35735,100);
     Libderiv->ABCD[107] = int_stack + 388896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+394896,int_stack+403219,int_stack+38735,100);
     Libderiv->ABCD[106] = int_stack + 394896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+400896,int_stack+44735,int_stack+41735,100);
     Libderiv->ABCD[105] = int_stack + 400896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+9000,int_stack+409819,int_stack+49235,100);
     Libderiv->ABCD[104] = int_stack + 9000;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+406896,int_stack+414319,int_stack+52235,100);
     Libderiv->ABCD[103] = int_stack + 406896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+15000,int_stack+418819,int_stack+55235,100);
     Libderiv->ABCD[95] = int_stack + 15000;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+412896,int_stack+58235,int_stack+202212,100);
     Libderiv->ABCD[94] = int_stack + 412896;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+21000,int_stack+62735,int_stack+221655,100);
     Libderiv->ABCD[93] = int_stack + 21000;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+27000,int_stack+67235,int_stack+239578,100);
     Libderiv->ABCD[92] = int_stack + 27000;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+33000,int_stack+71735,int_stack+256446,100);
     Libderiv->ABCD[91] = int_stack + 33000;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+39000,int_stack+76235,int_stack+273369,100);
     Libderiv->ABCD[90] = int_stack + 39000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+45000,int_stack+423319,int_stack+80735, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+198582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[47] = int_stack + 45000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+418896,int_stack+427819,int_stack+83735, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[46] = int_stack + 418896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+424896,int_stack+432319,int_stack+86735, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+235578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[45] = int_stack + 424896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+51000,int_stack+92735,int_stack+89735, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+253446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[44] = int_stack + 51000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+57000,int_stack+439819,int_stack+436819, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[43] = int_stack + 57000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+430896,int_stack+97235,int_stack+186632, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+287292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[42] = int_stack + 430896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+436896,int_stack+444319,int_stack+101735, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+305160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[38] = int_stack + 436896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+63000,int_stack+448819,int_stack+104735, 0.0, zero_stack, 1.0, int_stack+198582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[35] = int_stack + 63000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+442896,int_stack+110735,int_stack+107735, 0.0, zero_stack, 1.0, int_stack+218655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[34] = int_stack + 442896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+69000,int_stack+456319,int_stack+453319, 0.0, zero_stack, 1.0, int_stack+235578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[33] = int_stack + 69000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+448896,int_stack+460819,int_stack+115235, 0.0, zero_stack, 1.0, int_stack+253446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[32] = int_stack + 448896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+454896,int_stack+121235,int_stack+118235, 0.0, zero_stack, 1.0, int_stack+270369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[31] = int_stack + 454896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+75000,int_stack+465319,int_stack+191862, 0.0, zero_stack, 1.0, int_stack+287292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[30] = int_stack + 75000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+460896,int_stack+128735,int_stack+125735, 0.0, zero_stack, 1.0, int_stack+305160, 1.0, int_stack+320823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[26] = int_stack + 460896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+81000,int_stack+472819,int_stack+469819, 0.0, zero_stack, 2.0, int_stack+320823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[25] = int_stack + 81000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+466896,int_stack+136235,int_stack+133235, 1.0, int_stack+198582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[23] = int_stack + 466896;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+197862,int_stack+205212,int_stack+140735, 1.0, int_stack+218655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[22] = int_stack + 197862;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+87000,int_stack+143735,int_stack+194862, 1.0, int_stack+235578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[21] = int_stack + 87000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+93000,int_stack+148235,int_stack+224655, 1.0, int_stack+253446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[20] = int_stack + 93000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+99000,int_stack+152735,int_stack+249651, 1.0, int_stack+270369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[19] = int_stack + 99000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+105000,int_stack+157235,int_stack+323823, 1.0, int_stack+287292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[18] = int_stack + 105000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+111000,int_stack+161735,int_stack+337896, 1.0, int_stack+305160, 0.0, zero_stack, 1.0, int_stack+334896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[14] = int_stack + 111000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+117000,int_stack+308160,int_stack+301365, 1.0, int_stack+320823, 1.0, int_stack+334896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[13] = int_stack + 117000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+318660,int_stack+166235,int_stack+0, 2.0, int_stack+334896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[12] = int_stack + 318660;

}
