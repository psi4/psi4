#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ffff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (ff|ff) integrals */

void d12hrr_order_ffff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[6][6][11] = int_stack + 0;
 Libderiv->deriv_classes[6][6][10] = int_stack + 784;
 Libderiv->deriv_classes[6][6][9] = int_stack + 1568;
 Libderiv->deriv_classes[6][6][8] = int_stack + 2352;
 Libderiv->deriv_classes[6][6][7] = int_stack + 3136;
 Libderiv->dvrr_classes[6][5] = int_stack + 3920;
 Libderiv->deriv_classes[6][6][6] = int_stack + 4508;
 Libderiv->deriv_classes[6][6][2] = int_stack + 5292;
 Libderiv->deriv_classes[6][6][1] = int_stack + 6076;
 Libderiv->dvrr_classes[5][6] = int_stack + 6860;
 Libderiv->deriv_classes[6][6][0] = int_stack + 7448;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 8232;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 8332;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 8482;
 Libderiv->deriv2_classes[3][6][143] = int_stack + 8692;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 8972;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 9122;
 Libderiv->deriv2_classes[4][5][143] = int_stack + 9347;
 Libderiv->deriv2_classes[4][6][143] = int_stack + 9662;
 Libderiv->deriv2_classes[5][3][143] = int_stack + 10082;
 Libderiv->deriv2_classes[5][4][143] = int_stack + 10292;
 Libderiv->deriv2_classes[5][5][143] = int_stack + 10607;
 Libderiv->deriv2_classes[5][6][143] = int_stack + 11048;
 Libderiv->deriv2_classes[6][3][143] = int_stack + 11636;
 Libderiv->deriv2_classes[6][4][143] = int_stack + 11916;
 Libderiv->deriv2_classes[6][5][143] = int_stack + 12336;
 Libderiv->deriv2_classes[6][6][143] = int_stack + 12924;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 13708;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 13808;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 13958;
 Libderiv->deriv2_classes[3][6][131] = int_stack + 14168;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 14448;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 14598;
 Libderiv->deriv2_classes[4][5][131] = int_stack + 14823;
 Libderiv->deriv2_classes[4][6][131] = int_stack + 15138;
 Libderiv->deriv2_classes[5][3][131] = int_stack + 15558;
 Libderiv->deriv2_classes[5][4][131] = int_stack + 15768;
 Libderiv->deriv2_classes[5][5][131] = int_stack + 16083;
 Libderiv->deriv2_classes[5][6][131] = int_stack + 16524;
 Libderiv->deriv2_classes[6][3][131] = int_stack + 17112;
 Libderiv->deriv2_classes[6][4][131] = int_stack + 17392;
 Libderiv->deriv2_classes[6][5][131] = int_stack + 17812;
 Libderiv->deriv2_classes[6][6][131] = int_stack + 18400;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 19184;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 19284;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 19434;
 Libderiv->deriv2_classes[3][6][130] = int_stack + 19644;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 19924;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 20074;
 Libderiv->deriv2_classes[4][5][130] = int_stack + 20299;
 Libderiv->deriv2_classes[4][6][130] = int_stack + 20614;
 Libderiv->deriv2_classes[5][3][130] = int_stack + 21034;
 Libderiv->deriv2_classes[5][4][130] = int_stack + 21244;
 Libderiv->deriv2_classes[5][5][130] = int_stack + 21559;
 Libderiv->deriv2_classes[5][6][130] = int_stack + 22000;
 Libderiv->deriv2_classes[6][3][130] = int_stack + 22588;
 Libderiv->deriv2_classes[6][4][130] = int_stack + 22868;
 Libderiv->deriv2_classes[6][5][130] = int_stack + 23288;
 Libderiv->deriv2_classes[6][6][130] = int_stack + 23876;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 24660;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 24760;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 24910;
 Libderiv->deriv2_classes[3][6][119] = int_stack + 25120;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 25400;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 25550;
 Libderiv->deriv2_classes[4][5][119] = int_stack + 25775;
 Libderiv->deriv2_classes[4][6][119] = int_stack + 26090;
 Libderiv->deriv2_classes[5][3][119] = int_stack + 26510;
 Libderiv->deriv2_classes[5][4][119] = int_stack + 26720;
 Libderiv->deriv2_classes[5][5][119] = int_stack + 27035;
 Libderiv->deriv2_classes[5][6][119] = int_stack + 27476;
 Libderiv->deriv2_classes[6][3][119] = int_stack + 28064;
 Libderiv->deriv2_classes[6][4][119] = int_stack + 28344;
 Libderiv->deriv2_classes[6][5][119] = int_stack + 28764;
 Libderiv->deriv2_classes[6][6][119] = int_stack + 29352;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 30136;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 30236;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 30386;
 Libderiv->deriv2_classes[3][6][118] = int_stack + 30596;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 30876;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 31026;
 Libderiv->deriv2_classes[4][5][118] = int_stack + 31251;
 Libderiv->deriv2_classes[4][6][118] = int_stack + 31566;
 Libderiv->deriv2_classes[5][3][118] = int_stack + 31986;
 Libderiv->deriv2_classes[5][4][118] = int_stack + 32196;
 Libderiv->deriv2_classes[5][5][118] = int_stack + 32511;
 Libderiv->deriv2_classes[5][6][118] = int_stack + 32952;
 Libderiv->deriv2_classes[6][3][118] = int_stack + 33540;
 Libderiv->deriv2_classes[6][4][118] = int_stack + 33820;
 Libderiv->deriv2_classes[6][5][118] = int_stack + 34240;
 Libderiv->deriv2_classes[6][6][118] = int_stack + 34828;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 35612;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 35712;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 35862;
 Libderiv->deriv2_classes[3][6][117] = int_stack + 36072;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 36352;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 36502;
 Libderiv->deriv2_classes[4][5][117] = int_stack + 36727;
 Libderiv->deriv2_classes[4][6][117] = int_stack + 37042;
 Libderiv->deriv2_classes[5][3][117] = int_stack + 37462;
 Libderiv->deriv2_classes[5][4][117] = int_stack + 37672;
 Libderiv->deriv2_classes[5][5][117] = int_stack + 37987;
 Libderiv->deriv2_classes[5][6][117] = int_stack + 38428;
 Libderiv->deriv2_classes[6][3][117] = int_stack + 39016;
 Libderiv->deriv2_classes[6][4][117] = int_stack + 39296;
 Libderiv->deriv2_classes[6][5][117] = int_stack + 39716;
 Libderiv->deriv2_classes[6][6][117] = int_stack + 40304;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 41088;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 41188;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 41338;
 Libderiv->deriv2_classes[3][6][107] = int_stack + 41548;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 41828;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 41978;
 Libderiv->deriv2_classes[4][5][107] = int_stack + 42203;
 Libderiv->deriv2_classes[4][6][107] = int_stack + 42518;
 Libderiv->deriv2_classes[5][3][107] = int_stack + 42938;
 Libderiv->deriv2_classes[5][4][107] = int_stack + 43148;
 Libderiv->deriv2_classes[5][5][107] = int_stack + 43463;
 Libderiv->deriv2_classes[5][6][107] = int_stack + 43904;
 Libderiv->deriv2_classes[6][3][107] = int_stack + 44492;
 Libderiv->deriv2_classes[6][4][107] = int_stack + 44772;
 Libderiv->deriv2_classes[6][5][107] = int_stack + 45192;
 Libderiv->deriv2_classes[6][6][107] = int_stack + 45780;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 46564;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 46664;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 46814;
 Libderiv->deriv2_classes[3][6][106] = int_stack + 47024;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 47304;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 47454;
 Libderiv->deriv2_classes[4][5][106] = int_stack + 47679;
 Libderiv->deriv2_classes[4][6][106] = int_stack + 47994;
 Libderiv->deriv2_classes[5][3][106] = int_stack + 48414;
 Libderiv->deriv2_classes[5][4][106] = int_stack + 48624;
 Libderiv->deriv2_classes[5][5][106] = int_stack + 48939;
 Libderiv->deriv2_classes[5][6][106] = int_stack + 49380;
 Libderiv->deriv2_classes[6][3][106] = int_stack + 49968;
 Libderiv->deriv2_classes[6][4][106] = int_stack + 50248;
 Libderiv->deriv2_classes[6][5][106] = int_stack + 50668;
 Libderiv->deriv2_classes[6][6][106] = int_stack + 51256;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 52040;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 52140;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 52290;
 Libderiv->deriv2_classes[3][6][105] = int_stack + 52500;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 52780;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 52930;
 Libderiv->deriv2_classes[4][5][105] = int_stack + 53155;
 Libderiv->deriv2_classes[4][6][105] = int_stack + 53470;
 Libderiv->deriv2_classes[5][3][105] = int_stack + 53890;
 Libderiv->deriv2_classes[5][4][105] = int_stack + 54100;
 Libderiv->deriv2_classes[5][5][105] = int_stack + 54415;
 Libderiv->deriv2_classes[5][6][105] = int_stack + 54856;
 Libderiv->deriv2_classes[6][3][105] = int_stack + 55444;
 Libderiv->deriv2_classes[6][4][105] = int_stack + 55724;
 Libderiv->deriv2_classes[6][5][105] = int_stack + 56144;
 Libderiv->deriv2_classes[6][6][105] = int_stack + 56732;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 57516;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 57616;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 57766;
 Libderiv->deriv2_classes[3][6][104] = int_stack + 57976;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 58256;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 58406;
 Libderiv->deriv2_classes[4][5][104] = int_stack + 58631;
 Libderiv->deriv2_classes[4][6][104] = int_stack + 58946;
 Libderiv->deriv2_classes[5][3][104] = int_stack + 59366;
 Libderiv->deriv2_classes[5][4][104] = int_stack + 59576;
 Libderiv->deriv2_classes[5][5][104] = int_stack + 59891;
 Libderiv->deriv2_classes[5][6][104] = int_stack + 60332;
 Libderiv->deriv2_classes[6][3][104] = int_stack + 60920;
 Libderiv->deriv2_classes[6][4][104] = int_stack + 61200;
 Libderiv->deriv2_classes[6][5][104] = int_stack + 61620;
 Libderiv->deriv2_classes[6][6][104] = int_stack + 62208;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 62992;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 63092;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 63242;
 Libderiv->deriv2_classes[3][6][95] = int_stack + 63452;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 63732;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 63882;
 Libderiv->deriv2_classes[4][5][95] = int_stack + 64107;
 Libderiv->deriv2_classes[4][6][95] = int_stack + 64422;
 Libderiv->deriv2_classes[5][3][95] = int_stack + 64842;
 Libderiv->deriv2_classes[5][4][95] = int_stack + 65052;
 Libderiv->deriv2_classes[5][5][95] = int_stack + 65367;
 Libderiv->deriv2_classes[5][6][95] = int_stack + 65808;
 Libderiv->deriv2_classes[6][3][95] = int_stack + 66396;
 Libderiv->deriv2_classes[6][4][95] = int_stack + 66676;
 Libderiv->deriv2_classes[6][5][95] = int_stack + 67096;
 Libderiv->deriv2_classes[6][6][95] = int_stack + 67684;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 68468;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 68568;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 68718;
 Libderiv->deriv2_classes[3][6][94] = int_stack + 68928;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 69208;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 69358;
 Libderiv->deriv2_classes[4][5][94] = int_stack + 69583;
 Libderiv->deriv2_classes[4][6][94] = int_stack + 69898;
 Libderiv->deriv2_classes[5][3][94] = int_stack + 70318;
 Libderiv->deriv2_classes[5][4][94] = int_stack + 70528;
 Libderiv->deriv2_classes[5][5][94] = int_stack + 70843;
 Libderiv->deriv2_classes[5][6][94] = int_stack + 71284;
 Libderiv->deriv2_classes[6][3][94] = int_stack + 71872;
 Libderiv->deriv2_classes[6][4][94] = int_stack + 72152;
 Libderiv->deriv2_classes[6][5][94] = int_stack + 72572;
 Libderiv->deriv2_classes[6][6][94] = int_stack + 73160;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 73944;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 74044;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 74194;
 Libderiv->deriv2_classes[3][6][93] = int_stack + 74404;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 74684;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 74834;
 Libderiv->deriv2_classes[4][5][93] = int_stack + 75059;
 Libderiv->deriv2_classes[4][6][93] = int_stack + 75374;
 Libderiv->deriv2_classes[5][3][93] = int_stack + 75794;
 Libderiv->deriv2_classes[5][4][93] = int_stack + 76004;
 Libderiv->deriv2_classes[5][5][93] = int_stack + 76319;
 Libderiv->deriv2_classes[5][6][93] = int_stack + 76760;
 Libderiv->deriv2_classes[6][3][93] = int_stack + 77348;
 Libderiv->deriv2_classes[6][4][93] = int_stack + 77628;
 Libderiv->deriv2_classes[6][5][93] = int_stack + 78048;
 Libderiv->deriv2_classes[6][6][93] = int_stack + 78636;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 79420;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 79520;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 79670;
 Libderiv->deriv2_classes[3][6][92] = int_stack + 79880;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 80160;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 80310;
 Libderiv->deriv2_classes[4][5][92] = int_stack + 80535;
 Libderiv->deriv2_classes[4][6][92] = int_stack + 80850;
 Libderiv->deriv2_classes[5][3][92] = int_stack + 81270;
 Libderiv->deriv2_classes[5][4][92] = int_stack + 81480;
 Libderiv->deriv2_classes[5][5][92] = int_stack + 81795;
 Libderiv->deriv2_classes[5][6][92] = int_stack + 82236;
 Libderiv->deriv2_classes[6][3][92] = int_stack + 82824;
 Libderiv->deriv2_classes[6][4][92] = int_stack + 83104;
 Libderiv->deriv2_classes[6][5][92] = int_stack + 83524;
 Libderiv->deriv2_classes[6][6][92] = int_stack + 84112;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 84896;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 84996;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 85146;
 Libderiv->deriv2_classes[3][6][91] = int_stack + 85356;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 85636;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 85786;
 Libderiv->deriv2_classes[4][5][91] = int_stack + 86011;
 Libderiv->deriv2_classes[4][6][91] = int_stack + 86326;
 Libderiv->deriv2_classes[5][3][91] = int_stack + 86746;
 Libderiv->deriv2_classes[5][4][91] = int_stack + 86956;
 Libderiv->deriv2_classes[5][5][91] = int_stack + 87271;
 Libderiv->deriv2_classes[5][6][91] = int_stack + 87712;
 Libderiv->deriv2_classes[6][3][91] = int_stack + 88300;
 Libderiv->deriv2_classes[6][4][91] = int_stack + 88580;
 Libderiv->deriv2_classes[6][5][91] = int_stack + 89000;
 Libderiv->deriv2_classes[6][6][91] = int_stack + 89588;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 90372;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 90472;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 90622;
 Libderiv->deriv2_classes[3][6][83] = int_stack + 90832;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 91112;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 91262;
 Libderiv->deriv2_classes[4][5][83] = int_stack + 91487;
 Libderiv->deriv2_classes[4][6][83] = int_stack + 91802;
 Libderiv->deriv2_classes[5][3][83] = int_stack + 92222;
 Libderiv->deriv2_classes[5][4][83] = int_stack + 92432;
 Libderiv->deriv2_classes[5][5][83] = int_stack + 92747;
 Libderiv->deriv2_classes[5][6][83] = int_stack + 93188;
 Libderiv->deriv_classes[6][3][11] = int_stack + 93776;
 Libderiv->deriv2_classes[6][3][83] = int_stack + 94056;
 Libderiv->deriv_classes[6][4][11] = int_stack + 94336;
 Libderiv->deriv2_classes[6][4][83] = int_stack + 94756;
 Libderiv->deriv_classes[6][5][11] = int_stack + 95176;
 Libderiv->deriv2_classes[6][5][83] = int_stack + 95764;
 Libderiv->deriv2_classes[6][6][83] = int_stack + 96352;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 97136;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 97236;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 97386;
 Libderiv->deriv2_classes[3][6][82] = int_stack + 97596;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 97876;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 98026;
 Libderiv->deriv2_classes[4][5][82] = int_stack + 98251;
 Libderiv->deriv2_classes[4][6][82] = int_stack + 98566;
 Libderiv->deriv2_classes[5][3][82] = int_stack + 98986;
 Libderiv->deriv2_classes[5][4][82] = int_stack + 99196;
 Libderiv->deriv2_classes[5][5][82] = int_stack + 99511;
 Libderiv->deriv2_classes[5][6][82] = int_stack + 99952;
 Libderiv->deriv_classes[6][3][10] = int_stack + 100540;
 Libderiv->deriv2_classes[6][3][82] = int_stack + 100820;
 Libderiv->deriv_classes[6][4][10] = int_stack + 101100;
 Libderiv->deriv2_classes[6][4][82] = int_stack + 101520;
 Libderiv->deriv_classes[6][5][10] = int_stack + 101940;
 Libderiv->deriv2_classes[6][5][82] = int_stack + 102528;
 Libderiv->deriv2_classes[6][6][82] = int_stack + 103116;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 103900;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 104000;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 104150;
 Libderiv->deriv2_classes[3][6][81] = int_stack + 104360;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 104640;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 104790;
 Libderiv->deriv2_classes[4][5][81] = int_stack + 105015;
 Libderiv->deriv2_classes[4][6][81] = int_stack + 105330;
 Libderiv->deriv2_classes[5][3][81] = int_stack + 105750;
 Libderiv->deriv2_classes[5][4][81] = int_stack + 105960;
 Libderiv->deriv2_classes[5][5][81] = int_stack + 106275;
 Libderiv->deriv2_classes[5][6][81] = int_stack + 106716;
 Libderiv->deriv_classes[6][3][9] = int_stack + 107304;
 Libderiv->deriv2_classes[6][3][81] = int_stack + 107584;
 Libderiv->deriv_classes[6][4][9] = int_stack + 107864;
 Libderiv->deriv2_classes[6][4][81] = int_stack + 108284;
 Libderiv->deriv_classes[6][5][9] = int_stack + 108704;
 Libderiv->deriv2_classes[6][5][81] = int_stack + 109292;
 Libderiv->deriv2_classes[6][6][81] = int_stack + 109880;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 110664;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 110764;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 110914;
 Libderiv->deriv2_classes[3][6][80] = int_stack + 111124;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 111404;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 111554;
 Libderiv->deriv2_classes[4][5][80] = int_stack + 111779;
 Libderiv->deriv2_classes[4][6][80] = int_stack + 112094;
 Libderiv->deriv2_classes[5][3][80] = int_stack + 112514;
 Libderiv->deriv2_classes[5][4][80] = int_stack + 112724;
 Libderiv->deriv2_classes[5][5][80] = int_stack + 113039;
 Libderiv->deriv2_classes[5][6][80] = int_stack + 113480;
 Libderiv->deriv_classes[6][3][8] = int_stack + 114068;
 Libderiv->deriv2_classes[6][3][80] = int_stack + 114348;
 Libderiv->deriv_classes[6][4][8] = int_stack + 114628;
 Libderiv->deriv2_classes[6][4][80] = int_stack + 115048;
 Libderiv->deriv_classes[6][5][8] = int_stack + 115468;
 Libderiv->deriv2_classes[6][5][80] = int_stack + 116056;
 Libderiv->deriv2_classes[6][6][80] = int_stack + 116644;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 117428;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 117528;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 117678;
 Libderiv->deriv2_classes[3][6][79] = int_stack + 117888;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 118168;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 118318;
 Libderiv->deriv2_classes[4][5][79] = int_stack + 118543;
 Libderiv->deriv2_classes[4][6][79] = int_stack + 118858;
 Libderiv->deriv2_classes[5][3][79] = int_stack + 119278;
 Libderiv->deriv2_classes[5][4][79] = int_stack + 119488;
 Libderiv->deriv2_classes[5][5][79] = int_stack + 119803;
 Libderiv->deriv2_classes[5][6][79] = int_stack + 120244;
 Libderiv->deriv_classes[6][3][7] = int_stack + 120832;
 Libderiv->deriv2_classes[6][3][79] = int_stack + 121112;
 Libderiv->deriv_classes[6][4][7] = int_stack + 121392;
 Libderiv->deriv2_classes[6][4][79] = int_stack + 121812;
 Libderiv->deriv_classes[6][5][7] = int_stack + 122232;
 Libderiv->deriv2_classes[6][5][79] = int_stack + 122820;
 Libderiv->deriv2_classes[6][6][79] = int_stack + 123408;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 124192;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 124292;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 124442;
 Libderiv->deriv2_classes[3][6][78] = int_stack + 124652;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 124932;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 125082;
 Libderiv->deriv2_classes[4][5][78] = int_stack + 125307;
 Libderiv->deriv2_classes[4][6][78] = int_stack + 125622;
 Libderiv->deriv2_classes[5][3][78] = int_stack + 126042;
 Libderiv->deriv2_classes[5][4][78] = int_stack + 126252;
 Libderiv->deriv2_classes[5][5][78] = int_stack + 126567;
 Libderiv->deriv2_classes[5][6][78] = int_stack + 127008;
 Libderiv->dvrr_classes[6][3] = int_stack + 127596;
 Libderiv->deriv_classes[6][3][6] = int_stack + 127876;
 Libderiv->deriv2_classes[6][3][78] = int_stack + 128156;
 Libderiv->dvrr_classes[6][4] = int_stack + 128436;
 Libderiv->deriv_classes[6][4][6] = int_stack + 128856;
 Libderiv->deriv2_classes[6][4][78] = int_stack + 129276;
 Libderiv->deriv_classes[6][5][6] = int_stack + 129696;
 Libderiv->deriv2_classes[6][5][78] = int_stack + 130284;
 Libderiv->deriv2_classes[6][6][78] = int_stack + 130872;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 131656;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 131756;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 131906;
 Libderiv->deriv2_classes[3][6][35] = int_stack + 132116;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 132396;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 132546;
 Libderiv->deriv2_classes[4][5][35] = int_stack + 132771;
 Libderiv->deriv2_classes[4][6][35] = int_stack + 133086;
 Libderiv->deriv2_classes[5][3][35] = int_stack + 133506;
 Libderiv->deriv2_classes[5][4][35] = int_stack + 133716;
 Libderiv->deriv2_classes[5][5][35] = int_stack + 134031;
 Libderiv->deriv2_classes[5][6][35] = int_stack + 134472;
 Libderiv->deriv2_classes[6][3][35] = int_stack + 135060;
 Libderiv->deriv2_classes[6][4][35] = int_stack + 135340;
 Libderiv->deriv2_classes[6][5][35] = int_stack + 135760;
 Libderiv->deriv2_classes[6][6][35] = int_stack + 136348;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 137132;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 137232;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 137382;
 Libderiv->deriv2_classes[3][6][34] = int_stack + 137592;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 137872;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 138022;
 Libderiv->deriv2_classes[4][5][34] = int_stack + 138247;
 Libderiv->deriv2_classes[4][6][34] = int_stack + 138562;
 Libderiv->deriv2_classes[5][3][34] = int_stack + 138982;
 Libderiv->deriv2_classes[5][4][34] = int_stack + 139192;
 Libderiv->deriv2_classes[5][5][34] = int_stack + 139507;
 Libderiv->deriv2_classes[5][6][34] = int_stack + 139948;
 Libderiv->deriv2_classes[6][3][34] = int_stack + 140536;
 Libderiv->deriv2_classes[6][4][34] = int_stack + 140816;
 Libderiv->deriv2_classes[6][5][34] = int_stack + 141236;
 Libderiv->deriv2_classes[6][6][34] = int_stack + 141824;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 142608;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 142708;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 142858;
 Libderiv->deriv2_classes[3][6][33] = int_stack + 143068;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 143348;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 143498;
 Libderiv->deriv2_classes[4][5][33] = int_stack + 143723;
 Libderiv->deriv2_classes[4][6][33] = int_stack + 144038;
 Libderiv->deriv2_classes[5][3][33] = int_stack + 144458;
 Libderiv->deriv2_classes[5][4][33] = int_stack + 144668;
 Libderiv->deriv2_classes[5][5][33] = int_stack + 144983;
 Libderiv->deriv2_classes[5][6][33] = int_stack + 145424;
 Libderiv->deriv2_classes[6][3][33] = int_stack + 146012;
 Libderiv->deriv2_classes[6][4][33] = int_stack + 146292;
 Libderiv->deriv2_classes[6][5][33] = int_stack + 146712;
 Libderiv->deriv2_classes[6][6][33] = int_stack + 147300;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 148084;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 148184;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 148334;
 Libderiv->deriv2_classes[3][6][32] = int_stack + 148544;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 148824;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 148974;
 Libderiv->deriv2_classes[4][5][32] = int_stack + 149199;
 Libderiv->deriv2_classes[4][6][32] = int_stack + 149514;
 Libderiv->deriv2_classes[5][3][32] = int_stack + 149934;
 Libderiv->deriv2_classes[5][4][32] = int_stack + 150144;
 Libderiv->deriv2_classes[5][5][32] = int_stack + 150459;
 Libderiv->deriv2_classes[5][6][32] = int_stack + 150900;
 Libderiv->deriv2_classes[6][3][32] = int_stack + 151488;
 Libderiv->deriv2_classes[6][4][32] = int_stack + 151768;
 Libderiv->deriv2_classes[6][5][32] = int_stack + 152188;
 Libderiv->deriv2_classes[6][6][32] = int_stack + 152776;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 153560;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 153660;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 153810;
 Libderiv->deriv2_classes[3][6][31] = int_stack + 154020;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 154300;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 154450;
 Libderiv->deriv2_classes[4][5][31] = int_stack + 154675;
 Libderiv->deriv2_classes[4][6][31] = int_stack + 154990;
 Libderiv->deriv2_classes[5][3][31] = int_stack + 155410;
 Libderiv->deriv2_classes[5][4][31] = int_stack + 155620;
 Libderiv->deriv2_classes[5][5][31] = int_stack + 155935;
 Libderiv->deriv2_classes[5][6][31] = int_stack + 156376;
 Libderiv->deriv2_classes[6][3][31] = int_stack + 156964;
 Libderiv->deriv2_classes[6][4][31] = int_stack + 157244;
 Libderiv->deriv2_classes[6][5][31] = int_stack + 157664;
 Libderiv->deriv2_classes[6][6][31] = int_stack + 158252;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 159036;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 159136;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 159286;
 Libderiv->deriv2_classes[3][6][30] = int_stack + 159496;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 159776;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 159926;
 Libderiv->deriv2_classes[4][5][30] = int_stack + 160151;
 Libderiv->deriv2_classes[4][6][30] = int_stack + 160466;
 Libderiv->deriv2_classes[5][3][30] = int_stack + 160886;
 Libderiv->deriv2_classes[5][4][30] = int_stack + 161096;
 Libderiv->deriv2_classes[5][5][30] = int_stack + 161411;
 Libderiv->deriv2_classes[5][6][30] = int_stack + 161852;
 Libderiv->deriv_classes[6][3][2] = int_stack + 162440;
 Libderiv->deriv2_classes[6][3][30] = int_stack + 162720;
 Libderiv->deriv_classes[6][4][2] = int_stack + 163000;
 Libderiv->deriv2_classes[6][4][30] = int_stack + 163420;
 Libderiv->deriv_classes[6][5][2] = int_stack + 163840;
 Libderiv->deriv2_classes[6][5][30] = int_stack + 164428;
 Libderiv->deriv2_classes[6][6][30] = int_stack + 165016;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 165800;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 165900;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 166050;
 Libderiv->deriv2_classes[3][6][26] = int_stack + 166260;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 166540;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 166690;
 Libderiv->deriv2_classes[4][5][26] = int_stack + 166915;
 Libderiv->deriv2_classes[4][6][26] = int_stack + 167230;
 Libderiv->deriv2_classes[5][3][26] = int_stack + 167650;
 Libderiv->deriv2_classes[5][4][26] = int_stack + 167860;
 Libderiv->deriv2_classes[5][5][26] = int_stack + 168175;
 Libderiv->deriv2_classes[5][6][26] = int_stack + 168616;
 Libderiv->deriv2_classes[6][3][26] = int_stack + 169204;
 Libderiv->deriv2_classes[6][4][26] = int_stack + 169484;
 Libderiv->deriv2_classes[6][5][26] = int_stack + 169904;
 Libderiv->deriv2_classes[6][6][26] = int_stack + 170492;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 171276;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 171376;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 171526;
 Libderiv->deriv2_classes[3][6][23] = int_stack + 171736;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 172016;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 172166;
 Libderiv->deriv2_classes[4][5][23] = int_stack + 172391;
 Libderiv->deriv2_classes[4][6][23] = int_stack + 172706;
 Libderiv->deriv2_classes[5][3][23] = int_stack + 173126;
 Libderiv->deriv2_classes[5][4][23] = int_stack + 173336;
 Libderiv->deriv2_classes[5][5][23] = int_stack + 173651;
 Libderiv->deriv2_classes[5][6][23] = int_stack + 174092;
 Libderiv->deriv2_classes[6][3][23] = int_stack + 174680;
 Libderiv->deriv2_classes[6][4][23] = int_stack + 174960;
 Libderiv->deriv2_classes[6][5][23] = int_stack + 175380;
 Libderiv->deriv2_classes[6][6][23] = int_stack + 175968;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 176752;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 176852;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 177002;
 Libderiv->deriv2_classes[3][6][22] = int_stack + 177212;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 177492;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 177642;
 Libderiv->deriv2_classes[4][5][22] = int_stack + 177867;
 Libderiv->deriv2_classes[4][6][22] = int_stack + 178182;
 Libderiv->deriv2_classes[5][3][22] = int_stack + 178602;
 Libderiv->deriv2_classes[5][4][22] = int_stack + 178812;
 Libderiv->deriv2_classes[5][5][22] = int_stack + 179127;
 Libderiv->deriv2_classes[5][6][22] = int_stack + 179568;
 Libderiv->deriv2_classes[6][3][22] = int_stack + 180156;
 Libderiv->deriv2_classes[6][4][22] = int_stack + 180436;
 Libderiv->deriv2_classes[6][5][22] = int_stack + 180856;
 Libderiv->deriv2_classes[6][6][22] = int_stack + 181444;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 182228;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 182328;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 182478;
 Libderiv->deriv2_classes[3][6][21] = int_stack + 182688;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 182968;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 183118;
 Libderiv->deriv2_classes[4][5][21] = int_stack + 183343;
 Libderiv->deriv2_classes[4][6][21] = int_stack + 183658;
 Libderiv->deriv2_classes[5][3][21] = int_stack + 184078;
 Libderiv->deriv2_classes[5][4][21] = int_stack + 184288;
 Libderiv->deriv2_classes[5][5][21] = int_stack + 184603;
 Libderiv->deriv2_classes[5][6][21] = int_stack + 185044;
 Libderiv->deriv2_classes[6][3][21] = int_stack + 185632;
 Libderiv->deriv2_classes[6][4][21] = int_stack + 185912;
 Libderiv->deriv2_classes[6][5][21] = int_stack + 186332;
 Libderiv->deriv2_classes[6][6][21] = int_stack + 186920;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 187704;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 187804;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 187954;
 Libderiv->deriv2_classes[3][6][20] = int_stack + 188164;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 188444;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 188594;
 Libderiv->deriv2_classes[4][5][20] = int_stack + 188819;
 Libderiv->deriv2_classes[4][6][20] = int_stack + 189134;
 Libderiv->deriv2_classes[5][3][20] = int_stack + 189554;
 Libderiv->deriv2_classes[5][4][20] = int_stack + 189764;
 Libderiv->deriv2_classes[5][5][20] = int_stack + 190079;
 Libderiv->deriv2_classes[5][6][20] = int_stack + 190520;
 Libderiv->deriv2_classes[6][3][20] = int_stack + 191108;
 Libderiv->deriv2_classes[6][4][20] = int_stack + 191388;
 Libderiv->deriv2_classes[6][5][20] = int_stack + 191808;
 Libderiv->deriv2_classes[6][6][20] = int_stack + 192396;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 193180;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 193280;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 193430;
 Libderiv->deriv2_classes[3][6][19] = int_stack + 193640;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 193920;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 194070;
 Libderiv->deriv2_classes[4][5][19] = int_stack + 194295;
 Libderiv->deriv2_classes[4][6][19] = int_stack + 194610;
 Libderiv->deriv2_classes[5][3][19] = int_stack + 195030;
 Libderiv->deriv2_classes[5][4][19] = int_stack + 195240;
 Libderiv->deriv2_classes[5][5][19] = int_stack + 195555;
 Libderiv->deriv2_classes[5][6][19] = int_stack + 195996;
 Libderiv->deriv2_classes[6][3][19] = int_stack + 196584;
 Libderiv->deriv2_classes[6][4][19] = int_stack + 196864;
 Libderiv->deriv2_classes[6][5][19] = int_stack + 197284;
 Libderiv->deriv2_classes[6][6][19] = int_stack + 197872;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 198656;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 198756;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 198906;
 Libderiv->deriv2_classes[3][6][18] = int_stack + 199116;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 199396;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 199546;
 Libderiv->deriv2_classes[4][5][18] = int_stack + 199771;
 Libderiv->deriv2_classes[4][6][18] = int_stack + 200086;
 Libderiv->deriv2_classes[5][3][18] = int_stack + 200506;
 Libderiv->deriv2_classes[5][4][18] = int_stack + 200716;
 Libderiv->deriv2_classes[5][5][18] = int_stack + 201031;
 Libderiv->deriv2_classes[5][6][18] = int_stack + 201472;
 Libderiv->deriv_classes[6][3][1] = int_stack + 202060;
 Libderiv->deriv2_classes[6][3][18] = int_stack + 202340;
 Libderiv->deriv_classes[6][4][1] = int_stack + 202620;
 Libderiv->deriv2_classes[6][4][18] = int_stack + 203040;
 Libderiv->deriv_classes[6][5][1] = int_stack + 203460;
 Libderiv->deriv2_classes[6][5][18] = int_stack + 204048;
 Libderiv->deriv2_classes[6][6][18] = int_stack + 204636;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 205420;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 205520;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 205670;
 Libderiv->deriv2_classes[3][6][14] = int_stack + 205880;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 206160;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 206310;
 Libderiv->deriv2_classes[4][5][14] = int_stack + 206535;
 Libderiv->deriv2_classes[4][6][14] = int_stack + 206850;
 Libderiv->deriv2_classes[5][3][14] = int_stack + 207270;
 Libderiv->deriv2_classes[5][4][14] = int_stack + 207480;
 Libderiv->deriv2_classes[5][5][14] = int_stack + 207795;
 Libderiv->deriv2_classes[5][6][14] = int_stack + 208236;
 Libderiv->deriv2_classes[6][3][14] = int_stack + 208824;
 Libderiv->deriv2_classes[6][4][14] = int_stack + 209104;
 Libderiv->deriv2_classes[6][5][14] = int_stack + 209524;
 Libderiv->deriv2_classes[6][6][14] = int_stack + 210112;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 210896;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 210996;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 211146;
 Libderiv->deriv2_classes[3][6][13] = int_stack + 211356;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 211636;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 211786;
 Libderiv->deriv2_classes[4][5][13] = int_stack + 212011;
 Libderiv->deriv2_classes[4][6][13] = int_stack + 212326;
 Libderiv->deriv2_classes[5][3][13] = int_stack + 212746;
 Libderiv->deriv2_classes[5][4][13] = int_stack + 212956;
 Libderiv->deriv2_classes[5][5][13] = int_stack + 213271;
 Libderiv->deriv2_classes[5][6][13] = int_stack + 213712;
 Libderiv->deriv2_classes[6][3][13] = int_stack + 214300;
 Libderiv->deriv2_classes[6][4][13] = int_stack + 214580;
 Libderiv->deriv2_classes[6][5][13] = int_stack + 215000;
 Libderiv->deriv2_classes[6][6][13] = int_stack + 215588;
 Libderiv->deriv_classes[3][3][11] = int_stack + 216372;
 Libderiv->deriv_classes[3][4][11] = int_stack + 216472;
 Libderiv->deriv_classes[3][5][11] = int_stack + 216622;
 Libderiv->deriv_classes[3][6][11] = int_stack + 216832;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 217112;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 217212;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 217362;
 Libderiv->deriv2_classes[3][6][11] = int_stack + 217572;
 Libderiv->deriv_classes[4][3][11] = int_stack + 217852;
 Libderiv->deriv_classes[4][4][11] = int_stack + 218002;
 Libderiv->deriv_classes[4][5][11] = int_stack + 218227;
 Libderiv->deriv_classes[4][6][11] = int_stack + 218542;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 218962;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 219112;
 Libderiv->deriv2_classes[4][5][11] = int_stack + 219337;
 Libderiv->deriv2_classes[4][6][11] = int_stack + 219652;
 Libderiv->deriv_classes[5][3][11] = int_stack + 220072;
 Libderiv->deriv_classes[5][4][11] = int_stack + 220282;
 Libderiv->deriv_classes[5][5][11] = int_stack + 220597;
 Libderiv->deriv_classes[5][6][11] = int_stack + 221038;
 Libderiv->deriv2_classes[5][3][11] = int_stack + 221626;
 Libderiv->deriv2_classes[5][4][11] = int_stack + 221836;
 Libderiv->deriv2_classes[5][5][11] = int_stack + 222151;
 Libderiv->deriv2_classes[5][6][11] = int_stack + 222592;
 Libderiv->deriv2_classes[6][3][11] = int_stack + 223180;
 Libderiv->deriv2_classes[6][4][11] = int_stack + 223460;
 Libderiv->deriv2_classes[6][5][11] = int_stack + 223880;
 Libderiv->deriv2_classes[6][6][11] = int_stack + 224468;
 Libderiv->deriv_classes[3][3][10] = int_stack + 225252;
 Libderiv->deriv_classes[3][4][10] = int_stack + 225352;
 Libderiv->deriv_classes[3][5][10] = int_stack + 225502;
 Libderiv->deriv_classes[3][6][10] = int_stack + 225712;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 225992;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 226092;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 226242;
 Libderiv->deriv2_classes[3][6][10] = int_stack + 226452;
 Libderiv->deriv_classes[4][3][10] = int_stack + 226732;
 Libderiv->deriv_classes[4][4][10] = int_stack + 226882;
 Libderiv->deriv_classes[4][5][10] = int_stack + 227107;
 Libderiv->deriv_classes[4][6][10] = int_stack + 227422;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 227842;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 227992;
 Libderiv->deriv2_classes[4][5][10] = int_stack + 228217;
 Libderiv->deriv2_classes[4][6][10] = int_stack + 228532;
 Libderiv->deriv_classes[5][3][10] = int_stack + 228952;
 Libderiv->deriv_classes[5][4][10] = int_stack + 229162;
 Libderiv->deriv_classes[5][5][10] = int_stack + 229477;
 Libderiv->deriv_classes[5][6][10] = int_stack + 229918;
 Libderiv->deriv2_classes[5][3][10] = int_stack + 230506;
 Libderiv->deriv2_classes[5][4][10] = int_stack + 230716;
 Libderiv->deriv2_classes[5][5][10] = int_stack + 231031;
 Libderiv->deriv2_classes[5][6][10] = int_stack + 231472;
 Libderiv->deriv2_classes[6][3][10] = int_stack + 232060;
 Libderiv->deriv2_classes[6][4][10] = int_stack + 232340;
 Libderiv->deriv2_classes[6][5][10] = int_stack + 232760;
 Libderiv->deriv2_classes[6][6][10] = int_stack + 233348;
 Libderiv->deriv_classes[3][3][9] = int_stack + 234132;
 Libderiv->deriv_classes[3][4][9] = int_stack + 234232;
 Libderiv->deriv_classes[3][5][9] = int_stack + 234382;
 Libderiv->deriv_classes[3][6][9] = int_stack + 234592;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 234872;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 234972;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 235122;
 Libderiv->deriv2_classes[3][6][9] = int_stack + 235332;
 Libderiv->deriv_classes[4][3][9] = int_stack + 235612;
 Libderiv->deriv_classes[4][4][9] = int_stack + 235762;
 Libderiv->deriv_classes[4][5][9] = int_stack + 235987;
 Libderiv->deriv_classes[4][6][9] = int_stack + 236302;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 236722;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 236872;
 Libderiv->deriv2_classes[4][5][9] = int_stack + 237097;
 Libderiv->deriv2_classes[4][6][9] = int_stack + 237412;
 Libderiv->deriv_classes[5][3][9] = int_stack + 237832;
 Libderiv->deriv_classes[5][4][9] = int_stack + 238042;
 Libderiv->deriv_classes[5][5][9] = int_stack + 238357;
 Libderiv->deriv_classes[5][6][9] = int_stack + 238798;
 Libderiv->deriv2_classes[5][3][9] = int_stack + 239386;
 Libderiv->deriv2_classes[5][4][9] = int_stack + 239596;
 Libderiv->deriv2_classes[5][5][9] = int_stack + 239911;
 Libderiv->deriv2_classes[5][6][9] = int_stack + 240352;
 Libderiv->deriv2_classes[6][3][9] = int_stack + 240940;
 Libderiv->deriv2_classes[6][4][9] = int_stack + 241220;
 Libderiv->deriv2_classes[6][5][9] = int_stack + 241640;
 Libderiv->deriv2_classes[6][6][9] = int_stack + 242228;
 Libderiv->deriv_classes[3][3][8] = int_stack + 243012;
 Libderiv->deriv_classes[3][4][8] = int_stack + 243112;
 Libderiv->deriv_classes[3][5][8] = int_stack + 243262;
 Libderiv->deriv_classes[3][6][8] = int_stack + 243472;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 243752;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 243852;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 244002;
 Libderiv->deriv2_classes[3][6][8] = int_stack + 244212;
 Libderiv->deriv_classes[4][3][8] = int_stack + 244492;
 Libderiv->deriv_classes[4][4][8] = int_stack + 244642;
 Libderiv->deriv_classes[4][5][8] = int_stack + 244867;
 Libderiv->deriv_classes[4][6][8] = int_stack + 245182;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 245602;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 245752;
 Libderiv->deriv2_classes[4][5][8] = int_stack + 245977;
 Libderiv->deriv2_classes[4][6][8] = int_stack + 246292;
 Libderiv->deriv_classes[5][3][8] = int_stack + 246712;
 Libderiv->deriv_classes[5][4][8] = int_stack + 246922;
 Libderiv->deriv_classes[5][5][8] = int_stack + 247237;
 Libderiv->deriv_classes[5][6][8] = int_stack + 247678;
 Libderiv->deriv2_classes[5][3][8] = int_stack + 248266;
 Libderiv->deriv2_classes[5][4][8] = int_stack + 248476;
 Libderiv->deriv2_classes[5][5][8] = int_stack + 248791;
 Libderiv->deriv2_classes[5][6][8] = int_stack + 249232;
 Libderiv->deriv2_classes[6][3][8] = int_stack + 249820;
 Libderiv->deriv2_classes[6][4][8] = int_stack + 250100;
 Libderiv->deriv2_classes[6][5][8] = int_stack + 250520;
 Libderiv->deriv2_classes[6][6][8] = int_stack + 251108;
 Libderiv->deriv_classes[3][3][7] = int_stack + 251892;
 Libderiv->deriv_classes[3][4][7] = int_stack + 251992;
 Libderiv->deriv_classes[3][5][7] = int_stack + 252142;
 Libderiv->deriv_classes[3][6][7] = int_stack + 252352;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 252632;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 252732;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 252882;
 Libderiv->deriv2_classes[3][6][7] = int_stack + 253092;
 Libderiv->deriv_classes[4][3][7] = int_stack + 253372;
 Libderiv->deriv_classes[4][4][7] = int_stack + 253522;
 Libderiv->deriv_classes[4][5][7] = int_stack + 253747;
 Libderiv->deriv_classes[4][6][7] = int_stack + 254062;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 254482;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 254632;
 Libderiv->deriv2_classes[4][5][7] = int_stack + 254857;
 Libderiv->deriv2_classes[4][6][7] = int_stack + 255172;
 Libderiv->deriv_classes[5][3][7] = int_stack + 255592;
 Libderiv->deriv_classes[5][4][7] = int_stack + 255802;
 Libderiv->deriv_classes[5][5][7] = int_stack + 256117;
 Libderiv->deriv_classes[5][6][7] = int_stack + 256558;
 Libderiv->deriv2_classes[5][3][7] = int_stack + 257146;
 Libderiv->deriv2_classes[5][4][7] = int_stack + 257356;
 Libderiv->deriv2_classes[5][5][7] = int_stack + 257671;
 Libderiv->deriv2_classes[5][6][7] = int_stack + 258112;
 Libderiv->deriv2_classes[6][3][7] = int_stack + 258700;
 Libderiv->deriv2_classes[6][4][7] = int_stack + 258980;
 Libderiv->deriv2_classes[6][5][7] = int_stack + 259400;
 Libderiv->deriv2_classes[6][6][7] = int_stack + 259988;
 Libderiv->deriv_classes[3][3][6] = int_stack + 260772;
 Libderiv->deriv_classes[3][4][6] = int_stack + 260872;
 Libderiv->deriv_classes[3][5][6] = int_stack + 261022;
 Libderiv->deriv_classes[3][6][6] = int_stack + 261232;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 261512;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 261612;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 261762;
 Libderiv->deriv2_classes[3][6][6] = int_stack + 261972;
 Libderiv->deriv_classes[4][3][6] = int_stack + 262252;
 Libderiv->deriv_classes[4][4][6] = int_stack + 262402;
 Libderiv->deriv_classes[4][5][6] = int_stack + 262627;
 Libderiv->deriv_classes[4][6][6] = int_stack + 262942;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 263362;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 263512;
 Libderiv->deriv2_classes[4][5][6] = int_stack + 263737;
 Libderiv->deriv2_classes[4][6][6] = int_stack + 264052;
 Libderiv->dvrr_classes[5][3] = int_stack + 264472;
 Libderiv->deriv_classes[5][3][6] = int_stack + 264682;
 Libderiv->dvrr_classes[5][4] = int_stack + 264892;
 Libderiv->deriv_classes[5][4][6] = int_stack + 265207;
 Libderiv->dvrr_classes[5][5] = int_stack + 265522;
 Libderiv->deriv_classes[5][5][6] = int_stack + 265963;
 Libderiv->deriv_classes[5][6][6] = int_stack + 266404;
 Libderiv->deriv2_classes[5][3][6] = int_stack + 266992;
 Libderiv->deriv2_classes[5][4][6] = int_stack + 267202;
 Libderiv->deriv2_classes[5][5][6] = int_stack + 267517;
 Libderiv->deriv2_classes[5][6][6] = int_stack + 267958;
 Libderiv->deriv_classes[6][3][0] = int_stack + 268546;
 Libderiv->deriv2_classes[6][3][6] = int_stack + 268826;
 Libderiv->deriv_classes[6][4][0] = int_stack + 269106;
 Libderiv->deriv2_classes[6][4][6] = int_stack + 269526;
 Libderiv->deriv_classes[6][5][0] = int_stack + 269946;
 Libderiv->deriv2_classes[6][5][6] = int_stack + 270534;
 Libderiv->deriv2_classes[6][6][6] = int_stack + 271122;
 Libderiv->deriv_classes[3][3][2] = int_stack + 271906;
 Libderiv->deriv_classes[3][4][2] = int_stack + 272006;
 Libderiv->deriv_classes[3][5][2] = int_stack + 272156;
 Libderiv->deriv_classes[3][6][2] = int_stack + 272366;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 272646;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 272746;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 272896;
 Libderiv->deriv2_classes[3][6][2] = int_stack + 273106;
 Libderiv->deriv_classes[4][3][2] = int_stack + 273386;
 Libderiv->deriv_classes[4][4][2] = int_stack + 273536;
 Libderiv->deriv_classes[4][5][2] = int_stack + 273761;
 Libderiv->deriv_classes[4][6][2] = int_stack + 274076;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 274496;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 274646;
 Libderiv->deriv2_classes[4][5][2] = int_stack + 274871;
 Libderiv->deriv2_classes[4][6][2] = int_stack + 275186;
 Libderiv->deriv_classes[5][3][2] = int_stack + 275606;
 Libderiv->deriv_classes[5][4][2] = int_stack + 275816;
 Libderiv->deriv_classes[5][5][2] = int_stack + 276131;
 Libderiv->deriv_classes[5][6][2] = int_stack + 276572;
 Libderiv->deriv2_classes[5][3][2] = int_stack + 277160;
 Libderiv->deriv2_classes[5][4][2] = int_stack + 277370;
 Libderiv->deriv2_classes[5][5][2] = int_stack + 277685;
 Libderiv->deriv2_classes[5][6][2] = int_stack + 278126;
 Libderiv->deriv2_classes[6][3][2] = int_stack + 278714;
 Libderiv->deriv2_classes[6][4][2] = int_stack + 278994;
 Libderiv->deriv2_classes[6][5][2] = int_stack + 279414;
 Libderiv->deriv2_classes[6][6][2] = int_stack + 280002;
 Libderiv->deriv_classes[3][3][1] = int_stack + 280786;
 Libderiv->deriv_classes[3][4][1] = int_stack + 280886;
 Libderiv->deriv_classes[3][5][1] = int_stack + 281036;
 Libderiv->deriv_classes[3][6][1] = int_stack + 281246;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 281526;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 281626;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 281776;
 Libderiv->deriv2_classes[3][6][1] = int_stack + 281986;
 Libderiv->deriv_classes[4][3][1] = int_stack + 282266;
 Libderiv->deriv_classes[4][4][1] = int_stack + 282416;
 Libderiv->deriv_classes[4][5][1] = int_stack + 282641;
 Libderiv->deriv_classes[4][6][1] = int_stack + 282956;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 283376;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 283526;
 Libderiv->deriv2_classes[4][5][1] = int_stack + 283751;
 Libderiv->deriv2_classes[4][6][1] = int_stack + 284066;
 Libderiv->deriv_classes[5][3][1] = int_stack + 284486;
 Libderiv->deriv_classes[5][4][1] = int_stack + 284696;
 Libderiv->deriv_classes[5][5][1] = int_stack + 285011;
 Libderiv->deriv_classes[5][6][1] = int_stack + 285452;
 Libderiv->deriv2_classes[5][3][1] = int_stack + 286040;
 Libderiv->deriv2_classes[5][4][1] = int_stack + 286250;
 Libderiv->deriv2_classes[5][5][1] = int_stack + 286565;
 Libderiv->deriv2_classes[5][6][1] = int_stack + 287006;
 Libderiv->deriv2_classes[6][3][1] = int_stack + 287594;
 Libderiv->deriv2_classes[6][4][1] = int_stack + 287874;
 Libderiv->deriv2_classes[6][5][1] = int_stack + 288294;
 Libderiv->deriv2_classes[6][6][1] = int_stack + 288882;
 Libderiv->dvrr_classes[3][3] = int_stack + 289666;
 Libderiv->dvrr_classes[3][4] = int_stack + 289766;
 Libderiv->dvrr_classes[3][5] = int_stack + 289916;
 Libderiv->dvrr_classes[3][6] = int_stack + 290126;
 Libderiv->deriv_classes[3][3][0] = int_stack + 290406;
 Libderiv->deriv_classes[3][4][0] = int_stack + 290506;
 Libderiv->deriv_classes[3][5][0] = int_stack + 290656;
 Libderiv->deriv_classes[3][6][0] = int_stack + 290866;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 291146;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 291246;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 291396;
 Libderiv->deriv2_classes[3][6][0] = int_stack + 291606;
 Libderiv->dvrr_classes[4][3] = int_stack + 291886;
 Libderiv->dvrr_classes[4][4] = int_stack + 292036;
 Libderiv->dvrr_classes[4][5] = int_stack + 292261;
 Libderiv->dvrr_classes[4][6] = int_stack + 292576;
 Libderiv->deriv_classes[4][3][0] = int_stack + 292996;
 Libderiv->deriv_classes[4][4][0] = int_stack + 293146;
 Libderiv->deriv_classes[4][5][0] = int_stack + 293371;
 Libderiv->deriv_classes[4][6][0] = int_stack + 293686;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 294106;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 294256;
 Libderiv->deriv2_classes[4][5][0] = int_stack + 294481;
 Libderiv->deriv2_classes[4][6][0] = int_stack + 294796;
 Libderiv->deriv_classes[5][3][0] = int_stack + 295216;
 Libderiv->deriv_classes[5][4][0] = int_stack + 295426;
 Libderiv->deriv_classes[5][5][0] = int_stack + 295741;
 Libderiv->deriv_classes[5][6][0] = int_stack + 296182;
 Libderiv->deriv2_classes[5][3][0] = int_stack + 296770;
 Libderiv->deriv2_classes[5][4][0] = int_stack + 296980;
 Libderiv->deriv2_classes[5][5][0] = int_stack + 297295;
 Libderiv->deriv2_classes[5][6][0] = int_stack + 297736;
 Libderiv->deriv2_classes[6][3][0] = int_stack + 298324;
 Libderiv->deriv2_classes[6][4][0] = int_stack + 298604;
 Libderiv->deriv2_classes[6][5][0] = int_stack + 299024;
 Libderiv->deriv2_classes[6][6][0] = int_stack + 299612;
 memset(int_stack,0,2403168);

 Libderiv->dvrr_stack = int_stack + 1013746;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ffff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+300396,int_stack+289766,int_stack+289666,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+300696,int_stack+289916,int_stack+289766,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+301146,int_stack+300696,int_stack+300396,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+301746,int_stack+216472,int_stack+216372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289666,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+302046,int_stack+216622,int_stack+216472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289766,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+302496,int_stack+302046,int_stack+301746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+303096,int_stack+216832,int_stack+216622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289916,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+303726,int_stack+303096,int_stack+302046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300696,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+304626,int_stack+303726,int_stack+302496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301146,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+303096,int_stack+292036,int_stack+291886,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+303546,int_stack+292261,int_stack+292036,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+305626,int_stack+303546,int_stack+303096,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+306526,int_stack+218002,int_stack+217852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+291886,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+306976,int_stack+218227,int_stack+218002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292036,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+307651,int_stack+306976,int_stack+306526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+308551,int_stack+218542,int_stack+218227, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292261,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+309496,int_stack+308551,int_stack+306976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303546,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+310846,int_stack+309496,int_stack+307651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+312346,int_stack+310846,int_stack+304626,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+308551,int_stack+264892,int_stack+264472,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+309181,int_stack+265522,int_stack+264892,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+315346,int_stack+309181,int_stack+308551,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+310126,int_stack+220282,int_stack+220072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264472,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+316606,int_stack+220597,int_stack+220282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264892,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+317551,int_stack+316606,int_stack+310126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+318811,int_stack+221038,int_stack+220597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265522,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+320134,int_stack+318811,int_stack+316606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309181,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+322024,int_stack+320134,int_stack+317551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+324124,int_stack+322024,int_stack+310846,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+328624,int_stack+324124,int_stack+312346,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+318811,int_stack+128436,int_stack+127596,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+319651,int_stack+3920,int_stack+128436,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+334624,int_stack+319651,int_stack+318811,28);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+320911,int_stack+94336,int_stack+93776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127596,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+336304,int_stack+95176,int_stack+94336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128436,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+337564,int_stack+336304,int_stack+320911, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+318811,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+339244,int_stack+0,int_stack+95176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3920,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+341008,int_stack+339244,int_stack+336304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+319651,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+343528,int_stack+341008,int_stack+337564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+346328,int_stack+343528,int_stack+322024,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+352628,int_stack+346328,int_stack+324124,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+339244,int_stack+225352,int_stack+225252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289666, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+339544,int_stack+225502,int_stack+225352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289766, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+339994,int_stack+339544,int_stack+339244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+340594,int_stack+225712,int_stack+225502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289916, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+341224,int_stack+340594,int_stack+339544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300696, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+342124,int_stack+341224,int_stack+339994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301146, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+340594,int_stack+226882,int_stack+226732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+291886, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+341044,int_stack+227107,int_stack+226882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292036, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+343124,int_stack+341044,int_stack+340594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+344024,int_stack+227422,int_stack+227107, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292261, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+344969,int_stack+344024,int_stack+341044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303546, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+346319,int_stack+344969,int_stack+343124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+347819,int_stack+346319,int_stack+342124,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+344024,int_stack+229162,int_stack+228952, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264472, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+344654,int_stack+229477,int_stack+229162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264892, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+350819,int_stack+344654,int_stack+344024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+361628,int_stack+229918,int_stack+229477, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265522, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+362951,int_stack+361628,int_stack+344654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309181, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+364841,int_stack+362951,int_stack+350819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+366941,int_stack+364841,int_stack+346319,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+371441,int_stack+366941,int_stack+347819,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+361628,int_stack+101100,int_stack+100540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127596, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+362468,int_stack+101940,int_stack+101100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128436, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+377441,int_stack+362468,int_stack+361628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+318811, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+379121,int_stack+784,int_stack+101940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3920, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+380885,int_stack+379121,int_stack+362468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+319651, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+383405,int_stack+380885,int_stack+377441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+386205,int_stack+383405,int_stack+364841,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+392505,int_stack+386205,int_stack+366941,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+379121,int_stack+234232,int_stack+234132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289666, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+379421,int_stack+234382,int_stack+234232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289766, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+379871,int_stack+379421,int_stack+379121, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+380471,int_stack+234592,int_stack+234382, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289916, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+381101,int_stack+380471,int_stack+379421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300696, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+382001,int_stack+381101,int_stack+379871, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301146, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+380471,int_stack+235762,int_stack+235612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+291886, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+380921,int_stack+235987,int_stack+235762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292036, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+383001,int_stack+380921,int_stack+380471, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+383901,int_stack+236302,int_stack+235987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292261, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+384846,int_stack+383901,int_stack+380921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303546, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+386196,int_stack+384846,int_stack+383001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+387696,int_stack+386196,int_stack+382001,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+383901,int_stack+238042,int_stack+237832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264472, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+384531,int_stack+238357,int_stack+238042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264892, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+390696,int_stack+384531,int_stack+383901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+238798,int_stack+238357, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265522, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+0,int_stack+384531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309181, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+403395,int_stack+401505,int_stack+390696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+405495,int_stack+403395,int_stack+386196,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+409995,int_stack+405495,int_stack+387696,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+401505,int_stack+107864,int_stack+107304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127596, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+108704,int_stack+107864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128436, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+415995,int_stack+0,int_stack+401505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+318811, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+417675,int_stack+1568,int_stack+108704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3920, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+419439,int_stack+417675,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+319651, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+421959,int_stack+419439,int_stack+415995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+424759,int_stack+421959,int_stack+403395,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+431059,int_stack+424759,int_stack+405495,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+417675,int_stack+243112,int_stack+243012, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+417975,int_stack+243262,int_stack+243112, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+418425,int_stack+417975,int_stack+417675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+419025,int_stack+243472,int_stack+243262, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+289916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+419655,int_stack+419025,int_stack+417975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+420555,int_stack+419655,int_stack+418425, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+419025,int_stack+244642,int_stack+244492, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+291886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+419475,int_stack+244867,int_stack+244642, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+421555,int_stack+419475,int_stack+419025, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+422455,int_stack+245182,int_stack+244867, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+423400,int_stack+422455,int_stack+419475, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+424750,int_stack+423400,int_stack+421555, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+426250,int_stack+424750,int_stack+420555,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+422455,int_stack+246922,int_stack+246712, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+423085,int_stack+247237,int_stack+246922, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+429250,int_stack+423085,int_stack+422455, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+440059,int_stack+247678,int_stack+247237, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+441382,int_stack+440059,int_stack+423085, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309181, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+443272,int_stack+441382,int_stack+429250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+445372,int_stack+443272,int_stack+424750,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+449872,int_stack+445372,int_stack+426250,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+440059,int_stack+114628,int_stack+114068, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+440899,int_stack+115468,int_stack+114628, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+455872,int_stack+440899,int_stack+440059, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+457552,int_stack+2352,int_stack+115468, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+459316,int_stack+457552,int_stack+440899, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+319651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+461836,int_stack+459316,int_stack+455872, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+464636,int_stack+461836,int_stack+443272,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+470936,int_stack+464636,int_stack+445372,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+457552,int_stack+251992,int_stack+251892, 0.0, zero_stack, 1.0, int_stack+289666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+457852,int_stack+252142,int_stack+251992, 0.0, zero_stack, 1.0, int_stack+289766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+458302,int_stack+457852,int_stack+457552, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+458902,int_stack+252352,int_stack+252142, 0.0, zero_stack, 1.0, int_stack+289916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+459532,int_stack+458902,int_stack+457852, 0.0, zero_stack, 1.0, int_stack+300696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+460432,int_stack+459532,int_stack+458302, 0.0, zero_stack, 1.0, int_stack+301146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+458902,int_stack+253522,int_stack+253372, 0.0, zero_stack, 1.0, int_stack+291886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+459352,int_stack+253747,int_stack+253522, 0.0, zero_stack, 1.0, int_stack+292036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+461432,int_stack+459352,int_stack+458902, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+462332,int_stack+254062,int_stack+253747, 0.0, zero_stack, 1.0, int_stack+292261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+463277,int_stack+462332,int_stack+459352, 0.0, zero_stack, 1.0, int_stack+303546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+464627,int_stack+463277,int_stack+461432, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+466127,int_stack+464627,int_stack+460432,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+462332,int_stack+255802,int_stack+255592, 0.0, zero_stack, 1.0, int_stack+264472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+462962,int_stack+256117,int_stack+255802, 0.0, zero_stack, 1.0, int_stack+264892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+469127,int_stack+462962,int_stack+462332, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1260,int_stack+256558,int_stack+256117, 0.0, zero_stack, 1.0, int_stack+265522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+479936,int_stack+1260,int_stack+462962, 0.0, zero_stack, 1.0, int_stack+309181, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+481826,int_stack+479936,int_stack+469127, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+483926,int_stack+481826,int_stack+464627,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+488426,int_stack+483926,int_stack+466127,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+479936,int_stack+121392,int_stack+120832, 0.0, zero_stack, 1.0, int_stack+127596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+122232,int_stack+121392, 0.0, zero_stack, 1.0, int_stack+128436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+494426,int_stack+1260,int_stack+479936, 0.0, zero_stack, 1.0, int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+496106,int_stack+3136,int_stack+122232, 0.0, zero_stack, 1.0, int_stack+3920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+497870,int_stack+496106,int_stack+1260, 0.0, zero_stack, 1.0, int_stack+319651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+500390,int_stack+497870,int_stack+494426, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+503190,int_stack+500390,int_stack+481826,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+509490,int_stack+503190,int_stack+483926,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+496106,int_stack+260872,int_stack+260772, 1.0, int_stack+289666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+496406,int_stack+261022,int_stack+260872, 1.0, int_stack+289766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+496856,int_stack+496406,int_stack+496106, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+497456,int_stack+261232,int_stack+261022, 1.0, int_stack+289916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+498086,int_stack+497456,int_stack+496406, 1.0, int_stack+300696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+498986,int_stack+498086,int_stack+496856, 1.0, int_stack+301146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+497456,int_stack+262402,int_stack+262252, 1.0, int_stack+291886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+497906,int_stack+262627,int_stack+262402, 1.0, int_stack+292036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+499986,int_stack+497906,int_stack+497456, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+500886,int_stack+262942,int_stack+262627, 1.0, int_stack+292261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+501831,int_stack+500886,int_stack+497906, 1.0, int_stack+303546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+503181,int_stack+501831,int_stack+499986, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+504681,int_stack+503181,int_stack+498986,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+500886,int_stack+265207,int_stack+264682, 1.0, int_stack+264472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+501516,int_stack+265963,int_stack+265207, 1.0, int_stack+264892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+507681,int_stack+501516,int_stack+500886, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2520,int_stack+266404,int_stack+265963, 1.0, int_stack+265522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+518490,int_stack+2520,int_stack+501516, 1.0, int_stack+309181, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+520380,int_stack+518490,int_stack+507681, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+522480,int_stack+520380,int_stack+503181,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+526980,int_stack+522480,int_stack+504681,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+518490,int_stack+128856,int_stack+127876, 1.0, int_stack+127596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2520,int_stack+129696,int_stack+128856, 1.0, int_stack+128436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+2520,int_stack+518490, 1.0, int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+534660,int_stack+4508,int_stack+129696, 1.0, int_stack+3920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+536424,int_stack+534660,int_stack+2520, 1.0, int_stack+319651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+538944,int_stack+536424,int_stack+532980, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+541744,int_stack+538944,int_stack+520380,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+548044,int_stack+541744,int_stack+522480,100);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+308551,int_stack+290126,int_stack+289916,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+334624,int_stack+308551,int_stack+300696,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+534660,int_stack+334624,int_stack+301146,10);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+334624,int_stack+292576,int_stack+292261,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+300396,int_stack+334624,int_stack+303546,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+334624,int_stack+300396,int_stack+305626,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+535660,int_stack+334624,int_stack+534660,100);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+300396,int_stack+6860,int_stack+265522,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+538660,int_stack+300396,int_stack+309181,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+538660,int_stack+315346,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+538660,int_stack+318811,int_stack+334624,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+557044,int_stack+538660,int_stack+535660,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+315346,int_stack+272006,int_stack+271906,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+315646,int_stack+272156,int_stack+272006,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+543160,int_stack+315646,int_stack+315346,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+543760,int_stack+272366,int_stack+272156,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+305626,int_stack+543760,int_stack+315646,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+543760,int_stack+305626,int_stack+543160,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+273536,int_stack+273386,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+544760,int_stack+273761,int_stack+273536,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+545435,int_stack+544760,int_stack+305626,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+546335,int_stack+274076,int_stack+273761,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+300396,int_stack+546335,int_stack+544760,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+546335,int_stack+300396,int_stack+545435,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+563044,int_stack+546335,int_stack+543760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+534660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+300396,int_stack+275816,int_stack+275606,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+308551,int_stack+276131,int_stack+275816,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3780,int_stack+308551,int_stack+300396,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+303096,int_stack+276572,int_stack+276131,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+303096,int_stack+308551,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+567934,int_stack+566044,int_stack+3780,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+570034,int_stack+567934,int_stack+546335, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+574534,int_stack+570034,int_stack+563044, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+535660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+566044,int_stack+163000,int_stack+162440,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+303096,int_stack+163840,int_stack+163000,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+580534,int_stack+303096,int_stack+566044,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+582214,int_stack+5292,int_stack+163840,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+583978,int_stack+582214,int_stack+303096,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+586498,int_stack+583978,int_stack+580534,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+589298,int_stack+586498,int_stack+567934, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+595598,int_stack+589298,int_stack+570034, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+538660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+582214,int_stack+280886,int_stack+280786,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+306076,int_stack+281036,int_stack+280886,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+582514,int_stack+306076,int_stack+582214,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+309496,int_stack+281246,int_stack+281036,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+583114,int_stack+309496,int_stack+306076,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+584014,int_stack+583114,int_stack+582514,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+583114,int_stack+282416,int_stack+282266,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+585014,int_stack+282641,int_stack+282416,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+585689,int_stack+585014,int_stack+583114,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+586589,int_stack+282956,int_stack+282641,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+587534,int_stack+586589,int_stack+585014,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+588884,int_stack+587534,int_stack+585689,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+590384,int_stack+588884,int_stack+584014, 0.0, zero_stack, 1.0, int_stack+534660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+309496,int_stack+284696,int_stack+284486,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+586589,int_stack+285011,int_stack+284696,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+587534,int_stack+586589,int_stack+309496,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+593384,int_stack+285452,int_stack+285011,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+604598,int_stack+593384,int_stack+586589,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+593384,int_stack+604598,int_stack+587534,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+604598,int_stack+593384,int_stack+588884, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+609098,int_stack+604598,int_stack+590384, 0.0, zero_stack, 1.0, int_stack+535660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+615098,int_stack+202620,int_stack+202060,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+615938,int_stack+203460,int_stack+202620,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+617198,int_stack+615938,int_stack+615098,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+618878,int_stack+6076,int_stack+203460,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+620642,int_stack+618878,int_stack+615938,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+623162,int_stack+620642,int_stack+617198,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+625962,int_stack+623162,int_stack+593384, 0.0, zero_stack, 1.0, int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+632262,int_stack+625962,int_stack+604598, 0.0, zero_stack, 1.0, int_stack+538660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+618878,int_stack+290506,int_stack+290406,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+583564,int_stack+290656,int_stack+290506,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+619178,int_stack+583564,int_stack+618878,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+619778,int_stack+290866,int_stack+290656,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+620408,int_stack+619778,int_stack+583564,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+621308,int_stack+620408,int_stack+619178,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+619778,int_stack+293146,int_stack+292996,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+620228,int_stack+293371,int_stack+293146,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+622308,int_stack+620228,int_stack+619778,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+623208,int_stack+293686,int_stack+293371,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+624153,int_stack+623208,int_stack+620228,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+625503,int_stack+624153,int_stack+622308,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+627003,int_stack+625503,int_stack+621308, 1.0, int_stack+534660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+534660,int_stack+295426,int_stack+295216,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+623208,int_stack+295741,int_stack+295426,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+624153,int_stack+623208,int_stack+534660,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+630003,int_stack+296182,int_stack+295741,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5040,int_stack+630003,int_stack+623208,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+630003,int_stack+5040,int_stack+624153,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+641262,int_stack+630003,int_stack+625503, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+645762,int_stack+641262,int_stack+627003, 1.0, int_stack+535660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+334624,int_stack+269106,int_stack+268546,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5040,int_stack+269946,int_stack+269106,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+535290,int_stack+5040,int_stack+334624,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+651762,int_stack+7448,int_stack+269946,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+653526,int_stack+651762,int_stack+5040,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+656046,int_stack+653526,int_stack+535290,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+658846,int_stack+656046,int_stack+630003, 1.0, int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+665146,int_stack+658846,int_stack+641262, 1.0, int_stack+538660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+318811,int_stack+8332,int_stack+8232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+216372,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+319111,int_stack+8482,int_stack+8332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+216472,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+319561,int_stack+319111,int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+301746,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+320161,int_stack+8692,int_stack+8482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+216622,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+651762,int_stack+320161,int_stack+319111, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+302046,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+652662,int_stack+651762,int_stack+319561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+302496,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+9122,int_stack+8972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+217852,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+318811,int_stack+9347,int_stack+9122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+218002,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+319486,int_stack+318811,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+306526,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653662,int_stack+9662,int_stack+9347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+218227,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+654607,int_stack+653662,int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+306976,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+655957,int_stack+654607,int_stack+319486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+307651,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+657457,int_stack+655957,int_stack+652662,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+318811,int_stack+10292,int_stack+10082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+220072,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+319441,int_stack+10607,int_stack+10292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+220282,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+651762,int_stack+319441,int_stack+318811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+310126,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653022,int_stack+11048,int_stack+10607, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+220597,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660457,int_stack+653022,int_stack+319441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+316606,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+660457,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+317551,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+660457,int_stack+318811,int_stack+655957,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+536970,int_stack+660457,int_stack+657457,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+11916,int_stack+11636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+93776,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+651762,int_stack+12336,int_stack+11916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+94336,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653022,int_stack+651762,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+320911,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654702,int_stack+12924,int_stack+12336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+95176,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656466,int_stack+654702,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+336304,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6300,int_stack+656466,int_stack+653022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+337564,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+651762,int_stack+6300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+674146,int_stack+651762,int_stack+660457,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+13808,int_stack+13708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216372, 1.0, int_stack+225252,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+13958,int_stack+13808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216472, 1.0, int_stack+225352,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301746, 1.0, int_stack+339244,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+14168,int_stack+13958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216622, 1.0, int_stack+225502,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302046, 1.0, int_stack+339544,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302496, 1.0, int_stack+339994,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+14598,int_stack+14448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217852, 1.0, int_stack+226732,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+14823,int_stack+14598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218002, 1.0, int_stack+226882,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306526, 1.0, int_stack+340594,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+15138,int_stack+14823, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218227, 1.0, int_stack+227107,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306976, 1.0, int_stack+341044,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+307651, 1.0, int_stack+343124,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+15768,int_stack+15558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220072, 1.0, int_stack+228952,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+16083,int_stack+15768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220282, 1.0, int_stack+229162,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+310126, 1.0, int_stack+344024,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+16524,int_stack+16083, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220597, 1.0, int_stack+229477,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+316606, 1.0, int_stack+344654,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+317551, 1.0, int_stack+350819,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+6300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+17392,int_stack+17112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93776, 1.0, int_stack+100540,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+17812,int_stack+17392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94336, 1.0, int_stack+101100,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+320911, 1.0, int_stack+361628,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+18400,int_stack+17812, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95176, 1.0, int_stack+101940,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+336304, 1.0, int_stack+362468,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+337564, 1.0, int_stack+377441,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+12300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+683146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+19284,int_stack+19184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+225252, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+19434,int_stack+19284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+225352, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+339244, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+19644,int_stack+19434, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+225502, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+339544, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+339994, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+20074,int_stack+19924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+226732, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+20299,int_stack+20074, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+226882, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+340594, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+20614,int_stack+20299, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+227107, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+341044, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+343124, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+21244,int_stack+21034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+228952, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+21559,int_stack+21244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+229162, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+344024, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+22000,int_stack+21559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+229477, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+344654, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+350819, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+12300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+22868,int_stack+22588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+100540, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+23288,int_stack+22868, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+101100, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+361628, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+23876,int_stack+23288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+101940, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+362468, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+377441, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+18300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+692146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+24760,int_stack+24660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216372, 0.0, zero_stack, 1.0, int_stack+234132,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+24910,int_stack+24760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216472, 0.0, zero_stack, 1.0, int_stack+234232,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301746, 0.0, zero_stack, 1.0, int_stack+379121,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+25120,int_stack+24910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216622, 0.0, zero_stack, 1.0, int_stack+234382,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302046, 0.0, zero_stack, 1.0, int_stack+379421,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302496, 0.0, zero_stack, 1.0, int_stack+379871,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+25550,int_stack+25400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217852, 0.0, zero_stack, 1.0, int_stack+235612,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+25775,int_stack+25550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218002, 0.0, zero_stack, 1.0, int_stack+235762,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306526, 0.0, zero_stack, 1.0, int_stack+380471,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+26090,int_stack+25775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218227, 0.0, zero_stack, 1.0, int_stack+235987,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306976, 0.0, zero_stack, 1.0, int_stack+380921,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+307651, 0.0, zero_stack, 1.0, int_stack+383001,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+26720,int_stack+26510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220072, 0.0, zero_stack, 1.0, int_stack+237832,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+27035,int_stack+26720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220282, 0.0, zero_stack, 1.0, int_stack+238042,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+310126, 0.0, zero_stack, 1.0, int_stack+383901,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+27476,int_stack+27035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220597, 0.0, zero_stack, 1.0, int_stack+238357,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+316606, 0.0, zero_stack, 1.0, int_stack+384531,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+317551, 0.0, zero_stack, 1.0, int_stack+390696,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+18300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+28344,int_stack+28064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93776, 0.0, zero_stack, 1.0, int_stack+107304,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+28764,int_stack+28344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94336, 0.0, zero_stack, 1.0, int_stack+107864,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+320911, 0.0, zero_stack, 1.0, int_stack+401505,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+29352,int_stack+28764, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95176, 0.0, zero_stack, 1.0, int_stack+108704,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+336304, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+24300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+337564, 0.0, zero_stack, 1.0, int_stack+415995,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+24300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+701146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+30236,int_stack+30136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225252, 1.0, int_stack+234132, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+30386,int_stack+30236, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225352, 1.0, int_stack+234232, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+339244, 1.0, int_stack+379121, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+30596,int_stack+30386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225502, 1.0, int_stack+234382, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+339544, 1.0, int_stack+379421, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+339994, 1.0, int_stack+379871, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+31026,int_stack+30876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+226732, 1.0, int_stack+235612, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+31251,int_stack+31026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+226882, 1.0, int_stack+235762, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+340594, 1.0, int_stack+380471, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+31566,int_stack+31251, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+227107, 1.0, int_stack+235987, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+341044, 1.0, int_stack+380921, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+343124, 1.0, int_stack+383001, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+32196,int_stack+31986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+228952, 1.0, int_stack+237832, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+32511,int_stack+32196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+229162, 1.0, int_stack+238042, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+344024, 1.0, int_stack+383901, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+32952,int_stack+32511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+229477, 1.0, int_stack+238357, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+344654, 1.0, int_stack+384531, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+350819, 1.0, int_stack+390696, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+24300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+33820,int_stack+33540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100540, 1.0, int_stack+107304, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+34240,int_stack+33820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101100, 1.0, int_stack+107864, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+361628, 1.0, int_stack+401505, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+34828,int_stack+34240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101940, 1.0, int_stack+108704, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+362468, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+30300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+377441, 1.0, int_stack+415995, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+30300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+710146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+35712,int_stack+35612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+234132, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+35862,int_stack+35712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+234232, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+379121, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+36072,int_stack+35862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+234382, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+379421, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+379871, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+36502,int_stack+36352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+235612, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+36727,int_stack+36502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+235762, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+380471, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+37042,int_stack+36727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+235987, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+380921, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+383001, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+37672,int_stack+37462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+237832, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+37987,int_stack+37672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+238042, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+383901, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+38428,int_stack+37987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+238357, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+384531, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+390696, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+30300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+39296,int_stack+39016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+107304, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+39716,int_stack+39296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+107864, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+401505, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+40304,int_stack+39716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+108704, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+36300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+415995, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+36300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+719146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+41188,int_stack+41088, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216372, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243012,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+41338,int_stack+41188, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216472, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243112,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+301746, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+417675,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+41548,int_stack+41338, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216622, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+243262,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302046, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+417975,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+302496, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+418425,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+41978,int_stack+41828, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+217852, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+244492,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+42203,int_stack+41978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218002, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+244642,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306526, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+419025,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+42518,int_stack+42203, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218227, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+244867,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306976, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+419475,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+307651, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+421555,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+43148,int_stack+42938, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220072, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+246712,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+43463,int_stack+43148, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220282, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+246922,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+310126, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+422455,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+43904,int_stack+43463, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220597, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+247237,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+316606, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+423085,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+317551, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+429250,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+36300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+44772,int_stack+44492, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93776, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114068,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+45192,int_stack+44772, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94336, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114628,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+320911, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+440059,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+45780,int_stack+45192, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95176, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115468,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+336304, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+440899,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+42300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+337564, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+455872,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+42300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+728146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+46664,int_stack+46564, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225252, 0.0, zero_stack, 1.0, int_stack+243012, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+46814,int_stack+46664, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225352, 0.0, zero_stack, 1.0, int_stack+243112, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+339244, 0.0, zero_stack, 1.0, int_stack+417675, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+47024,int_stack+46814, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225502, 0.0, zero_stack, 1.0, int_stack+243262, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+339544, 0.0, zero_stack, 1.0, int_stack+417975, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+339994, 0.0, zero_stack, 1.0, int_stack+418425, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+47454,int_stack+47304, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+226732, 0.0, zero_stack, 1.0, int_stack+244492, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+47679,int_stack+47454, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+226882, 0.0, zero_stack, 1.0, int_stack+244642, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+340594, 0.0, zero_stack, 1.0, int_stack+419025, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+47994,int_stack+47679, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+227107, 0.0, zero_stack, 1.0, int_stack+244867, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+341044, 0.0, zero_stack, 1.0, int_stack+419475, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+343124, 0.0, zero_stack, 1.0, int_stack+421555, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+48624,int_stack+48414, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+228952, 0.0, zero_stack, 1.0, int_stack+246712, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+48939,int_stack+48624, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+229162, 0.0, zero_stack, 1.0, int_stack+246922, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+344024, 0.0, zero_stack, 1.0, int_stack+422455, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+49380,int_stack+48939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+229477, 0.0, zero_stack, 1.0, int_stack+247237, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+344654, 0.0, zero_stack, 1.0, int_stack+423085, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+350819, 0.0, zero_stack, 1.0, int_stack+429250, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+42300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+50248,int_stack+49968, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100540, 0.0, zero_stack, 1.0, int_stack+114068, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+50668,int_stack+50248, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101100, 0.0, zero_stack, 1.0, int_stack+114628, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+361628, 0.0, zero_stack, 1.0, int_stack+440059, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+51256,int_stack+50668, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101940, 0.0, zero_stack, 1.0, int_stack+115468, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+362468, 0.0, zero_stack, 1.0, int_stack+440899, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+48300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+377441, 0.0, zero_stack, 1.0, int_stack+455872, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+48300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+737146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+52140,int_stack+52040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234132, 1.0, int_stack+243012, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+52290,int_stack+52140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234232, 1.0, int_stack+243112, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+379121, 1.0, int_stack+417675, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+52500,int_stack+52290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+234382, 1.0, int_stack+243262, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+379421, 1.0, int_stack+417975, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+379871, 1.0, int_stack+418425, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+52930,int_stack+52780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+235612, 1.0, int_stack+244492, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+53155,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+235762, 1.0, int_stack+244642, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+380471, 1.0, int_stack+419025, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+53470,int_stack+53155, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+235987, 1.0, int_stack+244867, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+380921, 1.0, int_stack+419475, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+383001, 1.0, int_stack+421555, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+54100,int_stack+53890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+237832, 1.0, int_stack+246712, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+54415,int_stack+54100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+238042, 1.0, int_stack+246922, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+383901, 1.0, int_stack+422455, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+54856,int_stack+54415, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+238357, 1.0, int_stack+247237, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+384531, 1.0, int_stack+423085, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+390696, 1.0, int_stack+429250, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+48300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+55724,int_stack+55444, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107304, 1.0, int_stack+114068, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+56144,int_stack+55724, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107864, 1.0, int_stack+114628, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+401505, 1.0, int_stack+440059, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+56732,int_stack+56144, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108704, 1.0, int_stack+115468, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+440899, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+54300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+415995, 1.0, int_stack+455872, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+54300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+746146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+57616,int_stack+57516, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+243012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+57766,int_stack+57616, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+243112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+417675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+57976,int_stack+57766, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+243262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+417975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+418425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+58406,int_stack+58256, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+244492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+58631,int_stack+58406, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+244642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+419025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+58946,int_stack+58631, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+244867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+419475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+421555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+59576,int_stack+59366, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+246712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+59891,int_stack+59576, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+246922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+422455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+60332,int_stack+59891, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+247237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+423085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+429250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+54300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+61200,int_stack+60920, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+114068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+61620,int_stack+61200, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+114628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+440059, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+62208,int_stack+61620, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+115468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+440899, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+755146,int_stack+660966,int_stack+657522, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+455872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+755146,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+755146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+63092,int_stack+62992, 0.0, zero_stack, 1.0, int_stack+216372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+251892,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+63242,int_stack+63092, 0.0, zero_stack, 1.0, int_stack+216472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+251992,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+301746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+457552,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+63452,int_stack+63242, 0.0, zero_stack, 1.0, int_stack+216622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+252142,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 1.0, int_stack+302046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+457852,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 1.0, int_stack+302496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+458302,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+63882,int_stack+63732, 0.0, zero_stack, 1.0, int_stack+217852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+253372,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+64107,int_stack+63882, 0.0, zero_stack, 1.0, int_stack+218002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+253522,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+306526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+458902,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+64422,int_stack+64107, 0.0, zero_stack, 1.0, int_stack+218227, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+253747,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 1.0, int_stack+306976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+459352,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 1.0, int_stack+307651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+461432,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+65052,int_stack+64842, 0.0, zero_stack, 1.0, int_stack+220072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+255592,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+65367,int_stack+65052, 0.0, zero_stack, 1.0, int_stack+220282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+255802,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+310126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+462332,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+65808,int_stack+65367, 0.0, zero_stack, 1.0, int_stack+220597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+256117,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 1.0, int_stack+316606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+462962,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 1.0, int_stack+317551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+469127,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+60300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+66676,int_stack+66396, 0.0, zero_stack, 1.0, int_stack+93776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120832,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+67096,int_stack+66676, 0.0, zero_stack, 1.0, int_stack+94336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121392,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 1.0, int_stack+320911, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+479936,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+67684,int_stack+67096, 0.0, zero_stack, 1.0, int_stack+95176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122232,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 1.0, int_stack+336304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+764146,int_stack+660966,int_stack+657522, 0.0, zero_stack, 1.0, int_stack+337564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+494426,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+764146,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+764146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+68568,int_stack+68468, 0.0, zero_stack, 1.0, int_stack+225252, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+251892, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+68718,int_stack+68568, 0.0, zero_stack, 1.0, int_stack+225352, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+251992, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+339244, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+457552, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+68928,int_stack+68718, 0.0, zero_stack, 1.0, int_stack+225502, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+252142, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 1.0, int_stack+339544, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+457852, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 1.0, int_stack+339994, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+458302, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+69358,int_stack+69208, 0.0, zero_stack, 1.0, int_stack+226732, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+253372, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+69583,int_stack+69358, 0.0, zero_stack, 1.0, int_stack+226882, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+253522, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+340594, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+458902, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+69898,int_stack+69583, 0.0, zero_stack, 1.0, int_stack+227107, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+253747, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 1.0, int_stack+341044, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+459352, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 1.0, int_stack+343124, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+461432, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+70528,int_stack+70318, 0.0, zero_stack, 1.0, int_stack+228952, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+255592, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+70843,int_stack+70528, 0.0, zero_stack, 1.0, int_stack+229162, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+255802, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+344024, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+462332, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+71284,int_stack+70843, 0.0, zero_stack, 1.0, int_stack+229477, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+256117, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 1.0, int_stack+344654, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+462962, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 1.0, int_stack+350819, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+469127, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+773146,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+72152,int_stack+71872, 0.0, zero_stack, 1.0, int_stack+100540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120832, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+72572,int_stack+72152, 0.0, zero_stack, 1.0, int_stack+101100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121392, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 1.0, int_stack+361628, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+479936, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+73160,int_stack+72572, 0.0, zero_stack, 1.0, int_stack+101940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122232, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 1.0, int_stack+362468, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+66300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 1.0, int_stack+377441, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+494426, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+66300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+779146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+74044,int_stack+73944, 0.0, zero_stack, 1.0, int_stack+234132, 0.0, zero_stack, 1.0, int_stack+251892, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+74194,int_stack+74044, 0.0, zero_stack, 1.0, int_stack+234232, 0.0, zero_stack, 1.0, int_stack+251992, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+379121, 0.0, zero_stack, 1.0, int_stack+457552, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+74404,int_stack+74194, 0.0, zero_stack, 1.0, int_stack+234382, 0.0, zero_stack, 1.0, int_stack+252142, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 1.0, int_stack+379421, 0.0, zero_stack, 1.0, int_stack+457852, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 1.0, int_stack+379871, 0.0, zero_stack, 1.0, int_stack+458302, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+74834,int_stack+74684, 0.0, zero_stack, 1.0, int_stack+235612, 0.0, zero_stack, 1.0, int_stack+253372, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+75059,int_stack+74834, 0.0, zero_stack, 1.0, int_stack+235762, 0.0, zero_stack, 1.0, int_stack+253522, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+380471, 0.0, zero_stack, 1.0, int_stack+458902, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+75374,int_stack+75059, 0.0, zero_stack, 1.0, int_stack+235987, 0.0, zero_stack, 1.0, int_stack+253747, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 1.0, int_stack+380921, 0.0, zero_stack, 1.0, int_stack+459352, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 1.0, int_stack+383001, 0.0, zero_stack, 1.0, int_stack+461432, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+76004,int_stack+75794, 0.0, zero_stack, 1.0, int_stack+237832, 0.0, zero_stack, 1.0, int_stack+255592, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+76319,int_stack+76004, 0.0, zero_stack, 1.0, int_stack+238042, 0.0, zero_stack, 1.0, int_stack+255802, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+383901, 0.0, zero_stack, 1.0, int_stack+462332, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+76760,int_stack+76319, 0.0, zero_stack, 1.0, int_stack+238357, 0.0, zero_stack, 1.0, int_stack+256117, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 1.0, int_stack+384531, 0.0, zero_stack, 1.0, int_stack+462962, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 1.0, int_stack+390696, 0.0, zero_stack, 1.0, int_stack+469127, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+66300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+77628,int_stack+77348, 0.0, zero_stack, 1.0, int_stack+107304, 0.0, zero_stack, 1.0, int_stack+120832, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+78048,int_stack+77628, 0.0, zero_stack, 1.0, int_stack+107864, 0.0, zero_stack, 1.0, int_stack+121392, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 1.0, int_stack+401505, 0.0, zero_stack, 1.0, int_stack+479936, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+78636,int_stack+78048, 0.0, zero_stack, 1.0, int_stack+108704, 0.0, zero_stack, 1.0, int_stack+122232, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+72300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 1.0, int_stack+415995, 0.0, zero_stack, 1.0, int_stack+494426, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+72300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+788146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+79520,int_stack+79420, 0.0, zero_stack, 1.0, int_stack+243012, 1.0, int_stack+251892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+79670,int_stack+79520, 0.0, zero_stack, 1.0, int_stack+243112, 1.0, int_stack+251992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+417675, 1.0, int_stack+457552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+79880,int_stack+79670, 0.0, zero_stack, 1.0, int_stack+243262, 1.0, int_stack+252142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 1.0, int_stack+417975, 1.0, int_stack+457852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 1.0, int_stack+418425, 1.0, int_stack+458302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+80310,int_stack+80160, 0.0, zero_stack, 1.0, int_stack+244492, 1.0, int_stack+253372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+80535,int_stack+80310, 0.0, zero_stack, 1.0, int_stack+244642, 1.0, int_stack+253522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+419025, 1.0, int_stack+458902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+80850,int_stack+80535, 0.0, zero_stack, 1.0, int_stack+244867, 1.0, int_stack+253747, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 1.0, int_stack+419475, 1.0, int_stack+459352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 1.0, int_stack+421555, 1.0, int_stack+461432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+81480,int_stack+81270, 0.0, zero_stack, 1.0, int_stack+246712, 1.0, int_stack+255592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+81795,int_stack+81480, 0.0, zero_stack, 1.0, int_stack+246922, 1.0, int_stack+255802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 1.0, int_stack+422455, 1.0, int_stack+462332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+82236,int_stack+81795, 0.0, zero_stack, 1.0, int_stack+247237, 1.0, int_stack+256117, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 1.0, int_stack+423085, 1.0, int_stack+462962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 1.0, int_stack+429250, 1.0, int_stack+469127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+72300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+83104,int_stack+82824, 0.0, zero_stack, 1.0, int_stack+114068, 1.0, int_stack+120832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+83524,int_stack+83104, 0.0, zero_stack, 1.0, int_stack+114628, 1.0, int_stack+121392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 1.0, int_stack+440059, 1.0, int_stack+479936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+84112,int_stack+83524, 0.0, zero_stack, 1.0, int_stack+115468, 1.0, int_stack+122232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 1.0, int_stack+440899, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+78300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 1.0, int_stack+455872, 1.0, int_stack+494426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+78300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+797146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+84996,int_stack+84896, 0.0, zero_stack, 2.0, int_stack+251892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+85146,int_stack+84996, 0.0, zero_stack, 2.0, int_stack+251992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 0.0, zero_stack, 2.0, int_stack+457552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+85356,int_stack+85146, 0.0, zero_stack, 2.0, int_stack+252142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 0.0, zero_stack, 2.0, int_stack+457852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 0.0, zero_stack, 2.0, int_stack+458302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+85786,int_stack+85636, 0.0, zero_stack, 2.0, int_stack+253372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+86011,int_stack+85786, 0.0, zero_stack, 2.0, int_stack+253522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 0.0, zero_stack, 2.0, int_stack+458902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+86326,int_stack+86011, 0.0, zero_stack, 2.0, int_stack+253747, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 0.0, zero_stack, 2.0, int_stack+459352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 0.0, zero_stack, 2.0, int_stack+461432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+86956,int_stack+86746, 0.0, zero_stack, 2.0, int_stack+255592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+87271,int_stack+86956, 0.0, zero_stack, 2.0, int_stack+255802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+653337,int_stack+652392,int_stack+651762, 0.0, zero_stack, 2.0, int_stack+462332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+654597,int_stack+87712,int_stack+87271, 0.0, zero_stack, 2.0, int_stack+256117, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+655920,int_stack+654597,int_stack+652392, 0.0, zero_stack, 2.0, int_stack+462962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+655920,int_stack+653337, 0.0, zero_stack, 2.0, int_stack+469127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+78300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+88580,int_stack+88300, 0.0, zero_stack, 2.0, int_stack+120832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+89000,int_stack+88580, 0.0, zero_stack, 2.0, int_stack+121392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 0.0, zero_stack, 2.0, int_stack+479936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+89588,int_stack+89000, 0.0, zero_stack, 2.0, int_stack+122232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 0.0, zero_stack, 2.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+84300,int_stack+660966,int_stack+657522, 0.0, zero_stack, 2.0, int_stack+494426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+84300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+806146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+90472,int_stack+90372, 1.0, int_stack+216372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260772,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+90622,int_stack+90472, 1.0, int_stack+216472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260872,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 1.0, int_stack+301746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496106,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+90832,int_stack+90622, 1.0, int_stack+216622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+261022,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 1.0, int_stack+302046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496406,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 1.0, int_stack+302496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496856,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+91262,int_stack+91112, 1.0, int_stack+217852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262252,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+91487,int_stack+91262, 1.0, int_stack+218002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262402,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 1.0, int_stack+306526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+497456,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+91802,int_stack+91487, 1.0, int_stack+218227, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262627,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 1.0, int_stack+306976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+497906,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 1.0, int_stack+307651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+499986,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+306526,int_stack+92432,int_stack+92222, 1.0, int_stack+220072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264682,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+307156,int_stack+92747,int_stack+92432, 1.0, int_stack+220282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265207,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+651762,int_stack+307156,int_stack+306526, 1.0, int_stack+310126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500886,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653022,int_stack+93188,int_stack+92747, 1.0, int_stack+220597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265963,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+654345,int_stack+653022,int_stack+307156, 1.0, int_stack+316606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+501516,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+318811,int_stack+654345,int_stack+651762, 1.0, int_stack+317551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+507681,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+318811,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+84300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+94756,int_stack+94056, 1.0, int_stack+93776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127876,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+95764,int_stack+94756, 1.0, int_stack+94336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128856,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+335464, 1.0, int_stack+320911, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+518490,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+659202,int_stack+96352,int_stack+95764, 1.0, int_stack+95176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129696,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+660966,int_stack+659202,int_stack+656262, 1.0, int_stack+336304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2520,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+90300,int_stack+660966,int_stack+657522, 1.0, int_stack+337564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+532980,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+90300,int_stack+318811,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+815146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+97236,int_stack+97136, 1.0, int_stack+225252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260772, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+97386,int_stack+97236, 1.0, int_stack+225352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260872, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 1.0, int_stack+339244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496106, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+97596,int_stack+97386, 1.0, int_stack+225502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+261022, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 1.0, int_stack+339544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496406, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 1.0, int_stack+339994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496856, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+98026,int_stack+97876, 1.0, int_stack+226732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262252, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+98251,int_stack+98026, 1.0, int_stack+226882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262402, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 1.0, int_stack+340594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+497456, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+98566,int_stack+98251, 1.0, int_stack+227107, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262627, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 1.0, int_stack+341044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+497906, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 1.0, int_stack+343124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+499986, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+343124,int_stack+99196,int_stack+98986, 1.0, int_stack+228952, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264682, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+651762,int_stack+99511,int_stack+99196, 1.0, int_stack+229162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265207, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652707,int_stack+651762,int_stack+343124, 1.0, int_stack+344024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500886, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+343124,int_stack+99952,int_stack+99511, 1.0, int_stack+229477, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265963, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653967,int_stack+343124,int_stack+651762, 1.0, int_stack+344654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+501516, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+343124,int_stack+653967,int_stack+652707, 1.0, int_stack+350819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+507681, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+343124,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+90300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+350819,int_stack+101520,int_stack+100820, 1.0, int_stack+100540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127876, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+102528,int_stack+101520, 1.0, int_stack+101100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128856, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+350819, 1.0, int_stack+361628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+518490, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+350819,int_stack+103116,int_stack+102528, 1.0, int_stack+101940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129696, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+659202,int_stack+350819,int_stack+656262, 1.0, int_stack+362468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2520, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+361628,int_stack+659202,int_stack+657522, 1.0, int_stack+377441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+532980, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+656262,int_stack+361628,int_stack+343124,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+824146,int_stack+656262,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+104000,int_stack+103900, 1.0, int_stack+234132, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260772, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+104150,int_stack+104000, 1.0, int_stack+234232, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+260872, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 1.0, int_stack+379121, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496106, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+104360,int_stack+104150, 1.0, int_stack+234382, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+261022, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 1.0, int_stack+379421, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496406, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 1.0, int_stack+379871, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+496856, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+104790,int_stack+104640, 1.0, int_stack+235612, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262252, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+105015,int_stack+104790, 1.0, int_stack+235762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262402, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 1.0, int_stack+380471, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+497456, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+105330,int_stack+105015, 1.0, int_stack+235987, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262627, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 1.0, int_stack+380921, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+497906, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 1.0, int_stack+383001, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+499986, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+383001,int_stack+105960,int_stack+105750, 1.0, int_stack+237832, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+264682, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+651762,int_stack+106275,int_stack+105960, 1.0, int_stack+238042, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265207, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652707,int_stack+651762,int_stack+383001, 1.0, int_stack+383901, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500886, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+383001,int_stack+106716,int_stack+106275, 1.0, int_stack+238357, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+265963, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653967,int_stack+383001,int_stack+651762, 1.0, int_stack+384531, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+501516, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+383001,int_stack+653967,int_stack+652707, 1.0, int_stack+390696, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+507681, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+651762,int_stack+383001,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+96300,int_stack+651762,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+390696,int_stack+108284,int_stack+107584, 1.0, int_stack+107304, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127876, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+656262,int_stack+109292,int_stack+108284, 1.0, int_stack+107864, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128856, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+657522,int_stack+656262,int_stack+390696, 1.0, int_stack+401505, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+518490, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+401505,int_stack+109880,int_stack+109292, 1.0, int_stack+108704, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129696, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+659202,int_stack+401505,int_stack+656262, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+661722,int_stack+659202,int_stack+657522, 1.0, int_stack+415995, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+532980, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+102300,int_stack+661722,int_stack+383001,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+833146,int_stack+102300,int_stack+651762,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+110764,int_stack+110664, 1.0, int_stack+243012, 0.0, zero_stack, 1.0, int_stack+260772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652062,int_stack+110914,int_stack+110764, 1.0, int_stack+243112, 0.0, zero_stack, 1.0, int_stack+260872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652512,int_stack+652062,int_stack+651762, 1.0, int_stack+417675, 0.0, zero_stack, 1.0, int_stack+496106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+653112,int_stack+111124,int_stack+110914, 1.0, int_stack+243262, 0.0, zero_stack, 1.0, int_stack+261022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653742,int_stack+653112,int_stack+652062, 1.0, int_stack+417975, 0.0, zero_stack, 1.0, int_stack+496406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654642,int_stack+653742,int_stack+652512, 1.0, int_stack+418425, 0.0, zero_stack, 1.0, int_stack+496856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+111554,int_stack+111404, 1.0, int_stack+244492, 0.0, zero_stack, 1.0, int_stack+262252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652212,int_stack+111779,int_stack+111554, 1.0, int_stack+244642, 0.0, zero_stack, 1.0, int_stack+262402, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+652887,int_stack+652212,int_stack+651762, 1.0, int_stack+419025, 0.0, zero_stack, 1.0, int_stack+497456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+655642,int_stack+112094,int_stack+111779, 1.0, int_stack+244867, 0.0, zero_stack, 1.0, int_stack+262627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+656587,int_stack+655642,int_stack+652212, 1.0, int_stack+419475, 0.0, zero_stack, 1.0, int_stack+497906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+657937,int_stack+656587,int_stack+652887, 1.0, int_stack+421555, 0.0, zero_stack, 1.0, int_stack+499986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+659437,int_stack+657937,int_stack+654642,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+421555,int_stack+112724,int_stack+112514, 1.0, int_stack+246712, 0.0, zero_stack, 1.0, int_stack+264682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+651762,int_stack+113039,int_stack+112724, 1.0, int_stack+246922, 0.0, zero_stack, 1.0, int_stack+265207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+651762,int_stack+421555, 1.0, int_stack+422455, 0.0, zero_stack, 1.0, int_stack+500886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+421555,int_stack+113480,int_stack+113039, 1.0, int_stack+247237, 0.0, zero_stack, 1.0, int_stack+265963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+421555,int_stack+651762, 1.0, int_stack+423085, 0.0, zero_stack, 1.0, int_stack+501516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+651762,int_stack+401505,int_stack+0, 1.0, int_stack+429250, 0.0, zero_stack, 1.0, int_stack+507681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+651762,int_stack+657937,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+106800,int_stack+102300,int_stack+659437,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+429250,int_stack+115048,int_stack+114348, 1.0, int_stack+114068, 0.0, zero_stack, 1.0, int_stack+127876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+116056,int_stack+115048, 1.0, int_stack+114628, 0.0, zero_stack, 1.0, int_stack+128856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+401505,int_stack+0,int_stack+429250, 1.0, int_stack+440059, 0.0, zero_stack, 1.0, int_stack+518490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+429250,int_stack+116644,int_stack+116056, 1.0, int_stack+115468, 0.0, zero_stack, 1.0, int_stack+129696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653862,int_stack+429250,int_stack+0, 1.0, int_stack+440899, 0.0, zero_stack, 1.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+440059,int_stack+653862,int_stack+401505, 1.0, int_stack+455872, 0.0, zero_stack, 1.0, int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+653862,int_stack+440059,int_stack+651762,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+842146,int_stack+653862,int_stack+102300,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+117528,int_stack+117428, 1.0, int_stack+251892, 1.0, int_stack+260772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102600,int_stack+117678,int_stack+117528, 1.0, int_stack+251992, 1.0, int_stack+260872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103050,int_stack+102600,int_stack+102300, 1.0, int_stack+457552, 1.0, int_stack+496106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103650,int_stack+117888,int_stack+117678, 1.0, int_stack+252142, 1.0, int_stack+261022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+104280,int_stack+103650,int_stack+102600, 1.0, int_stack+457852, 1.0, int_stack+496406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+105180,int_stack+104280,int_stack+103050, 1.0, int_stack+458302, 1.0, int_stack+496856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+118318,int_stack+118168, 1.0, int_stack+253372, 1.0, int_stack+262252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+118543,int_stack+118318, 1.0, int_stack+253522, 1.0, int_stack+262402, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 1.0, int_stack+458902, 1.0, int_stack+497456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+651762,int_stack+118858,int_stack+118543, 1.0, int_stack+253747, 1.0, int_stack+262627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+652707,int_stack+651762,int_stack+102750, 1.0, int_stack+459352, 1.0, int_stack+497906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654057,int_stack+652707,int_stack+103425, 1.0, int_stack+461432, 1.0, int_stack+499986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+655557,int_stack+654057,int_stack+105180,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+461432,int_stack+119488,int_stack+119278, 1.0, int_stack+255592, 1.0, int_stack+264682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+651762,int_stack+119803,int_stack+119488, 1.0, int_stack+255802, 1.0, int_stack+265207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+651762,int_stack+461432, 1.0, int_stack+462332, 1.0, int_stack+500886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+461432,int_stack+120244,int_stack+119803, 1.0, int_stack+256117, 1.0, int_stack+265963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+461432,int_stack+651762, 1.0, int_stack+462962, 1.0, int_stack+501516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+651762,int_stack+401505,int_stack+0, 1.0, int_stack+469127, 1.0, int_stack+507681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+651762,int_stack+654057,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+658557,int_stack+102300,int_stack+655557,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+469127,int_stack+121812,int_stack+121112, 1.0, int_stack+120832, 1.0, int_stack+127876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+122820,int_stack+121812, 1.0, int_stack+121392, 1.0, int_stack+128856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+401505,int_stack+0,int_stack+469127, 1.0, int_stack+479936, 1.0, int_stack+518490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+479936,int_stack+123408,int_stack+122820, 1.0, int_stack+122232, 1.0, int_stack+129696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+653862,int_stack+479936,int_stack+0, 1.0, int_stack+1260, 1.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+461432,int_stack+653862,int_stack+401505, 1.0, int_stack+494426, 1.0, int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+112800,int_stack+461432,int_stack+651762,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+851146,int_stack+112800,int_stack+102300,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+124292,int_stack+124192, 2.0, int_stack+260772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102600,int_stack+124442,int_stack+124292, 2.0, int_stack+260872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103050,int_stack+102600,int_stack+102300, 2.0, int_stack+496106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103650,int_stack+124652,int_stack+124442, 2.0, int_stack+261022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+104280,int_stack+103650,int_stack+102600, 2.0, int_stack+496406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+105180,int_stack+104280,int_stack+103050, 2.0, int_stack+496856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+125082,int_stack+124932, 2.0, int_stack+262252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+125307,int_stack+125082, 2.0, int_stack+262402, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 2.0, int_stack+497456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+112800,int_stack+125622,int_stack+125307, 2.0, int_stack+262627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+113745,int_stack+112800,int_stack+102750, 2.0, int_stack+497906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+115095,int_stack+113745,int_stack+103425, 2.0, int_stack+499986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+116595,int_stack+115095,int_stack+105180,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+499986,int_stack+126252,int_stack+126042, 2.0, int_stack+264682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+112800,int_stack+126567,int_stack+126252, 2.0, int_stack+265207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+113745,int_stack+112800,int_stack+499986, 2.0, int_stack+500886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+499986,int_stack+127008,int_stack+126567, 2.0, int_stack+265963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+499986,int_stack+112800, 2.0, int_stack+501516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+499986,int_stack+401505,int_stack+113745, 2.0, int_stack+507681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+499986,int_stack+115095,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+119595,int_stack+102300,int_stack+116595,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+507681,int_stack+129276,int_stack+128156, 2.0, int_stack+127876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+401505,int_stack+130284,int_stack+129276, 2.0, int_stack+128856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+112800,int_stack+401505,int_stack+507681, 2.0, int_stack+518490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+518490,int_stack+130872,int_stack+130284, 2.0, int_stack+129696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+518490,int_stack+401505, 2.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+114480,int_stack+0,int_stack+112800, 2.0, int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+651762,int_stack+114480,int_stack+499986,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+860146,int_stack+651762,int_stack+102300,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+131756,int_stack+131656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+271906,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102600,int_stack+131906,int_stack+131756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272006,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103050,int_stack+102600,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103650,int_stack+132116,int_stack+131906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272156,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+104280,int_stack+103650,int_stack+102600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315646,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+105180,int_stack+104280,int_stack+103050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+543160,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+132546,int_stack+132396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273386,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+132771,int_stack+132546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273536,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+651762,int_stack+133086,int_stack+132771, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273761,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+652707,int_stack+651762,int_stack+102750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+544760,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+654057,int_stack+652707,int_stack+103425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+545435,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+655557,int_stack+654057,int_stack+105180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+304626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+651762,int_stack+133716,int_stack+133506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275606,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+652392,int_stack+134031,int_stack+133716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275816,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102300,int_stack+652392,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103560,int_stack+134472,int_stack+134031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276131,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+103560,int_stack+652392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+103560,int_stack+401505,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3780,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+112800,int_stack+103560,int_stack+654057, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+310846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+125595,int_stack+112800,int_stack+655557, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+312346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+135340,int_stack+135060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162440,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+401505,int_stack+135760,int_stack+135340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163000,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+401505,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+566044,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+651762,int_stack+136348,int_stack+135760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163840,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+651762,int_stack+401505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+651762,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+580534,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+335464,int_stack+651762,int_stack+103560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+322024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+869146,int_stack+335464,int_stack+112800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+112800,int_stack+137232,int_stack+137132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+271906, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113100,int_stack+137382,int_stack+137232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272006, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+113550,int_stack+113100,int_stack+112800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+114150,int_stack+137592,int_stack+137382, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272156, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+114780,int_stack+114150,int_stack+113100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315646, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+115680,int_stack+114780,int_stack+113550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+543160, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+112800,int_stack+138022,int_stack+137872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273386, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113250,int_stack+138247,int_stack+138022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273536, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+113925,int_stack+113250,int_stack+112800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+116680,int_stack+138562,int_stack+138247, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273761, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+117625,int_stack+116680,int_stack+113250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+544760, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+335464,int_stack+117625,int_stack+113925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+545435, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+336964,int_stack+335464,int_stack+115680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+342124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+112800,int_stack+139192,int_stack+138982, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275606, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113430,int_stack+139507,int_stack+139192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275816, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+114375,int_stack+113430,int_stack+112800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+115635,int_stack+139948,int_stack+139507, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276131, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+115635,int_stack+113430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+115635,int_stack+401505,int_stack+114375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3780, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+115635,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+346319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+651762,int_stack+102300,int_stack+336964, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+347819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+335464,int_stack+140816,int_stack+140536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162440, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+336304,int_stack+141236,int_stack+140816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163000, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+336304,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+566044, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+337564,int_stack+141824,int_stack+141236, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163840, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+337564,int_stack+336304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+335464,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+580534, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+131595,int_stack+335464,int_stack+115635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+364841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+878146,int_stack+131595,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+366941, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+142708,int_stack+142608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+271906, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102600,int_stack+142858,int_stack+142708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272006, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103050,int_stack+102600,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103650,int_stack+143068,int_stack+142858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272156, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+104280,int_stack+103650,int_stack+102600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315646, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+105180,int_stack+104280,int_stack+103050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+543160, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+143498,int_stack+143348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273386, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+143723,int_stack+143498, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273536, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131595,int_stack+144038,int_stack+143723, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273761, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+132540,int_stack+131595,int_stack+102750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+544760, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+133890,int_stack+132540,int_stack+103425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+545435, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+135390,int_stack+133890,int_stack+105180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+382001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+144668,int_stack+144458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275606, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+132225,int_stack+144983,int_stack+144668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275816, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102300,int_stack+132225,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103560,int_stack+145424,int_stack+144983, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276131, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+103560,int_stack+132225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+103560,int_stack+401505,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3780, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+138390,int_stack+103560,int_stack+133890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+386196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+335464,int_stack+138390,int_stack+135390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+387696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+146292,int_stack+146012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162440, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+401505,int_stack+146712,int_stack+146292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163000, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+401505,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+566044, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131595,int_stack+147300,int_stack+146712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163840, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+131595,int_stack+401505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+131595,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+580534, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+112800,int_stack+131595,int_stack+103560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+403395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+887146,int_stack+112800,int_stack+138390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+405495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+112800,int_stack+148184,int_stack+148084, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+271906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113100,int_stack+148334,int_stack+148184, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+113550,int_stack+113100,int_stack+112800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+114150,int_stack+148544,int_stack+148334, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+272156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+114780,int_stack+114150,int_stack+113100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+115680,int_stack+114780,int_stack+113550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+543160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+112800,int_stack+148974,int_stack+148824, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113250,int_stack+149199,int_stack+148974, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+113925,int_stack+113250,int_stack+112800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+116680,int_stack+149514,int_stack+149199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+273761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+117625,int_stack+116680,int_stack+113250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+544760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+131595,int_stack+117625,int_stack+113925, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+545435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+133095,int_stack+131595,int_stack+115680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+112800,int_stack+150144,int_stack+149934, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113430,int_stack+150459,int_stack+150144, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+275816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+114375,int_stack+113430,int_stack+112800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+115635,int_stack+150900,int_stack+150459, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+276131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+115635,int_stack+113430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+115635,int_stack+401505,int_stack+114375, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+115635,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+424750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+136095,int_stack+102300,int_stack+133095, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+426250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+151768,int_stack+151488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+132435,int_stack+152188,int_stack+151768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+132435,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+133695,int_stack+152776,int_stack+152188, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+133695,int_stack+132435, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+131595,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+580534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+142095,int_stack+131595,int_stack+115635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+443272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+896146,int_stack+142095,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+445372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+153660,int_stack+153560, 0.0, zero_stack, 1.0, int_stack+271906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102600,int_stack+153810,int_stack+153660, 0.0, zero_stack, 1.0, int_stack+272006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103050,int_stack+102600,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103650,int_stack+154020,int_stack+153810, 0.0, zero_stack, 1.0, int_stack+272156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+104280,int_stack+103650,int_stack+102600, 0.0, zero_stack, 1.0, int_stack+315646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+105180,int_stack+104280,int_stack+103050, 0.0, zero_stack, 1.0, int_stack+543160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+154450,int_stack+154300, 0.0, zero_stack, 1.0, int_stack+273386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+154675,int_stack+154450, 0.0, zero_stack, 1.0, int_stack+273536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+142095,int_stack+154990,int_stack+154675, 0.0, zero_stack, 1.0, int_stack+273761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+143040,int_stack+142095,int_stack+102750, 0.0, zero_stack, 1.0, int_stack+544760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+144390,int_stack+143040,int_stack+103425, 0.0, zero_stack, 1.0, int_stack+545435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+145890,int_stack+144390,int_stack+105180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+460432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+142095,int_stack+155620,int_stack+155410, 0.0, zero_stack, 1.0, int_stack+275606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+142725,int_stack+155935,int_stack+155620, 0.0, zero_stack, 1.0, int_stack+275816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102300,int_stack+142725,int_stack+142095, 0.0, zero_stack, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103560,int_stack+156376,int_stack+155935, 0.0, zero_stack, 1.0, int_stack+276131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+103560,int_stack+142725, 0.0, zero_stack, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+103560,int_stack+401505,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+3780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+103560,int_stack+144390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+464627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+148890,int_stack+131595,int_stack+145890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+466127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+157244,int_stack+156964, 0.0, zero_stack, 1.0, int_stack+162440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+401505,int_stack+157664,int_stack+157244, 0.0, zero_stack, 1.0, int_stack+163000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+401505,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+142095,int_stack+158252,int_stack+157664, 0.0, zero_stack, 1.0, int_stack+163840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+142095,int_stack+401505, 0.0, zero_stack, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+142095,int_stack+264472,int_stack+532980, 0.0, zero_stack, 1.0, int_stack+580534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+112800,int_stack+142095,int_stack+103560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+481826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+905146,int_stack+112800,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+483926, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+159136,int_stack+159036, 1.0, int_stack+271906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131895,int_stack+159286,int_stack+159136, 1.0, int_stack+272006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132345,int_stack+131895,int_stack+131595, 1.0, int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132945,int_stack+159496,int_stack+159286, 1.0, int_stack+272156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+133575,int_stack+132945,int_stack+131895, 1.0, int_stack+315646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+315346,int_stack+133575,int_stack+132345, 1.0, int_stack+543160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+316346,int_stack+159926,int_stack+159776, 1.0, int_stack+273386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+316796,int_stack+160151,int_stack+159926, 1.0, int_stack+273536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+317471,int_stack+316796,int_stack+316346, 1.0, int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+318371,int_stack+160466,int_stack+160151, 1.0, int_stack+273761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+319316,int_stack+318371,int_stack+316796, 1.0, int_stack+544760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+131595,int_stack+319316,int_stack+317471, 1.0, int_stack+545435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+133095,int_stack+131595,int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+498986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+315346,int_stack+161096,int_stack+160886, 1.0, int_stack+275606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+315976,int_stack+161411,int_stack+161096, 1.0, int_stack+275816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+316921,int_stack+315976,int_stack+315346, 1.0, int_stack+300396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+300396,int_stack+161852,int_stack+161411, 1.0, int_stack+276131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+401505,int_stack+300396,int_stack+315976, 1.0, int_stack+308551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+300396,int_stack+401505,int_stack+316921, 1.0, int_stack+3780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+300396,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+503181, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+315346,int_stack+102300,int_stack+133095, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+504681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+163420,int_stack+162720, 1.0, int_stack+162440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+132435,int_stack+164428,int_stack+163420, 1.0, int_stack+163000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+132435,int_stack+131595, 1.0, int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+566044,int_stack+165016,int_stack+164428, 1.0, int_stack+163840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+566044,int_stack+132435, 1.0, int_stack+303096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+131595,int_stack+264472,int_stack+532980, 1.0, int_stack+580534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+112800,int_stack+131595,int_stack+300396, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+520380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+154890,int_stack+112800,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+522480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+165900,int_stack+165800,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+305626,int_stack+166050,int_stack+165900,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+102600,int_stack+305626,int_stack+102300,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+103200,int_stack+166260,int_stack+166050,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+103830,int_stack+103200,int_stack+305626,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+104730,int_stack+103830,int_stack+102600,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+166690,int_stack+166540,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+166915,int_stack+166690,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+102975,int_stack+102300,int_stack+305626,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+105730,int_stack+167230,int_stack+166915,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+112800,int_stack+105730,int_stack+102300,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+114150,int_stack+112800,int_stack+102975,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+115650,int_stack+114150,int_stack+104730, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+543760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+112800,int_stack+167860,int_stack+167650,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+118650,int_stack+168175,int_stack+167860,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+102300,int_stack+118650,int_stack+112800,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+112800,int_stack+168616,int_stack+168175,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+112800,int_stack+118650,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+103560,int_stack+566044,int_stack+102300,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+103560,int_stack+114150, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+546335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+142095,int_stack+131595,int_stack+115650, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+563044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+169484,int_stack+169204,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+566044,int_stack+169904,int_stack+169484,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+580534,int_stack+566044,int_stack+102300,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+112800,int_stack+170492,int_stack+169904,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+112800,int_stack+566044,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+112800,int_stack+264472,int_stack+580534,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+163890,int_stack+112800,int_stack+103560, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+567934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+914146,int_stack+163890,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+570034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+171376,int_stack+171276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280786,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+305626,int_stack+171526,int_stack+171376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280886,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+131895,int_stack+305626,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582214,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132495,int_stack+171736,int_stack+171526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+281036,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+133125,int_stack+132495,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306076,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+134025,int_stack+133125,int_stack+131895, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582514,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+172166,int_stack+172016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282266,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+172391,int_stack+172166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282416,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132270,int_stack+131595,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583114,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+135025,int_stack+172706,int_stack+172391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282641,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+163890,int_stack+135025,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585014,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+165240,int_stack+163890,int_stack+132270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585689,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+166740,int_stack+165240,int_stack+134025, 0.0, zero_stack, 1.0, int_stack+304626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+163890,int_stack+173336,int_stack+173126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284486,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+173651,int_stack+173336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284696,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132540,int_stack+131595,int_stack+163890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309496,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+163890,int_stack+174092,int_stack+173651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+285011,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+163890,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+586589,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+133800,int_stack+566044,int_stack+132540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+587534,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+133800,int_stack+165240, 0.0, zero_stack, 1.0, int_stack+310846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+112800,int_stack+102300,int_stack+166740, 0.0, zero_stack, 1.0, int_stack+312346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+566044,int_stack+174960,int_stack+174680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202060,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+175380,int_stack+174960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202620,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+580534,int_stack+131595,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615098,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+566044,int_stack+175968,int_stack+175380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203460,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+566044,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615938,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+163890,int_stack+264472,int_stack+580534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+617198,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+166690,int_stack+163890,int_stack+133800, 0.0, zero_stack, 1.0, int_stack+322024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+923146,int_stack+166690,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+324124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+176852,int_stack+176752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280786, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+305626,int_stack+177002,int_stack+176852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280886, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102600,int_stack+305626,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582214, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103200,int_stack+177212,int_stack+177002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+281036, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103830,int_stack+103200,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306076, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104730,int_stack+103830,int_stack+102600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582514, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+177642,int_stack+177492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282266, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+177867,int_stack+177642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282416, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102975,int_stack+102300,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583114, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+105730,int_stack+178182,int_stack+177867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282641, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+163890,int_stack+105730,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585014, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+165240,int_stack+163890,int_stack+102975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585689, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+166740,int_stack+165240,int_stack+104730, 0.0, zero_stack, 1.0, int_stack+342124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+163890,int_stack+178812,int_stack+178602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284486, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+179127,int_stack+178812, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284696, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103245,int_stack+102300,int_stack+163890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309496, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+163890,int_stack+179568,int_stack+179127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+285011, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+163890,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+586589, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104505,int_stack+566044,int_stack+103245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+587534, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+104505,int_stack+165240, 0.0, zero_stack, 1.0, int_stack+346319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+169740,int_stack+131595,int_stack+166740, 0.0, zero_stack, 1.0, int_stack+347819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+566044,int_stack+180436,int_stack+180156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202060, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+180856,int_stack+180436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202620, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+580534,int_stack+102300,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615098, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+566044,int_stack+181444,int_stack+180856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203460, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+566044,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615938, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+163890,int_stack+264472,int_stack+580534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+617198, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+175740,int_stack+163890,int_stack+104505, 0.0, zero_stack, 1.0, int_stack+364841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+932146,int_stack+175740,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+366941, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+182328,int_stack+182228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280786, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+305626,int_stack+182478,int_stack+182328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280886, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+131895,int_stack+305626,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582214, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132495,int_stack+182688,int_stack+182478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+281036, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+133125,int_stack+132495,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306076, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+134025,int_stack+133125,int_stack+131895, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582514, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+183118,int_stack+182968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282266, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+183343,int_stack+183118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282416, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132270,int_stack+131595,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583114, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+135025,int_stack+183658,int_stack+183343, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282641, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+175740,int_stack+135025,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585014, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+177090,int_stack+175740,int_stack+132270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585689, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+178590,int_stack+177090,int_stack+134025, 0.0, zero_stack, 1.0, int_stack+382001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+175740,int_stack+184288,int_stack+184078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284486, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+184603,int_stack+184288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284696, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132540,int_stack+131595,int_stack+175740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309496, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+175740,int_stack+185044,int_stack+184603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+285011, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+175740,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+586589, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+133800,int_stack+566044,int_stack+132540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+587534, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+133800,int_stack+177090, 0.0, zero_stack, 1.0, int_stack+386196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+941146,int_stack+102300,int_stack+178590, 0.0, zero_stack, 1.0, int_stack+387696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+566044,int_stack+185912,int_stack+185632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202060, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+186332,int_stack+185912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202620, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+580534,int_stack+131595,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615098, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+566044,int_stack+186920,int_stack+186332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203460, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+566044,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615938, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+175740,int_stack+264472,int_stack+580534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+617198, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+178540,int_stack+175740,int_stack+133800, 0.0, zero_stack, 1.0, int_stack+403395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+947146,int_stack+178540,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+405495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+187804,int_stack+187704, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+305626,int_stack+187954,int_stack+187804, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+280886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102600,int_stack+305626,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103200,int_stack+188164,int_stack+187954, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+281036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103830,int_stack+103200,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+306076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104730,int_stack+103830,int_stack+102600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+582514, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+188594,int_stack+188444, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+188819,int_stack+188594, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102975,int_stack+102300,int_stack+305626, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+105730,int_stack+189134,int_stack+188819, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+282641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+175740,int_stack+105730,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+177090,int_stack+175740,int_stack+102975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+178590,int_stack+177090,int_stack+104730, 0.0, zero_stack, 1.0, int_stack+420555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+175740,int_stack+189764,int_stack+189554, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+190079,int_stack+189764, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+284696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103245,int_stack+102300,int_stack+175740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+309496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+175740,int_stack+190520,int_stack+190079, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+285011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+175740,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+586589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104505,int_stack+566044,int_stack+103245, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+587534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+104505,int_stack+177090, 0.0, zero_stack, 1.0, int_stack+424750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+181590,int_stack+131595,int_stack+178590, 0.0, zero_stack, 1.0, int_stack+426250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+566044,int_stack+191388,int_stack+191108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+191808,int_stack+191388, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+580534,int_stack+102300,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+566044,int_stack+192396,int_stack+191808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+566044,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+615938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+175740,int_stack+264472,int_stack+580534, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+617198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+956146,int_stack+175740,int_stack+104505, 0.0, zero_stack, 1.0, int_stack+443272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+962446,int_stack+956146,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+445372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+193280,int_stack+193180, 0.0, zero_stack, 1.0, int_stack+280786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+305626,int_stack+193430,int_stack+193280, 0.0, zero_stack, 1.0, int_stack+280886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+131895,int_stack+305626,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+582214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132495,int_stack+193640,int_stack+193430, 0.0, zero_stack, 1.0, int_stack+281036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+133125,int_stack+132495,int_stack+305626, 0.0, zero_stack, 1.0, int_stack+306076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+134025,int_stack+133125,int_stack+131895, 0.0, zero_stack, 1.0, int_stack+582514, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+194070,int_stack+193920, 0.0, zero_stack, 1.0, int_stack+282266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+194295,int_stack+194070, 0.0, zero_stack, 1.0, int_stack+282416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132270,int_stack+131595,int_stack+305626, 0.0, zero_stack, 1.0, int_stack+583114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+135025,int_stack+194610,int_stack+194295, 0.0, zero_stack, 1.0, int_stack+282641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+956146,int_stack+135025,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+585014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+957496,int_stack+956146,int_stack+132270, 0.0, zero_stack, 1.0, int_stack+585689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+958996,int_stack+957496,int_stack+134025, 0.0, zero_stack, 1.0, int_stack+460432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+956146,int_stack+195240,int_stack+195030, 0.0, zero_stack, 1.0, int_stack+284486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+195555,int_stack+195240, 0.0, zero_stack, 1.0, int_stack+284696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132540,int_stack+131595,int_stack+956146, 0.0, zero_stack, 1.0, int_stack+309496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+956146,int_stack+195996,int_stack+195555, 0.0, zero_stack, 1.0, int_stack+285011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+956146,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+586589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+133800,int_stack+566044,int_stack+132540, 0.0, zero_stack, 1.0, int_stack+587534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+133800,int_stack+957496, 0.0, zero_stack, 1.0, int_stack+464627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+187590,int_stack+102300,int_stack+958996, 0.0, zero_stack, 1.0, int_stack+466127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+566044,int_stack+196864,int_stack+196584, 0.0, zero_stack, 1.0, int_stack+202060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+197284,int_stack+196864, 0.0, zero_stack, 1.0, int_stack+202620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+580534,int_stack+131595,int_stack+566044, 0.0, zero_stack, 1.0, int_stack+615098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+566044,int_stack+197872,int_stack+197284, 0.0, zero_stack, 1.0, int_stack+203460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+566044,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+615938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+956146,int_stack+264472,int_stack+580534, 0.0, zero_stack, 1.0, int_stack+617198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+971446,int_stack+956146,int_stack+133800, 0.0, zero_stack, 1.0, int_stack+481826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+977746,int_stack+971446,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+483926, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+198756,int_stack+198656, 1.0, int_stack+280786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+305626,int_stack+198906,int_stack+198756, 1.0, int_stack+280886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102600,int_stack+305626,int_stack+102300, 1.0, int_stack+582214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103200,int_stack+199116,int_stack+198906, 1.0, int_stack+281036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103830,int_stack+103200,int_stack+305626, 1.0, int_stack+306076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+305626,int_stack+103830,int_stack+102600, 1.0, int_stack+582514, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+306626,int_stack+199546,int_stack+199396, 1.0, int_stack+282266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+307076,int_stack+199771,int_stack+199546, 1.0, int_stack+282416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+307751,int_stack+307076,int_stack+306626, 1.0, int_stack+583114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+102300,int_stack+200086,int_stack+199771, 1.0, int_stack+282641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103245,int_stack+102300,int_stack+307076, 1.0, int_stack+585014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104595,int_stack+103245,int_stack+307751, 1.0, int_stack+585689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+971446,int_stack+104595,int_stack+305626, 0.0, zero_stack, 1.0, int_stack+498986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+305626,int_stack+200716,int_stack+200506, 1.0, int_stack+284486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+306256,int_stack+201031,int_stack+200716, 1.0, int_stack+284696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+307201,int_stack+306256,int_stack+305626, 1.0, int_stack+309496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+308461,int_stack+201472,int_stack+201031, 1.0, int_stack+285011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+308461,int_stack+306256, 1.0, int_stack+586589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+308461,int_stack+566044,int_stack+307201, 1.0, int_stack+587534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+308461,int_stack+104595, 0.0, zero_stack, 1.0, int_stack+503181, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+956146,int_stack+131595,int_stack+971446, 0.0, zero_stack, 1.0, int_stack+504681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+971446,int_stack+203040,int_stack+202340, 1.0, int_stack+202060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+972286,int_stack+204048,int_stack+203040, 1.0, int_stack+202620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+972286,int_stack+971446, 1.0, int_stack+615098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+973546,int_stack+204636,int_stack+204048, 1.0, int_stack+203460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+973546,int_stack+972286, 1.0, int_stack+615938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+971446,int_stack+264472,int_stack+532980, 1.0, int_stack+617198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+193590,int_stack+971446,int_stack+308461, 0.0, zero_stack, 1.0, int_stack+520380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+986746,int_stack+193590,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+522480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+962146,int_stack+205520,int_stack+205420,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+205670,int_stack+205520,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+132045,int_stack+131595,int_stack+962146,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+132645,int_stack+205880,int_stack+205670,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+133275,int_stack+132645,int_stack+131595,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+134175,int_stack+133275,int_stack+132045,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+206310,int_stack+206160,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+132045,int_stack+206535,int_stack+206310,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+132720,int_stack+132045,int_stack+131595,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+193590,int_stack+206850,int_stack+206535,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+194535,int_stack+193590,int_stack+132045,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+195885,int_stack+194535,int_stack+132720,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+197385,int_stack+195885,int_stack+134175, 0.0, zero_stack, 1.0, int_stack+543760, 1.0, int_stack+584014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+193590,int_stack+207480,int_stack+207270,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+194220,int_stack+207795,int_stack+207480,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+131595,int_stack+194220,int_stack+193590,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+132855,int_stack+208236,int_stack+207795,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+132855,int_stack+194220,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+132855,int_stack+566044,int_stack+131595,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+132855,int_stack+195885, 0.0, zero_stack, 1.0, int_stack+546335, 1.0, int_stack+588884, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+200385,int_stack+102300,int_stack+197385, 0.0, zero_stack, 1.0, int_stack+563044, 1.0, int_stack+590384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+209104,int_stack+208824,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+566044,int_stack+209524,int_stack+209104,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+566044,int_stack+131595,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+193590,int_stack+210112,int_stack+209524,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+193590,int_stack+566044,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+193590,int_stack+264472,int_stack+532980,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+971446,int_stack+193590,int_stack+132855, 0.0, zero_stack, 1.0, int_stack+567934, 1.0, int_stack+593384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+995746,int_stack+971446,int_stack+102300, 0.0, zero_stack, 1.0, int_stack+570034, 1.0, int_stack+604598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+962146,int_stack+210996,int_stack+210896,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+211146,int_stack+210996,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+102750,int_stack+102300,int_stack+962146,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+103350,int_stack+211356,int_stack+211146,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+103980,int_stack+103350,int_stack+102300,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+104880,int_stack+103980,int_stack+102750,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+211786,int_stack+211636,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+212011,int_stack+211786,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+971446,int_stack+212326,int_stack+212011,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+972391,int_stack+971446,int_stack+102750,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+973741,int_stack+972391,int_stack+103425,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+193590,int_stack+973741,int_stack+104880, 0.0, zero_stack, 2.0, int_stack+584014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+971446,int_stack+212956,int_stack+212746,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+972076,int_stack+213271,int_stack+212956,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+102300,int_stack+972076,int_stack+971446,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+103560,int_stack+213712,int_stack+213271,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+103560,int_stack+972076,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+103560,int_stack+566044,int_stack+102300,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+103560,int_stack+973741, 0.0, zero_stack, 2.0, int_stack+588884, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+971446,int_stack+131595,int_stack+193590, 0.0, zero_stack, 2.0, int_stack+590384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+193590,int_stack+214580,int_stack+214300,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+215000,int_stack+214580,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+102300,int_stack+193590,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+193590,int_stack+215588,int_stack+215000,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+193590,int_stack+102300,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+193590,int_stack+264472,int_stack+532980,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+206385,int_stack+193590,int_stack+103560, 0.0, zero_stack, 2.0, int_stack+593384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+1004746,int_stack+206385,int_stack+131595, 0.0, zero_stack, 2.0, int_stack+604598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+977446,int_stack+217212,int_stack+217112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290406,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+217362,int_stack+217212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290506,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132045,int_stack+131595,int_stack+977446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+618878,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132645,int_stack+217572,int_stack+217362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290656,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+133275,int_stack+132645,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583564,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+134175,int_stack+133275,int_stack+132045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619178,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+219112,int_stack+218962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292996,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+132045,int_stack+219337,int_stack+219112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293146,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132720,int_stack+132045,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619778,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+206385,int_stack+219652,int_stack+219337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293371,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+207330,int_stack+206385,int_stack+132045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+620228,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+208680,int_stack+207330,int_stack+132720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+622308,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+210180,int_stack+208680,int_stack+134175, 1.0, int_stack+304626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+206385,int_stack+221836,int_stack+221626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295216,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+207015,int_stack+222151,int_stack+221836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295426,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+131595,int_stack+207015,int_stack+206385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+534660,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132855,int_stack+222592,int_stack+222151, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295741,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+132855,int_stack+207015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+623208,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+132855,int_stack+566044,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+624153,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+132855,int_stack+208680, 1.0, int_stack+310846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+213180,int_stack+102300,int_stack+210180, 1.0, int_stack+312346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+223460,int_stack+223180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+268546,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+566044,int_stack+223880,int_stack+223460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269106,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+566044,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+206385,int_stack+224468,int_stack+223880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269946,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+206385,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5040,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+206385,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+535290,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+219180,int_stack+206385,int_stack+132855, 1.0, int_stack+322024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+300396,int_stack+219180,int_stack+102300, 1.0, int_stack+324124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+977446,int_stack+226092,int_stack+225992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290406, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+226242,int_stack+226092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290506, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102750,int_stack+102300,int_stack+977446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+618878, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103350,int_stack+226452,int_stack+226242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290656, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103980,int_stack+103350,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583564, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104880,int_stack+103980,int_stack+102750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619178, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+227992,int_stack+227842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292996, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+228217,int_stack+227992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293146, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619778, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+219180,int_stack+228532,int_stack+228217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293371, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+220125,int_stack+219180,int_stack+102750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+620228, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+221475,int_stack+220125,int_stack+103425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+622308, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+222975,int_stack+221475,int_stack+104880, 1.0, int_stack+342124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+219180,int_stack+230716,int_stack+230506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295216, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+219810,int_stack+231031,int_stack+230716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295426, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102300,int_stack+219810,int_stack+219180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+534660, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103560,int_stack+231472,int_stack+231031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295741, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+103560,int_stack+219810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+623208, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+103560,int_stack+566044,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+624153, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+103560,int_stack+221475, 1.0, int_stack+346319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+225975,int_stack+131595,int_stack+222975, 1.0, int_stack+347819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+232340,int_stack+232060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+268546, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+566044,int_stack+232760,int_stack+232340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269106, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+566044,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+219180,int_stack+233348,int_stack+232760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269946, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+219180,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5040, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+219180,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+535290, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+206385,int_stack+219180,int_stack+103560, 1.0, int_stack+364841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+341464,int_stack+206385,int_stack+131595, 1.0, int_stack+366941, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+977446,int_stack+234972,int_stack+234872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290406, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+235122,int_stack+234972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290506, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132045,int_stack+131595,int_stack+977446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+618878, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132645,int_stack+235332,int_stack+235122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290656, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+133275,int_stack+132645,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583564, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+134175,int_stack+133275,int_stack+132045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619178, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+236872,int_stack+236722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292996, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+132045,int_stack+237097,int_stack+236872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293146, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132720,int_stack+132045,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619778, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+206385,int_stack+237412,int_stack+237097, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293371, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+207330,int_stack+206385,int_stack+132045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+620228, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+208680,int_stack+207330,int_stack+132720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+622308, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+210180,int_stack+208680,int_stack+134175, 1.0, int_stack+382001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+206385,int_stack+239596,int_stack+239386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295216, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+207015,int_stack+239911,int_stack+239596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295426, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+131595,int_stack+207015,int_stack+206385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+534660, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132855,int_stack+240352,int_stack+239911, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295741, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+132855,int_stack+207015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+623208, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+132855,int_stack+566044,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+624153, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+132855,int_stack+208680, 1.0, int_stack+386196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+219180,int_stack+102300,int_stack+210180, 1.0, int_stack+387696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+241220,int_stack+240940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+268546, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+566044,int_stack+241640,int_stack+241220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269106, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+566044,int_stack+131595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+206385,int_stack+242228,int_stack+241640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269946, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+206385,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5040, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+206385,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+535290, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+231975,int_stack+206385,int_stack+132855, 1.0, int_stack+403395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+361628,int_stack+231975,int_stack+102300, 1.0, int_stack+405495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+977446,int_stack+243852,int_stack+243752, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+244002,int_stack+243852, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102750,int_stack+102300,int_stack+977446, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+618878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103350,int_stack+244212,int_stack+244002, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+290656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103980,int_stack+103350,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+583564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104880,int_stack+103980,int_stack+102750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+245752,int_stack+245602, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+292996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+245977,int_stack+245752, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+619778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+231975,int_stack+246292,int_stack+245977, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+293371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+232920,int_stack+231975,int_stack+102750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+620228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234270,int_stack+232920,int_stack+103425, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+622308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+235770,int_stack+234270,int_stack+104880, 1.0, int_stack+420555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+231975,int_stack+248476,int_stack+248266, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+232605,int_stack+248791,int_stack+248476, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102300,int_stack+232605,int_stack+231975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+534660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103560,int_stack+249232,int_stack+248791, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295741, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+103560,int_stack+232605, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+623208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+103560,int_stack+566044,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+624153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+103560,int_stack+234270, 1.0, int_stack+424750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+238770,int_stack+131595,int_stack+235770, 1.0, int_stack+426250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+250100,int_stack+249820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+268546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+566044,int_stack+250520,int_stack+250100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+566044,int_stack+102300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+231975,int_stack+251108,int_stack+250520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+269946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+231975,int_stack+566044, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+231975,int_stack+264472,int_stack+532980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+535290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+244770,int_stack+231975,int_stack+103560, 1.0, int_stack+443272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+415995,int_stack+244770,int_stack+131595, 1.0, int_stack+445372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+977446,int_stack+252732,int_stack+252632, 0.0, zero_stack, 1.0, int_stack+290406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+252882,int_stack+252732, 0.0, zero_stack, 1.0, int_stack+290506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132045,int_stack+131595,int_stack+977446, 0.0, zero_stack, 1.0, int_stack+618878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132645,int_stack+253092,int_stack+252882, 0.0, zero_stack, 1.0, int_stack+290656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+133275,int_stack+132645,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+583564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+134175,int_stack+133275,int_stack+132045, 0.0, zero_stack, 1.0, int_stack+619178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+254632,int_stack+254482, 0.0, zero_stack, 1.0, int_stack+292996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+132045,int_stack+254857,int_stack+254632, 0.0, zero_stack, 1.0, int_stack+293146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132720,int_stack+132045,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+619778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+244770,int_stack+255172,int_stack+254857, 0.0, zero_stack, 1.0, int_stack+293371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+245715,int_stack+244770,int_stack+132045, 0.0, zero_stack, 1.0, int_stack+620228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+247065,int_stack+245715,int_stack+132720, 0.0, zero_stack, 1.0, int_stack+622308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+248565,int_stack+247065,int_stack+134175, 1.0, int_stack+460432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+244770,int_stack+257356,int_stack+257146, 0.0, zero_stack, 1.0, int_stack+295216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+245400,int_stack+257671,int_stack+257356, 0.0, zero_stack, 1.0, int_stack+295426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+131595,int_stack+245400,int_stack+244770, 0.0, zero_stack, 1.0, int_stack+534660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132855,int_stack+258112,int_stack+257671, 0.0, zero_stack, 1.0, int_stack+295741, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+132855,int_stack+245400, 0.0, zero_stack, 1.0, int_stack+623208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+132855,int_stack+566044,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+624153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+132855,int_stack+247065, 1.0, int_stack+464627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+251565,int_stack+102300,int_stack+248565, 1.0, int_stack+466127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+258980,int_stack+258700, 0.0, zero_stack, 1.0, int_stack+268546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+566044,int_stack+259400,int_stack+258980, 0.0, zero_stack, 1.0, int_stack+269106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+532980,int_stack+566044,int_stack+131595, 0.0, zero_stack, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+244770,int_stack+259988,int_stack+259400, 0.0, zero_stack, 1.0, int_stack+269946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+264472,int_stack+244770,int_stack+566044, 0.0, zero_stack, 1.0, int_stack+5040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+244770,int_stack+264472,int_stack+532980, 0.0, zero_stack, 1.0, int_stack+535290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+231975,int_stack+244770,int_stack+132855, 1.0, int_stack+481826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+440059,int_stack+231975,int_stack+102300, 1.0, int_stack+483926, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+977446,int_stack+261612,int_stack+261512, 1.0, int_stack+290406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102300,int_stack+261762,int_stack+261612, 1.0, int_stack+290506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+102750,int_stack+102300,int_stack+977446, 1.0, int_stack+618878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+103350,int_stack+261972,int_stack+261762, 1.0, int_stack+290656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103980,int_stack+103350,int_stack+102300, 1.0, int_stack+583564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104880,int_stack+103980,int_stack+102750, 1.0, int_stack+619178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+263512,int_stack+263362, 1.0, int_stack+292996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102750,int_stack+263737,int_stack+263512, 1.0, int_stack+293146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103425,int_stack+102750,int_stack+102300, 1.0, int_stack+619778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+231975,int_stack+264052,int_stack+263737, 1.0, int_stack+293371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+232920,int_stack+231975,int_stack+102750, 1.0, int_stack+620228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234270,int_stack+232920,int_stack+103425, 1.0, int_stack+622308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+235770,int_stack+234270,int_stack+104880, 1.0, int_stack+498986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+622308,int_stack+267202,int_stack+266992, 1.0, int_stack+295216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+231975,int_stack+267517,int_stack+267202, 1.0, int_stack+295426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+232920,int_stack+231975,int_stack+622308, 1.0, int_stack+534660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+102300,int_stack+267958,int_stack+267517, 1.0, int_stack+295741, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+102300,int_stack+231975, 1.0, int_stack+623208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+102300,int_stack+566044,int_stack+232920, 1.0, int_stack+624153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+131595,int_stack+102300,int_stack+234270, 1.0, int_stack+503181, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+244770,int_stack+131595,int_stack+235770, 1.0, int_stack+504681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+566044,int_stack+269526,int_stack+268826, 1.0, int_stack+268546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+231975,int_stack+270534,int_stack+269526, 1.0, int_stack+269106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+233235,int_stack+231975,int_stack+566044, 1.0, int_stack+334624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+566044,int_stack+271122,int_stack+270534, 1.0, int_stack+269946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234915,int_stack+566044,int_stack+231975, 1.0, int_stack+5040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+622308,int_stack+234915,int_stack+233235, 1.0, int_stack+535290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+0,int_stack+622308,int_stack+102300, 1.0, int_stack+520380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+257565,int_stack+0,int_stack+131595, 1.0, int_stack+522480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+977446,int_stack+272746,int_stack+272646,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+272896,int_stack+272746,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+132045,int_stack+131595,int_stack+977446,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+132645,int_stack+273106,int_stack+272896,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+133275,int_stack+132645,int_stack+131595,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+134175,int_stack+133275,int_stack+132045,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+131595,int_stack+274646,int_stack+274496,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+132045,int_stack+274871,int_stack+274646,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+132720,int_stack+132045,int_stack+131595,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+275186,int_stack+274871,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+945,int_stack+0,int_stack+132045,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+2295,int_stack+945,int_stack+132720,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+102300,int_stack+2295,int_stack+134175, 1.0, int_stack+543760, 0.0, zero_stack, 1.0, int_stack+621308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+277370,int_stack+277160,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+630,int_stack+277685,int_stack+277370,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+131595,int_stack+630,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+132855,int_stack+278126,int_stack+277685,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+566044,int_stack+132855,int_stack+630,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+132855,int_stack+566044,int_stack+131595,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+231975,int_stack+132855,int_stack+2295, 1.0, int_stack+546335, 0.0, zero_stack, 1.0, int_stack+625503, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+0,int_stack+231975,int_stack+102300, 1.0, int_stack+563044, 0.0, zero_stack, 1.0, int_stack+627003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+334624,int_stack+278994,int_stack+278714,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+131595,int_stack+279414,int_stack+278994,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+563044,int_stack+131595,int_stack+334624,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+564724,int_stack+280002,int_stack+279414,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+102300,int_stack+564724,int_stack+131595,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+564724,int_stack+102300,int_stack+563044,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+266565,int_stack+564724,int_stack+132855, 1.0, int_stack+567934, 0.0, zero_stack, 1.0, int_stack+630003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+494426,int_stack+266565,int_stack+231975, 1.0, int_stack+570034, 0.0, zero_stack, 1.0, int_stack+641262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6000,int_stack+281626,int_stack+281526,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+231975,int_stack+281776,int_stack+281626,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+232425,int_stack+231975,int_stack+6000,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+233025,int_stack+281986,int_stack+281776,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+233655,int_stack+233025,int_stack+231975,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+234555,int_stack+233655,int_stack+232425,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+231975,int_stack+283526,int_stack+283376,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+232425,int_stack+283751,int_stack+283526,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+233100,int_stack+232425,int_stack+231975,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+235555,int_stack+284066,int_stack+283751,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+236500,int_stack+235555,int_stack+232425,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+266565,int_stack+236500,int_stack+233100,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+235555,int_stack+266565,int_stack+234555, 1.0, int_stack+584014, 1.0, int_stack+621308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+231975,int_stack+286250,int_stack+286040,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+232605,int_stack+286565,int_stack+286250,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+233550,int_stack+232605,int_stack+231975,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+268065,int_stack+287006,int_stack+286565,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+269388,int_stack+268065,int_stack+232605,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+271278,int_stack+269388,int_stack+233550,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+102300,int_stack+271278,int_stack+266565, 1.0, int_stack+588884, 1.0, int_stack+625503, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+273378,int_stack+102300,int_stack+235555, 1.0, int_stack+590384, 1.0, int_stack+627003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+334624,int_stack+287874,int_stack+287594,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+266565,int_stack+288294,int_stack+287874,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+267825,int_stack+266565,int_stack+334624,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+269505,int_stack+288882,int_stack+288294,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+231975,int_stack+269505,int_stack+266565,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+234495,int_stack+231975,int_stack+267825,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+279378,int_stack+234495,int_stack+271278, 1.0, int_stack+593384, 1.0, int_stack+630003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+563044,int_stack+279378,int_stack+102300, 1.0, int_stack+604598, 1.0, int_stack+641262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6000,int_stack+291246,int_stack+291146,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+604598,int_stack+291396,int_stack+291246,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+605048,int_stack+604598,int_stack+6000,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+605648,int_stack+291606,int_stack+291396,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+606278,int_stack+605648,int_stack+604598,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+607178,int_stack+606278,int_stack+605048,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+604598,int_stack+294256,int_stack+294106,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+605048,int_stack+294481,int_stack+294256,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+605723,int_stack+605048,int_stack+604598,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+102300,int_stack+294796,int_stack+294481,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+103245,int_stack+102300,int_stack+605048,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+104595,int_stack+103245,int_stack+605723,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+279378,int_stack+104595,int_stack+607178, 2.0, int_stack+621308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+102300,int_stack+296980,int_stack+296770,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+102930,int_stack+297295,int_stack+296980,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+604598,int_stack+102930,int_stack+102300,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+605858,int_stack+297736,int_stack+297295,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+607181,int_stack+605858,int_stack+102930,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+102300,int_stack+607181,int_stack+604598,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+604598,int_stack+102300,int_stack+104595, 2.0, int_stack+625503, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+282378,int_stack+604598,int_stack+279378, 2.0, int_stack+627003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+334624,int_stack+298604,int_stack+298324,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+279378,int_stack+299024,int_stack+298604,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+280638,int_stack+279378,int_stack+334624,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+104400,int_stack+299612,int_stack+299024,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+288378,int_stack+104400,int_stack+279378,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+290898,int_stack+288378,int_stack+280638,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+293698,int_stack+290898,int_stack+102300, 2.0, int_stack+630003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+615098,int_stack+293698,int_stack+604598, 2.0, int_stack+641262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+288378,int_stack+352628,int_stack+328624,100);
     Libderiv->ABCD[11] = int_stack + 288378;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+580534,int_stack+392505,int_stack+371441,100);
     Libderiv->ABCD[10] = int_stack + 580534;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+350464,int_stack+431059,int_stack+409995,100);
     Libderiv->ABCD[9] = int_stack + 350464;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+455872,int_stack+470936,int_stack+449872,100);
     Libderiv->ABCD[8] = int_stack + 455872;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+465872,int_stack+509490,int_stack+488426,100);
     Libderiv->ABCD[7] = int_stack + 465872;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+503426,int_stack+548044,int_stack+526980,100);
     Libderiv->ABCD[6] = int_stack + 503426;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+513426,int_stack+595598,int_stack+574534, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+557044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 513426;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+590534,int_stack+632262,int_stack+609098, 0.0, zero_stack, 1.0, int_stack+557044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 590534;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+624098,int_stack+665146,int_stack+645762, 1.0, int_stack+557044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 624098;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+634098,int_stack+674146,int_stack+536970,100);
     Libderiv->ABCD[155] = int_stack + 634098;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+532980,int_stack+683146,int_stack+6300,100);
     Libderiv->ABCD[143] = int_stack + 532980;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+542980,int_stack+692146,int_stack+12300,100);
     Libderiv->ABCD[142] = int_stack + 542980;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+6000,int_stack+701146,int_stack+18300,100);
     Libderiv->ABCD[131] = int_stack + 6000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+552980,int_stack+710146,int_stack+24300,100);
     Libderiv->ABCD[130] = int_stack + 552980;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+16000,int_stack+719146,int_stack+30300,100);
     Libderiv->ABCD[129] = int_stack + 16000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+26000,int_stack+728146,int_stack+36300,100);
     Libderiv->ABCD[119] = int_stack + 26000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+664557,int_stack+737146,int_stack+42300,100);
     Libderiv->ABCD[118] = int_stack + 664557;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+36000,int_stack+746146,int_stack+48300,100);
     Libderiv->ABCD[117] = int_stack + 36000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+674557,int_stack+755146,int_stack+54300,100);
     Libderiv->ABCD[116] = int_stack + 674557;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+46000,int_stack+764146,int_stack+60300,100);
     Libderiv->ABCD[107] = int_stack + 46000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+56000,int_stack+779146,int_stack+773146,100);
     Libderiv->ABCD[106] = int_stack + 56000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+684557,int_stack+788146,int_stack+66300,100);
     Libderiv->ABCD[105] = int_stack + 684557;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+694557,int_stack+797146,int_stack+72300,100);
     Libderiv->ABCD[104] = int_stack + 694557;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+66000,int_stack+806146,int_stack+78300,100);
     Libderiv->ABCD[103] = int_stack + 66000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+704557,int_stack+815146,int_stack+84300,100);
     Libderiv->ABCD[95] = int_stack + 704557;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+76000,int_stack+824146,int_stack+90300,100);
     Libderiv->ABCD[94] = int_stack + 76000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+86000,int_stack+833146,int_stack+96300,100);
     Libderiv->ABCD[93] = int_stack + 86000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+96000,int_stack+842146,int_stack+106800,100);
     Libderiv->ABCD[92] = int_stack + 96000;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+714557,int_stack+851146,int_stack+658557,100);
     Libderiv->ABCD[91] = int_stack + 714557;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+724557,int_stack+860146,int_stack+119595,100);
     Libderiv->ABCD[90] = int_stack + 724557;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+734557,int_stack+869146,int_stack+125595, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+328624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[47] = int_stack + 734557;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+118800,int_stack+878146,int_stack+651762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+371441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[46] = int_stack + 118800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+651762,int_stack+887146,int_stack+335464, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+409995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[45] = int_stack + 651762;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+744557,int_stack+896146,int_stack+136095, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+449872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[44] = int_stack + 744557;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+128800,int_stack+905146,int_stack+148890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+488426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[43] = int_stack + 128800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+754557,int_stack+154890,int_stack+315346, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+526980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[42] = int_stack + 754557;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+309396,int_stack+914146,int_stack+142095, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+574534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[38] = int_stack + 309396;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+138800,int_stack+923146,int_stack+112800, 0.0, zero_stack, 1.0, int_stack+328624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[35] = int_stack + 138800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+106000,int_stack+932146,int_stack+169740, 0.0, zero_stack, 1.0, int_stack+371441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[34] = int_stack + 106000;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+148800,int_stack+947146,int_stack+941146, 0.0, zero_stack, 1.0, int_stack+409995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[33] = int_stack + 148800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+158800,int_stack+962446,int_stack+181590, 0.0, zero_stack, 1.0, int_stack+449872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[32] = int_stack + 158800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+168800,int_stack+977746,int_stack+187590, 0.0, zero_stack, 1.0, int_stack+488426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[31] = int_stack + 168800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+178800,int_stack+986746,int_stack+956146, 0.0, zero_stack, 1.0, int_stack+526980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[30] = int_stack + 178800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+977446,int_stack+995746,int_stack+200385, 0.0, zero_stack, 1.0, int_stack+574534, 1.0, int_stack+609098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[26] = int_stack + 977446;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+987446,int_stack+1004746,int_stack+971446, 0.0, zero_stack, 2.0, int_stack+609098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[25] = int_stack + 987446;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+997446,int_stack+300396,int_stack+213180, 1.0, int_stack+328624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[23] = int_stack + 997446;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+298378,int_stack+341464,int_stack+225975, 1.0, int_stack+371441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[22] = int_stack + 298378;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+225180,int_stack+361628,int_stack+219180, 1.0, int_stack+409995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[21] = int_stack + 225180;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+188800,int_stack+415995,int_stack+238770, 1.0, int_stack+449872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[20] = int_stack + 188800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+198800,int_stack+440059,int_stack+251565, 1.0, int_stack+488426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[19] = int_stack + 198800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+208800,int_stack+257565,int_stack+244770, 1.0, int_stack+526980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[18] = int_stack + 208800;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+235180,int_stack+494426,int_stack+0, 1.0, int_stack+574534, 0.0, zero_stack, 1.0, int_stack+645762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[14] = int_stack + 235180;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+245180,int_stack+563044,int_stack+273378, 1.0, int_stack+609098, 1.0, int_stack+645762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[13] = int_stack + 245180;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+562980,int_stack+615098,int_stack+282378, 2.0, int_stack+645762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[12] = int_stack + 562980;

}
