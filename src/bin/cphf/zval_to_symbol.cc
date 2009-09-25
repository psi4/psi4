/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
/*** ZVAL_TO_SYMBOL() return atom symbol ***/ 

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace psi { namespace cphf {

void zval_to_symbol(double zval, char *sym) {
  int z;
  z = (int) zval;

  if (z==0) strcpy(sym,"G");
  else if (z==1) strcpy(sym,"H"); 
  else if (z==2) strcpy(sym,"HE"); 
  else if (z==3) strcpy(sym,"LI"); 
  else if (z==4) strcpy(sym,"BE"); 
  else if (z==5) strcpy(sym,"B"); 
  else if (z==6) strcpy(sym,"C"); 
  else if (z==7) strcpy(sym,"N"); 
  else if (z==8) strcpy(sym,"O"); 
  else if (z==9) strcpy(sym,"F"); 
  else if (z==10) strcpy(sym,"NE"); 
  else if (z==11) strcpy(sym,"NA"); 
  else if (z==12) strcpy(sym,"MG"); 
  else if (z==13) strcpy(sym,"AL"); 
  else if (z==14) strcpy(sym,"SI"); 
  else if (z==15) strcpy(sym,"P"); 
  else if (z==16) strcpy(sym,"S"); 
  else if (z==17) strcpy(sym,"CL"); 
  else if (z==18) strcpy(sym,"AR"); 
  else if (z==19) strcpy(sym,"K"); 
  else if (z==20) strcpy(sym,"CA"); 
  else if (z==21) strcpy(sym,"SC"); 
  else if (z==22) strcpy(sym,"TI"); 
  else if (z==23) strcpy(sym,"V"); 
  else if (z==24) strcpy(sym,"CR"); 
  else if (z==25) strcpy(sym,"MN"); 
  else if (z==26) strcpy(sym,"FE"); 
  else if (z==27) strcpy(sym,"CO"); 
  else if (z==28) strcpy(sym,"NI"); 
  else if (z==29) strcpy(sym,"CU"); 
  else if (z==30) strcpy(sym,"ZN"); 
  else if (z==31) strcpy(sym,"GA"); 
  else if (z==32) strcpy(sym,"GE"); 
  else if (z==33) strcpy(sym,"AS"); 
  else if (z==34) strcpy(sym,"SE"); 
  else if (z==35) strcpy(sym,"BR"); 
  else if (z==36) strcpy(sym,"KR"); 
  else if (z==37) strcpy(sym,"RB"); 
  else if (z==38) strcpy(sym,"SR"); 
  else if (z==39) strcpy(sym,"Y"); 
  else if (z==40) strcpy(sym,"ZR"); 
  else if (z==41) strcpy(sym,"NB");
  else if (z==42) strcpy(sym,"MO");
  else if (z==43) strcpy(sym,"TC");
  else if (z==44) strcpy(sym,"RU");
  else if (z==45) strcpy(sym,"RH");
  else if (z==46) strcpy(sym,"PD");
  else if (z==47) strcpy(sym,"AG");
  else if (z==48) strcpy(sym,"CD");
  else if (z==49) strcpy(sym,"IN"); 
  else if (z==50) strcpy(sym,"SN");
  else if (z==51) strcpy(sym,"SB");
  else if (z==52) strcpy(sym,"TE");
  else if (z==53) strcpy(sym,"I");
  else if (z==54) strcpy(sym,"XE");
  else if (z==55) strcpy(sym,"CS");
  else if (z==56) strcpy(sym,"BA");
  else if (z==57) strcpy(sym,"LA");
  else if (z==58) strcpy(sym,"CE");
  else if (z==59) strcpy(sym,"PR");
  else if (z==60) strcpy(sym,"ND");
  else if (z==61) strcpy(sym,"PM");
  else if (z==62) strcpy(sym,"SM");
  else if (z==63) strcpy(sym,"EU");
  else if (z==64) strcpy(sym,"GD");
  else if (z==65) strcpy(sym,"TB");
  else if (z==66) strcpy(sym,"DY");
  else if (z==67) strcpy(sym,"HO");
  else if (z==68) strcpy(sym,"ER");
  else if (z==69) strcpy(sym,"TM");
  else if (z==70) strcpy(sym,"TY");
  else if (z==71) strcpy(sym,"LU");
  else if (z==72) strcpy(sym,"HF");
  else if (z==73) strcpy(sym,"TA");
  else if (z==74) strcpy(sym,"W");
  else if (z==75) strcpy(sym,"RE");
  else if (z==76) strcpy(sym,"OS");
  else if (z==77) strcpy(sym,"IR");
  else if (z==78) strcpy(sym,"PT");
  else if (z==79) strcpy(sym,"AU");
  else if (z==80) strcpy(sym,"HG");
  else if (z==81) strcpy(sym,"TL");
  else if (z==82) strcpy(sym,"PB");
  else if (z==83) strcpy(sym,"BI");
  else if (z==84) strcpy(sym,"PO");
  else if (z==85) strcpy(sym,"AT");
  else if (z==86) strcpy(sym,"RN");
  else if (z==87) strcpy(sym,"FR");
  else if (z==88) strcpy(sym,"RA");
  else if (z==89) strcpy(sym,"AC");
  else if (z==90) strcpy(sym,"TH");
  else if (z==91) strcpy(sym,"PA");
  else if (z==92) strcpy(sym,"U");
  else if (z==93) strcpy(sym,"NP");
  else if (z==94) strcpy(sym,"PU");
  else if (z==95) strcpy(sym,"AM");
  else if (z==96) strcpy(sym,"CM");
  else if (z==97) strcpy(sym,"BK");
  else if (z==98) strcpy(sym,"CF");
  else if (z==99) strcpy(sym,"ES");
  else if (z==100) strcpy(sym,"FM");
  else if (z==101) strcpy(sym,"MD");
  else if (z==102) strcpy(sym,"NO");
  else if (z==103) strcpy(sym,"UNQ");
  else if (z==104) strcpy(sym,"UNP");
  else if (z==105) strcpy(sym,"UNH");
  else if (z==106) strcpy(sym,"UNS");
  return ;
}

}} // namespace psi::cphf
