#include <cstdio>
#include <cstring>
#include <ccfiles.h>

/* A simple function to convert coupled cluster code file names to
** unit numbers, since I can only seem to remember the former.  
**
** -TDC, April 2001
*/ 
namespace psi{
  namespace tocprint{

  int cc2unit(char *cc)
  {

   if(!strcmp(cc,"INFO") || !strcmp(cc,"info")) return CC_INFO;
   else if(!strcmp(cc, "OEI") || !strcmp(cc, "oei")) return CC_OEI;
  
   else if (!strcmp(cc, "AINTS") || !strcmp(cc, "aints")) return CC_AINTS;
   else if (!strcmp(cc, "BINTS") || !strcmp(cc, "bints")) return CC_BINTS;
   else if (!strcmp(cc, "CINTS") || !strcmp(cc, "cints")) return CC_CINTS;
   else if (!strcmp(cc, "DINTS") || !strcmp(cc, "dints")) return CC_DINTS;
   else if (!strcmp(cc, "EINTS") || !strcmp(cc, "eints")) return CC_EINTS;
   else if (!strcmp(cc, "FINTS") || !strcmp(cc, "fints")) return CC_FINTS;
  
   else if (!strcmp(cc, "DENOM") || !strcmp(cc, "denom")) return CC_DENOM;
   else if (!strcmp(cc, "TAMPS") || !strcmp(cc, "tamps")) return CC_TAMPS;
   else if (!strcmp(cc, "LAMPS") || !strcmp(cc, "lamps")) return CC_LAMPS;
   else if (!strcmp(cc, "LAMBDA") || !strcmp(cc, "lambda")) return CC_LAMBDA;
   else if (!strcmp(cc, "RAMPS") || !strcmp(cc, "ramps")) return CC_RAMPS;
   else if (!strcmp(cc, "HBAR") || !strcmp(cc, "hbar")) return CC_HBAR;
   else if (!strcmp(cc, "GAMMA") || !strcmp(cc, "gamma")) return CC_GAMMA;
   else if (!strcmp(cc, "MISC") || !strcmp(cc, "misc")) return CC_MISC;
   else if (!strcmp(cc, "GLG") || !strcmp(cc, "glg")) return CC_GLG;
   else if (!strcmp(cc, "GL") || !strcmp(cc, "gl")) return CC_GL;
   else if (!strcmp(cc, "GR") || !strcmp(cc, "gr")) return CC_GR;
  
   else if (!strcmp(cc, "TMP") || !strcmp(cc, "tmp")) return CC_TMP;
   else if (!strcmp(cc, "TMP0") || !strcmp(cc, "tmp0")) return CC_TMP0;
   else if (!strcmp(cc, "TMP1") || !strcmp(cc, "tmp1")) return CC_TMP1;
   else if (!strcmp(cc, "TMP2") || !strcmp(cc, "tmp2")) return CC_TMP2;
   else if (!strcmp(cc, "TMP3") || !strcmp(cc, "tmp3")) return CC_TMP3;
   else if (!strcmp(cc, "TMP4") || !strcmp(cc, "tmp4")) return CC_TMP4;
   else if (!strcmp(cc, "TMP5") || !strcmp(cc, "tmp5")) return CC_TMP5;
   else if (!strcmp(cc, "TMP6") || !strcmp(cc, "tmp6")) return CC_TMP6;
   else if (!strcmp(cc, "TMP7") || !strcmp(cc, "tmp7")) return CC_TMP7;
   else if (!strcmp(cc, "TMP8") || !strcmp(cc, "tmp8")) return CC_TMP8;
   else if (!strcmp(cc, "TMP9") || !strcmp(cc, "tmp9")) return CC_TMP9;
   else if (!strcmp(cc, "TMP10") || !strcmp(cc, "tmp10")) return CC_TMP10;
   else if (!strcmp(cc, "TMP11") || !strcmp(cc, "tmp11")) return CC_TMP11;
  
   else { return 0; }
  }

  } /* Namespace tocprint */
} /* Namespace psi */
