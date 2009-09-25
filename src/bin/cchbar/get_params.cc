/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void get_params()
{
  int errcod;
  char *junk;

  fndcor(&(params.memory),infile,outfile);

  /* compute the Tamplitude equation matrix elements (usually 0) */
  params.Tamplitude = 0;
  errcod = ip_boolean("TAMPLITUDE", &(params.Tamplitude),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);

  errcod = ip_string("WFN", &(params.wfn), 0);

  params.dertype = 0;
  if(ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if(errcod != IPE_OK) params.dertype = 0;
    else if(!strcmp(junk,"NONE")) params.dertype = 0;
    else if(!strcmp(junk,"FIRST")) params.dertype = 1;
    else if(!strcmp(junk,"RESPONSE")) params.dertype = 3; /* linear response */
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }

  /* Should we use the minimal-disk algorithm for Wabei?  It's VERY slow! */
  params.wabei_lowdisk = 0;
  errcod = ip_boolean("WABEI_LOWDISK", &params.wabei_lowdisk, 0);
}

}} // namespace psi::cchbar
