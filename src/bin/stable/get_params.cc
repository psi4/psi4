/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

void get_params()
{
  int errcod, tol, ref;
  char *junk;

  params.print_lvl = 1;
  errcod = ip_data("PRINT","%d",&(params.print_lvl),0);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
  params.cachelev = 0;

  errcod = ip_string("REFERENCE", &(junk),0);
  /* if no reference is given, assume rhf */
  if (errcod != IPE_OK) {
    ref = 0;
  }
  else {
    if(!strcmp(junk, "RHF")) ref = 0;
    else if(!strcmp(junk, "ROHF")) ref = 1;
    else if(!strcmp(junk, "UHF")) ref = 2;
    else { 
      printf("Invalid value of input keyword REFERENCE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }

  /* Make sure the value of ref matches that from CC_INFO */
  if(params.ref != ref) {
    fprintf(outfile, "Value of REFERENCE from input.dat (%1d) and " 
            "CC_INFO (%1d) do not match!\n", ref, params.ref);
    fprintf(outfile, "Is this what you want to do?\n");
    params.ref = ref;
  }

  params.follow_instab = 0;
  errcod = ip_boolean("FOLLOW",&(params.follow_instab),0);
  if (params.follow_instab == 1 && params.ref != 2) {
    fprintf(outfile, "\nCan't follow instabilities unless REFERENCE is UHF\n");
    fprintf(outfile, "Instability following turned off.\n");
    params.follow_instab = 0;
  }

  params.num_evecs_print = 0;
  if (params.print_lvl > 2) params.num_evecs_print = 5;
  errcod = ip_data("NUM_EVECS_PRINT","%d",&(params.num_evecs_print),0);

  params.rotation_method = 0;
  errcod = ip_data("ROTATION_METHOD","%d",&(params.rotation_method),0);

  params.scale = 0.5;
  errcod = ip_data("SCALE","%lf",&(params.scale),0);

}


}} // namespace psi::stable
