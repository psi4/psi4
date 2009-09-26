/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

void get_params()
{
  int tol, ref;
  char *junk;

  params.print_lvl = options.get_int("PRINT");

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = options.get_int("CACHELEV"); 
  params.cachelev = 0;

  junk = options.get_str("REFERENCE");

  if(!strcmp(junk, "RHF")) ref = 0;
  else if(!strcmp(junk, "ROHF")) ref = 1;
  else if(!strcmp(junk, "UHF")) ref = 2;
  else { 
    printf("Invalid value of input keyword REFERENCE: %s\n", junk);
    throw PsiException("stable error", __FILE__, __LINE__);
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
  params.follow_instab = options.get_bool("FOLLOW");

  if (params.follow_instab == 1 && params.ref != 2) {
    fprintf(outfile, "\nCan't follow instabilities unless REFERENCE is UHF\n");
    fprintf(outfile, "Instability following turned off.\n");
    params.follow_instab = 0;
  }

  if (params.print_lvl > 2) params.num_evecs_print = 5;
  params.num_evecs_print = options.get_int("NUM_EVECS_PRINT");

  params.rotation_method = options.get_int("ROTATION_METHOD");

  params.scale = 0.5;
  params.scale = options.get_double("SCALE");

}


}} // namespace psi::stable
