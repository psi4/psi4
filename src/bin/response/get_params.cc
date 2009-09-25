/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

void get_params()
{
  int errcod, ref, count;
  char *junk, units[20];

  params.print = 1;
  errcod = ip_data("PRINT","%d",&(params.print),0);

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
    fprintf(outfile, "Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", 
	   ref, params.ref);
    fprintf(outfile, "Is this what you want to do?\n");
    params.ref = ref;
  }

  /* grab the field strength from input -- a few different units are converted to E_h */
  params.omega = 0.0; /* static polarizability by default */
  if(ip_exist("OMEGA",0)) {
    errcod = ip_count("OMEGA", &count, 0);

    if(errcod == IPE_NOT_AN_ARRAY)  /* assume Hartrees */
      errcod = ip_data("OMEGA", "%lf", &(params.omega), 0);

    else if(count == 2) {

      errcod = ip_data("OMEGA", "%lf", &(params.omega), 1, 0);
      errcod = ip_data("OMEGA", "%s", units, 1, 1);

      for(junk = units; *junk != '\0'; junk++)
	if(*junk>='a' && *junk <= 'z') *junk += 'A' - 'a';

      if(!strcmp(units, "HZ")) params.omega *= _h / _hartree2J;
      else if(!strcmp(units, "NM")) params.omega = (_c*_h/_hartree2J) / (params.omega * 1e-9);
      else if(!strcmp(units, "EV")) params.omega /= _hartree2ev;
    }
  }

  if(ip_exist("PROPERTY",0)) {
    errcod = ip_string("PROPERTY", &(params.prop), 0);
    if(strcmp(params.prop,"POLARIZABILITY") && strcmp(params.prop,"ROTATION")
       && strcmp(params.prop,"ROA") && strcmp(params.prop,"ALL")) {
      fprintf(outfile, "Invalid choice of resp. property: %s\n", params.prop);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.prop = strdup("POLARIZABILITY");

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tProperty        = %s\n", params.prop);
  fprintf(outfile, "\tReference wfn   =    %5s\n",
           (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tCache Level     =    %1d\n", params.cachelev);
  fprintf(outfile, "\tPrint Level     =    %1d\n",  params.print);
  fprintf(outfile, "\tApplied field   =    %5.3f E_h\n", params.omega);
}


}} // namespace psi::response
