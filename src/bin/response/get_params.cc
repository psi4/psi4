/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
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
  int ref, count;
  std::string units;
  std::string junk;

  params.print = options["PRINT"].to_integer();

  params.memory = module.get_memory();

  /*  Clearly, this isn't doing anything right now
  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
  */
  params.cachelev = 0;


  /* if no reference is given, assume rhf */
  ref = 0;
  junk = options.get_str("REFERENCE");
  if(junk == "RHF") ref = 0;
  else if(junk == "ROHF") ref = 1;
  else if(junk == "UHF") ref = 2;
  else {
      throw PsiException("Invalid value of the keyword REFERENCE", __FILE__, __LINE__);
  }

  /* Make sure the value of ref matches that from CC_INFO */
  if(params.ref != ref) {
    fprintf(outfile, "Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n",
           ref, params.ref);
    fprintf(outfile, "Is this what you want to do?\n");
    params.ref = ref;
  }

  /* grab the field frequencies from input -- a few different units are converted to E_h */
  count = options["OMEGA"].size();
  if(count == 1) { /* assume Hartrees and only one frequency */
    params.nomega = 1;
    params.omega = init_array(1);
    params.omega[0] = options["OMEGA"][0].to_double();
  }
  else if(count >= 2) {
    params.nomega = count-1;
    params.omega = init_array(params.nomega);

    units = options["OMEGA"][count-1].to_string();
    for(int i=0; i < count-1; i++) {
      params.omega[i] = options["OMEGA"][i].to_double();

      if(units == "HZ") params.omega[i] *= _h / _hartree2J;
      else if(units == "AU") 1; /* do nothing */
      else if(units == "NM") params.omega[i] = (_c*_h*1e9)/(params.omega[i]*_hartree2J);
      else if(units == "EV") params.omega[i] /= _hartree2ev;
      else {
        fprintf(outfile, "\n\tError in unit for input field frequencies.  Must use one of:\n");
        fprintf(outfile,   "\tau, hz, nm, or ev.\n");
        throw PsiException("Failure in response involving the OMEGA option.", __FILE__, __LINE__);
      }
    }
  }
  else {
    fprintf(outfile, "\n\tError reading input field frequencies.  Please use the format:\n");
    fprintf(outfile,   "\t  omega = (value1 value2 ... units)\n");
    fprintf(outfile,   "\twhere units = hartrees, hz, nm, or ev.\n");
    throw PsiException("Failure in response involving the OMEGA option.", __FILE__, __LINE__);
  }

  if(options["PROPERTY"].has_changed()) {
    params.prop = options.get_str("PROPERTY");
    if(!(params.prop == "POLARIZABILITY") && !(params.prop == "ROTATION")
       && !(params.prop == "ROA") && !(params.prop,"ALL")) {
      throw PsiException("Invalid choice of response property", __FILE__, __LINE__);
    }
  }

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tProperty        = ");
  fprintf(outfile, params.prop.c_str());
  fprintf(outfile, "\n\tReference wfn   =    %5s\n",
           (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tCache Level     =    %1d\n", params.cachelev);
  fprintf(outfile, "\tPrint Level     =    %1d\n",  params.print);
  fprintf(outfile, "\tApplied field   =    %5.3f E_h\n", params.omega);
}


}} // namespace psi::response
