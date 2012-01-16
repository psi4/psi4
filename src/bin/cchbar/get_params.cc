/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void get_params(Options &options)
{
  params.memory = Process::environment.get_memory();

  /* compute the Tamplitude equation matrix elements (usually 0) */
//  params.Tamplitude = 0;
//  errcod = ip_boolean("TAMPLITUDE", &(params.Tamplitude),0);
  params.Tamplitude = options.get_bool("T_AMPS");

//  params.cachelev = 2;
//  errcod = ip_data("CACHELEVEL", "%d", &(params.cachelev),0);
  params.cachelev = options.get_int("CACHELEVEL");

//  params.print = 0;
//  errcod = ip_data("PRINT", "%d", &(params.print),0);
  params.print = options.get_int("PRINT");

//  errcod = ip_string("WFN", &(params.wfn), 0);
  params.wfn = options.get_str("WFN");

//  params.dertype = 0;
  std::string junk = options.get_str("DERTYPE");
  if(junk == "NONE") params.dertype = 0;
  else if(junk == "FIRST") params.dertype = 1;
  else if(junk == "RESPONSE") params.dertype = 3; /* linear response */
  else {
//    printf("Invalid value of input keyword DERTYPE: %s\n", junk);
//    return PSI_RETURN_FAILURE;
      throw PsiException("CCHBAR: Invalid value of input keyword DERTYPE",__FILE__,__LINE__);
  }

  /* Should we use the minimal-disk algorithm for Wabei?  It's VERY slow! */
//  params.wabei_lowdisk = 0;
//  errcod = ip_boolean("WABEI_LOWDISK", &params.wabei_lowdisk, 0);
  params.wabei_lowdisk = options.get_bool("WABEI_LOWDISK");
}

}} // namespace psi::cchbar
