/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <liboptions/liboptions.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void get_params( Options& options)
{
  int errcod, tol;
  std::string junk;

  params.wfn  = options.get_str("WFN");

  junk = options.get_str("REFERENCE");
  if(junk == "RHF") params.ref = 0;
  else if(junk == "ROHF") params.ref = 1;
  else if(junk == "UHF") params.ref = 2;
  else {
    printf("Invalid value of input keyword REFERENCE: %s\n", junk.c_str());
    throw PsiException("ccdensity: error", __FILE__, __LINE__);
  }

  /* For EOM-CCSD Zeta calcs to use ROHF refs for now */
  if(params.wfn == "EOM_CCSD" && params.ref==0 && params.use_zeta) params.ref = 1;

  params.tolerance = 1e-14;
  tol = options.get_int("TOLERANCE");
  params.tolerance = 1.0*pow(10.0,(double) -tol);

  params.memory = Process::environment.get_memory();
  //fndcor(&(params.memory),infile,outfile);

//  params.cachelev = 2;
  params.cachelev = options.get_int("CACHELEV");

  //params.aobasis = 0;
  params.aobasis = options.get_bool("AO_BASIS");

  params.ael = 0;
  params.ael = options.get_bool("AEL");

  params.gauge = options.get_str("GAUGE");
  if( params.gauge != "LENGTH" && params.gauge != "VELOCITY") {
    printf("Invalid choice of gauge: %s\n", params.gauge.c_str());
    throw PsiException("ccdensity: error", __FILE__, __LINE__);
  }

  /*** determine DERTYPE from input */
  params.dertype = 0;
  junk = options.get_str("DERTYPE");
  if(junk == "NONE") params.dertype = 0;
  else if(junk == "FIRST") params.dertype = 1;
  else if(junk == "RESPONSE") params.dertype = 3;
  else {
    printf("Invalid value of input keyword DERTYPE: %s\n", junk.c_str() );
    throw PsiException("ccdensity: error", __FILE__, __LINE__);
  }

  if ( (params.dertype == 1) || params.wfn == "CCSD_MVD")
    params.relax_opdm = 1;  /* default for gradients, or MVD correction */
  else
    params.relax_opdm = 0;  /* otherwise, default is relax_opdm off */

  if(params.transition)
    params.relax_opdm = 0;

  params.relax_opdm = options.get_bool("RELAX_OPDM");
  if ( (params.onepdm) && (params.relax_opdm) ) { /* can't do relaxation without twopdm */
    fprintf(outfile,"\tTurning orbital relaxation off since only onepdm is requested.\n");
    params.relax_opdm = 0;
  }

  if ( params.wfn == "EOM_CCSD" && (params.dertype == 0) )
    params.connect_xi = 0;
  else
    params.connect_xi = 1;
  params.connect_xi = options.get_bool("CONNECT_XI");


  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function    = %6s\n", params.wfn.c_str() );
  fprintf(outfile, "\tReference wfn    = %5s\n", (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tTolerance        = %3.1e\n", params.tolerance);
  fprintf(outfile, "\tCache Level      = %1d\n", params.cachelev);
  fprintf(outfile, "\tAO Basis         = %s\n",
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tOPDM Only        = %s\n",
          params.onepdm ? "Yes" : "No");
  fprintf(outfile, "\tRelax OPDM       = %s\n",
          params.relax_opdm ? "Yes" : "No");
  fprintf(outfile, "\tCompute Xi       = %s\n",
          (params.calc_xi) ? "Yes" : "No");
  fprintf(outfile, "\tUse Zeta         = %s\n",
          (params.use_zeta) ? "Yes" : "No");
  fprintf(outfile, "\tXi connected     = %s\n",
          (params.connect_xi) ? "Yes" : "No");
  fflush(outfile);
}


}} // namespace psi::ccdensity
