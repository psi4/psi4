/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/psifiles.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psi4-dec.h"
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

  //junk = options.get_str("REFERENCE");
  //if(junk == "RHF") params.ref = 0;
  //else if(junk == "ROHF") params.ref = 1;
  //else if(junk == "UHF") params.ref = 2;
  //else {
  //  printf("Invalid value of input keyword REFERENCE: %s\n", junk.c_str());
  //  throw PsiException("ccdensity: error", __FILE__, __LINE__);
  //}
  
  psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *) &(params.ref),
                    sizeof(int));

  params.onepdm = options.get_bool("ONEPDM");
  params.onepdm_grid_dump = options.get_bool("ONEPDM_GRID_DUMP");

  params.calc_xi = options.get_bool("XI");
  if(params.calc_xi) {
    params.ground = 0;
    params.restart = 0;
  }

  params.use_zeta = options.get_bool("ZETA");
  if(params.use_zeta) {
    params.ground = 0;
    params.restart = 1;
  }

  /* For EOM-CCSD Zeta calcs to use ROHF refs for now */
  if(params.wfn == "EOM_CCSD" && params.ref==0 && params.use_zeta) params.ref = 1;

  params.tolerance = 1e-14;
  params.tolerance = options.get_double("INTS_TOLERANCE");

  params.memory = Process::environment.get_memory();
  //fndcor(&(params.memory),infile,outfile);

//  params.cachelev = 2;
  params.cachelev = options.get_int("CACHELEVEL");

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

  if(options["OPDM_RELAX"].has_changed())
    params.relax_opdm = options.get_bool("OPDM_RELAX");
  if ( (params.onepdm) && (params.relax_opdm) ) { /* can't do relaxation without twopdm */
    outfile->Printf("\tTurning orbital relaxation off since only onepdm is requested.\n");
    params.relax_opdm = 0;
  }

  if ( params.wfn == "EOM_CCSD" && (params.dertype == 0) )
    params.connect_xi = 0;
  else
    params.connect_xi = 1;
  if(options["XI_CONNECT"].has_changed())
    params.connect_xi = options.get_bool("XI_CONNECT");

  params.write_nos = options.get_bool("WRITE_NOS");

  outfile->Printf( "\n\tInput parameters:\n");
  outfile->Printf( "\t-----------------\n");
  outfile->Printf( "\tWave function    = %6s\n", params.wfn.c_str() );
  outfile->Printf( "\tReference wfn    = %5s\n", (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  outfile->Printf( "\tDertype          = %d\n", params.dertype);
  outfile->Printf( "\tTolerance        = %3.1e\n", params.tolerance);
  outfile->Printf( "\tCache Level      = %1d\n", params.cachelev);
  outfile->Printf( "\tAO Basis         = %s\n", params.aobasis ? "Yes" : "No");
  outfile->Printf( "\tOPDM Only        = %s\n", params.onepdm ? "Yes" : "No");
  outfile->Printf( "\tRelax OPDM       = %s\n", params.relax_opdm ? "Yes" : "No");
  outfile->Printf( "\tCompute Xi       = %s\n", (params.calc_xi) ? "Yes" : "No");
  outfile->Printf( "\tUse Zeta         = %s\n", (params.use_zeta) ? "Yes" : "No");
  outfile->Printf( "\tXi connected     = %s\n", (params.connect_xi) ? "Yes" : "No");
  outfile->Printf( "\tCompute NO       = %s\n", (params.write_nos) ? "Yes" : "No");
}


}} // namespace psi::ccdensity
