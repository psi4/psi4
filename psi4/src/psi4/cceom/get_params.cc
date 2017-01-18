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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psi4-dec.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void get_params(Options &options)
{
  params.memory = Process::environment.get_memory();

  params.wfn = options.get_str("WFN");
  if(params.wfn == "EOM_CC2") {
    psio_read_entry(PSIF_CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc), sizeof(double));
    outfile->Printf("\tCC2 energy          (file100) = %20.15f\n",moinfo.ecc);
  }
  else if(params.wfn == "EOM_CCSD") {
    psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc), sizeof(double));
    outfile->Printf("\tCCSD energy         (file100) = %20.15f\n",moinfo.ecc);
  }
  else if(params.wfn == "EOM_CC3") {
    psio_read_entry(PSIF_CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc), sizeof(double));
    outfile->Printf("\tCC3 energy          (file100) = %20.15f\n",moinfo.ecc);
  }

  params.semicanonical = 0;
  std::string read_ref = options.get_str("REFERENCE");
  if(read_ref == "RHF") params.ref = 0;
  else if(read_ref == "ROHF" && params.wfn == "EOM_CC3") {
    params.ref = 2;
    params.semicanonical = 1;
  }
  else if(read_ref == "ROHF") params.ref = 1;
  else if(read_ref == "UHF") params.ref = 2;

  std::string read_eom_ref = options.get_str("EOM_REFERENCE");
  if (params.ref == 0) { /* for RHF refs, allow CCEOM to do RHF, ROHF, UHF */
    if(read_eom_ref == "RHF") params.eom_ref = 0;
    else if(read_eom_ref == "ROHF") params.eom_ref = 1;
    else if(read_eom_ref == "UHF") params.eom_ref = 2;
    else params.eom_ref = 0;
  }
  else if (params.ref == 1) { /* for ROHF refs, allow CCEOM to do ROHF & UHF modes */
    if(read_eom_ref == "ROHF") params.eom_ref = 1;
    else if(read_eom_ref == "UHF") params.eom_ref = 2;
    else params.eom_ref = 1;
  }
  else params.eom_ref = 2; /* run in UHF mode - ignore EOM_REFERENCE */

  params.full_matrix = options["FULL_MATRIX"].to_integer();
  params.cachelev = options.get_int("CACHELEVEL");

  std::string cachetype = options.get_str("CACHETYPE");
  if(cachetype == "LOW") params.cachetype = 1;
  else if(cachetype == "LRU") params.cachetype = 0;
  if(params.ref == 2) /* No LRU cacheing yet for UHF references */
    params.cachetype = 0;

  params.nthreads = Process::environment.get_n_threads();
  if (options["CC_NUM_THREADS"].has_changed()){
     params.nthreads = options.get_int("CC_NUM_THREADS");
  }
  params.abcd = options.get_str("ABCD");
  params.t3_Ws_incore = options["T3_WS_INCORE"].to_integer();
  params.local = options["LOCAL"].to_integer();
  if(params.local) {
    local.cutoff = options.get_double("LOCAL_CUTOFF");
    local.method = options.get_str("LOCAL_METHOD");
    local.weakp = options.get_str("LOCAL_WEAKP");
    local.precon = options.get_str("LOCAL_PRECONDITIONER");
    local.ghost = options.get_int("LOCAL_GHOST");
    local.do_singles = options["LOCAL_DO_SINGLES"].to_integer();
    local.filter_singles = options["LOCAL_FILTER_SINGLES"].to_integer();
  }

  params.newtrips = options["NEW_TRIPLES"].to_integer();


  outfile->Printf( "\n\tInput parameters:\n");
  outfile->Printf( "\t-----------------\n");
  params.overlap = options["OVERLAP_CHECK"].to_integer();

  if(params.semicanonical)
    outfile->Printf( "\tReference wfn   = ROHF changed to UHF for Semicanonical Orbitals\n");
  else
    outfile->Printf( "\tReference wfn   =    %4s\n",
             (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  outfile->Printf( "\tReference EOM wfn=    %4s\n",
             (params.eom_ref == 0) ? "RHF" : ((params.eom_ref == 1) ? "ROHF" : "UHF"));
  outfile->Printf( "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  outfile->Printf( "\tABCD            =     %s\n", params.abcd.c_str());
  outfile->Printf( "\tCache Level     =    %1d\n", params.cachelev);
  outfile->Printf( "\tCache Type      =    %4s\n", params.cachetype ? "LOW" : "LRU");
  if (params.wfn == "EOM_CC3") outfile->Printf( "\tT3 Ws incore  =    %4s\n", params.t3_Ws_incore ? "Yes" : "No");
  outfile->Printf( "\tNum. of threads =     %d\n",params.nthreads);
  outfile->Printf( "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    outfile->Printf( "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    outfile->Printf( "\tLocal Method    =    %s\n", local.method.c_str());
    outfile->Printf( "\tWeak pairs      =    %s\n", local.weakp.c_str());
    outfile->Printf( "\tLocal precon.   =    %s\n", local.precon.c_str());
    outfile->Printf( "\tGhost atom      =    %d\n", local.ghost);
    outfile->Printf( "\tLocal guess     =    %s\n",
	    local.do_singles ? "HBAR_SS" : "UNIT VECTORS" );
    outfile->Printf( "\tFilter singles  =    %s\n", local.filter_singles ? "Yes" : "No");
  }
  outfile->Printf( "\n");
}


}} // namespace psi::cceom
