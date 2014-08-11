/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <psifiles.h>
#include <psi4-dec.h>
#define EXTERN
#include "globals.h"

namespace psi {
  namespace transqt2 {

void get_params(Options & options)
{
  int tol;

  params.wfn = options.get_str("WFN");

  /* NB: SCF wfns are allowed because, at present, ccsort is needed for
     RPA-type calculations */

  params.semicanonical = 0;
  std::string reference = options.get_str("REFERENCE");
  params.ref = 0; /* if no reference is given, assume rhf */

  if(reference == "RHF"){
    params.ref = 0;
  }
  else if((reference == "ROHF") and
     ((params.wfn == "MP2") or (params.wfn == "CCSD_T") or (params.wfn == "CCSD_AT") or
     (params.wfn == "CC3") or (params.wfn == "EOM_CC3") or
     (params.wfn == "CC2") or (params.wfn == "EOM_CC2"))) {
      params.ref = 2;
      params.semicanonical = 1;
    }
  else if(reference == "ROHF") params.ref = 1;
  else if(reference == "UHF") params.ref = 2;
  else if(reference == "TWOCON") params.ref = 1; /* Treat twocon as rhf */
  else
    throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);

  // Allow user to force semicanonical
  if(options["SEMICANONICAL"].has_changed()) {
   params.semicanonical = options.get_bool("SEMICANONICAL");
   params.ref = 2;
  }

  params.dertype = 0;
  std::string dertype = options.get_str("DERTYPE");
  if(dertype == "NONE") params.dertype = 0;
  else if(dertype == "FIRST") params.dertype = 1;
  else if(dertype == "SECOND") params.dertype = 2;
  else if(dertype == "RESPONSE") params.dertype = 3; /* linear response */

  params.print_lvl = options.get_int("PRINT");

  params.print_tei = 0;
  params.print_tei = options.get_bool("PRINT_TEI");

  params.tolerance = 1e-14;
  params.tolerance = options.get_double("INTS_TOLERANCE");

  params.memory = Process::environment.get_memory();

  params.cachelev = 2;
  params.cachelev = options.get_int("CACHELEVEL");

  params.delete_tei = 1;
  /* If AO-basis chosen, keep the SO_TEI file */
  std::string aobasis;
  aobasis = options.get_str("AO_BASIS");
  if (aobasis == "DISK" || options.get_bool("DELETE_TEI") == false) params.delete_tei = 0;

  // any MCSCF-type wavefunction needs multiple transforms so don't delete
  // the AO two-electron ints
  if ((params.wfn == "OOCCD") || (params.wfn == "DETCAS") ||
       (params.wfn == "CASSCF") || (params.wfn == "RASSCF"))
    params.delete_tei = 0;

  if(params.print_lvl) {
    outfile->Printf( "\n\tInput parameters:\n");
    outfile->Printf( "\t-----------------\n");
    outfile->Printf( "\tWave function   =\t%s\n", params.wfn.c_str());
    outfile->Printf( "\tBacktransform   =\t%s\n", params.backtr ? "Yes" : "No");
    outfile->Printf( "\tPrint Level     =\t%d\n", params.print_lvl);
    outfile->Printf( "\tPrint TEIs      =\t%s\n", params.print_tei ? "Yes" : "No");
    if(params.semicanonical) {
      outfile->Printf( "\tReference wfn   =\tROHF (using UHF for semicanonical orbitals)\n");
    }
    else {
      outfile->Printf( "\tReference wfn   =\t%s\n",
          (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
    }
    if(params.dertype == 0) outfile->Printf( "\tDerivative      =\tNone\n");
    else if(params.dertype == 1) outfile->Printf( "\tDerivative      =\tFirst\n");
    else if(params.dertype == 2) outfile->Printf( "\tDerivative      =\tSecond\n");
    else if(params.dertype == 3) outfile->Printf( "\tDerivative      =\tResponse\n");
    outfile->Printf( "\tDelete TEI File =\t%s\n", params.delete_tei ? "Yes" : "No");
    outfile->Printf( "\tMemory (Mbytes) =\t%.1f\n", params.memory/1e6);
    outfile->Printf( "\tCache Level     =\t%d\n", params.cachelev);
    outfile->Printf( "\tCache Type      =\t%s\n", "LRU");
    
  }
}

  } // namespace transqt2
} // namespace psi
