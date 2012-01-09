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
     ((params.wfn == "MP2") or (params.wfn == "CCSD_T") or
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
  params.cachelev = options.get_int("CACHELEV");

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
    fprintf(outfile, "\n\tInput parameters:\n");
    fprintf(outfile, "\t-----------------\n");
    fprintf(outfile, "\tWave function   =\t%s\n", params.wfn.c_str());
    fprintf(outfile, "\tBacktransform   =\t%s\n", params.backtr ? "Yes" : "No");
    fprintf(outfile, "\tPrint Level     =\t%d\n", params.print_lvl);
    fprintf(outfile, "\tPrint TEIs      =\t%s\n", params.print_tei ? "Yes" : "No");
    if(params.semicanonical) {
      fprintf(outfile, "\tReference wfn   =\tROHF (using UHF for semicanonical orbitals)\n");
    }
    else {
      fprintf(outfile, "\tReference wfn   =\t%s\n",
          (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
    }
    if(params.dertype == 0) fprintf(outfile, "\tDerivative      =\tNone\n");
    else if(params.dertype == 1) fprintf(outfile, "\tDerivative      =\tFirst\n");
    else if(params.dertype == 2) fprintf(outfile, "\tDerivative      =\tSecond\n");
    else if(params.dertype == 3) fprintf(outfile, "\tDerivative      =\tResponse\n");
    fprintf(outfile, "\tDelete TEI File =\t%s\n", params.delete_tei ? "Yes" : "No");
    fprintf(outfile, "\tMemory (Mbytes) =\t%.1f\n", params.memory/1e6);
    fprintf(outfile, "\tCache Level     =\t%d\n", params.cachelev);
    fprintf(outfile, "\tCache Type      =\t%s\n", "LRU");
    fflush(outfile);
  }
}

  } // namespace transqt2
} // namespace psi
