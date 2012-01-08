/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void get_eom_params(Options &options)
{
  // Number of excited states per irrep
  chkpt_init(PSIO_OPEN_OLD);
  if (chkpt_rd_override_occ()) eom_params.states_per_irrep = chkpt_rd_statespi();
  else if (options["STATES_PER_IRREP"].has_changed()) {
    if(options["STATES_PER_IRREP"].size() != moinfo.nirreps) 
      throw PsiException("STATES_PER_IRREP is wrong size. Should be number of irreps.", __FILE__, __LINE__);
    eom_params.states_per_irrep = new int[moinfo.nirreps];
    for (int h=0; h < moinfo.nirreps; ++h)
      eom_params.states_per_irrep[h] = options["STATES_PER_IRREP"][h].to_integer();
  }
  else throw PsiException("Must provide states_per_irrep vector in input.", __FILE__, __LINE__);
  chkpt_close();

  // Number of guess vectors per irrep
  eom_params.cs_per_irrep = new int[moinfo.nirreps];
  eom_params.number_of_states = 0;
  for (int state_irrep = 0; state_irrep < moinfo.nirreps; ++state_irrep) {
    eom_params.cs_per_irrep[state_irrep^moinfo.sym] = eom_params.states_per_irrep[state_irrep];
    eom_params.number_of_states += eom_params.states_per_irrep[state_irrep];
  }
  eom_params.state_energies = new double[eom_params.number_of_states];

  eom_params.max_iter = 80 * moinfo.nirreps;
  eom_params.max_iter = options.get_int("MAXITER");

  // Use prop_sym and prop_root only to determine what energy to write to chkpt
  if (options["PROP_SYM"].has_changed()) {
    eom_params.prop_sym = options.get_int("PROP_SYM");
    eom_params.prop_sym = (eom_params.prop_sym - 1)^moinfo.sym;
  }
  else {
    for (int i = 0;i < moinfo.nirreps; ++i)
      if (eom_params.states_per_irrep[i]) eom_params.prop_sym = i^moinfo.sym;
  }
  if (options["PROP_ROOT"].has_changed()) {
    eom_params.prop_root = options.get_int("PROP_ROOT");
    if (eom_params.prop_root > eom_params.states_per_irrep[eom_params.prop_sym^moinfo.sym])
      throw PsiException("Value of prop_root is too large.", __FILE__, __LINE__);
  }
  else eom_params.prop_root = eom_params.states_per_irrep[eom_params.prop_sym^moinfo.sym];
  --eom_params.prop_root;

  // by default for CC3 use follow_root to root-following based on CCSD soln
  if (params.wfn == "EOM_CC3" && eom_params.prop_root != 0) eom_params.follow_root = true;
  // allow user to explicitly turn off root-following which may help in bizarre cases
  // in this case prop_root is used to determine which residuals are used
  if (options["CC3_FOLLOW_ROOT"].has_changed()) eom_params.follow_root = options.get_bool("CC3_FOLLOW_ROOT");

  /* so far, all R's are always kept so this is not used */
  eom_params.save_all = 0;

  eom_params.mult = 1;
  eom_params.rhf_triplets = options.get_bool("RHF_TRIPLETS");
  if (eom_params.rhf_triplets) eom_params.mult = 3;

  eom_params.excitation_range = options.get_int("EXCITATION_RANGE");
  eom_params.print_singles = options["PRINT_SINGLES"].to_integer();
  eom_params.vectors_per_root_SS = options.get_int("VECS_PER_ROOT_SS");
  eom_params.vectors_per_root = options.get_int("VECS_PER_ROOT");
  eom_params.vectors_cc3 = eom_params.prop_root + 9;
  eom_params.vectors_cc3 = options.get_int("VECS_CC3");
  if (eom_params.vectors_cc3 > eom_params.vectors_per_root)
    eom_params.vectors_cc3 = eom_params.vectors_per_root;

  eom_params.collapse_with_last = options.get_bool("COLLAPSE_WITH_LAST");
  eom_params.complex_tol = options.get_double("COMPLEX_TOL");
  eom_params.residual_tol = options.get_double("RESIDUAL_TOL");
  eom_params.residual_tol_SS = options.get_double("RESIDUAL_TOL_SS");
  eom_params.eval_tol = options.get_double("EVAL_TOL");
  eom_params.eval_tol_SS = options.get_double("EVAL_TOL_SS");
  eom_params.amps_to_print = options.get_int("AMPS_TO_PRINT");
  eom_params.schmidt_add_residual_tol = options.get_double("SCHMIDT_ADD_RESIDUAL_TOL");
  eom_params.skip_diagSS = options["SKIP_DIAGSS"].to_integer();
  eom_params.restart_eom_cc3 = options["RESTART_EOM_CC3"].to_integer();
  eom_params.max_iter_SS = 500;
  eom_params.guess = options.get_str("EOM_GUESS");

  fprintf(outfile, "\n\tCCEOM parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tStates sought per irrep     =");
  for(int i = 0; i < moinfo.nirreps; ++i) 
    fprintf(outfile, " %s %d,", moinfo.irr_labs[i], eom_params.states_per_irrep[i]);

  fprintf(outfile,"\n");
  fprintf(outfile, "\tMax. number of iterations   = %5d\n", eom_params.max_iter);
  fprintf(outfile, "\tVectors stored per root     = %5d\n", eom_params.vectors_per_root);
  fprintf(outfile, "\tPrint HbarSS iterations?    = %5d\n", eom_params.print_singles);
  fprintf(outfile, "\tExcitation range for HBarSS = %5d\n", eom_params.excitation_range);
  fprintf(outfile, "\tEigenvalue tolerance        = %5.1e\n", eom_params.eval_tol);
  fprintf(outfile, "\tEigenvalue toleranceSS      = %5.1e\n", eom_params.eval_tol_SS);
  fprintf(outfile, "\tResidual vector tolerance   = %5.1e\n", eom_params.residual_tol);
  fprintf(outfile, "\tResidual vector toleranceSS = %5.1e\n", eom_params.residual_tol_SS);
  fprintf(outfile, "\tComplex tolerance           = %5.1e\n", eom_params.complex_tol);
  fprintf(outfile, "\tRoot for properties         = %5d\n", eom_params.prop_root + 1);
  fprintf(outfile, "\tSym of state for properties = %6s\n", moinfo.irr_labs[eom_params.prop_sym]);
  fprintf(outfile, "\tGuess vectors taken from    = %s\n", eom_params.guess.c_str());
  fprintf(outfile, "\tRestart EOM CC3             = %s\n", eom_params.restart_eom_cc3?"YES":"NO");
  fprintf(outfile, "\tCollapse with last vector   = %s\n", eom_params.collapse_with_last ? "YES":"NO");
  if (eom_params.follow_root) fprintf(outfile, "\tRoot following for CC3 turned on.\n");
  fprintf(outfile, "\n\n");
}

}} // namespace psi::cceom
