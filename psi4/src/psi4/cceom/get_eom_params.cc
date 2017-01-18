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
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/corrtab.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

namespace {

void map_irreps(int* array, std::shared_ptr<PointGroup> full, std::shared_ptr<PointGroup> sub)
{
    // If the parent symmetry hasn't been set, no displacements have been made
    if(!full) return;

    // If the point group between the full and sub are the same return
    if (full->symbol() == sub->symbol())
        return;

    // Build the correlation table between full, and subgroup
    CorrelationTable corrtab(full, sub);
    int nirreps = corrtab.n();
    int temp[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    for(int h = 0; h < nirreps; ++h){
        int target = corrtab.gamma(h, 0);
        temp[target] += array[h];
    }
    for(int h = 0; h < nirreps; ++h)
        array[h] = temp[h];
}

}

void get_eom_params(SharedWavefunction ref_wfn, Options &options)
{
  // Number of excited states per irrep
  if (options["ROOTS_PER_IRREP"].has_changed()) {
    // map the symmetry of the input ROOTS_PER_IRREP to account for displacements.
    std::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
    if (old_pg) {
        // This is one of a series of displacements;  check the dimension against the parent point group
        size_t full_nirreps = old_pg->char_table().nirrep();
        if(options["ROOTS_PER_IRREP"].size() != full_nirreps)
            throw PSIEXCEPTION("Input ROOTS_PER_IRREP array has the wrong dimensions");
        int *temp_docc = new int[full_nirreps];
        for(int h = 0; h < full_nirreps; ++h)
            temp_docc[h] = options["ROOTS_PER_IRREP"][h].to_integer();
        map_irreps(temp_docc, old_pg, ref_wfn->molecule()->point_group());
        eom_params.states_per_irrep = new int[moinfo.nirreps];
        for (int h=0; h<moinfo.nirreps; h++)
            eom_params.states_per_irrep[h] = temp_docc[h];
        delete[] temp_docc;
    }
    else {
      if(options["ROOTS_PER_IRREP"].size() != moinfo.nirreps)
        throw PsiException("ROOTS_PER_IRREP is wrong size. Should be number of irreps.", __FILE__, __LINE__);
      eom_params.states_per_irrep = new int[moinfo.nirreps];
      for (int h=0; h < moinfo.nirreps; ++h)
        eom_params.states_per_irrep[h] = options["ROOTS_PER_IRREP"][h].to_integer();
    }
  }
  else throw PsiException("Must provide roots_per_irrep vector in input.", __FILE__, __LINE__);

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

  // Use prop_sym and prop_root only to determine what energy to write
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
  eom_params.print_singles = options["SINGLES_PRINT"].to_integer();
  eom_params.vectors_per_root_SS = options.get_int("SS_VECS_PER_ROOT");
  eom_params.vectors_per_root = options.get_int("VECS_PER_ROOT");
  eom_params.vectors_cc3 = eom_params.prop_root + 9;
  eom_params.vectors_cc3 = options.get_int("VECS_CC3");
  if (eom_params.vectors_cc3 > eom_params.vectors_per_root)
    eom_params.vectors_cc3 = eom_params.vectors_per_root;

  eom_params.collapse_with_last = options.get_bool("COLLAPSE_WITH_LAST");
  eom_params.complex_tol = options.get_double("COMPLEX_TOLERANCE");
  eom_params.residual_tol = options.get_double("R_CONVERGENCE");
  eom_params.residual_tol_SS = options.get_double("SS_R_CONVERGENCE");
  eom_params.eval_tol = options.get_double("E_CONVERGENCE");
  eom_params.eval_tol_SS = options.get_double("SS_E_CONVERGENCE");
  eom_params.amps_to_print = options.get_int("NUM_AMPS_PRINT");
  eom_params.schmidt_add_residual_tol = options.get_double("SCHMIDT_ADD_RESIDUAL_TOLERANCE");
  eom_params.skip_diagSS = options["SS_SKIP_DIAG"].to_integer();
  eom_params.restart_eom_cc3 = options["RESTART_EOM_CC3"].to_integer();
  eom_params.max_iter_SS = 500;
  eom_params.guess = options.get_str("EOM_GUESS");

  outfile->Printf( "\n\tCCEOM parameters:\n");
  outfile->Printf( "\t-----------------\n");
  outfile->Printf( "\tStates sought per irrep     =");
  for(int i = 0; i < moinfo.nirreps; ++i)
    outfile->Printf( " %s %d,", moinfo.irr_labs[i], eom_params.states_per_irrep[i]);

  outfile->Printf("\n");
  outfile->Printf( "\tMax. number of iterations   = %5d\n", eom_params.max_iter);
  outfile->Printf( "\tVectors stored per root     = %5d\n", eom_params.vectors_per_root);
  outfile->Printf( "\tPrint HbarSS iterations?    = %5d\n", eom_params.print_singles);
  outfile->Printf( "\tExcitation range for HBarSS = %5d\n", eom_params.excitation_range);
  outfile->Printf( "\tEigenvalue tolerance        = %5.1e\n", eom_params.eval_tol);
  outfile->Printf( "\tEigenvalue toleranceSS      = %5.1e\n", eom_params.eval_tol_SS);
  outfile->Printf( "\tResidual vector tolerance   = %5.1e\n", eom_params.residual_tol);
  outfile->Printf( "\tResidual vector toleranceSS = %5.1e\n", eom_params.residual_tol_SS);
  outfile->Printf( "\tComplex tolerance           = %5.1e\n", eom_params.complex_tol);
  outfile->Printf( "\tRoot for properties         = %5d\n", eom_params.prop_root + 1);
  outfile->Printf( "\tSym of state for properties = %6s\n", moinfo.irr_labs[eom_params.prop_sym]);
  outfile->Printf( "\tGuess vectors taken from    = %s\n", eom_params.guess.c_str());
  outfile->Printf( "\tRestart EOM CC3             = %s\n", eom_params.restart_eom_cc3?"YES":"NO");
  outfile->Printf( "\tCollapse with last vector   = %s\n", eom_params.collapse_with_last ? "YES":"NO");
  if (eom_params.follow_root) outfile->Printf( "\tRoot following for CC3 turned on.\n");
  outfile->Printf( "\n\n");
}

}} // namespace psi::cceom
