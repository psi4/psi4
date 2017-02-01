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

#ifndef _psi_src_bin_cceom_params_h
#define _psi_src_bin_cceom_params_h

#include <string>

namespace psi { namespace cceom {

/* Input parameters for ccenergy */
struct Params {
  long int memory;
  int cachelev;
  int cachetype;
  int ref;
  int eom_ref;
  int local;
  std::string wfn;
  int semicanonical;
  int full_matrix; /* include reference rows/cols in diagonalization */
  std::string abcd;
  int t3_Ws_incore;
  int nthreads;
  int newtrips;
  int overlap; // check for overlaps between current wfn set and older set stored on disk
};

struct Eom_params {
  int max_iter;
  int vectors_per_root;
  int *states_per_irrep;
  int *cs_per_irrep;
  double *state_energies;
  int number_of_states;
  double eval_tol;
  double eval_tol_SS;
  double residual_tol;
  int prop_root;
  int prop_sym;
  int save_all;
  int print_singles;
  double complex_tol;
  double schmidt_add_residual_tol;
  int max_iter_SS;
  int vectors_per_root_SS;
  int excitation_range;
  double residual_tol_SS;
  std::string guess;
  int rhf_triplets;
  int mult;
  bool follow_root;
  int collapse_with_last;
  int skip_diagSS;
  int vectors_cc3;
  int restart_eom_cc3;
  int amps_to_print;

  /* compute overlap of normalized R with L (must run cclambda first) */
  int dot_with_L;
  double L0;
  int L_irr;
};


}} // namespace psi::cceom

#endif // _psi_src_bin_cceom_params_h