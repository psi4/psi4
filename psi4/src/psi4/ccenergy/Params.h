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
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_ccenergy_params_h
#define _psi_src_bin_ccenergy_params_h

#include <string>

namespace psi { namespace ccenergy {

/* Input parameters */
struct Params {
  int maxiter;
  double convergence;
  double e_convergence;
  int restart;
  long int memory;
  std::string aobasis;
  int cachelev;
  int cachetype;
  int ref;
  int diis;
  std::string wfn;
  int print;
  int local;
  int num_amps;
  int print_mp2_amps;
  int brueckner;
  double bconv;
  int analyze;
  int print_pair_energies;
  int spinadapt_energies;
  int semicanonical;
  int dertype;
  int t2_coupled;
  std::string prop;            /* user-selected property */
  int just_energy; /* just compute energy from T amplitudes on disk and quit */
	int just_residuals; /* just compute residuals from T amplitudes on disk and quit */
  std::string abcd;
  int t3_Ws_incore;
  int nthreads;
  int scs;
  int scsn;
  int scscc;
  double scsmp2_scale_os;
  double scsmp2_scale_ss;
  double scscc_scale_os;
  double scscc_scale_ss;
  int newtrips;
  int df;
};

}} // namespace psi::ccenergy

#endif // _psi_src_bin_ccenergy_params_h