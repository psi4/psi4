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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/

#include <string>

namespace psi { namespace cclambda {

/* Input parameters for cclambda */
struct Params {
  int maxiter;
  double convergence;
  int restart;
  long int memory;
  int cachelev;
  int aobasis;
  std::string wfn;
  int ref;
  int local;  /* boolean for using simulated local-CC framework */
  int nstates; /* total number of L vectors to compute */
  int zeta; /* boolean for solving zeta equations - implies excited state*/
  int print;
  int dertype;
  int diis;
  std::string abcd;
  int sekino;  /* Sekino-Bartlett size-extensive models */
	/* the following should be obseleted now or soon */
  int all; /* find Ls for all excited states plus ground state */
  int ground; /* find L for only ground state */
  int num_amps;
};

struct L_Params {
  int irrep; /* same as corresponding R */
  double R0;    /* same as corresponding R */
  double cceom_energy; /* same as corresponding R */
  int root; /* index of root within irrep */
  int ground; /* boolean, is this a ground state L ? */
  char L1A_lbl[32];
  char L1B_lbl[32];
  char L2AA_lbl[32];
  char L2BB_lbl[32];
  char L2AB_lbl[32];
  char L2RHF_lbl[32];
};


}} // namespace psi::cclambda