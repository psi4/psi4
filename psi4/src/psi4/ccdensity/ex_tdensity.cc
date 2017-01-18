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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void ex_tdensity_rohf(struct TD_Params S, struct TD_Params U);
void ex_tdensity_uhf(struct TD_Params S, struct TD_Params U);
void ex_sort_td_rohf(char hand, int Tirrep);
void ex_sort_td_uhf(char hand, int Tirrep);

void ex_tdensity(char hand, struct TD_Params S, struct TD_Params U) {
  // FYI: "Sort" needs to know L or R in order to
  //       put density in correct place (ltd or rtd)
  int Tirrep = S.irrep^U.irrep;
  if(params.ref == 0 || params.ref == 1) {
    ex_tdensity_rohf(S,U);
    outfile->Printf( "\t\t***...density has been built...\n");

    ex_sort_td_rohf(hand,Tirrep);
    outfile->Printf( "\t\t***...density has been sorted...\n");

  }
  else if(params.ref == 2) {
    ex_tdensity_uhf(S,U);
    outfile->Printf( "\t\t***...density has been built...\n");

    ex_sort_td_uhf(hand,Tirrep);
    outfile->Printf( "\t\t***...density has been sorted...\n");

  }

  return;
}

}} // namespace psi::ccdensity
