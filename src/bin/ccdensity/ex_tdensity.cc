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
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

//void ex_tdensity_rohf(struct TD_Params S, struct TD_Params U);
//void ex_tdensity_uhf(struct TD_Params S, struct TD_Params U);
//void ex_tdensity_intermediates(struct TD_Params S, struct TD_Params U);
void ex_tdensity_rohf(char hand, struct TD_Params S, struct TD_Params U);
void ex_tdensity_uhf(char hand, struct TD_Params S, struct TD_Params U);
void ex_tdensity_intermediates(char hand, struct TD_Params S, struct TD_Params U);
void ex_sort_td_rohf(char hand, struct TD_Params S);
void ex_sort_td_uhf(char hand, struct TD_Params S);

void ex_tdensity(char hand, struct TD_Params S, struct TD_Params U) {
  /*  "Density code" might need L or R for one reason:
   *  1) Which state's irrep to use for transition density
   *     being constructed (the higher state's irrep).
   *  "Sort" might need to know L or R for two reasons:
   *  1) Where to put density (ltd or rtd) -> DEFINITELY THIS
   *  2) Which state's irrep to use for sorting -> MAYBE THIS
   */
  if(params.ref == 0 || params.ref == 1) {
    //ex_tdensity_rohf(S,U);
    //ex_sort_td_rohf(hand,U);
    ex_tdensity_rohf(hand,S,U);
    fprintf(outfile, "    *** A density has been built.\n");
    fflush(outfile);
    if(hand=='l') ex_sort_td_rohf(hand,U);
    if(hand=='r') ex_sort_td_rohf(hand,S);
    fprintf(outfile, "    *** A density has been sorted.\n");
    fflush(outfile);
  }
  else if(params.ref == 2) {
    //ex_tdensity_uhf(S,U);
    //ex_sort_td_uhf(hand, U);
    ex_tdensity_uhf(hand,S,U);
    if(hand=='l') ex_sort_td_uhf(hand,U);
    if(hand=='r') ex_sort_td_uhf(hand,S);
  }

  return;
}

}} // namespace psi::ccdensity
