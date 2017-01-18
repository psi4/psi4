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

void ltdensity_rohf(struct TD_Params S);
void ltdensity_uhf(struct TD_Params S);
void ltdensity_intermediates(struct TD_Params S);
void sort_ltd_rohf(struct TD_Params S);
void sort_ltd_uhf(struct TD_Params S);
void rtdensity(struct TD_Params S);
void sort_rtd_rohf(struct TD_Params S);
void sort_rtd_uhf(struct TD_Params S);

void tdensity(struct TD_Params S) {

  if(params.ref == 0 || params.ref == 1) {
    ltdensity_rohf(S);
    sort_ltd_rohf(S);
    rtdensity(S);
    sort_rtd_rohf(S);
  }
  else if(params.ref == 2) {
    ltdensity_uhf(S);
    sort_ltd_uhf(S);
    rtdensity(S);
    sort_rtd_uhf(S);
  }

  return;
}


}} // namespace psi::ccdensity
