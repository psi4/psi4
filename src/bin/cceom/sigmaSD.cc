/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void FSD(int i, int C_irr);
void WamefSD(int i, int C_irr);
void WmnieSD(int i, int C_irr);

/* This function computes the H-bar singles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaSD(int i, int C_irr) {

#ifdef TIME_CCEOM
  timer_on("FSD");     FSD(i, C_irr);     timer_off("FSD");
  timer_on("WamefSD"); WamefSD(i, C_irr); timer_off("WamefSD");
  timer_on("WmnieSD"); WmnieSD(i, C_irr); timer_off("WmnieSD");
#else
  FSD(i, C_irr);
  WamefSD(i, C_irr);
  WmnieSD(i, C_irr);
#endif

  return;
}

}} // namespace psi::cceom