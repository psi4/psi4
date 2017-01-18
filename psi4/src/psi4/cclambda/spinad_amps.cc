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
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/* spinad_amps(): For RHF references, build the T2 AA and BB amplitudes from
** the existing T2 AB amplitudes and copy the existing T1 A amplitudes
** into B.
**
** T2(IJ,AB) = T2(ij,ab) = T2(Ij,Ab) - T2(Ij,Ba)
**
** T1(I,A) = T1(i,a)
**
*/

void spinad_amps(void)
{
  dpdfile2 T1;
  dpdbuf4 T2;

  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&T1, PSIF_CC_LAMBDA, 0, 0, 1, "LIA");
    global_dpd_->file2_copy(&T1, PSIF_CC_LAMBDA, "Lia");
    global_dpd_->file2_close(&T1);

    global_dpd_->buf4_init(&T2, PSIF_CC_LAMBDA, 0, 2, 7, 0, 5, 1, "LIjAb");
    global_dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "LIJAB");
    global_dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "Lijab");
    global_dpd_->buf4_close(&T2);

  }
}

}} // namespace psi::cclambda
