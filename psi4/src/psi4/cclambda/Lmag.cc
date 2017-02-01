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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void Lmag(int L_irr)
{
  dpdfile2 R1, L1, LIA, Lia, RIA, Ria;
  dpdbuf4 R2, L2, LIJAB, Lijab, LIjAb, RIJAB, Rijab, RIjAb;
  double norm;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    norm = global_dpd_->file2_dot_self(&LIA);
    norm += global_dpd_->file2_dot_self(&Lia);
    norm += global_dpd_->buf4_dot_self(&LIJAB);
    norm += global_dpd_->buf4_dot_self(&Lijab);
    norm += global_dpd_->buf4_dot_self(&LIjAb);
    outfile->Printf("size of L <L|L>     %15.10lf\n",norm);

    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&Lia);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_close(&LIjAb);
  }
}

}} // namespace psi::cclambda
