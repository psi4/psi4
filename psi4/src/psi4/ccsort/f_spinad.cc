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
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void f_spinad(void)
{
  dpdbuf4 F, F1;

  if(params.ref == 0) { /*** RHF ***/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_scmcopy(&F, PSIF_CC_FINTS, "F 2<ia|bc> - <ia|cb>", 2);
    global_dpd_->buf4_sort_ooc(&F, PSIF_CC_TMP0, pqsr, 10, 5, "F <ia|cb>");
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
    global_dpd_->buf4_init(&F1, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F <ia|cb>");
    global_dpd_->buf4_axpy(&F1, &F, -1);
    global_dpd_->buf4_close(&F1);
    global_dpd_->buf4_close(&F);

  }
}

}} // namespace psi::ccsort
