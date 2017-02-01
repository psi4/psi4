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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void denom1(dpdfile2 *X1, double omega);
void denom2(dpdbuf4 *X2, double omega);
void local_filter_T1(dpdfile2 *T1);
void local_filter_T2(dpdbuf4 *T2);

void init_X(const char *pert, int irrep, double omega)
{
  char lbl[32];
  dpdfile2 mu1, X1, FAE, FMI;
  dpdbuf4 X2, mu2;

  sprintf(lbl, "%sBAR_IA", pert);
  global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  if(!params.restart || !psio_tocscan(PSIF_CC_OEI, lbl)) {
    global_dpd_->file2_copy(&mu1, PSIF_CC_OEI, lbl);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    if(params.local && local.filter_singles) local_filter_T1(&X1);
    else denom1(&X1, omega);
    global_dpd_->file2_close(&X1);
  }
  else outfile->Printf( "\tUsing existing %s amplitudes.\n", lbl);
  global_dpd_->file2_close(&mu1);

  sprintf(lbl, "%sBAR_IjAb", pert);
  global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  if(!params.restart || !psio_tocscan(PSIF_CC_LR, lbl)) {
    global_dpd_->buf4_copy(&mu2, PSIF_CC_LR, lbl);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    if(params.local) local_filter_T2(&X2);
    else denom2(&X2, omega);
    global_dpd_->buf4_close(&X2);
  }
  else outfile->Printf( "\tUsing existing %s amplitudes.\n", lbl);
  global_dpd_->buf4_close(&mu2);
}

}} // namespace psi::ccresponse
