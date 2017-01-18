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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* spinad_amps(): For RHF references, build the T2 AA and BB amplitudes from
** the existing T2 AB amplitudes and copy the existing T1 A amplitudes
** into B.  Also build the remaining W and F spin cases.
**
** T2(IJ,AB) = T2(ij,ab) = T2(Ij,Ab) - T2(Ij,Ba)
**
** T1(I,A) = T1(i,a)
**
** WMbEj = WmBeJ
** WMbeJ = WmBEj
** WMBEJ = Wmbej = WMbeJ + WMbEj
**
** WMNIJ = Wmnij = WMnIj - WMnJi
**
** FMI = Fmi
** FAE = Fae
** FME = Fme
*/

void CCEnergyWavefunction::spinad_amps(void)
{
  dpdfile2 T1, F;
  dpdbuf4 T2AB1, T2AB2, T2, W, W1, W2;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_copy(&T1, PSIF_CC_OEI, "tia");
    global_dpd_->file2_close(&T1);

    global_dpd_->buf4_init(&T2AB1, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_copy(&T2AB1, PSIF_CC_TMP0, "tIjAb");
    global_dpd_->buf4_sort(&T2AB1, PSIF_CC_TMP0, pqsr, 0, 5, "tIjBa");
    global_dpd_->buf4_close(&T2AB1);

    global_dpd_->buf4_init(&T2AB1, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&T2AB2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjBa");
    global_dpd_->buf4_axpy(&T2AB2, &T2AB1, -1.0);
    global_dpd_->buf4_close(&T2AB2);
    global_dpd_->buf4_close(&T2AB1);

    global_dpd_->buf4_init(&T2AB1, PSIF_CC_TMP0, 0, 2, 7, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_copy(&T2AB1, PSIF_CC_TAMPS, "tIJAB");
    global_dpd_->buf4_copy(&T2AB1, PSIF_CC_TAMPS, "tijab");
    global_dpd_->buf4_close(&T2AB1);

    if(params_.wfn == "CC2" || params_.wfn == "EOM_CC2") {

      /*** Wmbej intermediates ***/
      global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
      global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "WmBEj");
      global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "WMBEJ");
      global_dpd_->buf4_close(&W);

      global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
      global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "WmBeJ");
      global_dpd_->buf4_close(&W);

      /* WMBEJ = WMbeJ + WMbEj */
      global_dpd_->buf4_init(&W1, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
      global_dpd_->buf4_init(&W2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
      global_dpd_->buf4_axpy(&W2, &W1, 1);
      global_dpd_->buf4_close(&W2);
      global_dpd_->buf4_close(&W1);

      global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
      global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "Wmbej");
      global_dpd_->buf4_close(&W);

      /*** Wmnij intermediates ***/

      global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
      global_dpd_->buf4_copy(&W, PSIF_CC_TMP0, "WMnIj");
      global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, pqsr, 0, 0, "WMnJi");
      global_dpd_->buf4_close(&W);

      global_dpd_->buf4_init(&W1, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "WMnIj");
      global_dpd_->buf4_init(&W2, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "WMnJi");
      global_dpd_->buf4_axpy(&W2, &W1, -1);
      global_dpd_->buf4_close(&W2);
      global_dpd_->buf4_close(&W1);

      global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 2, 0, 0, 0, "WMnIj");
      global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "WMNIJ");
      global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "Wmnij");
      global_dpd_->buf4_close(&W);
    }

    /*** FMI and FAE intermediates ***/

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_copy(&F, PSIF_CC_OEI, "Fmit");
    global_dpd_->file2_close(&F);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_copy(&F, PSIF_CC_OEI, "Faet");
    global_dpd_->file2_close(&F);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_copy(&F, PSIF_CC_OEI, "Fme");
    global_dpd_->file2_close(&F);

  }
}
}} // namespace psi::ccenergy
