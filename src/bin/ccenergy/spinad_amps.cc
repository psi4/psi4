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
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

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

void spinad_amps(void)
{
  dpdfile2 T1, F;
  dpdbuf4 T2AB1, T2AB2, T2, W, W1, W2;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_copy(&T1, PSIF_CC_OEI, "tia");
    dpd_file2_close(&T1);

    dpd_buf4_init(&T2AB1, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_copy(&T2AB1, PSIF_CC_TMP0, "tIjAb");
    dpd_buf4_sort(&T2AB1, PSIF_CC_TMP0, pqsr, 0, 5, "tIjBa");
    dpd_buf4_close(&T2AB1);

    dpd_buf4_init(&T2AB1, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&T2AB2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjBa");
    dpd_buf4_axpy(&T2AB2, &T2AB1, -1.0);
    dpd_buf4_close(&T2AB2);
    dpd_buf4_close(&T2AB1);

    dpd_buf4_init(&T2AB1, PSIF_CC_TMP0, 0, 2, 7, 0, 5, 0, "tIjAb");
    dpd_buf4_copy(&T2AB1, PSIF_CC_TAMPS, "tIJAB");
    dpd_buf4_copy(&T2AB1, PSIF_CC_TAMPS, "tijab");
    dpd_buf4_close(&T2AB1);

    if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {

      /*** Wmbej intermediates ***/
      dpd_buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
      dpd_buf4_copy(&W, PSIF_CC_HBAR, "WmBEj");
      dpd_buf4_copy(&W, PSIF_CC_HBAR, "WMBEJ");
      dpd_buf4_close(&W);

      dpd_buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
      dpd_buf4_copy(&W, PSIF_CC_HBAR, "WmBeJ");
      dpd_buf4_close(&W);

      /* WMBEJ = WMbeJ + WMbEj */
      dpd_buf4_init(&W1, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
      dpd_buf4_init(&W2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
      dpd_buf4_axpy(&W2, &W1, 1);
      dpd_buf4_close(&W2);
      dpd_buf4_close(&W1);

      dpd_buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
      dpd_buf4_copy(&W, PSIF_CC_HBAR, "Wmbej");
      dpd_buf4_close(&W);

      /*** Wmnij intermediates ***/

      dpd_buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
      dpd_buf4_copy(&W, PSIF_CC_TMP0, "WMnIj");
      dpd_buf4_sort(&W, PSIF_CC_TMP0, pqsr, 0, 0, "WMnJi");
      dpd_buf4_close(&W);

      dpd_buf4_init(&W1, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "WMnIj");
      dpd_buf4_init(&W2, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "WMnJi");
      dpd_buf4_axpy(&W2, &W1, -1);
      dpd_buf4_close(&W2);
      dpd_buf4_close(&W1);

      dpd_buf4_init(&W, PSIF_CC_TMP0, 0, 2, 2, 0, 0, 0, "WMnIj");
      dpd_buf4_copy(&W, PSIF_CC_HBAR, "WMNIJ");
      dpd_buf4_copy(&W, PSIF_CC_HBAR, "Wmnij");
      dpd_buf4_close(&W);
    }

    /*** FMI and FAE intermediates ***/

    dpd_file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_copy(&F, PSIF_CC_OEI, "Fmit");
    dpd_file2_close(&F);

    dpd_file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_copy(&F, PSIF_CC_OEI, "Faet");
    dpd_file2_close(&F);

    dpd_file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_file2_copy(&F, PSIF_CC_OEI, "Fme");
    dpd_file2_close(&F);

  }
}
}} // namespace psi::ccenergy
