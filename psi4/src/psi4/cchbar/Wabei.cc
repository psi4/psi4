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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <string>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* Wabei(): Computes all contributions to the Wabei HBAR matrix
** elements, whose spin-orbital definition is:
**
** Wabei = <ab||ei> - Fme t(mi,ab) + t(i,f) <ab||ef>
**            (I)         (II)             (IIIa)
**   - P(ab) t(i,f) t(m,b) <am||ef> + 1/2 t(i,f) t(mn,ab) <mn||ef>
**               (IIIb)                            (IIIc)
**   + 1/2 P(ab) t(i,f) t(m,a) t(n,b) <mn||ef>
**                   (IIId)
**   + 1/2 t(mn,ab) <mn||ei> + 1/2 P(ab) t(m,a) t(n,b) <mn||ei>
**             (IVa)                         (IVb)
**   + P(ab) t(mi,fb) <am||ef> - P(ab) t(m,a) <mb||ei>
**               (V)                    (VI)
**   - P(ab) t(m,a) t(ni,fb) <mn||ef>
**                 (VII)
**
**  [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
*/

void Wabei_RHF(void);
void Wabei_ROHF(void);
void WABEI_UHF(void);
void Wabei_UHF(void);
void WAbEi_UHF(void);
void WaBeI_UHF(void);
void Wabei_UHF_sort_ints();

void Wabei_build(void)
{
  if(params.ref == 0) Wabei_RHF();
  else if(params.ref == 1) Wabei_ROHF();
  else if(params.ref == 2) {
    Wabei_UHF_sort_ints();
    WABEI_UHF();
    Wabei_UHF();
    WAbEi_UHF();
    WaBeI_UHF();
  }
}

void Wabei_UHF_sort_ints(void){

  dpdbuf4 F, B;
  //required in WABEI
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 5, 21, 5, 1 , "F <AI|BC>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 5, 20, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_close(&F);

  // required in Wabei
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 15, 31, 15, 1 , "F <ai|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 15, 30, "F <ai||bc> (ab,ic)");
  global_dpd_->buf4_close(&F);

  // required in WAbEi
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0 , "F <Ia|Bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, qrps, 29,24, "F <Ia|Bc> (aB,Ic)");
  global_dpd_->buf4_close(&F);

  // required in WaBeI
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  global_dpd_->buf4_sort(&B,PSIF_CC_BINTS, qpsr,29,29,"B <aB|cD>");
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0 , "F <iA|bC>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, qrps, 28,27, "F <iA|bC> (Ab,iC)");
  global_dpd_->buf4_close(&F);
}

}} // namespace psi::cchbar
