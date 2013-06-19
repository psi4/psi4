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
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void Fkc_build(int irrep, int root, enum Spin spin)
{
  char lbl[32];
  dpdfile2 B, F, B_A, B_B;
  dpdbuf4 D;

  if(params.ref == 0) { /** RHF **/

    if(spin == singlet)
      sprintf(lbl, "BIA(%d)[%d] singlet", root, irrep);
    else 
      sprintf(lbl, "BIA(%d)[%d] triplet", root, irrep);

    dpd_->file2_init(&B, PSIF_CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "FKC(%d)[%d]", root, irrep);
    dpd_->file2_init(&F, PSIF_CC_MISC, irrep, 0, 1, lbl);

    if(spin == singlet) {
      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      dpd_->dot13(&B, &D, &F, 0, 0, 1, 0);
    }
    else {
      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_->dot14(&B, &D, &F, 0, 0, -1, 0);
    }

    dpd_->buf4_close(&D);

    dpd_->file2_close(&F);
    dpd_->file2_close(&B);
  }
  else if(params.ref == 2) { /** UHF **/

    sprintf(lbl, "BIA(%d)[%d]", root, irrep);
    dpd_->file2_init(&B_A, PSIF_CC_OEI, irrep, 0, 1, lbl);
    sprintf(lbl, "Bia(%d)[%d]", root, irrep);
    dpd_->file2_init(&B_B, PSIF_CC_OEI, irrep, 2, 3, lbl);

    sprintf(lbl, "FKC(%d)[%d]", root, irrep);
    dpd_->file2_init(&F, PSIF_CC_MISC, irrep, 0, 1, lbl);
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_->dot13(&B_A, &D, &F, 0, 0, 1, 0);
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->dot24(&B_B, &D, &F, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&F);

    sprintf(lbl, "Fkc(%d)[%d]", root, irrep);
    dpd_->file2_init(&F, PSIF_CC_MISC, irrep, 2, 3, lbl);
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_->dot13(&B_B, &D, &F, 0, 0, 1, 0);
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->dot13(&B_A, &D, &F, 0, 0, 1, 1);
    dpd_->buf4_close(&D);
    dpd_->file2_close(&F);

    dpd_->file2_close(&B_A);
    dpd_->file2_close(&B_B);
  }
}

}} // namespace psi::cis
