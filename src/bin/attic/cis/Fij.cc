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
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void Fij_build(void)
{
  dpdfile2 F;
  dpdbuf4 D, T2;

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&F, PSIF_CC_MISC, 0, 0, 0, "FIJ");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    dpd_->contract442(&T2, &D, &F, 0, 0, -1, 0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&D);
    dpd_->file2_close(&F);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_->file2_init(&F, PSIF_CC_MISC, 0, 0, 0, "FIJ");

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 7, 2, 7, 0, "MP2 tIJAB");
    dpd_->contract442(&T2, &D, &F, 0, 0, -1, 0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    dpd_->contract442(&T2, &D, &F, 0, 0, -1, 1);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&F);

    dpd_->file2_init(&F, PSIF_CC_MISC, 0, 2, 2, "Fij");

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 10, 17, 12, 17, 0, "MP2 tijab");
    dpd_->contract442(&T2, &D, &F, 0, 0, -1, 0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    dpd_->contract442(&T2, &D, &F, 1, 1, -1, 1);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&F);
  }
}

}} // namespace psi::cis
