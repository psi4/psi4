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
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void cc2_faeT2(void) {

  dpdfile2 fme, fME, Fae, FAE, fAE, tIA, tia;
  dpdbuf4 tIjAb, tIJAB, tijab, t2;
  dpdbuf4 newtIjAb, newtIJAB, newtijab;
  dpdbuf4 Zijab;

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&fAE, PSIF_CC_OEI, 0, 1, 1, "fAB");

    dpd_->buf4_init(&Zijab, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZIjAb");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->contract424(&tIjAb, &fAE, &Zijab, 3, 1, 0, 1, 0);
    dpd_->buf4_close(&tIjAb);
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->buf4_axpy(&Zijab, &newtIjAb, 1);
    dpd_->buf4_close(&newtIjAb);
    dpd_->buf4_sort_axpy(&Zijab, PSIF_CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    dpd_->buf4_close(&Zijab);

    dpd_->file2_close(&fAE);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&fME, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&FAE, PSIF_CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_->contract222(&tIA, &fME, &FAE, 1, 1, -1, 0);
    dpd_->file2_close(&FAE);  
    dpd_->file2_close(&fME);  
    dpd_->file2_close(&tIA);

    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->file2_init(&fme, PSIF_CC_OEI, 0, 0, 1, "fia");
    dpd_->file2_init(&Fae, PSIF_CC2_HET1, 0, 1, 1, "CC2 Fae");
    dpd_->contract222(&tia, &fme, &Fae, 1, 1, -1, 0);
    dpd_->file2_close(&Fae);
    dpd_->file2_close(&fme);
    dpd_->file2_close(&tia);

    /** F -> tijab **/
    dpd_->file2_init(&FAE, PSIF_CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_->file2_init(&Fae, PSIF_CC2_HET1, 0, 1, 1, "CC2 Fae");

    /*** AA ***/
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_->contract424(&tIJAB, &FAE, &t2, 3, 1, 0, 1, 0);
    dpd_->contract244(&FAE, &tIJAB, &t2, 1, 2, 1, 1, 1);
    dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    dpd_->buf4_close(&newtIJAB);
    dpd_->buf4_close(&t2);
    dpd_->buf4_close(&tIJAB);

    /*** BB ***/
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_->contract424(&tijab, &Fae, &t2, 3, 1, 0, 1, 0);
    dpd_->contract244(&Fae, &tijab, &t2, 1, 2, 1, 1, 1);
    dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    dpd_->buf4_axpy(&t2, &newtijab, 1);
    dpd_->buf4_close(&newtijab);
    dpd_->buf4_close(&t2);
    dpd_->buf4_close(&tijab);

    /*** AB ***/
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->contract424(&tIjAb, &Fae, &newtIjAb, 3, 1, 0, 1, 1);
    dpd_->contract244(&FAE, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);
    dpd_->buf4_close(&newtIjAb);
    dpd_->buf4_close(&tIjAb);

    dpd_->file2_close(&FAE);  
    dpd_->file2_close(&Fae);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&fME, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&FAE, PSIF_CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_->contract222(&tIA, &fME, &FAE, 1, 1, -1, 0);
    dpd_->file2_close(&FAE);  
    dpd_->file2_close(&fME);  
    dpd_->file2_close(&tIA);

    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&fme, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_->file2_init(&Fae, PSIF_CC2_HET1, 0, 3, 3, "CC2 Fae");
    dpd_->contract222(&tia, &fme, &Fae, 1, 1, -1, 0);
    dpd_->file2_close(&Fae);
    dpd_->file2_close(&fme);
    dpd_->file2_close(&tia);

    /** F -> tijab **/
    dpd_->file2_init(&FAE, PSIF_CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_->file2_init(&Fae, PSIF_CC2_HET1, 0, 3, 3, "CC2 Fae");

    /*** AA ***/
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_->contract424(&tIJAB, &FAE, &t2, 3, 1, 0, 1, 0);
    dpd_->contract244(&FAE, &tIJAB, &t2, 1, 2, 1, 1, 1);
    dpd_->buf4_close(&tIJAB);
    dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    dpd_->buf4_close(&newtIJAB);
    dpd_->buf4_close(&t2);

    /*** BB ***/
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_->contract424(&tijab, &Fae, &t2, 3, 1, 0, 1, 0);
    dpd_->contract244(&Fae, &tijab, &t2, 1, 2, 1, 1, 1);
    dpd_->buf4_close(&tijab);
    dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    dpd_->buf4_axpy(&t2, &newtijab, 1);
    dpd_->buf4_close(&newtijab);
    dpd_->buf4_close(&t2);

    /*** AB ***/
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_->contract424(&tIjAb, &Fae, &newtIjAb, 3, 1, 0, 1, 1);
    dpd_->contract244(&FAE, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);
    dpd_->buf4_close(&newtIjAb);
    dpd_->buf4_close(&tIjAb);

    dpd_->file2_close(&FAE);  
    dpd_->file2_close(&Fae);
  }
}
}} // namespace psi::ccenergy
