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

void cc2_fmiT2(void) {

  dpdfile2 fia, fIA, Fmi, FMI, fMI, tIA, tia;
  dpdbuf4 tIjAb, tIJAB, tijab, t2;
  dpdbuf4 newtIjAb, newtIJAB, newtijab;
  dpdbuf4 Zijab;

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&fMI, PSIF_CC_OEI, 0, 0, 0, "fIJ");

    dpd_->buf4_init(&Zijab, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZIjAb");
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->contract244(&fMI, &tIjAb, &Zijab, 0, 0, 0, -1, 0);
    dpd_->buf4_close(&tIjAb);
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->buf4_axpy(&Zijab, &newtIjAb, 1);
    dpd_->buf4_close(&newtIjAb);
    dpd_->buf4_sort_axpy(&Zijab, PSIF_CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    dpd_->buf4_close(&Zijab);

    dpd_->file2_close(&fMI);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&FMI, PSIF_CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 1, 0);
    dpd_->file2_close(&FMI); 
    dpd_->file2_close(&fIA); 
    dpd_->file2_close(&tIA);

    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    dpd_->file2_init(&Fmi, PSIF_CC2_HET1, 0, 0, 0, "CC2 Fmi");
    dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 1, 0);
    dpd_->file2_close(&Fmi);
    dpd_->file2_close(&fia); 
    dpd_->file2_close(&tia);

    /** F -> tijab **/
    dpd_->file2_init(&FMI, PSIF_CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_->file2_init(&Fmi, PSIF_CC2_HET1, 0, 0, 0, "CC2 Fmi");

    /*** AA ***/
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_->contract424(&tIJAB, &FMI, &t2, 1, 0, 1, -1, 0);
    dpd_->contract244(&FMI, &tIJAB, &t2, 0, 0, 0, -1, 1);
    dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    dpd_->buf4_close(&newtIJAB);
    dpd_->buf4_close(&t2);
    dpd_->buf4_close(&tIJAB);

    /*** BB ***/
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_->contract424(&tijab, &Fmi, &t2, 1, 0, 1, -1, 0);
    dpd_->contract244(&Fmi, &tijab, &t2, 0, 0, 0, -1, 1);
    dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    dpd_->buf4_axpy(&t2, &newtijab, 1);
    dpd_->buf4_close(&newtijab);
    dpd_->buf4_close(&t2);
    dpd_->buf4_close(&tijab);

    /*** AB ***/
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->contract424(&tIjAb, &Fmi, &newtIjAb, 1, 0, 1, -1, 1);
    dpd_->contract244(&FMI, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);
    dpd_->buf4_close(&newtIjAb);
    dpd_->buf4_close(&tIjAb);

    dpd_->file2_close(&FMI); 
    dpd_->file2_close(&Fmi);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_init(&FMI, PSIF_CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 1, 0);
    dpd_->file2_close(&FMI); 
    dpd_->file2_close(&fIA); 
    dpd_->file2_close(&tIA);

    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_->file2_init(&Fmi, PSIF_CC2_HET1, 0, 2, 2, "CC2 Fmi");
    dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 1, 0);
    dpd_->file2_close(&Fmi);
    dpd_->file2_close(&fia);
    dpd_->file2_close(&tia);

    /** tijab <- Fmi **/
    dpd_->file2_init(&FMI, PSIF_CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_->file2_init(&Fmi, PSIF_CC2_HET1, 0, 2, 2, "CC2 Fmi");

    /*** AA ***/
    dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_->contract424(&tIJAB, &FMI, &t2, 1, 0, 1, -1, 0);
    dpd_->contract244(&FMI, &tIJAB, &t2, 0, 0, 0, -1, 1);
    dpd_->buf4_close(&tIJAB);
    dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    dpd_->buf4_close(&newtIJAB);
    dpd_->buf4_close(&t2);

    /*** BB ***/
    dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_->contract424(&tijab, &Fmi, &t2, 1, 0, 1, -1, 0);
    dpd_->contract244(&Fmi, &tijab, &t2, 0, 0, 0, -1, 1);
    dpd_->buf4_close(&tijab);
    dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
    dpd_->buf4_axpy(&t2, &newtijab, 1);
    dpd_->buf4_close(&newtijab);
    dpd_->buf4_close(&t2);

    /*** AB ***/
    dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_->contract424(&tIjAb, &Fmi, &newtIjAb, 1, 0, 1, -1, 1);
    dpd_->contract244(&FMI, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);
    dpd_->buf4_close(&newtIjAb);
    dpd_->buf4_close(&tIjAb);

    dpd_->file2_close(&Fmi);
    dpd_->file2_close(&FMI); 
  }
}
}} // namespace psi::ccenergy
