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
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void local_filter_T2(dpdbuf4 *T2);

void init_amps(void)
{
  dpdfile2 tIA, tia, fIA, fia, dIA, dia;
  dpdbuf4 tIJAB, tijab, tIjAb, D, dIJAB, dijab, dIjAb;

  if(params.ref == 0) { /** RHF **/

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    if(!params.restart || !psio_tocscan(PSIF_CC_OEI, "tIA"))
      dpd_->file2_scm(&tIA, 0);
    else fprintf(outfile, "\tUsing old T1 amplitudes.\n");
    dpd_->file2_close(&tIA);

    if(!params.restart || !psio_tocscan(PSIF_CC_TAMPS, "tIjAb")) {
      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb");
      dpd_->buf4_close(&D);

      dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      if(params.local) {
	local_filter_T2(&tIjAb);
      }
      else {
	dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
	dpd_->buf4_dirprd(&dIjAb, &tIjAb);
	dpd_->buf4_close(&dIjAb);
      }
      dpd_->buf4_close(&tIjAb);
    }
    else fprintf(outfile, "\tUsing old T2 amplitudes.\n\n");
  }
  else if(params.ref == 1) { /*** ROHF ***/
    if(!params.restart || !psio_tocscan(PSIF_CC_OEI, "tIA") ||
       !psio_tocscan(PSIF_CC_OEI, "tia")) {

      dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      dpd_->file2_copy(&fIA, PSIF_CC_OEI, "tIA");
      dpd_->file2_close(&fIA);
      dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      /*  dpd_oe_scm(&tIA, 0);  */
      dpd_->file2_close(&tIA);
  
      dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
      dpd_->file2_copy(&fia, PSIF_CC_OEI, "tia");
      dpd_->file2_close(&fia);
      dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
      /*  dpd_oe_scm(&tia, 0); */
      dpd_->file2_close(&tia);

      dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
      dpd_->file2_dirprd(&dIA, &tIA);
      dpd_->file2_close(&tIA);
      dpd_->file2_close(&dIA);

      dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
      dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
      dpd_->file2_dirprd(&dia, &tia);
      dpd_->file2_close(&tia);
      dpd_->file2_close(&dia);
    }
    else fprintf(outfile, "\tUsing old T1 amplitudes.\n");

    if(!params.restart || !psio_tocscan(PSIF_CC_TAMPS, "tIjAb") || 
       !psio_tocscan(PSIF_CC_TAMPS, "tIJAB") || !psio_tocscan(PSIF_CC_TAMPS, "tijab")) {

      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
      dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIJAB");
      dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tijab");
      dpd_->buf4_close(&D);

      dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
      dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      dpd_->buf4_close(&tIJAB);
      dpd_->buf4_close(&dIJAB);

      dpd_->buf4_init(&dijab, PSIF_CC_DENOM, 0, 1, 6, 1, 6, 0, "dijab");
      dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
      dpd_->buf4_dirprd(&dijab, &tijab);
      dpd_->buf4_close(&tijab);
      dpd_->buf4_close(&dijab);

      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb");
      dpd_->buf4_close(&D);
  
      dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_->buf4_dirprd(&dIjAb, &tIjAb);
      dpd_->buf4_close(&tIjAb);
      dpd_->buf4_close(&dIjAb);
    }
    else 
      fprintf(outfile, "\tUsing old T2 amplitudes.\n");
  }
  else if(params.ref == 2) { /*** UHF ***/

    if(!params.restart || !psio_tocscan(PSIF_CC_OEI, "tIA") ||
       !psio_tocscan(PSIF_CC_OEI, "tia")) {

      dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      dpd_->file2_copy(&fIA, PSIF_CC_OEI, "tIA");
      dpd_->file2_close(&fIA);
      dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      /*    dpd_file2_scm(&tIA, 0); */
      dpd_->file2_close(&tIA);

      dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
      dpd_->file2_copy(&fia, PSIF_CC_OEI, "tia");
      dpd_->file2_close(&fia);
      dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
      /*    dpd_file2_scm(&tia, 0); */
      dpd_->file2_close(&tia);

      dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
      dpd_->file2_dirprd(&dIA, &tIA);
      dpd_->file2_close(&tIA);
      dpd_->file2_close(&dIA);

      dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 2, 3, "dia");
      dpd_->file2_dirprd(&dia, &tia);
      dpd_->file2_close(&tia);
      dpd_->file2_close(&dia);
    }
    else fprintf(outfile, "\tUsing old T1 amplitudes.\n");

    if(!params.restart || !psio_tocscan(PSIF_CC_TAMPS, "tIjAb") || 
       !psio_tocscan(PSIF_CC_TAMPS, "tIJAB") || !psio_tocscan(PSIF_CC_TAMPS, "tijab")) {

      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
      dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIJAB");
      dpd_->buf4_close(&D);
      dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
      dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      dpd_->buf4_close(&tIJAB);
      dpd_->buf4_close(&dIJAB);

      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
      dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tijab");
      dpd_->buf4_close(&D);
      dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
      dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
      dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      dpd_->buf4_close(&tIJAB);
      dpd_->buf4_close(&dIJAB);

      dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
      dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb");
      dpd_->buf4_close(&D);
      dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
      dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
      dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      dpd_->buf4_close(&tIJAB);
      dpd_->buf4_close(&dIJAB);
    }
    else 
      fprintf(outfile, "\tUsing old T2 amplitudes.\n");
  }
  else {  /*** RHF/ROHF ***/

  }
}
}} // namespace psi::ccenergy
