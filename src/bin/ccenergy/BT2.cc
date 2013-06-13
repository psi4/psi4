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
#include <string>
#include <cmath>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void BT2(void)
{
  int h;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 B_anti, B;
  dpdbuf4 tauIJAB, tauijab, tauIjAb;
  dpdbuf4 Z1,Z2;
  dpdbuf4 tau_a, tau_s, tau;
  dpdbuf4 B_a, B_s;
  dpdbuf4 S, A;
  double **B_diag, **tau_diag;
  int ij, Gc, C, c, cc;
  int nbuckets, rows_per_bucket, rows_left, m, row_start, ab, cd, dc, d;
  int nrows, ncols, nlinks;
  psio_address next;

  if(params.ref == 0) { /** RHF **/
    if(params.abcd == "OLD") {
#ifdef TIME_CCENERGY
      timer_on("ABCD:old");
#endif
      dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
      dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
      dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
      dpd_->contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
      dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_close(&B);
      dpd_->buf4_close(&tauIjAb);
#ifdef TIME_CCENERGY
      timer_off("ABCD:old");
#endif
    }
    else if(params.abcd == "NEW") {

#ifdef TIME_CCENERGY
      timer_on("ABCD:new");
#endif
      /* tau(-)(ij,ab) (i>j, a>b) = tau(ij,ab) - tau(ij,ba) */
      dpd_->buf4_init(&tau_a, PSIF_CC_TAMPS, 0, 4, 9, 0, 5, 1, "tauIjAb");
      dpd_->buf4_copy(&tau_a, PSIF_CC_TAMPS, "tau(-)(ij,ab)");
      dpd_->buf4_close(&tau_a);

      /* tau_s(+)(ij,ab) (i>=j, a>=b) = tau(ij,ab) + tau(ij,ba) */
      dpd_->buf4_init(&tau_a, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
      dpd_->buf4_copy(&tau_a, PSIF_CC_TMP0, "tau(+)(ij,ab)");
      dpd_->buf4_sort_axpy(&tau_a, PSIF_CC_TMP0, pqsr, 0, 5, "tau(+)(ij,ab)", 1);
      dpd_->buf4_close(&tau_a);
      dpd_->buf4_init(&tau_a, PSIF_CC_TMP0, 0, 3, 8, 0, 5, 0, "tau(+)(ij,ab)");
      dpd_->buf4_copy(&tau_a, PSIF_CC_TAMPS, "tau(+)(ij,ab)");
      dpd_->buf4_close(&tau_a);

#ifdef TIME_CCENERGY
      timer_on("ABCD:S");
#endif
      dpd_->buf4_init(&tau_s, PSIF_CC_TAMPS, 0, 3, 8, 3, 8, 0, "tau(+)(ij,ab)");
      dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      dpd_->buf4_init(&S, PSIF_CC_TMP0, 0, 8, 3, 8, 3, 0, "S(ab,ij)");
      dpd_->contract444(&B_s, &tau_s, &S, 0, 0, 0.5, 0);
      dpd_->buf4_close(&S);
      dpd_->buf4_close(&B_s);
      dpd_->buf4_close(&tau_s);
#ifdef TIME_CCENERGY
      timer_off("ABCD:S");
#endif

      /* tau_diag(ij,c)  = 2 * tau(ij,cc)*/
      dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 3, 8, 3, 8, 0, "tau(+)(ij,ab)");
      dpd_->buf4_mat_irrep_init(&tau, 0);
      dpd_->buf4_mat_irrep_rd(&tau, 0);
      tau_diag = dpd_->dpd_block_matrix(tau.params->rowtot[0], moinfo.nvirt);
      for(ij=0; ij < tau.params->rowtot[0]; ij++)
	for(Gc=0; Gc < moinfo.nirreps; Gc++)
	  for(C=0; C < moinfo.virtpi[Gc]; C++) {
	    c = C + moinfo.vir_off[Gc];
	    cc = tau.params->colidx[c][c];
	    tau_diag[ij][c] = tau.matrix[0][ij][cc];
	  }
      dpd_->buf4_mat_irrep_close(&tau, 0);

      dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      dpd_->buf4_init(&S, PSIF_CC_TMP0, 0, 8, 3, 8, 3, 0, "S(ab,ij)");
      dpd_->buf4_mat_irrep_init(&S, 0);
      dpd_->buf4_mat_irrep_rd(&S, 0);

      rows_per_bucket = dpd_->dpd_memfree()/(B_s.params->coltot[0] + moinfo.nvirt);
      if(rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
      nbuckets = (int) ceil((double) B_s.params->rowtot[0]/(double) rows_per_bucket);
      rows_left = B_s.params->rowtot[0] % rows_per_bucket;

      B_diag = dpd_->dpd_block_matrix(rows_per_bucket, moinfo.nvirt);
      next = PSIO_ZERO;
      ncols = tau.params->rowtot[0];
      nlinks = moinfo.nvirt;
      for(m=0; m < (rows_left ? nbuckets-1:nbuckets); m++) {
	row_start = m * rows_per_bucket;
	nrows = rows_per_bucket;
	if(nrows && ncols && nlinks) {
	  psio_read(PSIF_CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	  C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		  tau_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
	}

      }
      if(rows_left) {
	row_start = m * rows_per_bucket;
	nrows = rows_left;
	if(nrows && ncols && nlinks) {
	  psio_read(PSIF_CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	  C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		  tau_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
	}
      }
      dpd_->buf4_mat_irrep_wrt(&S, 0);
      dpd_->buf4_mat_irrep_close(&S, 0);
      dpd_->buf4_close(&S);
      dpd_->buf4_close(&B_s);
      dpd_->free_dpd_block(B_diag, rows_per_bucket, moinfo.nvirt);
      dpd_->free_dpd_block(tau_diag, tau.params->rowtot[0], moinfo.nvirt);
      dpd_->buf4_close(&tau);

#ifdef TIME_CCENERGY
      timer_on("ABCD:A");
#endif
      dpd_->buf4_init(&tau_a, PSIF_CC_TAMPS, 0, 4, 9, 4, 9, 0, "tau(-)(ij,ab)");
      dpd_->buf4_init(&B_a, PSIF_CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
      dpd_->buf4_init(&A, PSIF_CC_TMP0, 0, 9, 4, 9, 4, 0, "A(ab,ij)");
      dpd_->contract444(&B_a, &tau_a, &A, 0, 0, 0.5, 0);
      dpd_->buf4_close(&A);
      dpd_->buf4_close(&B_a);
      dpd_->buf4_close(&tau_a);
#ifdef TIME_CCENERGY
      timer_off("ABCD:A");
#endif

#ifdef TIME_CCENERGY
      timer_on("ABCD:axpy");
#endif
      dpd_->buf4_init(&S, PSIF_CC_TMP0, 0, 5, 0, 8, 3, 0, "S(ab,ij)");
      dpd_->buf4_sort_axpy(&S, PSIF_CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
      dpd_->buf4_close(&S);
      dpd_->buf4_init(&A, PSIF_CC_TMP0, 0, 5, 0, 9, 4, 0, "A(ab,ij)");
      dpd_->buf4_sort_axpy(&A, PSIF_CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
      dpd_->buf4_close(&A);
#ifdef TIME_CCENERGY
      timer_off("ABCD:axpy");
      timer_off("ABCD:new");
#endif
    }
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    dpd_->buf4_init(&B_anti, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");

    /* AA and BB terms */
    dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 7, 2, 7, 2, 0, "Z(ab,ij)");

    dpd_->contract444(&B_anti, &tauIJAB, &Z1, 0, 0, 1, 0);
    dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 2, 7, "Z(ij,ab)");
    dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(ij,ab)");
    dpd_->buf4_axpy(&Z2, &newtIJAB, 1);
    dpd_->buf4_close(&Z2);

    dpd_->contract444(&B_anti, &tauijab, &Z1, 0, 0, 1, 0);
    dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 2, 7, "Z(ij,ab)");
    dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(ij,ab)");
    dpd_->buf4_axpy(&Z2, &newtijab, 1);
    dpd_->buf4_close(&Z2);

    dpd_->buf4_close(&Z1);

    /* AB term */
    dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    dpd_->contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
    dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 0, 5, "Z(Ij,Ab)");
    dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_->buf4_axpy(&Z2, &newtIjAb, 1);
    dpd_->buf4_close(&Z2);
    dpd_->buf4_close(&Z1);

    dpd_->buf4_close(&B_anti);  
    dpd_->buf4_close(&B);

    dpd_->buf4_close(&tauIJAB);
    dpd_->buf4_close(&tauijab);
    dpd_->buf4_close(&tauIjAb);

    dpd_->buf4_close(&newtIJAB);
    dpd_->buf4_close(&newtijab);
    dpd_->buf4_close(&newtIjAb);

  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");

    dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
    dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 7, 2, 7, 2, 0, "Z(AB,IJ)");
    dpd_->contract444(&B, &tauIJAB, &Z1, 0, 0, 1, 0);
    dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 2, 7, "Z(IJ,AB)");
    dpd_->buf4_close(&Z1);
    dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(IJ,AB)");
    dpd_->buf4_axpy(&Z2, &newtIJAB, 1);
    dpd_->buf4_close(&Z2);
    dpd_->buf4_close(&B);

    dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
    dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 17, 12, 17, 12, 0, "Z(ab,ij)");
    dpd_->contract444(&B, &tauijab, &Z1, 0, 0, 1, 0);
    dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 12, 17, "Z(ij,ab)");
    dpd_->buf4_close(&Z1);
    dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 12, 17, 12, 17, 0, "Z(ij,ab)");
    dpd_->buf4_axpy(&Z2, &newtijab, 1);
    dpd_->buf4_close(&Z2);
    dpd_->buf4_close(&B);

    dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    dpd_->contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
    dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 22, 28, "Z(Ij,Ab)");
    dpd_->buf4_close(&Z1);
    dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
    dpd_->buf4_axpy(&Z2, &newtIjAb, 1);
    dpd_->buf4_close(&Z2);
    dpd_->buf4_close(&B);

    dpd_->buf4_close(&tauIJAB);
    dpd_->buf4_close(&tauijab);
    dpd_->buf4_close(&tauIjAb);

    dpd_->buf4_close(&newtIJAB);
    dpd_->buf4_close(&newtijab);
    dpd_->buf4_close(&newtIjAb);

  }

}

}} // namespace psi::ccenergy
