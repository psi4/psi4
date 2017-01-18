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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/* WefabL2(): Computes the contribution of the Wefab HBAR matrix
** elements to the Lambda double de-excitation amplitude equations.
** These contributions are given in spin-orbitals as:
**
** L_ij^ab = 1/2 L_ij^ef Wefab
**
** where W_efab is defined in spin orbitals as:
**
** Wefab = <ef||ab> - P(ef) t_m^f <em||ab> + 1/2 tau_mn^ef <mn||ab>
**
** and tau_mn^ef = t_mn^ef + t_m^e t_n^f - t_m^f t_n^e.
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).]
**
** NB: Wefab is not symmetric, so one must be careful when defining
** intermediate quantities for efficient contractions.  I use the
** following contraction steps for each spin case:
**
** Wefab term II spin cases:
**
**   L_IJ^AB <-- 1/2 ( -t_M^F <EM||AB> L_IJ^EF + t_M^E <FM||AB> L_IJ^EF )
**
**     Z(IJ,EM) = -t_M^F L_IJ^EF
**
**   L_IJ^AB <-- Z(IJ,EM) <EM||AB>
*******
**   L_ij^ab <-- 1/2 ( -t_m^f <em||ab> L_ij^ef + t_m^e <fm||ab> L_ij^ef )
**
**     Z(ij,em) = -t_m^f L_ij^ef
**
**   L_ij^ab <-- Z(ij,em) <em||ab>
*******
**   L_Ij^Ab <-- -t_m^f <Em|Ab> L_Ij^Ef - t_M^E <Mf|Ab>  L_Ij^Ef
**
**     Z(Ij,Em) = -t_m^f L_Ij^Ef
**     Z(Ij,Mf) = -t_M^E L_Ij^Ef
**
**   L_Ij^Ab <-- Z(Ij,Em) <Em|Ab> + Z(Ij,Mf) <Mf|Ab>
**
** Wefab term III:
**
**   L_IJ^AB <-- 1/4 tau_MN^EF <MN||AB> L_IJ^EF
**
**     Z(IJ,MN) = 1/2 tau_MN^EF L_IJ^EF
**
**   L_IJ^AB <-- 1/2 Z(IJ,MN) <MN||AB>
*******
**   L_ij^ab <-- 1/4 tau_mn^ef <mn||ab> L_ij^ef
**
**     Z(ij,mn) = 1/2 tau_mn^ef L_ij^ef
**
**   L_ij^ab <-- 1/2 Z(ij,mn) <mn||ab>
*******
**   L_Ij^Ab <-- tau_Mn^Ef <Mn|Ab> L_Ij^Ef
**
**     Z(Ij,Mn) = tau_Mn^Ef L_Ij^Ef
**
**   L_Ij^Ab <-- Z(Ij,Mn) <Mn|Ab>
*******
**
** TDC, July 2002
*/

void WefabL2(int L_irr)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 Tau, T2, Z, Z1, Z2, L, L2, B, D, F, Ltmp;
  dpdfile2 tIA, tia;
  dpdbuf4 tau_a, tau_s, tau;
  dpdbuf4 B_a, B_s;
  dpdbuf4 S, A;
  int h;
  double **B_diag, **tau_diag;
  int ij, Gc, C, c, cc;
  int nbuckets, rows_per_bucket, rows_left, m, row_start, ab, cd, dc, d;
  int nrows, ncols, nlinks;
  psio_address next;

  /* RHS += Wefab*Lijef  */
  if(params.ref == 0) { /** RHF **/

    if(params.abcd == "OLD") {
      global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 5, 0, 5, 0, 0, "ZAbIj");
      global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
      global_dpd_->contract444(&B, &LIjAb, &Z, 0, 0, 1, 0);
      global_dpd_->buf4_close(&B);
      global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
      global_dpd_->buf4_close(&Z);
      global_dpd_->buf4_close(&LIjAb);
    }
    else if(params.abcd == "NEW") {
      timer_on("ABCD:new");

      /* L_a(-)(ij,ab) (i>j, a>b) = L(ij,ab) - L(ij,ba) */
      global_dpd_->buf4_init(&tau_a, PSIF_CC_LAMBDA, L_irr, 4, 9, 0, 5, 1, "LIjAb");
      global_dpd_->buf4_copy(&tau_a, PSIF_CC_LAMBDA, "L(-)(ij,ab)");
      global_dpd_->buf4_close(&tau_a);

      /* L_s(+)(ij,ab) (i>=j, a>=b) = L(ij,ab) + L(ij,ba) */
      global_dpd_->buf4_init(&tau_a, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      global_dpd_->buf4_copy(&tau_a, PSIF_CC_TMP0, "L(+)(ij,ab)");
      global_dpd_->buf4_sort_axpy(&tau_a, PSIF_CC_TMP0, pqsr, 0, 5, "L(+)(ij,ab)", 1);
      global_dpd_->buf4_close(&tau_a);
      global_dpd_->buf4_init(&tau_a, PSIF_CC_TMP0, L_irr, 3, 8, 0, 5, 0, "L(+)(ij,ab)");
      global_dpd_->buf4_copy(&tau_a, PSIF_CC_LAMBDA, "L(+)(ij,ab)");
      global_dpd_->buf4_close(&tau_a);

      timer_on("ABCD:S");
      global_dpd_->buf4_init(&tau_s, PSIF_CC_LAMBDA, L_irr, 3, 8, 3, 8, 0, "L(+)(ij,ab)");
      global_dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      global_dpd_->buf4_init(&S, PSIF_CC_TMP0, L_irr, 8, 3, 8, 3, 0, "S(ab,ij)");
      global_dpd_->contract444(&B_s, &tau_s, &S, 0, 0, 0.5, 0);
      global_dpd_->buf4_close(&S);
      global_dpd_->buf4_close(&B_s);
      global_dpd_->buf4_close(&tau_s);
      timer_off("ABCD:S");

      /* L_diag(ij,c)  = 2 * L(ij,cc)*/

      /* NB: Gcc = 0, and B is totally symmetric, so Gab = 0 */
      /* But Gij = L_irr ^ Gab = L_irr */
      global_dpd_->buf4_init(&tau, PSIF_CC_LAMBDA, L_irr, 3, 8, 3, 8, 0, "L(+)(ij,ab)");
      global_dpd_->buf4_mat_irrep_init(&tau, L_irr);
      global_dpd_->buf4_mat_irrep_rd(&tau, L_irr);
      tau_diag = global_dpd_->dpd_block_matrix(tau.params->rowtot[L_irr], moinfo.nvirt);
      for(ij=0; ij < tau.params->rowtot[L_irr]; ij++)
	for(Gc=0; Gc < moinfo.nirreps; Gc++)
	  for(C=0; C < moinfo.virtpi[Gc]; C++) {
	    c = C + moinfo.vir_off[Gc];
	    cc = tau.params->colidx[c][c];
	    tau_diag[ij][c] = tau.matrix[L_irr][ij][cc];
	  }
      global_dpd_->buf4_mat_irrep_close(&tau, L_irr);

      global_dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      global_dpd_->buf4_init(&S, PSIF_CC_TMP0, L_irr, 8, 3, 8, 3, 0, "S(ab,ij)");
      global_dpd_->buf4_mat_irrep_init(&S, 0);
      global_dpd_->buf4_mat_irrep_rd(&S, 0);

      rows_per_bucket = dpd_memfree()/(B_s.params->coltot[0] + moinfo.nvirt);
      if(rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
      nbuckets = (int) ceil((double) B_s.params->rowtot[0]/(double) rows_per_bucket);
      rows_left = B_s.params->rowtot[0] % rows_per_bucket;

      B_diag = global_dpd_->dpd_block_matrix(rows_per_bucket, moinfo.nvirt);
      next = PSIO_ZERO;
      ncols = tau.params->rowtot[L_irr];
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
      global_dpd_->buf4_mat_irrep_wrt(&S, 0);
      global_dpd_->buf4_mat_irrep_close(&S, 0);
      global_dpd_->buf4_close(&S);
      global_dpd_->buf4_close(&B_s);
      global_dpd_->free_dpd_block(B_diag, rows_per_bucket, moinfo.nvirt);
      global_dpd_->free_dpd_block(tau_diag, tau.params->rowtot[L_irr], moinfo.nvirt);
      global_dpd_->buf4_close(&tau);

      timer_on("ABCD:A");
      global_dpd_->buf4_init(&tau_a, PSIF_CC_LAMBDA, L_irr, 4, 9, 4, 9, 0, "L(-)(ij,ab)");
      global_dpd_->buf4_init(&B_a, PSIF_CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
      global_dpd_->buf4_init(&A, PSIF_CC_TMP0, L_irr, 9, 4, 9, 4, 0, "A(ab,ij)");
      global_dpd_->contract444(&B_a, &tau_a, &A, 0, 0, 0.5, 0);
      global_dpd_->buf4_close(&A);
      global_dpd_->buf4_close(&B_a);
      global_dpd_->buf4_close(&tau_a);
      timer_off("ABCD:A");

      timer_on("ABCD:axpy");
      global_dpd_->buf4_init(&S, PSIF_CC_TMP0, L_irr, 5, 0, 8, 3, 0, "S(ab,ij)");
      global_dpd_->buf4_sort_axpy(&S, PSIF_CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
      global_dpd_->buf4_close(&S);
      global_dpd_->buf4_init(&A, PSIF_CC_TMP0, L_irr, 5, 0, 9, 4, 0, "A(ab,ij)");
      global_dpd_->buf4_sort_axpy(&A, PSIF_CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
      global_dpd_->buf4_close(&A);
      timer_off("ABCD:axpy");
      timer_off("ABCD:new");
    }

    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 10, 0, 10, 0, 0, "Z(Mf,Ij)");
    global_dpd_->contract244(&tIA, &LIjAb, &Z, 1, 2, 0, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&LIjAb);
    global_dpd_->file2_close(&tIA);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 10, 0, 10, 0, 0, "Z(Mf,Ij)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, L_irr, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    global_dpd_->contract444(&F, &Z, &Z1, 1, 1, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&newLIjAb);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_LAMBDA, srqp, 0, 5, "New LIjAb", 1);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 0, 0, 0, 0, 0, "Z(Ij,Mn)");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&Z, &D, &newLIjAb, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_close(&newLIjAb);

  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, L_irr, 7, 2, 7, 2, 0, "ZABIJ");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    global_dpd_->contract444(&B, &LIJAB, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, rspq, 2, 7, "New LIJAB", 1);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");

    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&Ltmp, PSIF_CC_TMP0, L_irr, 2, 10, 2, 10, 0, "Ltmp (I>J,MF)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->contract244(&tIA, &LIJAB, &Ltmp, 1, 2, 1, 1.0, 0.0);
    global_dpd_->contract444(&Ltmp, &F, &newLIJAB, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Ltmp);
    global_dpd_->buf4_close(&LIJAB);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 2, 2, 2, 2, 0, "Z(IJ,MN)");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract444(&Z, &D, &newLIJAB, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&newLIJAB);

    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, L_irr, 7, 2, 7, 2, 0, "Zabij");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    global_dpd_->contract444(&B, &Lijab, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, rspq, 2, 7, "New Lijab", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");

    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&Ltmp, PSIF_CC_TMP0, L_irr, 2, 10, 2, 10, 0, "Ltmp (i>j,mf)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->contract244(&tia, &Lijab, &Ltmp, 1, 2, 1, 1.0, 0.0);
    global_dpd_->contract444(&Ltmp, &F, &newLijab, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Ltmp);
    global_dpd_->buf4_close(&Lijab);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 2, 2, 2, 2, 0, "Z(ij,mn)");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract444(&Z, &D, &newLijab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&newLijab);


    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, L_irr, 5, 0, 5, 0, 0, "ZAbIj");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract444(&B, &LIjAb, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    global_dpd_->buf4_init(&Ltmp, PSIF_CC_TMP1, L_irr, 0, 11, 0, 11, 0, "Lt (Ij,Em)");
    global_dpd_->contract424(&LIjAb, &tia, &Ltmp, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_sort(&Ltmp, PSIF_CC_TMP2, pqsr, 0, 10, "Lt (Ij,mE)");
    global_dpd_->buf4_close(&Ltmp);
    global_dpd_->buf4_init(&Ltmp, PSIF_CC_TMP3, L_irr, 0, 10, 0, 10, 0, "Lt (Ij,Mf)");
    global_dpd_->contract244(&tIA, &LIjAb, &Ltmp, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Ltmp);

    global_dpd_->buf4_close(&LIjAb);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Ltmp, PSIF_CC_TMP3, L_irr, 0, 10, 0, 10, 0, "Lt (Ij,Mf)");
    global_dpd_->contract444(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&Ltmp);
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, pqsr, 10, 5, "<me|ab> (mE,Ab)");
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "<me|ab> (mE,Ab)");
    global_dpd_->buf4_init(&Ltmp, PSIF_CC_TMP2, L_irr, 0, 10, 0, 10, 0, "Lt (Ij,mE)");
    global_dpd_->contract444(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&Ltmp);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 0, 0, 0, 0, 0, "Z(Ij,Mn)");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&Z, &D, &newLIjAb, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");


    /** Z(AB,IJ) = <AB||CD> L(IJ,CD) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 7, 2, 7, 2, 0, "Z(AB,IJ)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
    global_dpd_->contract444(&B, &L2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&L2);
    /** Z(AB,IJ) --> New L(IJ,AB) **/
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, rspq, 2, 7, "New LIJAB", 1);
    global_dpd_->buf4_close(&Z);

    /** Z(IJ,EM) = -L(IJ,EFf) t(M,F) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 2, 21, 2, 21, 0, "Z(IJ,EM)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract424(&L2, &tIA, &Z, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&L2);
    /** New L(IJ,AB) <-- Z(IJ,EM) <EM||AB> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->contract444(&Z, &F, &L2, 0, 1, 1, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    /** Z(IJ,MN) = 1/2 L(IJ,EF) tau_MN^EF **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 2, 2, 2, 2, 0, "Z(IJ,MN)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&L2, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    /** New L(IJ,AB) <-- 1/2 Z(IJ,MN) <MN||AB> **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->contract444(&Z, &D, &L2, 0, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z);


    /** Z(ab,ij) = <ab||cd> L(ij,cd) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 17, 12, 17, 12, 0, "Z(ab,ij)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
    global_dpd_->contract444(&B, &L2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&L2);
    /** Z(ab,ij) --> New L(ij,ab) **/
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, rspq, 12, 17, "New Lijab", 1);
    global_dpd_->buf4_close(&Z);

    /** Z(ij,em) = -L(ij,ef) t(m,f) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 12, 31, 12, 31, 0, "Z(ij,em)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->contract424(&L2, &tia, &Z, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&L2);
    /** New L(ij,ab) <-- Z(ij,em) <em||ab> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    global_dpd_->contract444(&Z, &F, &L2, 0, 1, 1, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    /** Z(ij,mn) = 1/2 L(ij,ef) tau_mn^ef **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 12, 12, 12, 12, 0, "Z(ij,mn)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract444(&L2, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    /** New L(ij,ab) <-- 1/2 Z(ij,mn) <mn||ab> **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract444(&Z, &D, &L2, 0, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z);


    /** Z(Ab,Ij) = <Ab|Cd> L(Ij,Cd) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->contract444(&B, &L2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&L2);
    /** Z(Ab,Ij) --> New L(Ij,Ab) **/
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, rspq, 22, 28, "New LIjAb", 1);
    global_dpd_->buf4_close(&Z);

    /** Z(Ij,Em) = -L(Ij,Ef) t(m,f) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 22, 26, 22, 26, 0, "Z(Ij,Em)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract424(&L2, &tia, &Z, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&L2);
    /** New L(Ij,Ab) <-- Z(Ij,Em) <Em|Ab> **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    global_dpd_->contract444(&Z, &F, &L2, 0, 1, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z);

    /** Z(Ij,Mf) = -t(M,E) L(Ij,Ef) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 22, 24, 22, 24, 0, "Z(Ij,Mf)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract244(&tIA, &L2, &Z, 1, 2, 1, -1, 0);
    global_dpd_->buf4_close(&L2);
    /** New L(Ij,Ab) <-- Z(Ij,Mf) <Mf|Ab> **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract444(&Z, &F, &L2, 0, 1, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z);

    /** Z(Ij,Mn) = L(Ij,Ef) tau(Mn,Ef) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 22, 22, 22, 22, 0, "Z(Ij,Mn)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->contract444(&L2, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    /** New L(Ij,Ab) <-- Z(Ij,Mn) <Mn|Ab> **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract444(&Z, &D, &L2, 0, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  }

}

}} // namespace psi::cclambda
