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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void c_clean(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

/* This function computes the H-bar doubles-doubles block contribution
   from Wabef to a Sigma vector stored at Sigma plus 'i' */

void WabefDD(int i, int C_irr) {
  dpdfile2 tIA, tia, SIA, Sia;
  dpdbuf4 SIJAB, Sijab, SIjAb, B;
  dpdbuf4 CMNEF, Cmnef, CMnEf, X, F, tau, D, WM, WP, Z;
  char CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];
  char SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32], SIA_lbl[32], Sia_lbl[32];
  char lbl_a[32], lbl_s[32];
  dpdbuf4 tau_a, tau_s;
  dpdbuf4 B_a, B_s;
  dpdbuf4 S, A;
  double **B_diag, **tau_diag;
  int ij, Gc, C, c, cc;
  int nbuckets, rows_per_bucket, rows_left, m, row_start, ab, cd, dc, d;
  int nrows, ncols, nlinks;
  psio_address next;

  if (params.eom_ref == 0) { /* RHF */
    /* SIjAb += WAbEf*CIjEf */
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);

    /* SIjAb += <Ab|Ef> CIjEf -- allow out of core algorithm */

#ifdef TIME_CCEOM
    timer_on("WabefDD Z");
#endif

    if(params.abcd == "OLD") {
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
      global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 5, 0, 5, 0, 0, "WabefDD Z(Ab,Ij)");
      global_dpd_->buf4_init(&B, PSIF_CC_BINTS, H_IRR, 5, 5, 5, 5, 0, "B <ab|cd>");
      global_dpd_->contract444(&B, &CMnEf, &Z, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&B);
      global_dpd_->buf4_close(&CMnEf);
      global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, rspq, 0, 5, "WabefDD Z(Ij,Ab)");
      global_dpd_->buf4_close(&Z);

      global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
      global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabefDD Z(Ij,Ab)");
      global_dpd_->buf4_axpy(&Z, &SIjAb, 1);
      global_dpd_->buf4_close(&Z);
      global_dpd_->buf4_close(&SIjAb);
    }
    else if(params.abcd == "NEW") {

      sprintf(lbl_a, "CMnEf(-)(mn,ef) %d", i);
      sprintf(lbl_s, "CMnEf(+)(mn,ef) %d", i);

      /* L_a(-)(ij,ab) (i>j, a>b) = L(ij,ab) - L(ij,ba) */
      global_dpd_->buf4_init(&tau_a, PSIF_EOM_CMnEf, C_irr, 4, 9, 0, 5, 1, CMnEf_lbl);
      global_dpd_->buf4_copy(&tau_a, PSIF_EOM_CMnEf, lbl_a);
      global_dpd_->buf4_close(&tau_a);

      /* L_s(+)(ij,ab) (i>=j, a>=b) = L(ij,ab) + L(ij,ba) */
      global_dpd_->buf4_init(&tau_a, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
      global_dpd_->buf4_copy(&tau_a, PSIF_EOM_TMP, lbl_s);
      global_dpd_->buf4_sort_axpy(&tau_a, PSIF_EOM_TMP, pqsr, 0, 5, lbl_s, 1);
      global_dpd_->buf4_close(&tau_a);
      global_dpd_->buf4_init(&tau_a, PSIF_EOM_TMP, C_irr, 3, 8, 0, 5, 0, lbl_s);
      global_dpd_->buf4_copy(&tau_a, PSIF_EOM_CMnEf, lbl_s);
      global_dpd_->buf4_close(&tau_a);

      timer_on("ABCD:S");
      global_dpd_->buf4_init(&tau_s, PSIF_EOM_CMnEf, C_irr, 3, 8, 3, 8, 0, lbl_s);
      global_dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      global_dpd_->buf4_init(&S, PSIF_EOM_TMP, C_irr, 8, 3, 8, 3, 0, "S(ab,ij)");
      global_dpd_->contract444(&B_s, &tau_s, &S, 0, 0, 0.5, 0);
      global_dpd_->buf4_close(&S);
      global_dpd_->buf4_close(&B_s);
      global_dpd_->buf4_close(&tau_s);
      timer_off("ABCD:S");

      /* L_diag(ij,c)  = 2 * L(ij,cc)*/

      /* NB: Gcc = 0, and B is totally symmetric, so Gab = 0 */
      /* But Gij = L_irr ^ Gab = L_irr */
      global_dpd_->buf4_init(&tau, PSIF_EOM_CMnEf, C_irr, 3, 8, 3, 8, 0, lbl_s);
      global_dpd_->buf4_mat_irrep_init(&tau, C_irr);
      global_dpd_->buf4_mat_irrep_rd(&tau, C_irr);
      tau_diag = global_dpd_->dpd_block_matrix(tau.params->rowtot[C_irr], moinfo.nvirt);
      for(ij=0; ij < tau.params->rowtot[C_irr]; ij++)
	for(Gc=0; Gc < moinfo.nirreps; Gc++)
	  for(C=0; C < moinfo.virtpi[Gc]; C++) {
	    c = C + moinfo.vir_off[Gc];
	    cc = tau.params->colidx[c][c];
	    tau_diag[ij][c] = tau.matrix[C_irr][ij][cc];
	  }
      global_dpd_->buf4_mat_irrep_close(&tau, C_irr);

      global_dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      global_dpd_->buf4_init(&S, PSIF_EOM_TMP, C_irr, 8, 3, 8, 3, 0, "S(ab,ij)");
      global_dpd_->buf4_mat_irrep_init(&S, 0);
      global_dpd_->buf4_mat_irrep_rd(&S, 0);

      rows_per_bucket = dpd_memfree()/(B_s.params->coltot[0] + moinfo.nvirt);
      if(rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
      nbuckets = (int) ceil((double) B_s.params->rowtot[0]/(double) rows_per_bucket);
      rows_left = B_s.params->rowtot[0] % rows_per_bucket;

      B_diag = global_dpd_->dpd_block_matrix(rows_per_bucket, moinfo.nvirt);
      next = PSIO_ZERO;
      ncols = tau.params->rowtot[C_irr];
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
      global_dpd_->free_dpd_block(tau_diag, tau.params->rowtot[C_irr], moinfo.nvirt);
      global_dpd_->buf4_close(&tau);

      timer_on("ABCD:A");
      global_dpd_->buf4_init(&tau_a, PSIF_EOM_CMnEf, C_irr, 4, 9, 4, 9, 0, lbl_a);
      global_dpd_->buf4_init(&B_a, PSIF_CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
      global_dpd_->buf4_init(&A, PSIF_EOM_TMP, C_irr, 9, 4, 9, 4, 0, "A(ab,ij)");
      global_dpd_->contract444(&B_a, &tau_a, &A, 0, 0, 0.5, 0);
      global_dpd_->buf4_close(&A);
      global_dpd_->buf4_close(&B_a);
      global_dpd_->buf4_close(&tau_a);
      timer_off("ABCD:A");

      timer_on("ABCD:axpy");
      global_dpd_->buf4_init(&S, PSIF_EOM_TMP, C_irr, 5, 0, 8, 3, 0, "S(ab,ij)");
      global_dpd_->buf4_sort_axpy(&S, PSIF_EOM_SIjAb, rspq, 0, 5, SIjAb_lbl, 1);
      global_dpd_->buf4_close(&S);
      global_dpd_->buf4_init(&A, PSIF_EOM_TMP, C_irr, 5, 0, 9, 4, 0, "A(ab,ij)");
      global_dpd_->buf4_sort_axpy(&A, PSIF_EOM_SIjAb, rspq, 0, 5, SIjAb_lbl, 1);
      global_dpd_->buf4_close(&A);
      timer_off("ABCD:axpy");
    }

#ifdef TIME_CCEOM
    timer_off("WabefDD Z");
#endif


    /* construct XIjMb = CIjEf * <mb|ef> */
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 10, 0, 10, 0, 0, "WabefDD X(Mb,Ij)");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract444(&F, &CMnEf, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&CMnEf);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 5, 0, 5, 0, 0, "WabefDD Z(Ab,Ij)");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, H_IRR, 0, 1, "tIA");
    /* outfile->Printf("\n begin contract244 in WabefDD\n"); */
    global_dpd_->contract244(&tIA, &X, &Z, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&tIA);
    /* dpd_buf4_print(&Z,outfile,1); */
    global_dpd_->buf4_close(&X);

    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_SIjAb, rspq, 0, 5, SIjAb_lbl, -1);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_SIjAb, srqp, 0, 5, SIjAb_lbl, -1);

    /* SIjAb += tau_MnAb <Mn||ef> CIjEf */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, 0, 0, 0, "WabefDD XIjMn");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEF*CIJEF */
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <ab|cd>");
    global_dpd_->contract444(&CMNEF, &B, &SIJAB, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 2, 10, 2, 10, 0, "XIJMA");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->contract444(&CMNEF, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_M");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, H_IRR, 0, 1, "tIA");
    global_dpd_->contract244(&tIA, &X, &WM, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 2, 5, "WabefDD_P");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_P");
    global_dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract444(&CMNEF, &D, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&X, &tau, &SIJAB, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&X);

    /* Sijab += Wabef*Cijef */
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <ab|cd>");
    global_dpd_->contract444(&Cmnef, &B, &Sijab, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 2, 10, 2, 10, 0, "Xijma");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->contract444(&Cmnef, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_M");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, H_IRR, 0, 1, "tia");
    global_dpd_->contract244(&tia, &X, &WM, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 2, 5, "WabefDD_P");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_P");
    global_dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract444(&Cmnef, &D, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->contract444(&X, &tau, &Sijab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&X);

    /* SIjAb += WAbEf*CIjEf */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);

    /* make use of a more efficient algorithm */
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, H_IRR, 5, 5, 5, 5, 0, "B <ab|cd>");
    /*  dpd_contract444(&CMnEf, &B, &SIjAb, 0, 0, 1.0, 1.0); */
    global_dpd_->contract444(&B, &CMnEf, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, rspq, 0, 5, "Z(Ij,Ab)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_axpy(&Z, &SIjAb, 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 0, 10, 0, 10, 0, "XIjMa");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnfE");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_sort(&X, PSIF_EOM_TMP, pqsr, 0, 11, "XIjaM");
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 0, 11, 0, 11, 0, "XIjaM");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, H_IRR, 0, 1, "tia");
    global_dpd_->contract424(&X, &tia, &SIjAb, 3, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 0, 10, 0, 10, 0, "XIjMb");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, H_IRR, 0, 1, "tIA");
    global_dpd_->contract244(&tIA, &X, &SIjAb, 0, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, 0, 0, 0, "XIjMn");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&X);

    sprintf(SIA_lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, SIA_lbl);
    sprintf(Sia_lbl, "%s %d", "Sia", i);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 0, 1, Sia_lbl);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    c_clean(&SIA,&Sia,&SIJAB,&Sijab,&SIjAb);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) { /* UHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEF*CIJEF */
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <AB|CD>");
    global_dpd_->contract444(&CMNEF, &B, &SIJAB, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 2, 20, 2, 20, 0, "XIJMA");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 20, 7, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract444(&CMNEF, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WABEFDD_M");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, H_IRR, 0, 1, "tIA");
    global_dpd_->contract244(&tIA, &X, &WM, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 2, 5, "WABEFDD_P");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WABEFDD_P");
    global_dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->contract444(&CMNEF, &D, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&X, &tau, &SIJAB, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&X);

    /* Sijab += Wabef*Cijef */
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 12, 17, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, H_IRR, 17, 17, 15, 15, 1, "B <ab|cd>");
    global_dpd_->contract444(&Cmnef, &B, &Sijab, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 12, 30, 12, 30, 0, "Xijma");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 30, 17, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract444(&Cmnef, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WabefDD_M");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, H_IRR, 2, 3, "tia");
    global_dpd_->contract244(&tia, &X, &WM, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 12, 15, "WabefDD_P");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 12, 15, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WabefDD_P");
    global_dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 12, 12, 12, 12, 0, "Xijmn");
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract444(&Cmnef, &D, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 12, 17, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, H_IRR, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->contract444(&X, &tau, &Sijab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&X);

    /* SIjAb += WAbEf*CIjEf */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);

    /* make use of a more efficient algorithm */
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, H_IRR, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    /*  dpd_contract444(&CMnEf, &B, &SIjAb, 0, 0, 1.0, 1.0); */
    global_dpd_->contract444(&B, &CMnEf, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, rspq, 22, 28, "Z(Ij,Ab)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_axpy(&Z, &SIjAb, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&CMnEf);

    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 22, 27, 22, 27, 0, "XIjmA");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 22, 29, 22, 29, 0, "CMnfE");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_sort(&X, PSIF_EOM_TMP, pqsr, 22, 26, "XIjAm");
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 22, 26, 22, 26, 0, "XIjAm");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, H_IRR, 2, 3, "tia");
    global_dpd_->contract424(&X, &tia, &SIjAb, 3, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&X);

    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 22, 24, 22, 24, 0, "XIjMb");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, H_IRR, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, H_IRR, 0, 1, "tIA");
    global_dpd_->contract244(&tIA, &X, &SIjAb, 0, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&X);

    /* Sijab += tau_mneb <mn||ef> C_ijef */
    global_dpd_->buf4_init(&X, PSIF_EOM_TMP, C_irr, 22, 22, 22, 22, 0, "XIjMn");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&X);
  }

#ifdef EOM_DEBUG
  check_sum("WabefDD",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
