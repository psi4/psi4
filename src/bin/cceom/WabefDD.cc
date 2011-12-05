/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
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
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
      dpd_buf4_init(&Z, EOM_TMP, C_irr, 5, 0, 5, 0, 0, "WabefDD Z(Ab,Ij)");
      dpd_buf4_init(&B, CC_BINTS, H_IRR, 5, 5, 5, 5, 0, "B <ab|cd>");
      dpd_contract444(&B, &CMnEf, &Z, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&B);
      dpd_buf4_close(&CMnEf);
      dpd_buf4_sort(&Z, EOM_TMP, rspq, 0, 5, "WabefDD Z(Ij,Ab)");
      dpd_buf4_close(&Z);

      dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
      dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabefDD Z(Ij,Ab)");
      dpd_buf4_axpy(&Z, &SIjAb, 1);
      dpd_buf4_close(&Z);
      dpd_buf4_close(&SIjAb);
    }
    else if(params.abcd == "NEW") {

      sprintf(lbl_a, "CMnEf(-)(mn,ef) %d", i);
      sprintf(lbl_s, "CMnEf(+)(mn,ef) %d", i);

      /* L_a(-)(ij,ab) (i>j, a>b) = L(ij,ab) - L(ij,ba) */
      dpd_buf4_init(&tau_a, EOM_CMnEf, C_irr, 4, 9, 0, 5, 1, CMnEf_lbl);
      dpd_buf4_copy(&tau_a, EOM_CMnEf, lbl_a);
      dpd_buf4_close(&tau_a);

      /* L_s(+)(ij,ab) (i>=j, a>=b) = L(ij,ab) + L(ij,ba) */
      dpd_buf4_init(&tau_a, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
      dpd_buf4_copy(&tau_a, EOM_TMP, lbl_s);
      dpd_buf4_sort_axpy(&tau_a, EOM_TMP, pqsr, 0, 5, lbl_s, 1);
      dpd_buf4_close(&tau_a);
      dpd_buf4_init(&tau_a, EOM_TMP, C_irr, 3, 8, 0, 5, 0, lbl_s);
      dpd_buf4_copy(&tau_a, EOM_CMnEf, lbl_s);
      dpd_buf4_close(&tau_a);

      timer_on("ABCD:S");
      dpd_buf4_init(&tau_s, EOM_CMnEf, C_irr, 3, 8, 3, 8, 0, lbl_s);
      dpd_buf4_init(&B_s, CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      dpd_buf4_init(&S, EOM_TMP, C_irr, 8, 3, 8, 3, 0, "S(ab,ij)");
      dpd_contract444(&B_s, &tau_s, &S, 0, 0, 0.5, 0);
      dpd_buf4_close(&S);
      dpd_buf4_close(&B_s);
      dpd_buf4_close(&tau_s);
      timer_off("ABCD:S");

      /* L_diag(ij,c)  = 2 * L(ij,cc)*/

      /* NB: Gcc = 0, and B is totally symmetric, so Gab = 0 */
      /* But Gij = L_irr ^ Gab = L_irr */
      dpd_buf4_init(&tau, EOM_CMnEf, C_irr, 3, 8, 3, 8, 0, lbl_s);
      dpd_buf4_mat_irrep_init(&tau, C_irr);
      dpd_buf4_mat_irrep_rd(&tau, C_irr);
      tau_diag = dpd_block_matrix(tau.params->rowtot[C_irr], moinfo.nvirt);
      for(ij=0; ij < tau.params->rowtot[C_irr]; ij++)
	for(Gc=0; Gc < moinfo.nirreps; Gc++)
	  for(C=0; C < moinfo.virtpi[Gc]; C++) {
	    c = C + moinfo.vir_off[Gc];
	    cc = tau.params->colidx[c][c];
	    tau_diag[ij][c] = tau.matrix[C_irr][ij][cc];
	  }
      dpd_buf4_mat_irrep_close(&tau, C_irr);

      dpd_buf4_init(&B_s, CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
      dpd_buf4_init(&S, EOM_TMP, C_irr, 8, 3, 8, 3, 0, "S(ab,ij)");
      dpd_buf4_mat_irrep_init(&S, 0);
      dpd_buf4_mat_irrep_rd(&S, 0);

      rows_per_bucket = dpd_memfree()/(B_s.params->coltot[0] + moinfo.nvirt);
      if(rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
      nbuckets = (int) ceil((double) B_s.params->rowtot[0]/(double) rows_per_bucket);
      rows_left = B_s.params->rowtot[0] % rows_per_bucket;

      B_diag = dpd_block_matrix(rows_per_bucket, moinfo.nvirt);
      next = PSIO_ZERO;
      ncols = tau.params->rowtot[C_irr];
      nlinks = moinfo.nvirt;
      for(m=0; m < (rows_left ? nbuckets-1:nbuckets); m++) {
	row_start = m * rows_per_bucket;
	nrows = rows_per_bucket;
	if(nrows && ncols && nlinks) {
	  psio_read(CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	  C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		  tau_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
	}

      }
      if(rows_left) {
	row_start = m * rows_per_bucket;
	nrows = rows_left;
	if(nrows && ncols && nlinks) {
	  psio_read(CC_BINTS,"B(+) <ab|cc>",(char *) B_diag[0],nrows*nlinks*sizeof(double),next, &next);
	  C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks,
		  tau_diag[0], nlinks, 1, S.matrix[0][row_start], ncols);
	}
      }
      dpd_buf4_mat_irrep_wrt(&S, 0);
      dpd_buf4_mat_irrep_close(&S, 0);
      dpd_buf4_close(&S);
      dpd_buf4_close(&B_s);
      dpd_free_block(B_diag, rows_per_bucket, moinfo.nvirt);
      dpd_free_block(tau_diag, tau.params->rowtot[C_irr], moinfo.nvirt);
      dpd_buf4_close(&tau);

      timer_on("ABCD:A");
      dpd_buf4_init(&tau_a, EOM_CMnEf, C_irr, 4, 9, 4, 9, 0, lbl_a);
      dpd_buf4_init(&B_a, CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
      dpd_buf4_init(&A, EOM_TMP, C_irr, 9, 4, 9, 4, 0, "A(ab,ij)");
      dpd_contract444(&B_a, &tau_a, &A, 0, 0, 0.5, 0);
      dpd_buf4_close(&A);
      dpd_buf4_close(&B_a);
      dpd_buf4_close(&tau_a);
      timer_off("ABCD:A");

      timer_on("ABCD:axpy");
      dpd_buf4_init(&S, EOM_TMP, C_irr, 5, 0, 8, 3, 0, "S(ab,ij)");
      dpd_buf4_sort_axpy(&S, EOM_SIjAb, rspq, 0, 5, SIjAb_lbl, 1);
      dpd_buf4_close(&S);
      dpd_buf4_init(&A, EOM_TMP, C_irr, 5, 0, 9, 4, 0, "A(ab,ij)");
      dpd_buf4_sort_axpy(&A, EOM_SIjAb, rspq, 0, 5, SIjAb_lbl, 1);
      dpd_buf4_close(&A);
      timer_off("ABCD:axpy");
    }

#ifdef TIME_CCEOM
    timer_off("WabefDD Z");
#endif


    /* construct XIjMb = CIjEf * <mb|ef> */
    dpd_buf4_init(&X, EOM_TMP, C_irr, 10, 0, 10, 0, 0, "WabefDD X(Mb,Ij)");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract444(&F, &CMnEf, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 5, 0, 5, 0, 0, "WabefDD Z(Ab,Ij)");
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
    /* fprintf(outfile,"\n begin contract244 in WabefDD\n"); */
    dpd_contract244(&tIA, &X, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_file2_close(&tIA);
    /* dpd_buf4_print(&Z,outfile,1); */
    dpd_buf4_close(&X);

    dpd_buf4_sort_axpy(&Z, EOM_SIjAb, rspq, 0, 5, SIjAb_lbl, -1);
    dpd_buf4_sort_axpy(&Z, EOM_SIjAb, srqp, 0, 5, SIjAb_lbl, -1);

    /* SIjAb += tau_MnAb <Mn||ef> CIjEf */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 0, 0, 0, 0, "WabefDD XIjMn");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&X);
    dpd_buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEF*CIJEF */
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&CMNEF, &B, &SIJAB, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 10, 2, 10, 0, "XIJMA");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_contract444(&CMNEF, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_M");
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
    dpd_contract244(&tIA, &X, &WM, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tIA);
    dpd_buf4_sort(&WM, EOM_TMP, pqsr, 2, 5, "WabefDD_P");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_P");
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&CMNEF, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&X, &tau, &SIJAB, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&X);

    /* Sijab += Wabef*Cijef */
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&Cmnef, &B, &Sijab, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Sijab);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 10, 2, 10, 0, "Xijma");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_contract444(&Cmnef, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_M");
    dpd_file2_init(&tia, CC_OEI, H_IRR, 0, 1, "tia");
    dpd_contract244(&tia, &X, &WM, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tia);
    dpd_buf4_sort(&WM, EOM_TMP, pqsr, 2, 5, "WabefDD_P");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_P");
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&Sijab);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&Cmnef, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauijab");
    dpd_contract444(&X, &tau, &Sijab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&X);

    /* SIjAb += WAbEf*CIjEf */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);

    /* make use of a more efficient algorithm */
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 5, 5, 5, 5, 0, "B <ab|cd>");
    /*  dpd_contract444(&CMnEf, &B, &SIjAb, 0, 0, 1.0, 1.0); */
    dpd_contract444(&B, &CMnEf, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort(&Z, EOM_TMP, rspq, 0, 5, "Z(Ij,Ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_axpy(&Z, &SIjAb, 1);
    dpd_buf4_close(&Z);

    dpd_buf4_close(&CMnEf);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 10, 0, 10, 0, "XIjMa");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnfE");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_sort(&X, EOM_TMP, pqsr, 0, 11, "XIjaM");
    dpd_buf4_close(&X);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 11, 0, 11, 0, "XIjaM");
    dpd_file2_init(&tia, CC_OEI, H_IRR, 0, 1, "tia");
    dpd_contract424(&X, &tia, &SIjAb, 3, 0, 0, -1.0, 1.0);
    dpd_file2_close(&tia);
    dpd_buf4_close(&X);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 10, 0, 10, 0, "XIjMb");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
    dpd_contract244(&tIA, &X, &SIjAb, 0, 2, 1, -1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&X);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 0, 0, 0, 0, "XIjMn");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&SIjAb);
    dpd_buf4_close(&X);

    sprintf(SIA_lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, SIA_lbl);
    sprintf(Sia_lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, Sia_lbl);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    c_clean(&SIA,&Sia,&SIJAB,&Sijab,&SIjAb);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) { /* UHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEF*CIJEF */
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <AB|CD>");
    dpd_contract444(&CMNEF, &B, &SIJAB, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 20, 2, 20, 0, "XIJMA");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_contract444(&CMNEF, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WABEFDD_M");
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
    dpd_contract244(&tIA, &X, &WM, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tIA);
    dpd_buf4_sort(&WM, EOM_TMP, pqsr, 2, 5, "WABEFDD_P");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WABEFDD_P");
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_contract444(&CMNEF, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&X, &tau, &SIJAB, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&X);

    /* Sijab += Wabef*Cijef */
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, Cmnef_lbl);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 17, 12, 17, 0, Sijab_lbl);
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 17, 17, 15, 15, 1, "B <ab|cd>");
    dpd_contract444(&Cmnef, &B, &Sijab, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Sijab);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 12, 30, 12, 30, 0, "Xijma");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_contract444(&Cmnef, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WabefDD_M");
    dpd_file2_init(&tia, CC_OEI, H_IRR, 2, 3, "tia");
    dpd_contract244(&tia, &X, &WM, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tia);
    dpd_buf4_sort(&WM, EOM_TMP, pqsr, 12, 15, "WabefDD_P");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 15, 12, 17, 0, Sijab_lbl);
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WabefDD_P");
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&Sijab);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 12, 12, 12, 12, 0, "Xijmn");
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, Cmnef_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&Cmnef, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 17, 12, 17, 0, Sijab_lbl);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 12, 17, 12, 17, 0, "tauijab");
    dpd_contract444(&X, &tau, &Sijab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&X);

    /* SIjAb += WAbEf*CIjEf */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);

    /* make use of a more efficient algorithm */
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    /*  dpd_contract444(&CMnEf, &B, &SIjAb, 0, 0, 1.0, 1.0); */
    dpd_contract444(&B, &CMnEf, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort(&Z, EOM_TMP, rspq, 22, 28, "Z(Ij,Ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
    dpd_buf4_axpy(&Z, &SIjAb, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&CMnEf);

    dpd_buf4_init(&X, EOM_TMP, C_irr, 22, 27, 22, 27, 0, "XIjmA");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 22, 29, 22, 29, 0, "CMnfE");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_sort(&X, EOM_TMP, pqsr, 22, 26, "XIjAm");
    dpd_buf4_close(&X);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 22, 26, 22, 26, 0, "XIjAm");
    dpd_file2_init(&tia, CC_OEI, H_IRR, 2, 3, "tia");
    dpd_contract424(&X, &tia, &SIjAb, 3, 0, 0, -1.0, 1.0);
    dpd_file2_close(&tia);
    dpd_buf4_close(&X);

    dpd_buf4_init(&X, EOM_TMP, C_irr, 22, 24, 22, 24, 0, "XIjMb");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
    dpd_contract244(&tIA, &X, &SIjAb, 0, 2, 1, -1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&X);

    /* Sijab += tau_mneb <mn||ef> C_ijef */
    dpd_buf4_init(&X, EOM_TMP, C_irr, 22, 22, 22, 22, 0, "XIjMn");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&SIjAb);
    dpd_buf4_close(&X);
  }

#ifdef EOM_DEBUG
  check_sum("WabefDD",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
