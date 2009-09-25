/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar doubles-singles block contribution
   of Wbmfe to a Sigma vector stored at Sigma plus 'i' */

void WbmfeDS(int i, int C_irr) {
  dpdfile2 CME, Cme, XBF, Xbf;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF, WM, WP, W, Z;
  dpdbuf4 TIJAB, TIjAb, Tijab;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];
  int Gbm, Gfe, bm, b, m, Gb, Gm, Ge, Gf, B, M, f, e, fe, ef, nrows, ncols;
  dpdfile2 C;
  double *X;

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
/*     dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF"); */
/*     dpd_file2_scm(&XBF, 0.0); */
/*     dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl); */
/*     dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); */
/*     dpd_dot24(&CME, &W, &XBF, 0, 0, 1.0, 1.0); */
/*     dpd_buf4_close(&W); */
/*     dpd_file2_close(&CME); */

    /* OOC code below added 7/27/05, -TDC */
    /* X(b,f) = [ 2 Wbmfe - Wbmef ] * C(m,e) */
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_file2_scm(&XBF, 0.0);
    dpd_file2_mat_init(&XBF);
    dpd_file2_mat_rd(&XBF);
    dpd_file2_init(&C, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_mat_init(&C);
    dpd_file2_mat_rd(&C);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf");
    for(Gbm=0; Gbm < moinfo.nirreps; Gbm++) {
      Gfe = Gbm ^ H_IRR;
      dpd_buf4_mat_irrep_row_init(&W, Gbm);
      X = init_array(W.params->coltot[Gfe]);
      for(bm=0; bm < W.params->rowtot[Gbm]; bm++) {
	dpd_buf4_mat_irrep_row_rd(&W, Gbm, bm);

	b = W.params->roworb[Gbm][bm][0];
	m = W.params->roworb[Gbm][bm][1];
	Gb = W.params->psym[b];
	Gm = Gbm ^ Gb;
	Ge = Gm ^ C_irr;
	Gf = Ge ^ Gfe;
	B = b - moinfo.vir_off[Gb];
	M = m - moinfo.occ_off[Gm];

	zero_arr(X, W.params->coltot[Gfe]);

	for(fe=0; fe < W.params->coltot[Gfe]; fe++) {
	  f = W.params->colorb[Gfe][fe][0];
	  e = W.params->colorb[Gfe][fe][1];
	  ef = W.params->colidx[e][f];
	  X[fe] = 2.0 * W.matrix[Gbm][0][fe] - W.matrix[Gbm][0][ef];
	}

	nrows = moinfo.virtpi[Gf];
	ncols = moinfo.virtpi[Ge];

	if(nrows && ncols)
	  C_DGEMV('n',nrows,ncols,1,&X[W.col_offset[Gfe][Gf]],ncols,
		  C.matrix[Gm][M],1,1,XBF.matrix[Gb][B],1);

      }
      free(X);
      dpd_buf4_mat_irrep_row_close(&W, Gbm);
    }
    dpd_buf4_close(&W);
    dpd_file2_close(&C);
    dpd_file2_mat_wrt(&XBF);
    dpd_file2_mat_close(&XBF);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WbmfeDS Z(Ij,Ab)");
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract424(&TIjAb, &XBF, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&TIjAb);
    dpd_buf4_sort_axpy(&Z, EOM_SIjAb, qpsr, 0, 5, SIjAb_lbl, 1);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);
    dpd_file2_close(&XBF);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
    /* XBF = CME * WBMFE + Cme * WBmFe */
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_file2_scm(&XBF, 0.0);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WAMEF, CC_HBAR, H_IRR, 11, 5, 11, 7, 0, "WAMEF");
    dpd_dot24(&CME, &WAMEF, &XBF, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAMEF);
    dpd_file2_close(&CME);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_buf4_init(&WAmEf, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf");
    dpd_dot24(&Cme, &WAmEf, &XBF, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAmEf);
    dpd_file2_close(&Cme);
    dpd_file2_close(&XBF);

    /* Xbf = Cme * Wbmfe + CME * WbMfE */
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 1, 1, "Xbf");
    dpd_file2_scm(&Xbf, 0.0);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_buf4_init(&Wamef, CC_HBAR, H_IRR, 11, 5, 11, 7, 0, "Wamef");
    dpd_dot24(&Cme, &Wamef, &Xbf, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Wamef);
    dpd_file2_close(&Cme);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WaMeF, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WaMeF");
    dpd_dot24(&CME, &WaMeF, &Xbf, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WaMeF);
    dpd_file2_close(&CME);
    dpd_file2_close(&Xbf);

    /* SIJAB += XBF * TIJAF - XAF * TIJBF */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_buf4_init(&TIJAB, CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract424(&TIJAB, &XBF, &WP, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&TIJAB);
    dpd_file2_close(&XBF);
    dpd_buf4_sort(&WP, EOM_TMP, pqsr, 2, 5, "WbmfeDS_M"); 
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&SIJAB);

    /* Sijab += Xbf * Tijaf - Xaf * Tijbf */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 1, 1, "Xbf");
    dpd_buf4_init(&Tijab, CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tijab");
    dpd_contract424(&Tijab, &Xbf, &WP, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&Tijab);
    dpd_file2_close(&Xbf);
    dpd_buf4_sort(&WP, EOM_TMP, pqsr, 2, 5, "WbmfeDS_M");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&Sijab);

    /* SIjAb += Xbf * tIjAf + XAF * TIjbF */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 1, 1, "Xbf");
    dpd_contract424(&TIjAb, &Xbf, &SIjAb, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&Xbf);
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_contract244(&XBF, &TIjAb, &SIjAb, 1, 2, 1, 1.0, 1.0);
    dpd_file2_close(&XBF);
    dpd_buf4_close(&TIjAb);
    dpd_buf4_close(&SIjAb);
  }

  else { /* UHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
    /* XBF = CME * WBMFE + Cme * WBmFe */
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_file2_scm(&XBF, 0.0);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WAMEF, CC_HBAR, H_IRR, 21, 5, 21, 7, 0, "WAMEF");
    dpd_dot24(&CME, &WAMEF, &XBF, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAMEF);
    dpd_file2_close(&CME);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_buf4_init(&WAmEf, CC_HBAR, H_IRR, 26, 28, 26, 28, 0, "WAmEf");
    dpd_dot24(&Cme, &WAmEf, &XBF, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAmEf);
    dpd_file2_close(&Cme);
/*
fprintf(outfile,"XBF self dot %15.10lf\n", dpd_file2_dot_self(&XBF));
*/
    dpd_file2_close(&XBF);

    /* Xbf = Cme * Wbmfe + CME * WbMfE */
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 3, 3, "Xbf");
    dpd_file2_scm(&Xbf, 0.0);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_buf4_init(&Wamef, CC_HBAR, H_IRR, 31, 15, 31, 17, 0, "Wamef");
    dpd_dot24(&Cme, &Wamef, &Xbf, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Wamef);
    dpd_file2_close(&Cme);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WaMeF, CC_HBAR, H_IRR, 25, 29, 25, 29, 0, "WaMeF");
    dpd_dot24(&CME, &WaMeF, &Xbf, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WaMeF);
    dpd_file2_close(&CME);
/*
fprintf(outfile,"Xbf self dot %15.10lf\n", dpd_file2_dot_self(&Xbf));
*/
    dpd_file2_close(&Xbf);

    /* SIJAB += XBF * TIJAF - XAF * TIJBF */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_buf4_init(&TIJAB, CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract424(&TIJAB, &XBF, &WP, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&TIJAB);
    dpd_file2_close(&XBF);
    dpd_buf4_sort(&WP, EOM_TMP, pqsr, 2, 5, "WbmfeDS_M"); 
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&SIJAB);

    /* Sijab += Xbf * Tijaf - Xaf * Tijbf */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WbmfeDS_PB");
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 3, 3, "Xbf");
    dpd_buf4_init(&Tijab, CC_TAMPS, H_IRR, 12, 15, 12, 17, 0, "tijab");
    dpd_contract424(&Tijab, &Xbf, &WP, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&Tijab);
    dpd_file2_close(&Xbf);
    dpd_buf4_sort(&WP, EOM_TMP, pqsr, 12, 15, "WbmfeDS_MB");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 15, 12, 17, 0, Sijab_lbl);
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WbmfeDS_MB");
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&Sijab);

    /* SIjAb += Xbf * tIjAf + XAF * TIjbF */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 3, 3, "Xbf");
    dpd_contract424(&TIjAb, &Xbf, &SIjAb, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&Xbf);
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_contract244(&XBF, &TIjAb, &SIjAb, 1, 2, 1, 1.0, 1.0);
    dpd_file2_close(&XBF);
    dpd_buf4_close(&TIjAb);
    dpd_buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WbmfeDS",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
