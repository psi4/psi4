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
#include <cmath>
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar reference-singles block contribution
   to a Sigma vector stored at Sigma plus 'i' */

void sigma0S(int i, int C_irr) {
  dpdfile2 FME, Fme;
  dpdfile2 CME, Cme;
  dpdbuf4 W;
  char lbl[32];
  double S0, S0_old;

  /* S0 += FME*CME */
  if (params.eom_ref == 0) { /* RHF */
    if (C_irr == H_IRR) {
      sprintf(lbl, "%s %d", "CME", i);
      global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
      global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
      S0 = global_dpd_->file2_dot(&FME, &CME);
      global_dpd_->file2_close(&FME);
      global_dpd_->file2_close(&CME);
    }
    else {
      S0 = 0.0;
    }

    sprintf(lbl, "%s %d", "S0", i);
    psio_read_entry(PSIF_EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));
    S0 += S0_old;
    psio_write_entry(PSIF_EOM_SIA, lbl, (char *) &(S0), sizeof(double));
  }

#ifdef EOM_DEBUG
  check_sum("Sigma0S",i,C_irr);
#endif
  return;
}

void sigma00(int i, int C_irr) {
  char lbl[32];
        double C0, S0, S0_old, reference_expectation_value;

  psio_read_entry(PSIF_CC_HBAR, "Reference expectation value",
        (char *) &(reference_expectation_value), sizeof(double));

  sprintf(lbl, "%s %d", "C0", i);
  psio_read_entry(PSIF_EOM_CME, lbl, (char *) &(C0), sizeof(double));
  sprintf(lbl, "%s %d", "S0", i);
  psio_read_entry(PSIF_EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));
    S0_old += C0 * reference_expectation_value;
  psio_write_entry(PSIF_EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));

#ifdef EOM_DEBUG
  check_sum("Sigma00",i,C_irr);
#endif
  return;
}


/* This function computes the H-bar reference-doubles block contribution
   to a Sigma vector stored at Sigma plus 'i' */

void sigma0D(int i, int C_irr) {
  dpdbuf4 D;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CmNeF;
  char lbl[32];
  double S0, S0_old;

  /* S0 += Cijab <ij||ab> */
  if (params.eom_ref == 0) { /* RHF */
    if (C_irr == H_IRR) {
      sprintf(lbl, "%s %d", "CMnEf", i);
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      S0 = global_dpd_->buf4_dot(&D, &CMnEf);
      global_dpd_->buf4_close(&D);
      global_dpd_->buf4_close(&CMnEf);
    }
    else {
      S0 = 0.0;
    }

    sprintf(lbl, "%s %d", "S0", i);
    psio_read_entry(PSIF_EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));
    S0 += S0_old;
    psio_write_entry(PSIF_EOM_SIA, lbl, (char *) &(S0), sizeof(double));
  }
#ifdef EOM_DEBUG
  check_sum("Sigma0D",i,C_irr);
#endif
  return;
}


/* This function adds unconnected terms to the H-bar singles-singles
block contributions to a sigma vector stored at Sigma plus 'i' */

void sigmaSS_full(int i, int C_irr) {
  dpdfile2 CME, SIA;
  char lbl[32];
    double reference_expectation_value;

  psio_read_entry(PSIF_CC_HBAR, "Reference expectation value",
        (char *) &(reference_expectation_value), sizeof(double));

  if (params.eom_ref == 0) { /* RHF */
  /* SIA += RIA * <0|Hbar|0> */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_axpy(&CME, &SIA, reference_expectation_value,0);
    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&SIA);
  }

#ifdef EOM_DEBUG
  check_sum("SigmaSS_full",i,C_irr);
#endif
  return;
}


/* This function adds unconnected terms to the H-bar doubles-doubles
block contributions to a sigma vector stored at Sigma plus 'i' */

void sigmaDD_full(int i, int C_irr) {
  dpdbuf4 CMnEf, SIjAb;
  char lbl[32];
  double reference_expectation_value;

  psio_read_entry(PSIF_CC_HBAR, "Reference expectation value",
    (char *) &(reference_expectation_value), sizeof(double));

  if (params.eom_ref == 0) { /* RHF */
  /* Sijab += Rijab * <0|Hbar|0> */
    sprintf(lbl, "%s %d", "SIjAb", i);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_axpy(&CMnEf, &SIjAb, reference_expectation_value);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&SIjAb);
  }
#ifdef EOM_DEBUG
  check_sum("SigmaDD_full",i,C_irr);
#endif
  return;
}


/* This function computes the unconnected terms to the <S|Hbar|0> contributions
to the sigma vector stored at sigma i */
void sigmaS0(int i, int C_irr) {
  dpdfile2 FAI, SIA;
    double reference_expectation_value;
    char lbl[32];

  psio_read_entry(PSIF_CC_HBAR, "Reference expectation value",
    (char *) &(reference_expectation_value), sizeof(double));

  if (C_irr == H_IRR) {
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    global_dpd_->file2_init(&FAI, PSIF_CC_OEI, H_IRR, 0, 1, "FAI residual");
        global_dpd_->file2_axpy(&FAI, &SIA, reference_expectation_value, 0);
      global_dpd_->file2_close(&FAI);
      global_dpd_->file2_close(&SIA);
    }
#ifdef EOM_DEBUG
  check_sum("SigmaS0",i,C_irr);
#endif
    return;
}

/* This function computes the unconnected terms to the <D|Hbar|0> contributions
to the sigma vector stored at sigma i */
void sigmaD0(int i, int C_irr) {
  dpdbuf4 WAbIj, SIjAb;
    double reference_expectation_value;
    char lbl[32];

  psio_read_entry(PSIF_CC_HBAR, "Reference expectation value",
    (char *) &(reference_expectation_value), sizeof(double));

  if (C_irr == H_IRR) {
    sprintf(lbl, "%s %d", "SIjAb", i);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&WAbIj, PSIF_CC_HBAR, H_IRR, 0, 5, 0, 5, 0, "WAbIj residual");
        global_dpd_->buf4_axpy(&WAbIj, &SIjAb, reference_expectation_value);
      global_dpd_->buf4_close(&WAbIj);
      global_dpd_->buf4_close(&SIjAb);
    }
#ifdef EOM_DEBUG
  check_sum("SigmaD0",i,C_irr);
#endif
    return;
}

/* computes unconnected contributions of the <D|Hbar|S> block to
the sigma vector stored at sigma i */
/* Sijab += P(ij) P(ab) Cia Fbj, that is, [<Phi_j^b|Hbar|0>] */

void sigmaDS_full(int i_root, int C_irr) {
  dpdbuf4 SIjAb;
  dpdfile2 FBJ, CIA;
  int h, nirreps;
  int row,col;
  int i,j,a,b,I,J,A,B,Isym,Jsym,Asym,Bsym;
    char lbl[32];

  nirreps = moinfo.nirreps;

  /* SIjAb += CIA * Fjb + Cjb * FIA */
  if (params.ref == 0) {
    sprintf(lbl, "%s %d", "CME", i_root);
    global_dpd_->file2_init(&CIA, PSIF_EOM_CME, C_irr, 0, 1, lbl);
      global_dpd_->file2_mat_init(&CIA);
      global_dpd_->file2_mat_rd(&CIA);
    global_dpd_->file2_init(&FBJ, PSIF_CC_OEI, H_IRR, 0, 1, "FAI residual");
      global_dpd_->file2_mat_init(&FBJ);
      global_dpd_->file2_mat_rd(&FBJ);

    sprintf(lbl, "%s %d", "SIjAb", i_root);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&SIjAb, h);
      global_dpd_->buf4_mat_irrep_rd(&SIjAb, h);

      for(row=0; row < SIjAb.params->rowtot[h]; row++) {
        i = SIjAb.params->roworb[h][row][0];
        j = SIjAb.params->roworb[h][row][1];

        for(col=0; col < SIjAb.params->coltot[h^C_irr]; col++) {
          a = SIjAb.params->colorb[h^C_irr][col][0];
          b = SIjAb.params->colorb[h^C_irr][col][1];

          I = CIA.params->rowidx[i]; Isym = CIA.params->psym[i];
          J = FBJ.params->rowidx[j]; Jsym = FBJ.params->psym[j];
          A = CIA.params->colidx[a]; Asym = CIA.params->qsym[a];
          B = FBJ.params->colidx[b]; Bsym = FBJ.params->qsym[b];

          // ^ has lower precedence than ==; == will be evaluated first
          // original: if((Isym^Asym == C_irr) && (Jsym == Bsym))
          if(((Isym^Asym) == C_irr) && (Jsym == Bsym))
            SIjAb.matrix[h][row][col] += (CIA.matrix[Isym][I][A] * FBJ.matrix[Jsym][J][B]);

          // ^ has lower precedence than ==; == will be evaluated first
          //if((Isym == Asym) && (Jsym^Bsym == C_irr))
          if((Isym == Asym) && ((Jsym^Bsym) == C_irr))
            SIjAb.matrix[h][row][col] += (CIA.matrix[Jsym][J][B] * FBJ.matrix[Isym][I][A]);
        }
      }
      global_dpd_->buf4_mat_irrep_wrt(&SIjAb, h);
      global_dpd_->buf4_mat_irrep_close(&SIjAb, h);
    }
    global_dpd_->buf4_close(&SIjAb);

    global_dpd_->file2_mat_close(&FBJ);
    global_dpd_->file2_close(&FBJ);
    global_dpd_->file2_mat_close(&CIA);
    global_dpd_->file2_close(&CIA);
  }
#ifdef EOM_DEBUG
  check_sum("SigmaDS_full",i_root,C_irr);
#endif
    return;
}


}} // namespace psi::cceom
