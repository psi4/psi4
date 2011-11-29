/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <libpsio/psio.h>
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
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
      S0 = dpd_file2_dot(&FME, &CME);
      dpd_file2_close(&FME);
      dpd_file2_close(&CME);
    }
    else {
      S0 = 0.0;
    }

    sprintf(lbl, "%s %d", "S0", i);
    psio_read_entry(EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));
    S0 += S0_old;
    psio_write_entry(EOM_SIA, lbl, (char *) &(S0), sizeof(double));
  }

#ifdef EOM_DEBUG
  check_sum("Sigma0S",i,C_irr);
#endif
  return;
}

void sigma00(int i, int C_irr) {
  char lbl[32];
		double C0, S0, S0_old, reference_expectation_value;

  psio_read_entry(CC_HBAR, "Reference expectation value",
		(char *) &(reference_expectation_value), sizeof(double));

  sprintf(lbl, "%s %d", "C0", i);
  psio_read_entry(EOM_CME, lbl, (char *) &(C0), sizeof(double));
  sprintf(lbl, "%s %d", "S0", i);
  psio_read_entry(EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));
	S0_old += C0 * reference_expectation_value;
  psio_write_entry(EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));

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
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      S0 = dpd_buf4_dot(&D, &CMnEf);
      dpd_buf4_close(&D);
      dpd_buf4_close(&CMnEf);
    }
    else {
      S0 = 0.0;
    }

    sprintf(lbl, "%s %d", "S0", i);
    psio_read_entry(EOM_SIA, lbl, (char *) &(S0_old), sizeof(double));
    S0 += S0_old;
    psio_write_entry(EOM_SIA, lbl, (char *) &(S0), sizeof(double));
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

  psio_read_entry(CC_HBAR, "Reference expectation value",
		(char *) &(reference_expectation_value), sizeof(double));

  if (params.eom_ref == 0) { /* RHF */
  /* SIA += RIA * <0|Hbar|0> */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_axpy(&CME, &SIA, reference_expectation_value,0);
    dpd_file2_close(&CME);
    dpd_file2_close(&SIA);
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

  psio_read_entry(CC_HBAR, "Reference expectation value",
    (char *) &(reference_expectation_value), sizeof(double));

  if (params.eom_ref == 0) { /* RHF */
  /* Sijab += Rijab * <0|Hbar|0> */
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_axpy(&CMnEf, &SIjAb, reference_expectation_value);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&SIjAb);
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

  psio_read_entry(CC_HBAR, "Reference expectation value",
    (char *) &(reference_expectation_value), sizeof(double));

  if (C_irr == H_IRR) {
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    dpd_file2_init(&FAI, CC_OEI, H_IRR, 0, 1, "FAI residual");
		dpd_file2_axpy(&FAI, &SIA, reference_expectation_value, 0);
	  dpd_file2_close(&FAI);
	  dpd_file2_close(&SIA);
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

  psio_read_entry(CC_HBAR, "Reference expectation value",
    (char *) &(reference_expectation_value), sizeof(double));

  if (C_irr == H_IRR) {
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_init(&WAbIj, CC_HBAR, H_IRR, 0, 5, 0, 5, 0, "WAbIj residual");
		dpd_buf4_axpy(&WAbIj, &SIjAb, reference_expectation_value);
	  dpd_buf4_close(&WAbIj);
	  dpd_buf4_close(&SIjAb);
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
    dpd_file2_init(&CIA, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&CIA);
	  dpd_file2_mat_rd(&CIA);
    dpd_file2_init(&FBJ, CC_OEI, H_IRR, 0, 1, "FAI residual");
	  dpd_file2_mat_init(&FBJ);
	  dpd_file2_mat_rd(&FBJ);

    sprintf(lbl, "%s %d", "SIjAb", i_root);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&SIjAb, h);
      dpd_buf4_mat_irrep_rd(&SIjAb, h);

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

          if((Isym^Asym == C_irr) && (Jsym == Bsym))
            SIjAb.matrix[h][row][col] += (CIA.matrix[Isym][I][A] * FBJ.matrix[Jsym][J][B]);

          if((Isym == Asym) && (Jsym^Bsym == C_irr))
            SIjAb.matrix[h][row][col] += (CIA.matrix[Jsym][J][B] * FBJ.matrix[Isym][I][A]);
        }
      }
      dpd_buf4_mat_irrep_wrt(&SIjAb, h);
      dpd_buf4_mat_irrep_close(&SIjAb, h);
    }
    dpd_buf4_close(&SIjAb);

    dpd_file2_mat_close(&FBJ);
    dpd_file2_close(&FBJ);
    dpd_file2_mat_close(&CIA);
    dpd_file2_close(&CIA);
  }
#ifdef EOM_DEBUG
  check_sum("SigmaDS_full",i_root,C_irr);
#endif
	return;
}


}} // namespace psi::cceom
