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

/* This function orthogonalizes the r residual vector with the set of
numCs C vectors and adds the new vector to the C list if its norm is greater
than params.residual_tol */

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
  dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
  dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);

/* use for ROHF and UHF */
void schmidt_add(dpdfile2 *RIA, dpdfile2 *Ria,
  dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int *numCs, int irrep)
{
  double dotval;
  double norm;
  int i, I;
  dpdfile2 Cme, CME, Cme2, CME2;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CMNEF2, Cmnef2, CMnEf2;
  dpdbuf4 CMnEf_buf;
  char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];

  for (i=0; i<*numCs; i++) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);

    dpd_file2_init(&CME, EOM_CME, irrep, 0, 1, CME_lbl);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, irrep, 2, 7, 2, 7, 0, CMNEF_lbl);
    if (params.eom_ref == 1) {
      dpd_file2_init(&Cme, EOM_Cme, irrep, 0, 1, Cme_lbl);
      dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 2, 7, 2, 7, 0, Cmnef_lbl);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
    }
    else if (params.eom_ref == 2) {
      dpd_file2_init(&Cme, EOM_Cme, irrep, 2, 3, Cme_lbl);
      dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 12, 17, 12, 17, 0, Cmnef_lbl);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 22, 28, 22, 28, 0, CMnEf_lbl);
    }

    dotval  = dpd_file2_dot(RIA, &CME);
    dotval += dpd_file2_dot(Ria, &Cme);
  //fprintf(outfile, "OE Dotval for vector %d = %20.14f\n", i, dotval);
    dotval += dpd_buf4_dot(RIJAB, &CMNEF);
    dotval += dpd_buf4_dot(Rijab, &Cmnef);
    dotval += dpd_buf4_dot(RIjAb, &CMnEf);

  //fprintf(outfile, "Dotval for vector %d = %20.14f\n", i, dotval);

    dpd_file2_axpy(&CME, RIA, -1.0*dotval, 0);
    dpd_file2_axpy(&Cme, Ria, -1.0*dotval, 0);
    dpd_buf4_axpy(&CMNEF, RIJAB, -1.0*dotval);
    dpd_buf4_axpy(&Cmnef, Rijab, -1.0*dotval);
    dpd_buf4_axpy(&CMnEf, RIjAb, -1.0*dotval);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_close(&CMnEf);
  }

  norm = norm_C(RIA, Ria, RIJAB, Rijab, RIjAb);
  //fprintf(outfile, "Norm of residual (TDC) = %20.14f\n", norm);

  if (norm < eom_params.schmidt_add_residual_tol) {
    return;
  }
  else {
    scm_C(RIA, Ria, RIJAB, Rijab, RIjAb, 1.0/norm);
    sprintf(CME_lbl, "%s %d", "CME", *numCs);
    sprintf(Cme_lbl, "%s %d", "Cme", *numCs);
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", *numCs);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", *numCs);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", *numCs);

    dpd_file2_copy(RIA, EOM_CME, CME_lbl);
    dpd_file2_copy(Ria, EOM_Cme, Cme_lbl);
    dpd_buf4_copy(RIJAB, EOM_CMNEF, CMNEF_lbl);
    dpd_buf4_copy(Rijab, EOM_Cmnef, Cmnef_lbl);
    dpd_buf4_copy(RIjAb, EOM_CMnEf, CMnEf_lbl);

    ++(*numCs);
  }
  return;
}

void schmidt_add_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int *numCs, int irrep)
{
  double dotval, norm, R0, C0;
  int i, I;
  dpdfile2 CME;
  dpdbuf4 CMnEf, CAB1, CAB2;
  dpdfile2 R1;
  dpdbuf4 R2a, R2b;
  char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32], C0_lbl[32];

  if (params.full_matrix) psio_read_entry(EOM_R, "R0", (char *) &R0, sizeof(double));

  for (i=0; i<*numCs; i++) {
    /* Spin-adapt the residual */
    dpd_buf4_copy(RIjAb, EOM_TMP, "RIjAb");
    dpd_buf4_sort(RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");

    dpd_buf4_init(&R2a, EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_init(&R2b, EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");
    dpd_buf4_scm(&R2a, 2.0);
    dpd_buf4_axpy(&R2b, &R2a, -1.0);
    dpd_buf4_close(&R2b);

    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    dpd_file2_init(&CME, EOM_CME, irrep, 0, 1, CME_lbl);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
    dotval  = 2.0 * dpd_file2_dot(RIA, &CME);
 //fprintf(outfile, "OE Dotval for vector %d = %20.14f\n", i, dotval);
    dotval += dpd_buf4_dot(&R2a, &CMnEf);
    dpd_buf4_close(&R2a);
		if (params.full_matrix) {
      sprintf(C0_lbl, "%s %d", "C0", i);
			psio_read_entry(EOM_CME, C0_lbl, (char *) &C0, sizeof(double));
			dotval += C0 * R0;
		}

 //fprintf(outfile, "Dotval for vector %d = %20.14f\n", i, dotval);
		R0 = R0 - 1.0 * dotval * C0;
    dpd_file2_axpy(&CME, RIA, -1.0*dotval, 0);
    dpd_buf4_axpy(&CMnEf, RIjAb, -1.0*dotval);
    dpd_file2_close(&CME);
    dpd_buf4_close(&CMnEf);
  }

  dpd_buf4_sort(RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
  dpd_buf4_init(&R2b, EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");

  /* norm = norm_C_rhf(RIA, RIjAb, &R2b); */
  norm  = 2.0 * dpd_file2_dot_self(RIA);
  norm += 2.0 * dpd_buf4_dot_self(RIjAb);
  norm -= dpd_buf4_dot(RIjAb, &R2b);
	if (params.full_matrix)
	  norm += R0 * R0;
  norm = sqrt(norm);

  dpd_buf4_close(&R2b);

  //fprintf(outfile, "Norm of residual (TDC) = %20.14f\n", norm);

  if (norm < eom_params.schmidt_add_residual_tol) {
    return;
  }
  else {

    if (params.full_matrix) R0 *= 1.0/norm;
    dpd_file2_scm(RIA, 1.0/norm);
    dpd_buf4_scm(RIjAb, 1.0/norm);

#ifdef EOM_DEBUG
    dpd_buf4_sort(RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
    dpd_buf4_init(&R2b, EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");
    norm  = 2.0 * dpd_file2_dot_self(RIA);
    norm += 2.0 * dpd_buf4_dot_self(RIjAb);
    norm -= dpd_buf4_dot(RIjAb, &R2b);
		if (params.full_matrix) norm += R0 * R0;
    norm = sqrt(norm);
    fprintf(outfile,"Norm of final new C in schmidt_add(): %20.15lf\n", norm);
    dpd_buf4_close(&R2b);
#endif

    sprintf(CME_lbl, "%s %d", "CME", *numCs);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", *numCs);

    dpd_file2_copy(RIA, EOM_CME, CME_lbl);
    dpd_buf4_copy(RIjAb, EOM_CMnEf, CMnEf_lbl);

    /* Generate AA and BB C2 vectors from AB vector */
    /* C(IJ,AB) = C(ij,ab) = C(Ij,Ab) - C(Ij,bA) */
    dpd_buf4_copy(RIjAb, EOM_TMP, "CMnEf");
    dpd_buf4_sort(RIjAb, EOM_TMP, pqsr, 0, 5, "CMnfE");

    dpd_buf4_init(&CAB1, EOM_TMP, irrep, 0, 5, 0, 5, 0, "CMnEf");
    dpd_buf4_init(&CAB2, EOM_TMP, irrep, 0, 5, 0, 5, 0, "CMnfE");
    dpd_buf4_axpy(&CAB2, &CAB1, -1.0);
    dpd_buf4_close(&CAB2);
    dpd_buf4_close(&CAB1);

		if (params.full_matrix) {
      sprintf(C0_lbl, "%s %d", "C0", *numCs);
		  psio_write_entry(EOM_CME, C0_lbl, (char *) &R0, sizeof(double));
		}
    ++(*numCs);
  }
  return;
}

}} // namespace psi::cceom
