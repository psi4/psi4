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
#include <string>
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
            dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double dot_C(dpdfile2 *CME, dpdfile2 *Cme,
            dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);
extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
            dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);

/* this function determines R0, properly normalizes R, and checks orthogonality
 * with the ground state left eigenvector (1+lambda) */

/* for ROHF and UHF */
void rzero(int C_irr, int *converged) {
  double rzero=0.0, energy, norm, dotval;
  double dot_IA, dot_ia, dot_IJAB, dot_ijab, dot_IjAb;
  dpdfile2 RIA, Ria, RIA2, Ria2, FIA, Fia, LIA, Lia;
  dpdbuf4 RIJAB, Rijab, RIjAb, D, R2, LIJAB, Lijab, LIjAb;
  dpdbuf4 fRIJAB, fRijab, fRIjAb;
  int L_irr, i;
  int A_OCC, B_OCC, A_VIR, B_VIR;
  int AA_OCC, AA_VIR, BB_OCC, BB_VIR, AB_OCC, AB_VIR;
  char lbl[32], E_lbl[32], R1A_lbl[32], R1B_lbl[32];
  char R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];
  int R_index = -1;

  A_OCC = 0; A_VIR = 1;
  AA_OCC = 2; AA_VIR = 7;
  if (params.eom_ref <= 1) {
    B_OCC = 0; B_VIR = 1;
    BB_OCC = 2; BB_VIR = 7;
    AB_OCC = 0; AB_VIR = 5;
  }
  else if (params.eom_ref == 2) {
    B_OCC = 2; B_VIR = 3;
    BB_OCC = 12; BB_VIR = 17;
    AB_OCC = 22; AB_VIR = 28;
  }
  L_irr = eom_params.L_irr;

  for(i=0; i < eom_params.cs_per_irrep[C_irr]; i++) {
    if (!converged[i]) continue; /* this root did not converged */
    ++R_index;

    if(params.wfn == "EOM_CC2") {
      sprintf(E_lbl, "EOM CC2 Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(PSIF_CC_INFO, E_lbl) == NULL) {
        outfile->Printf("No EOM CC2 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(PSIF_CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      sprintf(E_lbl, "EOM CCSD Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(PSIF_CC_INFO, E_lbl) == NULL) {
        outfile->Printf("No EOM CCSD Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(PSIF_CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      sprintf(E_lbl, "EOM CC3 Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(PSIF_CC_INFO, E_lbl) == NULL) {
        outfile->Printf("No EOM CC3 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(PSIF_CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }

    sprintf(R1A_lbl, "RIA %d %d", C_irr, R_index);
    sprintf(R1B_lbl, "Ria %d %d", C_irr, R_index);
    sprintf(R2AA_lbl, "RIJAB %d %d", C_irr, R_index);
    sprintf(R2BB_lbl, "Rijab %d %d", C_irr, R_index);
    sprintf(R2AB_lbl, "RIjAb %d %d", C_irr, R_index);

    /* Calculate <0| Hbar R |0> */
    if (C_irr == H_IRR) {
      global_dpd_->file2_init(&FIA, PSIF_CC_OEI, H_IRR, A_OCC, A_VIR, "FME");
      global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, A_OCC, A_VIR, R1A_lbl);
      dot_IA = global_dpd_->file2_dot(&FIA, &RIA);
      global_dpd_->file2_close(&RIA);
      global_dpd_->file2_close(&FIA);

      global_dpd_->file2_init(&Fia, PSIF_CC_OEI, H_IRR, B_OCC, B_VIR, "Fme");
      global_dpd_->file2_init(&Ria, PSIF_CC_RAMPS, C_irr, B_OCC, B_VIR, R1B_lbl);
      dot_ia = global_dpd_->file2_dot(&Fia, &Ria);
      global_dpd_->file2_close(&Ria);
      global_dpd_->file2_close(&Fia);

      if (params.eom_ref == 1) {
        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
        global_dpd_->buf4_init(&RIJAB, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
        dot_IJAB = global_dpd_->buf4_dot(&D, &RIJAB);
        global_dpd_->buf4_close(&RIJAB);
        global_dpd_->buf4_init(&Rijab, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2BB_lbl);
        dot_ijab = global_dpd_->buf4_dot(&D, &Rijab);
        global_dpd_->buf4_close(&Rijab);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
        global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
        dot_IjAb = global_dpd_->buf4_dot(&D, &RIjAb);
        global_dpd_->buf4_close(&RIjAb);
        global_dpd_->buf4_close(&D);
      }
      else if (params.eom_ref == 2) {
        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
        global_dpd_->buf4_init(&RIJAB, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
        dot_IJAB = global_dpd_->buf4_dot(&D, &RIJAB);
        global_dpd_->buf4_close(&RIJAB);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
        global_dpd_->buf4_init(&Rijab, PSIF_CC_RAMPS, C_irr, 12, 17, 12, 17, 0, R2BB_lbl);
        dot_ijab = global_dpd_->buf4_dot(&D, &Rijab);
        global_dpd_->buf4_close(&Rijab);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 22, 28, 22, 28, 0, "D <Ij|Ab>");
        global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 22, 28, 22, 28, 0, R2AB_lbl);
        dot_IjAb = global_dpd_->buf4_dot(&D, &RIjAb);
        global_dpd_->buf4_close(&RIjAb);
        global_dpd_->buf4_close(&D);
      }
      rzero = (dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb)/energy;
    }
    else { /* C and H are different irreps */
      rzero = 0.0;
    }

    /* Now normalize so that <R|R> = 1 */
    global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, A_OCC, A_VIR, R1A_lbl);
    global_dpd_->file2_init(&Ria, PSIF_CC_RAMPS, C_irr, B_OCC, B_VIR, R1B_lbl);
    global_dpd_->buf4_init(&fRIJAB, PSIF_CC_RAMPS, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, R2AA_lbl);
    global_dpd_->buf4_init(&fRijab, PSIF_CC_RAMPS, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, R2BB_lbl);
    global_dpd_->buf4_init(&fRIjAb, PSIF_CC_RAMPS, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, R2AB_lbl);

    /* make R0 a positive number */
    /*
    if (rzero < 0.0) {
      rzero *= -1.0;
      dpd_file2_scm(&RIA,-1.0);
      dpd_file2_scm(&Ria,-1.0);
      dpd_buf4_scm(&fRIJAB,-1.0);
      dpd_buf4_scm(&fRijab,-1.0);
      dpd_buf4_scm(&fRIjAb,-1.0);
    }
    */

    /*
    norm = norm_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb);
    norm *= norm;
    */
    norm = dot_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb);
    norm += rzero * rzero;
    norm = sqrt(norm);
    rzero = rzero / norm;
    scm_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb, 1.0/norm);

    norm = dot_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb);
    norm += rzero * rzero;
    outfile->Printf("<R|R> = %20.16lf\n",norm);

/* just debugging with converged solutions - also my need a sort_C() */
/*
dpd_file2_copy(&RIA, EOM_CME, "CME 0");
dpd_file2_copy(&Ria, EOM_Cme, "Cme 0");
dpd_buf4_copy(&fRIJAB, EOM_CMNEF, "CMNEF 0");
dpd_buf4_copy(&fRijab, EOM_Cmnef, "Cmnef 0");
dpd_buf4_copy(&fRIjAb, EOM_CMnEf, "CMnEf 0");
*/
/* end debugging stuff */


    global_dpd_->file2_close(&RIA);
    global_dpd_->file2_close(&Ria);
    global_dpd_->buf4_close(&fRIJAB);
    global_dpd_->buf4_close(&fRijab);
    global_dpd_->buf4_close(&fRIjAb);

    if(params.wfn == "EOM_CC2") {
      outfile->Printf("EOM CC2 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC2 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      outfile->Printf("EOM CCSD R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CCSD R0 for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      outfile->Printf("EOM CC3 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC3 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }

    if (eom_params.dot_with_L) {
      /* evaluate check <R|L> == 0 */
      if (C_irr == L_irr ) {
        global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, A_OCC, A_VIR, R1A_lbl);
        global_dpd_->file2_init(&Ria, PSIF_CC_RAMPS, C_irr, B_OCC, B_VIR, R1B_lbl);
        global_dpd_->buf4_init(&RIJAB, PSIF_CC_RAMPS, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, R2AA_lbl);
        global_dpd_->buf4_init(&Rijab, PSIF_CC_RAMPS, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, R2BB_lbl);
        global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, R2AB_lbl);

        global_dpd_->file2_init(&LIA, PSIF_CC_OEI, L_irr, A_OCC, A_VIR, "LIA");
        global_dpd_->file2_init(&Lia, PSIF_CC_OEI, L_irr, B_OCC, B_VIR, "Lia");
        global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMPS, L_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, "LIJAB");
        global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMPS, L_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, "Lijab");
        global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMPS, L_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, "LIjAb");

        outfile->Printf("\nROHF orthogonality test\n");
        outfile->Printf("<L0|R0>            = %15.10lf\n", eom_params.L0 * rzero);
        dot_IA = global_dpd_->file2_dot(&LIA, &RIA);
        outfile->Printf("<LIA|RIA>          = %15.10lf\n", dot_IA);
        dot_ia = global_dpd_->file2_dot(&Lia, &Ria);
        outfile->Printf("<Lia|Ria>          = %15.10lf\n", dot_ia);
        dot_IJAB = global_dpd_->buf4_dot(&LIJAB, &RIJAB);
        outfile->Printf("<LIJAB|RIJAB>      = %15.10lf\n", dot_IJAB);
        dot_ijab = global_dpd_->buf4_dot(&Lijab, &Rijab);
        outfile->Printf("<Lijab|Rijab>      = %15.10lf\n", dot_ijab);
        dot_IjAb = global_dpd_->buf4_dot(&LIjAb, &RIjAb);
        outfile->Printf("<LIjAb|RIjAb>      = %15.10lf\n", dot_IjAb);
        outfile->Printf("<L|R>              = %15.10lf\n", (eom_params.L0 * rzero)
            + dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb);

        global_dpd_->file2_close(&LIA);
        global_dpd_->file2_close(&Lia);
        global_dpd_->buf4_close(&LIJAB);
        global_dpd_->buf4_close(&Lijab);
        global_dpd_->buf4_close(&LIjAb);
        global_dpd_->file2_close(&RIA);
        global_dpd_->file2_close(&Ria);
        global_dpd_->buf4_close(&RIJAB);
        global_dpd_->buf4_close(&Rijab);
        global_dpd_->buf4_close(&RIjAb);
      }
      else {
        outfile->Printf("\nOverlap <R|L> zero by symmetry\n");
      }
    } /* end dot_with_L loop */
  }
  return;
}

/* normalizes R and produces copies of R that are ROHF-like */
/* sort_amps then produces the others sorted versions of R */

void rzero_rhf(int C_irr, int *converged) {
  double r1, r2, rzero=0.0, energy, norm, dotval;
  double dot_IA, dot_ia, dot_IJAB, dot_ijab, dot_IjAb;
  dpdfile2 RIA, FIA, LIA, Lia, Ria;
  dpdbuf4 RIjAb, RIjbA, RIjAb1, RIjbA1, D, R2, LIjAb, RIJAB, Rijab;
  dpdbuf4 LIJAB, Lijab;
  int L_irr,i;
  char lbl[32], E_lbl[32], R1A_lbl[32], R1B_lbl[32], *blank;
  char R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];
  int R_index=-1;

  L_irr = eom_params.L_irr;

  for(i=0; i < eom_params.cs_per_irrep[C_irr]; i++) {
    if (!converged[i]) continue; /* this root did not converge */
    ++R_index;

    if(params.wfn == "EOM_CC2") {
      sprintf(E_lbl, "EOM CC2 Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(PSIF_CC_INFO, E_lbl) == NULL) {
        outfile->Printf("No EOM CC2 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(PSIF_CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      sprintf(E_lbl, "EOM CCSD Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(PSIF_CC_INFO, E_lbl) == NULL) {
        outfile->Printf("No EOM CCSD Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(PSIF_CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      sprintf(E_lbl, "EOM CC3 Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(PSIF_CC_INFO, E_lbl) == NULL) {
        outfile->Printf("No EOM CC3 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(PSIF_CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }

    sprintf(R1A_lbl, "RIA %d %d", C_irr, R_index);
    sprintf(R1B_lbl, "Ria %d %d", C_irr, R_index);
    sprintf(R2AB_lbl, "RIjAb %d %d", C_irr, R_index);
    sprintf(R2AA_lbl, "RIJAB %d %d", C_irr, R_index);
    sprintf(R2BB_lbl, "Rijab %d %d", C_irr, R_index);

    /* produce RIjbA and 2RIjAb-RIjbA copies - not yet normalized */
    global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    global_dpd_->buf4_sort(&RIjAb, PSIF_CC_TMP, pqsr, 0, 5, "RIjbA");
    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    global_dpd_->buf4_copy(&RIjAb, PSIF_CC_RAMPS, lbl);
    global_dpd_->buf4_close(&RIjAb);

    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_scm(&RIjAb, 2.0);
    global_dpd_->buf4_init(&RIjbA, PSIF_CC_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");
    global_dpd_->buf4_axpy(&RIjbA, &RIjAb, -1.0);
    global_dpd_->buf4_close(&RIjbA);
    global_dpd_->buf4_close(&RIjAb);

    /* calculate <0| hbar | 0> */
    if (!params.full_matrix) {
      if (C_irr == H_IRR) {
        global_dpd_->file2_init(&FIA, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
        global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, 0, 1, R1A_lbl);
        r1 = 2.0 * global_dpd_->file2_dot(&FIA, &RIA);
        global_dpd_->file2_close(&RIA);
        global_dpd_->file2_close(&FIA);

        sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
        global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
        r2 = global_dpd_->buf4_dot(&D, &RIjAb);
        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_close(&RIjAb);
        rzero = (r1 + r2)/energy;
      }
      else {
        rzero = 0.0;
      }
    }
    else { /* full matrix */
      sprintf(lbl, "%s %d %d", "R0", C_irr, R_index);
      psio_read_entry(PSIF_CC_RAMPS, lbl, (char *) &rzero, sizeof(double));
    }

    /* Now normalize R so that <R|R> = 1 */
    global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    global_dpd_->buf4_init(&RIjbA, PSIF_CC_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");

    /*
    if (rzero < 0.0) {
      rzero *= -1.0;
      dpd_file2_scm(&RIA,-1.0);
      dpd_buf4_scm(&RIjAb,-1.0);
      dpd_buf4_scm(&RIjbA,-1.0);
    }
    */

    norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
    norm *= norm;
    norm += rzero * rzero;
    norm = sqrt(norm);
    rzero = rzero / norm;
    global_dpd_->file2_scm(&RIA, 1.0/norm);
    global_dpd_->buf4_scm(&RIjAb, 1.0/norm);
    global_dpd_->buf4_scm(&RIjbA, 1.0/norm);

/* just for debugging cc3 put normalized vector back into C as well */
/*
dpd_file2_copy(&RIA, EOM_CME, "CME 0");
dpd_buf4_copy(&RIjAb, EOM_CMnEf, "CMnEf 0");
*/
/* end debugging stuff */

    global_dpd_->file2_close(&RIA);
    global_dpd_->buf4_close(&RIjAb);
    global_dpd_->buf4_close(&RIjbA);

    if(params.wfn == "EOM_CC2") {
      outfile->Printf("EOM CC2 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC2 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      outfile->Printf("EOM CCSD R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CCSD R0 for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      outfile->Printf("EOM CC3 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC3 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }

    /* produce ROHF like quantities and 2RIjAb-RIjbA */
    global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    global_dpd_->file2_copy(&RIA, PSIF_CC_RAMPS, R1B_lbl);
    global_dpd_->file2_close(&RIA);

    global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    global_dpd_->buf4_copy(&RIjAb, PSIF_CC_RAMPS, lbl);
    global_dpd_->buf4_close(&RIjAb);

    global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 2, 7, 0, 5, 1, R2AB_lbl);
    global_dpd_->buf4_copy(&RIjAb, PSIF_CC_RAMPS, R2AA_lbl);
    global_dpd_->buf4_copy(&RIjAb, PSIF_CC_RAMPS, R2BB_lbl);
    global_dpd_->buf4_close(&RIjAb);

    global_dpd_->buf4_init(&RIjbA, PSIF_CC_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");
    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_scm(&RIjAb, 2.0);
    global_dpd_->buf4_axpy(&RIjbA, &RIjAb, -1.0);
    global_dpd_->buf4_close(&RIjAb);
    global_dpd_->buf4_close(&RIjbA);

    /* test normalization with produced ROHF-like quantities */

    global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    norm = global_dpd_->file2_dot_self(&RIA);
    global_dpd_->file2_close(&RIA);
    global_dpd_->file2_init(&Ria, PSIF_CC_RAMPS, C_irr, 0, 1, R1B_lbl);
    norm += global_dpd_->file2_dot_self(&Ria);
    global_dpd_->file2_close(&Ria);
    global_dpd_->buf4_init(&RIJAB, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += global_dpd_->buf4_dot_self(&RIJAB);
    global_dpd_->buf4_close(&RIJAB);
    global_dpd_->buf4_init(&Rijab, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2BB_lbl);
    norm += global_dpd_->buf4_dot_self(&Rijab);
    global_dpd_->buf4_close(&Rijab);
    global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    norm += global_dpd_->buf4_dot_self(&RIjAb);
    global_dpd_->buf4_close(&RIjAb);
    norm += rzero * rzero;

#ifdef EOM_DEBUG
    outfile->Printf("Norm with produced ROHF-like quantities = %15.10lf\n", norm);
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    dpd_file2_print(&RIA, outfile);
    dpd_file2_close(&RIA);
#endif

    /* check orthogonality with L */
    if (eom_params.dot_with_L) {
      if (C_irr == L_irr) {
        global_dpd_->file2_init(&LIA, PSIF_CC_OEI, L_irr, 0, 1, "LIA");
        global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, 0, 1, R1A_lbl);
        r1 = 2.0 * global_dpd_->file2_dot(&LIA, &RIA);
        global_dpd_->file2_close(&RIA);
        global_dpd_->file2_close(&LIA);

        global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
        sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
        global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
        r2 = global_dpd_->buf4_dot(&LIjAb, &RIjAb);
        global_dpd_->buf4_close(&RIjAb);
        global_dpd_->buf4_close(&LIjAb);

        dotval = r1 + r2 + (eom_params.L0 * rzero);
        outfile->Printf("Performing RHF orthogonality test\n");
        outfile->Printf("<L0|R0>              = %15.10lf\n", eom_params.L0 * rzero);
        outfile->Printf("2*<LIA|RIA>          = %15.10lf\n", r1);
        outfile->Printf("<LIjAb|2RIjAb-RIjbA> = %15.10lf\n", r2);
        outfile->Printf("<L|R>                = %15.10lf\n", dotval);
      }
      else {
        outfile->Printf("<L|R> zero by symmetry\n");
      }
      if (C_irr == L_irr ) {
       /* double check orthogonality rohf-like */
        global_dpd_->file2_init(&RIA, PSIF_CC_RAMPS, C_irr, 0, 1, R1A_lbl);
        global_dpd_->file2_init(&Ria, PSIF_CC_RAMPS, C_irr, 0, 1, R1B_lbl);
        global_dpd_->buf4_init(&RIJAB, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
        global_dpd_->buf4_init(&Rijab, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2BB_lbl);
        global_dpd_->buf4_init(&RIjAb, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);

        global_dpd_->file2_init(&LIA, PSIF_CC_OEI, L_irr, 0, 1, "LIA");
        global_dpd_->file2_init(&Lia, PSIF_CC_OEI, L_irr, 0, 1, "Lia");
        global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "LIJAB");
        global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "Lijab");
        global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");

        dot_IA = global_dpd_->file2_dot(&LIA, &RIA);
        dot_ia = global_dpd_->file2_dot(&Lia, &Ria);
        dot_IJAB = global_dpd_->buf4_dot(&LIJAB, &RIJAB);
        dot_ijab = global_dpd_->buf4_dot(&Lijab, &Rijab);
        dot_IjAb = global_dpd_->buf4_dot(&LIjAb, &RIjAb);

        global_dpd_->file2_close(&RIA);
        global_dpd_->file2_close(&Ria);
        global_dpd_->buf4_close(&RIJAB);
        global_dpd_->buf4_close(&Rijab);
        global_dpd_->buf4_close(&RIjAb);

        global_dpd_->file2_close(&LIA);
        global_dpd_->file2_close(&Lia);
        global_dpd_->buf4_close(&LIJAB);
        global_dpd_->buf4_close(&Lijab);
        global_dpd_->buf4_close(&LIjAb);

        outfile->Printf("\nROHF-like orthogonality test\n");
        outfile->Printf("<L0|R0>              = %15.10lf\n", eom_params.L0 * rzero);
        outfile->Printf("<LIA|RIA>            = %15.10lf\n", dot_IA);
        outfile->Printf("<Lia|Ria>            = %15.10lf\n", dot_ia);
        outfile->Printf("<LIJAB|RIJAB>        = %15.10lf\n", dot_IJAB);
        outfile->Printf("<Lijab|Rijab>        = %15.10lf\n", dot_ijab);
        outfile->Printf("<LIjAb|RIjAb>        = %15.10lf\n", dot_IjAb);
        outfile->Printf("<L|R>                = %15.10lf\n", eom_params.L0 * rzero +
        dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb);
      }
    } /* end dot with L, <L|R> overlap checks */
  } /* end loop over Cs */
  return;
}

}} // namespace psi::cceom
