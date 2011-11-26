/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <string>
#include <libpsio/psio.h>
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
      if(psio_tocscan(CC_INFO, E_lbl) == NULL) {
        fprintf(outfile,"No EOM CC2 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      sprintf(E_lbl, "EOM CCSD Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(CC_INFO, E_lbl) == NULL) {
        fprintf(outfile,"No EOM CCSD Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      sprintf(E_lbl, "EOM CC3 Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(CC_INFO, E_lbl) == NULL) {
        fprintf(outfile,"No EOM CC3 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }

    sprintf(R1A_lbl, "RIA %d %d", C_irr, R_index);
    sprintf(R1B_lbl, "Ria %d %d", C_irr, R_index);
    sprintf(R2AA_lbl, "RIJAB %d %d", C_irr, R_index);
    sprintf(R2BB_lbl, "Rijab %d %d", C_irr, R_index);
    sprintf(R2AB_lbl, "RIjAb %d %d", C_irr, R_index);

    /* Calculate <0| Hbar R |0> */
    if (C_irr == H_IRR) {
      dpd_file2_init(&FIA, CC_OEI, H_IRR, A_OCC, A_VIR, "FME");
      dpd_file2_init(&RIA, CC_RAMPS, C_irr, A_OCC, A_VIR, R1A_lbl);
      dot_IA = dpd_file2_dot(&FIA, &RIA);
      dpd_file2_close(&RIA);
      dpd_file2_close(&FIA);
  
      dpd_file2_init(&Fia, CC_OEI, H_IRR, B_OCC, B_VIR, "Fme");
      dpd_file2_init(&Ria, CC_RAMPS, C_irr, B_OCC, B_VIR, R1B_lbl);
      dot_ia = dpd_file2_dot(&Fia, &Ria);
      dpd_file2_close(&Ria);
      dpd_file2_close(&Fia);

      if (params.eom_ref == 1) {
        dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
        dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
        dot_IJAB = dpd_buf4_dot(&D, &RIJAB);
        dpd_buf4_close(&RIJAB);
        dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2BB_lbl);
        dot_ijab = dpd_buf4_dot(&D, &Rijab);
        dpd_buf4_close(&Rijab);
        dpd_buf4_close(&D);
    
        dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
        dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
        dot_IjAb = dpd_buf4_dot(&D, &RIjAb);
        dpd_buf4_close(&RIjAb);
        dpd_buf4_close(&D);
      }
      else if (params.eom_ref == 2) {
        dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
        dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
        dot_IJAB = dpd_buf4_dot(&D, &RIJAB);
        dpd_buf4_close(&RIJAB);
        dpd_buf4_close(&D);
  
        dpd_buf4_init(&D, CC_DINTS, H_IRR, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
        dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 12, 17, 12, 17, 0, R2BB_lbl);
        dot_ijab = dpd_buf4_dot(&D, &Rijab);
        dpd_buf4_close(&Rijab);
        dpd_buf4_close(&D);
                           
        dpd_buf4_init(&D, CC_DINTS, H_IRR, 22, 28, 22, 28, 0, "D <Ij|Ab>");
        dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 22, 28, 22, 28, 0, R2AB_lbl);
        dot_IjAb = dpd_buf4_dot(&D, &RIjAb);
        dpd_buf4_close(&RIjAb);
        dpd_buf4_close(&D);
      }
      rzero = (dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb)/energy;
    }
    else { /* C and H are different irreps */
      rzero = 0.0;
    }

    /* Now normalize so that <R|R> = 1 */
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, A_OCC, A_VIR, R1A_lbl);
    dpd_file2_init(&Ria, CC_RAMPS, C_irr, B_OCC, B_VIR, R1B_lbl);
    dpd_buf4_init(&fRIJAB, CC_RAMPS, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, R2AA_lbl);
    dpd_buf4_init(&fRijab, CC_RAMPS, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, R2BB_lbl);
    dpd_buf4_init(&fRIjAb, CC_RAMPS, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, R2AB_lbl);

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
    fprintf(outfile,"<R|R> = %20.16lf\n",norm);

/* just debugging with converged solutions - also my need a sort_C() */
/*
dpd_file2_copy(&RIA, EOM_CME, "CME 0");
dpd_file2_copy(&Ria, EOM_Cme, "Cme 0");
dpd_buf4_copy(&fRIJAB, EOM_CMNEF, "CMNEF 0");
dpd_buf4_copy(&fRijab, EOM_Cmnef, "Cmnef 0");
dpd_buf4_copy(&fRIjAb, EOM_CMnEf, "CMnEf 0");
*/
/* end debugging stuff */ 


    dpd_file2_close(&RIA);
    dpd_file2_close(&Ria);
    dpd_buf4_close(&fRIJAB);
    dpd_buf4_close(&fRijab);
    dpd_buf4_close(&fRIjAb);

    if(params.wfn == "EOM_CC2") {
      fprintf(outfile,"EOM CC2 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC2 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      fprintf(outfile,"EOM CCSD R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CCSD R0 for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      fprintf(outfile,"EOM CC3 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC3 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }

    if (eom_params.dot_with_L) {
      /* evaluate check <R|L> == 0 */
      if (C_irr == L_irr ) {
        dpd_file2_init(&RIA, CC_RAMPS, C_irr, A_OCC, A_VIR, R1A_lbl);
        dpd_file2_init(&Ria, CC_RAMPS, C_irr, B_OCC, B_VIR, R1B_lbl);
        dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, R2AA_lbl);
        dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, R2BB_lbl);
        dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, R2AB_lbl);
  
        dpd_file2_init(&LIA, CC_OEI, L_irr, A_OCC, A_VIR, "LIA");
        dpd_file2_init(&Lia, CC_OEI, L_irr, B_OCC, B_VIR, "Lia");
        dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, "LIJAB");
        dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, "Lijab");
        dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, "LIjAb");
  
        fprintf(outfile,"\nROHF orthogonality test\n");
        fprintf(outfile,"<L0|R0>            = %15.10lf\n", eom_params.L0 * rzero);
        dot_IA = dpd_file2_dot(&LIA, &RIA);
        fprintf(outfile,"<LIA|RIA>          = %15.10lf\n", dot_IA);
        dot_ia = dpd_file2_dot(&Lia, &Ria);
        fprintf(outfile,"<Lia|Ria>          = %15.10lf\n", dot_ia);
        dot_IJAB = dpd_buf4_dot(&LIJAB, &RIJAB);
        fprintf(outfile,"<LIJAB|RIJAB>      = %15.10lf\n", dot_IJAB);
        dot_ijab = dpd_buf4_dot(&Lijab, &Rijab);
        fprintf(outfile,"<Lijab|Rijab>      = %15.10lf\n", dot_ijab);
        dot_IjAb = dpd_buf4_dot(&LIjAb, &RIjAb);
        fprintf(outfile,"<LIjAb|RIjAb>      = %15.10lf\n", dot_IjAb);
        fprintf(outfile,"<L|R>              = %15.10lf\n", (eom_params.L0 * rzero)
            + dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb);
  
        dpd_file2_close(&LIA);
        dpd_file2_close(&Lia);
        dpd_buf4_close(&LIJAB);
        dpd_buf4_close(&Lijab);
        dpd_buf4_close(&LIjAb);
        dpd_file2_close(&RIA);
        dpd_file2_close(&Ria);
        dpd_buf4_close(&RIJAB);
        dpd_buf4_close(&Rijab);
        dpd_buf4_close(&RIjAb);
      }
      else {
        fprintf(outfile,"\nOverlap <R|L> zero by symmetry\n");
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
      if(psio_tocscan(CC_INFO, E_lbl) == NULL) { 
        fprintf(outfile,"No EOM CC2 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      sprintf(E_lbl, "EOM CCSD Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(CC_INFO, E_lbl) == NULL) { 
        fprintf(outfile,"No EOM CCSD Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      sprintf(E_lbl, "EOM CC3 Energy for root %d %d", C_irr, R_index);
      if(psio_tocscan(CC_INFO, E_lbl) == NULL) { 
        fprintf(outfile,"No EOM CC3 Energy found in CC_INFO.  Not normalizing R.\n");
        return;
      }
      psio_read_entry(CC_INFO, E_lbl, (char *) &(energy), sizeof(double));
    }

    sprintf(R1A_lbl, "RIA %d %d", C_irr, R_index);
    sprintf(R1B_lbl, "Ria %d %d", C_irr, R_index);
    sprintf(R2AB_lbl, "RIjAb %d %d", C_irr, R_index);
    sprintf(R2AA_lbl, "RIJAB %d %d", C_irr, R_index);
    sprintf(R2BB_lbl, "Rijab %d %d", C_irr, R_index);

    /* produce RIjbA and 2RIjAb-RIjbA copies - not yet normalized */
    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    dpd_buf4_sort(&RIjAb, CC_TMP, pqsr, 0, 5, "RIjbA");
    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    dpd_buf4_copy(&RIjAb, CC_RAMPS, lbl);
    dpd_buf4_close(&RIjAb);

    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&RIjAb, 2.0);
    dpd_buf4_init(&RIjbA, CC_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");
    dpd_buf4_axpy(&RIjbA, &RIjAb, -1.0);
    dpd_buf4_close(&RIjbA);
    dpd_buf4_close(&RIjAb);

    /* calculate <0| hbar | 0> */
    if (!params.full_matrix) {
      if (C_irr == H_IRR) {
        dpd_file2_init(&FIA, CC_OEI, H_IRR, 0, 1, "FME");
        dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
        r1 = 2.0 * dpd_file2_dot(&FIA, &RIA);
        dpd_file2_close(&RIA);
        dpd_file2_close(&FIA);
    
        sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
        dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
        dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
        r2 = dpd_buf4_dot(&D, &RIjAb);
        dpd_buf4_close(&D);
        dpd_buf4_close(&RIjAb);
        rzero = (r1 + r2)/energy;
      }
      else {
        rzero = 0.0;
      }
    }
    else { /* full matrix */
      sprintf(lbl, "%s %d %d", "R0", C_irr, R_index);
      psio_read_entry(CC_RAMPS, lbl, (char *) &rzero, sizeof(double));
    }

    /* Now normalize R so that <R|R> = 1 */
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    dpd_buf4_init(&RIjbA, CC_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");

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
    dpd_file2_scm(&RIA, 1.0/norm);
    dpd_buf4_scm(&RIjAb, 1.0/norm);
    dpd_buf4_scm(&RIjbA, 1.0/norm);

/* just for debugging cc3 put normalized vector back into C as well */
/*
dpd_file2_copy(&RIA, EOM_CME, "CME 0");
dpd_buf4_copy(&RIjAb, EOM_CMnEf, "CMnEf 0");
*/
/* end debugging stuff */ 

    dpd_file2_close(&RIA);
    dpd_buf4_close(&RIjAb);
    dpd_buf4_close(&RIjbA);

    if(params.wfn == "EOM_CC2") {
      fprintf(outfile,"EOM CC2 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC2 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      fprintf(outfile,"EOM CCSD R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CCSD R0 for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      fprintf(outfile,"EOM CC3 R0 for root %d = %15.11lf\n", R_index, rzero);
      sprintf(lbl, "EOM CC3 R0 for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, lbl, (char *) &rzero, sizeof(double));
    }

    /* produce ROHF like quantities and 2RIjAb-RIjbA */
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    dpd_file2_copy(&RIA, CC_RAMPS, R1B_lbl);
    dpd_file2_close(&RIA);

    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    dpd_buf4_copy(&RIjAb, CC_RAMPS, lbl);
    dpd_buf4_close(&RIjAb);

    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 2, 7, 0, 5, 1, R2AB_lbl);
    dpd_buf4_copy(&RIjAb, CC_RAMPS, R2AA_lbl);
    dpd_buf4_copy(&RIjAb, CC_RAMPS, R2BB_lbl);
    dpd_buf4_close(&RIjAb);
      
    dpd_buf4_init(&RIjbA, CC_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");
    sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&RIjAb, 2.0);
    dpd_buf4_axpy(&RIjbA, &RIjAb, -1.0);
    dpd_buf4_close(&RIjAb);
    dpd_buf4_close(&RIjbA);
  
    /* test normalization with produced ROHF-like quantities */

    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    norm = dpd_file2_dot_self(&RIA);
    dpd_file2_close(&RIA);
    dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, R1B_lbl);
    norm += dpd_file2_dot_self(&Ria);
    dpd_file2_close(&Ria);
    dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += dpd_buf4_dot_self(&RIJAB);
    dpd_buf4_close(&RIJAB);
    dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2BB_lbl);
    norm += dpd_buf4_dot_self(&Rijab);
    dpd_buf4_close(&Rijab);
    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    norm += dpd_buf4_dot_self(&RIjAb);
    dpd_buf4_close(&RIjAb);
    norm += rzero * rzero;

#ifdef EOM_DEBUG
    fprintf(outfile,"Norm with produced ROHF-like quantities = %15.10lf\n", norm);
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
    dpd_file2_print(&RIA, outfile);
    dpd_file2_close(&RIA);
#endif

    /* check orthogonality with L */
    if (eom_params.dot_with_L) {
      if (C_irr == L_irr) {
        dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
        dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
        r1 = 2.0 * dpd_file2_dot(&LIA, &RIA);
        dpd_file2_close(&RIA);
        dpd_file2_close(&LIA);
  
        dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
        sprintf(lbl, "%s %d %d", "2RIjAb - RIjbA", C_irr, R_index);
        dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
        r2 = dpd_buf4_dot(&LIjAb, &RIjAb);
        dpd_buf4_close(&RIjAb);
        dpd_buf4_close(&LIjAb);
  
        dotval = r1 + r2 + (eom_params.L0 * rzero);
        fprintf(outfile,"Performing RHF orthogonality test\n");
        fprintf(outfile,"<L0|R0>              = %15.10lf\n", eom_params.L0 * rzero);
        fprintf(outfile,"2*<LIA|RIA>          = %15.10lf\n", r1);
        fprintf(outfile,"<LIjAb|2RIjAb-RIjbA> = %15.10lf\n", r2);
        fprintf(outfile,"<L|R>                = %15.10lf\n", dotval);
      }
      else {
        fprintf(outfile,"<L|R> zero by symmetry\n");
      }
      if (C_irr == L_irr ) {
       /* double check orthogonality rohf-like */
        dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, R1A_lbl);
        dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, R1B_lbl);
        dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2AA_lbl);
        dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, R2BB_lbl);
        dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, R2AB_lbl);
  
        dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
        dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");
        dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "LIJAB");
        dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "Lijab");
        dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  
        dot_IA = dpd_file2_dot(&LIA, &RIA);
        dot_ia = dpd_file2_dot(&Lia, &Ria);
        dot_IJAB = dpd_buf4_dot(&LIJAB, &RIJAB);
        dot_ijab = dpd_buf4_dot(&Lijab, &Rijab);
        dot_IjAb = dpd_buf4_dot(&LIjAb, &RIjAb);
  
        dpd_file2_close(&RIA);
        dpd_file2_close(&Ria);
        dpd_buf4_close(&RIJAB);
        dpd_buf4_close(&Rijab);
        dpd_buf4_close(&RIjAb);
  
        dpd_file2_close(&LIA);
        dpd_file2_close(&Lia);
        dpd_buf4_close(&LIJAB);
        dpd_buf4_close(&Lijab);
        dpd_buf4_close(&LIjAb);

        fprintf(outfile,"\nROHF-like orthogonality test\n");
        fprintf(outfile,"<L0|R0>              = %15.10lf\n", eom_params.L0 * rzero);
        fprintf(outfile,"<LIA|RIA>            = %15.10lf\n", dot_IA);
        fprintf(outfile,"<Lia|Ria>            = %15.10lf\n", dot_ia);
        fprintf(outfile,"<LIJAB|RIJAB>        = %15.10lf\n", dot_IJAB);
        fprintf(outfile,"<Lijab|Rijab>        = %15.10lf\n", dot_ijab);
        fprintf(outfile,"<LIjAb|RIjAb>        = %15.10lf\n", dot_IjAb);
        fprintf(outfile,"<L|R>                = %15.10lf\n", eom_params.L0 * rzero +
        dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb);
      }
    } /* end dot with L, <L|R> overlap checks */
  } /* end loop over Cs */
  return;
}

}} // namespace psi::cceom
