
/*
 *  mp2.cc
 *  
 *
 *  Created by M.Saitow on 10/07/31.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {

void mp2(void)
{
  int ij, ab;
  int i, j, I, J, isym, jsym;
  int iter, h, nirreps, row, col;
  double energy, conv, rms, value, sq_norm, Nij;
  dpdfile2 A, F;
  dpdbuf4 D, T2, newT2, Z;

  nirreps = moinfo.nirreps;

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_MISC, "MP2 tIjAb");
  dpd_buf4_close(&D);
  
  dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
  
  dpd_buf4_init(&D, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&D, &T2);
  dpd_buf4_close(&D);
  
  dpd_buf4_copy(&T2, CC_MISC, "New MP2 tIjAb");
  dpd_buf4_close(&T2);
  
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
  energy = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  
  fprintf(outfile, "\n\tSolving for MP2 wave function:\n");
  fprintf(outfile,   "\t-------------------------------\n");
  fprintf(outfile, "\titer = %d  MP2 Energy = %20.14f\n", 0, energy);
  
  conv = 0;
  for(iter=1; iter < params.maxiter; iter++) {
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_MISC, "New MP2 tIjAb Increment");
    dpd_buf4_close(&D);
    
    dpd_buf4_init(&newT2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb Increment");
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
    dpd_contract424(&T2, &F, &newT2, 1, 0, 1, -1, 1);
    dpd_contract244(&F, &T2, &newT2, 0, 0, 0, -1, 1);
    dpd_file2_close(&F);
    
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    dpd_contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
    dpd_contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
    dpd_file2_close(&F);
    
    dpd_buf4_close(&T2);
    
    dpd_buf4_init(&D, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    dpd_buf4_dirprd(&D, &newT2);
    dpd_buf4_close(&D);
    
    dpd_buf4_close(&newT2);
    
    dpd_buf4_init(&newT2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb Increment");
    dpd_buf4_axpy(&T2, &newT2, 1);
    dpd_buf4_close(&T2);
    
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    energy = dpd_buf4_dot(&D, &newT2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&newT2);
    
    dpd_buf4_init(&newT2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    rms = 0.0;
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&newT2, h);
      dpd_buf4_mat_irrep_rd(&newT2, h);
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++) {
	  value = newT2.matrix[h][row][col] - T2.matrix[h][row][col];
	  rms += value * value;
	}
      
      dpd_buf4_mat_irrep_close(&T2, h);
      dpd_buf4_mat_irrep_close(&newT2, h);
    }
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);
    rms = sqrt(rms);
    
    fprintf(outfile, "\titer = %d  MP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms);
    
    if(rms < params.criterion) {
      conv = 1;
      fprintf(outfile, "\n\tMP2 iterations converged.\n");
      break;
    }
    else {
      dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
      dpd_buf4_copy(&T2, CC_MISC, "MP2 tIjAb");
      dpd_buf4_close(&T2);
    }
  }
  
  if(!conv) {
    fprintf(outfile, "\n\tMP2 iterative procedure failed.\n");
    exit(PSI_RETURN_FAILURE);
  }
  
  dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
  dpd_buf4_sort(&T2, CC_TMP0, pqsr, 0, 5, "MP2 tIjbA");
  dpd_buf4_copy(&T2, CC_MISC, "MP2 2 tIjAb - tIjbA");
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
  dpd_buf4_scm(&T2, 2);
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "MP2 tIjbA");
  dpd_buf4_axpy(&Z, &T2, -1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
  sq_norm = 1 + dpd_buf4_dot(&Z, &T2);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&Z);
  fprintf(outfile, "\t[Squared norm of MP2 wavefunction = %10.7f]\n", sq_norm);
  fflush(outfile);
  
  if(params.corr_type == "PR"){
    /* 
     *  This bunch of codes is no so efficient because we can accomplish the partial
     *  renormalization of pair-correlation function at cost of merely O(N^4). So, 
     *  will be modified in near future !
     */
    dpd_file2_init(&A, CC_MISC, 0, 0, 0, "Aij");
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    dpd_buf4_init(&D,  CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
    dpd_contract442(&D, &T2, &A, 0, 0, 1, 0);
    dpd_buf4_close(&D);
    
    dpd_file2_mat_init(&A);
    dpd_file2_mat_rd(&A);
    dpd_buf4_init(&newT2, CC_MISC, 0, 0, 5, 0, 5, 0, "tilde tIjAb");
    
    for(h = 0;h < moinfo.nirreps;h++){
      dpd_buf4_mat_irrep_init(&T2, h);
      dpd_buf4_mat_irrep_init(&newT2, h);
      dpd_buf4_mat_irrep_rd(&T2, h);
      
      for(ij = 0;ij < T2.params->rowtot[h];ij++){
	i = T2.params->roworb[h][ij][0];
	j = T2.params->roworb[h][ij][1];
	isym = T2.params->psym[i];
	jsym = T2.params->qsym[j];
	
	I = i - moinfo.occ_off[isym];
	J = j - moinfo.occ_off[jsym];
	
	Nij = 1 + (A.matrix[isym][I][I] + A.matrix[jsym][J][J]) / 2;
	for(ab = 0;ab < T2.params->coltot[h];ab++)
	  newT2.matrix[h][ij][ab] = T2.matrix[h][ij][ab] / Nij;
      }
      dpd_buf4_mat_irrep_wrt(&newT2, h);
      dpd_buf4_mat_irrep_close(&newT2, h);
      dpd_buf4_mat_irrep_close(&T2, h);
    }
    dpd_file2_mat_close(&A);
    dpd_file2_close(&A);
    dpd_buf4_close(&T2);
    
    dpd_buf4_sort(&newT2, CC_MISC, pqsr, 0, 5, "tilde tIjbA");
    dpd_buf4_scmcopy(&newT2, CC_MISC, "tilde 2 tIjAb - tIjbA", 2.0);
    dpd_buf4_close(&newT2);
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "tilde 2 tIjAb - tIjbA");
    dpd_buf4_init(&Z,  CC_MISC, 0, 0, 5, 0, 5, 0, "tilde tIjbA");
    dpd_buf4_axpy(&Z, &T2, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    energy = dpd_buf4_dot(&T2, &Z);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    sq_norm = 1 + dpd_buf4_dot(&Z, &T2);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    fprintf(outfile, "\n");
    fprintf(outfile, "\tPR-MP2 energy = %20.14f\n", energy);
    fprintf(outfile, "\t[Squared norm of PR-MP2 wavefunction = %10.7f]\n", sq_norm);
    fflush(outfile);
  }
  
}

}}
