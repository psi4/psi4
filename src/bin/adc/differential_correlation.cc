
/*
 *  differential_correlation.cc
 *  
 *
 *  Created by M.Saitow on 11/07/14.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {
	
void renormalize(void)
{
	char lbl[32];
	int ai, bj, a, i, b, j, A, I, B, J;
	int irrep, Asym, Bsym, Isym, Jsym;
	dpdbuf4 V, K, Amat;
	dpdfile2 Aab, Aij, Xab, Xij;
      
	dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI, BJ)");
	dpd_buf4_copy(&Amat, CC_MISC, "dA(AI, BJ)");
	dpd_buf4_close(&Amat);
      
	if(params.corr_type == "MP2")
		dpd_buf4_init(&K, CC_MISC,  0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
	else if(params.corr_type == "CC")
		dpd_buf4_init(&K, CC_TAMPS, 0, 0, 5, 0, 5, 0,     "2 tIjAb - tIjbA");
	else if(params.corr_type == "PR")
		dpd_buf4_init(&K, CC_MISC,  0, 0, 5, 0, 5, 0, "tilde 2 tIjAb - tIjbA");

	dpd_buf4_init(&V, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      
	dpd_file2_init(&Xab, CC_TMP1, 0, 1, 1, "Xab");
	dpd_contract442(&K, &V, &Xab, 2, 2, 1, 0);
	dpd_file2_init(&Aab, CC_TMP1, 0, 1, 1, "Aab");
	dpd_file2_axpy(&Xab, &Aab, -0.5, 0);
	dpd_file2_axpy(&Xab, &Aab, -0.5, 1);
	dpd_file2_close(&Xab);
      
	dpd_file2_init(&Xij, CC_TMP0, 0, 0, 0, "Xij");
	dpd_contract442(&K, &V, &Xij, 0, 0, 1, 0);
	dpd_file2_init(&Aij, CC_TMP0, 0, 0, 0, "Aij");
	dpd_file2_axpy(&Xij, &Aij, -0.5, 0);
	dpd_file2_axpy(&Xij, &Aij, -0.5, 1);
	dpd_file2_close(&Xij);
      
	dpd_buf4_close(&V);
	dpd_buf4_close(&K);
      
	dpd_file2_mat_init(&Aab);
	dpd_file2_mat_rd(&Aab);
	dpd_file2_mat_init(&Aij);
	dpd_file2_mat_rd(&Aij);
      
	dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "dA(AI, BJ)");
	for(irrep = 0;irrep < moinfo.nirreps;irrep++){
		dpd_buf4_mat_irrep_init(&Amat, irrep);
		dpd_buf4_mat_irrep_rd(&Amat, irrep);
		for(ai = 0;ai < Amat.params->rowtot[irrep];ai++){
			a = Amat.params->roworb[irrep][ai][0];
			i = Amat.params->roworb[irrep][ai][1];
			A = Aab.params->rowidx[a];
			I = Aij.params->rowidx[i];
			Asym = Aab.params->psym[a];
			Isym = Aij.params->psym[i];
			for(bj = 0;bj < Amat.params->coltot[irrep];bj++){
				b = Amat.params->colorb[irrep][bj][0];
				j = Amat.params->colorb[irrep][bj][1];
				B = Aab.params->colidx[b];
				J = Aij.params->colidx[j];
				Bsym = Aab.params->qsym[b];
				Jsym = Aij.params->qsym[j];
				if((A==B) && (Isym==Jsym)) 
					Amat.matrix[irrep][ai][bj] += Aij.matrix[Isym][I][J];
				if((I==J) && (Asym==Bsym)) 
					Amat.matrix[irrep][ai][bj] += Aab.matrix[Asym][A][B];
			}
		}
		dpd_buf4_mat_irrep_wrt(&Amat, irrep);
		dpd_buf4_mat_irrep_close(&Amat, irrep);
	}
	dpd_buf4_close(&Amat);
      
	dpd_file2_mat_close(&Aab);
	dpd_file2_close(&Aab);
	dpd_file2_mat_close(&Aij);
	dpd_file2_close(&Aij);
}
    
}}
