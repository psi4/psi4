
/*
 *  primary_space.cc
 *  
 *
 *  Created by M.Saitow on 10/07/27.
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

void make_pp(void)
{
	int ai, bj, a, i, b, j, A, I, B, J;
	int Asym, Isym, Bsym, Jsym, h, nirreps;
	dpdbuf4 C, D, Amat;
	dpdfile2 fIJ, fAB;
      
	nirreps = moinfo.nirreps;
	if(params.spin == 0){
		dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
		dpd_buf4_sort(&D, CC_MISC, rpsq, 11, 11, "A(AI, BJ)");
		dpd_buf4_close(&D);
	}
	dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI, BJ)");
	dpd_buf4_scm(&Amat, 2);
	dpd_buf4_close(&Amat);
	dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_buf4_sort_axpy(&C, CC_MISC, qpsr, 11, 11, "A(AI, BJ)", -1);	
	dpd_buf4_close(&C);
      
	dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
	dpd_file2_mat_init(&fIJ);
	dpd_file2_mat_rd(&fIJ);
	dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
	dpd_file2_mat_init(&fAB);
	dpd_file2_mat_rd(&fAB);	
      
	dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI, BJ)");
	for(h = 0;h < nirreps;h++){
		dpd_buf4_mat_irrep_init(&Amat, h);
		dpd_buf4_mat_irrep_rd(&Amat, h);
		for(ai = 0;ai < Amat.params->rowtot[h];ai++){
			a = Amat.params->roworb[h][ai][0];
			i = Amat.params->roworb[h][ai][1];
			A = fAB.params->rowidx[a];
			I = fIJ.params->rowidx[i];
			Asym = fAB.params->psym[a];
			Isym = fIJ.params->psym[i];
			for(bj = 0;bj < Amat.params->coltot[h];bj++){
				b = Amat.params->colorb[h][bj][0];
				j = Amat.params->colorb[h][bj][1];
				B = fAB.params->colidx[b];
				J = fIJ.params->colidx[j];
				Bsym = fAB.params->qsym[b];
				Jsym = fIJ.params->qsym[j];
				if((A==B) && (Isym==Jsym)) Amat.matrix[h][ai][bj] -= fIJ.matrix[Isym][I][J];
				if((I==J) && (Asym==Bsym)) Amat.matrix[h][ai][bj] += fAB.matrix[Asym][A][B];
			}
		}
		dpd_buf4_mat_irrep_wrt(&Amat, h);
		dpd_buf4_mat_irrep_close(&Amat, h);
	}
	dpd_buf4_close(&Amat);
      
	dpd_file2_mat_close(&fAB);
	dpd_file2_close(&fAB);
	dpd_file2_mat_close(&fIJ);
	dpd_file2_close(&fIJ);
}
    
void make_pp_d(int iroot, int irrep)
{
	char lbl[32];
	dpdfile2 B, D, E, S;
	dpdbuf4 K, V;
      
	sprintf(lbl, "S^(%d)[%d]", iroot, irrep);
	dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
	sprintf(lbl, "B^(%d)[%d]", iroot, irrep);
	dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
      
	sprintf(lbl, "Dkc[%d]", irrep);
	dpd_file2_init(&D, CC_TMP0, irrep, 0, 1, lbl);
      
	if(params.spin == 0){
		dpd_buf4_init(&V, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
		dpd_dot13(&B, &V, &D, 0, 0, 1, 0);
	} else {
		dpd_buf4_init(&V, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
		dpd_dot13(&B, &V, &D, 0, 0, 1, 0);
	}
	dpd_buf4_close(&V);
      
	if(params.corr_type == "MP2")
		if(params.spin == 0)
			dpd_buf4_init(&K, CC_MISC,  0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
		else 
			dpd_buf4_init(&K, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
	else if(params.corr_type == "CC")
		if(params.spin == 0)
			dpd_buf4_init(&K, CC_TAMPS, 0, 0, 5, 0, 5, 0,     "2 tIjAb - tIjbA");
		else 
			dpd_buf4_init(&K, CC_TAMPS, 0, 0, 5, 0, 5, 0,     "tIjAb");
	else if (params.corr_type == "PR")
		if(params.spin == 0)
			dpd_buf4_init(&K, CC_MISC, 0, 0, 5, 0, 5, 0, "tilde 2 tIjAb - tIjbA");
		else
			dpd_buf4_init(&K, CC_MISC, 0, 0, 5, 0, 5, 0, "tilde tIjAb");
      
	dpd_dot24(&D, &K, &S, 0, 0,  0.5, 1);		
	dpd_buf4_close(&K);
	dpd_file2_close(&D);
      
	sprintf(lbl, "Ekc[%d]", irrep);
	dpd_file2_init(&E, CC_TMP0, irrep, 0, 1, lbl);
      
	if(params.corr_type == "MP2")
		if(params.spin == 0)
			dpd_buf4_init(&K, CC_MISC,  0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
		else
			dpd_buf4_init(&K, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
	else if (params.corr_type == "CC")
		if(params.spin == 0)
			dpd_buf4_init(&K, CC_TAMPS, 0, 0, 5, 0, 5, 0,     "2 tIjAb - tIjbA");
		else
			dpd_buf4_init(&K, CC_TAMPS, 0, 0, 5, 0, 5, 0,     "tIjAb");
	else if(params.corr_type == "PR")
		if(params.spin == 0)
			dpd_buf4_init(&K, CC_MISC,  0, 0, 5, 0, 5, 0, "tilde 2 tIjAb - tIjbA");
		else
			dpd_buf4_init(&K, CC_MISC, 0, 0, 5, 0, 5, 0, "tilde tIjAb");
	
	dpd_dot13(&B, &K, &E, 0, 0, 1, 0);
	dpd_buf4_close(&K);
      
	if(params.spin == 0)
		dpd_buf4_init(&V, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
	else
		dpd_buf4_init(&V, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      
	dpd_dot24(&E, &V, &S, 0, 0, 0.5, 1);
	dpd_buf4_close(&V);
	dpd_file2_close(&E);
      
	dpd_file2_close(&B);
	dpd_file2_close(&S);
}
    
}}
