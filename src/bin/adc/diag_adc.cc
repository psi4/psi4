
/*
 *  diag_adc.cc
 *  
 *
 *  Created by M.Saitow on 10/07/27.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {

void transpert(void);
void sort_pert(void);
double calc_ocss(char *, int, int, int, double);
void amp_write_T1(dpdfile2 *, int, FILE *);
    
void diag_adc(void)
{
	char lbl[32], *axis;
	int i, irrep, dim, init_dim;
	int root, nroot;
	int ai, bj, ck, c, k, C, K, Ksym;
	double value, *omega, **V;
	dpdfile2 B;
	dpdbuf4 A;
	FILE *iter_adc;
      
//	transpert();
//	sort_pert();
      
	dpd_buf4_init(&A, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI, BJ)");
	for(irrep = 0;irrep < moinfo.nirreps;irrep++){
		init_dim = params.init_dim[irrep];
		dim = A.params->rowtot[irrep];
		guess[irrep] = init_array(init_dim);
		params.ocss[irrep]     = init_array(params.ppi[irrep]);
		params.d_ocss[irrep]   = init_array(params.ppi[irrep]);
		params.mu_sqs[irrep]   = (double **)malloc(params.ppi[irrep]*sizeof(double *));
		params.d_mu_sqs[irrep] = (double **)malloc(params.ppi[irrep]*sizeof(double *));
		
		if(params.diag_method == "DAVIDSON"){				
			omega = init_array(init_dim);
			V = block_matrix(dim, init_dim);
		}
		else {
			omega = init_array(dim);
			V = block_matrix(dim, dim);
		}
		dpd_buf4_mat_irrep_init(&A, irrep);
		dpd_buf4_mat_irrep_rd(&A, irrep);
		if(params.diag_method == "DAVIDSON"){
			david(A.matrix[irrep], dim, init_dim, omega, V, 1e-14, 0);
		}
		else{
			sq_rsp(dim, dim, A.matrix[irrep], omega, 1, V, 1e-14);
		}
		for(root = 0;root < init_dim;root++){
	  
			for(i = 0;i < 3;i++){
				params.mu_sqs[irrep][root]   = init_array(3);
				params.d_mu_sqs[irrep][root] = init_array(3);
			}
	  
			sprintf(lbl, "BIA(%d)[%d] singlet", root, irrep);
			dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
			dpd_file2_mat_init(&B);
			for(ck = 0;ck < dim;ck++){
				c = A.params->roworb[irrep][ck][0];
				k = A.params->roworb[irrep][ck][1];
		
				K = B.params->rowidx[k];
				C = B.params->colidx[c];
				Ksym = B.params->psym[k];
				B.matrix[Ksym][K][C] = V[ck][root];
			}
			dpd_file2_mat_wrt(&B);
			dpd_file2_mat_close(&B);
			dpd_file2_close(&B);
	  
//			params.ocss[irrep][root] = calc_ocss(lbl, irrep, 0, root, omega[root]);
		}
		C_DCOPY(init_dim, omega, 1, guess[irrep], 1);
		free_block(V);
		
		ffile(&iter_adc, "iter.dat", 1);
		fprintf(iter_adc, "  irrep %s, dim=%d:\n", moinfo.labels[irrep], dim);
		for(i = 0;i < params.ppi[irrep];i++) 
			fprintf(iter_adc, "\t%8.5lf   %8.5lf\n", omega[i], omega[i]*_hartree2ev);
		fprintf(iter_adc, "\n");
		fclose(iter_adc);
			
		free(omega);
	}
	dpd_buf4_close(&A);
      
	axis = (char *)malloc(3*sizeof(char));
	axis[0] = 'X'; axis[1] = 'Y'; axis[2] = 'Z';
		
	fprintf(outfile, "\tRHF-ADC(1)/CIS");
	if(params.spin == 0) fprintf(outfile, " Singlet");
	else fprintf(outfile, " Triplet");
	fprintf(outfile, " Excitation Energies:\n");
	fprintf(outfile, "\t------------------------------------\n\n");
	fprintf(outfile, "\t  Root Irrep       Hartree          eV          cm-1  \n");
	fprintf(outfile, "\t  ---- ----- ------------------  ---------  ----------\n");
	for(irrep = 0;irrep < moinfo.nirreps;irrep++) {
		for(root = 0;root < params.ppi[irrep];root++) {
			value = guess[irrep][root];
			fprintf(outfile,   "\tState %4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[irrep],
							value, value*_hartree2ev, value*_hartree2wavenumbers);
			fprintf(outfile, "\n\tLargest components of singlet excited wave function #%d/#%d:\n",
							irrep, root);
			sprintf(lbl, "BIA(%d)[%d] singlet", root, irrep);
			dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
			amp_write_T1(&B, params.num_amps, outfile);
			dpd_file2_close(&B);
			fprintf(outfile, "\n");
//			fprintf(outfile, "\tComponents of squared dipole moment:\n");
//			fprintf(outfile, "\t");
//			for(i = 0;i < 3;i++)
//				fprintf(outfile, "%c: %10.7f ", axis[i], params.mu_sqs[irrep][root][i]);
//			fprintf(outfile, "\n\tOscillator strength is %10.7f.\n\n", params.ocss[irrep][root]);
			fflush(outfile);
		}
	}
		
}

}}
