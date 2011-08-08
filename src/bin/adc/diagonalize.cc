
/*
 *  diagonalize.cc
 *  
 *
 *  Created by M.Saitow on 10/07/31.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {namespace adc {
	
#define MAXIT 100
#define NORM_TOL 1e-5
    
    void denom(int, double);
    void make_pp_d(int, int);
    void make_screened_interaction(int, int, double);
    void shift_diag(double, int, int);
    
    double diagonalize(int irrep, int root, double omega)
    {
      char lbl[32];
      int I, J, k;
      int nexc, nirreps, nroots, converged, skip_check, length, init_dim;
      int *conv, dim, maxdim, iter, iroot;
      double **Alpha, **G, *eps, *lambda_o, *lambda, cutoff, sum, coeff;
      double norm, diff;
      dpdfile2 B, Bp, B_n, F, L, S, V;
      dpdbuf4 A;
      FILE *iter_adc;
      
      init_dim = params.init_dim[irrep];
      nirreps  = moinfo.nirreps;
      nexc     = moinfo.nexc[irrep];
      nroots   = params.ppi[irrep];
      maxdim   = (params.dim+params.width) * nroots;
      cutoff   = params.criterion;
      
      G        = block_matrix(maxdim, maxdim);
      Alpha    = block_matrix(maxdim, maxdim);
      eps      = init_array(nroots);
      lambda   = init_array(maxdim);
      lambda_o = init_array(maxdim);
      
      for(iroot = 0;iroot < init_dim;iroot++){
	sprintf(lbl, "BIA(%d)[%d] singlet", iroot, irrep);
	dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	sprintf(lbl, "B^(%d)[%d]", iroot, irrep);
	dpd_file2_copy(&B, CC_OEI, lbl);
	dpd_file2_close(&B);
      }
      C_DCOPY(init_dim, guess[irrep], 1, lambda_o, 1);
      
      length = init_dim;
      iter = 0;
      converged = 0;
      conv = init_int_array(nroots);
      
      denom(irrep, omega);
      
      ffile(&iter_adc, "iter.dat", 1);
      
      timer_on("Total_SEM");
      while (converged < nroots && iter < MAXIT){
	skip_check = 0;
	fprintf(iter_adc, "\niter = %d, dim = %d\n", iter, length);
	
	timer_on("Sigma");
	for(iroot = 0;iroot < length;iroot++){
	  sprintf(lbl, "B^(%d)[%d]", iroot, irrep);
	  dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	  sprintf(lbl, "S^(%d)[%d]", iroot, irrep);
	  dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
	  dpd_buf4_init(&A, CC_MISC, 0, 11, 11, 11, 11, 0, "dA(AI, BJ)");
	  dpd_contract422(&A, &B, &S, 1, 1, 1, 0);
	  dpd_buf4_close(&A);
	  dpd_file2_close(&B);
	  dpd_file2_close(&S);
	  
	  make_pp_d(iroot, irrep);
	  make_screened_interaction(iroot, irrep, omega);
	}
	for(I = 0;I < length;I++){
	  sprintf(lbl, "S^(%d)[%d]", I, irrep);
	  dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
	  for(J = 0;J <= I;J++){
	    sprintf(lbl, "B^(%d)[%d]", J, irrep);
	    dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	    sum = dpd_file2_dot(&B, &S);
	    if(I != J)
	      G[I][J] = G[J][I] = sum;
	    else
	      G[I][J] = sum;
	    dpd_file2_close(&B);
	  }
	  dpd_file2_close(&S);
	}
	timer_off("Sigma");
	sq_rsp(length, length, G, lambda, 1, Alpha, 1e-12);
	timer_on("SEM");
	for(k = 0;k < nroots;k++){
	  sprintf(lbl, "F^(%d)[%d]", k, irrep);
	  dpd_file2_init(&F, CC_OEI, irrep, 0, 1, lbl);
	  dpd_file2_scm(&F, 0.0);
	  for(I = 0;I < length;I++){
	    sprintf(lbl, "B^(%d)[%d]", I, irrep);
	    dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	    dpd_file2_axpy(&B, &F, -Alpha[I][k]*lambda[k], 0);
	    dpd_file2_close(&B);
	    sprintf(lbl, "S^(%d)[%d]", I, irrep);
	    dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
	    dpd_file2_axpy(&S, &F, Alpha[I][k], 0);
	    dpd_file2_close(&S);
	  }
	  shift_diag(lambda[k], k, irrep);
	  sprintf(lbl, "Lia^(%d)[%d]", k, irrep);
	  dpd_file2_init(&L, CC_MISC, irrep, 0, 1, lbl);
	  dpd_file2_dirprd(&L, &F);
	  dpd_file2_close(&L);
	  
	  norm = dpd_file2_dot_self(&F);
	  norm = sqrt(norm);
	  if(norm > 1e-6)
	    dpd_file2_scm(&F, 1/norm);
	  else
	    dpd_file2_scm(&F, 0.0);
	  
	  dpd_file2_close(&F);
	}
	
	for(k = 0;k < nroots;k++){
	  sprintf(lbl, "F^(%d)[%d]", k, irrep);
	  dpd_file2_init(&F, CC_OEI, irrep, 0, 1, lbl);
	  sprintf(lbl, "B^(%d)[%d]", length, irrep);
	  dpd_file2_copy(&F, CC_OEI, lbl);
	  for(I = 0;I < length;I++){
	    sprintf(lbl, "B^(%d)[%d]", I, irrep);
	    dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	    coeff = - dpd_file2_dot(&F, &B);
	    sprintf(lbl, "B^(%d)[%d]", length, irrep);
	    dpd_file2_init(&Bp, CC_OEI, irrep, 0, 1, lbl);
	    dpd_file2_axpy(&B, &Bp, coeff, 0);
	    dpd_file2_close(&B);
	    dpd_file2_close(&Bp);
	  }
	  sprintf(lbl, "B^(%d)[%d]", length, irrep);
	  dpd_file2_init(&Bp, CC_OEI, irrep, 0, 1, lbl);
	  norm = dpd_file2_dot_self(&Bp);
	  norm = sqrt(norm);
	  if(norm < NORM_TOL)
	    dpd_file2_scm(&Bp, 0.0);
	  else
	    dpd_file2_scm(&Bp, 1/norm);
	  dpd_file2_close(&Bp);
	  length++;
	}
	
	if(maxdim - length < nroots){
	  fprintf(iter_adc, "Subspace too large: maxdim = %d, L = %d\n", maxdim, length);
	  fprintf(iter_adc, "Collapsing eigenvectors.\n");
	  
	  for(k = 0;k < nroots;k++){
	    sprintf(lbl, "B_n^(%d)[%d]", k, irrep);
	    dpd_file2_init(&B_n, CC_OEI, irrep, 0, 1, lbl);
	    dpd_file2_scm(&B_n, 0.0);
	    for(I = 0;I < length;I++){
	      sprintf(lbl, "B^(%d)[%d]", I, irrep);
	      dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	      dpd_file2_axpy(&B, &B_n, Alpha[I][k], 0);
	      dpd_file2_close(&B);
	    }
	    dpd_file2_close(&B_n);
	  }
	  for(k = 0;k < nroots;k++){
	    sprintf(lbl, "B_n^(%d)[%d]", k, irrep);
	    dpd_file2_init(&B_n, CC_OEI, irrep, 0, 1, lbl);
	    sprintf(lbl, "B^(%d)[%d]", k, irrep);
	    dpd_file2_copy(&B_n, CC_OEI, lbl);
	    dpd_file2_close(&B_n);
	  }
	  
	  skip_check = 1;
	  length = nroots;
	}
	
	if(!skip_check){
	  converged = 0;
	  zero_int_array(conv, nroots);
	  fprintf(iter_adc, "Root      Eigenvalue       Delta  Converged?\n");
	  fprintf(iter_adc, "---- -------------------- ------- ----------\n");
	  
	  for(k = 0;k < nroots;k++){
	    diff = fabs(lambda[k]-lambda_o[k]);
	    if(diff < cutoff){
	      conv[k] = 1;
	      converged++;
	    }
	    lambda_o[k] = lambda[k];
	    fprintf(iter_adc, "%3d  %20.14f %4.3e    %1s\n", k, lambda[k], diff,
		    conv[k] == 1 ? "Y" : "N");
	    fflush(iter_adc);
	  }
	}
	if(conv[root]){
	  eps[root] = lambda[root];
	  sprintf(lbl, "V^(%d)[%d]", root, irrep);
	  dpd_file2_init(&V, CC_OEI, irrep, 0, 1, lbl);
	  dpd_file2_scm(&V, 0.0);
	  for(I = 0;I < length;I++){
	    sprintf(lbl, "B^(%d)[%d]", I, irrep);
	    dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	    dpd_file2_axpy(&B, &V, Alpha[I][root], 0);
	    dpd_file2_close(&B);
	  }
	  dpd_file2_close(&V);
	  fprintf(iter_adc, "Davidson algorithm converged in %d iterations for %dth root.\n", iter, root);
	  timer_off("SEM");
	  break;
	}
	timer_off("SEM");
	iter++;
      }
      timer_off("Total_SEM");
      
      fclose(iter_adc);
      free(conv);
      free_block(G);
      free_block(Alpha);
      free(lambda);
      free(lambda_o);
      
      return eps[root];
    }
    
}}

