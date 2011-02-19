/* This function calculates the Ind20,resp amplitudes */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "sapt2b.h"

namespace psi { namespace sapt {

void SAPT2B::cphf_induction()
{
  calc_info_.sA = uchf_ind(calc_info_.WBAR,calc_info_.evalsA,calc_info_.noccA,
    calc_info_.nvirA);
  if (workflow_.save_chf) {
    if (params_.print) {
      fprintf(outfile,"Solving CHF equations for monomer A\n\n");
      fflush(outfile);
    }
    calc_info_.CHFA = cphf_ind(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR RI Integrals", "RR RI Integrals", calc_info_.WBAR, calc_info_.sA, 
      calc_info_.evalsA, calc_info_.noccA, calc_info_.nvirA);
  }
  if (!workflow_.save_s)
    free_block(calc_info_.sA);

  calc_info_.sB = uchf_ind(calc_info_.WABS,calc_info_.evalsB,calc_info_.noccB, 
    calc_info_.nvirB);
  if (workflow_.save_chf) {
    if (params_.print) {
      fprintf(outfile,"Solving CHF equations for monomer B\n\n");
      fflush(outfile);
    }
    calc_info_.CHFB = cphf_ind(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS RI Integrals", "SS RI Integrals", calc_info_.WABS, calc_info_.sB,
      calc_info_.evalsB, calc_info_.noccB, calc_info_.nvirB);
  }
  if (!workflow_.save_s)
    free_block(calc_info_.sB);
}

double **SAPT::uchf_ind(double **W, double *evals, int nocc, int nvir)
{
  double **S = block_matrix(nocc,nvir);

  C_DCOPY(nocc*nvir,W[0],1,S[0],1);

  #pragma omp for schedule(static)
  for (int a=0; a < nocc; a++) {
    for (int r=0; r < nvir; r++) {
      double denom = evals[a] - evals[r+nocc];
      S[a][r] /= denom;
  }}

  return(S);
}

double **SAPT::cphf_ind(int dfnum, const char *OO, const char *OV, 
  const char *VV, double **W, double **S, double *evals, int nocc, int nvir)
{
  int a,r,ar;
  time_t start = time(NULL);
  time_t stop;
  int iter=0;
  double conv, tval, denom, E_old, E;
  double **C, **C_old;
  double **resC, **oldC;

  resC = block_matrix(nocc*nvir,params_.diisvec);
  oldC = block_matrix(nocc*nvir,params_.diisvec);

  C = block_matrix(nocc,nvir);
  C_old = block_matrix(nocc,nvir);

  C_DCOPY(nocc*nvir,S[0],1,C_old[0],1);

  E_old = 2.0*C_DDOT(nocc*nvir,&(C_old[0][0]),1,&(W[0][0]),1);

  conv = 1.0;

  if (params_.print)
    fprintf(outfile,"Iter      Energy (mH)         dE (mH)            RMS (mH)    Time (s)\n");

  do {

    if (iter > params_.diisvec-1) {
      diis_update(C_old, resC, oldC, nocc, nvir);
      }

  A_mat(dfnum, OO, OV, VV, C_old, C, nocc, nvir, iter);

  #pragma omp for private(a,r,denom) schedule(static)
  for (a=0; a<nocc; a++) {
    for (r=0; r<nvir; r++) {
      denom = evals[a] - evals[r+nocc];
      C[a][r] += W[a][r];
      C[a][r] /= denom;
  }}

    E = 2.0*C_DDOT(nocc*nvir,&(C[0][0]),1,&(W[0][0]),1);

    conv = 0.0;

    for (int a=0, i=0; a < nocc; a++) {
    for (int r=0; r < nvir; r++,i++) {
      resC[i][iter % params_.diisvec] = C[a][r] - C_old[a][r];
      conv += pow(C[a][r] - C_old[a][r],2);
      }}

    conv = sqrt(conv/(nocc*nvir));

    C_DCOPY(nocc*nvir,&(C[0][0]),1,&(C_old[0][0]),1);
    C_DCOPY(nocc*nvir,&(C[0][0]),1,&(oldC[0][iter % params_.diisvec]),
            params_.diisvec);

    iter++;
    stop = time(NULL);
    if (params_.print) {
      fprintf(outfile,"%4d %16.8lf %17.9lf %17.9lf    %10ld\n",iter,E*1000.0,
        (E_old-E)*1000.0,conv*1000.0,stop-start);
      fflush(outfile);
    }

    E_old = E;
    }
  while(conv > params_.d_conv && iter < params_.maxiter);

  if (conv <= params_.d_conv) {
    if (params_.print)
      fprintf(outfile,"\nCHF Iterations converged\n\n");
    }
  else {
    if (params_.print)
      fprintf(outfile,"\nCHF Iterations did not converge\n\n");
    }

  free_block(resC);
  free_block(oldC);
  free_block(C);

  fflush(outfile);

  return(C_old);
}

void SAPT::A_mat(int dfnum, const char *OO, const char *OV, const char *VV, 
  double **C_old, double **C_new, int nocc, int nvir, int iter)
{
  int nrio = ribasis_->nbf() + 3;

  double **B_p_AR = get_DF_ints(dfnum,OV,nocc*nvir);
  double *C_p = init_array(nrio);

  C_DGEMV('t',nocc*nvir,nrio,1.0,&(B_p_AR[0][0]),nrio,
          &(C_old[0][0]),1,0.0,C_p,1);
  C_DGEMV('n',nocc*nvir,nrio,4.0,&(B_p_AR[0][0]),nrio,
          C_p,1,0.0,&(C_new[0][0]),1);

  free(C_p);

  double **D_p_AA = block_matrix(nocc*nocc,nrio);

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','N',nocc,nrio,nvir,1.0,&(C_old[0][0]),nvir,
      &(B_p_AR[a*nvir][0]),nrio,0.0,&(D_p_AA[a][0]),
      nocc*nrio);
  }

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','T',nocc,nvir,nrio,-1.0,&(D_p_AA[a*nocc][0]),
      nrio,&(B_p_AR[a*nvir][0]),nrio,1.0,&(C_new[0][0]),
      nvir);
  }

  free_block(D_p_AA);
  memset(&(B_p_AR[0][0]),'\0',sizeof(double)*nocc*nvir*nrio);

  long int avail_mem = params_.memory;
  avail_mem -= 8*(nocc*nvir*(long int) nrio);

  long int temp_size = avail_mem / (8*nvir*(long int) nrio);
  
  if (temp_size > nvir)
    temp_size = nvir;

  int blocks = (nvir)/temp_size;
  if ((nvir)%temp_size) blocks++;
  
  if (temp_size < 1) {
    fprintf(outfile,"Not enough memory in A Matrix formation\n\n");
    exit(0);
  } 

  double **B_p_RR = block_matrix(temp_size*nvir,nrio);

  psio_address next_PSIF = PSIO_ZERO;
  for (int t_r=0; t_r<blocks; t_r++) {
    int r_start = temp_size*t_r;
    int r_stop = temp_size*(t_r+1);
    if (r_stop > nvir)
      r_stop = nvir;

    psio_->read(dfnum,VV,(char *) &(B_p_RR[0][0]),sizeof(double)*
      (r_stop-r_start)*nvir*(ULI) nrio,next_PSIF,&next_PSIF);
    for (int r=r_start; r<r_stop; r++) {
      C_DGEMM('N','N',nocc,nrio,nvir,1.0,&(C_old[0][0]),nvir,
        &(B_p_RR[(r-r_start)*nvir][0]),nrio,1.0,&(B_p_AR[r][0]),
        nvir*nrio);
    }

  }
  free_block(B_p_RR);

  double **B_p_AA = get_DF_ints(dfnum,OO,nocc*nocc);

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','T',nocc,nvir,nrio,-1.0,&(B_p_AA[a*nocc][0]),
      nrio,&(B_p_AR[a*nvir][0]),nrio,1.0,&(C_new[0][0]),
      nvir);
  }

  free_block(B_p_AA);
  free_block(B_p_AR);
}

void SAPT::diis_update(double **C, double **resC, double **oldC, int nocc, int nvir)
{
  int *ipiv;
  double *Cvec;
  double **Bmat;

  ipiv = init_int_array(params_.diisvec+1);
  Bmat = block_matrix(params_.diisvec+1,params_.diisvec+1);
  Cvec = (double *) malloc((params_.diisvec+1)*sizeof(double));

  for (int i=0; i<params_.diisvec; i++) {
    for (int j=0; j<=i; j++) {
      Bmat[i][j] = Bmat[j][i] = C_DDOT(nocc*nvir,&(resC[0][i]),params_.diisvec,
                          &(resC[0][j]),params_.diisvec);
      }}

  for (int i=0; i<params_.diisvec; i++) {
    Bmat[params_.diisvec][i] = -1.0;
    Bmat[i][params_.diisvec] = -1.0;
    Cvec[i] = 0.0;
    }

  Bmat[params_.diisvec][params_.diisvec] = 0.0;
  Cvec[params_.diisvec] = -1.0;

  C_DGESV(params_.diisvec+1,1,&(Bmat[0][0]),params_.diisvec+1,
          &(ipiv[0]),&(Cvec[0]),params_.diisvec+1);

  for (int i=0,m=0; i<nocc; i++) {
    for (int j=0; j<nvir; j++,m++) {
      C[i][j] = C_DDOT(params_.diisvec,&(Cvec[0]),1,&(oldC[m][0]),1);
        }}

  free(ipiv);
  free(Cvec);
  free_block(Bmat);
}

}}
