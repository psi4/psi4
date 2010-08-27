/* This function calculates the Disp20 energy */

#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <cstring>
#include <iostream>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libipv1/ip_data.gbl>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt0.h"

namespace psi { namespace sapt {

void SAPT0::df_disp20_chf()
{
  if (params_.print)
    fprintf(outfile,"Begining Disp20 (CHF) Calculation\n\n");

  double x_k[] = {
    0.0876494104789,
    0.462696328915,
    1.14105777483,
    2.12928364510,
    3.43708663389,
    5.07801861456,
    7.07033853501,
    9.43831433649,
    12.2142233686,
    15.4415273694,
    19.1801568554,
    23.5159056959,
    28.5787297412,
    34.5833987032,
    41.9404526474,
    51.7011603396 };

  double w_k[] = {
    0.225036314864,
    0.525836052762,
    0.831961391687,
    1.14609924096,
    1.47175131698,
    1.81313468736,
    2.17551751937,
    2.56576274964,
    2.99321508344,
    3.47123449709,
    4.02004410019,
    4.67251661146,
    5.48742063616,
    6.58536124449,
    8.27635797549,
    11.8242775489 };

  double **A_PQ, **B_PQ;
  double energy = 0.0, tval;

  for (int i=0; i<16; i++) {

    A_PQ = DF_FDDS(x_k[i],PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
      PSIF_SAPT_LRINTS,"A LR Integrals",calc_info_.evalsA,
      calc_info_.noccA,calc_info_.nvirA);

    B_PQ = DF_FDDS(x_k[i],PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
      PSIF_SAPT_LRINTS,"B LR Integrals",calc_info_.evalsB,
      calc_info_.noccB,calc_info_.nvirB);

    tval = -C_DDOT(calc_info_.nrio*calc_info_.nrio,&(A_PQ[0][0]),1,
      &(B_PQ[0][0]),1)/(8.0*M_PI);

    free_block(A_PQ);
    free_block(B_PQ);

    energy += w_k[i]*tval;

    if (params_.print) {
      fprintf(outfile,"  i=%2d CHF Dispersion Energy  = %18.12lf mH\n",i,
        tval*4000.0);
      fflush(outfile);
    }
  }

  results_.disp20chf = 4.0*energy;

  if (params_.print) {
    fprintf(outfile,"\nCHF Dispersion Energy  = %18.12lf mH\n\n",
      results_.disp20chf*1000.0);
    fflush(outfile);
  }
}

double **SAPT0::D_lambda_F(double omega, int file1, char *Dints,
  int file2, char *Fints, double *evals, int nocc, int nvir)
{
  double **lambda = block_matrix(nocc,nvir);

  for(int a1=0; a1<nocc; a1++) {
    for(int r1=0; r1<nvir; r1++) {
      double tval = evals[r1+nocc] - evals[a1];
      lambda[a1][r1] = -4.0*tval/(tval*tval + omega*omega);
  }}

  double **F_p_AR = get_DF_ints(file2,Fints,nocc*nvir);

  for(int a1=0,a1r1=0; a1<nocc; a1++) {
    for(int r1=0; r1<nvir; r1++,a1r1++) {
      C_DSCAL(calc_info_.nrio,lambda[a1][r1],F_p_AR[a1r1],1);
  }}

  free_block(lambda);

  double **D_p_AR = get_DF_ints(file1,Dints,nocc*nvir);
  double **C_pq = block_matrix(calc_info_.nrio,calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.nrio,calc_info_.nrio,nocc*nvir,1.0,
    &(D_p_AR[0][0]),calc_info_.nrio,&(F_p_AR[0][0]),calc_info_.nrio,0.0,
    &(C_pq[0][0]),calc_info_.nrio);

  free_block(F_p_AR);
  free_block(D_p_AR);

  return(C_pq);
}

double **SAPT0::DF_FDDS(double omega, int DFfile, char *AR_ints, 
  int LRfile, char *F_ints, double *evals, int nocc, int nvir)
{
  double **DyF = D_lambda_F(omega,DFfile,AR_ints,LRfile,F_ints,evals,
                            nocc,nvir);
  double **I = block_matrix(calc_info_.nrio,calc_info_.nrio);
  double **inv = block_matrix(calc_info_.nrio,calc_info_.nrio);

  for(int P=0; P<calc_info_.nrio; P++) {
    I[P][P] = 1.0;
  }

  C_DAXPY(calc_info_.nrio*calc_info_.nrio,-1.0,DyF[0],1,I[0],1);
  invert_matrix(I,inv,calc_info_.nrio,outfile);

  C_DGEMM('N','N',calc_info_.nrio,calc_info_.nrio,calc_info_.nrio,1.0,
    &(DyF[0][0]),calc_info_.nrio,&(inv[0][0]),calc_info_.nrio,0.0,
    &(I[0][0]),calc_info_.nrio);

  free_block(DyF);
  free_block(inv);

  double **DyD = D_lambda_F(omega,DFfile,AR_ints,DFfile,AR_ints,evals,
                            nocc,nvir);
  double **X = block_matrix(calc_info_.nrio,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.nrio,calc_info_.nrio,calc_info_.nrio,1.0,
    &(I[0][0]),calc_info_.nrio,&(DyD[0][0]),calc_info_.nrio,0.0,
    &(X[0][0]),calc_info_.nrio);

  C_DAXPY(calc_info_.nrio*calc_info_.nrio,1.0,DyD[0],1,X[0],1);

  free_block(DyD);
  free_block(I);

  return(X);
}

/*
void SAPT0::disp20_chf()
{
  fprintf(outfile,"Begining Disp20 (CHF) Calculation\n\n");

  double x_k[] = {
    0.0876494104789,
    0.462696328915,
    1.14105777483,
    2.12928364510,
    3.43708663389,
    5.07801861456,
    7.07033853501,
    9.43831433649,
    12.2142233686,
    15.4415273694,
    19.1801568554,
    23.5159056959,
    28.5787297412,
    34.5833987032,
    41.9404526474,
    51.7011603396 };

  double w_k[] = {
    0.225036314864,
    0.525836052762,
    0.831961391687,
    1.14609924096,
    1.47175131698,
    1.81313468736,
    2.17551751937,
    2.56576274964,
    2.99321508344,
    3.47123449709,
    4.02004410019,
    4.67251661146,
    5.48742063616,
    6.58536124449,
    8.27635797549,
    11.8242775489 };

  double **B_p_AR = get_AR_ints(0);
  double **B_p_BS = get_BS_ints(0);
  double **ARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,&(B_p_BS[0][0]),
    calc_info_.nrio,0.0,&(ARBS[0][0]),calc_info_.noccB*calc_info_.nvirB);

  free_block(B_p_AR);
  free_block(B_p_BS);

  double **C_ARAR, **C_BSBS;
  double **temp = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB);
  double **tARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB);

  double energy = 0.0, tval;

  for (int i=0; i<16; i++) {

    C_ARAR = FDDS(x_k[i],PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR RI Integrals","RR RI Integrals",calc_info_.evalsA,calc_info_.noccA,
      calc_info_.nvirA);

    C_BSBS = FDDS(x_k[i],PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS RI Integrals","SS RI Integrals",calc_info_.evalsB,calc_info_.noccB,
      calc_info_.nvirB);

    C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
      calc_info_.nvirB,calc_info_.noccA*calc_info_.nvirA,1.0,&(C_ARAR[0][0]),
      calc_info_.noccA*calc_info_.nvirA,&(ARBS[0][0]),calc_info_.noccB*
      calc_info_.nvirB,0.0,&(temp[0][0]),calc_info_.noccB*calc_info_.nvirB);

    C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
      calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB,-1.0/(8.0*M_PI),
      &(temp[0][0]),calc_info_.noccB*calc_info_.nvirB,&(C_BSBS[0][0]),
      calc_info_.noccB*calc_info_.nvirB,0.0,&(tARBS[0][0]),calc_info_.noccB*
      calc_info_.nvirB);

    free_block(C_ARAR);
    free_block(C_BSBS);

    tval = C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
                   &(ARBS[0][0]),1,&(tARBS[0][0]),1)*4.0;

    energy += w_k[i]*tval;

    fprintf(outfile,"  i=%2d CHF Dispersion Energy  = %18.12lf mH\n",i,
      tval*4000.0);
    fflush(outfile);
  }

  free_block(ARBS);
  free_block(tARBS);

  fprintf(outfile,"\nCHF Dispersion Energy  = %18.12lf mH\n\n",energy*4000.0);
  fflush(outfile);
}

double **SAPT0::FDDS(double omega, int DFfile, char *AA_ints, char *AR_ints,
  char *RR_ints, double *evals, int nocc, int nvir)
{
  double **B_p_AA = get_DF_ints(DFfile,AA_ints,nocc*nocc);
  double **B_p_AR = get_DF_ints(DFfile,AR_ints,nocc*nvir);
  double **B_p_RR = get_DF_ints(DFfile,RR_ints,nvir*nvir);

  double **ARAR = block_matrix(nocc*nvir,nocc*nvir);
  double **AARR = block_matrix(nocc*nocc,nvir*nvir);

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,&(B_p_AR[0][0]),calc_info_.nrio,0.0,
    &(ARAR[0][0]),nocc*nvir);

  C_DGEMM('N','T',nocc*nocc,nvir*nvir,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,&(B_p_RR[0][0]),calc_info_.nrio,0.0,
    &(AARR[0][0]),nvir*nvir);

  free_block(B_p_AA);
  free_block(B_p_AR);
  free_block(B_p_RR);

  double **H1 = block_matrix(nocc*nvir,nocc*nvir);
  double **H2 = block_matrix(nocc*nvir,nocc*nvir);

  for(int a1=0, a1r1=0; a1<nocc; a1++) {
    for(int r1=0; r1<nvir; r1++,a1r1++) {
      H1[a1r1][a1r1] += evals[r1+nocc] - evals[a1];
      H2[a1r1][a1r1] += evals[r1+nocc] - evals[a1];
      for(int a2=0, a2r2=0; a2<nocc; a2++) {
        for(int r2=0; r2<nvir; r2++,a2r2++) {
          int a1a2 = a1*nocc+a2;
          int r1r2 = r1*nvir+r2;
          int a1r2 = a1*nvir+r2;
          int a2r1 = a2*nvir+r1;
          H1[a1r1][a2r2] += 4.0*ARAR[a1r1][a2r2]; 
//                            - AARR[a1a2][r1r2] - ARAR[a1r2][a2r1];
//          H2[a1r1][a2r2] -= AARR[a1a2][r1r2] - ARAR[a1r2][a2r1];
  }}}}

  free_block(ARAR);
  free_block(AARR);

  double **temp = block_matrix(nocc*nvir,nocc*nvir);

  C_DGEMM('N','N',nocc*nvir,nocc*nvir,nocc*nvir,1.0,
    &(H2[0][0]),nocc*nvir,&(H1[0][0]),nocc*nvir,0.0,
    &(temp[0][0]),nocc*nvir);

  for(int a1=0, a1r1=0; a1<nocc; a1++) {
    for(int r1=0; r1<nvir; r1++,a1r1++) {
      temp[a1r1][a1r1] += omega*omega;
  }}

  invert_matrix(temp,H1,nocc*nvir,outfile);

  C_DGEMM('N','N',nocc*nvir,nocc*nvir,nocc*nvir,-2.0,
    &(H1[0][0]),nocc*nvir,&(H2[0][0]),nocc*nvir,0.0,
    &(temp[0][0]),nocc*nvir);

  free_block(H1);
  free_block(H2);

  return(temp);
}
*/
}}
