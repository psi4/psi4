/* This function calculates the Ind20,resp energy */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libipv1/ip_data.gbl>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "sapt2p3.h"

namespace psi { namespace sapt {

void SAPT2p3::ind_disp30()
{
  double e1,e2,e3;

  if (params_.print) 
    fprintf(outfile,"Begining Ind-Disp30 Calculation\n\n");

  ind_disp_ov("Ind-Disp30 AR Amplitudes","T(BS) AR",PSIF_SAPT_AA_DF_INTS,
    "AA RI Integrals","RR RI Integrals",calc_info_.sA,calc_info_.evalsA,
    calc_info_.noccA,calc_info_.nvirA);

  ind_disp_ov("Ind-Disp30 BS Amplitudes","T(AR) BS",PSIF_SAPT_BB_DF_INTS,
    "BB RI Integrals","SS RI Integrals",calc_info_.sB,calc_info_.evalsB,
    calc_info_.noccB,calc_info_.nvirB);

  ind_disp_ovov();

  double **X_AR = read_IJKL(PSIF_SAPT_AMPS,"Ind-Disp30 AR Amplitudes",
    calc_info_.noccA,calc_info_.nvirA);

  e1 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,X_AR[0],1,
    calc_info_.WBAR[0],1);

  double **X_BS = read_IJKL(PSIF_SAPT_AMPS,"Ind-Disp30 BS Amplitudes",
    calc_info_.noccB,calc_info_.nvirB);

  e2 = 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,X_BS[0],1,
    calc_info_.WABS[0],1);

  free_block(X_AR);
  free_block(X_BS);

  double **X_ARBS = read_IJKL(PSIF_SAPT_AMPS,"Ind-Disp30 ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);

  double **B_p_AR = get_AR_ints(1);
  double **B_p_BS = get_BS_ints(1);
  double **ARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nri,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(ARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  e3 = 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nvirB*
    calc_info_.noccB,&(ARBS[0][0]),1,&(X_ARBS[0][0]),1);

  free_block(X_ARBS);
  free_block(ARBS);

  results_.ind_disp30 = e1+e2+e3;
  if (params_.print) {
    fprintf(outfile,"Ind-Disp30         = %18.12lf  H\n\n",
      results_.ind_disp30);
    fflush(outfile);
  }
}

void SAPT2p3::ind_disp_ov(char *amp_out, char *theta_label, int AAfile, 
  char *AAlabel, char *RRlabel, double **sAR, double *evals, int nocc, 
  int nvir)
{
  double **iAR = block_matrix(nocc,nvir);

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,theta_label,nocc*nvir,
    calc_info_.nrio);
  double **B_p_RR = get_DF_ints(AAfile,RRlabel,nvir*nvir);

  C_DGEMM('N','T',nocc,nvir,nvir*calc_info_.nrio,2.0,&(T_p_AR[0][0]),
    nvir*calc_info_.nrio,&(B_p_RR[0][0]),nvir*calc_info_.nrio,0.0,&(iAR[0][0]),
    nvir);

  free_block(B_p_RR);

  double **B_p_AA = get_DF_ints(AAfile,AAlabel,nocc*nocc);

  for(int a=0; a<nocc; a++) {
    C_DGEMM('N','T',nocc,nvir,calc_info_.nrio,-2.0,&(B_p_AA[a*nocc][0]),
      calc_info_.nrio,&(T_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(iAR[0][0]),
      nvir);
  }

  free_block(B_p_AA);
  free_block(T_p_AR);

  for(int a=0; a<nocc; a++) {
    for(int r=0; r<nvir; r++) {
      double denom = evals[a] - evals[r+nocc];
      iAR[a][r] /= denom;
  }}

  write_IJKL(iAR,PSIF_SAPT_AMPS,amp_out,nocc,nvir);
}

void SAPT2p3::ind_disp_ovov()
{
  double **uARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);

  double **C_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);
  double **B_p_RR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",
    calc_info_.nvirA*calc_info_.nvirA);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.nvirA,1.0,calc_info_.sA[0],calc_info_.nvirA,&(B_p_RR[0][0]),
    calc_info_.nvirA*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio);

  free_block(B_p_RR);

  double **B_p_AA = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    calc_info_.noccA*calc_info_.noccA);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.noccA,-1.0,
      calc_info_.sA[0],calc_info_.nvirA,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,1.0,&(C_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio);
  }

  free_block(B_p_AA);

  double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(C_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(uARBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB);

  free_block(C_p_AR);
  free_block(B_p_BS);

  double **C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);
  double **B_p_SS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",
    calc_info_.nvirB*calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.nvirB,1.0,calc_info_.sB[0],calc_info_.nvirB,&(B_p_SS[0][0]),
    calc_info_.nvirB*calc_info_.nrio,0.0,&(C_p_BS[0][0]),calc_info_.nvirB*
    calc_info_.nrio);

  free_block(B_p_SS);

  double **B_p_BB = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    calc_info_.noccB*calc_info_.noccB);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccB,-1.0,
      calc_info_.sB[0],calc_info_.nvirB,&(B_p_BB[b*calc_info_.noccB][0]),
      calc_info_.nrio,1.0,&(C_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  free_block(B_p_BB);

  double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    calc_info_.noccA*calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(C_p_BS[0][0]),calc_info_.nrio,1.0,&(uARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(B_p_AR);
  free_block(C_p_BS);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA*calc_info_.nvirB*
    calc_info_.noccB,calc_info_.noccA,-1.0,&(calc_info_.WBAA[0][0]),
    calc_info_.noccA,&(tARBS[0][0]),calc_info_.nvirA*calc_info_.nvirB*
    calc_info_.noccB,1.0,&(uARBS[0][0]),calc_info_.nvirA*calc_info_.nvirB*
    calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.nvirA*calc_info_.noccA*calc_info_.noccB,
    calc_info_.nvirB,calc_info_.nvirB,1.0,&(tARBS[0][0]),calc_info_.nvirB,
    &(calc_info_.WASS[0][0]),calc_info_.nvirB,1.0,&(uARBS[0][0]),
    calc_info_.nvirB);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB,
      calc_info_.nvirA,1.0,&(calc_info_.WBRR[0][0]),calc_info_.nvirA,
      &(tARBS[a*calc_info_.nvirA][0]),calc_info_.noccB*calc_info_.nvirB,1.0,
      &(uARBS[a*calc_info_.nvirA][0]),calc_info_.noccB*calc_info_.nvirB);
  }

  for(int a=0, ar=0; a<calc_info_.noccA; a++) {
    for(int r=0; r<calc_info_.nvirA; r++,ar++) {
      C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccB,-1.0,
        &(calc_info_.WABB[0][0]),calc_info_.noccB,&(tARBS[ar][0]),
        calc_info_.nvirB,1.0,&(uARBS[ar][0]),calc_info_.nvirB);
  }}

  free_block(tARBS);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      for (int b=0,bs=0; b<calc_info_.noccB; b++) {
        for (int s=0; s<calc_info_.nvirB; s++,bs++) {
          double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
                  calc_info_.evalsA[r+calc_info_.noccA]-
                  calc_info_.evalsB[s+calc_info_.noccB];
          uARBS[ar][bs] /= denom;
        }}
    }}

  write_IJKL(uARBS,PSIF_SAPT_AMPS,"Ind-Disp30 ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
}

}}

