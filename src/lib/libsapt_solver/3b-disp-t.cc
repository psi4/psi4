/* This function calculates the Exch10 energy */

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
#include "structs.h"
#include "sapt3bn7.h"

namespace psi { namespace sapt {

void SAPT3BN7::disp211_T()
{
  results_.disp211t = disp211_T_0(PSIF_3B_SAPT_AA_DF_INTS,"AA RI Integrals",
    "AR RI Integrals","RR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BB RI Integrals","BS RI Integrals","SS RI Integrals",
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",calc_info_.evalsA,
    calc_info_.evalsB,calc_info_.evalsC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_211(T)         = %18.12lf  H\n",results_.disp211t);
  fflush(outfile);
  results_.disp121t = disp211_T_0(PSIF_3B_SAPT_CC_DF_INTS,"CC RI Integrals",
    "CT RI Integrals","TT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AA RI Integrals","AR RI Integrals","RR RI Integrals",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",calc_info_.evalsC,
    calc_info_.evalsA,calc_info_.evalsB,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"disp_121(T)         = %18.12lf  H\n",results_.disp121t);
  fflush(outfile);
  results_.disp112t = disp211_T_0(PSIF_3B_SAPT_BB_DF_INTS,"BB RI Integrals",
    "BS RI Integrals","SS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CC RI Integrals","CT RI Integrals","TT RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.evalsB,
    calc_info_.evalsC,calc_info_.evalsA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"disp_112(T)         = %18.12lf  H\n\n",results_.disp112t);
  fflush(outfile);
}

void SAPT3BN7::disp220_T()
{
  results_.disp220t = disp220_T_0(PSIF_3B_SAPT_CC_DF_INTS,"CC RI Integrals",
    "CT RI Integrals","TT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AA RI Integrals","AR RI Integrals","RR RI Integrals",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",calc_info_.evalsC,
    calc_info_.evalsA,calc_info_.evalsB,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"disp_220(T)         = %18.12lf  H\n",results_.disp220t);
  fflush(outfile);
  results_.disp202t = disp220_T_0(PSIF_3B_SAPT_AA_DF_INTS,"AA RI Integrals",
    "AR RI Integrals","RR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BB RI Integrals","BS RI Integrals","SS RI Integrals",
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",calc_info_.evalsA,
    calc_info_.evalsB,calc_info_.evalsC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_202(T)         = %18.12lf  H\n",results_.disp202t);
  fflush(outfile);
  results_.disp022t = disp220_T_0(PSIF_3B_SAPT_BB_DF_INTS,"BB RI Integrals",
    "BS RI Integrals","SS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CC RI Integrals","CT RI Integrals","TT RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.evalsB,
    calc_info_.evalsC,calc_info_.evalsA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"disp_022(T)         = %18.12lf  H\n\n",results_.disp022t);
  fflush(outfile);
}

double SAPT3BN7::disp211_T_0(int AAfile, char *AA_ints, char *AR_ints,
  char *RR_ints, int BBfile, char *BB_ints, char *BS_ints, char *SS_ints,
  int CCfile, char *CT_ints, double *e_A, double *e_B, double *e_C, int occA,
  int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **W_ARBS = block_matrix(occA*virA,occB*virB);
  double **Wt_BSAR = block_matrix(occB*virB,occA*virA);

  double **tARBS = IJKL_ints(AAfile,AR_ints,occA*virA,BBfile,BS_ints,
    occB*virB);
  double **tCTBS = IJKL_ints(CCfile,CT_ints,occC*virC,BBfile,BS_ints,
    occB*virB);

  for (int a=0,ar=0; a<occA; a++) {
  for (int r=0; r<virA; r++,ar++) {
    for (int b=0,bs=0; b<occB; b++) {
    for (int s=0; s<virB; s++,bs++) {
      double denom = e_A[a] + e_B[b] - e_A[r+occA] - e_B[s+occB];
      tARBS[ar][bs] /= denom;
  }}}}

  for (int c=0,ct=0; c<occC; c++) {
  for (int t=0; t<virC; t++,ct++) {
    for (int b=0,bs=0; b<occB; b++) {
    for (int s=0; s<virB; s++,bs++) {
      double denom = e_C[c] + e_B[b] - e_C[t+occC] - e_B[s+occB];
      tCTBS[ct][bs] /= denom;
  }}}}

  double **tBSAR = IJKL_ints(BBfile,BS_ints,occB*virB,AAfile,AR_ints,
    occA*virA);
  double **tCTAR = IJKL_ints(CCfile,CT_ints,occC*virC,AAfile,AR_ints,
    occA*virA);

  for (int b=0,bs=0; b<occB; b++) {
  for (int s=0; s<virB; s++,bs++) {
    for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      double denom = e_A[a] + e_B[b] - e_A[r+occA] - e_B[s+occB];
      tBSAR[bs][ar] /= denom;
  }}}}

  for (int c=0,ct=0; c<occC; c++) {
  for (int t=0; t<virC; t++,ct++) {
    for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      double denom = e_C[c] + e_A[a] - e_C[t+occC] - e_A[r+occA];
      tCTAR[ct][ar] /= denom;
  }}}}

  double **ARBB = IJKL_ints(AAfile,AR_ints,occA*virA,BBfile,BB_ints,
    occB*occB);
  double **CTBB = IJKL_ints(CCfile,CT_ints,occC*virC,BBfile,BB_ints,
    occB*occB);
  double **BSAA = IJKL_ints(BBfile,BS_ints,occB*virB,AAfile,AA_ints,
    occA*occA);
  double **CTAA = IJKL_ints(CCfile,CT_ints,occC*virC,AAfile,AA_ints,
    occA*occA);

  double **B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
  double **B_p_CT = get_DF_ints(CCfile,CT_ints,occC*virC);

  double **B_p_RR = get_DF_ints(AAfile,RR_ints,virA*virA);
  double **B_p_SS = get_DF_ints(BBfile,SS_ints,virB*virB);

  double **C_p_AR = block_matrix(occA*virA,calc_info_.nrio);
  double **C_p_BS = block_matrix(occB*virB,calc_info_.nrio);

  double **vRR = block_matrix(virA,virA);
  double **vSS = block_matrix(virB,virB);

  for (int c=0,ct=0; c<occC; c++) {
  for (int t=0; t<virC; t++,ct++) {

    C_DGEMV('n',virA*virA,calc_info_.nrio,1.0,&(B_p_RR[0][0]),calc_info_.nrio,
      &(B_p_CT[ct][0]),1,0.0,&(vRR[0][0]),1);

    C_DGEMV('n',virB*virB,calc_info_.nrio,1.0,&(B_p_SS[0][0]),calc_info_.nrio,
      &(B_p_CT[ct][0]),1,0.0,&(vSS[0][0]),1);

    C_DGEMM('N','N',occA*virA*occB,virB,virB,1.0,&(tARBS[0][0]),virB,
      &(vSS[0][0]),virB,0.0,&(W_ARBS[0][0]),virB);

    C_DGEMM('N','N',occA*virA*occB,virB,occB,-1.0,&(ARBB[0][0]),occB,
      &(tCTBS[ct][0]),virB,1.0,&(W_ARBS[0][0]),virB);

    C_DGEMM('N','N',occB*virB*occA,virA,virA,1.0,&(tBSAR[0][0]),virA,
      &(vRR[0][0]),virA,0.0,&(Wt_BSAR[0][0]),virA);

    C_DGEMM('N','N',occB*virB*occA,virA,occA,-1.0,&(BSAA[0][0]),occA,
      &(tCTAR[ct][0]),virA,1.0,&(Wt_BSAR[0][0]),virA);

    for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      C_DGEMM('N','N',occB,virB,occB,-1.0,&(CTBB[ct][0]),occB,&(tARBS[ar][0]),
        virB,1.0,&(W_ARBS[ar][0]),virB);
    }}

    for (int b=0,bs=0; b<occB; b++) {
    for (int s=0; s<virB; s++,bs++) {
      C_DGEMM('N','N',occA,virA,occA,-1.0,&(CTAA[ct][0]),occA,&(tBSAR[bs][0]),
        virA,1.0,&(Wt_BSAR[bs][0]),virA);
    }}

    C_DGEMM('N','N',occB,virB*calc_info_.nrio,virB,1.0,&(tCTBS[ct][0]),virB,
      &(B_p_SS[0][0]),virB*calc_info_.nrio,0.0,&(C_p_BS[0][0]),
      virB*calc_info_.nrio);

    C_DGEMM('N','T',occA*virA,occB*virB,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
      calc_info_.nrio,&(C_p_BS[0][0]),calc_info_.nrio,1.0,&(W_ARBS[0][0]),
      occB*virB);

    C_DGEMM('N','N',occA,virA*calc_info_.nrio,virA,1.0,&(tCTAR[ct][0]),virA,
      &(B_p_RR[0][0]),virA*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
      virA*calc_info_.nrio);

    C_DGEMM('N','T',occB*virB,occA*virA,calc_info_.nrio,1.0,&(B_p_BS[0][0]),
      calc_info_.nrio,&(C_p_AR[0][0]),calc_info_.nrio,1.0,&(Wt_BSAR[0][0]),
      occA*virA);

    for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      for (int b=0,bs=0; b<occB; b++) {
      for (int s=0; s<virB; s++,bs++) {
        double denom = e_A[a] + e_B[b] + e_C[c] -
          e_A[r+occA] - e_B[s+occB] - e_C[t+occC];
        energy += W_ARBS[ar][bs]*Wt_BSAR[bs][ar]/denom;
    }}}}

  }}

  free_block(W_ARBS);
  free_block(Wt_BSAR);
  free_block(tARBS);
  free_block(tBSAR);
  free_block(tCTBS);
  free_block(tCTAR);
  free_block(BSAA);
  free_block(CTAA);
  free_block(ARBB);
  free_block(CTBB);
  free_block(B_p_AR);
  free_block(B_p_BS);
  free_block(B_p_CT);
  free_block(B_p_RR);
  free_block(B_p_SS);
  free_block(C_p_AR);
  free_block(C_p_BS);
  free_block(vRR);
  free_block(vSS);

  return(16.0*energy);
}

double SAPT3BN7::disp220_T_0(int AAfile, char *AA_ints, char *AR_ints,
  char *RR_ints, int BBfile, char *BB_ints, char *BS_ints, char *SS_ints,
  int CCfile, char *CT_ints, double *e_A, double *e_B, double *e_C, int occA,
  int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **W_ARBS = block_matrix(occA*virA,occB*virB);

  double **tARBS = IJKL_ints(AAfile,AR_ints,occA*virA,BBfile,BS_ints,
    occB*virB);
  double **tCTBS = IJKL_ints(CCfile,CT_ints,occC*virC,BBfile,BS_ints,
    occB*virB);

  for (int a=0,ar=0; a<occA; a++) {
  for (int r=0; r<virA; r++,ar++) {
    for (int b=0,bs=0; b<occB; b++) {
    for (int s=0; s<virB; s++,bs++) {
      double denom = e_A[a] + e_B[b] - e_A[r+occA] - e_B[s+occB];
      tARBS[ar][bs] /= denom;
  }}}}

  for (int c=0,ct=0; c<occC; c++) {
  for (int t=0; t<virC; t++,ct++) {
    for (int b=0,bs=0; b<occB; b++) {
    for (int s=0; s<virB; s++,bs++) {
      double denom = e_C[c] + e_B[b] - e_C[t+occC] - e_B[s+occB];
      tCTBS[ct][bs] /= denom;
  }}}}

  double **ARBB = IJKL_ints(AAfile,AR_ints,occA*virA,BBfile,BB_ints,
    occB*occB);
  double **CTBB = IJKL_ints(CCfile,CT_ints,occC*virC,BBfile,BB_ints,
    occB*occB);

  double **B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  double **B_p_CT = get_DF_ints(CCfile,CT_ints,occC*virC);

  double **B_p_SS = get_DF_ints(BBfile,SS_ints,virB*virB);

  double **C_p_BS = block_matrix(occB*virB,calc_info_.nrio);

  double **vSS = block_matrix(virB,virB);

  for (int c=0,ct=0; c<occC; c++) {
  for (int t=0; t<virC; t++,ct++) {

    C_DGEMV('n',virB*virB,calc_info_.nrio,1.0,&(B_p_SS[0][0]),calc_info_.nrio,
      &(B_p_CT[ct][0]),1,0.0,&(vSS[0][0]),1);

    C_DGEMM('N','N',occA*virA*occB,virB,virB,1.0,&(tARBS[0][0]),virB,
      &(vSS[0][0]),virB,0.0,&(W_ARBS[0][0]),virB);

    C_DGEMM('N','N',occA*virA*occB,virB,occB,-1.0,&(ARBB[0][0]),occB,
      &(tCTBS[ct][0]),virB,1.0,&(W_ARBS[0][0]),virB);

    for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      C_DGEMM('N','N',occB,virB,occB,-1.0,&(CTBB[ct][0]),occB,&(tARBS[ar][0]),
        virB,1.0,&(W_ARBS[ar][0]),virB);
    }}

    C_DGEMM('N','N',occB,virB*calc_info_.nrio,virB,1.0,&(tCTBS[ct][0]),virB,
      &(B_p_SS[0][0]),virB*calc_info_.nrio,0.0,&(C_p_BS[0][0]),
      virB*calc_info_.nrio);

    C_DGEMM('N','T',occA*virA,occB*virB,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
      calc_info_.nrio,&(C_p_BS[0][0]),calc_info_.nrio,1.0,&(W_ARBS[0][0]),
      occB*virB);

    for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      for (int b=0,bs=0; b<occB; b++) {
      for (int s=0; s<virB; s++,bs++) {
        double denom = e_A[a] + e_B[b] + e_C[c] -
          e_A[r+occA] - e_B[s+occB] - e_C[t+occC];
        energy += W_ARBS[ar][bs]*W_ARBS[ar][bs]/denom;
    }}}}

  }}

  free_block(W_ARBS);
  free_block(tARBS);
  free_block(tCTBS);
  free_block(ARBB);
  free_block(CTBB);
  free_block(B_p_AR);
  free_block(B_p_CT);
  free_block(B_p_SS);
  free_block(C_p_BS);
  free_block(vSS);

  return(8.0*energy);
}

}}
