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
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::ind22()
{
  double ind220, ind202;

  if (params_.print) {
    fprintf(outfile,"Begining Ind22 Calculation\n\n");
    fflush(outfile);
  }

  ind220 = ind22_1(calc_info_.sA,calc_info_.WBAA,calc_info_.WBRR,
    "T ARAR Amplitudes",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    "AR RI Integrals","RR RI Integrals",calc_info_.evalsA,calc_info_.noccA,
    calc_info_.nvirA);
  ind220 += ind22_2(calc_info_.sA,calc_info_.WBAA,calc_info_.WBRR,
    "T AR Amplitudes",calc_info_.noccA,calc_info_.nvirA);
  ind220 += ind22_3(calc_info_.sA,calc_info_.WBAR,"AA MP2 OPDM","RR MP2 OPDM",
    calc_info_.noccA,calc_info_.nvirA);
  ind220 += ind22_4(calc_info_.sA,"Theta(AR) AR",PSIF_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccA,calc_info_.nvirA);
  ind220 += ind22_5(calc_info_.sA,"T2 ARAR Antisym Amplitudes",
    calc_info_.evalsA,calc_info_.noccA,calc_info_.nvirA);
  ind220 += ind22_6(calc_info_.sA,"T ARAR Antisym Amplitudes",
    PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals","RR RI Integrals",
    calc_info_.noccA,calc_info_.nvirA);
  ind220 += ind22_7(calc_info_.sB,"AA MP2 OPDM","RR MP2 OPDM",
    "T AR Amplitudes",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);

  if (params_.print) {
    fprintf(outfile,"ind220             = %18.12lf  H\n",ind220);
    fflush(outfile);
  }

  ind202 = ind22_1(calc_info_.sB,calc_info_.WABB,calc_info_.WASS,
    "T BSBS Amplitudes",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    "BS RI Integrals","SS RI Integrals",calc_info_.evalsB,calc_info_.noccB,
    calc_info_.nvirB);
  ind202 += ind22_2(calc_info_.sB,calc_info_.WABB,calc_info_.WASS,
    "T BS Amplitudes",calc_info_.noccB,calc_info_.nvirB);
  ind202 += ind22_3(calc_info_.sB,calc_info_.WABS,"BB MP2 OPDM","SS MP2 OPDM",
    calc_info_.noccB,calc_info_.nvirB);
  ind202 += ind22_4(calc_info_.sB,"Theta(BS) BS",PSIF_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccB,calc_info_.nvirB);
  ind202 += ind22_5(calc_info_.sB,"T2 BSBS Antisym Amplitudes",
    calc_info_.evalsB,calc_info_.noccB,calc_info_.nvirB);
  ind202 += ind22_6(calc_info_.sB,"T BSBS Antisym Amplitudes",
    PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals","SS RI Integrals",
    calc_info_.noccB,calc_info_.nvirB);
  ind202 += ind22_7(calc_info_.sA,"BB MP2 OPDM","SS MP2 OPDM",
    "T BS Amplitudes",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccA,calc_info_.nvirA);

  if (params_.print) {
    fprintf(outfile,"ind202             = %18.12lf  H\n\n",ind202);
    fflush(outfile);
  }

  results_.ind22 = ind220 + ind202;
}

double SAPT2::ind22_1(double **iAR, double **wAA, double **wRR, 
  const char *t_label, int dfnum, const char *OO_label, const char *OV_label, 
  const char *VV_label, double *evals, int nocc, int nvir)
{
  double energy = 0.0;

  double **xARAR = block_matrix(nocc*nvir,nocc*nvir);
  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,t_label,nocc*nvir,nocc*nvir);

  C_DGEMM('N','N',nocc,nvir*nocc*nvir,nocc,-1.0,&(wAA[0][0]),nocc,
    &(tARAR[0][0]),nvir*nocc*nvir,0.0,&(xARAR[0][0]),nvir*nocc*nvir);

  C_DGEMM('N','N',nocc*nvir*nocc,nvir,nvir,1.0,&(tARAR[0][0]),nvir,
    &(wRR[0][0]),nvir,1.0,&(xARAR[0][0]),nvir);

  free_block(tARAR);

  double **B_p_RR = get_DF_ints(dfnum,VV_label,nvir*nvir);
  double **C_p_AR = block_matrix(nocc*nvir,calc_info_.nrio);

  C_DGEMM('N','N',nocc,nvir*calc_info_.nrio,nvir,1.0,&(iAR[0][0]),nvir,
    &(B_p_RR[0][0]),nvir*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
    nvir*calc_info_.nrio);

  free_block(B_p_RR);

  double **B_p_AR = get_DF_ints(dfnum,OV_label,nocc*nvir);

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
    calc_info_.nrio,&(C_p_AR[0][0]),calc_info_.nrio,1.0,&(xARAR[0][0]),
    nocc*nvir);

  double **B_p_AA = get_DF_ints(dfnum,OO_label,nocc*nocc);

  for(int a=0; a<nocc; a++) {
    C_DGEMM('T','N',nvir,calc_info_.nrio,nocc,1.0,&(iAR[0][0]),nvir,
      &(B_p_AA[a*nocc][0]),calc_info_.nrio,0.0,&(C_p_AR[a*nvir][0]),
      calc_info_.nrio);
  }

  free_block(B_p_AA);

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,calc_info_.nrio,-1.0,&(B_p_AR[0][0]),
    calc_info_.nrio,&(C_p_AR[0][0]),calc_info_.nrio,1.0,&(xARAR[0][0]),
    nocc*nvir);

  free_block(B_p_AR);
  free_block(C_p_AR);

  for(int ar=0; ar<nocc*nvir; ar++) {
    for(int a1r1=0; a1r1<ar; a1r1++) {
      double tval = xARAR[ar][a1r1] + xARAR[a1r1][ar];
      xARAR[a1r1][ar] = tval;
      xARAR[ar][a1r1] = tval;
  }}

  C_DSCAL(nocc*nvir,2.0,&(xARAR[0][0]),nocc*nvir+1);

  for (int a=0, ar=0; a < nocc; a++) {
  for (int r=0; r < nvir; r++, ar++) {
    for (int aa=0, aarr=0; aa < nocc; aa++) {
    for (int rr=0; rr < nvir; rr++, aarr++) {
      int aar = aa*nvir + r;
      int arr = a*nvir + rr;
      double denom = evals[a]+evals[aa]-evals[r+nocc]-evals[rr+nocc];
      double tval = xARAR[ar][aarr];
      energy += (2.0*tval - xARAR[arr][aar])*tval/denom;
    }}
  }}

  free_block(xARAR);

  return(energy);
}

double SAPT2::ind22_2(double **iAR, double **wAA, double **wRR, 
  const char *t_label, int nocc, int nvir)
{
  double energy = 0.0;

  double **zAR = block_matrix(nocc,nvir);
  double **tAR = read_IJKL(PSIF_SAPT_AMPS,t_label,nocc,nvir);

  C_DGEMM('N','N',nocc,nvir,nvir,1.0,&(iAR[0][0]),nvir,&(wRR[0][0]),nvir,0.0,
    &(zAR[0][0]),nvir);

  C_DGEMM('N','N',nocc,nvir,nocc,-1.0,&(wAA[0][0]),nocc,&(iAR[0][0]),nvir,1.0,
    &(zAR[0][0]),nvir);

  energy = C_DDOT(nocc*nvir,&(tAR[0][0]),1,&(zAR[0][0]),1);

  free_block(tAR);
  free_block(zAR);

  return(4.0*energy);
}

double SAPT2::ind22_3(double **iAR, double **wAR, const char *OO_label, 
  const char *VV_label, int nocc, int nvir)
{
  double energy = 0.0;

  double **xAA = read_IJKL(PSIF_SAPT_AMPS,OO_label,nocc,nocc);
  double **yAA = block_matrix(nocc,nocc);

  C_DGEMM('N','T',nocc,nocc,nvir,1.0,&(iAR[0][0]),nvir,&(wAR[0][0]),nvir,0.0,
    &(yAA[0][0]),nocc);

  energy += C_DDOT(nocc*nocc,&(xAA[0][0]),1,&(yAA[0][0]),1);

  free_block(xAA);
  free_block(yAA);

  double **xRR = read_IJKL(PSIF_SAPT_AMPS,VV_label,nvir,nvir);
  double **yRR = block_matrix(nvir,nvir);

  C_DGEMM('T','N',nvir,nvir,nocc,1.0,&(iAR[0][0]),nvir,&(wAR[0][0]),nvir,0.0,
    &(yRR[0][0]),nvir);

  energy += C_DDOT(nvir*nvir,&(xRR[0][0]),1,&(yRR[0][0]),1);

  free_block(xRR);
  free_block(yRR);

  return(-2.0*energy);
}

double SAPT2::ind22_4(double **iAR, const char *theta_OV, int dfnum, 
  const char *OV_label, int nocc, int nvir)
{
  double energy = 0.0;

  double **B_p_AR = get_DF_ints(dfnum,OV_label,nocc*nvir);
  double **C_p_AR = block_matrix(nocc*nvir,calc_info_.nrio);

  double **yAA = block_matrix(nocc,nocc);

  C_DGEMM('N','T',nocc,nocc,nvir,1.0,&(iAR[0][0]),nvir,&(iAR[0][0]),nvir,0.0,
    &(yAA[0][0]),nocc);

  C_DGEMM('N','N',nocc,nvir*calc_info_.nrio,nocc,1.0,&(yAA[0][0]),nocc,
    &(B_p_AR[0][0]),nvir*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
    nvir*calc_info_.nrio);

  free_block(yAA); 

  double **yRR = block_matrix(nvir,nvir);

  C_DGEMM('T','N',nvir,nvir,nocc,1.0,&(iAR[0][0]),nvir,&(iAR[0][0]),nvir,0.0,
    &(yRR[0][0]),nvir);

  for(int a=0; a<nocc; a++) {
    C_DGEMM('N','N',nvir,calc_info_.nrio,nvir,1.0,&(yRR[0][0]),nvir,
      &(B_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(C_p_AR[a*nvir][0]),
      calc_info_.nrio);
  }

  free_block(yRR);
  free_block(B_p_AR);

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,theta_OV,nocc*nvir,
    calc_info_.nrio);

  energy = C_DDOT((long int) nocc*nvir*calc_info_.nrio,C_p_AR[0],1,
    T_p_AR[0],1);

  free_block(C_p_AR);
  free_block(T_p_AR);

  return(-2.0*energy);
}

double SAPT2::ind22_5(double **iAR, const char *t_label, double *evals,
  int nocc, int nvir)
{
  double energy = 0.0;

  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,t_label,nocc*nvir,nocc*nvir);

  for (int a=0, ar=0; a < nocc; a++) {
  for (int r=0; r < nvir; r++, ar++) {
    double tval = iAR[a][r];
    for (int aa=0, aarr=0; aa < nocc; aa++) {
    for (int rr=0; rr < nvir; rr++, aarr++) {
      double denom = evals[a]+evals[aa]-evals[r+nocc]-evals[rr+nocc];
      energy += tval*iAR[aa][rr]*tARAR[ar][aarr]*denom;
    }}
  }}

  free_block(tARAR);

  return(2.0*energy);
}

double SAPT2::ind22_6(double **iAR, const char *t_label, int dfnum, 
  const char *OO_label, const char *OV_label, const char *VV_label, 
  int nocc, int nvir)
{
  double energy = 0.0;

  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,t_label,nocc*nvir,nocc*nvir);
  double **tAR = block_matrix(nocc,nvir);

  C_DGEMV('n',nocc*nvir,nocc*nvir,1.0,&(tARAR[0][0]),nocc*nvir,
    &(iAR[0][0]),1,0.0,&(tAR[0][0]),1);

  free_block(tARAR);

  double *X = init_array(calc_info_.nrio);
  double *Y = init_array(calc_info_.nrio);

  double **B_p_AR = get_DF_ints(dfnum,OV_label,nocc*nvir);

  C_DGEMV('t',nocc*nvir,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(iAR[0][0]),1,0.0,X,1);

  C_DGEMV('t',nocc*nvir,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(tAR[0][0]),1,0.0,Y,1);

  energy += 2.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);
  free_block(B_p_AR);

  double **B_p_RR = get_DF_ints(dfnum,VV_label,nvir*nvir);
  double **C_p_AR = block_matrix(nocc*nvir,calc_info_.nrio);

  C_DGEMM('N','N',nocc,nvir*calc_info_.nrio,nvir,1.0,&(iAR[0][0]),nvir,
    &(B_p_RR[0][0]),nvir*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
    nvir*calc_info_.nrio);

  free_block(B_p_RR);
  double **C_p_AA = block_matrix(nocc*nocc,calc_info_.nrio);

  for(int a=0; a<nocc; a++) {
    C_DGEMM('N','N',nocc,calc_info_.nrio,nvir,1.0,&(tAR[0][0]),nvir,
      &(C_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(C_p_AA[a*nocc][0]),
      calc_info_.nrio);
  }

  free_block(tAR);
  free_block(C_p_AR);
  double **B_p_AA = get_DF_ints(dfnum,OO_label,nocc*nocc);

  energy -= C_DDOT((long int) nocc*nocc*calc_info_.nrio,B_p_AA[0],1,
    C_p_AA[0],1);

  free_block(B_p_AA);
  free_block(C_p_AA);

  return(-4.0*energy);
}

double SAPT2::ind22_7(double **iBS, const char *OOlabel, const char *VVlabel, 
  const char *t_label, int AAnum, const char *AA_label, const char *AR_label, 
  const char *RR_label, int BBnum, const char *BS_label, int noccA, int nvirA,
  int noccB, int nvirB)
{
  double energy = 0.0;

  double *X = init_array(calc_info_.nrio);

  double **B_p_BS = get_DF_ints(BBnum,BS_label,noccB*nvirB);

  C_DGEMV('t',noccB*nvirB,calc_info_.nrio,1.0,&(B_p_BS[0][0]),calc_info_.nrio,
    &(iBS[0][0]),1,0.0,X,1);

  free_block(B_p_BS);

  double *Y = init_array(calc_info_.nrio);

  double **xAA = read_IJKL(PSIF_SAPT_AMPS,OOlabel,noccA,noccA);
  double **B_p_AA = get_DF_ints(AAnum,AA_label,noccA*noccA);

  C_DGEMV('t',noccA*noccA,calc_info_.nrio,-1.0,&(B_p_AA[0][0]),calc_info_.nrio,
    &(xAA[0][0]),1,0.0,Y,1);

  free_block(xAA);
  free_block(B_p_AA);

  double **xRR = read_IJKL(PSIF_SAPT_AMPS,VVlabel,nvirA,nvirA);
  double **B_p_RR = get_DF_ints(AAnum,RR_label,nvirA*nvirA);

  C_DGEMV('t',nvirA*nvirA,calc_info_.nrio,1.0,&(B_p_RR[0][0]),calc_info_.nrio,
    &(xRR[0][0]),1,1.0,Y,1);

  free_block(xRR);
  free_block(B_p_RR);

  double **tAR = read_IJKL(PSIF_SAPT_AMPS,t_label,noccA,nvirA);
  double **B_p_AR = get_DF_ints(AAnum,AR_label,noccA*nvirA);

  C_DGEMV('t',noccA*nvirA,calc_info_.nrio,2.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(tAR[0][0]),1,1.0,Y,1);

  free_block(tAR);
  free_block(B_p_AR);

  energy = C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);

  return(8.0*energy);
}

}}
