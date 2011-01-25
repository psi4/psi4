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
#include "sapt2p.h"

namespace psi { namespace sapt {

void SAPT2p::disp22sdq()
{
  double d211,d220s,d202s,d220d,d202d,d220q,d202q;

  if (params_.print)
    fprintf(outfile,"Begining Disp22(SDQ) Calculation\n\n");

  d211 = disp211();

  if (params_.print) {
    fprintf(outfile,"disp211            = %18.12lf  H\n",d211);
    fflush(outfile);
  }

  d220s = disp22s("T AR Amplitudes","T(BS) AR",PSIF_SAPT_AA_DF_INTS,
    "AA RI Integrals","RR RI Integrals",calc_info_.noccA,calc_info_.nvirA);

  if (params_.print) {
    fprintf(outfile,"disp220s           = %18.12lf  H\n",d220s);
    fflush(outfile);
  }

  d202s = disp22s("T BS Amplitudes","T(AR) BS",PSIF_SAPT_BB_DF_INTS,
    "BB RI Integrals","SS RI Integrals",calc_info_.noccB,calc_info_.nvirB);

  if (params_.print) {
    fprintf(outfile,"disp202s           = %18.12lf  H\n",d202s);
    fflush(outfile);
  }

  d220d = disp220d();

  if (params_.print) {
    fprintf(outfile,"disp220d           = %18.12lf  H\n",d220d);
    fflush(outfile);
  }

  d202d = disp202d();

  if (params_.print) {
    fprintf(outfile,"disp202d           = %18.12lf  H\n",d202d);
    fflush(outfile);
  }

  d220q = disp22q_1("Theta(AR) AR","T(BS) AR","T ARAR Antisym Amplitudes",
    calc_info_.noccA,calc_info_.nvirA);
  d220q += disp22q_2("AA MP2 OPDM","RR MP2 OPDM","T(BS) AR",
    PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.noccA,calc_info_.nvirA);
  d220q += disp22q_3("T ARBS Amplitudes",'N','T',"g ARAR Integrals",
    "T ARAR Antisym Amplitudes",calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  d220q += disp22q_4("T ARAR Amplitudes","g ARAR Integrals",
    "T BSAR Amplitudes",calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  d220q += disp22q_5("T ARAR Amplitudes","g ARAR Integrals",
    "T ARBS Amplitudes",calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);

  if (params_.print) {
    fprintf(outfile,"disp220q           = %18.12lf  H\n",d220q);
    fflush(outfile);
  }

  d202q = disp22q_1("Theta(BS) BS","T(AR) BS","T BSBS Antisym Amplitudes",
    calc_info_.noccB,calc_info_.nvirB);
  d202q += disp22q_2("BB MP2 OPDM","SS MP2 OPDM","T(AR) BS",
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",calc_info_.noccB,calc_info_.nvirB);
  d202q += disp22q_3("T ARBS Amplitudes",'T','N',"g BSBS Integrals",
    "T BSBS Antisym Amplitudes",calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccA,calc_info_.nvirA);
  d202q += disp22q_4("T BSBS Amplitudes","g BSBS Integrals",
    "T ARBS Amplitudes",calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccA,calc_info_.nvirA);
  d202q += disp22q_5("T BSBS Amplitudes","g BSBS Integrals",
    "T BSAR Amplitudes",calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccA,calc_info_.nvirA);

  if (params_.print) {
    fprintf(outfile,"disp202q           = %18.12lf  H\n\n",d202q);
    fflush(outfile);
  }

  results_.disp22sdq = d211 + d220s + d202s + d220d + d202d + d220q + d202q;
}

double SAPT2p::disp211()
{
  double energy = 0.0;

  double **pARBS = read_IJKL(PSIF_SAPT_AMPS,"GxT(BS) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
  double **qARBS = read_IJKL(PSIF_SAPT_AMPS,"GxT(AR) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  double **B_p_AR = get_AR_ints(0);
  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"Theta(BS) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(T_p_BS[0][0]),calc_info_.nrio,1.0,&(pARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(B_p_AR);

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"Theta(AR) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **B_p_BS = get_BS_ints(0);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(T_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,1.0,&(qARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(B_p_BS);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int b=0, bs=0; b < calc_info_.noccB; b++) {
    for (int s=0; s < calc_info_.nvirB; s++, bs++) {
      double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
        calc_info_.evalsA[r+calc_info_.noccA]-
        calc_info_.evalsB[s+calc_info_.noccB];
      energy += 8.0*pARBS[ar][bs]*qARBS[ar][bs]/denom;
    }}
  }}

  free_block(qARBS);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(T_p_AR[0][0]),calc_info_.nrio,
    &(T_p_BS[0][0]),calc_info_.nrio,0.0,&(pARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(T_p_AR);
  free_block(T_p_BS);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  energy += 8.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.noccB*
    calc_info_.nvirB,tARBS[0],1,pARBS[0],1);

  free_block(pARBS);
  free_block(tARBS);

  return(energy);
}

double SAPT2p::disp22s(char *t_label, char *theta_label, int AAnum,
  char *AA_label, char *RR_label, int nocc, int nvir)
{
  double energy = 0.0;

  double **xAR = read_IJKL(PSIF_SAPT_AMPS,t_label,nocc,nvir);
  double **yAR = block_matrix(nocc,nvir);

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,theta_label,nocc*nvir,
    calc_info_.nrio);
  double **B_p_RR = get_DF_ints(AAnum,RR_label,nvir*nvir);

  C_DGEMM('N','T',nocc,nvir,nvir*calc_info_.nrio,1.0,&(T_p_AR[0][0]),
    nvir*calc_info_.nrio,&(B_p_RR[0][0]),nvir*calc_info_.nrio,0.0,&(yAR[0][0]),
    nvir);

  free_block(B_p_RR);

  double **B_p_AA = get_DF_ints(AAnum,AA_label,nocc*nocc);

  for(int a=0; a<nocc; a++) {
    C_DGEMM('N','T',nocc,nvir,calc_info_.nrio,-1.0,&(B_p_AA[a*nocc][0]),
      calc_info_.nrio,&(T_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(yAR[0][0]),
      nvir);
  }

  free_block(T_p_AR);
  free_block(B_p_AA);

  energy = 8.0*C_DDOT(nocc*nvir,xAR[0],1,yAR[0],1);

  free_block(xAR);
  free_block(yAR);

  return(energy);
}

double SAPT2p::disp220d()
{
  double energy = 0.0;

  double **zARAR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirA);
  double **B_p_AR = get_AR_ints(0);
  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(T_p_AR[0][0]),calc_info_.nrio,0.0,&(zARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA);

  free_block(B_p_AR);
  free_block(T_p_AR);

  for(int ar=0; ar<calc_info_.noccA*calc_info_.nvirA; ar++) {
    for(int a1r1=0; a1r1<ar; a1r1++) {
      double tval = zARAR[ar][a1r1] + zARAR[a1r1][ar];
      zARAR[a1r1][ar] = tval;
      zARAR[ar][a1r1] = tval;
  }}

  C_DSCAL(calc_info_.noccA*calc_info_.nvirA,2.0,&(zARAR[0][0]),
    calc_info_.noccA*calc_info_.nvirA+1);

  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,"T2 ARAR Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int aa=0, aarr=0; aa < calc_info_.noccA; aa++) {
    for (int rr=0; rr < calc_info_.nvirA; rr++, aarr++) {
      int aar = aa*calc_info_.nvirA + r;
      int arr = a*calc_info_.nvirA + rr;
      double tval = 2.0*zARAR[ar][aarr] - zARAR[aar][arr];
      energy += 4.0*tval*tARAR[ar][aarr];
    }}
  }}

  free_block(tARAR);
  free_block(zARAR);

  double **dARBS = read_IJKL(PSIF_SAPT_AMPS,"GxT(AR) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"Theta(AR) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **B_p_BS = get_BS_ints(0);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(T_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,1.0,&(dARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int b=0, bs=0; b < calc_info_.noccB; b++) {
    for (int s=0; s < calc_info_.nvirB; s++, bs++) {
      double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
        calc_info_.evalsA[r+calc_info_.noccA]-
        calc_info_.evalsB[s+calc_info_.noccB];
      double tval = dARBS[ar][bs];
      energy += 4.0*tval*tval/denom;
    }}
  }}

  free_block(dARBS);
  free_block(T_p_AR);
  free_block(B_p_BS);

  return(energy);
}

double SAPT2p::disp202d()
{
  double energy = 0.0;

  double **zBSBS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.nvirB);
  double **B_p_BS = get_BS_ints(0);
  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_BS[0][0]),calc_info_.nrio,
    &(T_p_BS[0][0]),calc_info_.nrio,0.0,&(zBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(B_p_BS);
  free_block(T_p_BS);

  for(int bs=0; bs<calc_info_.noccB*calc_info_.nvirB; bs++) {
    for(int b1s1=0; b1s1<bs; b1s1++) {
      double tval = zBSBS[bs][b1s1] + zBSBS[b1s1][bs];
      zBSBS[b1s1][bs] = tval;
      zBSBS[bs][b1s1] = tval;
  }}

  C_DSCAL(calc_info_.noccB*calc_info_.nvirB,2.0,&(zBSBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB+1);

  double **tBSBS = read_IJKL(PSIF_SAPT_AMPS,"T2 BSBS Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);

  for (int b=0, bs=0; b < calc_info_.noccB; b++) {
  for (int s=0; s < calc_info_.nvirB; s++, bs++) {
    for (int bb=0, bbss=0; bb < calc_info_.noccB; bb++) {
    for (int ss=0; ss < calc_info_.nvirB; ss++, bbss++) {
      int bbs = bb*calc_info_.nvirB + s;
      int bss = b*calc_info_.nvirB + ss;
      double tval = 2.0*zBSBS[bs][bbss] - zBSBS[bbs][bss];
      energy += 4.0*tval*tBSBS[bs][bbss];
    }}
  }}

  free_block(tBSBS);
  free_block(zBSBS);

  double **dARBS = read_IJKL(PSIF_SAPT_AMPS,"GxT(BS) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"Theta(BS) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);
  double **B_p_AR = get_AR_ints(0);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(T_p_BS[0][0]),calc_info_.nrio,1.0,&(dARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int b=0, bs=0; b < calc_info_.noccB; b++) {
    for (int s=0; s < calc_info_.nvirB; s++, bs++) {
      double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
        calc_info_.evalsA[r+calc_info_.noccA]-
        calc_info_.evalsB[s+calc_info_.noccB];
      double tval = dARBS[ar][bs];
      energy += 4.0*tval*tval/denom;
    }}
  }}

  free_block(dARBS);
  free_block(T_p_BS);
  free_block(B_p_AR);

  return(energy);
}

double SAPT2p::disp22q_1(char *label1, char *label2, char *label3, int nocc,
  int nvir)
{
  double energy = 0.0;

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,label1,nocc*nvir,calc_info_.nrio);
  double **T_q_AR = read_IJKL(PSIF_SAPT_AMPS,label2,nocc*nvir,calc_info_.nrio);
  double **xARAR = block_matrix(nocc*nvir,nocc*nvir);

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,calc_info_.nrio,1.0,&(T_p_AR[0][0]),
    calc_info_.nrio,&(T_q_AR[0][0]),calc_info_.nrio,1.0,&(xARAR[0][0]),
    nocc*nvir);

  free_block(T_p_AR);
  free_block(T_q_AR);

  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,label3,nocc*nvir,nocc*nvir);

  energy = 4.0*C_DDOT(nocc*nvir*nocc*nvir,xARAR[0],1,tARAR[0],1);

  free_block(xARAR);
  free_block(tARAR);

  return(energy);
}

double SAPT2p::disp22q_2(char *OO_opdm, char *VV_opdm, char *T_label,
  int DFnum, char *OV_label, int nocc, int nvir)
{
  double energy = 0.0;

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,T_label,nocc*nvir,
    calc_info_.nrio);
  double **B_p_AR = get_DF_ints(DFnum,OV_label,nocc*nvir);

  double **xAA = read_IJKL(PSIF_SAPT_AMPS,OO_opdm,nocc,nocc);
  double **yAA = block_matrix(nocc,nocc);

  C_DGEMM('N','T',nocc,nocc,nvir*calc_info_.nrio,1.0,&(T_p_AR[0][0]),
    nvir*calc_info_.nrio,&(B_p_AR[0][0]),nvir*calc_info_.nrio,0.0,&(yAA[0][0]),
    nocc);

  energy -= 4.0*C_DDOT(nocc*nocc,xAA[0],1,yAA[0],1);

  free_block(xAA);
  free_block(yAA);

  double **xRR = read_IJKL(PSIF_SAPT_AMPS,VV_opdm,nvir,nvir);
  double **yRR = block_matrix(nvir,nvir);

  for (int a=0; a < nocc; a++) {
    C_DGEMM('N','T',nvir,nvir,calc_info_.nrio,1.0,&(T_p_AR[a*nvir][0]),
      calc_info_.nrio,&(B_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(yRR[0][0]),
      nvir);
  }

  energy -= 4.0*C_DDOT(nvir*nvir,xRR[0],1,yRR[0],1);

  free_block(xRR);
  free_block(yRR);
  free_block(T_p_AR);
  free_block(B_p_AR);

  return(energy);
}

double SAPT2p::disp22q_3(char *t_label, char trans1, char trans2,
  char *g_label, char *th_label, int noccA, int nvirA, int noccB, int nvirB)
{
  double energy = 0.0;

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,t_label,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
  double **tARAR = block_matrix(noccA*nvirA,noccA*nvirA);

  C_DGEMM(trans1,trans2,noccA*nvirA,noccA*nvirA,noccB*nvirB,1.0,&(tARBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB,&(tARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,0.0,&(tARAR[0][0]),noccA*nvirA);

  free_block(tARBS);

  double **gARAR = read_IJKL(PSIF_SAPT_AMPS,g_label,noccA*nvirA,noccA*nvirA);
  double **xARAR = block_matrix(noccA*nvirA,noccA*nvirA);

  C_DGEMM('N','N',noccA*nvirA,noccA*nvirA,noccA*nvirA,1.0,&(tARAR[0][0]),
    noccA*nvirA,&(gARAR[0][0]),noccA*nvirA,0.0,&(xARAR[0][0]),noccA*nvirA);

  free_block(tARAR);
  free_block(gARAR);

  double **yARAR = read_IJKL(PSIF_SAPT_AMPS,th_label,noccA*nvirA,noccA*nvirA);

  energy = 4.0*C_DDOT(noccA*nvirA*noccA*nvirA,xARAR[0],1,yARAR[0],1);

  free_block(xARAR);
  free_block(yARAR);

  return(energy);
}

double SAPT2p::disp22q_4(char *t_label, char *g_label, char *d_label,
    int noccA, int nvirA, int noccB, int nvirB)
{
  double energy = 0.0;

  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,t_label,noccA*nvirA,noccA*nvirA);
  double **gARAR = read_IJKL(PSIF_SAPT_AMPS,g_label,noccA*nvirA,noccA*nvirA);
  double **xRR = block_matrix(nvirA,nvirA);

  C_DGEMM('T','N',nvirA,nvirA,noccA*nvirA*noccA,1.0,&(tARAR[0][0]),nvirA,
    &(gARAR[0][0]),nvirA,0.0,&(xRR[0][0]),nvirA);

  free_block(tARAR);
  free_block(gARAR);

  double **tBSAR = read_IJKL(PSIF_SAPT_AMPS,d_label,noccB*nvirB,noccA*nvirA);
  double **yRR = block_matrix(nvirA,nvirA);

  C_DGEMM('T','N',nvirA,nvirA,noccB*nvirB*noccA,1.0,&(tBSAR[0][0]),nvirA,
    &(tBSAR[0][0]),nvirA,0.0,&(yRR[0][0]),nvirA);

  energy = -4.0*C_DDOT(nvirA*nvirA,xRR[0],1,yRR[0],1);

  free_block(tBSAR);
  free_block(xRR);
  free_block(yRR);

  return(energy);
}

double SAPT2p::disp22q_5(char *t_label, char *g_label, char *d_label,
    int noccA, int nvirA, int noccB, int nvirB)
{
  double energy = 0.0;

  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,t_label,noccA*nvirA,noccA*nvirA);
  double **gARAR = read_IJKL(PSIF_SAPT_AMPS,g_label,noccA*nvirA,noccA*nvirA);
  double **xAA = block_matrix(noccA,noccA);

  C_DGEMM('N','T',noccA,noccA,nvirA*noccA*nvirA,1.0,&(tARAR[0][0]),
    nvirA*noccA*nvirA,&(gARAR[0][0]),nvirA*noccA*nvirA,0.0,&(xAA[0][0]),noccA);

  free_block(tARAR);
  free_block(gARAR);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,d_label,noccA*nvirA,noccB*nvirB);
  double **yAA = block_matrix(noccA,noccA);

  C_DGEMM('N','T',noccA,noccA,nvirA*noccB*nvirB,1.0,&(tARBS[0][0]),
    nvirA*noccB*nvirB,&(tARBS[0][0]),nvirA*noccB*nvirB,0.0,&(yAA[0][0]),noccA);

  energy = -4.0*C_DDOT(noccA*noccA,xAA[0],1,yAA[0],1);

  free_block(tARBS);
  free_block(xAA);
  free_block(yAA);

  return(energy);
}

}}
