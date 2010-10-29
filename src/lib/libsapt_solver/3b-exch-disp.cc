/* This function calculates the Exch10 energy */

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
#include "structs.h"
#include "sapt3bn5.h"

namespace psi { namespace sapt {

void SAPT3BN5::exch_disp200_s2()
{ 
  results_.exch_disp200_s2 = exch_disp200_s2_0(calc_info_.S_AC,calc_info_.S_BC,
    calc_info_.WABS,calc_info_.WBAR,"T(BS) AR Intermediates",
    "T(AR) BS Intermediates",PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals","T2 ARBS Amplitudes",'N',
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC);
  fprintf(outfile,"exch_disp_200 S^2   = %18.12lf  H\n",
    results_.exch_disp200_s2);
  fflush(outfile);
  results_.exch_disp020_s2 = exch_disp200_s2_0(calc_info_.S_CB,calc_info_.S_AB,
    calc_info_.WCAR,calc_info_.WACT,"T(AR) CT Intermediates",
    "T(CT) AR Intermediates",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals","T2 ARCT Amplitudes",'T',
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB);
  fprintf(outfile,"exch_disp_020 S^2   = %18.12lf  H\n",
    results_.exch_disp020_s2);
  fflush(outfile);
  results_.exch_disp002_s2 = exch_disp200_s2_0(calc_info_.S_BA,calc_info_.S_CA,
    calc_info_.WBCT,calc_info_.WCBS,"T(CT) BS Intermediates",
    "T(BS) CT Intermediates",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals","T2 BSCT Amplitudes",'N',
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA);
  fprintf(outfile,"exch_disp_002 S^2   = %18.12lf  H\n\n",
    results_.exch_disp002_s2);
  fflush(outfile);
}

void SAPT3BN5::exch_disp110_s2()
{ 
  results_.exch_disp110_s2 = exch_disp110_s2_0(calc_info_.S_CA,calc_info_.S_CB,
    calc_info_.S_AB,calc_info_.S_AC,calc_info_.S_BA,calc_info_.WCBS,
    calc_info_.WBCT,"T2 ARCT Amplitudes",'N',"T2 ARBS Amplitudes",'N',
    "T(CT) AR Intermediates","T(BS) AR Intermediates",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"exch_disp_110 S^2   = %18.12lf  H\n",
    results_.exch_disp110_s2);
  fflush(outfile);
  results_.exch_disp101_s2 = exch_disp110_s2_0(calc_info_.S_AB,calc_info_.S_AC,
    calc_info_.S_BC,calc_info_.S_BA,calc_info_.S_CB,calc_info_.WACT,
    calc_info_.WCAR,"T2 ARBS Amplitudes",'T',"T2 BSCT Amplitudes",'N',
    "T(AR) BS Intermediates","T(CT) BS Intermediates",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"exch_disp_101 S^2   = %18.12lf  H\n",
    results_.exch_disp101_s2); 
  fflush(outfile);
  results_.exch_disp011_s2 = exch_disp110_s2_0(calc_info_.S_BC,
    calc_info_.S_BA,calc_info_.S_CA,calc_info_.S_CB,calc_info_.S_AC,
    calc_info_.WBAR,calc_info_.WABS,"T2 BSCT Amplitudes",'T',
    "T2 ARCT Amplitudes",'T',"T(BS) CT Intermediates","T(AR) CT Intermediates",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"exch_disp_011 S^2   = %18.12lf  H\n\n",
    results_.exch_disp011_s2); 
  fflush(outfile);
}

void SAPT3BN5::exch_disp200_s3()
{
  results_.exch_disp200_s3 = exch_disp200_s3_0(calc_info_.S_AB,calc_info_.S_BC,
    calc_info_.S_CA,calc_info_.S_AC,calc_info_.S_BA,calc_info_.S_CB,
    calc_info_.WABS,calc_info_.WBAR,"T(AR) BS Intermediates",
    "T(BS) AR Intermediates",PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals","T2 ARBS Amplitudes",'n','t',
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"exch_disp_200 S^3   = %18.12lf  H\n",
    results_.exch_disp200_s3);
  fflush(outfile);
  results_.exch_disp020_s3 = exch_disp200_s3_0(calc_info_.S_CA,calc_info_.S_AB,
    calc_info_.S_BC,calc_info_.S_CB,calc_info_.S_AC,calc_info_.S_BA,
    calc_info_.WCAR,calc_info_.WACT,"T(CT) AR Intermediates",
    "T(AR) CT Intermediates",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals","T2 ARCT Amplitudes",'t','n',
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"exch_disp_020 S^3   = %18.12lf  H\n",
    results_.exch_disp020_s3);
  fflush(outfile);
  results_.exch_disp002_s3 = exch_disp200_s3_0(calc_info_.S_BC,calc_info_.S_CA,
    calc_info_.S_AB,calc_info_.S_BA,calc_info_.S_CB,calc_info_.S_AC,
    calc_info_.WBCT,calc_info_.WCBS,"T(BS) CT Intermediates",
    "T(CT) BS Intermediates",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals","T2 BSCT Amplitudes",'n','t',
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"exch_disp_002 S^3   = %18.12lf  H\n\n",
    results_.exch_disp002_s3);
  fflush(outfile);
}

void SAPT3BN5::exch_disp110_s3()
{
  results_.exch_disp110_s3 = exch_disp110_s3_0(calc_info_.S_CA,calc_info_.S_AB,
    calc_info_.S_BC,calc_info_.S_CB,calc_info_.S_AC,calc_info_.S_BA,
    calc_info_.WCAR,calc_info_.WBAR,calc_info_.WACT,calc_info_.WABS,
    "T(AR) BS Intermediates","T(AR) CT Intermediates",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals","T2 ARBS Amplitudes",'n',
    "T2 ARCT Amplitudes",'n',calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"exch_disp_110 S^3   = %18.12lf  H\n",
    results_.exch_disp110_s3);
  fflush(outfile);
  results_.exch_disp101_s3 = exch_disp110_s3_0(calc_info_.S_AB,calc_info_.S_BC,
    calc_info_.S_CA,calc_info_.S_AC,calc_info_.S_BA,calc_info_.S_CB,
    calc_info_.WABS,calc_info_.WCBS,calc_info_.WBAR,calc_info_.WBCT,
    "T(BS) CT Intermediates","T(BS) AR Intermediates",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals","T2 BSCT Amplitudes",'n',
    "T2 ARBS Amplitudes",'t',calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"exch_disp_101 S^3   = %18.12lf  H\n",
    results_.exch_disp101_s3);
  fflush(outfile);
  results_.exch_disp011_s3 = exch_disp110_s3_0(calc_info_.S_BC,calc_info_.S_CA,
    calc_info_.S_AB,calc_info_.S_BA,calc_info_.S_CB,calc_info_.S_AC,
    calc_info_.WBCT,calc_info_.WACT,calc_info_.WCBS,calc_info_.WCAR,
    "T(CT) AR Intermediates","T(CT) BS Intermediates",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals","T2 ARCT Amplitudes",'t',
    "T2 BSCT Amplitudes",'t',calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"exch_disp_011 S^3   = %18.12lf  H\n\n",
    results_.exch_disp011_s3);
  fflush(outfile);
}

double SAPT3BN5::exch_disp200_s2_0(double **SAC, double **SBC, double **WABS,
  double **WBAR, char *TBS_AR, char *TAR_BS, int AAfile, char *AR_ints,
  int BBfile, char *BS_ints, char *T2label, char trans, int occA, int virA,
  int occB, int virB, int occC)
{
  double energy = 0.0;
  energy += exch_disp200_s2_1(SAC,TBS_AR,AAfile,AR_ints,occA,virA,occC);
  energy += exch_disp200_s2_1(SBC,TAR_BS,BBfile,BS_ints,occB,virB,occC);
  energy += exch_disp200_s2_2(SAC,SBC,WABS,WBAR,T2label,trans,occA,virA,
    occB,virB,occC);
  return(energy);
}

double SAPT3BN5::exch_disp200_s2_1(double **SAC, char *TBS_AR, int AAfile, 
  char *AR_ints, int occA, int virA, int occC)
{
  double energy = 0.0;

  double **X_AA = block_matrix(occA,occA);
  double **X_RR = block_matrix(virA,virA);

  C_DGEMM('N','T',occA,occA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(X_AA[0][0]),occA);

  C_DGEMM('N','T',virA,virA,occC,1.0,&(SAC[occA][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(X_RR[0][0]),virA);

  double **B_p_AR = block_matrix(occA*virA,calc_info_.nrio);
  double **T_AR = read_IJKL(PSIF_3B_SAPT_AMPS,TBS_AR,occA*virA,
    calc_info_.nrio);

  C_DGEMM('N','N',occA,virA*calc_info_.nrio,occA,1.0,&(X_AA[0][0]),
          occA,&(T_AR[0][0]),virA*calc_info_.nrio,0.0,&(B_p_AR[0][0]),
          virA*calc_info_.nrio);

  for (int i=0; i<occA; i++) {
    C_DGEMM('T','N',virA,calc_info_.nrio,virA,-1.0,&(X_RR[0][0]),
            virA,&(T_AR[i*virA][0]),calc_info_.nrio,1.0,&(B_p_AR[i*virA][0]),
            calc_info_.nrio);
  }

  free_block(X_AA);
  free_block(X_RR);
  free_block(T_AR);

  double **B_q_AR = get_DF_ints(AAfile,AR_ints,occA*virA);

  energy = 4.0*C_DDOT(occA*virA*calc_info_.nrio,&(B_p_AR[0][0]),1,
                      &(B_q_AR[0][0]),1);

  free_block(B_p_AR);
  free_block(B_q_AR);

  return(energy);
}

double SAPT3BN5::exch_disp200_s2_2(double **SAC, double **SBC, double **WABS,
  double **WBAR, char *T2label, char trans, int occA, int virA, int occB, 
  int virB, int occC)
{
  double energy = 0.0;

  double **X_AR = block_matrix(occA,virA);
  double **X_BS = block_matrix(occB,virB);
  double **Y_AR = block_matrix(occA,virA);
  double **Y_BS = block_matrix(occB,virB);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(X_AR[0][0]),virA);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  if (trans == 'N') {
    double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virA*occA,occB*virB);
    C_DGEMV('n',occA*virA,occB*virB,1.0,&(tARBS[0][0]),occB*virB,
            &(WABS[0][0]),1,0.0,&(Y_AR[0][0]),1);
    C_DGEMV('t',occA*virA,occB*virB,1.0,&(tARBS[0][0]),occB*virB,
            &(WBAR[0][0]),1,0.0,&(Y_BS[0][0]),1);
    energy -= 4.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(Y_AR[0][0]),1);
    energy -= 4.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(Y_BS[0][0]),1);
    free_block(tARBS);
  }
  else {
    double **tBSAR = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occA*virA);
    C_DGEMV('t',occB*virB,occA*virA,1.0,&(tBSAR[0][0]),occA*virA,
            &(WABS[0][0]),1,0.0,&(Y_AR[0][0]),1);
    C_DGEMV('n',occB*virB,occA*virA,1.0,&(tBSAR[0][0]),occA*virA,
            &(WBAR[0][0]),1,0.0,&(Y_BS[0][0]),1);
    energy -= 4.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(Y_AR[0][0]),1);
    energy -= 4.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(Y_BS[0][0]),1);
    free_block(tBSAR);
  }

  free_block(X_AR);
  free_block(X_BS);
  free_block(Y_AR);
  free_block(Y_BS);
  return(energy);
}

double SAPT3BN5::exch_disp110_s2_0(double **SAB, double **SAC, double **SBC, 
  double **SBA, double **SCB, double **WACT, double **WCAR, char *T2ARBS, 
  char trans1, char *T2BSCT, char trans2, char *TAR_BS, char *TCT_BS, 
  int AAfile, char *AR_ints, int CCfile, char *CT_ints, int occA, int virA, 
  int occB, int virB, int occC, int virC)
{
  double energy = 0.0;
  energy += exch_disp110_s2_1(SAB,TCT_BS,AAfile,AR_ints,occA,virA,occB,virB);
  energy += exch_disp110_s2_1(SCB,TAR_BS,CCfile,CT_ints,occC,virC,occB,virB);
  energy += exch_disp110_s2_2(SAB,SBC,WACT,T2BSCT,trans2,occA,virA,occB,virB,
    occC,virC);
  energy += exch_disp110_s2_2(SCB,SBA,WCAR,T2ARBS,trans1,occC,virC,occB,virB,
    occA,virA);
  return(energy);
}

double SAPT3BN5::exch_disp110_s2_1(double **SAB, char *TCT_BS, int AAfile, 
  char *AR_ints, int occA, int virA, int occB, int virB)
{
  double energy = 0.0;

  double **B_p_AS = block_matrix(occA*virB,calc_info_.nrio);
  double **TCT_B = read_IJKL(PSIF_3B_SAPT_AMPS,TCT_BS,occB*virB,
    calc_info_.nrio);

  C_DGEMM('N','N',occA,virB*calc_info_.nrio,occB,1.0,&(SAB[0][0]),
    calc_info_.nmo,&(TCT_B[0][0]),virB*calc_info_.nrio,0.0,&(B_p_AS[0][0]),
    virB*calc_info_.nrio);

  free_block(TCT_B);

  double **B_p_AR = block_matrix(occA*virA,calc_info_.nrio);

  for (int i=0; i<occA; i++) {
    C_DGEMM('N','N',virA,calc_info_.nrio,virB,1.0,&(SAB[occA][occB]),
            calc_info_.nmo,&(B_p_AS[i*virB][0]),calc_info_.nrio,0.0,
            &(B_p_AR[i*virA][0]),calc_info_.nrio);
  }

  free_block(B_p_AS);
  double **B_q_AR = get_DF_ints(AAfile,AR_ints,occA*virA);

  energy -= 4.0*C_DDOT(occA*virA*calc_info_.nrio,&(B_p_AR[0][0]),1,
                       &(B_q_AR[0][0]),1);

  free_block(B_p_AR);
  free_block(B_q_AR);

  return(energy);
}

double SAPT3BN5::exch_disp110_s2_2(double **SAB, double **SBC, double **WACT,
  char *T2label, char trans, int occA, int virA, int occB, int virB, 
  int occC, int virC)
{
  double energy = 0.0;

  double **X_BS = block_matrix(occB,virB);
  double **Y_BS = block_matrix(occB,virB);

  C_DGEMM('T','N',occB,virB,occA,1.0,&(SAB[0][0]),calc_info_.nmo,
    &(SAB[0][occB]),calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  double **X_BT = block_matrix(occB,virC);
  double **X_CS = block_matrix(occC,virB);
  double **Y_BT = block_matrix(occB,virC);
  double **Y_CS = block_matrix(occC,virB);

  C_DGEMM('N','T',occC,virB,virC,1.0,&(WACT[0][0]),virC,
    &(SBC[occB][occC]),calc_info_.nmo,0.0,&(X_CS[0][0]),virB);

  C_DGEMM('N','N',occB,virC,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(WACT[0][0]),virC,0.0,&(X_BT[0][0]),virC);

  if (trans == 'N') {
    double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occC*virC);
    C_DGEMV('n',occB*virB,occC*virC,1.0,&(tBSCT[0][0]),occC*virC,
            &(WACT[0][0]),1,0.0,&(Y_BS[0][0]),1);
    for (int b=0;b<occB;b++) {
      for (int s=0;s<virB;s++) {
        for (int c=0;c<occC;c++) {
          for (int t=0;t<virC;t++) {
            int bs = b*virB + s;
            int ct = c*virC + t;
            double tval = tBSCT[bs][ct];
            Y_BT[b][t] += tval*SBC[s+occB][c];
            Y_CS[c][s] += tval*SBC[b][t+occC];
    }}}}
    free_block(tBSCT);
  }
  else {
    double **tCTBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virC*occC,occB*virB);
    C_DGEMV('t',occC*virC,occB*virB,1.0,&(tCTBS[0][0]),occB*virB,
            &(WACT[0][0]),1,0.0,&(Y_BS[0][0]),1);
    for (int c=0;c<occC;c++) {
      for (int t=0;t<virC;t++) {
        for (int b=0;b<occB;b++) {
          for (int s=0;s<virB;s++) {
            int bs = b*virB + s;
            int ct = c*virC + t;
            double tval = tCTBS[ct][bs];
            Y_BT[b][t] += tval*SBC[s+occB][c];
            Y_CS[c][s] += tval*SBC[b][t+occC];
    }}}}
    free_block(tCTBS);
  }

  energy -= 4.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(Y_BS[0][0]),1);
  energy -= 2.0*C_DDOT(occC*virB,&(X_CS[0][0]),1,&(Y_CS[0][0]),1);
  energy += 2.0*C_DDOT(occB*virC,&(X_BT[0][0]),1,&(Y_BT[0][0]),1);

  free_block(X_BS);
  free_block(Y_BS);
  free_block(X_BT);
  free_block(X_CS);
  free_block(Y_BT);
  free_block(Y_CS);

  return(energy);
}

double SAPT3BN5::exch_disp200_s3_0(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, double **WABS, double **WBAR,
  char *TAR_BS, char *TBS_AR, int AAfile, char *AR_ints, int BBfile,
  char *BS_ints, char *T2label, char trans1, char trans2, int occA, int virA,
  int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  energy += exch_disp200_s3_1(SAB,SBC,SCA,SAC,SBA,SCB,TAR_BS,BBfile,BS_ints,
    occA,virA,occB,virB,occC,virC);
  energy += exch_disp200_s3_1(SBA,SAC,SCB,SBC,SAB,SCA,TBS_AR,AAfile,AR_ints,
    occB,virB,occA,virA,occC,virC);
  energy += exch_disp200_s3_2(SAB,SBC,SCA,SAC,SBA,SCB,AAfile,AR_ints,
    BBfile,BS_ints,T2label,trans1,occA,virA,occB,virB,occC,virC);
  energy += exch_disp200_s3_3(SAB,SBC,SCA,SAC,SBA,SCB,AAfile,AR_ints,BBfile,
    BS_ints,T2label,trans1,occA,virA,occB,virB,occC,virC);
  energy += exch_disp200_s3_3(SBA,SAC,SCB,SBC,SAB,SCA,BBfile,BS_ints,AAfile,
    AR_ints,T2label,trans2,occB,virB,occA,virA,occC,virC);
  energy += exch_disp200_s3_4(SAB,SBC,SCA,SAC,SBA,SCB,WABS,T2label,trans1,occA,
    virA,occB,virB,occC,virC);
  energy += exch_disp200_s3_4(SBA,SAC,SCB,SBC,SAB,SCA,WBAR,T2label,trans2,occB,
    virB,occA,virA,occC,virC);

  return(energy);
}

double SAPT3BN5::exch_disp200_s3_1(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, char *TAR_BS, int BBfile,
  char *BS_ints, int occA, int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **T_BS = read_IJKL(PSIF_3B_SAPT_AMPS,TAR_BS,occB*virB,
    calc_info_.nrio);
  double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);

  double **X_BB = block_matrix(occB,occB);
  double **X_SS = block_matrix(virB,virB);

  C_DGEMM('N','T',occB,occB,virB*calc_info_.nrio,1.0,&(B_p_BS[0][0]),
    virB*calc_info_.nrio,&(T_BS[0][0]),virB*calc_info_.nrio,0.0,&(X_BB[0][0]),
    occB);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','T',virB,virB,calc_info_.nrio,1.0,&(B_p_BS[b*virB][0]),
      calc_info_.nrio,&(T_BS[b*virB][0]),calc_info_.nrio,1.0,&(X_SS[0][0]),
      virB);
  }

  free_block(T_BS);
  free_block(B_p_BS);

  double **Y_SS = block_matrix(virB,virB);
  double **Y_SC = block_matrix(virB,occC);
  double **Y_SA = block_matrix(virB,occA);

  double **Y_BB = block_matrix(occB,occB);
  double **Y_BA = block_matrix(occB,occA);
  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',virB,occC,occA,1.0,&(SBA[occB][0]),calc_info_.nmo,
    &(SAC[0][0]),calc_info_.nmo,0.0,&(Y_SC[0][0]),occC);

  C_DGEMM('N','N',virB,occA,occC,1.0,&(SBC[occB][0]),calc_info_.nmo,
    &(SCA[0][0]),calc_info_.nmo,0.0,&(Y_SA[0][0]),occA);

  C_DGEMM('N','N',occB,occA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SCA[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  C_DGEMM('N','N',occB,occC,occA,1.0,&(SBA[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  C_DGEMM('N','N',virB,virB,occC,1.0,&(Y_SC[0][0]),occC,&(SCB[0][occB]),
    calc_info_.nmo,0.0,&(Y_SS[0][0]),virB);

  C_DGEMM('N','N',virB,virB,occA,1.0,&(Y_SA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,1.0,&(Y_SS[0][0]),virB);

  C_DGEMM('N','N',occB,occB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  C_DGEMM('N','N',occB,occB,occC,1.0,&(Y_BC[0][0]),occC,&(SCB[0][0]),
    calc_info_.nmo,1.0,&(Y_BB[0][0]),occB);

  free_block(Y_SC);
  free_block(Y_SA);
  free_block(Y_BA);
  free_block(Y_BC);

  energy += 4.0*C_DDOT(virB*virB,&(X_SS[0][0]),1,&(Y_SS[0][0]),1);
  energy -= 4.0*C_DDOT(occB*occB,&(X_BB[0][0]),1,&(Y_BB[0][0]),1);

  free_block(X_BB);
  free_block(Y_BB);
  free_block(X_SS);
  free_block(Y_SS);

  return(energy);
}

double SAPT3BN5::exch_disp200_s3_2(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, int AAfile, char *AR_ints,
  int BBfile, char *BS_ints, char *T2label, char trans, int occA, int virA,
  int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **X_AB = block_matrix(occA,occB);
  double **Y_AB = block_matrix(occA,occB);

  for (int a=0; a<occA; a++) {
    C_DCOPY(occB,&(SAB[a][0]),1,&(X_AB[a][0]),1);
  }

  C_DGEMM('N','N',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SCB[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
  double **C_p_AS = block_matrix(occA*virB,calc_info_.nrio);
  double **D_p_AS = block_matrix(occA*virB,calc_info_.nrio);

  C_DGEMM('N','N',occA,virB*calc_info_.nrio,occB,1.0,&(X_AB[0][0]),occB,
    &(B_p_BS[0][0]),virB*calc_info_.nrio,0.0,&(C_p_AS[0][0]),
    virB*calc_info_.nrio);

  C_DGEMM('N','N',occA,virB*calc_info_.nrio,occB,1.0,&(Y_AB[0][0]),occB,
    &(B_p_BS[0][0]),virB*calc_info_.nrio,0.0,&(D_p_AS[0][0]),
    virB*calc_info_.nrio);

  free_block(B_p_BS);

  double **B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  double **C_p_BR = block_matrix(occB*virA,calc_info_.nrio);
  double **D_p_BR = block_matrix(occB*virA,calc_info_.nrio);

  C_DGEMM('T','N',occB,virA*calc_info_.nrio,occA,1.0,&(X_AB[0][0]),occB,
    &(B_p_AR[0][0]),virA*calc_info_.nrio,0.0,&(C_p_BR[0][0]),
    virA*calc_info_.nrio);

  C_DGEMM('T','N',occB,virA*calc_info_.nrio,occA,1.0,&(Y_AB[0][0]),occB,
    &(B_p_AR[0][0]),virA*calc_info_.nrio,0.0,&(D_p_BR[0][0]),
    virA*calc_info_.nrio);

  free_block(B_p_AR);
  free_block(X_AB);
  free_block(Y_AB);

  if (trans == 'n') {
    double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virA*occA,occB*virB);
    double **Y_RS = block_matrix(virA,virB);
    for (int a=0; a<occA; a++) {
      for (int b=0; b<occB; b++) {
        C_DGEMM('N','T',virA,virB,calc_info_.nrio,1.0,&(C_p_BR[b*virA][0]),
          calc_info_.nrio,&(D_p_AS[a*virB][0]),calc_info_.nrio,0.0,
          &(Y_RS[0][0]),virB);
        C_DGEMM('N','T',virA,virB,calc_info_.nrio,1.0,&(D_p_BR[b*virA][0]),
          calc_info_.nrio,&(C_p_AS[a*virB][0]),calc_info_.nrio,1.0,
          &(Y_RS[0][0]),virB);
        for (int r=0; r<virA; r++) {
          int ar = a*virA+r;
          energy += 2.0*C_DDOT(virB,&(tARBS[ar][b*virB]),1,&(Y_RS[r][0]),1);
        }
    }}
    free_block(tARBS);
    free_block(Y_RS);
  }
  else {
    double **tBSAR = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occA*virA);
    double **Y_SR = block_matrix(virB,virA);
    for (int b=0; b<occB; b++) {
      for (int a=0; a<occA; a++) {
        C_DGEMM('N','T',virB,virA,calc_info_.nrio,1.0,&(C_p_AS[a*virB][0]),
          calc_info_.nrio,&(D_p_BR[b*virA][0]),calc_info_.nrio,0.0,
          &(Y_SR[0][0]),virA);
        C_DGEMM('N','T',virB,virA,calc_info_.nrio,1.0,&(D_p_AS[a*virB][0]),
          calc_info_.nrio,&(C_p_BR[b*virA][0]),calc_info_.nrio,1.0,
          &(Y_SR[0][0]),virA);
        for (int s=0; s<virB; s++) {
          int bs = b*virB+s;
          energy += 2.0*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,&(Y_SR[s][0]),1);
        }
    }}
    free_block(tBSAR);
    free_block(Y_SR);
  }

  free_block(C_p_BR);
  free_block(C_p_AS);
  free_block(D_p_BR);
  free_block(D_p_AS);

  double **X_RS = block_matrix(virA,virB);
  double **Y_RS = block_matrix(virA,virB);

  for (int r=0; r<virA; r++) {
    C_DCOPY(virB,&(SAB[r+occA][occB]),1,&(X_RS[r][0]),1);
  }

  C_DGEMM('N','N',virA,virB,occC,1.0,&(SAC[occA][0]),calc_info_.nmo,
    &(SCB[0][occB]),calc_info_.nmo,0.0,&(Y_RS[0][0]),virB);

  B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
  C_p_BR = block_matrix(occB*virA,calc_info_.nrio);
  D_p_BR = block_matrix(occB*virA,calc_info_.nrio);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','N',virA,calc_info_.nrio,virB,1.0,&(X_RS[0][0]),virB,
      &(B_p_BS[b*virB][0]),calc_info_.nrio,0.0,&(C_p_BR[b*virA][0]),
      calc_info_.nrio);
    C_DGEMM('N','N',virA,calc_info_.nrio,virB,1.0,&(Y_RS[0][0]),virB,
      &(B_p_BS[b*virB][0]),calc_info_.nrio,0.0,&(D_p_BR[b*virA][0]),
      calc_info_.nrio);
  }

  free_block(B_p_BS);

  B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  C_p_AS = block_matrix(occA*virB,calc_info_.nrio);
  D_p_AS = block_matrix(occA*virB,calc_info_.nrio);

  for (int a=0; a<occA; a++) {
    C_DGEMM('T','N',virB,calc_info_.nrio,virA,1.0,&(X_RS[0][0]),virB,
      &(B_p_AR[a*virA][0]),calc_info_.nrio,0.0,&(C_p_AS[a*virB][0]),
      calc_info_.nrio);
    C_DGEMM('T','N',virB,calc_info_.nrio,virA,1.0,&(Y_RS[0][0]),virB,
      &(B_p_AR[a*virA][0]),calc_info_.nrio,0.0,&(D_p_AS[a*virB][0]),
      calc_info_.nrio);
  }

  free_block(B_p_AR);
  free_block(X_RS);
  free_block(Y_RS);

  if (trans == 'n') {
    double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virA*occA,occB*virB);
    double **Y_RS = block_matrix(virA,virB);
    for (int a=0; a<occA; a++) {
      for (int b=0; b<occB; b++) {
        C_DGEMM('N','T',virA,virB,calc_info_.nrio,1.0,&(C_p_BR[b*virA][0]),
          calc_info_.nrio,&(D_p_AS[a*virB][0]),calc_info_.nrio,0.0,
          &(Y_RS[0][0]),virB);
        C_DGEMM('N','T',virA,virB,calc_info_.nrio,1.0,&(D_p_BR[b*virA][0]),
          calc_info_.nrio,&(C_p_AS[a*virB][0]),calc_info_.nrio,1.0,
          &(Y_RS[0][0]),virB);
        for (int r=0; r<virA; r++) {
          int ar = a*virA+r;
          energy += 2.0*C_DDOT(virB,&(tARBS[ar][b*virB]),1,&(Y_RS[r][0]),1);
        }
    }}
    free_block(tARBS);
    free_block(Y_RS);
  }
  else {
    double **tBSAR = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occA*virA);
    double **Y_SR = block_matrix(virB,virA);
    for (int b=0; b<occB; b++) {
      for (int a=0; a<occA; a++) {
        C_DGEMM('N','T',virB,virA,calc_info_.nrio,1.0,&(C_p_AS[a*virB][0]),
          calc_info_.nrio,&(D_p_BR[b*virA][0]),calc_info_.nrio,0.0,
          &(Y_SR[0][0]),virA);
        C_DGEMM('N','T',virB,virA,calc_info_.nrio,1.0,&(D_p_AS[a*virB][0]),
          calc_info_.nrio,&(C_p_BR[b*virA][0]),calc_info_.nrio,1.0,
          &(Y_SR[0][0]),virA);
        for (int s=0; s<virB; s++) {
          int bs = b*virB+s;
          energy += 2.0*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,&(Y_SR[s][0]),1);
        }
    }}
    free_block(tBSAR);
    free_block(Y_SR);
  }

  free_block(C_p_BR);
  free_block(C_p_AS);
  free_block(D_p_BR);
  free_block(D_p_AS);

  return(energy);
}

double SAPT3BN5::exch_disp200_s3_3(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, int AAfile, char *AR_ints,
  int BBfile, char *BS_ints, char *T2label, char trans, int occA, int virA,
  int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **X_AS = block_matrix(occA,virB);
  double **Y_AS = block_matrix(occA,virB);

  for (int a=0; a<occA; a++) {
    C_DCOPY(virB,&(SAB[a][occB]),1,&(X_AS[a][0]),1);
  }

  C_DGEMM('N','N',occA,virB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SCB[0][occB]),calc_info_.nmo,0.0,&(Y_AS[0][0]),virB);

    double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
    double **C_p_AB = block_matrix(occA*occB,calc_info_.nrio);
    double **D_p_AB = block_matrix(occA*occB,calc_info_.nrio);

    for (int b=0; b<occB; b++) {
      C_DGEMM('N','N',occA,calc_info_.nrio,virB,1.0,&(X_AS[0][0]),virB,
        &(B_p_BS[b*virB][0]),calc_info_.nrio,0.0,&(C_p_AB[b][0]),
        occB*calc_info_.nrio);
      C_DGEMM('N','N',occA,calc_info_.nrio,virB,1.0,&(Y_AS[0][0]),virB,
        &(B_p_BS[b*virB][0]),calc_info_.nrio,0.0,&(D_p_AB[b][0]),
        occB*calc_info_.nrio);
    }

    free_block(B_p_BS);

    double **B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
    double **X_RB = block_matrix(virA,occB);
    double **Y_RB = block_matrix(virA,occB);

    for (int a=0; a<occA; a++) {
      C_DGEMM('N','T',virA,occB,calc_info_.nrio,1.0,&(B_p_AR[a*virA][0]),
        calc_info_.nrio,&(C_p_AB[a*occB][0]),calc_info_.nrio,1.0,
        &(X_RB[0][0]),occB);
      C_DGEMM('N','T',virA,occB,calc_info_.nrio,1.0,&(B_p_AR[a*virA][0]),
        calc_info_.nrio,&(D_p_AB[a*occB][0]),calc_info_.nrio,1.0,
        &(Y_RB[0][0]),occB);
    }

    free_block(B_p_AR);
    free_block(C_p_AB);
    free_block(D_p_AB);

  if (trans == 'n') {
    double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virA*occA,occB*virB);

    for (int a=0,ar=0; a<occA; a++) {
      for (int r=0; r<virA; r++,ar++) {
        for (int b=0; b<occB; b++) {
        energy -= 2.0*X_RB[r][b]*C_DDOT(virB,&(tARBS[ar][b*virB]),1,
                                        &(Y_AS[a][0]),1);
        energy -= 2.0*Y_RB[r][b]*C_DDOT(virB,&(tARBS[ar][b*virB]),1,
                                        &(X_AS[a][0]),1);
    }}}

    free_block(tARBS);
  }
  else {
    double **tBSAR = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occA*virA);

    for (int b=0,bs=0; b<occB; b++) {
      for (int s=0; s<virB; s++,bs++) {
        for (int a=0; a<occA; a++) {
        energy -= 2.0*X_AS[a][s]*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,
                                        &(Y_RB[0][b]),occB);
        energy -= 2.0*Y_AS[a][s]*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,
                                        &(X_RB[0][b]),occB);
    }}}

    free_block(tBSAR);
  }

  free_block(X_RB);
  free_block(Y_RB);
  free_block(X_AS);
  free_block(Y_AS);

  return(energy);
}

double SAPT3BN5::exch_disp200_s3_4(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, double **WABS, char *T2label,
  char trans, int occA, int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **X_AR = block_matrix(occA,virA);

  if (trans == 'n') {
    double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virA*occA,occB*virB);
    C_DGEMV('n',occA*virA,occB*virB,1.0,&(tARBS[0][0]),occB*virB,
      &(WABS[0][0]),1,0.0,&(X_AR[0][0]),1);
    free_block(tARBS);
  }
  else {
    double **tBSAR = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occA*virA);
    C_DGEMV('t',occB*virB,occA*virA,1.0,&(tBSAR[0][0]),occA*virA,
      &(WABS[0][0]),1,0.0,&(X_AR[0][0]),1);
    free_block(tBSAR);
  }

  double **Y_AB = block_matrix(occA,occB);
  double **Y_AC = block_matrix(occA,occC);
  double **Y_AR = block_matrix(occA,virA);

  C_DGEMM('N','N',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SCB[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  C_DGEMM('N','N',occA,occC,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  C_DGEMM('N','N',occA,virA,occB,1.0,&(Y_AB[0][0]),occB,&(SBA[0][occA]),
    calc_info_.nmo,0.0,&(Y_AR[0][0]),virA);

  C_DGEMM('N','N',occA,virA,occC,1.0,&(Y_AC[0][0]),occC,&(SCA[0][occA]),
    calc_info_.nmo,1.0,&(Y_AR[0][0]),virA);

  energy += 4.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(Y_AR[0][0]),1);

  free_block(X_AR);
  free_block(Y_AR);
  free_block(Y_AB);
  free_block(Y_AC);

  double **A_RB = block_matrix(virA,occB);
  double **A_AS = block_matrix(occA,virB);

  C_DGEMM('N','T',virA,occB,virB,1.0,&(SAB[occA][occB]),calc_info_.nmo,
    &(WABS[0][0]),virB,0.0,&(A_RB[0][0]),occB);

  C_DGEMM('N','N',occA,virB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SCB[0][occB]),calc_info_.nmo,0.0,&(A_AS[0][0]),virB);

  double **B_RB = block_matrix(virA,occB);
  double **B_AS = block_matrix(occA,virB);
  double **B_RS = block_matrix(virA,virB);

  C_DGEMM('N','N',virA,virB,occC,1.0,&(SAC[occA][0]),calc_info_.nmo,
    &(SCB[0][occB]),calc_info_.nmo,0.0,&(B_RS[0][0]),virB);

  C_DGEMM('N','T',virA,occB,virB,1.0,&(B_RS[0][0]),virB,&(WABS[0][0]),virB,
    0.0,&(B_RB[0][0]),occB);

  for (int a=0; a<occA; a++) {
    C_DCOPY(virB,&(SAB[a][occB]),1,&(B_AS[a][0]),1);
  }

  free_block(B_RS);

  double **C_RB = block_matrix(virA,occB);
  double **C_AS = block_matrix(occA,virB);

  C_DGEMM('N','N',virA,occB,occC,1.0,&(SAC[occA][0]),calc_info_.nmo,
    &(SCB[0][0]),calc_info_.nmo,0.0,&(C_RB[0][0]),occB);

  C_DGEMM('N','N',occA,virB,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(WABS[0][0]),
    virB,0.0,&(C_AS[0][0]),virB);

  double **D_RB = block_matrix(virA,occB);
  double **D_AS = block_matrix(occA,virB);
  double **D_AB = block_matrix(occA,occB);

  C_DGEMM('N','N',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SCB[0][0]),
    calc_info_.nmo,0.0,&(D_AB[0][0]),occB);

  C_DGEMM('N','N',occA,virB,occB,1.0,&(D_AB[0][0]),occB,&(WABS[0][0]),virB,
    0.0,&(D_AS[0][0]),virB);

  for (int r=0; r<virA; r++) {
    C_DCOPY(occB,&(SAB[r+occA][0]),1,&(D_RB[r][0]),1);
  }

  free_block(D_AB);

  if (trans == 'n') {
    double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virA*occA,occB*virB);

    for (int a=0,ar=0; a<occA; a++) {
      for (int r=0; r<virA; r++,ar++) {
        for (int b=0; b<occB; b++) {
          energy += 2.0*A_RB[r][b]*C_DDOT(virB,&(tARBS[ar][b*virB]),1,
                                          &(A_AS[a][0]),1);
          energy += 2.0*B_RB[r][b]*C_DDOT(virB,&(tARBS[ar][b*virB]),1,
                                          &(B_AS[a][0]),1);
          energy -= 2.0*C_RB[r][b]*C_DDOT(virB,&(tARBS[ar][b*virB]),1,
                                          &(C_AS[a][0]),1);
          energy -= 2.0*D_RB[r][b]*C_DDOT(virB,&(tARBS[ar][b*virB]),1,
                                          &(D_AS[a][0]),1);
    }}}

    free_block(tARBS);
  }
  else {
    double **tBSAR = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occA*virA);

    for (int b=0,bs=0; b<occB; b++) {
      for (int s=0; s<virB; s++,bs++) {
        for (int a=0; a<occA; a++) {
          energy += 2.0*A_AS[a][s]*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,
                                          &(A_RB[0][b]),occB);
          energy += 2.0*B_AS[a][s]*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,
                                          &(B_RB[0][b]),occB);
          energy -= 2.0*C_AS[a][s]*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,
                                          &(C_RB[0][b]),occB);
          energy -= 2.0*D_AS[a][s]*C_DDOT(virA,&(tBSAR[bs][a*virA]),1,
                                          &(D_RB[0][b]),occB);
    }}}

    free_block(tBSAR);
  }

  free_block(A_RB);
  free_block(B_RB);
  free_block(C_RB);
  free_block(D_RB);
  free_block(A_AS);
  free_block(B_AS);
  free_block(C_AS);
  free_block(D_AS);

  return(energy);
}

double SAPT3BN5::exch_disp110_s3_0(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, double **WABS, double **WCBS,
  double **WBAR, double **WBCT, char *TBS_CT, char *TBS_AR,
  int AAfile, char *AR_ints, int BBfile, char *BS_ints, int CCfile,
  char *CT_ints, char *tBSCT, char trans1, char *tARBS, char trans2, int occA,
  int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  energy += exch_disp110_s3_1(SAB,SBC,SCA,SAC,SBA,SCB,TBS_CT,AAfile,AR_ints,
    occA,virA,occB,virB,occC,virC);
  energy += exch_disp110_s3_1(SCB,SBA,SAC,SCA,SBC,SAB,TBS_AR,CCfile,CT_ints,
    occC,virC,occB,virB,occA,virA);
  energy += exch_disp110_s3_2(SAB,SBC,SCA,SAC,SBA,SCB,AAfile,AR_ints,BBfile,
    BS_ints,tBSCT,trans1,occA,virA,occB,virB,occC,virC);
  energy += exch_disp110_s3_2(SCB,SBA,SAC,SCA,SBC,SAB,CCfile,CT_ints,BBfile,
    BS_ints,tARBS,trans2,occC,virC,occB,virB,occA,virA);
  energy += exch_disp110_s3_3(SAB,SBC,SCA,SAC,SBA,SCB,WABS,tBSCT,trans1,
    occA,virA,occB,virB,occC,virC);
  energy += exch_disp110_s3_3(SCB,SBA,SAC,SCA,SBC,SAB,WCBS,tARBS,trans2,
    occC,virC,occB,virB,occA,virA);
  energy += exch_disp110_s3_4(SAB,SBC,SCA,SAC,SBA,SCB,WBAR,tBSCT,trans1,
    occA,virA,occB,virB,occC,virC);
  energy += exch_disp110_s3_4(SCB,SBA,SAC,SCA,SBC,SAB,WBCT,tARBS,trans2,
    occC,virC,occB,virB,occA,virA);

  return(energy);
}

double SAPT3BN5::exch_disp110_s3_1(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, char *TBS_CT, int AAfile,
  char *AR_ints, int occA, int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **X_AC = block_matrix(occA,occC);
  double **X_RT = block_matrix(virA,virC);

  for (int r=0; r<virA; r++) {
    C_DCOPY(virC,&(SAC[r+occA][occC]),1,&(X_RT[r][0]),1);
  }

  C_DGEMM('N','N',occA,occC,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(X_AC[0][0]),occC);

  double **Y_AC = block_matrix(occA,occC);
  double **Y_RT = block_matrix(virA,virC);

  for (int a=0; a<occA; a++) {
    C_DCOPY(occC,&(SAC[a][0]),1,&(Y_AC[a][0]),1);
  }

  C_DGEMM('N','N',virA,virC,occB,1.0,&(SAB[occA][0]),calc_info_.nmo,
    &(SBC[0][occC]),calc_info_.nmo,0.0,&(Y_RT[0][0]),virC);

  double **T_CT = read_IJKL(PSIF_3B_SAPT_AMPS,TBS_CT,occC*virC,
    calc_info_.nrio);
  double **C_p_CR = block_matrix(occC*virA,calc_info_.nrio);
  double **C_p_AR = block_matrix(occA*virA,calc_info_.nrio);

  for (int c=0; c<occC; c++) {
    C_DGEMM('N','N',virA,calc_info_.nrio,virC,1.0,&(X_RT[0][0]),virC,
      &(T_CT[c*virC][0]),calc_info_.nrio,0.0,&(C_p_CR[c*virA][0]),
      calc_info_.nrio);
  }

  C_DGEMM('N','N',occA,virA*calc_info_.nrio,occC,1.0,&(X_AC[0][0]),occC,
    &(C_p_CR[0][0]),virA*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
    virA*calc_info_.nrio);

  for (int c=0; c<occC; c++) {
    C_DGEMM('N','N',virA,calc_info_.nrio,virC,1.0,&(Y_RT[0][0]),virC,
      &(T_CT[c*virC][0]),calc_info_.nrio,0.0,&(C_p_CR[c*virA][0]),
      calc_info_.nrio);
  }

  C_DGEMM('N','N',occA,virA*calc_info_.nrio,occC,1.0,&(Y_AC[0][0]),occC,
    &(C_p_CR[0][0]),virA*calc_info_.nrio,1.0,&(C_p_AR[0][0]),
    virA*calc_info_.nrio);

  free_block(X_AC);
  free_block(X_RT);
  free_block(Y_AC);
  free_block(Y_RT);
  free_block(T_CT);
  free_block(C_p_CR);

  double **B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);

  energy = 4.0*C_DDOT(occA*virA*calc_info_.nrio,&(B_p_AR[0][0]),1,
    &(C_p_AR[0][0]),1);

  free_block(B_p_AR);
  free_block(C_p_AR);

  return(energy);
}

double SAPT3BN5::exch_disp110_s3_2(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, int AAfile, char *AR_ints,
  int BBfile, char *BS_ints, char *T2label, char trans, int occA, int virA,
  int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
  double **C_p_AB = block_matrix(occA*occB,calc_info_.nrio);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','N',occA,calc_info_.nrio,virB,1.0,&(SAB[0][occB]),
      calc_info_.nmo,&(B_p_BS[b*virB][0]),calc_info_.nrio,0.0,&(C_p_AB[b][0]),
      occB*calc_info_.nrio);
  }

  free_block(B_p_BS);

  double **B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  double **C_p_AT = block_matrix(occA*virC,calc_info_.nrio);

  for (int a=0; a<occA; a++) {
    C_DGEMM('N','N',virC,calc_info_.nrio,virA,1.0,&(SCA[occC][occA]),
      calc_info_.nmo,&(B_p_AR[a*virA][0]),calc_info_.nrio,0.0,
      &(C_p_AT[a*virC][0]),calc_info_.nrio);
  }

  free_block(B_p_AR);

  double **X_BT = block_matrix(occB,virC);

  for (int a=0; a<occA; a++) {
    C_DGEMM('N','T',occB,virC,calc_info_.nrio,1.0,&(C_p_AB[a*occB][0]),
      calc_info_.nrio,&(C_p_AT[a*virC][0]),calc_info_.nrio,1.0,&(X_BT[0][0]),
      virC);
  }

  free_block(C_p_AB);
  free_block(C_p_AT);

  B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  double **C_p_CR = block_matrix(occC*virA,calc_info_.nrio);

  C_DGEMM('N','N',occC,virA*calc_info_.nrio,occA,1.0,&(SCA[0][0]),
    calc_info_.nmo,&(B_p_AR[0][0]),virA*calc_info_.nrio,0.0,&(C_p_CR[0][0]),
    virA*calc_info_.nrio);

  free_block(B_p_AR);

  double **C_p_BC = block_matrix(occB*occC,calc_info_.nrio);

  for (int c=0; c<occC; c++) {
    C_DGEMM('N','N',occB,calc_info_.nrio,virA,1.0,&(SBA[0][occA]),
      calc_info_.nmo,&(C_p_CR[c*virA][0]),calc_info_.nrio,0.0,&(C_p_BC[c][0]),
      occC*calc_info_.nrio);
  }

  free_block(C_p_CR);

  B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
  double **X_SC = block_matrix(virB,occC);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','T',virB,occC,calc_info_.nrio,1.0,&(B_p_BS[b*virB][0]),
      calc_info_.nrio,&(C_p_BC[b*occC][0]),calc_info_.nrio,1.0,&(X_SC[0][0]),
      occC);
  }

  free_block(B_p_BS);
  free_block(C_p_BC);

  if (trans == 'n') {
    double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occC*virC);

    for (int b=0,bs=0; b<occB; b++) {
      for (int s=0; s<virB; s++,bs++) {
        for (int c=0; c<occC; c++) {
        energy += 2.0*SBC[s+occB][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                        &(X_BT[b][0]),1);
        energy -= 2.0*X_SC[s][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                        &(SBC[b][occC]),1);
    }}}

    free_block(tBSCT);
  }
  else {
    double **tCTBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virC*occC,occB*virB);

    for (int c=0,ct=0; c<occC; c++) {
      for (int t=0; t<virC; t++,ct++) {
        for (int b=0; b<occB; b++) {
        energy += 2.0*X_BT[b][t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                        &(SCB[c][occB]),1);
        energy -= 2.0*SBC[b][occC+t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                        &(X_SC[0][c]),occC);
    }}}

    free_block(tCTBS);
  }

  free_block(X_BT);
  free_block(X_SC);

  B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  double **C_p_AS = block_matrix(occA*virB,calc_info_.nrio);

 for (int a=0; a<occA; a++) {
    C_DGEMM('N','N',virB,calc_info_.nrio,virA,1.0,&(SBA[occB][occA]),
      calc_info_.nmo,&(B_p_AR[a*virA][0]),calc_info_.nrio,0.0,
      &(C_p_AS[a*virB][0]),calc_info_.nrio);
  }

  free_block(B_p_AR);

  double **C_p_CS = block_matrix(occC*virB,calc_info_.nrio);

  C_DGEMM('N','N',occC,virB*calc_info_.nrio,occA,1.0,&(SCA[0][0]),
    calc_info_.nmo,&(C_p_AS[0][0]),virB*calc_info_.nrio,0.0,&(C_p_CS[0][0]),
    virB*calc_info_.nrio);

  free_block(C_p_AS);

  B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
  double **C_p_BT = block_matrix(occB*virC,calc_info_.nrio);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','N',virC,calc_info_.nrio,virB,1.0,&(SCB[occC][occB]),
      calc_info_.nmo,&(B_p_BS[b*virB][0]),calc_info_.nrio,0.0,
      &(C_p_BT[b*virC][0]),calc_info_.nrio);
  }

  free_block(B_p_BS);

  if (trans == 'n') {
    double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occC*virC);
    double **X_ST = block_matrix(virB,virC);

    for (int c=0; c<occC; c++) {
      for (int b=0,bs=0; b<occB; b++) {
        C_DGEMM('N','T',virB,virC,calc_info_.nrio,1.0,&(C_p_CS[c*virB][0]),
          calc_info_.nrio,&(C_p_BT[b*virC][0]),calc_info_.nrio,0.0,
          &(X_ST[0][0]),virC);
        for (int s=0; s<virB; s++,bs++) {
          energy += 2.0*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,&(X_ST[s][0]),1);
    }}}

    free_block(X_ST);
    free_block(tBSCT);
  }
  else {
    double **tCTBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virC*occC,occB*virB);
    double **X_TS = block_matrix(virC,virB);

    for (int b=0; b<occB; b++) {
      for (int c=0,ct=0; c<occC; c++) {
        C_DGEMM('N','T',virC,virB,calc_info_.nrio,1.0,&(C_p_BT[b*virC][0]),
          calc_info_.nrio,&(C_p_CS[c*virB][0]),calc_info_.nrio,0.0,
          &(X_TS[0][0]),virB);
        for (int t=0; t<virC; t++,ct++) {
          energy += 2.0*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,&(X_TS[t][0]),1);
    }}}

    free_block(X_TS);
    free_block(tCTBS);
  }

  free_block(C_p_BT);
  free_block(C_p_CS);

  B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);
  C_p_AT = block_matrix(occA*virC,calc_info_.nrio);

  for (int a=0; a<occA; a++) {
    C_DGEMM('N','N',virC,calc_info_.nrio,virA,1.0,&(SCA[occC][occA]),
      calc_info_.nmo,&(B_p_AR[a*virA][0]),calc_info_.nrio,0.0,
      &(C_p_AT[a*virC][0]),calc_info_.nrio);
  }

  free_block(B_p_AR);

  C_p_BT = block_matrix(occB*virC,calc_info_.nrio);

  C_DGEMM('N','N',occB,virC*calc_info_.nrio,occA,1.0,&(SBA[0][0]),
    calc_info_.nmo,&(C_p_AT[0][0]),virC*calc_info_.nrio,0.0,&(C_p_BT[0][0]),
    virC*calc_info_.nrio);

  free_block(C_p_AT);

  B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);
  C_p_CS = block_matrix(occC*virB,calc_info_.nrio);

  C_DGEMM('N','N',occC,virB*calc_info_.nrio,occB,1.0,&(SCB[0][0]),
    calc_info_.nmo,&(B_p_BS[0][0]),virB*calc_info_.nrio,0.0,&(C_p_CS[0][0]),
    virB*calc_info_.nrio);

  free_block(B_p_BS);

  if (trans == 'n') {
    double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occC*virC);
    double **X_ST = block_matrix(virB,virC);

    for (int c=0; c<occC; c++) {
      for (int b=0,bs=0; b<occB; b++) {
        C_DGEMM('N','T',virB,virC,calc_info_.nrio,1.0,&(C_p_CS[c*virB][0]),
          calc_info_.nrio,&(C_p_BT[b*virC][0]),calc_info_.nrio,0.0,
          &(X_ST[0][0]),virC);
        for (int s=0; s<virB; s++,bs++) {
          energy -= 2.0*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,&(X_ST[s][0]),1);
    }}}

    free_block(X_ST);
    free_block(tBSCT);
  }
  else {
    double **tCTBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virC*occC,occB*virB);
    double **X_TS = block_matrix(virC,virB);

    for (int b=0; b<occB; b++) {
      for (int c=0,ct=0; c<occC; c++) {
        C_DGEMM('N','T',virC,virB,calc_info_.nrio,1.0,&(C_p_BT[b*virC][0]),
          calc_info_.nrio,&(C_p_CS[c*virB][0]),calc_info_.nrio,0.0,
          &(X_TS[0][0]),virB);
        for (int t=0; t<virC; t++,ct++) {
          energy -= 2.0*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,&(X_TS[t][0]),1);
    }}}

    free_block(X_TS);
    free_block(tCTBS);
  }

  free_block(C_p_BT);
  free_block(C_p_CS);

  return(energy);
}

double SAPT3BN5::exch_disp110_s3_3(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, double **WABS, char *T2label,
  char trans, int occA, int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **X_CT = block_matrix(occC,virC);

  if (trans == 'n') {
    double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occC*virC);
    C_DGEMV('t',virB*occB,occC*virC,1.0,&(tBSCT[0][0]),occC*virC,
      &(WABS[0][0]),1,0.0,&(X_CT[0][0]),1);
    free_block(tBSCT);
  }
  else {
    double **tCTBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virC*occC,occB*virB);
    C_DGEMV('n',virC*occC,occB*virB,1.0,&(tCTBS[0][0]),occB*virB,
      &(WABS[0][0]),1,0.0,&(X_CT[0][0]),1);
    free_block(tCTBS);
  }

  double **Y_CA = block_matrix(occC,occA);
  double **Y_CB = block_matrix(occC,occB);
  double **Y_CT = block_matrix(occC,virC);

  C_DGEMM('N','N',occC,occA,occB,1.0,&(SCB[0][0]),calc_info_.nmo,&(SBA[0][0]),
    calc_info_.nmo,0.0,&(Y_CA[0][0]),occA);

  C_DGEMM('N','N',occC,occB,occA,1.0,&(SCA[0][0]),calc_info_.nmo,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_CB[0][0]),occB);

  C_DGEMM('N','N',occC,virC,occA,1.0,&(Y_CA[0][0]),occA,&(SAC[0][occC]),
    calc_info_.nmo,0.0,&(Y_CT[0][0]),virC);

  C_DGEMM('N','N',occC,virC,occB,1.0,&(Y_CB[0][0]),occB,&(SBC[0][occC]),
    calc_info_.nmo,1.0,&(Y_CT[0][0]),virC);

  energy += 4.0*C_DDOT(occC*virC,&(X_CT[0][0]),1,&(Y_CT[0][0]),1);

  free_block(X_CT);
  free_block(Y_CT);
  free_block(Y_CA);
  free_block(Y_CB);

  double **A_BA = block_matrix(occB,occA);
  double **A_BT = block_matrix(occB,virC);
  double **A_SC = block_matrix(virB,occC);

  C_DGEMM('N','N',occB,occA,virB,1.0,&(WABS[0][0]),virB,
    &(SBA[occB][0]),calc_info_.nmo,0.0,&(A_BA[0][0]),occA);

  C_DGEMM('N','N',occB,virC,occA,1.0,&(A_BA[0][0]),occA,&(SAC[0][occC]),
    calc_info_.nmo,0.0,&(A_BT[0][0]),virC);

  for (int s=0; s<virB; s++) {
    C_DCOPY(occC,&(SBC[s+occB][0]),1,&(A_SC[s][0]),1);
  }

  free_block(A_BA);

  double **B_BT = block_matrix(occB,virC);
  double **B_SC = block_matrix(virB,occC);

  C_DGEMM('N','N',occB,virC,virB,1.0,&(WABS[0][0]),virB,&(SBC[occB][occC]),
    calc_info_.nmo,0.0,&(B_BT[0][0]),virC);

  C_DGEMM('N','N',virB,occC,occA,1.0,&(SBA[occB][0]),calc_info_.nmo,
    &(SAC[0][0]),calc_info_.nmo,0.0,&(B_SC[0][0]),occC);

  double **C_BT = block_matrix(occB,virC);
  double **C_SC = block_matrix(virB,occC);

  C_DGEMM('N','N',occB,virC,occA,1.0,&(SBA[0][0]),calc_info_.nmo,
    &(SAC[0][occC]),calc_info_.nmo,0.0,&(C_BT[0][0]),virC);

  C_DGEMM('T','N',virB,occC,occB,1.0,&(WABS[0][0]),virB,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(C_SC[0][0]),occC);

  double **D_BC = block_matrix(occB,occC);
  double **D_BT = block_matrix(occB,virC);
  double **D_SC = block_matrix(virB,occC);

  C_DGEMM('N','N',occB,occC,occA,1.0,&(SBA[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(D_BC[0][0]),occC);

  C_DGEMM('T','N',virB,occC,occB,1.0,&(WABS[0][0]),virB,&(D_BC[0][0]),occC,
    0.0,&(D_SC[0][0]),occC);

  for (int b=0; b<occB; b++) {
    C_DCOPY(virC,&(SBC[b][occC]),1,&(D_BT[b][0]),1);
  }

  free_block(D_BC);

  if (trans == 'n') {
    double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occC*virC);

    for (int b=0,bs=0; b<occB; b++) {
      for (int s=0; s<virB; s++,bs++) {
        for (int c=0; c<occC; c++) {
          energy += 2.0*A_SC[s][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                          &(A_BT[b][0]),1);
          energy += 2.0*B_SC[s][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                          &(B_BT[b][0]),1);
          energy -= 2.0*C_SC[s][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                          &(C_BT[b][0]),1);
          energy -= 2.0*D_SC[s][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                          &(D_BT[b][0]),1);
    }}}

    free_block(tBSCT);
  }
  else {
    double **tCTBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virC*occC,occB*virB);

    for (int c=0,ct=0; c<occC; c++) {
      for (int t=0; t<virC; t++,ct++) {
        for (int b=0; b<occB; b++) {
          energy += 2.0*A_BT[b][t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                          &(A_SC[0][c]),occC);
          energy += 2.0*B_BT[b][t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                          &(B_SC[0][c]),occC);
          energy -= 2.0*C_BT[b][t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                          &(C_SC[0][c]),occC);
          energy -= 2.0*D_BT[b][t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                          &(D_SC[0][c]),occC);
    }}}

    free_block(tCTBS);
  }

  free_block(A_BT);
  free_block(B_BT);
  free_block(C_BT);
  free_block(D_BT);
  free_block(A_SC);
  free_block(B_SC);
  free_block(C_SC);
  free_block(D_SC);

  return(energy);
}

double SAPT3BN5::exch_disp110_s3_4(double **SAB, double **SBC, double **SCA,
  double **SAC, double **SBA, double **SCB, double **WBAR, char *T2label,
  char trans, int occA, int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **A_AT = block_matrix(occA,virC);
  double **A_BT = block_matrix(occB,virC);
  double **A_SC = block_matrix(virB,occC);

  C_DGEMM('N','N',occA,virC,virA,1.0,&(WBAR[0][0]),virA,
    &(SAC[occA][occC]),calc_info_.nmo,0.0,&(A_AT[0][0]),virC);

  C_DGEMM('N','N',occB,virC,occA,1.0,&(SBA[0][0]),calc_info_.nmo,&(A_AT[0][0]),
    virC,0.0,&(A_BT[0][0]),virC);

  for (int s=0; s<virB; s++) {
    C_DCOPY(occC,&(SBC[s+occB][0]),1,&(A_SC[s][0]),1);
  }

  free_block(A_AT);

  double **B_SA = block_matrix(virB,occA);
  double **B_BT = block_matrix(occB,virC);
  double **B_SC = block_matrix(virB,occC);

  C_DGEMM('N','T',virB,occA,virA,1.0,&(SBA[occB][occA]),calc_info_.nmo,
    &(WBAR[0][0]),virA,0.0,&(B_SA[0][0]),occA);

  C_DGEMM('N','N',virB,occC,occA,1.0,&(B_SA[0][0]),occA,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(B_SC[0][0]),occC);

  for (int b=0; b<occB; b++) {
    C_DCOPY(virC,&(SBC[b][occC]),1,&(B_BT[b][0]),1);
  }

  free_block(B_SA);

  if (trans == 'n') {
    double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occC*virC);

    for (int b=0,bs=0; b<occB; b++) {
      for (int s=0; s<virB; s++,bs++) {
        for (int c=0; c<occC; c++) {
          energy += 2.0*A_SC[s][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                          &(A_BT[b][0]),1);
          energy += 2.0*B_SC[s][c]*C_DDOT(virC,&(tBSCT[bs][c*virC]),1,
                                          &(B_BT[b][0]),1);
    }}}

    free_block(tBSCT);
  }
  else {
    double **tCTBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virC*occC,occB*virB);

    for (int c=0,ct=0; c<occC; c++) {
      for (int t=0; t<virC; t++,ct++) {
        for (int b=0; b<occB; b++) {
          energy += 2.0*A_BT[b][t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                          &(A_SC[0][c]),occC);
          energy += 2.0*B_BT[b][t]*C_DDOT(virB,&(tCTBS[ct][b*virB]),1,
                                          &(B_SC[0][c]),occC);
    }}}

    free_block(tCTBS);
  }

  free_block(A_BT);
  free_block(B_BT);
  free_block(A_SC);
  free_block(B_SC);

  return(energy);
}

}}
