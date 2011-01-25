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
#include "sapt3bn5.h"

namespace psi { namespace sapt {

void SAPT3BN5::exch_ind200_s2()
{
  results_.exch_ind200_s2 = exch_ind200_s2_0(calc_info_.S_AC,calc_info_.S_BC,
    calc_info_.CHFB_A,calc_info_.CHFA_B,calc_info_.WBAR,calc_info_.WABS,
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC);
  fprintf(outfile,"exch_ind_200,r S^2  = %18.12lf  H\n",
    results_.exch_ind200_s2);
  fflush(outfile);
  results_.exch_ind020_s2 = exch_ind200_s2_0(calc_info_.S_CB,calc_info_.S_AB,
    calc_info_.CHFA_C,calc_info_.CHFC_A,calc_info_.WACT,calc_info_.WCAR,
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB);
  fprintf(outfile,"exch_ind_020,r S^2  = %18.12lf  H\n",
    results_.exch_ind020_s2);
  fflush(outfile);
  results_.exch_ind002_s2 = exch_ind200_s2_0(calc_info_.S_BA,calc_info_.S_CA,
    calc_info_.CHFC_B,calc_info_.CHFB_C,calc_info_.WCBS,calc_info_.WBCT,
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC,calc_info_.noccA);
  fprintf(outfile,"exch_ind_002,r S^2  = %18.12lf  H\n\n",
    results_.exch_ind002_s2);
  fflush(outfile);
}

void SAPT3BN5::exch_ind110_s2()
{
  results_.exch_ind110_s2 = exch_ind110_s2_0(calc_info_.S_CA,calc_info_.S_AC,
    calc_info_.S_CB,calc_info_.S_BC,calc_info_.S_BA,calc_info_.S_AB,
    calc_info_.CHFC_A,calc_info_.CHFA_C,calc_info_.CHFA_B,calc_info_.CHFB_A,
    calc_info_.WCAR,calc_info_.WACT,calc_info_.WABS,calc_info_.WBAR,
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"exch_ind_110,r S^2  = %18.12lf  H\n",
    results_.exch_ind110_s2);
  fflush(outfile);
  results_.exch_ind101_s2 = exch_ind110_s2_0(calc_info_.S_AB,calc_info_.S_BA,
    calc_info_.S_AC,calc_info_.S_CA,calc_info_.S_CB,calc_info_.S_BC,
    calc_info_.CHFA_B,calc_info_.CHFB_A,calc_info_.CHFB_C,calc_info_.CHFC_B,
    calc_info_.WABS,calc_info_.WBAR,calc_info_.WBCT,calc_info_.WCBS,
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"exch_ind_101,r S^2  = %18.12lf  H\n",
    results_.exch_ind101_s2);
  fflush(outfile);
  results_.exch_ind011_s2 = exch_ind110_s2_0(calc_info_.S_BC,calc_info_.S_CB,
    calc_info_.S_BA,calc_info_.S_AB,calc_info_.S_AC,calc_info_.S_CA,
    calc_info_.CHFB_C,calc_info_.CHFC_B,calc_info_.CHFC_A,calc_info_.CHFA_C,
    calc_info_.WBCT,calc_info_.WCBS,calc_info_.WCAR,calc_info_.WACT,
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"exch_ind_011,r S^2  = %18.12lf  H\n\n",
    results_.exch_ind011_s2);
  fflush(outfile);
}

void SAPT3BN5::exch_ind200_s3()
{
  results_.exch_ind200_s3 = exch_ind200_s3_0(calc_info_.S_AB,calc_info_.S_AC,
    calc_info_.S_BA,calc_info_.S_BC,calc_info_.CHFA_B,calc_info_.CHFB_A,
    calc_info_.WABS,calc_info_.WBAR,PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC);
  fprintf(outfile,"exch_ind_200,r S^3  = %18.12lf  H\n",
    results_.exch_ind200_s3);
  fflush(outfile);
  results_.exch_ind020_s3 = exch_ind200_s3_0(calc_info_.S_CA,calc_info_.S_CB,
    calc_info_.S_AC,calc_info_.S_AB,calc_info_.CHFC_A,calc_info_.CHFA_C,
    calc_info_.WCAR,calc_info_.WACT,PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.noccC,
    calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB);
  fprintf(outfile,"exch_ind_020,r S^3  = %18.12lf  H\n",
    results_.exch_ind020_s3);
  fflush(outfile);
  results_.exch_ind002_s3 = exch_ind200_s3_0(calc_info_.S_BC,calc_info_.S_BA,
    calc_info_.S_CB,calc_info_.S_CA,calc_info_.CHFB_C,calc_info_.CHFC_B,
    calc_info_.WBCT,calc_info_.WCBS,PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA);
  fprintf(outfile,"exch_ind_002,r S^3  = %18.12lf  H\n\n",
    results_.exch_ind002_s3);
  fflush(outfile);
}

void SAPT3BN5::exch_ind110_s3()
{
  results_.exch_ind110_s3 = exch_ind110_s3_0(calc_info_.S_CA,calc_info_.S_AC,
    calc_info_.S_CB,calc_info_.S_BC,calc_info_.S_BA,calc_info_.S_AB,
    calc_info_.CHFC_A,calc_info_.CHFA_C,calc_info_.CHFA_B,calc_info_.CHFB_A,
    calc_info_.WCAR,calc_info_.WACT,calc_info_.WABS,calc_info_.WBAR,
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"exch_ind_110,r S^3  = %18.12lf  H\n",
    results_.exch_ind110_s3);
  fflush(outfile);
  results_.exch_ind101_s3 = exch_ind110_s3_0(calc_info_.S_AB,calc_info_.S_BA,
    calc_info_.S_AC,calc_info_.S_CA,calc_info_.S_CB,calc_info_.S_BC,
    calc_info_.CHFA_B,calc_info_.CHFB_A,calc_info_.CHFB_C,calc_info_.CHFC_B,
    calc_info_.WABS,calc_info_.WBAR,calc_info_.WBCT,calc_info_.WCBS,
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"exch_ind_101,r S^3  = %18.12lf  H\n",
    results_.exch_ind101_s3);
  fflush(outfile);
  results_.exch_ind011_s3 = exch_ind110_s3_0(calc_info_.S_BC,calc_info_.S_CB,
    calc_info_.S_BA,calc_info_.S_AB,calc_info_.S_AC,calc_info_.S_CA,
    calc_info_.CHFB_C,calc_info_.CHFC_B,calc_info_.CHFC_A,calc_info_.CHFA_C,
    calc_info_.WBCT,calc_info_.WCBS,calc_info_.WCAR,calc_info_.WACT,
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"exch_ind_011,r S^3  = %18.12lf  H\n\n",
    results_.exch_ind011_s3);
  fflush(outfile);
}

double SAPT3BN5::exch_ind200_s2_0(double **SAC, double **SBC, double **tA,
  double **tB, double **WB, double **WA, int AAfile, char *AA_ints, int BBfile,
  char *BB_ints, int occA, int virA, int occB, int virB, int occC)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);
  double **Y_BB = block_matrix(occB,occB);

  C_DGEMM('N','T',occB,occB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  C_DGEMM('N','N',occB,virB,occB,1.0,&(Y_BB[0][0]),occB,&(tA[0][0]),
    virB,0.0,&(X_BS[0][0]),virB);

  free_block(Y_BB);

  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,virB,1.0,&(tA[0][0]),virB,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  C_DGEMM('N','T',occB,virB,occC,-1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BC);

  double **X_AR = block_matrix(occA,virA);
  double **Y_AA = block_matrix(occA,occA);

  C_DGEMM('N','T',occA,occA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_AA[0][0]),occA);

  C_DGEMM('N','N',occA,virA,occA,1.0,&(Y_AA[0][0]),occA,&(tB[0][0]),
    virA,0.0,&(X_AR[0][0]),virA);

  free_block(Y_AA);

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,virA,1.0,&(tB[0][0]),virA,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  C_DGEMM('N','T',occA,virA,occC,-1.0,&(Y_AC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AC);

  energy += 2.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(WB[0][0]),1);
  energy += 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WA[0][0]),1);

  free_block(X_AR);
  free_block(X_BS);

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);
  X_AR = block_matrix(occA,virA);
  X_BS = block_matrix(occB,virB);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(X_AR[0][0]),virA);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      for (int b=0; b<occB; b++) {
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          double tval = X_AR[a][r]*tA[b][s];
          tval += X_BS[b][s]*tB[a][r];
          energy -= 4.0*ARBS[ar][bs]*tval;
  }}}}

  free_block(X_AR);
  free_block(X_BS);
  free_block(ARBS);

  return(energy);
}

double SAPT3BN5::exch_ind110_s2_0(double **SAB, double **SBA, double **SAC,
  double **SCA, double **SCB, double **SBC, double **CA_B, double **CB_A,
  double **CB_C, double **CC_B, double **WAB, double **WBA, double **WBC,
  double **WCB, int AAfile, char *AA_ints, int BBfile, char *BB_ints,
  int CCfile, char *CC_ints, int occA, int virA, int occB, int virB, int occC,
  int virC)
{
  double energy = 0.0;

  energy += exch_ind110_s2_1(SAB,SAC,CB_C,CC_B,WBA,occA,virA,occB,virB,occC,
    virC);
  energy += exch_ind110_s2_1(SCA,SCB,CA_B,CB_A,WBC,occC,virC,occA,virA,occB,
    virB);
  energy += exch_ind110_s2_2(SAB,SBC,CB_C,CC_B,WAB,occA,virA,occB,virB,occC,
    virC);
  energy += exch_ind110_s2_2(SCB,SBA,CB_A,CA_B,WCB,occC,virC,occB,virB,occA,
    virA);
  energy += exch_ind110_s2_3(SAB,SAC,CB_C,AAfile,AA_ints,BBfile,BB_ints,
    occA,virA,occB,virB,occC);
  energy += exch_ind110_s2_3(SCB,SCA,CB_A,CCfile,CC_ints,BBfile,BB_ints,
    occC,virC,occB,virB,occA);

  return(energy);
}

double SAPT3BN5::exch_ind110_s2_1(double **SAB, double **SAC, double **CB_C,
  double **CC_B, double **WBA, int occA, int virA, int occB, int virB,
  int occC, int virC)
{
  double energy = 0.0;
  double **X_AR = block_matrix(occA,virA);
  double **Y_AS = block_matrix(occA,virB);

  C_DGEMM('N','N',occA,virB,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(CB_C[0][0]),
    virB,0.0,&(Y_AS[0][0]),virB);

  C_DGEMM('N','T',occA,virA,virB,1.0,&(Y_AS[0][0]),virB,&(SAB[occA][occB]),
    calc_info_.nmo,0.0,&(X_AR[0][0]),virA);

  free_block(Y_AS);

  double **Y_AT = block_matrix(occA,virC);

  C_DGEMM('N','N',occA,virC,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(CC_B[0][0]),
    virC,0.0,&(Y_AT[0][0]),virC);

  C_DGEMM('N','T',occA,virA,virC,1.0,&(Y_AT[0][0]),virC,&(SAC[occA][occC]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AT);

  energy -= 2.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(WBA[0][0]),1);

  free_block(X_AR);

  return(energy);
}

double SAPT3BN5::exch_ind110_s2_2(double **SAB, double **SBC, double **CB_C,
  double **CC_B, double **WAB, int occA, int virA, int occB, int virB,
  int occC, int virC)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);
  double **Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,virB,1.0,&(CB_C[0][0]),virB,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  C_DGEMM('N','N',occB,virB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  free_block(Y_BA);

  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,virB,1.0,&(CB_C[0][0]),virB,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BC);

  double **Y_BB = block_matrix(occB,occB);

  C_DGEMM('T','N',occB,occB,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  C_DGEMM('N','T',occB,occB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,1.0,&(Y_BB[0][0]),occB);

  C_DGEMM('N','N',occB,virB,occB,-1.0,&(Y_BB[0][0]),occB,&(CB_C[0][0]),
    virB,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BB);

  double **Y_BT = block_matrix(occB,virC);

  C_DGEMM('N','N',occB,virC,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(CC_B[0][0]),
    virC,0.0,&(Y_BT[0][0]),virC);

  C_DGEMM('N','T',occB,virB,virC,1.0,&(Y_BT[0][0]),virC,&(SBC[occB][occC]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BT);

  energy -= 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WAB[0][0]),1);

  free_block(X_BS);

  return(energy);
}

double SAPT3BN5::exch_ind110_s2_3(double **SAB, double **SAC, double **CB_C,
  int AAfile, char *AA_ints, int BBfile, char *BB_ints, int occA, int virA,
  int occB, int virB, int occC)
{
  double energy = 0.0;

  double **A_AS = block_matrix(occA,virB);

  C_DGEMM('N','N',occA,virB,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(CB_C[0][0]),
    virB,0.0,&(A_AS[0][0]),virB);

  double **B_AR = block_matrix(occA,virA);

  C_DGEMM('N','T',occA,virA,occB,1.0,&(SAB[0][0]),calc_info_.nmo,
    &(SAB[occA][0]),calc_info_.nmo,0.0,&(B_AR[0][0]),virA);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,1.0,&(B_AR[0][0]),virA);

  double **C_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occB,virA,virB,1.0,&(CB_C[0][0]),virB,&(SAB[occA][occB]),
    calc_info_.nmo,0.0,&(C_BR[0][0]),virA);

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      for (int b=0; b<occB; b++) {
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          double tval = A_AS[a][s]*SAB[r+occA][b];
          tval -= 2.0*B_AR[a][r]*CB_C[b][s];
          tval -= C_BR[b][r]*SAB[a][s+occB];
          energy += 2.0*ARBS[ar][bs]*tval;
  }}}}

  free_block(A_AS);
  free_block(B_AR);
  free_block(C_BR);
  free_block(ARBS);

  return(energy);
}

double SAPT3BN5::exch_ind200_s3_0(double **SAB, double **SAC, double **SBA,
  double **SBC, double **CA_B, double **CB_A, double **WAB, double **WBA,
  int AAfile, char *AA_ints, int BBfile, char *BB_ints, int occA, int virA,
  int occB, int virB, int occC)
{
  double energy = 0.0;

  energy += exch_ind200_s3_1(SAB,SAC,SBC,CA_B,CB_A,WAB,occA,virA,occB,virB,
    occC);
  energy += exch_ind200_s3_1(SBA,SBC,SAC,CB_A,CA_B,WBA,occB,virB,occA,virA,
    occC);
  energy += exch_ind200_s3_2(SAB,SAC,SBC,CA_B,AAfile,AA_ints,BBfile,BB_ints,
    occA,virA,occB,virB,occC);
  energy += exch_ind200_s3_2(SBA,SBC,SAC,CB_A,BBfile,BB_ints,AAfile,AA_ints,
    occB,virB,occA,virA,occC);

  return(energy);
}

double SAPT3BN5::exch_ind200_s3_1(double **SAB, double **SAC, double **SBC,
  double **CA_B, double **CB_A, double **WAB, int occA, int virA, int occB,
  int virB, int occC)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);
  double **Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','T',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  double **Y_BR = block_matrix(occB,virA);

  C_DGEMM('T','N',occB,virA,occA,1.0,&(Y_AB[0][0]),occB,&(CA_B[0][0]),
    virA,0.0,&(Y_BR[0][0]),virA);

  free_block(Y_AB);

  C_DGEMM('N','N',occB,virB,virA,1.0,&(Y_BR[0][0]),virA,&(SAB[occA][occB]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  free_block(Y_BR);

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,virA,1.0,&(CA_B[0][0]),virA,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('T','N',occB,occC,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(Y_AC[0][0]),
    occC,0.0,&(Y_BC[0][0]),occC);

  free_block(Y_AC);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BC);

  Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,virB,1.0,&(CB_A[0][0]),virB,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  double **Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,occC,1.0,&(Y_BC[0][0]),occC,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  free_block(Y_BC);

  C_DGEMM('N','N',occB,virB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BA);

  Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,virB,1.0,&(CB_A[0][0]),virB,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,occA,1.0,&(Y_BA[0][0]),occA,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  free_block(Y_BA);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BC);

  Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  double **Y_BB = block_matrix(occB,occB);

  C_DGEMM('N','N',occB,occB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  free_block(Y_BA);

  C_DGEMM('N','N',occB,virB,occB,-1.0,&(Y_BB[0][0]),occB,&(CB_A[0][0]),
    virB,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BB);

  Y_BC = block_matrix(occB,occC);

  C_DGEMM('T','N',occB,occC,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  Y_BB = block_matrix(occB,occB);

  C_DGEMM('N','T',occB,occB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  free_block(Y_BC);

  C_DGEMM('N','N',occB,virB,occB,-1.0,&(Y_BB[0][0]),occB,&(CB_A[0][0]),
    virB,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BB);

  energy += 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WAB[0][0]),1);

  free_block(X_BS);

  return(energy);
}

double SAPT3BN5::exch_ind200_s3_2(double **SAB, double **SAC, double **SBC,
  double **CA_B, int AAfile, char *AA_ints, int BBfile, char *BB_ints,
  int occA, int virA, int occB, int virB, int occC)
{
  double energy = 0.0;

  double **Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  double **A_BS = block_matrix(occB,virB);

  C_DGEMM('N','N',occB,virB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(A_BS[0][0]),virB);

  free_block(Y_BA);

  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('T','N',occB,occC,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,1.0,&(A_BS[0][0]),virB);

  free_block(Y_BC);

  double **B_AS = block_matrix(occA,virB);
  double **B_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occB,virA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(B_BR[0][0]),virA);

  C_DGEMM('N','N',occA,virB,virA,1.0,&(CA_B[0][0]),virA,&(SAB[occA][occB]),
    calc_info_.nmo,0.0,&(B_AS[0][0]),virB);

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,virA,1.0,&(CA_B[0][0]),virA,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  double **C_AS = block_matrix(occA,virB);

  C_DGEMM('N','T',occA,virB,occC,1.0,&(Y_AC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(C_AS[0][0]),virB);

  free_block(Y_AC);

  Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  double **D_BR = block_matrix(occB,virA);

  C_DGEMM('N','N',occB,virA,occA,1.0,&(Y_BA[0][0]),occA,&(CA_B[0][0]),
    virA,0.0,&(D_BR[0][0]),virA);

  double **E_AS = block_matrix(occA,virB);
  double **E_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occA,virB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(E_AS[0][0]),virB);

  C_DGEMM('T','N',occB,virA,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(CA_B[0][0]),
    virA,0.0,&(E_BR[0][0]),virA);

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      for (int b=0; b<occB; b++) {
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          double tval = 2.0*A_BS[b][s]*CA_B[a][r];
          tval += B_AS[a][s]*B_BR[b][r];
          tval += C_AS[a][s]*SAB[r+occA][b];
          tval -= SAB[a][s+occB]*D_BR[b][r];
          tval -= E_AS[a][s]*E_BR[b][r];
          energy += 2.0*ARBS[ar][bs]*tval;
  }}}}

  free_block(A_BS);
  free_block(B_AS);
  free_block(B_BR);
  free_block(C_AS);
  free_block(D_BR);
  free_block(E_AS);
  free_block(E_BR);
  free_block(ARBS);

  return(energy);
}

double SAPT3BN5::exch_ind110_s3_0(double **SAB, double **SBA, double **SAC,
  double **SCA, double **SCB, double **SBC, double **CA_B, double **CB_A,
  double **CB_C, double **CC_B, double **WAB, double **WBA, double **WBC,
  double **WCB, int AAfile, char *AA_ints, int BBfile, char *BB_ints,
  int CCfile, char *CC_ints, int occA, int virA, int occB, int virB, int occC,
  int virC)
{
  double energy = 0.0;

  energy += exch_ind110_s3_1(SAB,SAC,SBC,CB_C,CC_B,WAB,occA,virA,occB,virB,
    occC,virC);
  energy += exch_ind110_s3_1(SCB,SCA,SBA,CB_A,CA_B,WCB,occC,virC,occB,virB,
    occA,virA);
  energy += exch_ind110_s3_2(SAB,SAC,SBC,CB_C,CC_B,WBA,occA,virA,occB,virB,
    occC,virC);
  energy += exch_ind110_s3_2(SCB,SCA,SBA,CB_A,CA_B,WBC,occC,virC,occB,virB,
    occA,virA);
  energy += exch_ind110_s3_3(SAB,SAC,SBC,CB_C,CC_B,AAfile,AA_ints,BBfile,
    BB_ints,occA,virA,occB,virB,occC,virC);
  energy += exch_ind110_s3_3(SCB,SCA,SBA,CB_A,CA_B,CCfile,CC_ints,BBfile,
    BB_ints,occC,virC,occB,virB,occA,virA);

  return(energy);
}

double SAPT3BN5::exch_ind110_s3_1(double **SAB, double **SAC, double **SBC,
  double **CB_C, double **CC_B, double **WAB, int occA, int virA, int occB,
  int virB, int occC, int virC)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);
  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,virB,1.0,&(CB_C[0][0]),virB,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  double **Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,occC,1.0,&(Y_BC[0][0]),occC,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  free_block(Y_BC);

  C_DGEMM('N','N',occB,virB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  free_block(Y_BA);

  Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,virB,1.0,&(CB_C[0][0]),virB,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,occA,1.0,&(Y_BA[0][0]),occA,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  free_block(Y_BA);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BC);

  Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  double **Y_BB = block_matrix(occB,occB);

  C_DGEMM('N','N',occB,occB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  free_block(Y_BA);

  C_DGEMM('N','N',occB,virB,occB,-1.0,&(Y_BB[0][0]),occB,&(CB_C[0][0]),
    virB,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BB);

  Y_BC = block_matrix(occB,occC);

  C_DGEMM('T','N',occB,occC,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  Y_BB = block_matrix(occB,occB);

  C_DGEMM('N','T',occB,occB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  free_block(Y_BC);

  C_DGEMM('N','N',occB,virB,occB,-1.0,&(Y_BB[0][0]),occB,&(CB_C[0][0]),
    virB,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BB);

  double **Y_CA = block_matrix(occC,occA);

  C_DGEMM('N','T',occC,occA,virC,1.0,&(CC_B[0][0]),virC,&(SAC[0][occC]),
    calc_info_.nmo,0.0,&(Y_CA[0][0]),occA);

  Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','N',occB,occA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(Y_CA[0][0]),
    occA,0.0,&(Y_BA[0][0]),occA);

  free_block(Y_CA);

  C_DGEMM('N','N',occB,virB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BA);

  double **Y_AT = block_matrix(occA,virC);

  C_DGEMM('N','N',occA,virC,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(CC_B[0][0]),
    virC,0.0,&(Y_AT[0][0]),virC);

  double **Y_BT = block_matrix(occB,virC);

  C_DGEMM('T','N',occB,virC,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(Y_AT[0][0]),
    virC,0.0,&(Y_BT[0][0]),virC);

  free_block(Y_AT);

  C_DGEMM('N','T',occB,virB,virC,1.0,&(Y_BT[0][0]),virC,&(SBC[occB][occC]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BT);

  energy += 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WAB[0][0]),1);

  free_block(X_BS);

  return(energy);
}

double SAPT3BN5::exch_ind110_s3_2(double **SAB, double **SAC, double **SBC,
  double **CB_C, double **CC_B, double **WBA, int occA, int virA, int occB,
  int virB, int occC, int virC)
{
  double energy = 0.0;
  double **X_AR = block_matrix(occA,virA);
  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,virB,1.0,&(CB_C[0][0]),virB,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(Y_BC[0][0]),
    occC,0.0,&(Y_AC[0][0]),occC);

  free_block(Y_BC);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(Y_AC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(X_AR[0][0]),virA);

  free_block(Y_AC);

  double **Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','T',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  double **Y_AS = block_matrix(occA,virB);

  C_DGEMM('N','N',occA,virB,occB,1.0,&(Y_AB[0][0]),occB,&(CB_C[0][0]),
    virB,0.0,&(Y_AS[0][0]),virB);

  free_block(Y_AB);

  C_DGEMM('N','T',occA,virA,virB,1.0,&(Y_AS[0][0]),virB,&(SAB[occA][occB]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AS);

  Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  double **Y_AT = block_matrix(occA,virC);

  C_DGEMM('N','N',occA,virC,occC,1.0,&(Y_AC[0][0]),occC,&(CC_B[0][0]),
    virC,0.0,&(Y_AT[0][0]),virC);

  free_block(Y_AC);

  C_DGEMM('N','T',occA,virA,virC,1.0,&(Y_AT[0][0]),virC,&(SAC[occA][occC]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AT);

  double **Y_CB = block_matrix(occC,occB);

  C_DGEMM('N','T',occC,occB,virC,1.0,&(CC_B[0][0]),virC,&(SBC[0][occC]),
    calc_info_.nmo,0.0,&(Y_CB[0][0]),occB);

  Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','N',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(Y_CB[0][0]),
    occB,0.0,&(Y_AB[0][0]),occB);

  free_block(Y_CB);

  C_DGEMM('N','T',occA,virA,occB,1.0,&(Y_AB[0][0]),occB,&(SAB[occA][0]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AB);

  energy += 2.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(WBA[0][0]),1);

  free_block(X_AR);

  return(energy);
}

double SAPT3BN5::exch_ind110_s3_3(double **SAB, double **SAC, double **SBC,
  double **CB_C, double **CC_B, int AAfile, char *AA_ints, int BBfile,
  char *BB_ints, int occA, int virA, int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  double **A_AR = block_matrix(occA,virA);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(Y_AC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(A_AR[0][0]),virA);

  free_block(Y_AC);

  double **Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','T',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  C_DGEMM('N','T',occA,virA,occB,1.0,&(Y_AB[0][0]),occB,&(SAB[occA][0]),
    calc_info_.nmo,1.0,&(A_AR[0][0]),virA);

  free_block(Y_AB);

  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,virB,1.0,&(CB_C[0][0]),virB,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  double **B_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occB,virA,occC,1.0,&(Y_BC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(B_BR[0][0]),virA);

  free_block(Y_BC);

  double **C_AS = block_matrix(occA,virB);
  double **C_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occA,virB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(C_AS[0][0]),virB);

  C_DGEMM('N','T',occB,virA,virB,1.0,&(CB_C[0][0]),virB,&(SAB[occA][occB]),
    calc_info_.nmo,0.0,&(C_BR[0][0]),virA);

  double **Y_BT = block_matrix(occB,virC);

  C_DGEMM('N','N',occB,virC,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(CC_B[0][0]),
    virC,0.0,&(Y_BT[0][0]),virC);

  double **D_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occB,virA,virC,1.0,&(Y_BT[0][0]),virC,&(SAC[occA][occC]),
    calc_info_.nmo,0.0,&(D_BR[0][0]),virA);

  free_block(Y_BT);

  double **Y_AT = block_matrix(occA,virC);

  C_DGEMM('N','N',occA,virC,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(CC_B[0][0]),
    virC,0.0,&(Y_AT[0][0]),virC);

  double **E_AS = block_matrix(occA,virB);

  C_DGEMM('N','T',occA,virB,virC,1.0,&(Y_AT[0][0]),virC,&(SBC[occB][occC]),
    calc_info_.nmo,0.0,&(E_AS[0][0]),virB);

  free_block(Y_AT);

  double **F_AS = block_matrix(occA,virB);
  double **F_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occB,virA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(F_BR[0][0]),virA);

  C_DGEMM('N','N',occA,virB,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(CB_C[0][0]),
    virB,0.0,&(F_AS[0][0]),virB);

  Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','T',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  double **G_AS = block_matrix(occA,virB);

  C_DGEMM('N','N',occA,virB,occB,1.0,&(Y_AB[0][0]),occB,&(CB_C[0][0]),
    virB,0.0,&(G_AS[0][0]),virB);

  free_block(Y_AB);

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      for (int b=0; b<occB; b++) {
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          double tval = 2.0*A_AR[a][r]*CB_C[b][s];
          tval += B_BR[b][r]*SAB[a][s+occB];
          tval += C_BR[b][r]*C_AS[a][s];
          tval += D_BR[b][r]*SAB[a][s+occB];
          tval += SAB[r+occA][b]*E_AS[a][s];
          tval -= F_BR[b][r]*F_AS[a][s];
          tval -= SAB[r+occA][b]*G_AS[a][s];
          energy += 2.0*ARBS[ar][bs]*tval;
  }}}}

  free_block(A_AR);
  free_block(B_BR);
  free_block(C_BR);
  free_block(C_AS);
  free_block(D_BR);
  free_block(E_AS);
  free_block(F_BR);
  free_block(F_AS);
  free_block(G_AS);
  free_block(ARBS);

  return(energy);
}

}}
