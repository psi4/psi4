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

void SAPT3BN5::exch100_s2()
{
  results_.exch100_s2 = exch_s2(calc_info_.S_AC,calc_info_.S_BC,
    calc_info_.WBAR,calc_info_.WABS,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC);
  fprintf(outfile,"exch_100 S^2        = %18.12lf  H\n",results_.exch100_s2);
  fflush(outfile);
  results_.exch010_s2 = exch_s2(calc_info_.S_CB,calc_info_.S_AB,
    calc_info_.WACT,calc_info_.WCAR,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB);
  fprintf(outfile,"exch_010 S^2        = %18.12lf  H\n",results_.exch010_s2);
  fflush(outfile);
  results_.exch001_s2 = exch_s2(calc_info_.S_BA,calc_info_.S_CA,
    calc_info_.WCBS,calc_info_.WBCT,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA);
  fprintf(outfile,"exch_001 S^2        = %18.12lf  H\n\n",results_.exch001_s2);
  fflush(outfile);
}

void SAPT3BN5::exch100_s3()
{
  results_.exch100_s3 = exch_s3(calc_info_.S_AB,calc_info_.S_AC,
    calc_info_.S_BC,calc_info_.WBAR,calc_info_.WABS,PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB,
    calc_info_.noccC);
  fprintf(outfile,"exch_100 S^3        = %18.12lf  H\n",results_.exch100_s3);
  fflush(outfile);
  results_.exch010_s3 = exch_s3(calc_info_.S_CA,calc_info_.S_CB,
    calc_info_.S_AB,calc_info_.WACT,calc_info_.WCAR,PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB);
  fprintf(outfile,"exch_010 S^3        = %18.12lf  H\n",results_.exch010_s3);
  fflush(outfile);
  results_.exch001_s3 = exch_s3(calc_info_.S_BC,calc_info_.S_BA,
    calc_info_.S_CA,calc_info_.WCBS,calc_info_.WBCT,PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC,
    calc_info_.noccA);
  fprintf(outfile,"exch_001 S^3        = %18.12lf  H\n\n",results_.exch001_s3);
  fflush(outfile);
}

void SAPT3BN5::exch100_s4()
{
  results_.exch100_s4 = exch_s4(calc_info_.S_AB,calc_info_.S_BA,
    calc_info_.S_AC,calc_info_.S_BC,calc_info_.WBAR,calc_info_.WABS,
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC);
  fprintf(outfile,"exch_100 S^4        = %18.12lf  H\n",results_.exch100_s4);
  fflush(outfile);
  results_.exch010_s4 = exch_s4(calc_info_.S_CA,calc_info_.S_AC,
    calc_info_.S_CB,calc_info_.S_AB,calc_info_.WACT,calc_info_.WCAR,
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB);
  fprintf(outfile,"exch_010 S^4        = %18.12lf  H\n",results_.exch010_s4);
  fflush(outfile);
  results_.exch001_s4 = exch_s4(calc_info_.S_BC,calc_info_.S_CB,
    calc_info_.S_BA,calc_info_.S_CA,calc_info_.WCBS,calc_info_.WBCT,
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC,calc_info_.noccA);
  fprintf(outfile,"exch_001 S^4        = %18.12lf  H\n\n",results_.exch001_s4);
  fflush(outfile);
}

double SAPT3BN5::exch_s2(double **SAC, double **SBC, double **WB, double **WA,
  int occA, int virA, int occB, int virB, int occC)
{
  double energy = 0.0;
  double **X_AR = block_matrix(occA,virA);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(X_AR[0][0]),virA);

  double **X_BS = block_matrix(occB,virB);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  energy -= 2.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(WB[0][0]),1);
  energy -= 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WA[0][0]),1);

  free_block(X_AR);
  free_block(X_BS);

  return(energy);
}

double SAPT3BN5::exch_s3(double **SAB, double **SAC, double **SBC, double **WB,
  double **WA, int AAfile, char *AA_ints, int BBfile, char *BB_ints, int occA,
  int virA, int occB, int virB, int occC)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);
  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('T','N',occB,occC,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  free_block(Y_BC);

  double **Y_BA = block_matrix(occB,occA);

  C_DGEMM('N','T',occB,occA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BA[0][0]),occA);

  C_DGEMM('N','N',occB,virB,occA,1.0,&(Y_BA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,1.0,&(X_BS[0][0]),virB);

  free_block(Y_BA);

  double **X_AR = block_matrix(occA,virA);
  double **Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','T',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  C_DGEMM('N','T',occA,virA,occB,1.0,&(Y_AB[0][0]),occB,&(SAB[occA][0]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AB);

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(Y_AC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AC);

  energy += 2.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(WB[0][0]),1);
  energy += 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WA[0][0]),1);

  free_block(X_AR);
  free_block(X_BS);

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);
  double **X_AS = block_matrix(occA,virB);
  double **X_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occA,virB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(X_AS[0][0]),virB);

  C_DGEMM('N','T',occB,virA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(X_BR[0][0]),virA);

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      for (int b=0; b<occB; b++) {
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          double tval = X_AS[a][s]*SAB[r+occA][b];
          tval += X_BR[b][r]*SAB[a][s+occB];
          energy += 2.0*ARBS[ar][bs]*tval;
  }}}}

  free_block(X_AS);
  free_block(X_BR);
  free_block(ARBS);

  return(energy);
}

double SAPT3BN5::exch_s4(double **SAB, double **SBA, double **SAC,
  double **SBC, double **WB, double **WA, int AAfile, char *AA_ints,
  int BBfile, char *BB_ints, int occA, int virA, int occB, int virB, int occC)
{
  double energy = 0.0;

  energy += exch_s4_1(SBC,WA,occB,virB,occC);
  energy += exch_s4_1(SAC,WB,occA,virA,occC);
  energy += exch_s4_2(SAB,SAC,WB,occA,virA,occB,occC);
  energy += exch_s4_2(SBA,SBC,WA,occB,virB,occA,occC);
  energy += exch_s4_3(SAB,SAC,WA,occA,occB,virB,occC);
  energy += exch_s4_3(SBA,SBC,WB,occB,occA,virA,occC);
  energy += exch_s4_4(SAC,SBC,WB,occA,virA,occB,occC);
  energy += exch_s4_4(SBC,SAC,WA,occB,virB,occA,occC);
  energy += exch_s4_5(SAB,SAC,AAfile,AA_ints,BBfile,BB_ints,occA,virA,occB,
                      virB,occC);
  energy += exch_s4_5(SBA,SBC,BBfile,BB_ints,AAfile,AA_ints,occB,virB,occA,
                      virA,occC);
  energy += exch_s4_6(SAC,SBC,AAfile,AA_ints,BBfile,BB_ints,occA,virA,occB,
                      virB,occC);

  return(energy);
}


double SAPT3BN5::exch_s4_1(double **SBC, double **WA, int occB, int virB,
  int occC)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);
  double **Y_BB = block_matrix(occB,occB);

  C_DGEMM('N','T',occB,occB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_BB[0][0]),occB);

  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('N','N',occB,occC,occB,1.0,&(Y_BB[0][0]),occB,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  free_block(Y_BB);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(Y_BC[0][0]),occC,&(SBC[occB][0]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  free_block(Y_BC);

  energy -= 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WA[0][0]),1);

  free_block(X_BS);

  return(energy);
}

double SAPT3BN5::exch_s4_2(double **SAB, double **SAC, double **WB, int occA,
  int virA, int occB, int occC)
{
  double energy = 0.0;
  double **X_AR = block_matrix(occA,virA);
  double **Y_AA = block_matrix(occA,occA);

  C_DGEMM('N','T',occA,occA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_AA[0][0]),occA);

  double **Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','N',occA,occB,occA,1.0,&(Y_AA[0][0]),occA,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  free_block(Y_AA);

  C_DGEMM('N','T',occA,virA,occB,1.0,&(Y_AB[0][0]),occB,&(SAB[occA][0]),
    calc_info_.nmo,0.0,&(X_AR[0][0]),virA);

  free_block(Y_AB);

  Y_AA = block_matrix(occA,occA);

  C_DGEMM('N','T',occA,occA,occB,1.0,&(SAB[0][0]),calc_info_.nmo,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_AA[0][0]),occA);

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,occA,1.0,&(Y_AA[0][0]),occA,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  free_block(Y_AA);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(Y_AC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,1.0,&(X_AR[0][0]),virA);

  free_block(Y_AC);

  energy -= 2.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(WB[0][0]),1);

  free_block(X_AR);

  return(energy);
}

double SAPT3BN5::exch_s4_3(double **SAB, double **SAC, double **WA, int occA,
  int occB, int virB, int occC)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);
  double **Y_AA = block_matrix(occA,occA);

  C_DGEMM('N','T',occA,occA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_AA[0][0]),occA);

  double **Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','N',occA,occB,occA,1.0,&(Y_AA[0][0]),occA,&(SAB[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  free_block(Y_AA);

  C_DGEMM('T','N',occB,virB,occA,1.0,&(Y_AB[0][0]),occB,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),virB);

  free_block(Y_AB);

  energy -= 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(WA[0][0]),1);

  free_block(X_BS);

  return(energy);
}

double SAPT3BN5::exch_s4_4(double **SAC, double **SBC, double **WB, int occA,
  int virA, int occB, int occC)
{
  double energy = 0.0;
  double **X_AR = block_matrix(occA,virA);
  double **Y_AB = block_matrix(occA,occB);

  C_DGEMM('N','T',occA,occB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AB[0][0]),occB);

  double **Y_AC = block_matrix(occA,occC);

  C_DGEMM('N','N',occA,occC,occB,1.0,&(Y_AB[0][0]),occB,&(SBC[0][0]),
    calc_info_.nmo,0.0,&(Y_AC[0][0]),occC);

  free_block(Y_AB);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(Y_AC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(X_AR[0][0]),virA);

  free_block(Y_AC);

  energy -= 2.0*C_DDOT(occA*virA,&(X_AR[0][0]),1,&(WB[0][0]),1);

  free_block(X_AR);

  return(energy);
}

double SAPT3BN5::exch_s4_5(double **SAB, double **SAC, int AAfile,
  char *AA_ints, int BBfile, char *BB_ints, int occA, int virA, int occB,
  int virB, int occC)
{
  double energy = 0.0;
  double **D_AR = block_matrix(occA,virA);
  double **D_BS = block_matrix(occB,virB);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(D_AR[0][0]),virA);

  C_DGEMM('T','N',occB,virB,occA,1.0,&(SAB[0][0]),calc_info_.nmo,
    &(SAB[0][occB]),calc_info_.nmo,0.0,&(D_BS[0][0]),virB);

  double **Y_AA = block_matrix(occA,occA);

  C_DGEMM('N','T',occA,occA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_AA[0][0]),occA);

  double **E_AS = block_matrix(occA,virB);

  C_DGEMM('N','N',occA,virB,occA,1.0,&(Y_AA[0][0]),occA,&(SAB[0][occB]),
    calc_info_.nmo,0.0,&(E_AS[0][0]),virB);

  free_block(Y_AA);

  double **Y_BC = block_matrix(occB,occC);

  C_DGEMM('T','N',occB,occC,occA,1.0,&(SAB[0][0]),calc_info_.nmo,&(SAC[0][0]),
    calc_info_.nmo,0.0,&(Y_BC[0][0]),occC);

  double **F_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occB,virA,occC,1.0,&(Y_BC[0][0]),occC,&(SAC[occA][0]),
    calc_info_.nmo,0.0,&(F_BR[0][0]),virA);

  free_block(Y_BC);

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      for (int b=0; b<occB; b++) {
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          double tval = 2.0*D_AR[a][r]*D_BS[b][s];
          tval -= E_AS[a][s]*SAB[r+occA][b];
          tval -= F_BR[b][r]*SAB[a][s+occB];
          energy += 2.0*ARBS[ar][bs]*tval;
  }}}}

  free_block(D_AR);
  free_block(D_BS);
  free_block(E_AS);
  free_block(F_BR);
  free_block(ARBS);

  return(energy);
}

double SAPT3BN5::exch_s4_6(double **SAC, double **SBC, int AAfile,
  char *AA_ints, int BBfile, char *BB_ints, int occA, int virA, int occB,
  int virB, int occC)
{
  double energy = 0.0;
  double **D_AR = block_matrix(occA,virA);
  double **D_BS = block_matrix(occB,virB);

  C_DGEMM('N','T',occA,virA,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(D_AR[0][0]),virA);

  C_DGEMM('N','T',occB,virB,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(D_BS[0][0]),virB);

  double **E_AS = block_matrix(occA,virB);
  double **E_BR = block_matrix(occB,virA);

  C_DGEMM('N','T',occA,virB,occC,1.0,&(SAC[0][0]),calc_info_.nmo,
    &(SBC[occB][0]),calc_info_.nmo,0.0,&(E_AS[0][0]),virB);

  C_DGEMM('N','T',occB,virA,occC,1.0,&(SBC[0][0]),calc_info_.nmo,
    &(SAC[occA][0]),calc_info_.nmo,0.0,&(E_BR[0][0]),virA);

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      for (int b=0; b<occB; b++) {
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          double tval = 2.0*D_AR[a][r]*D_BS[b][s];
          tval -= E_AS[a][s]*E_BR[b][r];
          energy += 2.0*ARBS[ar][bs]*tval;
  }}}}

  free_block(D_AR);
  free_block(D_BS);
  free_block(E_AS);
  free_block(E_BR);
  free_block(ARBS);

  return(energy);
}

}}
