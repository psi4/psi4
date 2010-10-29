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

void SAPT3BN5::ind110()
{
  results_.ind110 = ind110_1(calc_info_.CHFA_C,calc_info_.CHFA_B,
    calc_info_.WCAR,calc_info_.WBAR,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"ind_110,r           = %18.12lf  H\n",results_.ind110);
  fflush(outfile);
  results_.ind101 = ind110_1(calc_info_.CHFB_A,calc_info_.CHFB_C,
    calc_info_.WABS,calc_info_.WCBS,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"ind_101,r           = %18.12lf  H\n",results_.ind101);
  fflush(outfile);
  results_.ind011 = ind110_1(calc_info_.CHFC_B,calc_info_.CHFC_A,
    calc_info_.WBCT,calc_info_.WACT,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"ind_011,r           = %18.12lf  H\n\n",results_.ind011);
  fflush(outfile);
}

void SAPT3BN5::ind111()
{
  results_.ind111 = ind111_1(calc_info_.sA_C,calc_info_.sB_C,
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB);
  results_.ind111 += ind111_1(calc_info_.sB_A,calc_info_.sC_A,
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC);
  results_.ind111 += ind111_1(calc_info_.sC_B,calc_info_.sA_B,
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,
    calc_info_.nvirA);
  fprintf(outfile,"ind_111             = %18.12lf  H\n\n",results_.ind111);
  fflush(outfile);
}

void SAPT3BN5::ind210()
{
  results_.ind210 = ind210_0(calc_info_.sB_A,calc_info_.sA_C,calc_info_.sA_B,
    calc_info_.WBRR,calc_info_.WBAA,calc_info_.WCRR,calc_info_.WCAA,
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,
    calc_info_.nvirA);
  fprintf(outfile,"ind_210             = %18.12lf  H\n",results_.ind210);
  fflush(outfile);
  results_.ind201 = ind210_0(calc_info_.sA_B,calc_info_.sB_C,calc_info_.sB_A,
    calc_info_.WASS,calc_info_.WABB,calc_info_.WCSS,calc_info_.WCBB,
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB);
  fprintf(outfile,"ind_201             = %18.12lf  H\n",results_.ind201);
  fflush(outfile);
  results_.ind120 = ind210_0(calc_info_.sC_A,calc_info_.sA_B,calc_info_.sA_C,
    calc_info_.WCRR,calc_info_.WCAA,calc_info_.WBRR,calc_info_.WBAA,
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,
    calc_info_.nvirA);
  fprintf(outfile,"ind_120             = %18.12lf  H\n",results_.ind120);
  fflush(outfile);
  results_.ind021 = ind210_0(calc_info_.sA_C,calc_info_.sC_B,calc_info_.sC_A,
    calc_info_.WATT,calc_info_.WACC,calc_info_.WBTT,calc_info_.WBCC,
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccC,
    calc_info_.nvirC);
  fprintf(outfile,"ind_021             = %18.12lf  H\n",results_.ind021);
  fflush(outfile);
  results_.ind102 = ind210_0(calc_info_.sC_B,calc_info_.sB_A,calc_info_.sB_C,
    calc_info_.WCSS,calc_info_.WCBB,calc_info_.WASS,calc_info_.WABB,
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccC,calc_info_.nvirC,calc_info_.noccB,
    calc_info_.nvirB);
  fprintf(outfile,"ind_102             = %18.12lf  H\n",results_.ind102);
  fflush(outfile);
  results_.ind012 = ind210_0(calc_info_.sB_C,calc_info_.sC_A,calc_info_.sC_B,
    calc_info_.WBTT,calc_info_.WBCC,calc_info_.WATT,calc_info_.WACC,
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC);
  fprintf(outfile,"ind_012             = %18.12lf  H\n\n",results_.ind012);
  fflush(outfile);
}

double SAPT3BN5::ind110_1(double **t_A, double **t_C, double **W_A, 
  double **W_C, int occB, int virB)
{
  double energy = 0.0;

  energy += 2.0*C_DDOT(occB*virB,&(t_A[0][0]),1,&(W_C[0][0]),1);
  energy += 2.0*C_DDOT(occB*virB,&(t_C[0][0]),1,&(W_A[0][0]),1);

  return(energy);
}

double SAPT3BN5::ind111_1(double **CA_C, double **CB_C, int AAfile, 
  char *AA_ints, int BBfile, char *BB_ints, int occA, int virA, int occB, 
  int virB)
{
  double energy = 0.0;

  double **ARBS = IJKL_ints(AAfile,AA_ints,occA*virA,BBfile,BB_ints,occB*virB);

  for (int a=0; a<occA; a++) { 
    for (int r=0; r<virA; r++) { 
      for (int b=0; b<occB; b++) { 
        for (int s=0; s<virB; s++) {
          int ar = a*virA + r;
          int bs = b*virB + s;
          energy += 16.0*ARBS[ar][bs]*CA_C[a][r]*CB_C[b][s];
  }}}}
  
  free_block(ARBS);
  
  return(energy);
}

double SAPT3BN5::ind210_0(double **CA_B, double **CB_C, double **CB_A, 
  double **WASS, double **WABB, double **WCSS, double **WCBB, int AAfile, 
  char *AA_ints, int BBfile, char *BB_ints, int occA, int virA, int occB, 
  int virB)
{
  double energy = ind111_1(CA_B,CB_C,AAfile,AA_ints,BBfile,BB_ints,
                           occA,virA,occB,virB);
  energy += ind210_1(CB_A, CB_C, WASS, WABB, occB, virB);
  energy += ind210_2(CB_A, WCSS, WCBB, occB, virB);
  return(energy);
}

double SAPT3BN5::ind210_1(double **CB_A, double **CB_C, double **WASS, 
  double **WABB, int occB, int virB)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);

  C_DGEMM('N','N',occB,virB,virB,1.0,&(CB_C[0][0]),virB,&(WASS[0][0]),
    virB,0.0,&(X_BS[0][0]),virB);

  C_DGEMM('N','N',occB,virB,occB,-1.0,&(WABB[0][0]),occB,&(CB_C[0][0]),
    virB,1.0,&(X_BS[0][0]),virB);

  energy += 4.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(CB_A[0][0]),1);

  free_block(X_BS);

  return(energy);
}

double SAPT3BN5::ind210_2(double **CB_A, double **WCSS, double **WCBB, 
  int occB, int virB)
{
  double energy = 0.0;
  double **X_BS = block_matrix(occB,virB);

  C_DGEMM('N','N',occB,virB,virB,1.0,&(CB_A[0][0]),virB,&(WCSS[0][0]),
    virB,0.0,&(X_BS[0][0]),virB);

  C_DGEMM('N','N',occB,virB,occB,-1.0,&(WCBB[0][0]),occB,&(CB_A[0][0]),
    virB,1.0,&(X_BS[0][0]),virB);

  energy += 2.0*C_DDOT(occB*virB,&(X_BS[0][0]),1,&(CB_A[0][0]),1);

  free_block(X_BS);

  return(energy);
}

}}
