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

void SAPT3BN5::ind_disp210()
{
  results_.ind_disp210 = ind_disp210_0(calc_info_.sA_C,
    "T(BS) AR Intermediates",PSIF_3B_SAPT_AA_DF_INTS,"AA RI Integrals",
    "RR RI Integrals","T2 ARBS Amplitudes",'T',calc_info_.WCRR,calc_info_.WCAA,
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"ind_disp_210        = %18.12lf  H\n",results_.ind_disp210);
  fflush(outfile);
  results_.ind_disp201 = ind_disp210_0(calc_info_.sB_C,
    "T(AR) BS Intermediates",PSIF_3B_SAPT_BB_DF_INTS,"BB RI Integrals",
    "SS RI Integrals","T2 ARBS Amplitudes",'N',calc_info_.WCSS,calc_info_.WCBB,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"ind_disp_201        = %18.12lf  H\n",results_.ind_disp201);
  fflush(outfile);
  results_.ind_disp120 = ind_disp210_0(calc_info_.sA_B,
    "T(CT) AR Intermediates",PSIF_3B_SAPT_AA_DF_INTS,"AA RI Integrals",
    "RR RI Integrals","T2 ARCT Amplitudes",'T',calc_info_.WBRR,calc_info_.WBAA,
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"ind_disp_120        = %18.12lf  H\n",results_.ind_disp120);
  fflush(outfile);
  results_.ind_disp021 = ind_disp210_0(calc_info_.sC_B,
    "T(AR) CT Intermediates",PSIF_3B_SAPT_CC_DF_INTS,"CC RI Integrals",
    "TT RI Integrals","T2 ARCT Amplitudes",'N',calc_info_.WBTT,calc_info_.WBCC,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"ind_disp_021        = %18.12lf  H\n",results_.ind_disp021);
  fflush(outfile);
  results_.ind_disp102 = ind_disp210_0(calc_info_.sB_A,
    "T(CT) BS Intermediates",PSIF_3B_SAPT_BB_DF_INTS,"BB RI Integrals",
    "SS RI Integrals","T2 BSCT Amplitudes",'T',calc_info_.WASS,calc_info_.WABB,
    calc_info_.noccC,calc_info_.nvirC,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"ind_disp_102        = %18.12lf  H\n",results_.ind_disp102);
  fflush(outfile);
  results_.ind_disp012 = ind_disp210_0(calc_info_.sC_A,
    "T(BS) CT Intermediates",PSIF_3B_SAPT_CC_DF_INTS,"CC RI Integrals",
    "TT RI Integrals","T2 BSCT Amplitudes",'N',calc_info_.WATT,calc_info_.WACC,
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"ind_disp_012        = %18.12lf  H\n\n",
    results_.ind_disp012);
  fflush(outfile);
}

double SAPT3BN5::ind_disp210_0(double **CB_C, char *TAR_BS, int BBfile,
  char *BB_ints, char *SS_ints, char *T2label, char trans, double **WCSS,
  double **WCBB, int occA, int virA, int occB, int virB)
{
  double energy = 0.0;
  energy += ind_disp210_1(CB_C,TAR_BS,BBfile,BB_ints,SS_ints,occB,virB);
  energy += ind_disp210_2(T2label,trans,WCSS,WCBB,occA,virA,occB,virB);
  return(energy);
}

double SAPT3BN5::ind_disp210_1(double **CB_C, char *TAR_BS, int BBfile,
  char *BB_ints, char *SS_ints, int occB, int virB)
{
  double energy = 0.0;
  double **B_p_BS = block_matrix(occB*virB,calc_info_.nrio);
  double **B_p_SS = get_DF_ints(BBfile,SS_ints,virB*virB);

  C_DGEMM('N','N',occB,virB*calc_info_.nrio,virB,1.0,&(CB_C[0][0]),
          virB,&(B_p_SS[0][0]),virB*calc_info_.nrio,0.0,&(B_p_BS[0][0]),
          virB*calc_info_.nrio);

  free_block(B_p_SS);

  double **B_p_BB = get_DF_ints(BBfile,BB_ints,occB*occB);

  for (int i=0; i<occB; i++) {
    C_DGEMM('T','N',virB,calc_info_.nrio,occB,-1.0,&(CB_C[0][0]),
            virB,&(B_p_BB[i*occB][0]),calc_info_.nrio,1.0,
            &(B_p_BS[i*virB][0]),calc_info_.nrio);
  }

  free_block(B_p_BB);

  double **T_BS = read_IJKL(PSIF_3B_SAPT_AMPS,TAR_BS,occB*virB,
    calc_info_.nrio);

  energy = 8.0*C_DDOT(occB*virB*calc_info_.nrio,&(T_BS[0][0]),1,
                      &(B_p_BS[0][0]),1);

  free_block(T_BS);
  free_block(B_p_BS);

  return(energy);
}

double SAPT3BN5::ind_disp210_2(char *T2label, char trans, double **WCSS,
  double **WCBB, int occA, int virA, int occB, int virB)
{
  double energy = 0.0;

  if (trans == 'N') {
    double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virA*occA,occB*virB);
    double **xARBS = block_matrix(virA*occA,occB*virB);
    C_DGEMM('N','N',occA*virA*occB,virB,virB,1.0,&(tARBS[0][0]),virB,
            &(WCSS[0][0]),virB,0.0,&(xARBS[0][0]),virB);
    for (int i=0,ij=0;i<occA;i++) {
      for (int j=0;j<virA;j++,ij++) {
        C_DGEMM('N','N',occB,virB,occB,-1.0,&(WCBB[0][0]),occB,
                &(tARBS[ij][0]),virB,1.0,&(xARBS[ij][0]),virB);
    }}
    energy += 4.0*C_DDOT(occA*virA*occB*virB,&(tARBS[0][0]),1,&(xARBS[0][0]),
                         1);
    free_block(tARBS);
    free_block(xARBS);
  }
  else {
    double **tBSAR = read_IJKL(PSIF_3B_SAPT_AMPS,T2label,virB*occB,occA*virA);
    double **xBSAR = block_matrix(virB*occB,occA*virA);
    C_DGEMM('N','N',occB,virB*occA*virA,occB,-1.0,&(WCBB[0][0]),occB,
            &(tBSAR[0][0]),virB*occA*virA,0.0,&(xBSAR[0][0]),virB*occA*virA);
    for (int i=0,ij=0;i<occB;i++) {
      C_DGEMM('N','N',virB,occA*virA,virB,1.0,&(WCSS[0][0]),virB,
            &(tBSAR[i*virB][0]),occA*virA,1.0,&(xBSAR[i*virB][0]),occA*virA);
    }
    energy += 4.0*C_DDOT(occA*virA*occB*virB,&(tBSAR[0][0]),1,&(xBSAR[0][0]),
                         1);
    free_block(tBSAR);
    free_block(xBSAR);
  }

  return(energy);
}

}}
