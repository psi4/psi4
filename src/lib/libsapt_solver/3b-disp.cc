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
#include "sapt3bn6.h"

namespace psi { namespace sapt {

void SAPT3BN5::disp111()
{
  results_.disp111 = disp111_1("T(BS) AR Intermediates",
    "T(CT) AR Intermediates",calc_info_.noccA,calc_info_.nvirA);
  results_.disp111 += disp111_1("T(CT) BS Intermediates",
    "T(AR) BS Intermediates",calc_info_.noccB,calc_info_.nvirB);
  results_.disp111 += disp111_1("T(AR) CT Intermediates",
    "T(BS) CT Intermediates",calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_111            = %18.12lf  H\n\n",results_.disp111);
  fflush(outfile);
}

void SAPT3BN6::disp3100()
{ 
  results_.disp3100 = disp3100_1("T(CT) AR Intermediates",
    "T(CT) BS Intermediates","Theta(AR) AR Intermediates",
    "T2 ARBS Amplitudes",'N',"G ARAR",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.evalsA,calc_info_.evalsB,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  results_.disp3100 += disp3100_1("T(BS) AR Intermediates",
    "T(BS) CT Intermediates","Theta(AR) AR Intermediates",
    "T2 ARCT Amplitudes",'N',"G ARAR",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    calc_info_.evalsA,calc_info_.evalsC,calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccC,calc_info_.nvirC);
  results_.disp3100 += disp3100_2("Theta ARAR Amplitudes",
    "T(BS) AR Intermediates","T(CT) AR Intermediates",calc_info_.noccA,
    calc_info_.nvirA);
  results_.disp3100 += disp3100_3("T2 ARBS Amplitudes",'N',
    "T2 BSCT Amplitudes",'N',"T2 ARCT Amplitudes",'N',PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  results_.disp3100 += disp3100_4("T2 ARBS Amplitudes",
    "Theta(AR) AR Intermediates","T(CT) BS Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  results_.disp3100 += disp3100_4("T2 ARCT Amplitudes",
    "Theta(AR) AR Intermediates","T(BS) CT Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_111_100        = %18.12lf  H\n",results_.disp3100);
  fflush(outfile);

  results_.disp3010 = disp3100_1("T(AR) BS Intermediates",
    "T(AR) CT Intermediates","Theta(BS) BS Intermediates","T2 BSCT Amplitudes",
    'N',"G BSBS",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",calc_info_.evalsB,
    calc_info_.evalsC,calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC);
  results_.disp3010 += disp3100_1("T(CT) BS Intermediates",
    "T(CT) AR Intermediates","Theta(BS) BS Intermediates","T2 ARBS Amplitudes",
    'T',"G BSBS",PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.evalsB,
    calc_info_.evalsA,calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,
    calc_info_.nvirA);
  results_.disp3010 += disp3100_2("Theta BSBS Amplitudes",
    "T(CT) BS Intermediates","T(AR) BS Intermediates",calc_info_.noccB,
    calc_info_.nvirB);
  results_.disp3010 += disp3100_3("T2 BSCT Amplitudes",'N',
    "T2 ARCT Amplitudes",'T',"T2 ARBS Amplitudes",'T',PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC,calc_info_.noccA,calc_info_.nvirA);
  results_.disp3010 += disp3100_4("T2 BSCT Amplitudes",
    "Theta(BS) BS Intermediates","T(AR) CT Intermediates",calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  results_.disp3010 += disp3100_4("T2 ARBS Amplitudes",
    "T(CT) AR Intermediates","Theta(BS) BS Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"disp_111_010        = %18.12lf  H\n",results_.disp3010);
  fflush(outfile);

  results_.disp3001 = disp3100_1("T(BS) CT Intermediates",
    "T(BS) AR Intermediates","Theta(CT) CT Intermediates","T2 ARCT Amplitudes",
    'T',"G CTCT",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    PSIF_3B_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.evalsC,
    calc_info_.evalsA,calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,
    calc_info_.nvirA);
  results_.disp3001 += disp3100_1("T(AR) CT Intermediates",
    "T(AR) BS Intermediates","Theta(CT) CT Intermediates","T2 BSCT Amplitudes",
    'T',"G CTCT",PSIF_3B_SAPT_CC_DF_INTS,"CT RI Integrals",
    PSIF_3B_SAPT_BB_DF_INTS,"BS RI Integrals",calc_info_.evalsC,
    calc_info_.evalsB,calc_info_.noccC,calc_info_.nvirC,calc_info_.noccB,
    calc_info_.nvirB);
  results_.disp3001 += disp3100_2("Theta CTCT Amplitudes",
    "T(AR) CT Intermediates","T(BS) CT Intermediates",calc_info_.noccC,
    calc_info_.nvirC);
  results_.disp3001 += disp3100_3("T2 ARCT Amplitudes",'T',
    "T2 ARBS Amplitudes",'N',"T2 BSCT Amplitudes",'T',PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccC,calc_info_.nvirC,calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  results_.disp3001 += disp3100_4("T2 ARCT Amplitudes",
    "T(BS) AR Intermediates","Theta(CT) CT Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccC,calc_info_.nvirC);
  results_.disp3001 += disp3100_4("T2 BSCT Amplitudes",
    "T(AR) BS Intermediates","Theta(CT) CT Intermediates",calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_111_001        = %18.12lf  H\n\n",results_.disp3001);
  fflush(outfile);
}

void SAPT3BN6::disp211_D() 
{
  results_.disp211d = disp211_D_1("K ARBS Amplitudes",
    "K tilde ARBS Amplitudes",calc_info_.evalsA,calc_info_.evalsB,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"disp_211(D)         = %18.12lf  H\n",results_.disp211d);
  fflush(outfile);
  results_.disp121d = disp211_D_1("K ARCT Amplitudes",
    "K tilde ARCT Amplitudes",calc_info_.evalsA,calc_info_.evalsC,
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_121(D)         = %18.12lf  H\n",results_.disp121d);
  fflush(outfile);
  results_.disp112d = disp211_D_1("K BSCT Amplitudes",
    "K tilde BSCT Amplitudes",calc_info_.evalsB,calc_info_.evalsC,
    calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_112(D)         = %18.12lf  H\n\n",results_.disp112d);
  fflush(outfile);
}

void SAPT3BN6::disp220_S()
{
  results_.disp220s = disp220_S_1("S(BS) AR Amplitudes","S(CT) AR Amplitudes",
    calc_info_.evalsA,calc_info_.noccA,calc_info_.nvirA);
  fprintf(outfile,"disp_220(S)         = %18.12lf  H\n",results_.disp220s);
  fflush(outfile);
  results_.disp202s = disp220_S_1("S(AR) BS Amplitudes","S(CT) BS Amplitudes",
    calc_info_.evalsB,calc_info_.noccB,calc_info_.nvirB);
  fprintf(outfile,"disp_202(S)         = %18.12lf  H\n",results_.disp202s);
  fflush(outfile);
  results_.disp022s = disp220_S_1("S(AR) CT Amplitudes","S(BS) CT Amplitudes",
    calc_info_.evalsC,calc_info_.noccC,calc_info_.nvirC);
  fprintf(outfile,"disp_022(S)         = %18.12lf  H\n\n",results_.disp022s);
  fflush(outfile);
}

void SAPT3BN6::disp220_D()
{
  results_.disp220d = disp220_D_1("K(BS) ARAR Amplitudes",
    "K(CT) ARAR Amplitudes",calc_info_.evalsA,calc_info_.noccA,
    calc_info_.nvirA);
  fprintf(outfile,"disp_220(D)         = %18.12lf  H\n",results_.disp220d);
  fflush(outfile);
  results_.disp202d = disp220_D_1("K(AR) BSBS Amplitudes",
    "K(CT) BSBS Amplitudes",calc_info_.evalsB,calc_info_.noccB,
    calc_info_.nvirB);
  fprintf(outfile,"disp_202(D)         = %18.12lf  H\n",results_.disp202d);
  fflush(outfile);
  results_.disp022d = disp220_D_1("K(AR) CTCT Amplitudes",
    "K(BS) CTCT Amplitudes",calc_info_.evalsC,calc_info_.noccC,
    calc_info_.nvirC);
  fprintf(outfile,"disp_022(D)         = %18.12lf  H\n\n",results_.disp022d);
  fflush(outfile);
}

void SAPT3BN6::disp220_Q()
{
  results_.disp220q = disp220_Q_1("T2 ARCT Amplitudes",
    "T(BS) AR Intermediates","T(AR) CT Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccC,calc_info_.nvirC);
  results_.disp220q += disp220_Q_1("T2 ARBS Amplitudes",
    "T(CT) AR Intermediates","T(AR) BS Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  results_.disp220q += disp220_Q_3("T2 ARCT Amplitudes",
    "T(BS) AR Intermediates",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccC,
    calc_info_.nvirC);
  results_.disp220q += disp220_Q_3("T2 ARBS Amplitudes",
    "T(CT) AR Intermediates",PSIF_3B_SAPT_AA_DF_INTS,
    "AR RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB);
  fprintf(outfile,"disp_220(Q)         = %18.12lf  H\n",results_.disp220q);
  fflush(outfile);

  results_.disp202q = disp220_Q_1("T2 ARBS Amplitudes",
    "T(BS) AR Intermediates","T(CT) BS Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);
  results_.disp202q += disp220_Q_1("T2 BSCT Amplitudes",
    "T(AR) BS Intermediates","T(BS) CT Intermediates",calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  results_.disp202q += disp220_Q_2("T2 ARBS Amplitudes",
    "T(CT) BS Intermediates",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB);
  results_.disp202q += disp220_Q_3("T2 BSCT Amplitudes",
    "T(AR) BS Intermediates",PSIF_3B_SAPT_BB_DF_INTS,
    "BS RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC);
  fprintf(outfile,"disp_202(Q)         = %18.12lf  H\n",results_.disp202q);
  fflush(outfile);

  results_.disp022q = disp220_Q_1("T2 BSCT Amplitudes",
    "T(CT) BS Intermediates","T(AR) CT Intermediates",calc_info_.noccB,
    calc_info_.nvirB,calc_info_.noccC,calc_info_.nvirC);
  results_.disp022q += disp220_Q_1("T2 ARCT Amplitudes",
    "T(CT) AR Intermediates","T(BS) CT Intermediates",calc_info_.noccA,
    calc_info_.nvirA,calc_info_.noccC,calc_info_.nvirC);
  results_.disp022q += disp220_Q_2("T2 BSCT Amplitudes",
    "T(AR) CT Intermediates",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccB,calc_info_.nvirB,calc_info_.noccC,
    calc_info_.nvirC);
  results_.disp022q += disp220_Q_2("T2 ARCT Amplitudes",
    "T(BS) CT Intermediates",PSIF_3B_SAPT_CC_DF_INTS,
    "CT RI Integrals",calc_info_.noccA,calc_info_.nvirA,calc_info_.noccC,
    calc_info_.nvirC);
  fprintf(outfile,"disp_022(Q)         = %18.12lf  H\n\n",results_.disp022q);
  fflush(outfile);
}

double SAPT3BN5::disp111_1(char *AA_ints, char *BB_ints, int occA, int virA)
{
  double energy = 0.0;

  double **T_BS = read_IJKL(PSIF_3B_SAPT_AMPS,AA_ints,occA*virA,
    calc_info_.nrio);
  double **T_CT = read_IJKL(PSIF_3B_SAPT_AMPS,BB_ints,occA*virA,
    calc_info_.nrio);

  energy = 16.0*C_DDOT(occA*virA*calc_info_.nrio,&(T_BS[0][0]),1,
    &(T_CT[0][0]),1);

  free_block(T_BS);
  free_block(T_CT);

  return(energy);
}

double SAPT3BN6::disp3100_1(char *T_1, char *T_2, char *Theta, char *t2, 
  char trans, char *g2, int AAfile, char *AR_ints, int BBfile, char *BS_ints,
  double *e_A, double *e_B, int occA, int virA, int occB, int virB)
{
  double energy = 0.0;

  double **Y_ARBS = block_matrix(occA*virA,occB*virB);
  double **Th_AR = read_IJKL(PSIF_3B_SAPT_AMPS,Theta,occA*virA,
    calc_info_.nrio);
  double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);

  C_DGEMM('N','T',occA*virA,occB*virB,calc_info_.nrio,1.0,&(Th_AR[0][0]),
          calc_info_.nrio,&(B_p_BS[0][0]),calc_info_.nrio,0.0,&(Y_ARBS[0][0]),
          occB*virB);
  
  free_block(Th_AR);
  free_block(B_p_BS);
  
  double **gARAR = read_IJKL(PSIF_3B_SAPT_AMPS,g2,occA*virA,occA*virA);
  double **tARBS;
  
  if (trans == 'N') {
    tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2,occA*virA,occB*virB);
    C_DGEMM('N','N',occA*virA,occB*virB,occA*virA,1.0,&(gARAR[0][0]),
            occA*virA,&(tARBS[0][0]),occB*virB,1.0,&(Y_ARBS[0][0]),occB*virB);
  }
  else {
    tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2,occB*virB,occA*virA);
    C_DGEMM('N','T',occA*virA,occB*virB,occA*virA,1.0,&(gARAR[0][0]),
            occA*virA,&(tARBS[0][0]),occA*virA,1.0,&(Y_ARBS[0][0]),occB*virB);
  }
  
  free_block(tARBS);
  free_block(gARAR);
  
  double **X_ARBS = block_matrix(occA*virA,occB*virB);
  double **T_AR = read_IJKL(PSIF_3B_SAPT_AMPS,T_1,occA*virA,calc_info_.nrio);
  B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);

  C_DGEMM('N','T',occA*virA,occB*virB,calc_info_.nrio,1.0,&(T_AR[0][0]),
          calc_info_.nrio,&(B_p_BS[0][0]),calc_info_.nrio,0.0,&(X_ARBS[0][0]),
          occB*virB);

  free_block(T_AR);
  free_block(B_p_BS);

  double **T_BS = read_IJKL(PSIF_3B_SAPT_AMPS,T_2,occB*virB,calc_info_.nrio);
  double **B_p_AR = get_DF_ints(AAfile,AR_ints,occA*virA);

  C_DGEMM('N','T',occA*virA,occB*virB,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
          calc_info_.nrio,&(T_BS[0][0]),calc_info_.nrio,1.0,&(X_ARBS[0][0]),
          occB*virB);

  free_block(T_BS);
  free_block(B_p_AR);

  double denom,tval;

  for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      for (int b=0,bs=0; b<occB; b++) {
        for (int s=0; s<virB; s++,bs++) {
          denom = e_A[a] + e_B[b] - e_A[r+occA] - e_B[s+occB];
          tval = X_ARBS[ar][bs]*Y_ARBS[ar][bs];
          energy += tval/denom;
  }}}}

  free_block(X_ARBS);
  free_block(Y_ARBS);

  return(16.0*energy);
}

double SAPT3BN6::disp3100_2(char *theta, char *T_1, char *T_2, int occA, 
  int virA)
{
  double energy = 0.0;

  double **X_ARAR = block_matrix(occA*virA,occA*virA);
  double **T1_AR = read_IJKL(PSIF_3B_SAPT_AMPS,T_1,occA*virA,calc_info_.nrio);
  double **T2_AR = read_IJKL(PSIF_3B_SAPT_AMPS,T_2,occA*virA,calc_info_.nrio);

  C_DGEMM('N','T',occA*virA,occA*virA,calc_info_.nrio,1.0,&(T1_AR[0][0]),
          calc_info_.nrio,&(T2_AR[0][0]),calc_info_.nrio,0.0,&(X_ARAR[0][0]),
          occA*virA);

  free_block(T1_AR);
  free_block(T2_AR);

  double **thetaARAR = read_IJKL(PSIF_3B_SAPT_AMPS,theta,occA*virA,occA*virA);

  energy = 8.0*C_DDOT(occA*virA*occA*virA,&(X_ARAR[0][0]),1,
                      &(thetaARAR[0][0]),1);

  free_block(X_ARAR);
  free_block(thetaARAR);

  return(energy);
}

double SAPT3BN6::disp3100_3(char *t2a, char trans1, char *t2b, char trans2, 
  char *t2c, char trans3, int AAfile, char *AR_ints, int occA, int virA, 
  int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **X_ARCT = block_matrix(occA*virA,occC*virC);
  double **T_ARBS;
  double **T_BSCT;

  if (trans1 == 'N' && trans2 == 'N') {
    T_ARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2a,virA*occA,occB*virB);
    T_BSCT = read_IJKL(PSIF_3B_SAPT_AMPS,t2b,virB*occB,occC*virC);
    C_DGEMM('N','N',occA*virA,occC*virC,occB*virB,1.0,&(T_ARBS[0][0]),
            occB*virB,&(T_BSCT[0][0]),occC*virC,0.0,&(X_ARCT[0][0]),
            occC*virC);
  }
  else if (trans1 == 'N' && trans2 == 'T') {
    T_ARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2a,virA*occA,occB*virB);
    T_BSCT = read_IJKL(PSIF_3B_SAPT_AMPS,t2b,virC*occC,occB*virB);
    C_DGEMM('N','T',occA*virA,occC*virC,occB*virB,1.0,&(T_ARBS[0][0]),
            occB*virB,&(T_BSCT[0][0]),occB*virB,0.0,&(X_ARCT[0][0]),
            occC*virC);
  }
  else if (trans1 == 'T' && trans2 == 'N') {
    T_ARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2a,virB*occB,occA*virA);
    T_BSCT = read_IJKL(PSIF_3B_SAPT_AMPS,t2b,virB*occB,occC*virC);
    C_DGEMM('T','N',occA*virA,occC*virC,occB*virB,1.0,&(T_ARBS[0][0]),
            occA*virA,&(T_BSCT[0][0]),occC*virC,0.0,&(X_ARCT[0][0]),
            occC*virC);
  }
  else {
    printf("You fucked up\n");
    exit(1);
  }

  free_block(T_ARBS);
  free_block(T_BSCT);

  double **X_ARAR = block_matrix(occA*virA,occA*virA);
  double **T_ARCT;

  if (trans3 == 'N') {
    T_ARCT = read_IJKL(PSIF_3B_SAPT_AMPS,t2c,virA*occA,occC*virC);
    C_DGEMM('N','T',occA*virA,occA*virA,occC*virC,1.0,&(X_ARCT[0][0]),
            occC*virC,&(T_ARCT[0][0]),occC*virC,0.0,&(X_ARAR[0][0]),
            occA*virA);
  }
  else {
    T_ARCT = read_IJKL(PSIF_3B_SAPT_AMPS,t2c,virC*occC,occA*virA);
    C_DGEMM('N','N',occA*virA,occA*virA,occC*virC,1.0,&(X_ARCT[0][0]),
            occC*virC,&(T_ARCT[0][0]),occA*virA,0.0,&(X_ARAR[0][0]),
            occA*virA);
  }

  free_block(X_ARCT);
  free_block(T_ARCT);

  double **ARAR = IJIJ_ints(AAfile,AR_ints,occA*virA);

  double tval;

  for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      for (int ap=0,aprp=0; ap<occA; ap++) {
        for (int rp=0; rp<virA; rp++,aprp++) {
          int arp = a*virA+rp;
          int apr = ap*virA+r;
          tval = 2.0*ARAR[ar][aprp] - ARAR[arp][apr];
          tval *= 8.0*X_ARAR[ar][aprp];
          energy += tval;
        }}
    }}

  free_block(X_ARAR);
  free_block(ARAR);

  return(energy);
}

double SAPT3BN6::disp3100_4(char *t2, char *T_1, char *T_2, int occA, int virA,
  int occB, int virB)
{
  double energy = 0.0;

  double **X_ARBS = block_matrix(occA*virA,occB*virB);
  double **T1_AR = read_IJKL(PSIF_3B_SAPT_AMPS,T_1,occA*virA,calc_info_.nrio);
  double **T2_BS = read_IJKL(PSIF_3B_SAPT_AMPS,T_2,occB*virB,calc_info_.nrio);

  C_DGEMM('N','T',occA*virA,occB*virB,calc_info_.nrio,1.0,&(T1_AR[0][0]),
          calc_info_.nrio,&(T2_BS[0][0]),calc_info_.nrio,0.0,&(X_ARBS[0][0]),
          occB*virB);

  free_block(T1_AR);
  free_block(T2_BS);

  double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2,occA*virA,occB*virB);

  energy = 8.0*C_DDOT(occA*virA*occB*virB,&(X_ARBS[0][0]),1,
                      &(tARBS[0][0]),1);

  free_block(X_ARBS);
  free_block(tARBS);

  return(energy);
}

double SAPT3BN6::disp211_D_1(char *K_1, char *K_2, double *e_A, double *e_B,
  int occA, int virA, int occB, int virB)
{
  double energy = 0.0;
  
  double **K1_ARBS = read_IJKL(PSIF_3B_SAPT_AMPS,K_1,occA*virA,occB*virB);
  double **K2_ARBS = read_IJKL(PSIF_3B_SAPT_AMPS,K_2,occA*virA,occB*virB);
  
  double denom,tval;
  
  for (int a=0,ar=0; a<occA; a++) { 
    for (int r=0; r<virA; r++,ar++) { 
      for (int b=0,bs=0; b<occB; b++) { 
        for (int s=0; s<virB; s++,bs++) {
          denom = e_A[a] + e_B[b] - e_A[r+occA] - e_B[s+occB];
          tval = K1_ARBS[ar][bs]*K2_ARBS[ar][bs];
          tval += K2_ARBS[ar][bs]*K2_ARBS[ar][bs];
          energy += tval/denom;
  }}}}
  
  free_block(K1_ARBS);
  free_block(K2_ARBS);
  
  return(16.0*energy);
}

double SAPT3BN6::disp220_S_1(char *BS_S, char *CT_S, double *e_A, int occA, 
  int virA)
{
  double energy = 0.0;

  double **S_BS = read_IJKL(PSIF_3B_SAPT_AMPS,BS_S,occA,virA);
  double **S_CT = read_IJKL(PSIF_3B_SAPT_AMPS,CT_S,occA,virA);

  double denom;

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      denom = e_A[a] - e_A[r+occA];
      energy += S_BS[a][r]*S_CT[a][r]/denom;
  }}

  free_block(S_BS);
  free_block(S_CT);

  return(16.0*energy);
}

double SAPT3BN6::disp220_D_1(char *K_1, char *K_2, double *e_A, int occA, 
  int virA)
{
  double energy = 0.0;

  double **K1_ARBS = read_IJKL(PSIF_3B_SAPT_AMPS,K_1,occA*virA,occA*virA);
  double **K2_ARBS = read_IJKL(PSIF_3B_SAPT_AMPS,K_2,occA*virA,occA*virA);

  double denom,tval;

  for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      for (int ap=0,aprp=0; ap<occA; ap++) {
        for (int rp=0; rp<virA; rp++,aprp++) {
          denom = e_A[a] + e_A[ap] - e_A[r+occA] - e_A[rp+occA];
          energy += K1_ARBS[ar][aprp]*K2_ARBS[ar][aprp]/denom;
  }}}}

  energy *= 16.0;

  for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      for (int ap=0,aprp=0; ap<occA; ap++) {
        for (int rp=0; rp<virA; rp++,aprp++) {
          int arp = a*virA+rp;
          int apr = ap*virA+r;
          denom = e_A[a] + e_A[ap] - e_A[r+occA] - e_A[rp+occA];
          tval = K1_ARBS[ar][aprp]*K2_ARBS[apr][arp];
          tval += K1_ARBS[apr][arp]*K2_ARBS[ar][aprp];
          energy -= 4.0*tval/denom;
  }}}}
    
  free_block(K1_ARBS);
  free_block(K2_ARBS);

  return(energy);
}

double SAPT3BN6::disp220_Q_1(char *t2, char *T_1, char *T_2, int occA, int virA,
  int occB, int virB)
{
  double energy = 0.0;

  double **X_ARBS = block_matrix(occA*virA,occB*virB);
  double **T_AR = read_IJKL(PSIF_3B_SAPT_AMPS,T_1,occA*virA,calc_info_.nrio);
  double **T_BS = read_IJKL(PSIF_3B_SAPT_AMPS,T_2,occB*virB,calc_info_.nrio);

  C_DGEMM('N','T',occA*virA,occB*virB,calc_info_.nrio,1.0,&(T_AR[0][0]),
          calc_info_.nrio,&(T_BS[0][0]),calc_info_.nrio,0.0,&(X_ARBS[0][0]),
          occB*virB);

  free_block(T_AR);
  free_block(T_BS);

  double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2,occA*virA,occB*virB);

  energy = 16.0*C_DDOT(occA*virA*occB*virB,&(X_ARBS[0][0]),1,&(tARBS[0][0]),1);

  free_block(X_ARBS);
  free_block(tARBS);

  return(energy);
}

double SAPT3BN6::disp220_Q_2(char *t2, char *T_1, int BBfile, char *BS_ints, 
  int occA, int virA, int occB, int virB)
{
  double energy = 0.0;

  double **X_BB = block_matrix(occB,occB);
  double **X_SS = block_matrix(virB,virB);

  double **T_BS = read_IJKL(PSIF_3B_SAPT_AMPS,T_1,occB*virB,calc_info_.nrio);
  double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);

  C_DGEMM('N','T',occB,occB,virB*calc_info_.nrio,1.0,&(T_BS[0][0]),
          virB*calc_info_.nrio,&(B_p_BS[0][0]),virB*calc_info_.nrio,0.0,
          &(X_BB[0][0]),occB);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','T',virB,virB,calc_info_.nrio,1.0,&(B_p_BS[b*virB][0]),
            calc_info_.nrio,&(T_BS[b*virB][0]),calc_info_.nrio,1.0,
            &(X_SS[0][0]),virB);
  }

  free_block(T_BS);
  free_block(B_p_BS);

  double **tARBS = read_IJKL(PSIF_3B_SAPT_AMPS,t2,occA*virA,occB*virB);
  double **X_ARBS = block_matrix(occA*virA,occB*virB);

  for (int a=0,ar=0; a<occA; a++) {
    for (int r=0; r<virA; r++,ar++) {
      C_DGEMM('N','N',occB,virB,occB,-1.0,&(X_BB[0][0]),occB,&(tARBS[ar][0]),
              virB,1.0,&(X_ARBS[ar][0]),virB);
  }}

  C_DGEMM('N','N',occA*virA*occB,virB,virB,-1.0,&(tARBS[0][0]),virB,
          &(X_SS[0][0]),virB,1.0,&(X_ARBS[0][0]),virB);

  energy = 8.0*C_DDOT(occA*virA*occB*virB,&(X_ARBS[0][0]),1,&(tARBS[0][0]),1);

  free_block(X_BB);
  free_block(X_SS);
  free_block(X_ARBS);
  free_block(tARBS);

  return(energy);
}

double SAPT3BN6::disp220_Q_3(char *t2, char *T_1, int BBfile, char *BS_ints, 
  int occB, int virB, int occC, int virC)
{
  double energy = 0.0;

  double **X_BB = block_matrix(occB,occB);
  double **X_SS = block_matrix(virB,virB);

  double **T_BS = read_IJKL(PSIF_3B_SAPT_AMPS,T_1,occB*virB,calc_info_.nrio);
  double **B_p_BS = get_DF_ints(BBfile,BS_ints,occB*virB);

  C_DGEMM('N','T',occB,occB,virB*calc_info_.nrio,1.0,&(T_BS[0][0]),
          virB*calc_info_.nrio,&(B_p_BS[0][0]),virB*calc_info_.nrio,0.0,
          &(X_BB[0][0]),occB);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','T',virB,virB,calc_info_.nrio,1.0,&(T_BS[b*virB][0]),
            calc_info_.nrio,&(B_p_BS[b*virB][0]),calc_info_.nrio,1.0,
            &(X_SS[0][0]),virB);
  }

  free_block(T_BS);
  free_block(B_p_BS);

  double **tBSCT = read_IJKL(PSIF_3B_SAPT_AMPS,t2,occB*virB,occC*virC);
  double **X_BSCT = block_matrix(occB*virB,occC*virC);

  C_DGEMM('N','N',occB,virB*occC*virC,occB,-1.0,&(X_BB[0][0]),occB,
          &(tBSCT[0][0]),virB*occC*virC,0.0,&(X_BSCT[0][0]),virB*occC*virC);

  for (int b=0; b<occB; b++) {
    C_DGEMM('N','N',virB,occC*virC,virB,-1.0,&(X_SS[0][0]),virB,
            &(tBSCT[b*virB][0]),occC*virC,1.0,&(X_BSCT[b*virB][0]),occC*virC);
  }

  energy = 8.0*C_DDOT(occB*virB*occC*virC,&(X_BSCT[0][0]),1,&(tBSCT[0][0]),1);

  free_block(X_BB);
  free_block(X_SS);
  free_block(X_BSCT);
  free_block(tBSCT);

  return(energy);
}

}}
