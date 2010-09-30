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

void SAPT2p3::exch_ind_disp30()
{
  double e10,e01,e11,e21,e12;

  if (params_.print) 
    fprintf(outfile,"Begining Exch-Ind-Disp30 Calculation\n\n");

  double **X_AR = read_IJKL(PSIF_SAPT_AMPS,"Ind-Disp30 AR Amplitudes",
    calc_info_.noccA,calc_info_.nvirA);
  double **V_AR = read_IJKL(PSIF_SAPT_AMPS,"AR Exch-Ind Integrals",
    calc_info_.noccA,calc_info_.nvirA);

  e10 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,X_AR[0],1,V_AR[0],1);

  free_block(X_AR);
  free_block(V_AR);

  double **X_BS = read_IJKL(PSIF_SAPT_AMPS,"Ind-Disp30 BS Amplitudes",
    calc_info_.noccB,calc_info_.nvirB);
  double **V_BS = read_IJKL(PSIF_SAPT_AMPS,"BS Exch-Ind Integrals",
    calc_info_.noccB,calc_info_.nvirB);

  e01 = -2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,X_BS[0],1,V_BS[0],1);

  free_block(X_BS);
  free_block(V_BS);

  double **X_ARBS = read_IJKL(PSIF_SAPT_AMPS,"Ind-Disp30 ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);
  double **V_ARBS = read_IJKL(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  e11 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nvirB*
    calc_info_.noccB,&(V_ARBS[0][0]),1,&(X_ARBS[0][0]),1);

  free_block(X_ARBS);
  free_block(V_ARBS);

  e21 = exch_ind_disp30_21();
  e12 = exch_ind_disp30_12();

  results_.exch_ind_disp30 = e10+e01+e11;
  if (params_.print) {
    fprintf(outfile,"Exch-Ind-Disp30    = %18.12lf  H\n\n",
      results_.exch_ind_disp30);
    fflush(outfile);
  }
}

double SAPT2p3::exch_ind_disp30_21()
{
  double energy = 0.0;

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);

  double **tAS_RB = block_matrix(calc_info_.nvirA,calc_info_.noccB);
  double **tRB_AS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      for (int b=0,bs=0; b<calc_info_.noccB; b++) {
        for (int s=0; s<calc_info_.nvirB; s++,bs++) {
          tAS_RB[r][b] += tARBS[ar][bs]*calc_info_.S_AB[a][s+calc_info_.noccB];
          tRB_AS[a][s] += tARBS[ar][bs]*calc_info_.S_AB[r+calc_info_.noccA][b];
  }}}}

  double **B_p_AR = get_AR_ints(1);
  double **B_p_RB = get_RB_ints(1);

  double **xRS = block_matrix(calc_info_.nvirA,calc_info_.nvirB);
  double **C_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(xRS[0][0]),
    calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(xRS[0][0]),calc_info_.nvirB,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  double **xRBS = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(B_p_RB[0][0]),calc_info_.nrio,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,&(xRBS[0][0]),
      calc_info_.nvirB);
    energy += C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      tARBS[a*calc_info_.nvirA],1,xRBS[0],1);
  }

  free_block(xRBS);
  free_block(C_p_AS);

  double **B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(tAS_RB[0][0]),calc_info_.noccB,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(B_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  energy += C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  free_block(B_p_AB);
  free_block(C_p_AB);

  double *X = init_array(calc_info_.nrio);
  double *Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.sA[0],1,0.0,X,1);
  
  C_DGEMV('t',calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_RB[0][0]),calc_info_.nrio,tAS_RB[0],1,0.0,Y,1);

  energy -= 2.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);

  double **C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(xRS[0][0]),calc_info_.nvirB,&(B_p_RB[b][0]),
      calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio);
  }

  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  energy -= 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    C_p_BS[0],1,T_p_BS[0],1);

  free_block(xRS);
  free_block(B_p_RB);
  free_block(C_p_BS);
  free_block(T_p_BS);

  double **B_p_BS = get_BS_ints(1);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  C_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);
  double **C_p_RB = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.nrio);

  for (int r=0; r<calc_info_.nvirA; r++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(xAB[0][0]),calc_info_.noccB,&(B_p_AR[r][0]),
      calc_info_.nvirA*calc_info_.nrio,0.0,&(C_p_RB[r*calc_info_.noccB][0]),
      calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_BS[0][0]),calc_info_.nvirB*calc_info_.nrio,0.0,&(C_p_AS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  xRBS = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,&(xRBS[0][0]),
      calc_info_.nvirB);
    energy -= C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      tARBS[a*calc_info_.nvirA],1,xRBS[0],1);
  }

  free_block(xAB);
  free_block(xRBS);
  free_block(C_p_AS);
  free_block(C_p_RB);

  double **xRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(xRB[0][0]),calc_info_.noccB);

  B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);
  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(xRB[0][0]),calc_info_.noccB,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(tRB_AS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(C_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }
    
  energy -= C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);
 
  free_block(xRB); 
  free_block(B_p_AB);
  free_block(C_p_AB);

  double **xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAR[0][0]),
    calc_info_.nvirA);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.nvirA,1.0,
    &(xAR[0][0]),calc_info_.nvirA,&(calc_info_.sA[0][0]),calc_info_.nvirA,
    0.0,&(xAA[0][0]),calc_info_.noccA);

  double **xRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nvirA,calc_info_.noccA,1.0,
    &(xAR[0][0]),calc_info_.nvirA,&(calc_info_.sA[0][0]),calc_info_.nvirA,
    0.0,&(xRR[0][0]),calc_info_.nvirA);
    
  double **C_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(xAA[0][0]),calc_info_.noccA,&(B_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(xRR[0][0]),calc_info_.nvirA,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,1.0,&(C_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio);
  }

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    C_p_AR[0],1,T_p_AR[0],1);

  free_block(xAR);
  free_block(xAA);
  free_block(xRR);
  free_block(C_p_AR);
  free_block(T_p_AR);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    0.0,&(xBB[0][0]),calc_info_.noccB);

  T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(xBB[0][0]),calc_info_.noccB,&(B_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio,0.0,&(C_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    C_p_BS[0],1,T_p_BS[0],1);

  free_block(xAB);
  free_block(xBB);
  free_block(C_p_BS);
  free_block(T_p_BS);

  double **xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(tRB_AS[0][0]),calc_info_.nvirB,
    0.0,&(xBS[0][0]),calc_info_.nvirB);

  X = init_array(calc_info_.nrio);
  Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.sA[0],1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,xBS[0],1,0.0,Y,1);

  energy += 2.0*C_DDOT(calc_info_.nrio,X,1,Y,1);
  
  free(X);
  free(Y);
  free_block(xBS);
  free_block(B_p_BS);

  double **B_p_BB = get_BB_ints(1);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(xAB[0][0]),calc_info_.noccB,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(B_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  xRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.noccB,calc_info_.nrio,
      1.0,&(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(xRB[0][0]),
      calc_info_.noccB);
  }

  energy -= C_DDOT(calc_info_.nvirA*calc_info_.noccB,xRB[0],1,tAS_RB[0],1);

  free_block(xAB);
  free_block(xRB);
  free_block(B_p_AB);

  xRS = block_matrix(calc_info_.nvirA,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(xRS[0][0]),
    calc_info_.nvirB);

  C_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio);
  C_p_RB = block_matrix(calc_info_.noccB*calc_info_.nvirA,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(xRS[0][0]),calc_info_.nvirB,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,&(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,
    &(C_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio);

  xRBS = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,&(xRBS[0][0]),
      calc_info_.nvirB);
    energy -= C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      tARBS[a*calc_info_.nvirA],1,xRBS[0],1);
  }

  free_block(xRS);
  free_block(xRBS);
  free_block(C_p_AS);
  free_block(C_p_RB);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(xBS[0][0]),calc_info_.nvirB);

  C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccB,1.0,
      &(xBS[0][0]),calc_info_.nvirB,&(B_p_BB[b*calc_info_.noccB][0]),
      calc_info_.nrio,0.0,&(C_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio);
  }

 T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
   calc_info_.nvirB,calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    C_p_BS[0],1,T_p_BS[0],1);

  free_block(xAB);
  free_block(xBS);
  free_block(C_p_BS);
  free_block(T_p_BS);

  X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.sA[0],1,0.0,X,1);

  xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,X,1,0.0,xBB[0],1);

  xRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(xBB[0][0]),calc_info_.noccB,0.0,&(xRB[0][0]),calc_info_.noccB);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      for (int b=0; b<calc_info_.noccB; b++) {
        energy += 2.0*xRB[r][b]*C_DDOT(calc_info_.nvirB,
          &(tARBS[ar][b*calc_info_.nvirB]),1,
          &(calc_info_.S_AB[a][calc_info_.noccB]),1);
  }}}

  free(X);
  free_block(xBB);
  free_block(xRB);

  xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,xAR[0],1);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  xRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccA,1.0,
    &(xAR[0][0]),calc_info_.nvirA,&(xAB[0][0]),calc_info_.noccB,0.0,
    &(xRB[0][0]),calc_info_.noccB);

  energy += 2.0*C_DDOT(calc_info_.nvirA*calc_info_.noccB,xRB[0],1,tAS_RB[0],1);

  free_block(xAB);
  free_block(xRB);

  xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,&(xAR[0][0]),calc_info_.nvirA,
    0.0,&(xAA[0][0]),calc_info_.noccA);

  double **xAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(xAA[0][0]),calc_info_.noccA,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(xAS[0][0]),calc_info_.nvirB);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirB,xAS[0],1,tRB_AS[0],1);

  free_block(xAR);
  free_block(xAA);
  free_block(xAS);
  free_block(B_p_BB);

  xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAR[0][0]),
    calc_info_.nvirA);

 T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
   calc_info_.nvirA,calc_info_.nrio);

  X = init_array(calc_info_.nrio);
  Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.sA[0],1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(T_p_AR[0][0]),calc_info_.nrio,xAR[0],1,0.0,Y,1);

  energy -= 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xAR);
  free_block(T_p_AR);

  xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,xAR[0],1);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(xBS[0][0]),calc_info_.nvirB);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      energy -= 4.0*xAR[a][r]*C_DDOT(calc_info_.noccB*calc_info_.nvirB,
        &(tARBS[ar][0]),1,&(xBS[0][0]),1);
  }}

  free_block(xAR);
  free_block(xAB);
  free_block(xBS);
  free_block(B_p_AR);

  return(2.0*energy);
}

double SAPT2p3::exch_ind_disp30_12()
{
  double energy = 0.0;

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);

  double **tAS_RB = block_matrix(calc_info_.nvirA,calc_info_.noccB);
  double **tRB_AS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      for (int b=0,bs=0; b<calc_info_.noccB; b++) {
        for (int s=0; s<calc_info_.nvirB; s++,bs++) {
          tAS_RB[r][b] += tARBS[ar][bs]*calc_info_.S_AB[a][s+calc_info_.noccB];
          tRB_AS[a][s] += tARBS[ar][bs]*calc_info_.S_AB[r+calc_info_.noccA][b];
  }}}}

  double **B_p_BS = get_BS_ints(1);
  double **B_p_AS = get_AS_ints(1);

  double **xRS = block_matrix(calc_info_.nvirA,calc_info_.nvirB);
  double **C_p_RB = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nrio);
  
  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xRS[0][0]),
    calc_info_.nvirB);
    
  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xRS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(C_p_RB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }

  double **xRBS = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,
      &(B_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,&(xRBS[0][0]),
      calc_info_.nvirB);
    energy += C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      tARBS[a*calc_info_.nvirA],1,xRBS[0],1);
  }

  free_block(xRBS);
  free_block(C_p_RB);

  double **B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio); 

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(tRB_AS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.sB[0][0]),calc_info_.nvirB,&(B_p_AS[a*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  energy += C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  free_block(B_p_AB);
  free_block(C_p_AB);

  double *X = init_array(calc_info_.nrio);
  double *Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.sB[0],1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_AS[0][0]),calc_info_.nrio,tRB_AS[0],1,0.0,Y,1);

  energy -= 2.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);

  double **C_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xRS[0][0]),calc_info_.nvirB,&(B_p_AS[a*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(C_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio);
  }

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    C_p_AR[0],1,T_p_AR[0],1);

  free_block(xRS);
  free_block(B_p_AS);
  free_block(C_p_AR);
  free_block(T_p_AR);

  double **B_p_AR = get_AR_ints(1);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  double **C_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);
  C_p_RB = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(xAB[0][0]),calc_info_.noccB,&(B_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio,0.0,&(C_p_AS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  for (int r=0; r<calc_info_.nvirA; r++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_AR[r][0]),
      calc_info_.nvirA*calc_info_.nrio,0.0,&(C_p_RB[r*calc_info_.noccB][0]),
      calc_info_.nrio);
  }

  xRBS = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,&(xRBS[0][0]),
      calc_info_.nvirB);
    energy -= C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      tARBS[a*calc_info_.nvirA],1,xRBS[0],1);
  }

  free_block(xAB);
  free_block(xRBS);
  free_block(C_p_AS);
  free_block(C_p_RB);

  double **xAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.sB[0][0]),
    calc_info_.nvirB,0.0,&(xAS[0][0]),calc_info_.nvirB);

  B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);
  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xAS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }
  
  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(tAS_RB[0][0]),calc_info_.noccB,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  } 

  energy -= C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);
    
  free_block(xAS); 
  free_block(B_p_AB);
  free_block(C_p_AB);

  double **xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(xBS[0][0]),
    calc_info_.nvirB);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(xBS[0][0]),calc_info_.nvirB,&(calc_info_.sB[0][0]),calc_info_.nvirB,
    0.0,&(xBB[0][0]),calc_info_.noccB);

  double **xSS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(xBS[0][0]),calc_info_.nvirB,&(calc_info_.sB[0][0]),calc_info_.nvirB,
    0.0,&(xSS[0][0]),calc_info_.nvirB);

  double **C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(xBB[0][0]),calc_info_.noccB,&(B_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio,0.0,&(C_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xSS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,1.0,&(C_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    C_p_BS[0],1,T_p_BS[0],1);

  free_block(xBS);
  free_block(xBB);
  free_block(xSS);
  free_block(C_p_BS);
  free_block(T_p_BS);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    0.0,&(xAA[0][0]),calc_info_.noccA);

  T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);

  C_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(xAA[0][0]),calc_info_.noccA,&(B_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    C_p_AR[0],1,T_p_AR[0],1);

  free_block(xAB);
  free_block(xAA);
  free_block(C_p_AR);
  free_block(T_p_AR);

  double **xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(tAS_RB[0][0]),calc_info_.noccB,
    0.0,&(xAR[0][0]),calc_info_.nvirA);

  X = init_array(calc_info_.nrio);
  Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.sB[0],1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,xAR[0],1,0.0,Y,1);

  energy += 2.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xAR);
  free_block(B_p_AR);

  double **B_p_AA = get_AA_ints(1);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  double **B_p_BA = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(xAB[0][0]),calc_info_.noccB,&(B_p_AA[0][0]),
    calc_info_.noccA*calc_info_.nrio,0.0,&(B_p_BA[0][0]),
    calc_info_.noccA*calc_info_.nrio);

  xAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirB,calc_info_.nrio,
      1.0,&(B_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(xAS[0][0]),
      calc_info_.nvirB);
  }

  energy -= C_DDOT(calc_info_.noccA*calc_info_.nvirB,xAS[0],1,tRB_AS[0],1);

  free_block(xAB);
  free_block(xAS);
  free_block(B_p_BA);

  xRS = block_matrix(calc_info_.nvirA,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xRS[0][0]),
    calc_info_.nvirB);

  C_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio);
  C_p_RB = block_matrix(calc_info_.noccB*calc_info_.nvirA,calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xRS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(C_p_RB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,0.0,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  xRBS = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,&(xRBS[0][0]),
      calc_info_.nvirB);
    energy -= C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      tARBS[a*calc_info_.nvirA],1,xRBS[0],1);
  }

  free_block(xRS);
  free_block(xRBS);
  free_block(C_p_AS);
  free_block(C_p_RB);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(xAR[0][0]),calc_info_.nvirA);

  C_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.noccA,1.0,
      &(xAR[0][0]),calc_info_.nvirA,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(C_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio);
  }

 T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
   calc_info_.nvirA,calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    C_p_AR[0],1,T_p_AR[0],1);

  free_block(xAB);
  free_block(xAR);
  free_block(C_p_AR);
  free_block(T_p_AR);

  X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.sB[0],1,0.0,X,1);

  xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,X,1,0.0,xAA[0],1);

  xAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(xAA[0][0]),calc_info_.noccA,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(xAS[0][0]),calc_info_.nvirB);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      for (int b=0; b<calc_info_.noccB; b++) {
        energy += 2.0*calc_info_.S_AB[r+calc_info_.noccA][b]*
          C_DDOT(calc_info_.nvirB,&(tARBS[ar][b*calc_info_.nvirB]),1,
          &(xAS[a][0]),1);
  }}}

  free(X);
  free_block(xAA);
  free_block(xAS);

  xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,xBS[0],1);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  xAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(xBS[0][0]),calc_info_.nvirB,0.0,
    &(xAS[0][0]),calc_info_.nvirB);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirB,xAS[0],1,tRB_AS[0],1);

  free_block(xAB);
  free_block(xAS);

  xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,&(xBS[0][0]),calc_info_.nvirB,
    0.0,&(xBB[0][0]),calc_info_.noccB);

  double **xRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,&(xBB[0][0]),
    calc_info_.noccB,0.0,&(xRB[0][0]),calc_info_.noccB);

  energy += 2.0*C_DDOT(calc_info_.nvirA*calc_info_.noccB,xRB[0],1,tAS_RB[0],1);

  free_block(xBS);
  free_block(xBB);
  free_block(xRB);
  free_block(B_p_AA);

  xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(xBS[0][0]),
    calc_info_.nvirB);

 T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
   calc_info_.nvirB,calc_info_.nrio);

  X = init_array(calc_info_.nrio);
  Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.sB[0],1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(T_p_BS[0][0]),calc_info_.nrio,xBS[0],1,0.0,Y,1);

  energy -= 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xBS);
  free_block(T_p_BS);

  xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,xBS[0],1);

  xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xAB[0][0]),calc_info_.noccB);

  xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(xAR[0][0]),calc_info_.nvirA);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      energy -= 4.0*xAR[a][r]*C_DDOT(calc_info_.noccB*calc_info_.nvirB,
        &(tARBS[ar][0]),1,&(xBS[0][0]),1);
  }}

  free_block(xAR);
  free_block(xAB);
  free_block(xBS);
  free_block(B_p_BS);

  return(2.0*energy);
}

}}

