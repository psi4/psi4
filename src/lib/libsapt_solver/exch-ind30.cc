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

void SAPT2p3::exch_ind30()
{
  double e10,e01,e11,e20,e02;

  if (params_.print) 
    fprintf(outfile,"Begining Exch-Ind30 Calculation\n\n");

  double **X_AR = read_IJKL(PSIF_SAPT_AMPS,"Ind30 AR Amplitudes",
    calc_info_.noccA,calc_info_.nvirA);
  double **V_AR = read_IJKL(PSIF_SAPT_AMPS,"AR Exch-Ind Integrals",
    calc_info_.noccA,calc_info_.nvirA);

  e10 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,X_AR[0],1,V_AR[0],1);

  free_block(V_AR);

  double **X_BS = read_IJKL(PSIF_SAPT_AMPS,"Ind30 BS Amplitudes",
    calc_info_.noccB,calc_info_.nvirB);
  double **V_BS = read_IJKL(PSIF_SAPT_AMPS,"BS Exch-Ind Integrals",
    calc_info_.noccB,calc_info_.nvirB);

  e01 = -2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,X_BS[0],1,V_BS[0],1);

  free_block(V_BS);

  e11 = 0.0;

  double **vARBS = read_IJKL(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      e11 -= 2.0*calc_info_.sA[a][r]*C_DDOT(calc_info_.noccB*calc_info_.nvirB,
        vARBS[ar],1,calc_info_.sB[0],1);
  }}

  free_block(vARBS);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(X_BS[0][0]),
    calc_info_.nvirB);

  double **B_p_AR = get_AR_ints(1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(X_AR[0][0]),1);

  free_block(B_p_AR);

  double t1 = C_DDOT(calc_info_.noccA*calc_info_.nvirA,X_AR[0],1,
    calc_info_.sA[0],1);
  double t2 = C_DDOT(calc_info_.noccB*calc_info_.nvirB,X_BS[0],1,
    calc_info_.sB[0],1);

  e11 += 8.0*t1*t2;

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(X_AR[0][0]),
    calc_info_.nvirA);

  double **B_p_BS = get_BS_ints(1);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(X_BS[0][0]),1);

  free_block(B_p_BS);

  t1 = C_DDOT(calc_info_.noccA*calc_info_.nvirA,X_AR[0],1,
    calc_info_.sA[0],1);
  t2 = C_DDOT(calc_info_.noccB*calc_info_.nvirB,X_BS[0],1,
    calc_info_.sB[0],1);
  
  e11 += 8.0*t1*t2;

  free_block(X_AR);
  free_block(X_BS);

  e20 = exch_ind30_20();
  e02 = exch_ind30_02();

  results_.exch_ind30 = e10+e01+e11+e20+e02;
  if (params_.print) {
    fprintf(outfile,"Exch-Ind30         = %18.12lf  H\n\n",
      results_.exch_ind30);
    fflush(outfile);
  }
}

double SAPT2p3::exch_ind30_20()
{
  double energy = 0.0;

  double **B_p_AR = get_AR_ints(1);
  double **B_p_RB = get_RB_ints(1);

  double **xRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(xRB[0][0]),calc_info_.noccB);

  double **B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(xRB[0][0]),calc_info_.noccB,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(B_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  double *X = init_array(calc_info_.nrio);
  double *Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
	  &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.sA[0],1,0.0,X,1);

  C_DGEMV('t',calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_RB[0][0]),calc_info_.nrio,xRB[0],1,0.0,Y,1);

  energy -= 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free_block(B_p_RB);

  double **xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);
  double **yAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAR[0][0]),
    calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,yAR[0],1);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);
  double **yAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.nvirA,1.0,
    &(xAR[0][0]),calc_info_.nvirA,&(calc_info_.sA[0][0]),calc_info_.nvirA,
    0.0,&(xAA[0][0]),calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,&(yAR[0][0]),calc_info_.nvirA,
    0.0,&(yAA[0][0]),calc_info_.noccA);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,yAA[0],1);

  free_block(xAR);
  free_block(yAR);
  free_block(xAA);
  free_block(yAA);

  double **B_p_BB = get_BB_ints(1);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.sA[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(xAB[0][0]),calc_info_.noccB,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(xAB[0][0]),calc_info_.noccB,
    0.0,&(xBB[0][0]),calc_info_.noccB);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,xBB[0],1,0.0,Y,1);

  energy += 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xAB);
  free_block(xBB);
  free_block(xRB);
  free_block(B_p_AB);
  free_block(C_p_AB);
  free_block(B_p_BB);
  free_block(B_p_AR);

  return(energy);
}

double SAPT2p3::exch_ind30_02()
{
  double energy = 0.0;

  double **B_p_BS = get_BS_ints(1);
  double **B_p_AS = get_AS_ints(1);

  double **xAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.sB[0][0]),
    calc_info_.nvirB,0.0,&(xAS[0][0]),calc_info_.nvirB);

  double **B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xAS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,
      calc_info_.nvirB,1.0,&(calc_info_.sB[0][0]),calc_info_.nvirB,
      &(B_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  double *X = init_array(calc_info_.nrio);
  double *Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
	  &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.sB[0],1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_AS[0][0]),calc_info_.nrio,xAS[0],1,0.0,Y,1);

  energy -= 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free_block(B_p_AS);

  double **xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);
  double **yBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(xBS[0][0]),
    calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,yBS[0],1);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);
  double **yBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(xBS[0][0]),calc_info_.nvirB,&(calc_info_.sB[0][0]),calc_info_.nvirB,
    0.0,&(xBB[0][0]),calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,&(yBS[0][0]),calc_info_.nvirB,
    0.0,&(yBB[0][0]),calc_info_.noccB);

  energy += 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,yBB[0],1);

  free_block(xBS);
  free_block(yBS);
  free_block(xBB);
  free_block(yBB);

  double **B_p_AA = get_AA_ints(1);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.sB[0][0]),calc_info_.nvirB,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {  
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(xAB[0][0]),calc_info_.noccB,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(xAB[0][0]),calc_info_.noccB,
    0.0,&(xAA[0][0]),calc_info_.noccA);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,xAA[0],1,0.0,Y,1);

  energy += 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xAB);
  free_block(xAA);
  free_block(xAS);
  free_block(B_p_AB);
  free_block(C_p_AB);
  free_block(B_p_AA);
  free_block(B_p_BS);

  return(energy);
}

}}

