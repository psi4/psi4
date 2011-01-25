/* This function calculates the Ind20,resp energy */

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
#include "sapt2p3.h"

namespace psi { namespace sapt {

void SAPT2p3::exch_disp30()
{
  double e11,e20,e02,e22;
  if (params_.print)
    fprintf(outfile,"Begining Exch-Disp30 Calculation\n\n");

  double **X_ARBS = read_IJKL(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);
  double **Y_ARBS = read_IJKL(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  e11 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nvirB*
    calc_info_.noccB,&(X_ARBS[0][0]),1,&(Y_ARBS[0][0]),1);

  free_block(X_ARBS);
  free_block(Y_ARBS);

  e20 = exch_disp30_20();
  e02 = exch_disp30_02();
  e22 = exch_disp30_22();
  results_.exch_disp30 = e11 + e20 + e02 + e22;

  if (params_.print) {
    fprintf(outfile,"Exch-Disp30        = %18.12lf  H\n\n",
      results_.exch_disp30);
    fflush(outfile);
  }
}

double SAPT2p3::exch_disp30_20()
{
  double energy = 0.0;

  double **uARAR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirA);
  double **B_p_AR = get_AR_ints(1);
  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nri,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(T_p_AR[0][0]),calc_info_.nrio,0.0,&(uARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA);

  free_block(T_p_AR);

  for(int ar=0; ar<calc_info_.noccA*calc_info_.nvirA; ar++) {
    for(int a1r1=0; a1r1<ar; a1r1++) {
      double tval = uARAR[ar][a1r1] + uARAR[a1r1][ar];
      uARAR[a1r1][ar] = tval;
      uARAR[ar][a1r1] = tval;
  }}

  C_DSCAL(calc_info_.noccA*calc_info_.nvirA,2.0,&(uARAR[0][0]),
    calc_info_.noccA*calc_info_.nvirA+1);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int aa=0, aarr=0; aa < calc_info_.noccA; aa++) {
    for (int rr=0; rr < calc_info_.nvirA; rr++, aarr++) {
      double denom = calc_info_.evalsA[a]+calc_info_.evalsA[aa]-
        calc_info_.evalsA[r+calc_info_.noccA]-
        calc_info_.evalsA[rr+calc_info_.noccA];
      uARAR[ar][aarr] /= denom;
    }}
  }}

  double **U_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,
    calc_info_.noccA*calc_info_.nvirA,1.0,&(uARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA,&(B_p_AR[0][0]),calc_info_.nrio,0.0,&(U_p_AR[0][0]),
    calc_info_.nrio);

  double *X = init_array(calc_info_.nvirA);

  for(int a=0; a<calc_info_.noccA; a++) {
  for(int a1=0; a1<=a; a1++) {
    for(int r=0; r<calc_info_.nvirA; r++) {
      int ar = a*calc_info_.nvirA+r;
      int a1r = a1*calc_info_.nvirA+r;
      C_DCOPY(calc_info_.nvirA,&(uARAR[ar][a1*calc_info_.nvirA]),1,X,1);
      C_DCOPY(calc_info_.nvirA,&(uARAR[a1r][a*calc_info_.nvirA]),1,
        &(uARAR[ar][a1*calc_info_.nvirA]),1);
      C_DCOPY(calc_info_.nvirA,X,1,&(uARAR[a1r][a*calc_info_.nvirA]),1);
  }}}

  free(X);

  double **U_p_ApR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,
    calc_info_.noccA*calc_info_.nvirA,1.0,&(uARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA,&(B_p_AR[0][0]),calc_info_.nrio,0.0,&(U_p_ApR[0][0]),
    calc_info_.nrio);

  free_block(B_p_AR);
  free_block(uARAR);

  double **B_p_RB = get_RB_ints(1);
  double **X_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  for(int r=0; r<calc_info_.nvirA; r++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
      &(B_p_RB[r*calc_info_.noccB][0]),calc_info_.nrio,0.0,
      &(X_p_AR[r][0]),calc_info_.nvirA*calc_info_.nrio);
  }

  energy = C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    &(U_p_ApR[0][0]),1,&(X_p_AR[0][0]),1);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    &(U_p_AR[0][0]),1,&(X_p_AR[0][0]),1);

  free_block(B_p_RB);
  free_block(X_p_AR);

  double **xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);
  double **yAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAR[0][0]),
    calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(U_p_ApR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,yAR[0],1);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,&(xAR[0][0]),1,
    &(yAR[0][0]),1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(U_p_AR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,yAR[0],1);

  energy -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,&(xAR[0][0]),1,
    &(yAR[0][0]),1);

  free_block(xAR);
  free_block(yAR);

  double **A_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **A_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,
      1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(U_p_ApR[a*calc_info_.nvirA][0]),calc_info_.nrio,0.0,
      &(A_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,-1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(A_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,
    &(A_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,
      1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(U_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,0.0,
      &(A_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,2.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(A_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,1.0,
    &(A_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio);

  double **B_p_BB = get_BB_ints(1);

  energy += C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    &(A_p_BB[0][0]),1,&(B_p_BB[0][0]),1);

  free_block(A_p_AB);
  free_block(A_p_BB);
  free_block(U_p_AR);
  free_block(U_p_ApR);
  free_block(B_p_BB);

  return(4.0*energy);
}

double SAPT2p3::exch_disp30_02()
{
  double energy = 0.0;

  double **uBSBS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.nvirB);
  double **B_p_BS = get_BS_ints(1);
  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nri,1.0,&(B_p_BS[0][0]),calc_info_.nrio,
    &(T_p_BS[0][0]),calc_info_.nrio,0.0,&(uBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(T_p_BS);

  for(int bs=0; bs<calc_info_.noccB*calc_info_.nvirB; bs++) {
    for(int b1s1=0; b1s1<bs; b1s1++) {
      double tval = uBSBS[bs][b1s1] + uBSBS[b1s1][bs];
      uBSBS[b1s1][bs] = tval;
      uBSBS[bs][b1s1] = tval;
  }}

  C_DSCAL(calc_info_.noccB*calc_info_.nvirB,2.0,&(uBSBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB+1);

  for (int b=0, bs=0; b < calc_info_.noccB; b++) {
  for (int s=0; s < calc_info_.nvirB; s++, bs++) {
    for (int bb=0, bbss=0; bb < calc_info_.noccB; bb++) {
    for (int ss=0; ss < calc_info_.nvirB; ss++, bbss++) {
      double denom = calc_info_.evalsB[b]+calc_info_.evalsB[bb]-
        calc_info_.evalsB[s+calc_info_.noccB]-
        calc_info_.evalsB[ss+calc_info_.noccB];
      uBSBS[bs][bbss] /= denom;
    }}
  }}

  double **U_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,
    calc_info_.noccB*calc_info_.nvirB,1.0,&(uBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,&(B_p_BS[0][0]),calc_info_.nrio,0.0,&(U_p_BS[0][0]),
    calc_info_.nrio);

  double *X = init_array(calc_info_.nvirB);

  for(int b=0; b<calc_info_.noccB; b++) {
  for(int b1=0; b1<=b; b1++) {
    for(int s=0; s<calc_info_.nvirB; s++) {
      int bs = b*calc_info_.nvirB+s;
      int b1s = b1*calc_info_.nvirB+s;
      C_DCOPY(calc_info_.nvirB,&(uBSBS[bs][b1*calc_info_.nvirB]),1,X,1);
      C_DCOPY(calc_info_.nvirB,&(uBSBS[b1s][b*calc_info_.nvirB]),1,
        &(uBSBS[bs][b1*calc_info_.nvirB]),1);
      C_DCOPY(calc_info_.nvirB,X,1,&(uBSBS[b1s][b*calc_info_.nvirB]),1);
  }}}

  free(X);

  double **U_p_BpS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,
    calc_info_.noccB*calc_info_.nvirB,1.0,&(uBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,&(B_p_BS[0][0]),calc_info_.nrio,0.0,&(U_p_BpS[0][0]),
    calc_info_.nrio);

  free_block(B_p_BS);
  free_block(uBSBS);

  double **B_p_AS = get_AS_ints(1);
  double **X_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AS[0][0]),calc_info_.nvirB*calc_info_.nrio,0.0,
    &(X_p_BS[0][0]),calc_info_.nvirB*calc_info_.nrio);

  energy = C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    &(U_p_BpS[0][0]),1,&(X_p_BS[0][0]),1);

  energy -= 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    &(U_p_BS[0][0]),1,&(X_p_BS[0][0]),1);

  free_block(B_p_AS);
  free_block(X_p_BS);

  double **xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);
  double **yBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(xBS[0][0]),
    calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(U_p_BpS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,yBS[0],1);

  energy += 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,&(xBS[0][0]),1,
    &(yBS[0][0]),1);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(U_p_BS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,yBS[0],1);

  energy -= 4.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,&(xBS[0][0]),1,
    &(yBS[0][0]),1);

  free_block(xBS);
  free_block(yBS);

  double **A_p_BA = block_matrix(calc_info_.noccB*calc_info_.noccA,
    calc_info_.nrio);
  double **A_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,
      1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(U_p_BpS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(A_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccB,-1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(A_p_BA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,
    &(A_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,
      1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(U_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(A_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccB,2.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(A_p_BA[0][0]),calc_info_.noccA*calc_info_.nrio,1.0,
    &(A_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio);

  double **B_p_AA = get_AA_ints(1);

  energy += C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    &(A_p_AA[0][0]),1,&(B_p_AA[0][0]),1);

  free_block(A_p_BA);
  free_block(A_p_AA);
  free_block(U_p_BS);
  free_block(U_p_BpS);
  free_block(B_p_AA);

  return(4.0*energy);
}

double SAPT2p3::exch_disp30_22()
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

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  double **B_p_AR = get_AR_ints(0);
  double **B_p_BS = get_BS_ints(0);

  double **X_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);
  double **xSS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(tAS_RB[0][0]),calc_info_.noccB,0.0,&(xBB[0][0]),calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccB,2.0,&(xBB[0][0]),calc_info_.noccB,&(B_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio,0.0,&(X_p_BS[0][0]),calc_info_.nvirB*
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,&(tRB_AS[0][0]),
    calc_info_.nvirB,0.0,&(xSS[0][0]),calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.nvirB,2.0,
      &(xSS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,1.0,&(X_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  energy += C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    &(T_p_BS[0][0]),1,&(X_p_BS[0][0]),1);

  free_block(xBB);
  free_block(xSS);
  free_block(X_p_BS);

  double **X_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);
  double **xRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(tRB_AS[0][0]),calc_info_.nvirB,0.0,&(xAA[0][0]),calc_info_.noccA);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.noccA,2.0,&(xAA[0][0]),calc_info_.noccA,&(B_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio,0.0,&(X_p_AR[0][0]),calc_info_.nvirA*
    calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,&(tAS_RB[0][0]),
    calc_info_.noccB,0.0,&(xRR[0][0]),calc_info_.nvirA);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.nvirA,2.0,
      &(xRR[0][0]),calc_info_.nvirA,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,1.0,&(X_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio);
  }

  energy += C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    &(T_p_AR[0][0]),1,&(X_p_AR[0][0]),1);

  free_block(xAA);
  free_block(xRR);
  free_block(X_p_AR);

  double **A_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **B_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,
      1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(T_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,0.0,
      &(A_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,
      1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(T_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }

  energy -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    &(A_p_AB[0][0]),1,&(B_p_AB[0][0]),1);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,
      1.0,&(tAS_RB[0][0]),calc_info_.noccB,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(A_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,
      1.0,&(tRB_AS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio);
  }

  energy -= C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    &(A_p_AB[0][0]),1,&(B_p_AB[0][0]),1);

  free_block(A_p_AB);
  free_block(B_p_AB);
  free_block(T_p_AR);
  free_block(T_p_BS);

  double **tABRS = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nvirA*calc_info_.nvirB);

  for (int a=0,ar=0; a<calc_info_.noccA; a++) {
    for (int r=0; r<calc_info_.nvirA; r++,ar++) {
      for (int b=0; b<calc_info_.noccB; b++) {
        int ab = a*calc_info_.noccB+b;
        C_DCOPY(calc_info_.nvirB,&(tARBS[ar][b*calc_info_.nvirB]),1,
          &(tABRS[ab][r*calc_info_.nvirB]),1);
  }}}

  free_block(tARBS);

  double **vABRS = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nvirA*calc_info_.nvirB);

  for(int a=0,ab=0; a<calc_info_.noccA; a++) {
    for(int b=0; b<calc_info_.noccB; b++,ab++) {
      C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nvirB,calc_info_.nrio,
        1.0,&(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,
        &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
        &(vABRS[ab][0]),calc_info_.nvirB);
  }}

  free_block(B_p_AR);
  free_block(B_p_BS);

  double **xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);
  double **ABAB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.noccA*calc_info_.noccB);

  for(int a=0,ab=0; a<calc_info_.noccA; a++) {
    for(int b=0; b<calc_info_.noccB; b++,ab++) {
      C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.nvirB,1.0,
        &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
        &(tABRS[ab][0]),calc_info_.nvirB,0.0,&(xAR[0][0]),calc_info_.nvirA);
      C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
        &(xAR[0][0]),calc_info_.nvirA,&(calc_info_.S_AB[calc_info_.noccA][0]),
        calc_info_.nmo,0.0,&(ABAB[ab][0]),calc_info_.noccB);
  }}

  free_block(xAR);

  double **xABRS = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nvirA*calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccA*calc_info_.noccB,calc_info_.nvirA*
    calc_info_.nvirB,calc_info_.noccA*calc_info_.noccB,1.0,&(ABAB[0][0]),
    calc_info_.noccA*calc_info_.noccB,&(vABRS[0][0]),calc_info_.nvirA*
    calc_info_.nvirB,0.0,&(xABRS[0][0]),calc_info_.nvirA*calc_info_.nvirB);

  energy -= C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nvirA*
    calc_info_.nvirB,&(tABRS[0][0]),1,&(xABRS[0][0]),1);

  free_block(tABRS);
  free_block(ABAB);
  free_block(vABRS);
  free_block(xABRS);

  return(2.0*energy);
}

}}
