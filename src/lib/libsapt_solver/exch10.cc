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
#include "sapt0.h"

namespace psi { namespace sapt {

void SAPT0::exch10_s2()
{
  double ex1, ex2, ex3, ex4, ex5, ex6;

  if (params_.print)
    fprintf(outfile,"Begining Exch10 Calculation\n\n");

  double **B_p_AB = get_AB_ints(1);
  double **B_q_AB = get_AB_ints(2);

  ex1 = C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    &(B_p_AB[0][0]),1,&(B_q_AB[0][0]),1);

  double **B_p_AA = get_AA_ints(1);
  double **B_p_BB = get_BB_ints(1);

  double **X_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++)
    C_DCOPY(calc_info_.noccB,&(calc_info_.S_AB[a][0]),1,&(X_AB[a][0]),1);

  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++){
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,1.0,
      &(X_AB[0][0]),calc_info_.noccB,&(B_q_AB[a*calc_info_.noccB][0]),
      calc_info_.nrio,0.0,&(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  double *Ap_diag = init_array(calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++){
    int aa = a*calc_info_.noccA+a;
    C_DAXPY(calc_info_.nrio,1.0,&(C_p_AA[aa][0]),1,&(Ap_diag[0]),1);
  }

  ex2 = 2.0*C_DDOT(calc_info_.nrio,calc_info_.diagAA,1,Ap_diag,1); 
  ex2 -= C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    &(B_p_AA[0][0]),1,&(C_p_AA[0][0]),1);

  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(X_AB[0][0]),calc_info_.noccB,&(B_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  double *Bp_diag = init_array(calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++){
    int bb = b*calc_info_.noccB+b;
    C_DAXPY(calc_info_.nrio,1.0,&(C_p_BB[bb][0]),1,&(Bp_diag[0]),1);
  }

  ex3 = 2.0*C_DDOT(calc_info_.nrio,calc_info_.diagBB,1,Bp_diag,1);
  ex3 -= C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    &(B_p_BB[0][0]),1,&(C_p_BB[0][0]),1);

  free_block(C_p_AA);
  free_block(C_p_BB);

  double **X_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
          &(X_AB[0][0]),calc_info_.noccB,&(X_AB[0][0]),calc_info_.noccB,
          0.0,&(X_AA[0][0]),calc_info_.noccA);
  
  double **X_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);
  
  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
          &(X_AB[0][0]),calc_info_.noccB,&(X_AB[0][0]),calc_info_.noccB,
          0.0,&(X_BB[0][0]),calc_info_.noccB);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,&(X_BB[0][0]),1,0.0,Bp_diag,1);

  ex4 = 2.0*C_DDOT(calc_info_.nrio,calc_info_.diagAA,1,Bp_diag,1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,&(X_AA[0][0]),1,0.0,Ap_diag,1);

  ex5 = 2.0*C_DDOT(calc_info_.nrio,calc_info_.diagBB,1,Ap_diag,1);

  free(Ap_diag);
  free(Bp_diag);
  free_block(X_AA);
  free_block(X_BB);

  for(int a=0; a<calc_info_.noccA; a++){
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(X_AB[0][0]),calc_info_.noccB,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(X_AB[0][0]),calc_info_.noccB,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(B_q_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  ex6 = C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    &(B_p_AB[0][0]),1,&(B_q_AB[0][0]),1);

  free_block(X_AB);
  free_block(B_p_AA);
  free_block(B_p_BB);
  free_block(B_p_AB);
  free_block(B_q_AB);

  results_.exch10_s2 = -2.0*(ex1+ex2+ex3-ex4-ex5+ex6);

  if (params_.print) {
    fprintf(outfile,"exch10(s^2) Energy = %18.12lf  H\n",results_.exch10_s2);
    fflush(outfile);
  }
}

void SAPT0::exch10()
{
  double ex1=0, ex2=0, ex3=0, ex4=0, ex5=0, ex6=0, ex7=0, ex8=0, ex9=0;

  double **P = block_matrix(calc_info_.noccA+calc_info_.noccB,
    calc_info_.noccA+calc_info_.noccB);
  double **Q = block_matrix(calc_info_.noccA+calc_info_.noccB,
    calc_info_.noccA+calc_info_.noccB);
  
  for (int i=0; i<calc_info_.noccA+calc_info_.noccB; i++)
    Q[i][i] = 1.0;
  
  for (int a=0; a<calc_info_.noccA; a++) { 
    for (int b=0; b<calc_info_.noccB; b++) {
      Q[a][b+calc_info_.noccA] = calc_info_.S_AB[a][b];
      Q[b+calc_info_.noccA][a] = calc_info_.S_AB[a][b];
  }}
  
  invert_matrix(Q,P,calc_info_.noccA+calc_info_.noccB,outfile);
  
  free_block(Q);
  
  for (int i=0; i<calc_info_.noccA+calc_info_.noccB; i++)
    P[i][i] -= 1.0;
  
  double **pAA = block_matrix(calc_info_.noccA,calc_info_.noccA);
  double **pBB = block_matrix(calc_info_.noccB,calc_info_.noccB);
  double **pAB = block_matrix(calc_info_.noccA,calc_info_.noccB);
  double **pBA = block_matrix(calc_info_.noccB,calc_info_.noccA);
  
  for (int a1=0; a1<calc_info_.noccA; a1++) { 
    for (int a2=0; a2<calc_info_.noccA; a2++) {
      pAA[a1][a2] = P[a1][a2];
  }}
  
  for (int b1=0; b1<calc_info_.noccB; b1++) { 
    for (int b2=0; b2<calc_info_.noccB; b2++) {
      pBB[b1][b2] = P[b1+calc_info_.noccA][b2+calc_info_.noccA];
  }}
  
  for (int a=0; a<calc_info_.noccA; a++) { 
    for (int b=0; b<calc_info_.noccB; b++) {
      pAB[a][b] = P[a][b+calc_info_.noccA];
  }}
  
  for (int b=0; b<calc_info_.noccB; b++) { 
    for (int a=0; a<calc_info_.noccA; a++) {
      pBA[b][a] = P[b+calc_info_.noccA][a];
  }}

  free_block(P);

  double **B_p_AB = get_AB_ints(1);
  double **A_p_AB = get_AB_ints(2);
  double **B_p_AA = get_AA_ints(1);
  double **A_p_BB = get_BB_ints(1);

  ex1 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    &(B_p_AB[0][0]),1,&(A_p_AB[0][0]),1);

  double *X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,&(pAA[0][0]),1,0.0,X,1);

  ex2 = 4.0*C_DDOT(calc_info_.nrio,calc_info_.diagBB,1,X,1);

  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,pAA[0],calc_info_.noccA,B_p_AB[0],calc_info_.noccB*
    calc_info_.nrio,0.0,C_p_AB[0],calc_info_.noccB*calc_info_.nrio);

  ex2 -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    A_p_AB[0],1,C_p_AB[0],1);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(A_p_BB[0][0]),calc_info_.nrio,&(pBB[0][0]),1,0.0,X,1);

  ex3 = 4.0*C_DDOT(calc_info_.nrio,calc_info_.diagAA,1,X,1);

  for (int a1=0; a1<calc_info_.noccA; a1++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccB,1.0,
      pBB[0],calc_info_.noccB,A_p_AB[a1*calc_info_.noccB],calc_info_.nrio,0.0,
      C_p_AB[a1*calc_info_.noccB],calc_info_.nrio);
  }

  ex3 -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(A_p_AB[0][0]),calc_info_.nrio,&(pAB[0][0]),1,0.0,X,1);

  ex4 = 4.0*C_DDOT(calc_info_.nrio,calc_info_.diagAA,1,X,1);

  for (int a1=0; a1<calc_info_.noccA; a1++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      pAB[0],calc_info_.noccB,B_p_AA[a1*calc_info_.noccA],calc_info_.nrio,0.0,
      C_p_AB[a1*calc_info_.noccB],calc_info_.nrio);
  }

  ex4 -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    A_p_AB[0],1,C_p_AB[0],1);

  free_block(C_p_AB);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,&(pAB[0][0]),1,0.0,X,1);

  ex5 = 4.0*C_DDOT(calc_info_.nrio,calc_info_.diagBB,1,X,1);

  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,pAB[0],calc_info_.noccB,B_p_AB[0],calc_info_.noccB*
    calc_info_.nrio,0.0,C_p_BB[0],calc_info_.noccB*calc_info_.nrio);

  ex5 -= 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    A_p_BB[0],1,C_p_BB[0],1);

  free_block(C_p_BB);

  double *Y = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,&(pAB[0][0]),1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(A_p_BB[0][0]),calc_info_.nrio,&(pBB[0][0]),1,0.0,Y,1);

  ex6 = 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  double **D_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **E_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,pAB[0],calc_info_.noccB,A_p_BB[0],calc_info_.noccB*
    calc_info_.nrio,0.0,D_p_AB[0],calc_info_.noccB*calc_info_.nrio);

  for (int a1=0; a1<calc_info_.noccA; a1++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccB,1.0,
      pBB[0],calc_info_.noccB,D_p_AB[a1*calc_info_.noccB],calc_info_.nrio,0.0,
      E_p_AB[a1*calc_info_.noccB],calc_info_.nrio);
  }

  ex6 -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,E_p_AB[0],1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(A_p_AB[0][0]),calc_info_.nrio,&(pAB[0][0]),1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,&(pAA[0][0]),1,0.0,Y,1);

  ex7 = 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  for (int a1=0; a1<calc_info_.noccA; a1++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      pAB[0],calc_info_.noccB,B_p_AA[a1*calc_info_.noccA],calc_info_.nrio,0.0,
      D_p_AB[a1*calc_info_.noccB],calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,pAA[0],calc_info_.noccA,D_p_AB[0],calc_info_.noccB*
    calc_info_.nrio,0.0,E_p_AB[0],calc_info_.noccB*calc_info_.nrio);

  ex7 -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    A_p_AB[0],1,E_p_AB[0],1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,&(pAA[0][0]),1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(A_p_BB[0][0]),calc_info_.nrio,&(pBB[0][0]),1,0.0,Y,1);

  ex8 = 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,pAA[0],calc_info_.noccA,B_p_AB[0],calc_info_.noccB*
    calc_info_.nrio,0.0,D_p_AB[0],calc_info_.noccB*calc_info_.nrio);

  for (int a1=0; a1<calc_info_.noccA; a1++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccB,1.0,
      pBB[0],calc_info_.noccB,D_p_AB[a1*calc_info_.noccB],calc_info_.nrio,0.0,
      E_p_AB[a1*calc_info_.noccB],calc_info_.nrio);
  }

  ex8 -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    A_p_AB[0],1,E_p_AB[0],1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(A_p_AB[0][0]),calc_info_.nrio,&(pAB[0][0]),1,0.0,X,1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,&(pAB[0][0]),1,0.0,Y,1);

  ex9 = 4.0*C_DDOT(calc_info_.nrio,X,1,Y,1);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,pAB[0],calc_info_.noccB,A_p_BB[0],calc_info_.noccB*
    calc_info_.nrio,0.0,D_p_AB[0],calc_info_.noccB*calc_info_.nrio);

  for (int a1=0; a1<calc_info_.noccA; a1++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      pAB[0],calc_info_.noccB,B_p_AA[a1*calc_info_.noccA],calc_info_.nrio,0.0,
      E_p_AB[a1*calc_info_.noccB],calc_info_.nrio);
  }

  ex9 -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    D_p_AB[0],1,E_p_AB[0],1);

  free(X);
  free(Y);
  free_block(D_p_AB);
  free_block(E_p_AB);

  free_block(pAA);
  free_block(pBB);
  free_block(pAB);
  free_block(pBA);

  free_block(B_p_AA);
  free_block(A_p_BB);
  free_block(A_p_AB);
  free_block(B_p_AB);

  results_.exch10 = ex1+ex2+ex3+ex4+ex5+ex6+ex7+ex8+ex9;

  if (params_.print) {
    fprintf(outfile,"exch10      Energy = %18.12lf  H\n\n",results_.exch10);
    fflush(outfile);
  }
}

}}
