/* This function calculates the Disp20 energy */

#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <cstring>
#include <iostream>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libipv1/ip_data.gbl>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::exch12()
{
  double ex111, ex120_k2u, ex102_k2u, ex120_k2f, ex102_k2f, ex120_k11u, 
    ex102_k11u;

  if (params_.print) {
    fprintf(outfile,"Begining Exch12 Calculation\n\n");
    fflush(outfile);
  }

  ex111 = exch111();

  if (params_.print) {
    fprintf(outfile,"exch111            = %18.12lf  H\n",ex111);
    fflush(outfile);
  }

  ex120_k2u = exch110("Theta(2)(AR) AR");

  if (params_.print) {
    fprintf(outfile,"exch120_k2u        = %18.12lf  H\n",ex120_k2u);
    fflush(outfile);
  }

  ex102_k2u = exch101("Theta(2)(BS) BS");

  if (params_.print) {
    fprintf(outfile,"exch102_k2u        = %18.12lf  H\n", ex102_k2u);
    fflush(outfile);
  }

  ex120_k2f = exch120_k2f();
  
  if (params_.print) {
    fprintf(outfile,"exch120_k2f        = %18.12lf  H\n",ex120_k2f);
    fflush(outfile);
  }

  ex102_k2f = exch102_k2f();
  
  if (params_.print) {
    fprintf(outfile,"exch102_k2f        = %18.12lf  H\n",ex102_k2f);
    fflush(outfile);
  }

  timer_on("ex120_k11u     ");
    ex120_k11u = exch120_k11u_1();
    ex120_k11u += exch120_k11u_2();
    ex120_k11u += exch120_k11u_3();
    ex120_k11u += exch120_k11u_4();
    ex120_k11u += exch120_k11u_5();
    ex120_k11u += exch120_k11u_6();
  timer_off("ex120_k11u     ");

  if (params_.print) {
    fprintf(outfile,"exch120_k11u       = %18.12lf  H\n",ex120_k11u);
    fflush(outfile);
  }

  timer_on("ex102_k11u     ");
    ex102_k11u = exch102_k11u_1();
    ex102_k11u += exch102_k11u_2();
    ex102_k11u += exch102_k11u_3();
    ex102_k11u += exch102_k11u_4();
    ex102_k11u += exch102_k11u_5();
    ex102_k11u += exch102_k11u_6();
  timer_off("ex102_k11u     ");

  if (params_.print) {
    fprintf(outfile,"exch102_k11u       = %18.12lf  H\n\n",ex102_k11u);
    fflush(outfile);
  }

  results_.exch12 = ex111 + ex120_k2u + ex102_k2u + ex120_k2f + ex102_k2f + 
    ex120_k11u + ex102_k11u;
}

double SAPT2::exch111()
{
  double e1, e2;

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"Theta(AR) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"Theta(BS) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  double **C_p_BA = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(T_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,0.0,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(T_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(C_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio);
  }

  e1 = 0.0;

  for (int a=0, ab=0; a<calc_info_.noccA; a++) {
    for (int b=0; b<calc_info_.noccB; b++, ab++) {
      int ba = b*calc_info_.noccA+a;
      e1 -= 4.0*C_DDOT(calc_info_.nrio,C_p_AB[ab],1,C_p_BA[ba],1);
  }}

  free_block(C_p_AB);
  free_block(C_p_BA);

  double **C_p_BR = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.S_AB[calc_info_.noccA][calc_info_.noccB]),calc_info_.nmo,
      &(T_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(C_p_BR[b*calc_info_.nvirA][0]),calc_info_.nrio);
  }

  free_block(T_p_BS);

  double **C_p_AR = block_matrix(calc_info_.nvirA*calc_info_.noccA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_BR[0][0]),calc_info_.nvirA*calc_info_.nrio,0.0,
    &(C_p_AR[0][0]),calc_info_.nvirA*calc_info_.nrio);

  e2 = -4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    &(T_p_AR[0][0]),1,&(C_p_AR[0][0]),1);

  free_block(T_p_AR);
  free_block(C_p_BR);
  free_block(C_p_AR);
 
  return(e1+e2);
}

double SAPT2::exch120_k2f()
{
  double e1,e2;

  double **T_AR = read_IJKL(PSIF_SAPT_AMPS,"T AR Amplitudes",calc_info_.noccA,
    calc_info_.nvirA);
  double **v_AR = read_IJKL(PSIF_SAPT_AMPS,"AR Exch12 K2f Integrals",
    calc_info_.noccA,calc_info_.nvirA);

  e1 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,T_AR[0],1,v_AR[0],1);

  free_block(v_AR);

  double **B_p_AB = get_AB_ints(2);
  double **B_p_RB = get_RB_ints(2);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(T_AR[0][0]),calc_info_.nvirA,&(B_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  e2 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  free_block(B_p_RB);

  double **X_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(C_p_AB[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(X_AB[0][0]),1);

  for (int a=0; a<calc_info_.noccA; a++) {
    e2 -= 4.0*C_DDOT(calc_info_.noccB,X_AB[a],1,calc_info_.S_AB[a],1);
  }

  double **B_p_BB = get_BB_ints(1);
  double **D_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(D_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  e2 += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    C_p_AB[0],1,D_p_AB[0],1);

  free_block(B_p_BB);
  free_block(C_p_AB);
  free_block(D_p_AB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(T_AR[0][0]),calc_info_.nvirA,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(X_AB[0][0]),calc_info_.noccB);

  double **Y_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(Y_AB[0][0]),1);

  e2 -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,X_AB[0],1,Y_AB[0],1);

  free_block(Y_AB);

  double **B_p_AA = get_AA_ints(1);
  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(X_AB[0][0]),calc_info_.noccB,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }  

  e2 += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  free_block(X_AB);
  free_block(B_p_AA);
  free_block(C_p_AB);

  double **B_p_AR = get_AR_ints(1);
  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(T_AR[0][0]),calc_info_.nvirA,&(B_p_AR[a*calc_info_.nvirA][0]),
      calc_info_.nrio,0.0,&(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  double **C_p_BA = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,&(C_p_BA[0][0]),
    calc_info_.noccA*calc_info_.nrio);

  free_block(C_p_AA);
  free_block(T_AR);

  for (int a=0,ab=0; a<calc_info_.noccA; a++) {
    for (int b=0; b<calc_info_.noccB; b++,ab++) {
      int ba = b*calc_info_.noccA+a;
      e2 += 2.0*C_DDOT(calc_info_.nrio,B_p_AB[ab],1,C_p_BA[ba],1);
  }}

  free_block(B_p_AB);
  free_block(C_p_BA);

  return(e1+e2);
}

double SAPT2::exch102_k2f()
{
  double e1,e2;

  double **T_BS = read_IJKL(PSIF_SAPT_AMPS,"T BS Amplitudes",calc_info_.noccB,
    calc_info_.nvirB);
  double **v_BS = read_IJKL(PSIF_SAPT_AMPS,"BS Exch12 K2f Integrals",
    calc_info_.noccB,calc_info_.nvirB);

  e1 = -2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,T_BS[0],1,v_BS[0],1);

  free_block(v_BS);

  double **B_p_AS = get_AS_ints(2);
  double **B_p_AB = get_AB_ints(1);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(T_BS[0][0]),calc_info_.nvirB,&(B_p_AS[a*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  free_block(B_p_AS);

  e2 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,B_p_AB[0],
    1,C_p_AB[0],1);

  double **X_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(C_p_AB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(X_AB[0][0]),1);

  for (int a=0; a<calc_info_.noccA; a++) {
    e2 -= 4.0*C_DDOT(calc_info_.noccB,calc_info_.S_AB[a],1,X_AB[a],1);
  }

  double **B_p_AA = get_AA_ints(1);
  double **D_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(D_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  e2 += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,C_p_AB[0],
    1,D_p_AB[0],1);

  free_block(B_p_AA);
  free_block(C_p_AB);
  free_block(D_p_AB);

  double **B_p_BS = get_BS_ints(1);
  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(T_BS[0][0]),calc_info_.nvirB,&(B_p_BS[b*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(C_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio);
  }

  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,
    &(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio);

  e2 += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,B_p_AB[0],
    1,C_p_AB[0],1);

  free_block(B_p_BS);
  free_block(C_p_BB);
  free_block(C_p_AB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,
    1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,&(T_BS[0][0]),
    calc_info_.nvirB,0.0,&(X_AB[0][0]),calc_info_.noccB);

  double **Y_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(Y_AB[0][0]),1);

  e2 -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,X_AB[0],1,Y_AB[0],1);

  free_block(Y_AB);

  double **B_p_BB = get_BB_ints(1);
  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(X_AB[0][0]),calc_info_.noccB,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  e2 += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,B_p_AB[0],
    1,C_p_AB[0],1);

  free_block(B_p_AB);
  free_block(B_p_BB);
  free_block(C_p_AB);

  free_block(T_BS);

  return(e1+e2);
}

double SAPT2::exch120_k11u_1()
{
  double energy=0.0;

  double **xRR = read_IJKL(PSIF_SAPT_AMPS,"RR MP2 OPDM",calc_info_.nvirA,
    calc_info_.nvirA);

  double **B_p_RB = get_RB_ints(1);
  double **C_p_RB = get_RB_ints(2);

  double **yRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nrio,1.0,&(B_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio,
    &(C_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(yRR[0][0]),
    calc_info_.nvirA);

  energy += 2.0*C_DDOT(calc_info_.nvirA*calc_info_.nvirA,xRR[0],1,yRR[0],1);

  free_block(yRR);

  double **D_p_RB = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(xRR[0][0]),calc_info_.nvirA,&(B_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(D_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(B_p_RB);

  double **E_p_RB = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(xRR[0][0]),calc_info_.nvirA,&(C_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(E_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(C_p_RB);

  double **B_p_AR = get_AR_ints(1);
  double **D_p_BR = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AR[0][0]),calc_info_.nvirA*calc_info_.nrio,0.0,&(D_p_BR[0][0]),
    calc_info_.nvirA*calc_info_.nrio);

  for (int b=0,br=0; b<calc_info_.noccB; b++) {
    for (int r=0; r<calc_info_.nvirA; r++,br++) {
      int rb = r*calc_info_.noccB+b;
      energy -= 2.0*C_DDOT(calc_info_.nrio,D_p_BR[br],1,D_p_RB[rb],1);
  }}

  double **xRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(D_p_RB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(xRB[0][0]),1);

  free_block(D_p_RB);

  for (int r=0; r<calc_info_.nvirA; r++) {
    energy += 4.0*C_DDOT(calc_info_.noccB,calc_info_.S_AB[r+calc_info_.noccA],
      1,xRB[r],1);
  }

  C_DGEMV('n',calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(E_p_RB[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(xRB[0][0]),1);

  for (int r=0; r<calc_info_.nvirA; r++) {
    energy += 4.0*C_DDOT(calc_info_.noccB,calc_info_.S_AB[r+calc_info_.noccA],
      1,xRB[r],1);
  }

  free_block(xRB);

  double **B_p_BB = get_BB_ints(1);
  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(E_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(E_p_RB);

  energy -= 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    B_p_BB[0],1,C_p_BB[0],1);

  free_block(C_p_BB);

  double **B_p_AB = get_AB_ints(2);
  double **yRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.noccB,calc_info_.nrio,1.0,
      B_p_AR[a*calc_info_.nvirA],calc_info_.nrio,B_p_AB[a*calc_info_.noccB],
      calc_info_.nrio,1.0,yRB[0],calc_info_.noccB);
  }

  free_block(B_p_AR);

  double **zRB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirA,1.0,
    xRR[0],calc_info_.nvirA,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,zRB[0],calc_info_.noccB);

  energy -= 2.0*C_DDOT(calc_info_.nvirA*calc_info_.noccB,yRB[0],1,zRB[0],1);

  free_block(yRB);

  double **xBR = block_matrix(calc_info_.noccB,calc_info_.nvirA);

  C_DGEMV('n',calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(D_p_BR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(xBR[0][0]),1);

  for (int b=0; b<calc_info_.noccB; b++) {
    for (int r=0; r<calc_info_.nvirA; r++) {
      energy -= 8.0*xBR[b][r]*zRB[r][b];
  }}  

  free_block(xBR);

  double **D_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      zRB[0],calc_info_.noccB,D_p_BR[b*calc_info_.nvirA],calc_info_.nrio,0.0,
      D_p_BB[b*calc_info_.noccB],calc_info_.nrio);
  }

  free_block(D_p_BR);

  energy += 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    &(B_p_BB[0][0]),1,&(D_p_BB[0][0]),1);

  free_block(D_p_BB);

  double **zBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.nvirA,1.0,
    calc_info_.S_AB[calc_info_.noccA],calc_info_.nmo,zRB[0],
    calc_info_.noccB,0.0,zBB[0],calc_info_.noccB);

  double **wBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(wBB[0][0]),1);

  energy -= 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,wBB[0],1,zBB[0],1);

  free_block(wBB);
  free_block(zBB);
  free_block(zRB);

  double **B_p_RR = get_RR_ints(1);

  double *X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.nvirA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_RR[0][0]),calc_info_.nrio,xRR[0],1,0.0,X,1);

  free_block(xRR);
  free_block(B_p_RR);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,X,1,0.0,&(xAB[0][0]),1);

  for (int a=0; a<calc_info_.noccA; a++) {
    energy += 4.0*C_DDOT(calc_info_.noccB,calc_info_.S_AB[a],1,xAB[a],1);
  }

  free_block(xAB);
  free_block(B_p_AB);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);
  double **yBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,X,1,0.0,&(xBB[0][0]),1);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    calc_info_.S_AB[0],calc_info_.nmo,calc_info_.S_AB[0],calc_info_.nmo,0.0,
    yBB[0],calc_info_.noccB);

  energy -= 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,yBB[0],1);

  free(X);
  free_block(xBB);
  free_block(yBB);
  free_block(B_p_BB);

  return(-energy);
}

double SAPT2::exch102_k11u_1()
{
  double energy=0.0;

  double **xSS = read_IJKL(PSIF_SAPT_AMPS,"SS MP2 OPDM",calc_info_.nvirB,
    calc_info_.nvirB);

  double **B_p_AS = get_AS_ints(1);
  double **C_p_AS = get_AS_ints(2);

  double **ySS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirB,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(B_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(ySS[0][0]),
      calc_info_.nvirB);
  }

  energy += 2.0*C_DDOT(calc_info_.nvirB*calc_info_.nvirB,xSS[0],1,ySS[0],1);

  free_block(ySS);

  double **D_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xSS[0][0]),calc_info_.nvirB,&(B_p_AS[a*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(D_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  free_block(B_p_AS);

  double **E_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(xSS[0][0]),calc_info_.nvirB,&(C_p_AS[a*calc_info_.nvirB][0]),
      calc_info_.nrio,0.0,&(E_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  free_block(C_p_AS);

  double **B_p_BS = get_BS_ints(1);
  double **F_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_BS[0][0]),calc_info_.nvirB*calc_info_.nrio,0.0,&(F_p_AS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirB*calc_info_.nrio,
    D_p_AS[0],1,F_p_AS[0],1);

  double **xAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(D_p_AS[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(xAS[0][0]),1);

  free_block(D_p_AS);

  for (int a=0; a<calc_info_.noccA; a++) {
    energy += 4.0*C_DDOT(calc_info_.nvirB,
      &(calc_info_.S_AB[a][calc_info_.noccB]),1,xAS[a],1);
  }

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(E_p_AS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(xAS[0][0]),1);

  for (int a=0; a<calc_info_.noccA; a++) {
    energy += 4.0*C_DDOT(calc_info_.nvirB,
      &(calc_info_.S_AB[a][calc_info_.noccB]),1,xAS[a],1);
  }

  double **B_p_AA = get_AA_ints(1);
  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(E_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  free_block(E_p_AS);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    B_p_AA[0],1,C_p_AA[0],1);

  free_block(C_p_AA);

  double **B_p_AB = get_AB_ints(1);
  double **yAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirB,calc_info_.nrio,1.0,
      B_p_AB[b],calc_info_.noccB*calc_info_.nrio,B_p_BS[b*calc_info_.nvirB],
      calc_info_.nrio,1.0,yAS[0],calc_info_.nvirB);
  }

  free_block(B_p_BS);

  double **zAS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,xSS[0],
    calc_info_.nvirB,0.0,zAS[0],calc_info_.nvirB);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirB,yAS[0],1,zAS[0],1);

  free_block(yAS);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(F_p_AS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(xAS[0][0]),1);

  energy -= 8.0*C_DDOT(calc_info_.noccA*calc_info_.nvirB,xAS[0],1,zAS[0],1);

  free_block(xAS);

  double **D_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      zAS[0],calc_info_.nvirB,F_p_AS[a*calc_info_.nvirB],calc_info_.nrio,0.0,
      D_p_AA[a*calc_info_.noccA],calc_info_.nrio);
  }

  free_block(F_p_AS);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    &(B_p_AA[0][0]),1,&(D_p_AA[0][0]),1);

  free_block(D_p_AA);

  double **zAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.nvirB,1.0,
    zAS[0],calc_info_.nvirB,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,zAA[0],calc_info_.noccA);

  double **wAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(wAA[0][0]),1);

  energy -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,wAA[0],1,zAA[0],1);

  free_block(wAA);
  free_block(zAA);
  free_block(zAS);

  double **B_p_SS = get_SS_ints(1);

  double *X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.nvirB*calc_info_.nvirB,calc_info_.nrio,1.0,
    B_p_SS[0],calc_info_.nrio,xSS[0],1,0.0,X,1);

  free_block(xSS);
  free_block(B_p_SS);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,X,1,0.0,&(xAB[0][0]),1);

  for (int a=0; a<calc_info_.noccA; a++) {
    energy += 4.0*C_DDOT(calc_info_.noccB,calc_info_.S_AB[a],1,xAB[a],1);
  }

  free_block(xAB);
  free_block(B_p_AB);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);
  double **yAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,X,1,0.0,&(xAA[0][0]),1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    calc_info_.S_AB[0],calc_info_.nmo,calc_info_.S_AB[0],calc_info_.nmo,0.0,
    yAA[0],calc_info_.noccA);

  energy -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,yAA[0],1);

  free(X);
  free_block(xAA);
  free_block(yAA);
  free_block(B_p_AA);

  return(-energy);
}

double SAPT2::exch120_k11u_2()
{
  double energy=0.0;

  double **xAA = read_IJKL(PSIF_SAPT_AMPS,"AA MP2 OPDM",calc_info_.noccA,
    calc_info_.noccA);

  double **B_p_AB = get_AB_ints(1);
  double **C_p_AB = get_AB_ints(2);

  double **yAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB*
    calc_info_.nrio,1.0,&(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,
    &(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(yAA[0][0]),
    calc_info_.noccA);

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,yAA[0],1);

  double **sAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DCOPY(calc_info_.noccB,calc_info_.S_AB[a],1,sAB[a],1);
  }

  double **B_p_AA = get_AA_ints(1);
  double *X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    C_p_AB[0],calc_info_.nrio,sAB[0],1,0.0,X,1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,X,1,0.0,&(yAA[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,yAA[0],1);

  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,1.0,
      &(sAB[0][0]),calc_info_.noccB,&(C_p_AB[a*calc_info_.noccB][0]),
      calc_info_.nrio,0.0,&(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccA*
    calc_info_.nrio,1.0,&(B_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio,
    &(C_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,&(yAA[0][0]),
    calc_info_.noccA);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,yAA[0],1);

  free_block(C_p_AA);

  double **tAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.noccA,1.0,
    &(xAA[0][0]),calc_info_.noccA,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(tAB[0][0]),calc_info_.noccB);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(C_p_AB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(xAB[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  memset(&(xAB[0][0]),'\0',sizeof(double)*calc_info_.noccA*calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(xAB[0][0]),
      calc_info_.noccB);
  }

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  free_block(C_p_AB);

  double **B_p_BB = get_BB_ints(1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(xAB[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.noccB*
    calc_info_.nrio,1.0,&(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,
    &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(xAB[0][0]),
    calc_info_.noccB);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  free_block(B_p_AB);

  double **sBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(sAB[0][0]),calc_info_.noccB,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(sBB[0][0]),calc_info_.noccB);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    B_p_BB[0],calc_info_.nrio,sBB[0],1,0.0,X,1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,X,1,0.0,&(yAA[0][0]),1);

  energy -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,yAA[0],1);

  free(X);
  free_block(sBB);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    B_p_BB[0],calc_info_.nrio,calc_info_.diagAA,1,0.0,xBB[0],1);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.noccB,1.0,
    &(sAB[0][0]),calc_info_.noccB,&(xBB[0][0]),calc_info_.noccB,0.0,
    &(xAB[0][0]),calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(yAA[0][0]),calc_info_.noccA);

  energy -= 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,yAA[0],1);

  free_block(xAA);
  free_block(xAB);
  free_block(xBB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(tAB[0][0]),calc_info_.noccB,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(yAA[0][0]),calc_info_.noccA);

  double **zAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(zAA[0][0]),1);

  energy -= 8.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,yAA[0],1,zAA[0],1);

  free_block(yAA);
  free_block(zAA);

  double **C_p_BA = block_matrix(calc_info_.noccB*calc_info_.noccA,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(tAB[0][0]),calc_info_.noccB,&(B_p_AA[0][0]),
    calc_info_.noccA*calc_info_.nrio,0.0,&(C_p_BA[0][0]),
    calc_info_.noccA*calc_info_.nrio);

  free_block(tAB);
  free_block(B_p_AA);

  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(sAB[0][0]),calc_info_.noccB,&(C_p_BA[b*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(C_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio);
  }

  free_block(sAB);
  free_block(C_p_BA);

  energy += 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    B_p_BB[0],1,C_p_BB[0],1);

  free_block(B_p_BB);
  free_block(C_p_BB);

  return(energy);
}

double SAPT2::exch102_k11u_2()
{
  double energy=0.0;

  double **xBB = read_IJKL(PSIF_SAPT_AMPS,"BB MP2 OPDM",calc_info_.noccB,
    calc_info_.noccB);

  double **B_p_AB = get_AB_ints(2);
  double **C_p_AB = get_AB_ints(1);

  double **yBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(yBB[0][0]),
      calc_info_.noccB);
  }

  energy += 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,yBB[0],1);

  double **sAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DCOPY(calc_info_.noccB,calc_info_.S_AB[a],1,sAB[a],1);
  }

  double **B_p_BB = get_BB_ints(1);
  double *X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    C_p_AB[0],calc_info_.nrio,sAB[0],1,0.0,X,1);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,X,1,0.0,&(yBB[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,yBB[0],1);

  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(sAB[0][0]),calc_info_.noccB,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  memset(&(yBB[0][0]),'\0',sizeof(double)*calc_info_.noccB*calc_info_.noccB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,
      &(C_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(yBB[0][0]),
      calc_info_.noccB);
  }

  energy -= 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,yBB[0],1);

  free_block(C_p_BB);

  double **tAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.noccB,1.0,
    &(sAB[0][0]),calc_info_.noccB,&(xBB[0][0]),calc_info_.noccB,0.0,
    &(tAB[0][0]),calc_info_.noccB);

  double **xAB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(C_p_AB[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(xAB[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,
    calc_info_.noccB*calc_info_.nrio,1.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(xAB[0][0]),calc_info_.noccB);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  free_block(C_p_AB);

  double **B_p_AA = get_AA_ints(1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(xAB[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  memset(&(xAB[0][0]),'\0',sizeof(double)*calc_info_.noccA*calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(xAB[0][0]),
      calc_info_.noccB);
  }

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,tAB[0],1,xAB[0],1);

  free_block(B_p_AB);

  double **sAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(sAB[0][0]),calc_info_.noccB,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(sAA[0][0]),calc_info_.noccA);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    B_p_AA[0],calc_info_.nrio,sAA[0],1,0.0,X,1);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,X,1,0.0,&(yBB[0][0]),1);

  energy -= 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,yBB[0],1);

  free(X);
  free_block(sAA);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    B_p_AA[0],calc_info_.nrio,calc_info_.diagBB,1,0.0,xAA[0],1);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.noccA,1.0,
    &(xAA[0][0]),calc_info_.noccA,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(xAB[0][0]),calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(xAB[0][0]),calc_info_.noccB,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(yBB[0][0]),calc_info_.noccB);

  energy -= 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,yBB[0],1);

  free_block(xBB);
  free_block(xAB);
  free_block(xAA);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(tAB[0][0]),calc_info_.noccB,&(sAB[0][0]),calc_info_.noccB,0.0,
    &(yBB[0][0]),calc_info_.noccB);

  double **zBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(zBB[0][0]),1);

  energy -= 8.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,yBB[0],1,zBB[0],1);

  free_block(yBB);
  free_block(zBB);

  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(tAB[0][0]),calc_info_.noccB,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(tAB);
  free_block(B_p_BB);

  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,1.0,
      &(sAB[0][0]),calc_info_.noccB,&(C_p_AB[a*calc_info_.noccB][0]),
      calc_info_.nrio,0.0,&(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  free_block(sAB);
  free_block(C_p_AB);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    B_p_AA[0],1,C_p_AA[0],1);

  free_block(B_p_AA);
  free_block(C_p_AA);

  return(energy);
}

double SAPT2::exch120_k11u_3()
{
  double energy=0.0;

  double **temp_thetaARAR = read_IJKL(PSIF_SAPT_AMPS,
    "T ARAR Antisym Amplitudes",calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirA);

  double **thetaRRAA = block_matrix(calc_info_.nvirA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.noccA);

  for(int a1=0,a1r1=0; a1<calc_info_.noccA; a1++) {
    for(int r1=0; r1<calc_info_.nvirA; r1++,a1r1++) {
      for(int a2=0,a2r2=0; a2<calc_info_.noccA; a2++) {
        for(int r2=0; r2<calc_info_.nvirA; r2++,a2r2++) {
          int a1a2 = a1*calc_info_.noccA+a2;
          int r1r2 = r1*calc_info_.nvirA+r2;
          thetaRRAA[r1r2][a1a2] = temp_thetaARAR[a1r1][a2r2];
  }}}}

  free_block(temp_thetaARAR);

  double **thetaRBAA = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.noccA*calc_info_.noccA);

  for(int r1=0; r1<calc_info_.nvirA; r1++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.noccA,
      calc_info_.nvirA,1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),
      calc_info_.nmo,&(thetaRRAA[r1*calc_info_.nvirA][0]),
      calc_info_.noccA*calc_info_.noccA,0.0,
      &(thetaRBAA[r1*calc_info_.noccB][0]),calc_info_.noccA*calc_info_.noccA);
  }

  free_block(thetaRRAA);

  double **temp_tARAR = read_IJKL(PSIF_SAPT_AMPS,"T ARAR Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

  double **tRRAA = block_matrix(calc_info_.nvirA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.noccA);

  for(int a1=0,a1r1=0; a1<calc_info_.noccA; a1++) {
    for(int r1=0; r1<calc_info_.nvirA; r1++,a1r1++) {
      for(int a2=0,a2r2=0; a2<calc_info_.noccA; a2++) {
        for(int r2=0; r2<calc_info_.nvirA; r2++,a2r2++) {
          int a1a2 = a1*calc_info_.noccA+a2;
          int r1r2 = r1*calc_info_.nvirA+r2;
          tRRAA[r1r2][a1a2] = temp_tARAR[a1r1][a2r2];
  }}}}

  free_block(temp_tARAR);

  double **B_RB_p = get_RB_ints(1);
  double **B_p_RB = block_matrix(calc_info_.nrio,calc_info_.nvirA*
    calc_info_.noccB);

  for(int r=0,rb=0; r<calc_info_.nvirA; r++) {
    for(int b=0; b<calc_info_.noccB; b++,rb++) {
      for(int P=0; P<calc_info_.nrio; P++) {
        B_p_RB[P][rb] = B_RB_p[rb][P];
  }}}

  free_block(B_RB_p);

  double **B_p_RR = get_RR_ints(1);

  double **yRB = block_matrix(calc_info_.nvirA,calc_info_.nvirA*
    calc_info_.noccB);
  double **C_p = block_matrix(calc_info_.nvirA,calc_info_.nrio);

  for (int r1=0; r1<calc_info_.nvirA; r1++) {
    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.nvirA,calc_info_.noccB,
      calc_info_.noccA*calc_info_.noccA,1.0,&(tRRAA[0][0]),
      calc_info_.noccA*calc_info_.noccA,&(thetaRBAA[r1*calc_info_.noccB][0]),
      calc_info_.noccA*calc_info_.noccA,0.0,&(yRB[0][0]),
      calc_info_.noccB);
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nrio,calc_info_.nvirA*
      calc_info_.noccB,1.0,&(yRB[0][0]),calc_info_.nvirA*calc_info_.noccB,
      &(B_p_RB[0][0]),calc_info_.nvirA*calc_info_.noccB,0.0,&(C_p[0][0]),
      calc_info_.nrio);
    energy += 2.0*C_DDOT(calc_info_.nvirA*calc_info_.nrio,
      B_p_RR[r1*calc_info_.nvirA],1,C_p[0],1);
  }

  free_block(yRB);
  free_block(C_p);
  free_block(B_p_RB);

  double **tRBAA = block_matrix(calc_info_.nvirA*calc_info_.noccB,
    calc_info_.noccA*calc_info_.noccA);

  for(int r1=0; r1<calc_info_.nvirA; r1++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.noccA,
      calc_info_.nvirA,1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),
      calc_info_.nmo,&(tRRAA[r1*calc_info_.nvirA][0]),
      calc_info_.noccA*calc_info_.noccA,0.0,&(tRBAA[r1*calc_info_.noccB][0]),
      calc_info_.noccA*calc_info_.noccA);
  }

  free_block(tRRAA);

  double **xRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);
  double **yRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nvirA,calc_info_.noccA*
    calc_info_.noccA*calc_info_.noccB,1.0,&(tRBAA[0][0]),calc_info_.noccA*
    calc_info_.noccA*calc_info_.noccB,&(thetaRBAA[0][0]),calc_info_.noccA*
    calc_info_.noccA*calc_info_.noccB,0.0,&(xRR[0][0]),calc_info_.nvirA);

  C_DGEMV('n',calc_info_.nvirA*calc_info_.nvirA,calc_info_.nrio,1.0,
    B_p_RR[0],calc_info_.nrio,calc_info_.diagBB,1,0.0,&(yRR[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.nvirA*calc_info_.nvirA,xRR[0],1,yRR[0],1);

  free_block(xRR);
  free_block(yRR);

  double **B_BB_p = get_BB_ints(1);
  double **B_p_BB = block_matrix(calc_info_.nrio,calc_info_.noccB*
    calc_info_.noccB);
    
  for(int b=0,bbp=0; b<calc_info_.noccB; b++) {
    for(int bp=0; bp<calc_info_.noccB; bp++,bbp++) {
      for(int P=0; P<calc_info_.nrio; P++) {
        B_p_BB[P][bbp] = B_BB_p[bbp][P];
  }}}
  
  free_block(B_BB_p);

  double **yRBB = block_matrix(calc_info_.nvirA,calc_info_.noccB*
    calc_info_.noccB);
  double **D_p = block_matrix(calc_info_.nvirA,calc_info_.nrio);

  for(int r1=0; r1<calc_info_.nvirA; r1++) {
      C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.noccB,
        calc_info_.noccA*calc_info_.noccA,1.0,&(tRBAA[0][0]),
        calc_info_.noccA*calc_info_.noccA,&(thetaRBAA[r1*calc_info_.noccB][0]),
        calc_info_.noccA*calc_info_.noccA,0.0,&(yRBB[0][0]),calc_info_.noccB);
      C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nrio,calc_info_.noccB*
        calc_info_.noccB,1.0,&(yRBB[0][0]),calc_info_.noccB*calc_info_.noccB,
        &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.noccB,0.0,&(D_p[0][0]),
        calc_info_.nrio);
      energy -= 2.0*C_DDOT(calc_info_.nvirA*calc_info_.nrio,
        B_p_RR[r1*calc_info_.nvirA],1,D_p[0],1);
  }

  free_block(tRBAA);
  free_block(thetaRBAA);
  free_block(D_p);
  free_block(B_p_BB);
  free_block(B_p_RR);
  free_block(yRBB);

  return(-energy);
}

double SAPT2::exch102_k11u_3()
{
  double energy=0.0;

  double **temp_thetaBSBS = read_IJKL(PSIF_SAPT_AMPS,
    "T BSBS Antisym Amplitudes",calc_info_.noccB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.nvirB);

  double **thetaSSBB = block_matrix(calc_info_.nvirB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.noccB);

  for(int b1=0,b1s1=0; b1<calc_info_.noccB; b1++) {
    for(int s1=0; s1<calc_info_.nvirB; s1++,b1s1++) {
      for(int b2=0,b2s2=0; b2<calc_info_.noccB; b2++) {
        for(int s2=0; s2<calc_info_.nvirB; s2++,b2s2++) {
          int b1b2 = b1*calc_info_.noccB+b2;
          int s1s2 = s1*calc_info_.nvirB+s2;
          thetaSSBB[s1s2][b1b2] = temp_thetaBSBS[b1s1][b2s2];
  }}}}

  free_block(temp_thetaBSBS);

  double **thetaSABB = block_matrix(calc_info_.nvirB*calc_info_.noccA,
    calc_info_.noccB*calc_info_.noccB);

  for(int s1=0; s1<calc_info_.nvirB; s1++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.noccB,
      calc_info_.nvirB,1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),
      calc_info_.nmo,&(thetaSSBB[s1*calc_info_.nvirB][0]),
      calc_info_.noccB*calc_info_.noccB,0.0,
      &(thetaSABB[s1*calc_info_.noccA][0]),calc_info_.noccB*calc_info_.noccB);
  }

  free_block(thetaSSBB);

  double **temp_tBSBS = read_IJKL(PSIF_SAPT_AMPS,"T BSBS Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);

  double **tSSBB = block_matrix(calc_info_.nvirB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.noccB);

  for(int b1=0,b1s1=0; b1<calc_info_.noccB; b1++) {
    for(int s1=0; s1<calc_info_.nvirB; s1++,b1s1++) {
      for(int b2=0,b2s2=0; b2<calc_info_.noccB; b2++) {
        for(int s2=0; s2<calc_info_.nvirB; s2++,b2s2++) {
          int b1b2 = b1*calc_info_.noccB+b2;
          int s1s2 = s1*calc_info_.nvirB+s2;
          tSSBB[s1s2][b1b2] = temp_tBSBS[b1s1][b2s2];
  }}}}

  free_block(temp_tBSBS);

  double **B_AS_p = get_AS_ints(1);
  double **B_p_AS = block_matrix(calc_info_.nrio,calc_info_.noccA*
    calc_info_.nvirB);

  for(int a=0,as=0; a<calc_info_.noccA; a++) {
    for(int s=0; s<calc_info_.nvirB; s++,as++) {
      for(int P=0; P<calc_info_.nrio; P++) {
        B_p_AS[P][as] = B_AS_p[as][P];
  }}}

  free_block(B_AS_p);

  double **B_p_SS = get_SS_ints(1);

  double **yAS = block_matrix(calc_info_.nvirB,calc_info_.noccA*
    calc_info_.nvirB);
  double **C_p = block_matrix(calc_info_.nvirB,calc_info_.nrio);

  for (int s3=0; s3<calc_info_.nvirB; s3++) {
    C_DGEMM('N','T',calc_info_.nvirB*calc_info_.noccA,calc_info_.nvirB,
      calc_info_.noccB*calc_info_.noccB,1.0,&(thetaSABB[0][0]),
      calc_info_.noccB*calc_info_.noccB,&(tSSBB[s3*calc_info_.nvirB][0]),
      calc_info_.noccB*calc_info_.noccB,0.0,&(yAS[0][0]),calc_info_.nvirB);
    C_DGEMM('N','T',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccA*
      calc_info_.nvirB,1.0,&(yAS[0][0]),calc_info_.noccA*calc_info_.nvirB,
      &(B_p_AS[0][0]),calc_info_.noccA*calc_info_.nvirB,0.0,&(C_p[0][0]),
      calc_info_.nrio);
    energy += 2.0*C_DDOT(calc_info_.nvirB*calc_info_.nrio,
      B_p_SS[s3*calc_info_.nvirB],1,C_p[0],1);
  }

  free_block(C_p);
  free_block(yAS);
  free_block(B_p_AS);

  double **tSABB = block_matrix(calc_info_.nvirB*calc_info_.noccA,
    calc_info_.noccB*calc_info_.noccB);

  for(int s1=0; s1<calc_info_.nvirB; s1++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.noccB,
      calc_info_.nvirB,1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),
      calc_info_.nmo,&(tSSBB[s1*calc_info_.nvirB][0]),
      calc_info_.noccB*calc_info_.noccB,0.0,
      &(tSABB[s1*calc_info_.noccA][0]),calc_info_.noccB*calc_info_.noccB);
  }

  free_block(tSSBB);

  double **xSS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);
  double **ySS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.nvirB,calc_info_.nvirB,calc_info_.noccB*
    calc_info_.noccB*calc_info_.noccA,1.0,&(tSABB[0][0]),calc_info_.noccB*
    calc_info_.noccB*calc_info_.noccA,&(thetaSABB[0][0]),calc_info_.noccB*
    calc_info_.noccB*calc_info_.noccA,0.0,&(xSS[0][0]),calc_info_.nvirB);

  C_DGEMV('n',calc_info_.nvirB*calc_info_.nvirB,calc_info_.nrio,1.0,
    B_p_SS[0],calc_info_.nrio,calc_info_.diagAA,1,0.0,&(ySS[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.nvirB*calc_info_.nvirB,xSS[0],1,ySS[0],1);

  free_block(xSS);
  free_block(ySS);

  double **B_AA_p = get_AA_ints(1);
  double **B_p_AA = block_matrix(calc_info_.nrio,calc_info_.noccA*
    calc_info_.noccA);

  for(int a=0,aap=0; a<calc_info_.noccA; a++) {
    for(int ap=0; ap<calc_info_.noccA; ap++,aap++) {
      for(int P=0; P<calc_info_.nrio; P++) {
        B_p_AA[P][aap] = B_AA_p[aap][P];
  }}}
  
  free_block(B_AA_p);

  double **ySAA = block_matrix(calc_info_.nvirB,calc_info_.noccA*
    calc_info_.noccA);
  double **D_p = block_matrix(calc_info_.nvirB,calc_info_.nrio);

  for(int s1=0; s1<calc_info_.nvirB; s1++) {
    C_DGEMM('N','T',calc_info_.nvirB*calc_info_.noccA,calc_info_.noccA,
      calc_info_.noccB*calc_info_.noccB,1.0,&(tSABB[0][0]),
      calc_info_.noccB*calc_info_.noccB,&(thetaSABB[s1*calc_info_.noccA][0]),
      calc_info_.noccB*calc_info_.noccB,0.0,&(ySAA[0][0]),calc_info_.noccA);
    C_DGEMM('N','T',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccA*
      calc_info_.noccA,1.0,&(ySAA[0][0]),calc_info_.noccA*calc_info_.noccA,
      &(B_p_AA[0][0]),calc_info_.noccA*calc_info_.noccA,0.0,&(D_p[0][0]),
      calc_info_.nrio);
    energy -= 2.0*C_DDOT(calc_info_.nvirB*calc_info_.nrio,
      B_p_SS[s1*calc_info_.nvirB],1,D_p[0],1);
  }

  free_block(tSABB);
  free_block(thetaSABB);
  free_block(D_p);
  free_block(B_p_AA);
  free_block(B_p_SS);
  free_block(ySAA);

  return(-energy);
}

double SAPT2::exch120_k11u_4()
{
  double energy=0.0;

  double *tARAR = init_array(calc_info_.noccA*calc_info_.nvirA*
    calc_info_.noccA*calc_info_.nvirA);
  double *thetaARAR = init_array(calc_info_.noccA*calc_info_.nvirA*
    calc_info_.noccA*calc_info_.nvirA);

  psio_->read_entry(PSIF_SAPT_AMPS,"T ARAR Amplitudes",(char *) &(tARAR[0]),
    sizeof(double)*calc_info_.noccA*calc_info_.nvirA*calc_info_.noccA*
    calc_info_.nvirA);

  psio_->read_entry(PSIF_SAPT_AMPS,"T ARAR Antisym Amplitudes",
    (char *) &(thetaARAR[0]),sizeof(double)*calc_info_.noccA*calc_info_.nvirA*
    calc_info_.noccA*calc_info_.nvirA);

  double **xRA = block_matrix(calc_info_.nvirA,calc_info_.noccA);

  for(int a=0; a<calc_info_.noccA; a++) {
    for(int r=0; r<calc_info_.nvirA; r++) {
      C_DCOPY(calc_info_.noccA*calc_info_.nvirA,&(tARAR[a*calc_info_.nvirA*
        calc_info_.noccA*calc_info_.nvirA+r]),calc_info_.nvirA,xRA[0],1);
      for(int a1=0; a1<calc_info_.noccA; a1++) {
        int aa1 = a*calc_info_.noccA+a1;
        int r1r = r;
        int aa1r1r = aa1*calc_info_.nvirA*calc_info_.nvirA+r1r;
        C_DCOPY(calc_info_.nvirA,&(xRA[0][a1]),calc_info_.noccA,
          &(tARAR[aa1r1r]),calc_info_.nvirA);
  }}}

  for(int a=0; a<calc_info_.noccA; a++) {
    for(int r=0; r<calc_info_.nvirA; r++) {
      C_DCOPY(calc_info_.noccA*calc_info_.nvirA,&(thetaARAR[a*calc_info_.nvirA*
        calc_info_.noccA*calc_info_.nvirA+r]),calc_info_.nvirA,xRA[0],1);
      for(int a1=0; a1<calc_info_.noccA; a1++) {
        int aa1 = a*calc_info_.noccA+a1;
        int r1r = r;
        int aa1r1r = aa1*calc_info_.nvirA*calc_info_.nvirA+r1r;
        C_DCOPY(calc_info_.nvirA,&(xRA[0][a1]),calc_info_.noccA,
          &(thetaARAR[aa1r1r]),calc_info_.nvirA);
  }}}

  free_block(xRA);

  double *AAAA = init_array(calc_info_.noccA*calc_info_.noccA*calc_info_.noccA*
    calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.noccA,calc_info_.noccA*
    calc_info_.noccA,calc_info_.nvirA*calc_info_.nvirA,1.0,thetaARAR,
    calc_info_.nvirA*calc_info_.nvirA,tARAR,calc_info_.nvirA*calc_info_.nvirA,
    0.0,AAAA,calc_info_.noccA*calc_info_.noccA);

  free(tARAR);
  free(thetaARAR);

  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  for(int a=0; a<calc_info_.noccA; a++) {
    for(int aa=0; aa<calc_info_.noccA; aa++) {
      C_DCOPY(calc_info_.noccA*calc_info_.noccA,&(AAAA[a*calc_info_.noccA*
        calc_info_.noccA*calc_info_.noccA+aa]),calc_info_.noccA,xAA[0],1);
      for(int a1=0; a1<calc_info_.noccA; a1++) {
        int aa1 = a*calc_info_.noccA+a1;
        int aa1a = aa1*calc_info_.noccA*calc_info_.noccA+aa;
        C_DCOPY(calc_info_.noccA,&(xAA[0][a1]),calc_info_.noccA,
          &(AAAA[aa1a]),calc_info_.noccA);
  }}}

  free_block(xAA);

  double **B_p_AA = get_AA_ints(1);
  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,
    calc_info_.noccA*calc_info_.noccA,1.0,AAAA,calc_info_.noccA*
    calc_info_.noccA,&(B_p_AA[0][0]),calc_info_.nrio,0.0,
    &(C_p_AA[0][0]),calc_info_.nrio);

  free(AAAA);
  free_block(B_p_AA);

  double **B_p_AB = get_AB_ints(2);
  double **D_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_AB[a*calc_info_.noccB][0]),
      calc_info_.nrio,0.0,&(D_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  energy += 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    C_p_AA[0],1,D_p_AA[0],1);

  free_block(B_p_AB);
  free_block(D_p_AA);

  double *X = init_array(calc_info_.nrio);
  double **sAA = block_matrix(calc_info_.noccA,calc_info_.noccA);
  
  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(sAA[0][0]),calc_info_.noccA);
  
  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    C_p_AA[0],calc_info_.nrio,sAA[0],1,0.0,X,1);

  energy += 4.0*C_DDOT(calc_info_.nrio,X,1,calc_info_.diagBB,1);

  double **B_p_BB = get_BB_ints(1);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(B_p_BB);

  double **E_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(C_p_AB[a*calc_info_.noccB][0]),
      calc_info_.nrio,0.0,&(E_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    C_p_AA[0],1,E_p_AA[0],1);

  free_block(C_p_AB);
  free_block(C_p_AA);
  free_block(E_p_AA);

  return(-energy);
}

double SAPT2::exch102_k11u_4()
{
  double energy=0.0;

  double *tBSBS = init_array(calc_info_.noccB*calc_info_.nvirB*
    calc_info_.noccB*calc_info_.nvirB);
  double *thetaBSBS = init_array(calc_info_.noccB*calc_info_.nvirB*
    calc_info_.noccB*calc_info_.nvirB);

  psio_->read_entry(PSIF_SAPT_AMPS,"T BSBS Amplitudes",(char *) &(tBSBS[0]),
    sizeof(double)*calc_info_.noccB*calc_info_.nvirB*calc_info_.noccB*
    calc_info_.nvirB);

  psio_->read_entry(PSIF_SAPT_AMPS,"T BSBS Antisym Amplitudes",
    (char *) &(thetaBSBS[0]),sizeof(double)*calc_info_.noccB*calc_info_.nvirB*
    calc_info_.noccB*calc_info_.nvirB);

  double **xSB = block_matrix(calc_info_.nvirB,calc_info_.noccB);

  for(int b=0; b<calc_info_.noccB; b++) {
    for(int s=0; s<calc_info_.nvirB; s++) {
      C_DCOPY(calc_info_.noccB*calc_info_.nvirB,&(tBSBS[b*calc_info_.nvirB*
        calc_info_.noccB*calc_info_.nvirB+s]),calc_info_.nvirB,xSB[0],1);
      for(int b1=0; b1<calc_info_.noccB; b1++) {
        int bb1 = b*calc_info_.noccB+b1;
        int s1s = s;
        int bb1s1s = bb1*calc_info_.nvirB*calc_info_.nvirB+s1s;
        C_DCOPY(calc_info_.nvirB,&(xSB[0][b1]),calc_info_.noccB,
          &(tBSBS[bb1s1s]),calc_info_.nvirB);
  }}}

  for(int b=0; b<calc_info_.noccB; b++) {
    for(int s=0; s<calc_info_.nvirB; s++) {
      C_DCOPY(calc_info_.noccB*calc_info_.nvirB,&(thetaBSBS[b*calc_info_.nvirB*
        calc_info_.noccB*calc_info_.nvirB+s]),calc_info_.nvirB,xSB[0],1);
      for(int b1=0; b1<calc_info_.noccB; b1++) {
        int bb1 = b*calc_info_.noccB+b1;
        int s1s = s;
        int bb1s1s = bb1*calc_info_.nvirB*calc_info_.nvirB+s1s;
        C_DCOPY(calc_info_.nvirB,&(xSB[0][b1]),calc_info_.noccB,
          &(thetaBSBS[bb1s1s]),calc_info_.nvirB);
  }}}

  free_block(xSB);

  double *BBBB = init_array(calc_info_.noccB*calc_info_.noccB*calc_info_.noccB*
    calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccB*calc_info_.noccB,calc_info_.noccB*
    calc_info_.noccB,calc_info_.nvirB*calc_info_.nvirB,1.0,thetaBSBS,
    calc_info_.nvirB*calc_info_.nvirB,tBSBS,calc_info_.nvirB*calc_info_.nvirB,
    0.0,BBBB,calc_info_.noccB*calc_info_.noccB);

  free(tBSBS);
  free(thetaBSBS);

  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  for(int b=0; b<calc_info_.noccB; b++) {
    for(int bb=0; bb<calc_info_.noccB; bb++) {
      C_DCOPY(calc_info_.noccB*calc_info_.noccB,&(BBBB[b*calc_info_.noccB*
        calc_info_.noccB*calc_info_.noccB+bb]),calc_info_.noccB,xBB[0],1);
      for(int b1=0; b1<calc_info_.noccB; b1++) {
        int bb1 = b*calc_info_.noccB+b1;
        int bb1b = bb1*calc_info_.noccB*calc_info_.noccB+bb;
        C_DCOPY(calc_info_.noccB,&(xBB[0][b1]),calc_info_.noccB,
          &(BBBB[bb1b]),calc_info_.noccB);
  }}}

  free_block(xBB);

  double **B_p_BB = get_BB_ints(1);
  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,
    calc_info_.noccB*calc_info_.noccB,1.0,BBBB,calc_info_.noccB*
    calc_info_.noccB,&(B_p_BB[0][0]),calc_info_.nrio,0.0,
    &(C_p_BB[0][0]),calc_info_.nrio);

  free(BBBB);
  free_block(B_p_BB);

  double **B_p_AB = get_AB_ints(1);
  double **D_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(D_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  energy += 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    C_p_BB[0],1,D_p_BB[0],1);

  free_block(B_p_AB);
  free_block(D_p_BB);

  double *X = init_array(calc_info_.nrio);
  double **sBB = block_matrix(calc_info_.noccB,calc_info_.noccB);
  
  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(sBB[0][0]),calc_info_.noccB);
  
  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    C_p_BB[0],calc_info_.nrio,sBB[0],1,0.0,X,1);

  energy += 4.0*C_DDOT(calc_info_.nrio,X,1,calc_info_.diagAA,1);

  double **B_p_AA = get_AA_ints(1);
  double **C_p_BA = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,&(C_p_BA[0][0]),
    calc_info_.noccA*calc_info_.nrio);

  free_block(B_p_AA);

  double **E_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(C_p_BA[b*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(E_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio);
  }

  energy -= 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    C_p_BB[0],1,E_p_BB[0],1);

  free_block(C_p_BA);
  free_block(C_p_BB);
  free_block(E_p_BB);

  return(-energy);
}

double SAPT2::exch120_k11u_5()
{
  double energy=0.0;

  double **theta_p_AR = read_IJKL(PSIF_SAPT_AMPS,"Theta(AR) AR",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio);
  double **thetaARAR = read_IJKL(PSIF_SAPT_AMPS,"T ARAR Antisym Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

  double **T_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,
    calc_info_.noccA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA,&(theta_p_AR[0][0]),calc_info_.nrio,0.0,
    &(T_p_AR[0][0]),calc_info_.nrio);

  free_block(theta_p_AR);
  free_block(thetaARAR);

  double **C_p_BR = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(T_p_AR[0][0]),calc_info_.nvirA*calc_info_.nrio,0.0,&(C_p_BR[0][0]),
    calc_info_.nvirA*calc_info_.nrio);

  double **B_p_RB = get_RB_ints(1);

  for (int r=0,rb=0; r<calc_info_.nvirA; r++){
    for (int b=0; b<calc_info_.noccB; b++,rb++){
      int br = b*calc_info_.nvirA+r;
      energy += C_DDOT(calc_info_.nrio,C_p_BR[br],1,B_p_RB[rb],1);
  }}

  free_block(B_p_RB);
  free_block(C_p_BR);

  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(T_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,0.0,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  double **B_p_AB = get_AB_ints(2);

  energy += C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    B_p_AB[0],1,C_p_AB[0],1);

  free_block(B_p_AB);

  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,
    &(C_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio);

  free_block(C_p_AB);

  double **B_p_BB = get_BB_ints(1);

  energy -= 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    B_p_BB[0],1,C_p_BB[0],1);

  free_block(B_p_BB);
  free_block(C_p_BB);

  double **xAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);
  double **yAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(xAR[0][0]),
    calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    T_p_AR[0],calc_info_.nrio,calc_info_.diagBB,1,0.0,&(yAR[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,xAR[0],1,yAR[0],1);

  free_block(xAR);
  free_block(yAR);
  free_block(T_p_AR);

  return(-2.0*energy);
}

double SAPT2::exch102_k11u_5()
{
  double energy=0.0;

  double **theta_p_BS = read_IJKL(PSIF_SAPT_AMPS,"Theta(BS) BS",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio);
  double **thetaBSBS = read_IJKL(PSIF_SAPT_AMPS,"T BSBS Antisym Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);

  double **T_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,
    calc_info_.noccB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,&(theta_p_BS[0][0]),calc_info_.nrio,0.0,
    &(T_p_BS[0][0]),calc_info_.nrio);

  free_block(theta_p_BS);
  free_block(thetaBSBS);

  double **C_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(T_p_BS[0][0]),calc_info_.nvirB*calc_info_.nrio,0.0,&(C_p_AS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  double **B_p_AS = get_AS_ints(1);

  energy += C_DDOT(calc_info_.noccA*calc_info_.nvirB*calc_info_.nrio,
    C_p_AS[0],1,B_p_AS[0],1);

  free_block(B_p_AS);
  free_block(C_p_AS);

  double **C_p_BA = block_matrix(calc_info_.noccB*calc_info_.noccA,
    calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(T_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(C_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio);
  }

  double **B_p_AB = get_AB_ints(1);

  for(int a=0,ab=0; a<calc_info_.noccA; a++) {
    for(int b=0; b<calc_info_.noccB; b++,ab++) {
      int ba = b*calc_info_.noccA+a;
      energy += C_DDOT(calc_info_.nrio,B_p_AB[ab],1,C_p_BA[ba],1);
  }}

  free_block(B_p_AB);

  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_BA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,
    &(C_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio);

  free_block(C_p_BA);

  double **B_p_AA = get_AA_ints(1);

  energy -= 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    B_p_AA[0],1,C_p_AA[0],1);

  free_block(B_p_AA);
  free_block(C_p_AA);

  double **xBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);
  double **yBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(xBS[0][0]),
    calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    T_p_BS[0],calc_info_.nrio,calc_info_.diagAA,1,0.0,&(yBS[0][0]),1);

  energy += 4.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,xBS[0],1,yBS[0],1);

  free_block(xBS);
  free_block(yBS);
  free_block(T_p_BS);

  return(-2.0*energy);
}

double SAPT2::exch120_k11u_6()
{
  double energy=0.0;

  double *T_ARAR = init_array(calc_info_.noccA*calc_info_.nvirA*
    calc_info_.noccA*calc_info_.nvirA);
  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,"T ARAR Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA,3.0,&(tARAR[0][0]),
    calc_info_.noccA*calc_info_.nvirA,&(tARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA,0.0,T_ARAR,calc_info_.noccA*calc_info_.nvirA);

  free_block(tARAR);

  double **thetaARAR = read_IJKL(PSIF_SAPT_AMPS,"T ARAR Antisym Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

  double *X = init_array(calc_info_.nvirA);

  for(int a=0; a<calc_info_.noccA; a++) {
  for(int a1=0; a1<=a; a1++) {
    for(int r=0; r<calc_info_.nvirA; r++) {
      int ar = a*calc_info_.nvirA+r;
      int a1r = a1*calc_info_.nvirA+r;
      C_DCOPY(calc_info_.nvirA,&(thetaARAR[ar][a1*calc_info_.nvirA]),1,X,1);
      C_DCOPY(calc_info_.nvirA,&(thetaARAR[a1r][a*calc_info_.nvirA]),1,
        &(thetaARAR[ar][a1*calc_info_.nvirA]),1);
      C_DCOPY(calc_info_.nvirA,X,1,&(thetaARAR[a1r][a*calc_info_.nvirA]),1);
  }}}

  free(X);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),
    calc_info_.noccA*calc_info_.nvirA,&(thetaARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA,1.0,T_ARAR,calc_info_.noccA*calc_info_.nvirA);

  free_block(thetaARAR);

  double **xRA = block_matrix(calc_info_.nvirA,calc_info_.noccA);
  
  for(int a=0; a<calc_info_.noccA; a++) {
    for(int r=0; r<calc_info_.nvirA; r++) {
      C_DCOPY(calc_info_.noccA*calc_info_.nvirA,&(T_ARAR[a*calc_info_.nvirA*
        calc_info_.noccA*calc_info_.nvirA+r]),calc_info_.nvirA,xRA[0],1);
      for(int a1=0; a1<calc_info_.noccA; a1++) {
        int aa1 = a*calc_info_.noccA+a1;
        int r1r = r;
        int aa1r1r = aa1*calc_info_.nvirA*calc_info_.nvirA+r1r;
        C_DCOPY(calc_info_.nvirA,&(xRA[0][a1]),calc_info_.noccA,
          &(T_ARAR[aa1r1r]),calc_info_.nvirA);
  }}}

  free_block(xRA);

  double **B_p_RR = get_RR_ints(1);
  double **T_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,
    calc_info_.nvirA*calc_info_.nvirA,1.0,T_ARAR,calc_info_.nvirA*
    calc_info_.nvirA,&(B_p_RR[0][0]),calc_info_.nrio,0.0,&(T_p_AA[0][0]),
    calc_info_.nrio);

  free_block(B_p_RR);

  double **B_p_AA = get_AA_ints(1);
  double **T_p_RR = block_matrix(calc_info_.nvirA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.nvirA*calc_info_.nvirA,calc_info_.nrio,
    calc_info_.noccA*calc_info_.noccA,1.0,T_ARAR,calc_info_.nvirA*
    calc_info_.nvirA,&(B_p_AA[0][0]),calc_info_.nrio,0.0,&(T_p_RR[0][0]),
    calc_info_.nrio);

  free(T_ARAR);
  free_block(B_p_AA);

  double **B_p_AB = get_AB_ints(2);
  double **C_p_BA = block_matrix(calc_info_.noccB*calc_info_.noccA,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(T_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,&(C_p_BA[0][0]),
    calc_info_.noccA*calc_info_.nrio);

  for(int a=0, ab=0; a<calc_info_.noccA; a++) {
    for(int b=0; b<calc_info_.noccB; b++, ab++) {
      int ba = b*calc_info_.noccA+a;
      energy -= C_DDOT(calc_info_.nrio,&(C_p_BA[ba][0]),1,&(B_p_AB[ab][0]),1);
  }}

  free_block(B_p_AB);

  double **B_p_BB = get_BB_ints(1);
  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(C_p_BA[b*calc_info_.noccA][0]),
      calc_info_.nrio,0.0,&(C_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio);
  }

  energy += C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    &(C_p_BB[0][0]),1,&(B_p_BB[0][0]),1);

  free_block(C_p_BB);
  free_block(C_p_BA);

  double **B_p_RB = get_RB_ints(1);
  double **C_p_BR = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,&(T_p_RR[0][0]),calc_info_.nvirA*calc_info_.nrio,0.0,
    &(C_p_BR[0][0]),calc_info_.nvirA*calc_info_.nrio);

  for(int r=0, rb=0; r<calc_info_.nvirA; r++) {
    for(int b=0; b<calc_info_.noccB; b++, rb++) {
      int br = b*calc_info_.nvirA+r;
      energy -= C_DDOT(calc_info_.nrio,B_p_RB[rb],1,C_p_BR[br],1);
  }}

  free_block(B_p_RB);

  double **D_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(C_p_BR[b*calc_info_.nvirA][0]),calc_info_.nrio,0.0,
      &(D_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio);
  }

  energy += C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,
    &(D_p_BB[0][0]),1,&(B_p_BB[0][0]),1);

  free_block(B_p_BB);
  free_block(C_p_BR);
  free_block(D_p_BB);

  double **sAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(sAA[0][0]),calc_info_.noccA);

  double **sRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(sRR[0][0]),
    calc_info_.nvirA);

  X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    T_p_AA[0],calc_info_.nrio,&(sAA[0][0]),1,0.0,X,1);

  energy -= 2.0*C_DDOT(calc_info_.nrio,X,1,calc_info_.diagBB,1);
  
  C_DGEMV('t',calc_info_.nvirA*calc_info_.nvirA,calc_info_.nrio,1.0,
    T_p_RR[0],calc_info_.nrio,&(sRR[0][0]),1,0.0,X,1);

  energy -= 2.0*C_DDOT(calc_info_.nrio,X,1,calc_info_.diagBB,1);

  free(X);
  free_block(sAA);
  free_block(sRR);
  free_block(T_p_AA);
  free_block(T_p_RR);

  return(-energy);
}

double SAPT2::exch102_k11u_6()
{
  double energy=0.0;

  double *T_BSBS = init_array(calc_info_.noccB*calc_info_.nvirB*
    calc_info_.noccB*calc_info_.nvirB);
  double **tBSBS = read_IJKL(PSIF_SAPT_AMPS,"T BSBS Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB,3.0,&(tBSBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB,&(tBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,0.0,T_BSBS,calc_info_.noccB*calc_info_.nvirB);

  free_block(tBSBS);

  double **thetaBSBS = read_IJKL(PSIF_SAPT_AMPS,"T BSBS Antisym Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);

  double *X = init_array(calc_info_.nvirB);

  for(int b=0; b<calc_info_.noccB; b++) {
  for(int b1=0; b1<=b; b1++) {
    for(int s=0; s<calc_info_.nvirB; s++) {
      int bs = b*calc_info_.nvirB+s;
      int b1s = b1*calc_info_.nvirB+s;
      C_DCOPY(calc_info_.nvirB,&(thetaBSBS[bs][b1*calc_info_.nvirB]),1,X,1);
      C_DCOPY(calc_info_.nvirB,&(thetaBSBS[b1s][b*calc_info_.nvirB]),1,
        &(thetaBSBS[bs][b1*calc_info_.nvirB]),1);
      C_DCOPY(calc_info_.nvirB,X,1,&(thetaBSBS[b1s][b*calc_info_.nvirB]),1);
  }}}

  free(X);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB,&(thetaBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,1.0,T_BSBS,calc_info_.noccB*calc_info_.nvirB);

  free_block(thetaBSBS);

  double **xSB = block_matrix(calc_info_.nvirB,calc_info_.noccB);
  
  for(int b=0; b<calc_info_.noccB; b++) {
    for(int s=0; s<calc_info_.nvirB; s++) {
      C_DCOPY(calc_info_.noccB*calc_info_.nvirB,&(T_BSBS[b*calc_info_.nvirB*
        calc_info_.noccB*calc_info_.nvirB+s]),calc_info_.nvirB,xSB[0],1);
      for(int b1=0; b1<calc_info_.noccB; b1++) {
        int bb1 = b*calc_info_.noccB+b1;
        int s1s = s;
        int bb1s1s = bb1*calc_info_.nvirB*calc_info_.nvirB+s1s;
        C_DCOPY(calc_info_.nvirB,&(xSB[0][b1]),calc_info_.noccB,
          &(T_BSBS[bb1s1s]),calc_info_.nvirB);
  }}}

  free_block(xSB);

  double **B_p_SS = get_SS_ints(1);
  double **T_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,
    calc_info_.nvirB*calc_info_.nvirB,1.0,T_BSBS,calc_info_.nvirB*
    calc_info_.nvirB,&(B_p_SS[0][0]),calc_info_.nrio,0.0,&(T_p_BB[0][0]),
    calc_info_.nrio);

  free_block(B_p_SS);

  double **B_p_BB = get_BB_ints(1);
  double **T_p_SS = block_matrix(calc_info_.nvirB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.nvirB*calc_info_.nvirB,calc_info_.nrio,
    calc_info_.noccB*calc_info_.noccB,1.0,T_BSBS,calc_info_.nvirB*
    calc_info_.nvirB,&(B_p_BB[0][0]),calc_info_.nrio,0.0,&(T_p_SS[0][0]),
    calc_info_.nrio);

  free(T_BSBS);
  free_block(B_p_BB);

  double **B_p_AB = get_AB_ints(1);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(T_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  energy -= C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,
    &(C_p_AB[0][0]),1,&(B_p_AB[0][0]),1);

  free_block(B_p_AB);

  double **B_p_AA = get_AA_ints(1);
  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(C_p_AB[a*calc_info_.noccB][0]),
      calc_info_.nrio,0.0,&(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  energy += C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    &(C_p_AA[0][0]),1,&(B_p_AA[0][0]),1);

  free_block(C_p_AA);
  free_block(C_p_AB);

  double **B_p_AS = get_AS_ints(1);
  double **C_p_AS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.nvirB,1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,&(T_p_SS[0][0]),calc_info_.nvirB*calc_info_.nrio,0.0,
    &(C_p_AS[0][0]),calc_info_.nvirB*calc_info_.nrio);

  energy -= C_DDOT(calc_info_.noccA*calc_info_.nvirB*calc_info_.nrio,
    B_p_AS[0],1,C_p_AS[0],1);

  free_block(B_p_AS);

  double **D_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(C_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(D_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  energy += C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,
    &(D_p_AA[0][0]),1,&(B_p_AA[0][0]),1);

  free_block(B_p_AA);
  free_block(C_p_AS);
  free_block(D_p_AA);

  double **sBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(sBB[0][0]),calc_info_.noccB);

  double **sSS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(sSS[0][0]),
    calc_info_.nvirB);

  X = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    T_p_BB[0],calc_info_.nrio,&(sBB[0][0]),1,0.0,X,1);

  energy -= 2.0*C_DDOT(calc_info_.nrio,X,1,calc_info_.diagAA,1);
  
  C_DGEMV('t',calc_info_.nvirB*calc_info_.nvirB,calc_info_.nrio,1.0,
    T_p_SS[0],calc_info_.nrio,&(sSS[0][0]),1,0.0,X,1);

  energy -= 2.0*C_DDOT(calc_info_.nrio,X,1,calc_info_.diagAA,1);

  free(X);
  free_block(sBB);
  free_block(sSS);
  free_block(T_p_BB);
  free_block(T_p_SS);

  return(-energy);
}

}}
