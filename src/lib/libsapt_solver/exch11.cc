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
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::exch11()
{
  double ex110, ex101;

  if (params_.print) {
    fprintf(outfile,"Begining Exch11 Calculation\n\n");
    fflush(outfile);
  }

  ex110 = exch110("Theta(AR) AR");

  if (params_.print) {
    fprintf(outfile,"exch110            = %18.12lf  H\n",ex110);
    fflush(outfile);
  }

  ex101 = exch101("Theta(BS) BS");

  if (params_.print) {
    fprintf(outfile,"exch101            = %18.12lf  H\n\n",ex101);
    fflush(outfile);
  }

  results_.exch11 = ex110 + ex101;
}

double SAPT2::exch110(char *theta_label)
{
  double e1, e2, e3, e4;

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,theta_label,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(T_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,0.0,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  double **B_p_AB = get_AB_ints(2);

  e1 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,B_p_AB[0],
    1,C_p_AB[0],1);

  free_block(B_p_AB);

  double **B_p_RB = get_RB_ints(1);
  double **C_p_BR = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(T_p_AR[0][0]),calc_info_.nvirA*calc_info_.nrio,0.0,
    &(C_p_BR[0][0]),calc_info_.nvirA*calc_info_.nrio);

  e2 = 0.0;

  for (int r=0, rb=0; r<calc_info_.nvirA; r++) {
    for (int b=0; b<calc_info_.noccB; b++, rb++) {
      int br = b*calc_info_.nvirA+r;
      e2 -= 2.0*C_DDOT(calc_info_.nrio,B_p_RB[rb],1,C_p_BR[br],1);
  }}

  free_block(B_p_RB);
  free_block(C_p_BR);

  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,
    &(C_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio);

  free_block(C_p_AB);

  double **B_p_BB = get_BB_ints(1);

  e3 = 4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB*calc_info_.nrio,B_p_BB[0],
    1,C_p_BB[0],1);

  free_block(B_p_BB);
  free_block(C_p_BB);

  double **X_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(X_AR[0][0]),
    calc_info_.nvirA);

  double *C_p = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(T_p_AR[0][0]),calc_info_.nrio,X_AR[0],1,0.0,C_p,1);

  e4 = -8.0*C_DDOT(calc_info_.nrio,calc_info_.diagBB,1,C_p,1);

  free(C_p);
  free_block(X_AR);
  free_block(T_p_AR);

  return(e1+e2+e3+e4);
}

double SAPT2::exch101(char *theta_label)
{
  double e1, e2, e3, e4;

  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,theta_label,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);
  double **C_p_BA = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(T_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,0.0,
      &(C_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio);
  }

  double **B_p_AB = get_AB_ints(1);

  e1 = 0.0;

  for (int a=0, ab=0; a<calc_info_.noccA; a++) {
    for (int b=0; b<calc_info_.noccB; b++, ab++) {
      int ba = b*calc_info_.noccA+a;
      e1 -= 2.0*C_DDOT(calc_info_.nrio,B_p_AB[ab],1,C_p_BA[ba],1);
  }}

  free_block(B_p_AB);

  double **B_p_AS = get_AS_ints(1);
  double **C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AS[0][0]),calc_info_.nvirB*calc_info_.nrio,0.0,
    &(C_p_BS[0][0]),calc_info_.nvirB*calc_info_.nrio);

  e2 = -2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,
    T_p_BS[0],1,C_p_BS[0],1);

  free_block(B_p_AS);
  free_block(C_p_BS);

  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_BA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,
    &(C_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio);

  free_block(C_p_BA);

  double **B_p_AA = get_AA_ints(1);

  e3 = 4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA*calc_info_.nrio,B_p_AA[0],
    1,C_p_AA[0],1);

  free_block(B_p_AA);
  free_block(C_p_AA);

  double **X_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(X_BS[0][0]),
    calc_info_.nvirB);

  double *C_p = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(T_p_BS[0][0]),calc_info_.nrio,X_BS[0],1,0.0,C_p,1);

  e4 = -8.0*C_DDOT(calc_info_.nrio,calc_info_.diagAA,1,C_p,1);

  free(C_p);
  free_block(X_BS);
  free_block(T_p_BS);

  return(e1+e2+e3+e4);
}

}}
