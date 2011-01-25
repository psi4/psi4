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

void SAPT0::exch10()
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

  results_.exch10 = -2.0*(ex1+ex2+ex3-ex4-ex5+ex6);

  if (params_.print) {
    fprintf(outfile,"exch10      Energy = %18.12lf  H\n\n",results_.exch10);
    fflush(outfile);
  }
}

}}
