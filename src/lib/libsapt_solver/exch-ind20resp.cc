/* This function calculates the Exch-Ind20 energy */

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
#include "sapt0.h"
#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT0::exch_ind20()
{
  double exindA_B, exindB_A;

  if (params_.print) {
    fprintf(outfile,"Begining Exch-Ind20 Calculation\n\n");
    fflush(outfile);
  }

  exindA_B = exch_ind20respA_B();

  if (params_.print) {
    fprintf(outfile,"exch_indrA_B       = %18.12lf  H\n",exindA_B);
    fflush(outfile);
  }

  exindB_A = exch_ind20respB_A();

  if (params_.print) {
    fprintf(outfile,"exch_indrB_A       = %18.12lf  H\n\n",exindB_A);
    fflush(outfile);
  }

  results_.exch_ind20 = exindA_B + exindB_A;
}

double SAPT0::exch_ind20respA_B()
{
  double e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13;

  if (params_.print) {
    fprintf(outfile,"Begining Exch-Ind20,resp (A<-B) Calculation\n\n");
    fflush(outfile);
  }

  double **B_p_AB = get_AB_ints(1);
  double **B_p_RB = get_RB_ints(1);
  double **Y_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB*calc_info_.nrio,
    1.0,&(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,&(B_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(Y_AR[0][0]),calc_info_.nvirA);

  e1 = C_DDOT(calc_info_.noccA*calc_info_.nvirA,&(calc_info_.CHFA[0][0]),1,
         &(Y_AR[0][0]),1);

  free_block(B_p_AB);
 
  if (params_.print) {
    fprintf(outfile,"EXCH1       Energy = %18.12lf  H\n",e1);
    fflush(outfile);
  }

  double **X_RB = block_matrix(calc_info_.nvirA,calc_info_.noccB);
  double **Y_RB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.CHFA[0][0]),calc_info_.nvirA,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(X_RB[0][0]),calc_info_.noccB);

  C_DGEMV('n',calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_RB[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,0.0,&(Y_RB[0][0]),1);

  e2 = 2.0*C_DDOT(calc_info_.nvirA*calc_info_.noccB,&(X_RB[0][0]),1,
             &(Y_RB[0][0]),1);

  if (params_.print) {
    fprintf(outfile,"EXCH2       Energy = %18.12lf  H\n",e2);
    fflush(outfile);
  }

  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,calc_info_.nvirA,
    1.0,&(calc_info_.CHFA[0][0]),calc_info_.nvirA,&(B_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(B_p_RB);

  double **Y_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);
  double **B_p_AA = get_AA_ints(1);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(Y_AB[0][0]),
      calc_info_.noccB);
  }

  e5 = 0.0;

  for (int a=0; a<calc_info_.noccA; a++) {
    e5 += C_DDOT(calc_info_.noccB,&(calc_info_.S_AB[a][0]),1,&(Y_AB[a][0]),1);
  }

  free_block(C_p_AB);
  free_block(B_p_AA);

  if (params_.print) {
    fprintf(outfile,"EXCH5       Energy = %18.12lf  H\n",e5);
    fflush(outfile);
  }

  double **B_p_AR = get_AR_ints(1);
  B_p_AB = get_AB_ints(2);

  double *C_p = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
    calc_info_.nrio,&(calc_info_.CHFA[0][0]),1,0.0,C_p,1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_AB[0][0]),
    calc_info_.nrio,C_p,1,0.0,&(Y_AB[0][0]),1);

  e3 = 0.0;

  for (int a=0; a<calc_info_.noccA; a++) {
    e3 += 2.0*C_DDOT(calc_info_.noccB,&(calc_info_.S_AB[a][0]),1,&(Y_AB[a][0]),1);
  }

  if (params_.print) {
    fprintf(outfile,"EXCH3       Energy = %18.12lf  H\n",e3);
    fflush(outfile);
  }

  memset(&(Y_RB[0][0]),'\0',sizeof(double)*calc_info_.nvirA*calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(Y_RB[0][0]),
      calc_info_.noccB);
  }

  e4 = C_DDOT(calc_info_.nvirA*calc_info_.noccB,&(X_RB[0][0]),1,&(Y_RB[0][0]),1);

  free_block(B_p_AB);

  if (params_.print) {
    fprintf(outfile,"EXCH4       Energy = %18.12lf  H\n",e4);
    fflush(outfile);
  }

  double **B_p_BB = get_BB_ints(1);

  double **X_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(X_BB[0][0]),calc_info_.noccB);

  double **Y_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_BB[0][0]),
    calc_info_.nrio,C_p,1,0.0,&(Y_BB[0][0]),1);

  free(C_p);

  e9 = 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,&(X_BB[0][0]),1,
             &(Y_BB[0][0]),1);

  if (params_.print) {
    fprintf(outfile,"EXCH9       Energy = %18.12lf  H\n",e9);
    fflush(outfile);
  }

  double **X_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(X_AA[0][0]),calc_info_.noccA);

  double **X_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccA,1.0,
    &(X_AA[0][0]),calc_info_.noccA,&(calc_info_.CHFA[0][0]),
    calc_info_.nvirA,0.0,&(X_AR[0][0]),calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,0.0,&(Y_AR[0][0]),1);

  e11 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,&(X_AR[0][0]),1,
             &(Y_AR[0][0]),1); 

  free_block(X_AR);
  free_block(Y_AR);

  if (params_.print) {
    fprintf(outfile,"EXCH11      Energy = %18.12lf  H\n",e11);
    fflush(outfile);
  }

  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,calc_info_.noccB,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  memset(&(Y_RB[0][0]),'\0',sizeof(double)*calc_info_.nvirA*calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(Y_RB[0][0]),
      calc_info_.noccB);
  }

  e13 = C_DDOT(calc_info_.nvirA*calc_info_.noccB,&(X_RB[0][0]),1,&(Y_RB[0][0]),1);

  free_block(X_RB);
  free_block(Y_RB);
	free_block(C_p_AB);
	free_block(B_p_AR);

  if (params_.print) {
    fprintf(outfile,"EXCH13      Energy = %18.12lf  H\n",e13);
    fflush(outfile);
  }

  B_p_AB = get_AB_ints(1);

  double **X_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirA,1.0,
    &(calc_info_.CHFA[0][0]),calc_info_.nvirA,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(X_AB[0][0]),calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_AB[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,0.0,&(Y_AB[0][0]),1);

  e6 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,&(X_AB[0][0]),1,
             &(Y_AB[0][0]),1);

  if (params_.print) {
    fprintf(outfile,"EXCH6       Energy = %18.12lf  H\n",e6);
    fflush(outfile);
  }

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    1.0,&(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(Y_AB[0][0]),calc_info_.noccB);

  e7 = C_DDOT(calc_info_.noccA*calc_info_.noccB,&(X_AB[0][0]),1,&(Y_AB[0][0]),1);

  free_block(Y_AB);
  free_block(B_p_AB);

  if (params_.print) {
    fprintf(outfile,"EXCH7       Energy = %18.12lf  H\n",e7);
    fflush(outfile);
  }

  B_p_AA = get_AA_ints(1);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(X_AB[0][0]),
    calc_info_.noccB,0.0,&(X_BB[0][0]),calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_BB[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,0.0,&(Y_BB[0][0]),1);

  e8 = 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,&(X_BB[0][0]),1,
             &(Y_BB[0][0]),1);

  free_block(X_BB);
  free_block(Y_BB);

  if (params_.print) {
    fprintf(outfile,"EXCH8       Energy = %18.12lf  H\n",e8);
    fflush(outfile);
  }

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(X_AB[0][0]),
    calc_info_.noccB,0.0,&(X_AA[0][0]),calc_info_.noccA);

  double **Y_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,&(B_p_AA[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,0.0,&(Y_AA[0][0]),1);

  e10 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,&(X_AA[0][0]),1,
              &(Y_AA[0][0]),1);

  free_block(X_AA);
  free_block(Y_AA);

  if (params_.print) {
    fprintf(outfile,"EXCH10      Energy = %18.12lf  H\n",e10);
    fflush(outfile);
  }

  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,1.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

	free_block(B_p_AA);
  double **D_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,calc_info_.noccB,
    1.0,&(X_AB[0][0]),calc_info_.noccB,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(D_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  e12 = C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,&(C_p_AB[0][0]),1,
          &(D_p_AB[0][0]),1);

  free_block(X_AB);
	free_block(B_p_BB);
	free_block(C_p_AB);
	free_block(D_p_AB);

  if (params_.print) 
    fprintf(outfile,"EXCH12      Energy = %18.12lf  H\n\n",e12);

  return(-2.0*(e1+e2+e3-e4-e5+e6-e7-e8-e9-e10-e11+e12+e13));
}

double SAPT0::exch_ind20respB_A()
{
  double e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13;

  if (params_.print) {
    fprintf(outfile,"Begining Exch-Ind20,resp (B<-A) Calculation\n\n");
    fflush(outfile);
  }

  double **B_p_AS = get_AS_ints(1);
  double **B_p_AB = get_AB_ints(2);
  double **Y_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,
      &(B_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,
      1.0,&(Y_BS[0][0]),calc_info_.nvirB);
  }

  e1 = C_DDOT(calc_info_.noccB*calc_info_.nvirB,&(calc_info_.CHFB[0][0]),1,
         &(Y_BS[0][0]),1);

  free_block(B_p_AB);

  if (params_.print) {
    fprintf(outfile,"EXCH1       Energy = %18.12lf  H\n",e1);
    fflush(outfile);
  }

  double **X_AS = block_matrix(calc_info_.noccA,calc_info_.nvirB);
  double **Y_AS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.CHFB[0][0]),
    calc_info_.nvirB,0.0,&(X_AS[0][0]),calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_AS[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,0.0,&(Y_AS[0][0]),1);

  e2 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirB,&(X_AS[0][0]),1,
             &(Y_AS[0][0]),1);

  if (params_.print) {
    fprintf(outfile,"EXCH2       Energy = %18.12lf  H\n",e2);
    fflush(outfile);
  }

  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(calc_info_.CHFB[0][0]),calc_info_.nvirB,&(B_p_AS[a*calc_info_.nvirB][0]),
      calc_info_.nrio,1.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  free_block(B_p_AS);

  double **Y_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);
  double **B_p_BB = get_BB_ints(1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    1.0,&(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(Y_AB[0][0]),calc_info_.noccB);

  e5 = 0.0;

  for (int a=0; a<calc_info_.noccA; a++) {
    e5 += C_DDOT(calc_info_.noccB,&(calc_info_.S_AB[a][0]),1,&(Y_AB[a][0]),1);
  }

  free_block(C_p_AB);
  free_block(B_p_BB);

  if (params_.print) {
    fprintf(outfile,"EXCH5       Energy = %18.12lf  H\n",e5);
    fflush(outfile);
  }

  double **B_p_BS = get_BS_ints(1);
  B_p_AB = get_AB_ints(1);

  double *C_p = init_array(calc_info_.nrio);

  C_DGEMV('t',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_BS[0][0]),
    calc_info_.nrio,&(calc_info_.CHFB[0][0]),1,0.0,C_p,1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_AB[0][0]),
    calc_info_.nrio,C_p,1,0.0,&(Y_AB[0][0]),1);

  e3 = 0.0;

  for (int a=0; a<calc_info_.noccA; a++) {
    e3 += 2.0*C_DDOT(calc_info_.noccB,&(calc_info_.S_AB[a][0]),1,&(Y_AB[a][0]),1);
  }

  if (params_.print) {
    fprintf(outfile,"EXCH3       Energy = %18.12lf  H\n",e3);
    fflush(outfile);
  }

  memset(&(Y_AS[0][0]),'\0',sizeof(double)*calc_info_.noccA*calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio,
      &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(Y_AS[0][0]),
      calc_info_.nvirB);
  }

  e4 = C_DDOT(calc_info_.noccA*calc_info_.nvirB,&(X_AS[0][0]),1,&(Y_AS[0][0]),1);

  free_block(B_p_AB);

  if (params_.print) {
    fprintf(outfile,"EXCH4       Energy = %18.12lf  H\n",e4);
    fflush(outfile);
  }

  double **B_p_AA = get_AA_ints(1);

  double **X_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(X_AA[0][0]),calc_info_.noccA);

  double **Y_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,&(B_p_AA[0][0]),
    calc_info_.nrio,C_p,1,0.0,&(Y_AA[0][0]),1);

  free(C_p);

  e9 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,&(X_AA[0][0]),1,
             &(Y_AA[0][0]),1);

  if (params_.print) {
    fprintf(outfile,"EXCH9       Energy = %18.12lf  H\n",e9);
    fflush(outfile);
  }

  double **X_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(X_BB[0][0]),calc_info_.noccB);

  double **X_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(X_BB[0][0]),calc_info_.noccB,&(calc_info_.CHFB[0][0]),
    calc_info_.nvirB,0.0,&(X_BS[0][0]),calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_BS[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,0.0,&(Y_BS[0][0]),1);

  e11 = 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,&(X_BS[0][0]),1,
             &(Y_BS[0][0]),1);

  free_block(X_BS);
  free_block(Y_BS);

  if (params_.print) {
    fprintf(outfile,"EXCH11      Energy = %18.12lf  H\n",e11);
    fflush(outfile);
  }

  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.nrio,calc_info_.noccA,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_AA[0][0]),
    calc_info_.noccA*calc_info_.nrio,0.0,&(C_p_AB[0][0]),
    calc_info_.noccA*calc_info_.nrio);

  memset(&(Y_AS[0][0]),'\0',sizeof(double)*calc_info_.noccA*calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(C_p_AB[b*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(Y_AS[0][0]),
      calc_info_.nvirB);
  }

  e13 = C_DDOT(calc_info_.noccA*calc_info_.nvirB,&(X_AS[0][0]),1,&(Y_AS[0][0]),1);

  free_block(X_AS);
  free_block(Y_AS);
  free_block(C_p_AB);
  free_block(B_p_BS);

  if (params_.print) {
    fprintf(outfile,"EXCH13      Energy = %18.12lf  H\n",e13);
    fflush(outfile);
  }

  B_p_AB = get_AB_ints(2);

  double **X_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nvirB,1.0,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,&(calc_info_.CHFB[0][0]),
    calc_info_.nvirB,0.0,&(X_AB[0][0]),calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_AB[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,0.0,&(Y_AB[0][0]),1);

  e6 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccB,&(X_AB[0][0]),1,
             &(Y_AB[0][0]),1);

  if (params_.print) {
    fprintf(outfile,"EXCH6       Energy = %18.12lf  H\n",e6);
    fflush(outfile);
  }

  memset(&(Y_AB[0][0]),'\0',sizeof(double)*calc_info_.noccA*calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(Y_AB[0][0]),
      calc_info_.noccB);
  }

  e7 = C_DDOT(calc_info_.noccA*calc_info_.noccB,&(X_AB[0][0]),1,&(Y_AB[0][0]),1);

  free_block(Y_AB);
  free_block(B_p_AB);

  if (params_.print) {
    fprintf(outfile,"EXCH7       Energy = %18.12lf  H\n",e7);
    fflush(outfile);
  }

  B_p_BB = get_BB_ints(1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(X_AB[0][0]),
    calc_info_.noccB,0.0,&(X_AA[0][0]),calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,&(B_p_AA[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,0.0,&(Y_AA[0][0]),1);

  e8 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,&(X_AA[0][0]),1,
             &(Y_AA[0][0]),1);

  free_block(X_AA);
  free_block(Y_AA);

  if (params_.print) {
    fprintf(outfile,"EXCH8       Energy = %18.12lf  H\n",e8);
    fflush(outfile);
  }

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(X_AB[0][0]),
    calc_info_.noccB,0.0,&(X_BB[0][0]),calc_info_.noccB);

  double **Y_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_BB[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,0.0,&(Y_BB[0][0]),1);

  e10 = 2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,&(X_BB[0][0]),1,
              &(Y_BB[0][0]),1);

  free_block(X_BB);
  free_block(Y_BB);

  if (params_.print) {
    fprintf(outfile,"EXCH10      Energy = %18.12lf  H\n",e10);
    fflush(outfile);
  }

  C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,1.0,
      &(X_AB[0][0]),calc_info_.noccB,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,1.0,&(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  free_block(B_p_AA);
  double **D_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,calc_info_.noccB,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio,0.0,&(D_p_AB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  e12 = C_DDOT(calc_info_.noccA*calc_info_.noccB*calc_info_.nrio,&(C_p_AB[0][0]),1,
          &(D_p_AB[0][0]),1);

  free_block(X_AB);
  free_block(B_p_BB);
  free_block(C_p_AB);
  free_block(D_p_AB);

  if (params_.print) 
    fprintf(outfile,"EXCH12      Energy = %18.12lf  H\n\n",e12);

  return(-2.0*(e1+e2+e3-e4-e5+e6-e7-e8-e9-e10-e11+e12+e13));
}

double SAPT2::exch_ind20respA_B()
{
  double energy;

  double **uAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);
  double **vAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);
  double **tempAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  double **B_p_AB = get_AB_ints(1);
  double **B_p_RB = get_RB_ints(1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nrio,1.0,&(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,
    &(B_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(tempAR[0][0]),
    calc_info_.nvirA);

  free_block(B_p_AB);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,1.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,1.0,tempAR[0],1,vAR[0],1);

  double **C_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  for (int r=0; r<calc_info_.nvirA; r++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
      &(B_p_RB[r*calc_info_.noccB][0]),calc_info_.nrio,0.0,&(C_p_AR[r][0]),
      calc_info_.nvirA*calc_info_.nrio);
  }

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(C_p_AR[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(tempAR[0][0]),1);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,2.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,2.0,tempAR[0],1,vAR[0],1);

  double **B_p_AA = get_AA_ints(1);

  memset(&(tempAR[0][0]),'\0',sizeof(double)*calc_info_.noccA*
    calc_info_.nvirA);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.nrio,1.0,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(C_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,1.0,&(tempAR[0][0]),
      calc_info_.nvirA);
  }

  free_block(C_p_AR);
  free_block(B_p_RB);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-1.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-1.0,tempAR[0],1,vAR[0],1);

  B_p_AB = get_AB_ints(2);
  double **C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,
    calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,0.0,
      &(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio);
  }

  free_block(B_p_AB);

  double *C_p = init_array(calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DAXPY(calc_info_.nrio,1.0,&(C_p_AA[a*calc_info_.noccA+a][0]),1,C_p,1);
  }

  double **B_p_AR = get_AR_ints(1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,C_p,1,0.0,&(tempAR[0][0]),1);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,4.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,2.0,tempAR[0],1,vAR[0],1);

  memset(&(tempAR[0][0]),'\0',sizeof(double)*calc_info_.noccA*
    calc_info_.nvirA);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.nrio,1.0,
      &(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,1.0,&(tempAR[0][0]),
      calc_info_.nvirA);
  }

  free_block(C_p_AA);
  free_block(B_p_AR);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-1.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-1.0,tempAR[0],1,vAR[0],1);

  B_p_AB = get_AB_ints(1);

  double **X_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(X_AB[0][0]),1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(X_AB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(tempAR[0][0]),calc_info_.nvirA);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,2.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,2.0,tempAR[0],1,vAR[0],1);

  double **B_p_BB = get_BB_ints(1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.noccB*
    calc_info_.nrio,1.0,&(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,
    &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(X_AB[0][0]),
    calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(X_AB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(tempAR[0][0]),calc_info_.nvirA);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-1.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-1.0,tempAR[0],1,vAR[0],1);

  double **X_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(X_BB[0][0]),1);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(X_BB[0][0]),calc_info_.noccB,
    0.0,&(X_AB[0][0]),calc_info_.noccB);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(X_AB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(tempAR[0][0]),calc_info_.nvirA);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-4.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-2.0,tempAR[0],1,vAR[0],1);

  free_block(X_AB);
  free_block(X_BB);

  double **S_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(S_AA[0][0]),calc_info_.noccA);

  double **S_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(S_BB[0][0]),calc_info_.noccB);

  double **S_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,0.0,&(S_AR[0][0]),
    calc_info_.nvirA);

  double **X_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(X_AA[0][0]),1);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccA,1.0,
    &(X_AA[0][0]),calc_info_.noccA,&(S_AR[0][0]),calc_info_.nvirA,0.0,
    &(tempAR[0][0]),calc_info_.nvirA);

  free_block(X_AA);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-4.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-2.0,tempAR[0],1,vAR[0],1);

  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);
  double **C_p_RB = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,0.0,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  free_block(B_p_AA);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,&(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,
    &(C_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nrio,1.0,&(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,
    &(C_p_RB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(tempAR[0][0]),
    calc_info_.nvirA);

  free_block(C_p_AB);
  free_block(C_p_RB);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,2.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,1.0,tempAR[0],1,vAR[0],1);

  B_p_AR = get_AR_ints(1);

  C_DGEMV('t',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,&(S_BB[0][0]),1,0.0,C_p,1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,C_p,1,0.0,&(tempAR[0][0]),1);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-4.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-2.0,tempAR[0],1,vAR[0],1);
 
  free(C_p);

  double **X_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,
    &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(X_AR[0][0]),1);
 
  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccA,1.0,
    &(S_AA[0][0]),calc_info_.noccA,&(X_AR[0][0]),calc_info_.nvirA,0.0,
    &(tempAR[0][0]),calc_info_.nvirA);

  free_block(X_AR);

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-4.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,-2.0,tempAR[0],1,vAR[0],1);

  double **C_p_BA = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  C_p_AA = block_matrix(calc_info_.noccA*calc_info_.noccA,calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
      &(B_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,0.0,
      &(C_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio);
  }

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccB,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_BA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,
    &(C_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio);

  memset(&(tempAR[0][0]),'\0',sizeof(double)*calc_info_.noccA*
    calc_info_.nvirA);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.nrio,1.0,
      &(C_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,1.0,&(tempAR[0][0]),
      calc_info_.nvirA);
  }

  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,2.0,tempAR[0],1,uAR[0],1);
  C_DAXPY(calc_info_.noccA*calc_info_.nvirA,1.0,tempAR[0],1,vAR[0],1);

  energy = -2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,calc_info_.CHFA[0],1,
    vAR[0],1);

  free_block(tempAR);
  free_block(C_p_BA);
  free_block(B_p_BB);
  free_block(B_p_AR);
  free_block(S_AA);
  free_block(S_BB);
  free_block(S_AR);

  write_IJKL(uAR,PSIF_SAPT_AMPS,"AR Exch12 K2f Integrals",calc_info_.noccA,
    calc_info_.nvirA);
  write_IJKL(vAR,PSIF_SAPT_AMPS,"AR Exch-Ind Integrals",calc_info_.noccA,
    calc_info_.nvirA);

  return(energy);
}

double SAPT2::exch_ind20respB_A()
{
  double energy;

  double **uBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);
  double **vBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);
  double **tempBS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  double **B_p_AS = get_AS_ints(1);
  double **B_p_AB = get_AB_ints(2);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,
      &(B_p_AS[a*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(tempBS[0][0]),
      calc_info_.nvirB);
  }

  free_block(B_p_AB);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,1.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,1.0,tempBS[0],1,vBS[0],1);

  double **C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AS[0][0]),calc_info_.nvirB*calc_info_.nrio,0.0,&(C_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(C_p_BS[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(tempBS[0][0]),1);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,2.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,2.0,tempBS[0],1,vBS[0],1);

  free_block(B_p_AS);
  double **B_p_BB = get_BB_ints(1);

  memset(&(tempBS[0][0]),'\0',sizeof(double)*calc_info_.noccB*
    calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(B_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,
      &(C_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(tempBS[0][0]),
      calc_info_.nvirB);
  }

  free_block(C_p_BS);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-1.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-1.0,tempBS[0],1,vBS[0],1);

  B_p_AB = get_AB_ints(1);
  double **C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_BB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(B_p_AB);

  double *C_p = init_array(calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DAXPY(calc_info_.nrio,1.0,&(C_p_BB[b*calc_info_.noccB+b][0]),1,C_p,1);
  }

  double **B_p_BS = get_BS_ints(1);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,C_p,1,0.0,&(tempBS[0][0]),1);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,4.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,2.0,tempBS[0],1,vBS[0],1);

  memset(&(tempBS[0][0]),'\0',sizeof(double)*calc_info_.noccB*
    calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(C_p_BB[b][0]),calc_info_.noccB*calc_info_.nrio,
      &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(tempBS[0][0]),
      calc_info_.nvirB);
  }

  free_block(C_p_BB);
  free_block(B_p_BS);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-1.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-1.0,tempBS[0],1,vBS[0],1);

  B_p_AB = get_AB_ints(2);

  double **X_AB = block_matrix(calc_info_.noccA,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_AB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(X_AB[0][0]),1);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(X_AB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(tempBS[0][0]),calc_info_.nvirB);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,2.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,2.0,tempBS[0],1,vBS[0],1);

  double **B_p_AA = get_AA_ints(1);

  memset(&(X_AB[0][0]),'\0',sizeof(double)*calc_info_.noccA*
    calc_info_.noccB);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(X_AB[0][0]),
      calc_info_.noccB);
  }

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(X_AB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(tempBS[0][0]),calc_info_.nvirB);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-1.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-1.0,tempBS[0],1,vBS[0],1);

  double **X_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(X_AA[0][0]),1);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB,calc_info_.noccA,1.0,
    &(X_AA[0][0]),calc_info_.noccA,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    0.0,&(X_AB[0][0]),calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(X_AB[0][0]),calc_info_.noccB,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(tempBS[0][0]),calc_info_.nvirB);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-4.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-2.0,tempBS[0],1,vBS[0],1);

  free_block(X_AB);
  free_block(X_AA);

  double **S_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(S_AA[0][0]),calc_info_.noccA);

  double **S_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(S_BB[0][0]),calc_info_.noccB);

  double **S_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,0.0,&(S_BS[0][0]),
    calc_info_.nvirB);

  double **X_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,
    &(B_p_BB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(X_BB[0][0]),1);

  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(X_BB[0][0]),calc_info_.noccB,&(S_BS[0][0]),calc_info_.nvirB,0.0,
    &(tempBS[0][0]),calc_info_.nvirB);

  free_block(X_BB);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-4.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-2.0,tempBS[0],1,vBS[0],1);

  double **C_p_BA = block_matrix(calc_info_.noccB*calc_info_.noccA,
    calc_info_.nrio);
  double **C_p_SA = block_matrix(calc_info_.nvirB*calc_info_.noccA,
    calc_info_.nrio);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','N',calc_info_.noccA,calc_info_.nrio,calc_info_.noccB,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
      &(B_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,0.0,
      &(C_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio);
  }

  free_block(B_p_BB);

  C_DGEMM('T','N',calc_info_.nvirB,calc_info_.noccA*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,&(B_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,
    &(C_p_SA[0][0]),calc_info_.noccA*calc_info_.nrio);

  C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA*
    calc_info_.nrio,1.0,&(C_p_BA[0][0]),calc_info_.noccA*calc_info_.nrio,
    &(C_p_SA[0][0]),calc_info_.noccA*calc_info_.nrio,0.0,&(tempBS[0][0]),
    calc_info_.nvirB);

  free_block(C_p_BA);
  free_block(C_p_SA);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,2.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,1.0,tempBS[0],1,vBS[0],1);

  B_p_BS = get_BS_ints(1);

  C_DGEMV('t',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,
    &(B_p_AA[0][0]),calc_info_.nrio,&(S_AA[0][0]),1,0.0,C_p,1);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,C_p,1,0.0,&(tempBS[0][0]),1);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-4.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-2.0,tempBS[0],1,vBS[0],1);
 
  free(C_p);

  double **X_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,
    &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(X_BS[0][0]),1);
 
  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccB,1.0,
    &(S_BB[0][0]),calc_info_.noccB,&(X_BS[0][0]),calc_info_.nvirB,0.0,
    &(tempBS[0][0]),calc_info_.nvirB);

  free_block(X_BS);

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-4.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,-2.0,tempBS[0],1,vBS[0],1);

  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nrio);
  C_p_BB = block_matrix(calc_info_.noccB*calc_info_.noccB,calc_info_.nrio);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,0.0,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio);
  }

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB*calc_info_.nrio,
    calc_info_.noccA,1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,
    &(C_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio);

  memset(&(tempBS[0][0]),'\0',sizeof(double)*calc_info_.noccB*
    calc_info_.nvirB);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(C_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,
      &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(tempBS[0][0]),
      calc_info_.nvirB);
  }

  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,2.0,tempBS[0],1,uBS[0],1);
  C_DAXPY(calc_info_.noccB*calc_info_.nvirB,1.0,tempBS[0],1,vBS[0],1);

  energy = -2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,calc_info_.CHFB[0],1,
    vBS[0],1);

  free_block(tempBS);
  free_block(C_p_AB);
  free_block(B_p_AA);
  free_block(B_p_BS);
  free_block(S_AA);
  free_block(S_BB);
  free_block(S_BS);

  write_IJKL(uBS,PSIF_SAPT_AMPS,"BS Exch12 K2f Integrals",calc_info_.noccB,
    calc_info_.nvirB);
  write_IJKL(vBS,PSIF_SAPT_AMPS,"BS Exch-Ind Integrals",calc_info_.noccB,
    calc_info_.nvirB);

  return(energy);
}

}}
