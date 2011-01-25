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

void SAPT2p3::ind30()
{
  double e1,e2;

  if (params_.print)
    fprintf(outfile,"Begining Ind30 Calculation\n\n");

  ind30_amps("Ind30 AR Amplitudes","T ARBS Amplitudes",PSIF_SAPT_AA_DF_INTS,
    "AR RI Integrals",PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",calc_info_.sA,
    calc_info_.sB,calc_info_.WBAA,calc_info_.WBRR,calc_info_.WABS,
    calc_info_.evalsA,calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,
    calc_info_.nvirB);

  double **X_AR = read_IJKL(PSIF_SAPT_AMPS,"Ind30 AR Amplitudes",
    calc_info_.noccA,calc_info_.nvirA);

  e1 = 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,X_AR[0],1,
    calc_info_.WBAR[0],1);

  free_block(X_AR);

  ind30_amps("Ind30 BS Amplitudes","T BSAR Amplitudes",PSIF_SAPT_BB_DF_INTS,
    "BS RI Integrals",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",calc_info_.sB,
    calc_info_.sA,calc_info_.WABB,calc_info_.WASS,calc_info_.WBAR,
    calc_info_.evalsB,calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,
    calc_info_.nvirA);

  double **X_BS = read_IJKL(PSIF_SAPT_AMPS,"Ind30 BS Amplitudes",
    calc_info_.noccB,calc_info_.nvirB);

  e2 = 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,X_BS[0],1,
    calc_info_.WABS[0],1);

  free_block(X_BS);

  results_.ind30 = e1+e2;
  if (params_.print) {
    fprintf(outfile,"Ind30              = %18.12lf  H\n\n",results_.ind30);
    fflush(outfile);
  }
}

void SAPT2p3::ind30_amps(char *ind_out, char *T_amps, int AAfile,
  char *AR_label, int BBfile, char *BS_label, double **sAR, double **sBS,
  double **WBAA, double **WBRR, double **WABS, double *evalsA, int noccA,
  int nvirA, int noccB, int nvirB)
{
  double **uAR = block_matrix(noccA,nvirA);

  C_DGEMM('N','N',noccA,nvirA,nvirA,1.0,&(sAR[0][0]),nvirA,&(WBRR[0][0]),
    nvirA,0.0,&(uAR[0][0]),nvirA);

  C_DGEMM('N','N',noccA,nvirA,noccA,-1.0,&(WBAA[0][0]),noccA,&(sAR[0][0]),
    nvirA,1.0,&(uAR[0][0]),nvirA);

  double **B_p_BS = get_DF_ints(BBfile,BS_label,noccB*nvirB);
  double *C_p = init_array(calc_info_.nrio);

  C_DGEMV('t',nvirB*noccB,calc_info_.nrio,1.0,&(B_p_BS[0][0]),calc_info_.nrio,
    &(sBS[0][0]),1,0.0,C_p,1);

  free_block(B_p_BS);

  double **B_p_AR = get_DF_ints(AAfile,AR_label,noccA*nvirA);

  C_DGEMV('n',noccA*nvirA,calc_info_.nrio,2.0,&(B_p_AR[0][0]),calc_info_.nrio,
    C_p,1,1.0,&(uAR[0][0]),1);

  free(C_p);
  free_block(B_p_AR);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,T_amps,noccA*nvirA,nvirB*noccB);

  C_DGEMV('n',noccA*nvirA,nvirB*noccB,2.0,&(tARBS[0][0]),nvirB*noccB,
    &(WABS[0][0]),1,1.0,&(uAR[0][0]),1);

  free_block(tARBS);

  for(int a=0; a<noccA; a++) {
  for(int r=0; r<nvirA; r++) {
      double denom = evalsA[a] - evalsA[r+noccA];
      uAR[a][r] = uAR[a][r]/denom;
  }}

  write_IJKL(uAR,PSIF_SAPT_AMPS,ind_out,noccA,nvirA);
}

}}

