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

void SAPT2p3::disp30()
{
  if (params_.print)
    fprintf(outfile,"Begining Disp30 Calculation\n\n");

  double **tARBS;

  if (params_.nat_orbs && 0) {
    psio_->open(PSIF_SAPT_TEMP,0);
    natural_orbitalify_disp30();

    double **wARBS = disp30_amps(PSIF_SAPT_TEMP,"T ARBS Amplitudes",PSIF_SAPT_TEMP,"AA RI Integrals","RR RI Integrals",PSIF_SAPT_TEMP,"BB RI Integrals","SS RI Integrals",
      no_info_.evalsA,no_info_.evalsB,calc_info_.noccA,no_info_.nvirA,
      params_.foccA,calc_info_.noccB,no_info_.nvirB,params_.foccB);

    double **xRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);
    for (int r=0; r < calc_info_.nvirA; r++)
      xRR[r][r] = 1.0;
    double **S_RR = block_matrix(calc_info_.nvirA,no_info_.nvirA);
    C_DGEMM('N','N',calc_info_.nvirA,no_info_.nvirA,calc_info_.nvirA,1.0,
      &(xRR[0][0]),calc_info_.nvirA,&(no_info_.CA[calc_info_.noccA]
      [calc_info_.noccA]),calc_info_.noccA+no_info_.nvirA,0.0,&(S_RR[0][0]),
      no_info_.nvirA);
    free_block(xRR);

    double **xSS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);
    for (int s=0; s < calc_info_.nvirB; s++)
      xSS[s][s] = 1.0;
    double **S_SS = block_matrix(calc_info_.nvirB,no_info_.nvirB);
    C_DGEMM('N','N',calc_info_.nvirB,no_info_.nvirB,calc_info_.nvirB,1.0,
      &(xSS[0][0]),calc_info_.nvirB,&(no_info_.CB[calc_info_.noccB]
      [calc_info_.noccB]),calc_info_.noccB+no_info_.nvirB,0.0,&(S_SS[0][0]),
      no_info_.nvirB);
    free_block(xSS);

    double **xARBS = block_matrix((calc_info_.noccA-params_.foccA)*
      calc_info_.nvirA,(calc_info_.noccB-params_.foccB)*no_info_.nvirB);

    for (int a=0; a < calc_info_.noccA-params_.foccA; a++) {
      C_DGEMM('N','N',calc_info_.nvirA,(calc_info_.noccB-params_.foccB)*
        no_info_.nvirB,no_info_.nvirA,1.0,&(S_RR[0][0]),no_info_.nvirA,
        &(wARBS[a*no_info_.nvirA][0]),(calc_info_.noccB-params_.foccB)*
        no_info_.nvirB,0.0,&(xARBS[a*calc_info_.nvirA][0]),
        (calc_info_.noccB-params_.foccB)*no_info_.nvirB);
    }

    free_block(wARBS);

    double **yARBS = block_matrix((calc_info_.noccA-params_.foccA)*
      calc_info_.nvirA,(calc_info_.noccB-params_.foccB)*calc_info_.nvirB);

    C_DGEMM('N','T',(calc_info_.noccA-params_.foccA)*calc_info_.nvirA*
      (calc_info_.noccB-params_.foccB),calc_info_.nvirB,no_info_.nvirB,1.0,
      &(xARBS[0][0]),no_info_.nvirB,&(S_SS[0][0]),no_info_.nvirB,0.0,
      &(yARBS[0][0]),calc_info_.nvirB);

    free_block(xARBS);
    free_block(S_RR);
    free_block(S_SS);

    tARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
      calc_info_.nvirB);
    for (int a=params_.foccA, ar=0; a < calc_info_.noccA; a++) {
      for (int r=0; r < calc_info_.nvirA; r++,ar++) {
        int ar_ = a*calc_info_.nvirA+r;
        C_DCOPY((calc_info_.noccB-params_.foccB)*calc_info_.nvirB,
          &(yARBS[ar][0]),1,&(tARBS[ar_][params_.foccB*calc_info_.nvirB]),1);
    }}
    free_block(yARBS);
    psio_->close(PSIF_SAPT_TEMP,0);
  }
  else if ((params_.foccA || params_.foccB) && 0) {
    psio_->open(PSIF_SAPT_TEMP,0);
    frzn_disp30_prep();
    double **xARBS = disp30_amps(PSIF_SAPT_TEMP,"T ARBS Amplitudes",PSIF_SAPT_TEMP,"AA RI Integrals","RR RI Integrals",PSIF_SAPT_TEMP,"BB RI Integrals","SS RI Integrals",
      calc_info_.evalsA,calc_info_.evalsB,calc_info_.noccA,calc_info_.nvirA,
      params_.foccA,calc_info_.noccB,calc_info_.nvirB,params_.foccB);
    tARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
      calc_info_.nvirB);
    for (int a=params_.foccA, ar=0; a < calc_info_.noccA; a++) {
      for (int r=0; r < calc_info_.nvirA; r++,ar++) {
        int ar_ = a*calc_info_.nvirA+r;
        C_DCOPY((calc_info_.noccB-params_.foccB)*calc_info_.nvirB,
          &(xARBS[ar][0]),1,&(tARBS[ar_][params_.foccB*calc_info_.nvirB]),1);
    }}
    free_block(xARBS);
    psio_->close(PSIF_SAPT_TEMP,0);
  }
  else {
    tARBS = disp30_amps(PSIF_SAPT_AMPS,"T ARBS Amplitudes",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","RR RI Integrals",
      PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","SS RI Integrals",calc_info_.evalsA,calc_info_.evalsB,
      calc_info_.noccA,calc_info_.nvirA,0,calc_info_.noccB,
      calc_info_.nvirB,0);
  }

  double **B_p_AR = get_AR_ints(0);
  double **B_p_BS = get_BS_ints(0);
  double **ARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nri,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(ARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(B_p_AR);
  free_block(B_p_BS);

  results_.disp30 = 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*
    calc_info_.nvirB*calc_info_.noccB,&(ARBS[0][0]),1,&(tARBS[0][0]),1);

  free_block(ARBS);

  write_IJKL(tARBS,PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);

  if (params_.print) {
    fprintf(outfile,"Disp30             = %18.12lf  H\n\n",results_.disp30);
    fflush(outfile);
  }
}

void SAPT2p3::frzn_disp30_prep()
{
  psio_address next_psio;
  double **B_p_AA = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    calc_info_.noccA*calc_info_.noccA);

  next_psio = PSIO_ZERO;
  for(int a=params_.foccA; a<calc_info_.noccA; a++) {
    int aa = a*calc_info_.noccA + params_.foccA;
    psio_->write(PSIF_SAPT_TEMP,"AA RI Integrals",(char *) &(B_p_AA[aa][0]),
      (calc_info_.noccA-params_.foccA)*calc_info_.nrio*(ULI) sizeof(double),
      next_psio,&next_psio);
  }

  free_block(B_p_AA);

  double **B_p_BB = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    calc_info_.noccB*calc_info_.noccB);

  next_psio = PSIO_ZERO;
  for(int b=params_.foccB; b<calc_info_.noccB; b++) {
    int bb = b*calc_info_.noccB + params_.foccB;
    psio_->write(PSIF_SAPT_TEMP,"BB RI Integrals",(char *) &(B_p_BB[bb][0]),
      (calc_info_.noccB-params_.foccB)*calc_info_.nrio*(ULI) sizeof(double),
      next_psio,&next_psio);
  }

  free_block(B_p_BB);

  double **B_p_RR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",
    calc_info_.nvirA*calc_info_.nvirA);

  psio_->write_entry(PSIF_SAPT_TEMP,"RR RI Integrals",(char *) &(B_p_RR[0][0]),
      calc_info_.nvirA*calc_info_.nvirA*calc_info_.nrio*(ULI) sizeof(double));

  free_block(B_p_RR);

  double **B_p_SS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",
    calc_info_.nvirB*calc_info_.nvirB);
  
  psio_->write_entry(PSIF_SAPT_TEMP,"SS RI Integrals",(char *) &(B_p_SS[0][0]),
      calc_info_.nvirB*calc_info_.nvirB*calc_info_.nrio*(ULI) sizeof(double));

  free_block(B_p_SS);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  next_psio = PSIO_ZERO;
  for (int a=params_.foccA; a < calc_info_.noccA; a++) {
    for (int r=0; r < calc_info_.nvirA; r++) {
      int ar = a*calc_info_.nvirA+r;
      for (int b=params_.foccB; b < calc_info_.noccB; b++) {
      psio_->write(PSIF_SAPT_TEMP,"T ARBS Amplitudes",(char *) 
        &(tARBS[ar][b*calc_info_.nvirB]),calc_info_.nvirB*(ULI) sizeof(double),
        next_psio,&next_psio);
  }}}

  free_block(tARBS);
}

void SAPT2p3::natural_orbitalify_disp30()
{
  psio_address next_psio;
  double **B_p_AA = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    calc_info_.noccA*calc_info_.noccA);
    
  next_psio = PSIO_ZERO;
  for(int a=params_.foccA; a<calc_info_.noccA; a++) {
    int aa = a*calc_info_.noccA + params_.foccA;
    psio_->write(PSIF_SAPT_TEMP,"AA RI Integrals",(char *) &(B_p_AA[aa][0]),
      (calc_info_.noccA-params_.foccA)*calc_info_.nrio*(ULI) sizeof(double),
      next_psio,&next_psio);
  }

  free_block(B_p_AA);

  double **B_p_BB = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    calc_info_.noccB*calc_info_.noccB);

  next_psio = PSIO_ZERO;
  for(int b=params_.foccB; b<calc_info_.noccB; b++) {
    int bb = b*calc_info_.noccB + params_.foccB;
    psio_->write(PSIF_SAPT_TEMP,"BB RI Integrals",(char *) &(B_p_BB[bb][0]),
      (calc_info_.noccB-params_.foccB)*calc_info_.nrio*(ULI) sizeof(double),
      next_psio,&next_psio);
  }

  free_block(B_p_BB);

  double **B_p_RR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",
    calc_info_.nvirA*calc_info_.nvirA);

  double **C_p_RR = block_matrix(no_info_.nvirA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('T','N',no_info_.nvirA,calc_info_.nvirA*calc_info_.nrio,
    calc_info_.nvirA,1.0,&(no_info_.CA[calc_info_.noccA][calc_info_.noccA]),
    calc_info_.noccA+no_info_.nvirA,B_p_RR[0],calc_info_.nvirA*calc_info_.nrio,
    0.0,C_p_RR[0],calc_info_.nvirA*calc_info_.nrio);

  free_block(B_p_RR);
  double **D_p_RR = block_matrix(no_info_.nvirA*no_info_.nvirA,
    calc_info_.nrio);

  for(int r=0; r<no_info_.nvirA; r++) {
    C_DGEMM('T','N',no_info_.nvirA,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(no_info_.CA[calc_info_.noccA][calc_info_.noccA]),calc_info_.noccA+
      no_info_.nvirA,C_p_RR[r*calc_info_.nvirA],calc_info_.nrio,0.0,
      D_p_RR[r*no_info_.nvirA],calc_info_.nrio);
  }
    
  psio_->write_entry(PSIF_SAPT_TEMP,"RR RI Integrals",(char *) &(D_p_RR[0][0]),
      no_info_.nvirA*no_info_.nvirA*calc_info_.nrio*(ULI) sizeof(double));

  free_block(C_p_RR);
  free_block(D_p_RR);

  double **B_p_SS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",    
    calc_info_.nvirB*calc_info_.nvirB);
  
  double **C_p_SS = block_matrix(no_info_.nvirB*calc_info_.nvirB,
    calc_info_.nrio);
  
  C_DGEMM('T','N',no_info_.nvirB,calc_info_.nvirB*calc_info_.nrio,
    calc_info_.nvirB,1.0,&(no_info_.CB[calc_info_.noccB][calc_info_.noccB]),
    calc_info_.noccB+no_info_.nvirB,B_p_SS[0],calc_info_.nvirB*calc_info_.nrio,
    0.0,C_p_SS[0],calc_info_.nvirB*calc_info_.nrio);

  free_block(B_p_SS);
  double **D_p_SS = block_matrix(no_info_.nvirB*no_info_.nvirB,
    calc_info_.nrio);

  for(int s=0; s<no_info_.nvirB; s++) {
    C_DGEMM('T','N',no_info_.nvirB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(no_info_.CB[calc_info_.noccB][calc_info_.noccB]),calc_info_.noccB+
      no_info_.nvirB,C_p_SS[s*calc_info_.nvirB],calc_info_.nrio,0.0,
      D_p_SS[s*no_info_.nvirB],calc_info_.nrio);
  }

  psio_->write_entry(PSIF_SAPT_TEMP,"SS RI Integrals",(char *) &(D_p_SS[0][0]),
      no_info_.nvirB*no_info_.nvirB*calc_info_.nrio*(ULI) sizeof(double));

  free_block(C_p_SS);
  free_block(D_p_SS);

  double **B_p_AR = block_matrix((calc_info_.noccA-params_.foccA)*
    calc_info_.nvirA,calc_info_.nrio);

  next_psio = psio_get_address(PSIO_ZERO,params_.foccA*calc_info_.nvirA*
    calc_info_.nrio*(ULI) sizeof(double));
  psio_->read(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *) &(B_p_AR[0][0]),
    (calc_info_.noccA-params_.foccA)*calc_info_.nvirA*calc_info_.nrio*
    (ULI) sizeof(double),next_psio,&next_psio);

  double **C_p_AR = block_matrix((calc_info_.noccA-params_.foccA)*
    no_info_.nvirA,calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA-params_.foccA; a++) {
    C_DGEMM('T','N',no_info_.nvirA,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(no_info_.CA[calc_info_.noccA][calc_info_.noccA]),calc_info_.noccA+
      no_info_.nvirA,B_p_AR[a*calc_info_.nvirA],calc_info_.nrio,0.0,
      C_p_AR[a*no_info_.nvirA],calc_info_.nrio);
  }

  free_block(B_p_AR);

  double **B_p_BS = block_matrix((calc_info_.noccB-params_.foccB)*
    calc_info_.nvirB,calc_info_.nrio);

  next_psio = psio_get_address(PSIO_ZERO,params_.foccB*calc_info_.nvirB*
    calc_info_.nrio*(ULI) sizeof(double));
  psio_->read(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *) &(B_p_BS[0][0]),
    (calc_info_.noccB-params_.foccB)*calc_info_.nvirB*calc_info_.nrio*
    (ULI) sizeof(double),next_psio,&next_psio);

  double **C_p_BS = block_matrix((calc_info_.noccB-params_.foccB)*
    no_info_.nvirB,calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB-params_.foccB; b++) {
    C_DGEMM('T','N',no_info_.nvirB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(no_info_.CB[calc_info_.noccB][calc_info_.noccB]),calc_info_.noccB+
      no_info_.nvirB,B_p_BS[b*calc_info_.nvirB],calc_info_.nrio,0.0,
      C_p_BS[b*no_info_.nvirB],calc_info_.nrio);
  }

  free_block(B_p_BS);
 
  double **tARBS = block_matrix((calc_info_.noccA-params_.foccA)*
    no_info_.nvirA,(calc_info_.noccB-params_.foccB)*no_info_.nvirB);

  C_DGEMM('N','T',(calc_info_.noccA-params_.foccA)*no_info_.nvirA,
    (calc_info_.noccB-params_.foccB)*no_info_.nvirB,calc_info_.nrio,
    1.0,&(C_p_AR[0][0]),calc_info_.nrio,&(C_p_BS[0][0]),calc_info_.nrio,0.0,
    &(tARBS[0][0]),(calc_info_.noccB-params_.foccB)*no_info_.nvirB);

  free_block(C_p_AR);
  free_block(C_p_BS);

  for (int a=0, ar=0; a < calc_info_.noccA-params_.foccA; a++) {
  for (int r=0; r < no_info_.nvirA; r++, ar++) {
    for (int b=0, bs=0; b < calc_info_.noccB-params_.foccB; b++) {
    for (int s=0; s < no_info_.nvirB; s++, bs++) {
      double denom = no_info_.evalsA[a+params_.foccA]+
        no_info_.evalsB[b+params_.foccB]-
        no_info_.evalsA[r+calc_info_.noccA]-
        no_info_.evalsB[s+calc_info_.noccB];
      tARBS[ar][bs] /= denom;
    }}
  }}

  psio_->write_entry(PSIF_SAPT_TEMP,"T ARBS Amplitudes",(char *)
    &(tARBS[0][0]),(calc_info_.noccA-params_.foccA)*no_info_.nvirA*
    (calc_info_.noccB-params_.foccB)*no_info_.nvirB*(ULI) sizeof(double));

  free_block(tARBS);
}

}}

