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

  disp30_amps();

  double **X_ARBS = read_IJKL(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.nvirB*calc_info_.noccB);

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
    calc_info_.nvirB*calc_info_.noccB,&(ARBS[0][0]),1,&(X_ARBS[0][0]),1);

  free_block(X_ARBS);
  free_block(ARBS);

  if (params_.print) {
    fprintf(outfile,"Disp30             = %18.12lf  H\n\n",results_.disp30);
    fflush(outfile);
  }
}

void SAPT2p3::disp30_amps()
{
  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
  double **tABRS = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nvirA*calc_info_.nvirB);

    for (int a=0, ar=0; a < calc_info_.noccA; a++) {
    for (int r=0; r < calc_info_.nvirA; r++, ar++) {
      for (int b=0, bs=0; b < calc_info_.noccB; b++) {
      for (int s=0; s < calc_info_.nvirB; s++, bs++) {
        int ab = a*calc_info_.noccB + b;
        int rs = r*calc_info_.nvirB + s;
        tABRS[ab][rs] = tARBS[ar][bs];
      }}
    }}

  free_block(tARBS);

  double **t2ABRS = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.nvirA*calc_info_.nvirB);

  double **B_p_RR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",
    calc_info_.nvirA*calc_info_.nvirA);
  double **B_p_SS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",
    calc_info_.nvirB*calc_info_.nvirB);

  double **X_RS = block_matrix(calc_info_.nvirA,calc_info_.nvirB*
    calc_info_.nvirB);

  for (int r=0; r < calc_info_.nvirA; r++) {
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.nvirB*calc_info_.nvirB,
      calc_info_.nrio,1.0,&(B_p_RR[r*calc_info_.nvirA][0]),calc_info_.nrio,
      &(B_p_SS[0][0]),calc_info_.nrio,0.0,&(X_RS[0][0]),calc_info_.nvirB*
      calc_info_.nvirB);
    C_DGEMM('N','T',calc_info_.noccA*calc_info_.noccB,
      calc_info_.nvirA*calc_info_.nvirB,calc_info_.nvirB,1.0,
      &(tABRS[0][r*calc_info_.nvirB]),calc_info_.nvirA*calc_info_.nvirB,
      &(X_RS[0][0]),calc_info_.nvirB,1.0,&(t2ABRS[0][0]),calc_info_.nvirA*
      calc_info_.nvirB);
  }

  free_block(B_p_RR);
  free_block(B_p_SS);
  free_block(X_RS);

  double **B_p_AA = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    calc_info_.noccA*calc_info_.noccA);
  double **B_p_BB = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    calc_info_.noccB*calc_info_.noccB);

  double **ABAB = block_matrix(calc_info_.noccA*calc_info_.noccB,
    calc_info_.noccA*calc_info_.noccB);

  for (int a=0, ab=0; a < calc_info_.noccA; a++) {
    for (int b=0; b < calc_info_.noccB; b++,ab++) {
      C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccB,calc_info_.nrio,1.0,
        &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
        &(B_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,0.0,&(ABAB[ab][0]),
        calc_info_.noccB);
  }}

  free_block(B_p_AA);
  free_block(B_p_BB);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.noccB,calc_info_.nvirA*
    calc_info_.nvirB,calc_info_.noccA*calc_info_.noccB,1.0,&(ABAB[0][0]),
    calc_info_.noccA*calc_info_.noccB,&(tABRS[0][0]),calc_info_.nvirA*
    calc_info_.nvirB,1.0,&(t2ABRS[0][0]),calc_info_.nvirA*calc_info_.nvirB);

  free_block(ABAB);

  double **tBRAS = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirB);

    for (int a=0, ab=0; a < calc_info_.noccA; a++) {
    for (int b=0; b < calc_info_.noccB; b++, ab++) {
      for (int r=0, rs=0; r < calc_info_.nvirA; r++) {
      for (int s=0; s < calc_info_.nvirB; s++, rs++) {
        int br = b*calc_info_.nvirA + r;
        int as = a*calc_info_.nvirB + s;
        tBRAS[br][as] = tABRS[ab][rs];
      }}
    }}

  free_block(tABRS);

  double **t2BRAS = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirB);

    for (int a=0, ab=0; a < calc_info_.noccA; a++) {
    for (int b=0; b < calc_info_.noccB; b++, ab++) {
      for (int r=0, rs=0; r < calc_info_.nvirA; r++) {
      for (int s=0; s < calc_info_.nvirB; s++, rs++) {
        int br = b*calc_info_.nvirA + r;
        int as = a*calc_info_.nvirB + s;
        t2BRAS[br][as] = t2ABRS[ab][rs];
      }}
    }}

  free_block(t2ABRS);

  B_p_BB = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    calc_info_.noccB*calc_info_.noccB);
  B_p_RR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",
    calc_info_.nvirA*calc_info_.nvirA);

  double **BRBR = block_matrix(calc_info_.noccB*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirA);

  for (int b=0, br=0; b < calc_info_.noccB; b++) {
    for (int r=0; r < calc_info_.nvirA; r++, br++) {
      C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirA,calc_info_.nrio,1.0,
        &(B_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,
        &(B_p_RR[r*calc_info_.nvirA][0]),calc_info_.nrio,0.0,&(BRBR[br][0]),
        calc_info_.nvirA);
  }}

  free_block(B_p_BB);
  free_block(B_p_RR);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirA,-1.0,&(BRBR[0][0]),
    calc_info_.noccB*calc_info_.nvirA,&(tBRAS[0][0]),calc_info_.noccA*
    calc_info_.nvirB,1.0,&(t2BRAS[0][0]),calc_info_.noccA*calc_info_.nvirB);

  free_block(BRBR);

  B_p_AA = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    calc_info_.noccA*calc_info_.noccA);
  B_p_SS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",
    calc_info_.nvirB*calc_info_.nvirB);

  double **ASAS = block_matrix(calc_info_.noccA*calc_info_.nvirB,
    calc_info_.noccA*calc_info_.nvirB);

  for (int a=0, as=0; a < calc_info_.noccA; a++) {
    for (int s=0; s < calc_info_.nvirB; s++, as++) {
      C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirB,calc_info_.nrio,1.0,
        &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
        &(B_p_SS[s*calc_info_.nvirB][0]),calc_info_.nrio,0.0,&(ASAS[as][0]),
        calc_info_.nvirB);
  }}

  free_block(B_p_AA);
  free_block(B_p_SS);

  C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirB,calc_info_.noccA*calc_info_.nvirB,-1.0,&(tBRAS[0][0]),
    calc_info_.noccA*calc_info_.nvirB,&(ASAS[0][0]),calc_info_.noccA*
    calc_info_.nvirB,1.0,&(t2BRAS[0][0]),calc_info_.noccA*calc_info_.nvirB);

  free_block(ASAS);
  free_block(tBRAS);

  tARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB);

  for (int a=0,ar=0; a < calc_info_.noccA; a++) {
    for (int r=0; r < calc_info_.nvirA; r++,ar++) {
      for (int b=0,bs=0; b < calc_info_.noccB; b++) {
        for (int s=0; s < calc_info_.nvirB; s++,bs++) {
          int br = b*calc_info_.nvirA + r;
          int as = a*calc_info_.nvirB + s;
          double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
                  calc_info_.evalsA[r+calc_info_.noccA]-
                  calc_info_.evalsB[s+calc_info_.noccB];
          tARBS[ar][bs] = t2BRAS[br][as]/denom;
  }}}}

  free_block(t2BRAS);

  write_IJKL(tARBS,PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
}

}}

