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
#include "sapt2b.h"
#include "sapt2p3.h"

namespace psi { namespace sapt {

void SAPT2p3::disp30()
{
  if (params_.print)
    fprintf(outfile,"Begining Disp30 Calculation\n\n");

  double e1 = disp30_1(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",PSIF_SAPT_BB_DF_INTS,
    "SS RI Integrals",calc_info_.noccA,calc_info_.nvirA,
    calc_info_.noccB,calc_info_.nvirB);
  double e2 = disp30_2(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","RR RI Integrals",
    PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","SS RI Integrals",
    calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,calc_info_.nvirB);

  results_.disp30 = e1+e2;

  if (params_.print) {
    fprintf(outfile,"Disp30             = %18.12lf  H\n\n",results_.disp30);
    fflush(outfile);
  }
}

double SAPT2p3::disp30_1(int ampfile, const char *amplabel, int AAintfile, 
  const char *RRlabel, int BBintfile, const char *SSlabel, int noccA, 
  int nvirA, int noccB, int nvirB)
{
  double **tARBS = read_IJKL(ampfile,amplabel,noccA*nvirA,noccB*nvirB);
  double **tRSAB = block_matrix(nvirA*nvirB,noccA*noccB);

    for (int a=0, ar=0; a < noccA; a++) {
    for (int r=0; r < nvirA; r++, ar++) {
      for (int b=0, bs=0; b < noccB; b++) {
      for (int s=0; s < nvirB; s++, bs++) {
        int ab = a*noccB + b;
        int rs = r*nvirB + s;
        int sr = s*nvirA + r;
        tRSAB[rs][ab] = tARBS[ar][bs];
      }}
    }}

  free_block(tARBS);

  double energy = 0.0;

  psio_address next_DF_RR = PSIO_ZERO;
  psio_address next_DF_SS = PSIO_ZERO;

  double **B_p_RR = block_matrix(nvirA*(nvirA+1)/2,ribasis_->nbf()+3);
  double **B_p_SS = block_matrix(nvirB*(nvirB+1)/2,ribasis_->nbf()+3);

  for (int r1=0,r1r2=0; r1 < nvirA; r1++) {
  for (int r2=0; r2 <= r1; r2++,r1r2++) {
    next_DF_RR = psio_get_address(PSIO_ZERO,(r1*nvirA+r2)*(ribasis_->nbf()+3)*
      (ULI) sizeof(double));
    psio_->read(AAintfile,RRlabel,(char *) &(B_p_RR[r1r2][0]),
      (ribasis_->nbf()+3)*(ULI) sizeof(double),next_DF_RR,&next_DF_RR);
    if (r1 != r2) C_DSCAL(ribasis_->nbf()+3,2.0,B_p_RR[r1r2],1);
  }}

  for (int s1=0,s1s2=0; s1 < nvirB; s1++) {
  for (int s2=0; s2 <= s1; s2++,s1s2++) {
    next_DF_SS = psio_get_address(PSIO_ZERO,(s1*nvirB+s2)*(ribasis_->nbf()+3)*
      (ULI) sizeof(double));
    psio_->read(BBintfile,SSlabel,(char *) &(B_p_SS[s1s2][0]),
      (ribasis_->nbf()+3)*(ULI) sizeof(double),next_DF_SS,&next_DF_SS);
    if (s1 != s2) C_DSCAL(ribasis_->nbf()+3,2.0,B_p_SS[s1s2],1);
  }}

  double **xRS = block_matrix(nvirA,nvirB*nvirB);
  double **yRS = block_matrix(nvirA,nvirB*(nvirB+1)/2);
  double *zSS = init_array(nvirB*(nvirB+1)/2);

  for (int r1=0; r1 < nvirA; r1++) {
    C_DGEMM('N','T',(r1+1)*nvirB,nvirB,noccA*noccB,1.0,tRSAB[0],noccA*noccB,
      tRSAB[r1*nvirB],noccA*noccB,0.0,xRS[0],nvirB); 
    C_DGEMM('N','T',(r1+1),nvirB*(nvirB+1)/2,ribasis_->nbf()+3,1.0,
      B_p_RR[calc_info_.ioff[r1]],ribasis_->nbf()+3,B_p_SS[0],
      ribasis_->nbf()+3,0.0,yRS[0],nvirB*(nvirB+1)/2);
    for (int r2=0; r2 <= r1; r2++) {
      for (int s1=0,s1s2=0; s1 < nvirB; s1++) {
      for (int s2=0; s2 <= s1; s2++,s1s2++) {
        zSS[s1s2] = xRS[r2][s1*nvirB+s2];
        zSS[s1s2] += xRS[r2][s2*nvirB+s1];
      }}
      energy += 2.0*C_DDOT(nvirB*(nvirB+1)/2,zSS,1,yRS[r2],1);
  }}

  free_block(B_p_RR);
  free_block(B_p_SS);
  free_block(xRS);
  free_block(yRS);
  free(zSS);
  free_block(tRSAB);

  return(energy);
}

double SAPT2p3::disp30_2(int ampfile, const char *amplabel, int AAintfile,
  const char *AAlabel, const char *RRlabel, int BBintfile, const char *BBlabel,
  const char *SSlabel, int noccA, int nvirA, int noccB, int nvirB)
{
  double **tARBS = read_IJKL(ampfile,amplabel,noccA*nvirA,noccB*nvirB);
  double **tABRS = block_matrix(noccA*noccB,nvirA*nvirB);

    for (int a=0, ar=0; a < noccA; a++) {
    for (int r=0; r < nvirA; r++, ar++) {
      for (int b=0, bs=0; b < noccB; b++) {
      for (int s=0; s < nvirB; s++, bs++) {
        int ab = a*noccB + b;
        int rs = r*nvirB + s;
        tABRS[ab][rs] = tARBS[ar][bs];
      }}
    }}

  free_block(tARBS);

  double **t2ABRS = block_matrix(noccA*noccB,nvirA*nvirB);

  double **B_p_AA = get_DF_ints(AAintfile,AAlabel,noccA*noccA);
  double **B_p_BB = get_DF_ints(BBintfile,BBlabel,noccB*noccB);

  double **ABAB = block_matrix(noccA*noccB,noccA*noccB);

  for (int a=0, ab=0; a < noccA; a++) {
    for (int b=0; b < noccB; b++,ab++) {
      C_DGEMM('N','T',noccA,noccB,ribasis_->nbf()+3,1.0,&(B_p_AA[a*noccA][0]),
        ribasis_->nbf()+3,&(B_p_BB[b*noccB][0]),ribasis_->nbf()+3,0.0,
        &(ABAB[ab][0]),noccB);
  }}

  free_block(B_p_AA);
  free_block(B_p_BB);

  C_DGEMM('N','N',noccA*noccB,nvirA*nvirB,noccA*noccB,1.0,&(ABAB[0][0]),
    noccA*noccB,&(tABRS[0][0]),nvirA*nvirB,1.0,&(t2ABRS[0][0]),nvirA*nvirB);

  free_block(ABAB);

  double **tBRAS = block_matrix(noccB*nvirA,noccA*nvirB);

    for (int a=0, ab=0; a < noccA; a++) {
    for (int b=0; b < noccB; b++, ab++) {
      for (int r=0, rs=0; r < nvirA; r++) {
      for (int s=0; s < nvirB; s++, rs++) {
        int br = b*nvirA + r;
        int as = a*nvirB + s;
        tBRAS[br][as] = tABRS[ab][rs];
      }}
    }}

  free_block(tABRS);

  double **t2BRAS = block_matrix(noccB*nvirA,noccA*nvirB);

    for (int a=0, ab=0; a < noccA; a++) {
    for (int b=0; b < noccB; b++, ab++) {
      for (int r=0, rs=0; r < nvirA; r++) {
      for (int s=0; s < nvirB; s++, rs++) {
        int br = b*nvirA + r;
        int as = a*nvirB + s;
        t2BRAS[br][as] = t2ABRS[ab][rs];
      }}
    }}

  free_block(t2ABRS);

  B_p_BB = get_DF_ints(BBintfile,BBlabel,noccB*noccB);
  double **B_p_RR = get_DF_ints(AAintfile,RRlabel,nvirA*nvirA);

  double **BRBR = block_matrix(noccB*nvirA,noccB*nvirA);

  for (int b=0, br=0; b < noccB; b++) {
    for (int r=0; r < nvirA; r++, br++) {
      C_DGEMM('N','T',noccB,nvirA,ribasis_->nbf()+3,1.0,&(B_p_BB[b*noccB][0]),
        ribasis_->nbf()+3,&(B_p_RR[r*nvirA][0]),ribasis_->nbf()+3,0.0,
        &(BRBR[br][0]),nvirA);
  }}

  free_block(B_p_BB);
  free_block(B_p_RR);

  C_DGEMM('N','N',noccB*nvirA,noccA*nvirB,noccB*nvirA,-1.0,&(BRBR[0][0]),
    noccB*nvirA,&(tBRAS[0][0]),noccA*nvirB,1.0,&(t2BRAS[0][0]),noccA*nvirB);

  free_block(BRBR);

  B_p_AA = get_DF_ints(AAintfile,AAlabel,noccA*noccA);
  double **B_p_SS = get_DF_ints(BBintfile,SSlabel,nvirB*nvirB);

  double **ASAS = block_matrix(noccA*nvirB,noccA*nvirB);

  for (int a=0, as=0; a < noccA; a++) {
    for (int s=0; s < nvirB; s++, as++) {
      C_DGEMM('N','T',noccA,nvirB,ribasis_->nbf()+3,1.0,&(B_p_AA[a*noccA][0]),
        ribasis_->nbf()+3,&(B_p_SS[s*nvirB][0]),ribasis_->nbf()+3,0.0,
        &(ASAS[as][0]),nvirB);
  }}

  free_block(B_p_AA);
  free_block(B_p_SS);

  C_DGEMM('N','N',noccB*nvirA,noccA*nvirB,noccA*nvirB,-1.0,&(tBRAS[0][0]),
    noccA*nvirB,&(ASAS[0][0]),noccA*nvirB,1.0,&(t2BRAS[0][0]),noccA*nvirB);

  free_block(ASAS);

  double energy = 4.0*C_DDOT(noccA*noccB*nvirA*nvirB,&(tBRAS[0][0]),1,
    &(t2BRAS[0][0]),1);

  free_block(tBRAS);
  free_block(t2BRAS);

  return(energy);
}

}}
