/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "sapt2p3.h"

namespace psi { namespace sapt {

void SAPT2p3::disp30()
{
  if (third_order_) {
    double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
      foccA_,noccA_,0,nvirA_);
    double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
      foccB_,noccB_,0,nvirB_);
  
    double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
    double **vARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
    psio_->read_entry(PSIF_SAPT_AMPS,"Disp30 uARBS Amplitudes",(char *)
      tARBS[0],sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);
  
    C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_+3,1.0,B_p_AR[0],ndf_+3,
      B_p_BS[0],ndf_+3,0.0,vARBS[0],aoccB_*nvirB_);
  
    e_disp30_ = 4.0*C_DDOT((long int) aoccA_*nvirA_*aoccB_*nvirB_,
      vARBS[0],1,tARBS[0],1);
  
    free_block(B_p_AR);
    free_block(B_p_BS);
    free_block(vARBS);
    free_block(tARBS);
  }
  else {
    double e1 = disp30_1(PSIF_SAPT_AMPS,"tARBS Amplitudes",
      PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",PSIF_SAPT_BB_DF_INTS,
      "SS RI Integrals",foccA_,noccA_,nvirA_,foccB_,noccB_,nvirB_);
    double e2 = disp30_2(PSIF_SAPT_AMPS,"tARBS Amplitudes",
      PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","RR RI Integrals",
      PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","SS RI Integrals",
      foccA_,noccA_,nvirA_,foccB_,noccB_,nvirB_);
  
    e_disp30_ = e1+e2;
  }

  if (print_) {
    outfile->Printf("    Disp30              = %18.12lf [Eh]\n",e_disp30_);
    
  }
}

double SAPT2p3::disp30_1(int ampfile, const char *amplabel, int AAintfile,
  const char *RRlabel, int BBintfile, const char *SSlabel, int foccA, 
  int noccA, int nvirA, int foccB, int noccB, int nvirB)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);
  double **tRSAB = block_matrix(nvirA*nvirB,aoccA*aoccB);

    for (int a=0, ar=0; a < aoccA; a++) {
    for (int r=0; r < nvirA; r++, ar++) {
      for (int b=0, bs=0; b < aoccB; b++) {
      for (int s=0; s < nvirB; s++, bs++) {
        int ab = a*aoccB + b;
        int rs = r*nvirB + s;
        int sr = s*nvirA + r;
        tRSAB[rs][ab] = tARBS[ar][bs];
      }}
    }}

  free_block(tARBS);

  double energy = 0.0;

  psio_address next_DF_RR = PSIO_ZERO;
  psio_address next_DF_SS = PSIO_ZERO;

  double **B_p_RR = block_matrix(nvirA*(nvirA+1)/2,ndf_+3);
  double **B_p_SS = block_matrix(nvirB*(nvirB+1)/2,ndf_+3);

  for (int r1=0,r1r2=0; r1 < nvirA; r1++) {
  for (int r2=0; r2 <= r1; r2++,r1r2++) {
    next_DF_RR = psio_get_address(PSIO_ZERO,sizeof(double)*(r1*nvirA+r2)*
      (ndf_+3));
    psio_->read(AAintfile,RRlabel,(char *) &(B_p_RR[r1r2][0]),sizeof(double)*
      (ndf_+3),next_DF_RR,&next_DF_RR);
    if (r1 != r2) C_DSCAL(ndf_+3,2.0,B_p_RR[r1r2],1);
  }}

  for (int s1=0,s1s2=0; s1 < nvirB; s1++) {
  for (int s2=0; s2 <= s1; s2++,s1s2++) {
    next_DF_SS = psio_get_address(PSIO_ZERO,sizeof(double)*(s1*nvirB+s2)*
      (ndf_+3));
    psio_->read(BBintfile,SSlabel,(char *) &(B_p_SS[s1s2][0]),sizeof(double)*
      (ndf_+3),next_DF_SS,&next_DF_SS);
    if (s1 != s2) C_DSCAL(ndf_+3,2.0,B_p_SS[s1s2],1);
  }}

  double **xRS = block_matrix(nvirA,nvirB*nvirB);
  double **yRS = block_matrix(nvirA,nvirB*(nvirB+1)/2);
  double *zSS = init_array(nvirB*(nvirB+1)/2);

  for (int r1=0; r1 < nvirA; r1++) {
    C_DGEMM('N','T',(r1+1)*nvirB,nvirB,aoccA*aoccB,1.0,tRSAB[0],aoccA*aoccB,
      tRSAB[r1*nvirB],aoccA*aoccB,0.0,xRS[0],nvirB);
    C_DGEMM('N','T',(r1+1),nvirB*(nvirB+1)/2,ndf_+3,1.0,
      B_p_RR[ioff_[r1]],ndf_+3,B_p_SS[0],
      ndf_+3,0.0,yRS[0],nvirB*(nvirB+1)/2);
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
  const char *SSlabel, int foccA, int noccA, int nvirA, int foccB, int noccB, 
  int nvirB)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);
  double **tABRS = block_matrix(aoccA*aoccB,nvirA*nvirB);

    for (int a=0, ar=0; a < aoccA; a++) {
    for (int r=0; r < nvirA; r++, ar++) {
      for (int b=0, bs=0; b < aoccB; b++) {
      for (int s=0; s < nvirB; s++, bs++) {
        int ab = a*aoccB + b;
        int rs = r*nvirB + s;
        tABRS[ab][rs] = tARBS[ar][bs];
      }}
    }}

  free_block(tARBS);

  double **t2ABRS = block_matrix(aoccA*aoccB,nvirA*nvirB);

  double **B_p_AA = get_DF_ints(AAintfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_BB = get_DF_ints(BBintfile,BBlabel,foccB,noccB,foccB,noccB);

  double **ABAB = block_matrix(aoccA*aoccB,aoccA*aoccB);

  for (int a=0, ab=0; a < aoccA; a++) {
    for (int b=0; b < aoccB; b++,ab++) {
      C_DGEMM('N','T',aoccA,aoccB,ndf_+3,1.0,&(B_p_AA[a*aoccA][0]),
        ndf_+3,&(B_p_BB[b*aoccB][0]),ndf_+3,0.0,
        &(ABAB[ab][0]),aoccB);
  }}

  free_block(B_p_AA);
  free_block(B_p_BB);

  C_DGEMM('N','N',aoccA*aoccB,nvirA*nvirB,aoccA*aoccB,1.0,&(ABAB[0][0]),
    aoccA*aoccB,&(tABRS[0][0]),nvirA*nvirB,1.0,&(t2ABRS[0][0]),nvirA*nvirB);

  free_block(ABAB);

  double **tBRAS = block_matrix(aoccB*nvirA,aoccA*nvirB);

    for (int a=0, ab=0; a < aoccA; a++) {
    for (int b=0; b < aoccB; b++, ab++) {
      for (int r=0, rs=0; r < nvirA; r++) {
      for (int s=0; s < nvirB; s++, rs++) {
        int br = b*nvirA + r;
        int as = a*nvirB + s;
        tBRAS[br][as] = tABRS[ab][rs];
      }}
    }}

  free_block(tABRS);

  double **t2BRAS = block_matrix(aoccB*nvirA,aoccA*nvirB);

    for (int a=0, ab=0; a < aoccA; a++) {
    for (int b=0; b < aoccB; b++, ab++) {
      for (int r=0, rs=0; r < nvirA; r++) {
      for (int s=0; s < nvirB; s++, rs++) {
        int br = b*nvirA + r;
        int as = a*nvirB + s;
        t2BRAS[br][as] = t2ABRS[ab][rs];
      }}
    }}

  free_block(t2ABRS);

  B_p_BB = get_DF_ints(BBintfile,BBlabel,foccB,noccB,foccB,noccB);
  double **B_p_RR = get_DF_ints(AAintfile,RRlabel,0,nvirA,0,nvirA);

  double **BRBR = block_matrix(aoccB*nvirA,aoccB*nvirA);

  for (int b=0, br=0; b < aoccB; b++) {
    for (int r=0; r < nvirA; r++, br++) {
      C_DGEMM('N','T',aoccB,nvirA,ndf_+3,1.0,&(B_p_BB[b*aoccB][0]),
        ndf_+3,&(B_p_RR[r*nvirA][0]),ndf_+3,0.0,
        &(BRBR[br][0]),nvirA);
  }}

  free_block(B_p_BB);
  free_block(B_p_RR);

  C_DGEMM('N','N',aoccB*nvirA,aoccA*nvirB,aoccB*nvirA,-1.0,&(BRBR[0][0]),
    aoccB*nvirA,&(tBRAS[0][0]),aoccA*nvirB,1.0,&(t2BRAS[0][0]),aoccA*nvirB);

  free_block(BRBR);

  B_p_AA = get_DF_ints(AAintfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_SS = get_DF_ints(BBintfile,SSlabel,0,nvirB,0,nvirB);

  double **ASAS = block_matrix(aoccA*nvirB,aoccA*nvirB);

  for (int a=0, as=0; a < aoccA; a++) {
    for (int s=0; s < nvirB; s++, as++) {
      C_DGEMM('N','T',aoccA,nvirB,ndf_+3,1.0,&(B_p_AA[a*aoccA][0]),
        ndf_+3,&(B_p_SS[s*nvirB][0]),ndf_+3,0.0,
        &(ASAS[as][0]),nvirB);
  }}

  free_block(B_p_AA);
  free_block(B_p_SS);

  C_DGEMM('N','N',aoccB*nvirA,aoccA*nvirB,aoccA*nvirB,-1.0,&(tBRAS[0][0]),
    aoccA*nvirB,&(ASAS[0][0]),aoccA*nvirB,1.0,&(t2BRAS[0][0]),aoccA*nvirB);

  free_block(ASAS);

  double energy = 4.0*C_DDOT((long int) aoccA*aoccB*nvirA*nvirB,
    &(tBRAS[0][0]),1,&(t2BRAS[0][0]),1);

  free_block(tBRAS);
  free_block(t2BRAS);

  return(energy);
}

}}
