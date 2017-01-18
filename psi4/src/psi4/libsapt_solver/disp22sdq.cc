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

#include "sapt2p.h"

namespace psi { namespace sapt {

void SAPT2p::disp22sdq()
{
  double e_disp211 = disp211();

  if (debug_) {
    outfile->Printf("    Disp211             = %18.12lf [Eh]\n",e_disp211);
    
  }

  double e_disp220s = disp220s(PSIF_SAPT_AMPS,"T2 AR Amplitudes",
    "T AR Intermediates",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    "RR RI Integrals",foccA_,noccA_,nvirA_);

  if (debug_) {
    outfile->Printf("    Disp220 (S)         = %18.12lf [Eh]\n",e_disp220s);
    
  }

  double e_disp202s = disp220s(PSIF_SAPT_AMPS,"T2 BS Amplitudes",
    "T BS Intermediates",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    "SS RI Integrals",foccB_,noccB_,nvirB_);

  if (debug_) {
    outfile->Printf("    Disp202 (S)         = %18.12lf [Eh]\n",e_disp202s);
    
  }

  double e_disp220d = disp220d_1(PSIF_SAPT_AMPS,"t2ARAR Amplitudes",
    "T AR Intermediates",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccA_,noccA_,nvirA_);
  e_disp220d += disp220d_2(PSIF_SAPT_AMPS,"gARAR x tARBS",
    "Theta AR Intermediates",PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccA_,noccA_,nvirA_,foccB_,noccB_,nvirB_,evalsA_,evalsB_,'N');

  if (debug_) {
    outfile->Printf("    Disp220 (D)         = %18.12lf [Eh]\n",e_disp220d);
    
  }

  double e_disp202d = disp220d_1(PSIF_SAPT_AMPS,"t2BSBS Amplitudes",
    "T BS Intermediates",PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccB_,noccB_,nvirB_);
  e_disp202d += disp220d_2(PSIF_SAPT_AMPS,"gBSBS x tARBS",
    "Theta BS Intermediates",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccB_,noccB_,nvirB_,foccA_,noccA_,nvirA_,evalsB_,evalsA_,'T');

  if (debug_) {
    outfile->Printf("    Disp202 (D)         = %18.12lf [Eh]\n",e_disp202d);
    
  }

  double e_disp220q = disp220q_1(PSIF_SAPT_AMPS,"tARAR Amplitudes",
    "T AR Intermediates","Theta AR Intermediates",aoccA_,nvirA_);
  e_disp220q += disp220q_2(PSIF_SAPT_AMPS,"pAA Density Matrix",
    "pRR Density Matrix","T AR Intermediates",PSIF_SAPT_AA_DF_INTS,
    "AR RI Integrals",foccA_,noccA_,nvirA_);
  e_disp220q += disp220q_3(PSIF_SAPT_AMPS,"tARAR Amplitudes",
    "tARBS Amplitudes",'N',PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccA_,noccA_,nvirA_,foccB_,noccB_,nvirB_);
  e_disp220q += disp220q_4(PSIF_SAPT_AMPS,"tARAR Amplitudes",
    "tARBS Amplitudes",'N',PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccA_,noccA_,nvirA_,foccB_,noccB_,nvirB_);

  if (debug_) {
    outfile->Printf("    Disp220 (Q)         = %18.12lf [Eh]\n",e_disp220q);
    
  }

  double e_disp202q = disp220q_1(PSIF_SAPT_AMPS,"tBSBS Amplitudes",
    "T BS Intermediates","Theta BS Intermediates",aoccB_,nvirB_);
  e_disp202q += disp220q_2(PSIF_SAPT_AMPS,"pBB Density Matrix",
    "pSS Density Matrix","T BS Intermediates",PSIF_SAPT_BB_DF_INTS,
    "BS RI Integrals",foccB_,noccB_,nvirB_);
  e_disp202q += disp220q_3(PSIF_SAPT_AMPS,"tBSBS Amplitudes",
    "tARBS Amplitudes",'T',PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccB_,noccB_,nvirB_,foccA_,noccA_,nvirA_);
  e_disp202q += disp220q_4(PSIF_SAPT_AMPS,"tBSBS Amplitudes",
    "tARBS Amplitudes",'T',PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccB_,noccB_,nvirB_,foccA_,noccA_,nvirA_);

  if (debug_) {
    outfile->Printf("    Disp202 (Q)         = %18.12lf [Eh]\n\n",e_disp202q);
    
  }

  e_disp22sdq_ = e_disp211 + e_disp220s + e_disp202s + e_disp220d +
    e_disp202d + e_disp220q + e_disp202q;

  if (print_) {
    outfile->Printf("    Disp22 (SDQ)        = %18.12lf [Eh]\n",e_disp22sdq_);
    
  }
}

double SAPT2p::disp211()
{
  double energy = 0.0;

  double **xARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  double **yARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"gBSBS x tARBS",(char *) xARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"gARAR x tARBS",(char *) yARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccA_,noccA_,0,nvirA_);
  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"Theta BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    T_p_BS[0],ndf_+3,1.0,xARBS[0],aoccB_*nvirB_);

  free_block(B_p_AR);

  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"Theta AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));
  double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccB_,noccB_,0,nvirB_);

  C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_+3,1.0,T_p_AR[0],ndf_+3,
    B_p_BS[0],ndf_+3,1.0,yARBS[0],aoccB_*nvirB_);

  free_block(B_p_BS);

  for (int a=0, ar=0; a<aoccA_; a++){
    for (int r=0; r<nvirA_; r++, ar++){
      for (int b=0, bs=0; b<aoccB_; b++){
        for (int s=0; s<nvirB_; s++, bs++){
          xARBS[ar][bs] /= evalsA_[a+foccA_]+evalsB_[b+foccB_]
            -evalsA_[r+noccA_]-evalsB_[s+noccB_];
  }}}}

  energy = 8.0*C_DDOT(aoccA_*nvirA_*aoccB_*nvirB_,xARBS[0],1,yARBS[0],1);

  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *) xARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_+3,1.0,T_p_AR[0],ndf_+3,
    T_p_BS[0],ndf_+3,0.0,yARBS[0],aoccB_*nvirB_);

  energy += 8.0*C_DDOT(aoccA_*nvirA_*aoccB_*nvirB_,xARBS[0],1,yARBS[0],1);

  free_block(xARBS);
  free_block(yARBS);
  free_block(T_p_AR);
  free_block(T_p_BS);

  return(energy);
}

double SAPT2p::disp220s(int ampfile, const char *tlabel, 
  const char *thetalabel, int intfile, const char *AAlabel, 
  const char *RRlabel, int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **yAR = block_matrix(aoccA,nvirA);

  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),1.0,T_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),0.0,yAR[0],nvirA);

  free_block(B_p_RR);

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);

  for (int a=0; a<aoccA; a++){
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-1.0,B_p_AA[a*aoccA],ndf_+3,
      T_p_AR[a*nvirA],ndf_+3,1.0,yAR[0],nvirA);
  }

  free_block(B_p_AA);
  free_block(T_p_AR);

  double **tAR = block_matrix(aoccA,nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tAR[0], 
    sizeof(double)*aoccA*nvirA);

  double energy = 8.0*C_DDOT(aoccA*nvirA,tAR[0],1,yAR[0],1);

  free_block(tAR);
  free_block(yAR);

  return(energy);
}

double SAPT2p::disp220d_1(int ampfile, const char *tlabel,
  const char *thetalabel, int intfile, const char *ARlabel,
  int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double *xARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,T_p_AR[0],ndf_+3,
    B_p_AR[0],ndf_+3,0.0,xARAR,aoccA*nvirA);

  symmetrize(xARAR,aoccA,nvirA);
  antisym(xARAR,aoccA,nvirA);

  free_block(B_p_AR);
  free_block(T_p_AR);

  double *t2ARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) t2ARAR,
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);

  double energy = 4.0*C_DDOT((long int) aoccA*nvirA*aoccA*nvirA,xARAR,1,
    t2ARAR,1);

  free(t2ARAR);
  free(xARAR);

  if (debug_) {
    outfile->Printf("\n    Disp22d_1           = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

double SAPT2p::disp220d_2(int ampfile, const char *glabel,
  const char *thetalabel, int intfile, const char *BSlabel,
  int foccA, int noccA, int nvirA, int foccB, int noccB, int nvirB,
  double *evalsA, double *evalsB, const char trans)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double energy = 0.0;

  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));
  double **B_p_BS = get_DF_ints(intfile,BSlabel,foccB,noccB,0,nvirB);

  if (trans =='n' || trans == 'N') {
    double **yARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
    psio_->read_entry(ampfile,glabel,(char *) yARBS[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    C_DGEMM('N','T',aoccA*nvirA,aoccB*nvirB,ndf_+3,1.0,T_p_AR[0],ndf_+3,
      B_p_BS[0],ndf_+3,1.0,yARBS[0],aoccB*nvirB);

    for (int a=0, ar=0; a<aoccA; a++){
      for (int r=0; r<nvirA; r++, ar++){
        for (int b=0, bs=0; b<aoccB; b++){
          for (int s=0; s<nvirB; s++, bs++){
            double tval = yARBS[ar][bs];
            tval *= tval;
            energy += 4.0*tval/(evalsA[a+foccA]+evalsB[b+foccB]
              -evalsA[r+noccA]-evalsB[s+noccB]);
    }}}}

    free_block(yARBS);
  }
  else if (trans =='t' || trans == 'T') {
    double **yBSAR = block_matrix(aoccB*nvirB,aoccA*nvirA);
    psio_->read_entry(ampfile,glabel,(char *) yBSAR[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    C_DGEMM('N','T',aoccB*nvirB,aoccA*nvirA,ndf_+3,1.0,B_p_BS[0],ndf_+3,
      T_p_AR[0],ndf_+3,1.0,yBSAR[0],aoccA*nvirA);
  
    for (int b=0, bs=0; b<aoccB; b++){
      for (int s=0; s<nvirB; s++, bs++){
        for (int a=0, ar=0; a<aoccA; a++){
          for (int r=0; r<nvirA; r++, ar++){
            double tval = yBSAR[bs][ar];
            tval *= tval;
            energy += 4.0*tval/(evalsA[a+foccA]+evalsB[b+foccB]
              -evalsA[r+noccA]-evalsB[s+noccB]);
    }}}}
  
    free_block(yBSAR);
  }
  else
    throw PsiException("You want me to do what to that matrix?",
       __FILE__,__LINE__);

  free_block(T_p_AR);
  free_block(B_p_BS);

  if (debug_) {
    outfile->Printf("    Disp22d_2           = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

double SAPT2p::disp220q_1(int ampfile, const char *tlabel, 
  const char *Tlabel, const char *thetalabel, int aoccA, int nvirA)
{
  double energy = 0.0;

  double **thetaARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) thetaARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  antisym(thetaARAR,aoccA,nvirA);

  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,Tlabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double **theta_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) theta_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double **xARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,T_p_AR[0],ndf_+3,
    theta_p_AR[0],ndf_+3,0.0,xARAR[0],aoccA*nvirA);

  energy = 4.0*C_DDOT((long int) aoccA*nvirA*aoccA*nvirA,xARAR[0],1,
    thetaARAR[0],1);

  free_block(T_p_AR);
  free_block(theta_p_AR);
  free_block(thetaARAR);
  free_block(xARAR);

  if (debug_) {
    outfile->Printf("\n    Disp22q_1           = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

double SAPT2p::disp220q_2(int ampfile, const char *pAAlabel,
  const char *pRRlabel, const char *Tlabel, int intfile, const char *ARlabel,
  int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **pAA = block_matrix(aoccA,aoccA);
  double **pRR = block_matrix(nvirA,nvirA);

  psio_->read_entry(ampfile,pAAlabel,(char *) pAA[0],
    sizeof(double)*aoccA*aoccA);
  psio_->read_entry(ampfile,pRRlabel,(char *) pRR[0],
    sizeof(double)*nvirA*nvirA);

  double **qAA = block_matrix(aoccA,aoccA);
  double **qRR = block_matrix(nvirA,nvirA);

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,Tlabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  C_DGEMM('N','T',aoccA,aoccA,nvirA*(ndf_+3),1.0,B_p_AR[0],nvirA*(ndf_+3),
    T_p_AR[0],nvirA*(ndf_+3),0.0,qAA[0],aoccA);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',nvirA,nvirA,ndf_+3,1.0,B_p_AR[a*nvirA],ndf_+3,
      T_p_AR[a*nvirA],ndf_+3,1.0,qRR[0],nvirA);
  }

  free_block(B_p_AR);
  free_block(T_p_AR);

  double energy = -4.0*C_DDOT(aoccA*aoccA,pAA[0],1,qAA[0],1);
  energy -= 4.0*C_DDOT(nvirA*nvirA,pRR[0],1,qRR[0],1);

  free_block(pAA);
  free_block(pRR);
  free_block(qAA);
  free_block(qRR);

  if (debug_) {
    outfile->Printf("    Disp22q_2           = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

double SAPT2p::disp220q_3(int ampfile, const char *tARARlabel,
  const char *tARBSlabel, const char trans, int intfile, const char *ARlabel, 
  int foccA, int noccA, int nvirA, int foccB, int noccB, int nvirB)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **xARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);

  if (trans =='n' || trans == 'N') {
    double **tARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
    psio_->read_entry(ampfile,tARBSlabel,(char *) tARBS[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccB*nvirB,1.0,tARBS[0],
      aoccB*nvirB,tARBS[0],aoccB*nvirB,0.0,xARAR[0],aoccA*nvirA);

    free_block(tARBS);
  }
  else if (trans =='t' || trans == 'T') {
    double **tBSAR = block_matrix(aoccB*nvirB,aoccA*nvirA);
    psio_->read_entry(ampfile,tARBSlabel,(char *) tBSAR[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    C_DGEMM('T','N',aoccA*nvirA,aoccA*nvirA,aoccB*nvirB,1.0,tBSAR[0],
      aoccA*nvirA,tBSAR[0],aoccA*nvirA,0.0,xARAR[0],aoccA*nvirA);

    free_block(tBSAR);
  }
  else
    throw PsiException("You want me to do what to that matrix?",
       __FILE__,__LINE__);

  double **tARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,tARARlabel,(char *) tARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  antisym(tARAR,aoccA,nvirA);

  double **yARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,1.0,xARAR[0],
    aoccA*nvirA,tARAR[0],aoccA*nvirA,0.0,yARAR[0],aoccA*nvirA);

  free_block(tARAR);

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    B_p_AR[0],ndf_+3,0.0,xARAR[0],aoccA*nvirA);
  antisym(xARAR,aoccA,nvirA);

  double energy = 4.0*C_DDOT((long int) aoccA*nvirA*aoccA*nvirA,xARAR[0],1,
    yARAR[0],1);

  free_block(xARAR);
  free_block(yARAR);
  free_block(B_p_AR);

  if (debug_) {
    outfile->Printf("    Disp22q_3           = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

double SAPT2p::disp220q_4(int ampfile, const char *tARARlabel,
  const char *tARBSlabel, const char trans, int intfile, const char *ARlabel,
  int foccA, int noccA, int nvirA, int foccB, int noccB, int nvirB)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **rAA = block_matrix(aoccA,aoccA);
  double **rRR = block_matrix(nvirA,nvirA);

  double **tARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,tARARlabel,(char *) tARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);

  double **gARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    B_p_AR[0],ndf_+3,0.0,gARAR[0],aoccA*nvirA);
  antisym(gARAR,aoccA,nvirA);

  C_DGEMM('N','T',aoccA,aoccA,nvirA*aoccA*nvirA,1.0,tARAR[0],nvirA*aoccA*nvirA,
    gARAR[0],nvirA*aoccA*nvirA,0.0,rAA[0],aoccA);

  C_DGEMM('T','N',nvirA,nvirA,aoccA*nvirA*aoccA,1.0,tARAR[0],nvirA,
    gARAR[0],nvirA,0.0,rRR[0],nvirA);

  free_block(gARAR);
  free_block(tARAR);
  free_block(B_p_AR);

  double **sAA = block_matrix(aoccA,aoccA);
  double **sRR = block_matrix(nvirA,nvirA);

  if (trans =='n' || trans == 'N') {
    double **tARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
    psio_->read_entry(ampfile,tARBSlabel,(char *) tARBS[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    C_DGEMM('N','T',aoccA,aoccA,nvirA*aoccB*nvirB,1.0,tARBS[0],
      nvirA*aoccB*nvirB,tARBS[0],nvirA*aoccB*nvirB,0.0,sAA[0],aoccA);

    for (int a=0; a<aoccA; a++) {
      C_DGEMM('N','T',nvirA,nvirA,aoccB*nvirB,1.0,tARBS[a*nvirA],aoccB*nvirB,
        tARBS[a*nvirA],aoccB*nvirB,1.0,sRR[0],nvirA);
    }

    free_block(tARBS);
  }
  else if (trans =='t' || trans == 'T') {
    double **tBSAR = block_matrix(aoccB*nvirB,aoccA*nvirA);
    psio_->read_entry(ampfile,tARBSlabel,(char *) tBSAR[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    for (int b=0, bs=0; b<aoccB; b++) {
      for (int s=0; s<nvirB; s++, bs++) {
        C_DGEMM('N','T',aoccA,aoccA,nvirA,1.0,tBSAR[bs],nvirA,tBSAR[bs],nvirA,
          1.0,sAA[0],aoccA);
    }}

    C_DGEMM('T','N',nvirA,nvirA,aoccA*aoccB*nvirB,1.0,tBSAR[0],nvirA,
      tBSAR[0],nvirA,0.0,sRR[0],nvirA);

    free_block(tBSAR);
  }
  else
    throw PsiException("You want me to do what to that matrix?",
       __FILE__,__LINE__);

  double energy = -4.0*C_DDOT(aoccA*aoccA,rAA[0],1,sAA[0],1);
  energy -= 4.0*C_DDOT(nvirA*nvirA,rRR[0],1,sRR[0],1);

  free_block(rAA);
  free_block(rRR);
  free_block(sAA);
  free_block(sRR);

  if (debug_) {
    outfile->Printf("    Disp22q_4           = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

}}
