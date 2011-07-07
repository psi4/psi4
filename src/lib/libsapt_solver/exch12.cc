#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::exch12()
{
  double e_exch111 = exch111();

  if (debug_) {
    fprintf(outfile,"    Exch111             = %18.12lf H\n",e_exch111);
  }

  double e_exch120_k2f = exch120_k2f();

  if (debug_) {
    fprintf(outfile,"    Exch120 K2f         = %18.12lf H\n",e_exch120_k2f);
  }

  double e_exch102_k2f = exch102_k2f();

  if (debug_) {
    fprintf(outfile,"    Exch102 K2f         = %18.12lf H\n",e_exch102_k2f);
  }

  e_exch12_ = e_exch111 + e_exch120_k2f + e_exch102_k2f; 

  if (print_) {
    fprintf(outfile,"    Exch12              = %18.12lf H\n",e_exch12_);
    fflush(outfile);
  }
}

double SAPT2::exch111()
{
  double e1 = 0.0, e2 = 0.0;

  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"Theta AR Intermediates",(char *) T_p_AR[0],
     sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"Theta BS Intermediates",(char *) T_p_BS[0],
     sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  double **C_p_AB = block_matrix(aoccA_*aoccB_,ndf_+3);
  double **D_p_AB = block_matrix(aoccA_*aoccB_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',aoccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][foccB_]),nmo_,
      T_p_AR[a*nvirA_],ndf_+3,0.0,C_p_AB[a*aoccB_],ndf_+3);
  }

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,nvirB_,1.0,&(sAB_[foccA_][noccB_]),nmo_,
      T_p_BS[b*nvirB_],ndf_+3,0.0,D_p_AB[b],aoccB_*(ndf_+3));
  }

  e1 -= 4.0*C_DDOT(aoccA_*aoccB_*(ndf_+3),C_p_AB[0],1,D_p_AB[0],1);

  free_block(C_p_AB);
  free_block(D_p_AB);

  double **C_p_AS = block_matrix(aoccA_*nvirB_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',nvirB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][noccB_]),nmo_,
      T_p_AR[a*nvirA_],ndf_+3,0.0,C_p_AS[a*nvirB_],ndf_+3);
  }

  free_block(T_p_AR);
  double **C_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('T','N',aoccB_,nvirB_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][foccB_]),
    nmo_,C_p_AS[0],nvirB_*(ndf_+3),0.0,C_p_BS[0],nvirB_*(ndf_+3));

  e2 -= 4.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),T_p_BS[0],1,C_p_BS[0],1);

  free_block(T_p_BS);
  free_block(C_p_AS);
  free_block(C_p_BS);

  if (debug_) {
    fprintf(outfile,"\n    Exch111_1           = %18.12lf H\n",e1);
    fprintf(outfile,"    Exch111_2           = %18.12lf H\n",e2);

  }

  return(e1+e2);
}

double SAPT2::exch120_k2f()
{
  double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0;
  double e5 = 0.0, e6 = 0.0, e7 = 0.0;

  double **tAR = block_matrix(aoccA_,nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"T2 AR Amplitudes",(char *) tAR[0],
     sizeof(double)*aoccA_*nvirA_);

  double **vAR = block_matrix(noccA_,nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"AR Exch12 K2f Integrals",(char *) vAR[0],
     sizeof(double)*noccA_*nvirA_);

  e1 -= 2.0*C_DDOT(aoccA_*nvirA_,tAR[0],1,vAR[foccA_],1);

  free_block(vAR);

  double **B_p_RB = get_RB_ints(2);
  double **B_p_AB = get_AB_ints(2);

  double **C_p_AB = block_matrix(aoccA_*noccB_,ndf_+3);

  C_DGEMM('N','N',aoccA_,noccB_*(ndf_+3),nvirA_,1.0,tAR[0],nvirA_,B_p_RB[0],
    noccB_*(ndf_+3),0.0,C_p_AB[0],noccB_*(ndf_+3));

  free_block(B_p_RB);

  e2 -= 2.0*C_DDOT(aoccA_*noccB_*(ndf_+3),B_p_AB[foccA_*noccB_],1,C_p_AB[0],1);

  double **C_p_BB = block_matrix(noccB_*noccB_,ndf_+3);

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][0]),nmo_,
    C_p_AB[0],noccB_*(ndf_+3),0.0,C_p_BB[0],noccB_*(ndf_+3));

  double **B_p_BB = get_BB_ints(1);

  e3 += 2.0*C_DDOT(noccB_*noccB_*(ndf_+3),B_p_BB[0],1,C_p_BB[0],1);

  free_block(B_p_BB);
  free_block(C_p_BB);

  double **xAB = block_matrix(aoccA_,noccB_);

  C_DGEMV('n',aoccA_*noccB_,ndf_+3,1.0,C_p_AB[0],ndf_+3,diagBB_,1,0.0,
    xAB[0],1);

  free_block(C_p_AB);

  for (int a=0; a<aoccA_; a++) {
    e4 -= 4.0*C_DDOT(noccB_,xAB[a],1,sAB_[a+foccA_],1);
  }

  C_DGEMV('n',aoccA_*noccB_,ndf_+3,1.0,B_p_AB[foccA_*noccB_],ndf_+3,
    diagAA_,1,0.0,xAB[0],1);

  double **yAB = block_matrix(aoccA_,noccB_);

  C_DGEMM('N','N',aoccA_,noccB_,nvirA_,1.0,tAR[0],nvirA_,&(sAB_[noccA_][0]),
    nmo_,0.0,yAB[0],noccB_);

  e5 -= 4.0*C_DDOT(aoccA_*noccB_,xAB[0],1,yAB[0],1);

  free_block(xAB);

  double **B_p_AA = get_AA_ints(1);
  double **D_p_AB = block_matrix(noccA_*noccB_,ndf_+3);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,aoccA_,1.0,yAB[0],noccB_,
      B_p_AA[a*noccA_+foccA_],ndf_+3,0.0,D_p_AB[a*noccB_],ndf_+3);
  }

  e6 += 2.0*C_DDOT(noccA_*noccB_*(ndf_+3),B_p_AB[0],1,D_p_AB[0],1);

  free_block(yAB);
  free_block(B_p_AA);
  free_block(D_p_AB);

  double **B_p_AR = get_AR_ints(1);
  double **C_p_AA = block_matrix(aoccA_*noccA_,ndf_+3);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,nvirA_,1.0,tAR[0],nvirA_,B_p_AR[a*nvirA_],
      ndf_+3,0.0,C_p_AA[a],noccA_*(ndf_+3));
  }

  free_block(B_p_AR);

  double **D_p_AA = block_matrix(aoccA_*noccA_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',noccA_,ndf_+3,noccB_,1.0,&(sAB_[0][0]),nmo_,
      B_p_AB[(a+foccA_)*noccB_],ndf_+3,0.0,D_p_AA[a*noccA_],ndf_+3);
  }

  e7 += 2.0*C_DDOT(aoccA_*noccA_*(ndf_+3),C_p_AA[0],1,D_p_AA[0],1);

  free_block(B_p_AB);
  free_block(C_p_AA);
  free_block(D_p_AA);
  free_block(tAR);

  if (debug_) {
    fprintf(outfile,"\n    Exch12_k2f_1        = %18.12lf H\n",e1);
    fprintf(outfile,"    Exch12_k2f_2        = %18.12lf H\n",e2);
    fprintf(outfile,"    Exch12_k2f_3        = %18.12lf H\n",e3);
    fprintf(outfile,"    Exch12_k2f_4        = %18.12lf H\n",e4);
    fprintf(outfile,"    Exch12_k2f_5        = %18.12lf H\n",e5);
    fprintf(outfile,"    Exch12_k2f_6        = %18.12lf H\n",e6);
    fprintf(outfile,"    Exch12_k2f_7        = %18.12lf H\n",e7);
  }

  return(e1+e2+e3+e4+e5+e6+e7);
}

double SAPT2::exch102_k2f()
{
  double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0;
  double e5 = 0.0, e6 = 0.0, e7 = 0.0;

  double **tBS = block_matrix(aoccB_,nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"T2 BS Amplitudes",(char *) tBS[0],
     sizeof(double)*aoccB_*nvirB_);

  double **vBS = block_matrix(noccB_,nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"BS Exch12 K2f Integrals",(char *) vBS[0],
     sizeof(double)*noccB_*nvirB_);

  e1 -= 2.0*C_DDOT(aoccB_*nvirB_,tBS[0],1,vBS[foccB_],1);

  free_block(vBS);

  double **B_p_AS = get_AS_ints(2);
  double **B_p_AB = get_AB_ints(1);

  double **C_p_AB = block_matrix(noccA_*aoccB_,ndf_+3);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',aoccB_,ndf_+3,nvirB_,1.0,tBS[0],nvirB_,B_p_AS[a*nvirB_],
      ndf_+3,0.0,C_p_AB[a*aoccB_],ndf_+3);
  }

  free_block(B_p_AS);

  for (int a=0; a<noccA_; a++) {
    e2 -= 2.0*C_DDOT(aoccB_*(ndf_+3),B_p_AB[a*noccB_+foccB_],1,
      C_p_AB[a*aoccB_],1);
  }

  double **C_p_AA = block_matrix(noccA_*noccA_,ndf_+3);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',noccA_,ndf_+3,aoccB_,1.0,&(sAB_[0][foccB_]),nmo_,
      C_p_AB[a*aoccB_],ndf_+3,0.0,C_p_AA[a*noccA_],ndf_+3);
  }

  double **B_p_AA = get_AA_ints(1);

  e3 += 2.0*C_DDOT(noccA_*noccA_*(ndf_+3),B_p_AA[0],1,C_p_AA[0],1);

  free_block(B_p_AA);
  free_block(C_p_AA);

  double **xAB = block_matrix(noccA_,aoccB_);

  C_DGEMV('n',noccA_*aoccB_,ndf_+3,1.0,C_p_AB[0],ndf_+3,diagAA_,1,0.0,
    xAB[0],1);

  free_block(C_p_AB);

  for (int a=0; a<noccA_; a++) {
    e4 -= 4.0*C_DDOT(aoccB_,xAB[a],1,&(sAB_[a][foccB_]),1);
  }

  for (int a=0; a<noccA_; a++) {
    C_DGEMV('n',aoccB_,ndf_+3,1.0,B_p_AB[a*noccB_+foccB_],ndf_+3,
      diagBB_,1,0.0,xAB[a],1);
  }

  double **yAB = block_matrix(noccA_,aoccB_);

  C_DGEMM('N','T',noccA_,aoccB_,nvirB_,1.0,&(sAB_[0][noccB_]),nmo_,tBS[0],
    nvirB_,0.0,yAB[0],aoccB_);

  e5 -= 4.0*C_DDOT(noccA_*aoccB_,xAB[0],1,yAB[0],1);

  free_block(xAB);

  double **B_p_BB = get_BB_ints(1);
  double **D_p_AB = block_matrix(noccA_*noccB_,ndf_+3);

  C_DGEMM('N','N',noccA_,noccB_*(ndf_+3),aoccB_,1.0,yAB[0],aoccB_,
    B_p_BB[foccB_*noccB_],noccB_*(ndf_+3),0.0,D_p_AB[0],noccB_*(ndf_+3));

  e6 += 2.0*C_DDOT(noccA_*noccB_*(ndf_+3),B_p_AB[0],1,D_p_AB[0],1);

  free_block(yAB);
  free_block(B_p_BB);
  free_block(D_p_AB);

  double **B_p_BS = get_BS_ints(1);
  double **C_p_BB = block_matrix(aoccB_*noccB_,ndf_+3);

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('N','N',aoccB_,ndf_+3,nvirB_,1.0,tBS[0],nvirB_,B_p_BS[b*nvirB_],
      ndf_+3,0.0,C_p_BB[b],noccB_*(ndf_+3));
  }

  free_block(B_p_BS);

  double **D_p_BB = block_matrix(aoccB_*noccB_,ndf_+3);

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('T','N',noccB_,ndf_+3,noccA_,1.0,&(sAB_[0][0]),nmo_,
      B_p_AB[b+foccB_],noccB_*(ndf_+3),0.0,D_p_BB[b*noccB_],ndf_+3);
  }

  e7 += 2.0*C_DDOT(aoccB_*noccB_*(ndf_+3),C_p_BB[0],1,D_p_BB[0],1);

  free_block(B_p_AB);
  free_block(C_p_BB);
  free_block(D_p_BB);
  free_block(tBS);

  if (debug_) {
    fprintf(outfile,"\n    Exch12_k2f_1        = %18.12lf H\n",e1);
    fprintf(outfile,"    Exch12_k2f_2        = %18.12lf H\n",e2);
    fprintf(outfile,"    Exch12_k2f_3        = %18.12lf H\n",e3);
    fprintf(outfile,"    Exch12_k2f_4        = %18.12lf H\n",e4);
    fprintf(outfile,"    Exch12_k2f_5        = %18.12lf H\n",e5);
    fprintf(outfile,"    Exch12_k2f_6        = %18.12lf H\n",e6);
    fprintf(outfile,"    Exch12_k2f_7        = %18.12lf H\n",e7);
  }

  return(e1+e2+e3+e4+e5+e6+e7);
}

}}
