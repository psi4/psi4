#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::exch12()
{
  double e_exch111 = exch111();

  if (debug_) {
    fprintf(outfile,"    Exch111             = %18.12lf H\n",e_exch111);
    fflush(outfile);
  }

  double e_exch120_k2u = exch110(PSIF_SAPT_AMPS,"Theta 2 AR Intermediates");

  if (debug_) {
    fprintf(outfile,"    Exch120 K2u         = %18.12lf H\n",e_exch120_k2u);
    fflush(outfile);
  }

  double e_exch102_k2u = exch101(PSIF_SAPT_AMPS,"Theta 2 BS Intermediates");

  if (debug_) {
    fprintf(outfile,"    Exch102 K2u         = %18.12lf H\n",e_exch102_k2u);
    fflush(outfile);
  }

  double e_exch120_k2f = exch120_k2f();

  if (debug_) {
    fprintf(outfile,"    Exch120 K2f         = %18.12lf H\n",e_exch120_k2f);
    fflush(outfile);
  }

  double e_exch102_k2f = exch102_k2f();

  if (debug_) {
    fprintf(outfile,"    Exch102 K2f         = %18.12lf H\n",e_exch102_k2f);
    fflush(outfile);
  }

  double e_exch120_k11u = exch120_k11u_1();
  e_exch120_k11u += exch120_k11u_2();
  e_exch120_k11u += exch120_k11u_3();
  e_exch120_k11u += exch120_k11u_4();
  e_exch120_k11u += exch120_k11u_5();
  e_exch120_k11u += exch120_k11u_6();

  if (debug_) {
    fprintf(outfile,"    Exch120 K11u        = %18.12lf H\n",e_exch120_k11u);
    fflush(outfile);
  }

  double e_exch102_k11u = exch102_k11u_1();
  e_exch102_k11u += exch102_k11u_2();
  e_exch102_k11u += exch102_k11u_3();
  e_exch102_k11u += exch102_k11u_4();
  e_exch102_k11u += exch102_k11u_5();
  e_exch102_k11u += exch102_k11u_6();

  if (debug_) {
    fprintf(outfile,"    Exch102 K11u        = %18.12lf H\n\n",e_exch102_k11u);
    fflush(outfile);
  }

  e_exch12_ = e_exch111 + e_exch120_k2f + e_exch102_k2f + e_exch120_k2u +
    e_exch102_k2u + e_exch120_k11u + e_exch102_k11u; 

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
    C_DGEMM('T','N',aoccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][foccB_]),nmoB_,
      T_p_AR[a*nvirA_],ndf_+3,0.0,C_p_AB[a*aoccB_],ndf_+3);
  }

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,nvirB_,1.0,&(sAB_[foccA_][noccB_]),nmoB_,
      T_p_BS[b*nvirB_],ndf_+3,0.0,D_p_AB[b],aoccB_*(ndf_+3));
  }

  e1 -= 4.0*C_DDOT(aoccA_*aoccB_*(ndf_+3),C_p_AB[0],1,D_p_AB[0],1);

  free_block(C_p_AB);
  free_block(D_p_AB);

  double **C_p_AS = block_matrix(aoccA_*nvirB_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',nvirB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][noccB_]),nmoB_,
      T_p_AR[a*nvirA_],ndf_+3,0.0,C_p_AS[a*nvirB_],ndf_+3);
  }

  free_block(T_p_AR);
  double **C_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('T','N',aoccB_,nvirB_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][foccB_]),
    nmoB_,C_p_AS[0],nvirB_*(ndf_+3),0.0,C_p_BS[0],nvirB_*(ndf_+3));

  e2 -= 4.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),T_p_BS[0],1,C_p_BS[0],1);

  free_block(T_p_BS);
  free_block(C_p_AS);
  free_block(C_p_BS);

  if (debug_) {
    fprintf(outfile,"\n    Exch111_1           = %18.12lf H\n",e1);
    fprintf(outfile,"    Exch111_2           = %18.12lf H\n",e2);
    fflush(outfile);

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

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][0]),nmoB_,
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
    nmoB_,0.0,yAB[0],noccB_);

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
    C_DGEMM('N','N',noccA_,ndf_+3,noccB_,1.0,&(sAB_[0][0]),nmoB_,
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
    fflush(outfile);
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
    C_DGEMM('N','N',noccA_,ndf_+3,aoccB_,1.0,&(sAB_[0][foccB_]),nmoB_,
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

  C_DGEMM('N','T',noccA_,aoccB_,nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,tBS[0],
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
    C_DGEMM('T','N',noccB_,ndf_+3,noccA_,1.0,&(sAB_[0][0]),nmoB_,
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
    fflush(outfile);
  }

  return(e1+e2+e3+e4+e5+e6+e7);
}

double SAPT2::exch120_k11u_1()
{
  double energy=0.0;

  double **pRR = block_matrix(nvirA_,nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"pRR Density Matrix",(char *) pRR[0],
    sizeof(double)*nvirA_*nvirA_);

  double **B_p_RB = get_RB_ints(1);
  double **C_p_RB = get_RB_ints(2);

  double **yRR = block_matrix(nvirA_,nvirA_);

  C_DGEMM('N','T',nvirA_,nvirA_,noccB_*(ndf_+3),1.0,&(B_p_RB[0][0]),
    noccB_*(ndf_+3),&(C_p_RB[0][0]),noccB_*(ndf_+3),0.0,&(yRR[0][0]),
    nvirA_);

  energy += 2.0*C_DDOT(nvirA_*nvirA_,pRR[0],1,yRR[0],1);

  free_block(yRR);

  double **D_p_RB = block_matrix(nvirA_*noccB_,(ndf_+3));

  C_DGEMM('N','N',nvirA_,noccB_*(ndf_+3),nvirA_,1.0,&(pRR[0][0]),nvirA_,
    &(B_p_RB[0][0]),noccB_*(ndf_+3),0.0,&(D_p_RB[0][0]),noccB_*(ndf_+3));

  free_block(B_p_RB);

  double **E_p_RB = block_matrix(nvirA_*noccB_,(ndf_+3));

  C_DGEMM('N','N',nvirA_,noccB_*(ndf_+3),nvirA_,1.0,&(pRR[0][0]),nvirA_,
    &(C_p_RB[0][0]),noccB_*(ndf_+3),0.0,&(E_p_RB[0][0]),noccB_*(ndf_+3));

  free_block(C_p_RB);

  double **B_p_AR = get_AR_ints(1);
  double **D_p_BR = block_matrix(nvirA_*noccB_,(ndf_+3));

  C_DGEMM('T','N',noccB_,nvirA_*(ndf_+3),noccA_,1.0,&(sAB_[0][0]),nmoB_,
    &(B_p_AR[0][0]),nvirA_*(ndf_+3),0.0,&(D_p_BR[0][0]),nvirA_*(ndf_+3));

  for (int b=0,br=0; b<noccB_; b++) {
    for (int r=0; r<nvirA_; r++,br++) {
      int rb = r*noccB_+b;
      energy -= 2.0*C_DDOT((ndf_+3),D_p_BR[br],1,D_p_RB[rb],1);
  }}

  double **xRB = block_matrix(nvirA_,noccB_);

  C_DGEMV('n',nvirA_*noccB_,(ndf_+3),1.0,&(D_p_RB[0][0]),(ndf_+3),diagAA_,1,
    0.0,&(xRB[0][0]),1);

  free_block(D_p_RB);

  for (int r=0; r<nvirA_; r++) {
    energy += 4.0*C_DDOT(noccB_,sAB_[r+noccA_],1,xRB[r],1);
  }

  C_DGEMV('n',nvirA_*noccB_,(ndf_+3),1.0,&(E_p_RB[0][0]),(ndf_+3),diagBB_,1,
    0.0,&(xRB[0][0]),1);

  for (int r=0; r<nvirA_; r++) {
    energy += 4.0*C_DDOT(noccB_,sAB_[r+noccA_],1,xRB[r],1);
  }

  free_block(xRB);

  double **B_p_BB = get_BB_ints(1);
  double **C_p_BB = block_matrix(noccB_*noccB_,(ndf_+3));

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),nvirA_,1.0,&(sAB_[noccA_][0]),nmoB_,
    &(E_p_RB[0][0]),noccB_*(ndf_+3),0.0,&(C_p_BB[0][0]),noccB_*(ndf_+3));

  free_block(E_p_RB);

  energy -= 2.0*C_DDOT((long int) noccB_*noccB_*(ndf_+3),B_p_BB[0],1,
    C_p_BB[0],1);

  free_block(C_p_BB);

  double **B_p_AB = get_AB_ints(2);
  double **yRB = block_matrix(nvirA_,noccB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirA_,noccB_,(ndf_+3),1.0,B_p_AR[a*nvirA_],(ndf_+3),
      B_p_AB[a*noccB_],(ndf_+3),1.0,yRB[0],noccB_);
  }

  free_block(B_p_AR);

  double **zRB = block_matrix(nvirA_,noccB_);

  C_DGEMM('N','N',nvirA_,noccB_,nvirA_,1.0,pRR[0],nvirA_,&(sAB_[noccA_][0]),
    nmoB_,0.0,zRB[0],noccB_);

  energy -= 2.0*C_DDOT(nvirA_*noccB_,yRB[0],1,zRB[0],1);

  free_block(yRB);

  double **xBR = block_matrix(noccB_,nvirA_);

  C_DGEMV('n',nvirA_*noccB_,(ndf_+3),1.0,&(D_p_BR[0][0]),(ndf_+3),diagBB_,1,
    0.0,&(xBR[0][0]),1);

  for (int b=0; b<noccB_; b++) {
    for (int r=0; r<nvirA_; r++) {
      energy -= 8.0*xBR[b][r]*zRB[r][b];
  }}

  free_block(xBR);

  double **D_p_BB = block_matrix(noccB_*noccB_,(ndf_+3));

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('T','N',noccB_,(ndf_+3),nvirA_,1.0,zRB[0],noccB_,
      D_p_BR[b*nvirA_],(ndf_+3),0.0,D_p_BB[b*noccB_],(ndf_+3));
  }

  free_block(D_p_BR);

  energy += 4.0*C_DDOT((long int) noccB_*noccB_*(ndf_+3),&(B_p_BB[0][0]),1,
    &(D_p_BB[0][0]),1);

  free_block(D_p_BB);

  double **zBB = block_matrix(noccB_,noccB_);

  C_DGEMM('T','N',noccB_,noccB_,nvirA_,1.0,sAB_[noccA_],nmoB_,zRB[0],
    noccB_,0.0,zBB[0],noccB_);

  double **wBB = block_matrix(noccB_,noccB_);

  C_DGEMV('n',noccB_*noccB_,(ndf_+3),1.0,&(B_p_BB[0][0]),(ndf_+3),diagAA_,1,
    0.0,&(wBB[0][0]),1);

  energy -= 4.0*C_DDOT(noccB_*noccB_,wBB[0],1,zBB[0],1);

  free_block(wBB);
  free_block(zBB);
  free_block(zRB);

  double **B_p_RR = get_RR_ints(1);

  double *X = init_array((ndf_+3));

  C_DGEMV('t',nvirA_*nvirA_,(ndf_+3),1.0,&(B_p_RR[0][0]),(ndf_+3),pRR[0],1,
    0.0,X,1);

  free_block(pRR);
  free_block(B_p_RR);

  double **xAB = block_matrix(noccA_,noccB_);

  C_DGEMV('n',noccA_*noccB_,(ndf_+3),1.0,&(B_p_AB[0][0]),(ndf_+3),X,1,0.0,
    &(xAB[0][0]),1);

  for (int a=0; a<noccA_; a++) {
    energy += 4.0*C_DDOT(noccB_,sAB_[a],1,xAB[a],1);
  }

  free_block(xAB);
  free_block(B_p_AB);

  double **xBB = block_matrix(noccB_,noccB_);
  double **yBB = block_matrix(noccB_,noccB_);

  C_DGEMV('n',noccB_*noccB_,(ndf_+3),1.0,&(B_p_BB[0][0]),(ndf_+3),X,1,0.0,
    &(xBB[0][0]),1);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,sAB_[0],nmoB_,sAB_[0],nmoB_,0.0,
    yBB[0],noccB_);

  energy -= 4.0*C_DDOT(noccB_*noccB_,xBB[0],1,yBB[0],1);

  free(X);
  free_block(xBB);
  free_block(yBB);
  free_block(B_p_BB);

  if (debug_) {
    fprintf(outfile,"\n    Exch12_k11u_1       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

double SAPT2::exch102_k11u_1()
{
  double energy=0.0;

  double **pSS = block_matrix(nvirB_,nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"pSS Density Matrix",(char *) pSS[0],
    sizeof(double)*nvirB_*nvirB_);

  double **B_p_AS = get_AS_ints(1);
  double **C_p_AS = get_AS_ints(2);

  double **ySS = block_matrix(nvirB_,nvirB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirB_,nvirB_,(ndf_+3),1.0,&(B_p_AS[a*nvirB_][0]),(ndf_+3),
      &(C_p_AS[a*nvirB_][0]),(ndf_+3),1.0,&(ySS[0][0]),nvirB_);
  }

  energy += 2.0*C_DDOT(nvirB_*nvirB_,pSS[0],1,ySS[0],1);

  free_block(ySS);

  double **D_p_AS = block_matrix(noccA_*nvirB_,(ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',nvirB_,(ndf_+3),nvirB_,1.0,&(pSS[0][0]),nvirB_,
      &(B_p_AS[a*nvirB_][0]),(ndf_+3),0.0,&(D_p_AS[a*nvirB_][0]),(ndf_+3));
  }

  free_block(B_p_AS);

  double **E_p_AS = block_matrix(noccA_*nvirB_,(ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',nvirB_,(ndf_+3),nvirB_,1.0,&(pSS[0][0]),nvirB_,
      &(C_p_AS[a*nvirB_][0]),(ndf_+3),0.0,&(E_p_AS[a*nvirB_][0]),(ndf_+3));
  }

  free_block(C_p_AS);

  double **B_p_BS = get_BS_ints(1);
  double **F_p_AS = block_matrix(noccA_*nvirB_,(ndf_+3));

  C_DGEMM('N','N',noccA_,nvirB_*(ndf_+3),noccB_,1.0,&(sAB_[0][0]),nmoB_,
    &(B_p_BS[0][0]),nvirB_*(ndf_+3),0.0,&(F_p_AS[0][0]),nvirB_*(ndf_+3));

  energy -= 2.0*C_DDOT((long int) noccA_*nvirB_*(ndf_+3),D_p_AS[0],1,
    F_p_AS[0],1);

  double **xAS = block_matrix(noccA_,nvirB_);

  C_DGEMV('n',noccA_*nvirB_,(ndf_+3),1.0,&(D_p_AS[0][0]),(ndf_+3),diagBB_,1,
    0.0,&(xAS[0][0]),1);

  free_block(D_p_AS);

  for (int a=0; a<noccA_; a++) {
    energy += 4.0*C_DDOT(nvirB_,&(sAB_[a][noccB_]),1,xAS[a],1);
  }

  C_DGEMV('n',noccA_*nvirB_,(ndf_+3),1.0,&(E_p_AS[0][0]),(ndf_+3),diagAA_,1,
    0.0,&(xAS[0][0]),1);

  for (int a=0; a<noccA_; a++) {
    energy += 4.0*C_DDOT(nvirB_,&(sAB_[a][noccB_]),1,xAS[a],1);
  }

  double **B_p_AA = get_AA_ints(1);
  double **C_p_AA = block_matrix(noccA_*noccA_,(ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',noccA_,(ndf_+3),nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
      &(E_p_AS[a*nvirB_][0]),(ndf_+3),0.0,&(C_p_AA[a*noccA_][0]),(ndf_+3));
  }

  free_block(E_p_AS);

  energy -= 2.0*C_DDOT((long int) noccA_*noccA_*(ndf_+3),B_p_AA[0],1,
    C_p_AA[0],1);

  free_block(C_p_AA);

  double **B_p_AB = get_AB_ints(1);
  double **yAS = block_matrix(noccA_,nvirB_);

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('N','T',noccA_,nvirB_,(ndf_+3),1.0,B_p_AB[b],noccB_*(ndf_+3),
      B_p_BS[b*nvirB_],(ndf_+3),1.0,yAS[0],nvirB_);
  }

  free_block(B_p_BS);

  double **zAS = block_matrix(noccA_,nvirB_);

  C_DGEMM('N','N',noccA_,nvirB_,nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,pSS[0],
    nvirB_,0.0,zAS[0],nvirB_);

  energy -= 2.0*C_DDOT(noccA_*nvirB_,yAS[0],1,zAS[0],1);

  free_block(yAS);

  C_DGEMV('n',noccA_*nvirB_,(ndf_+3),1.0,&(F_p_AS[0][0]),(ndf_+3),diagAA_,1,
    0.0,&(xAS[0][0]),1);

  energy -= 8.0*C_DDOT(noccA_*nvirB_,xAS[0],1,zAS[0],1);

  free_block(xAS);

  double **D_p_AA = block_matrix(noccA_*noccA_,(ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',noccA_,(ndf_+3),nvirB_,1.0,zAS[0],nvirB_,
      F_p_AS[a*nvirB_],(ndf_+3),0.0,D_p_AA[a*noccA_],(ndf_+3));
  }

  free_block(F_p_AS);

  energy += 4.0*C_DDOT((long int) noccA_*noccA_*(ndf_+3),&(B_p_AA[0][0]),1,
    &(D_p_AA[0][0]),1);

  free_block(D_p_AA);

  double **zAA = block_matrix(noccA_,noccA_);

  C_DGEMM('N','T',noccA_,noccA_,nvirB_,1.0,zAS[0],nvirB_,&(sAB_[0][noccB_]),
    nmoB_,0.0,zAA[0],noccA_);

  double **wAA = block_matrix(noccA_,noccA_);

  C_DGEMV('n',noccA_*noccA_,(ndf_+3),1.0,&(B_p_AA[0][0]),(ndf_+3),diagBB_,1,
    0.0,&(wAA[0][0]),1);

  energy -= 4.0*C_DDOT(noccA_*noccA_,wAA[0],1,zAA[0],1);

  free_block(wAA);
  free_block(zAA);
  free_block(zAS);

  double **B_p_SS = get_SS_ints(1);

  double *X = init_array((ndf_+3));

  C_DGEMV('t',nvirB_*nvirB_,(ndf_+3),1.0,B_p_SS[0],(ndf_+3),pSS[0],1,0.0,X,1);

  free_block(pSS);
  free_block(B_p_SS);

  double **xAB = block_matrix(noccA_,noccB_);

  C_DGEMV('n',noccA_*noccB_,(ndf_+3),1.0,&(B_p_AB[0][0]),(ndf_+3),X,1,0.0,
    &(xAB[0][0]),1);

  for (int a=0; a<noccA_; a++) {
    energy += 4.0*C_DDOT(noccB_,sAB_[a],1,xAB[a],1);
  }

  free_block(xAB);
  free_block(B_p_AB);

  double **xAA = block_matrix(noccA_,noccA_);
  double **yAA = block_matrix(noccA_,noccA_);

  C_DGEMV('n',noccA_*noccA_,(ndf_+3),1.0,&(B_p_AA[0][0]),(ndf_+3),X,1,
    0.0,&(xAA[0][0]),1);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,sAB_[0],nmoB_,sAB_[0],nmoB_,0.0,
    yAA[0],noccA_);

  energy -= 4.0*C_DDOT(noccA_*noccA_,xAA[0],1,yAA[0],1);

  free(X);
  free_block(xAA);
  free_block(yAA);
  free_block(B_p_AA);

  if (debug_) {
    fprintf(outfile,"\n    Exch12_k11u_1       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

double SAPT2::exch120_k11u_2()
{
  double energy=0.0;

  double **paa = block_matrix(aoccA_,aoccA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"pAA Density Matrix",(char *) paa[0],
    sizeof(double)*aoccA_*aoccA_);

  double **A_p_aB = get_AB_ints(1,foccA_,0);

  double **B_p_AB = get_AB_ints(2);
  double **B_p_aB = get_AB_ints(2,foccA_,0);

  double **A_p_aA = get_AA_ints(1,foccA_,0);
  double **A_p_Aa = get_AA_ints(1,0,foccA_);
  double **A_p_aa = get_AA_ints(1,foccA_,foccA_);

  double **B_p_BB = get_BB_ints(1);

  double **sAB = block_matrix(noccA_,noccB_);
  double **saB = block_matrix(aoccA_,noccB_);

  for (int a=0; a<noccA_; a++) {
    C_DCOPY(noccB_,sAB_[a],1,sAB[a],1);
  }

  for (int a=0; a<aoccA_; a++) {
    C_DCOPY(noccB_,sAB_[a+foccA_],1,saB[a],1);
  }

  double **C_p_aB = block_matrix(aoccA_*noccB_,ndf_+3);
  double **C_p_aA = block_matrix(aoccA_*noccA_,ndf_+3);
  double **xaa = block_matrix(aoccA_,aoccA_);
  double **xaA = block_matrix(aoccA_,noccA_);
  double **yaA = block_matrix(aoccA_,noccA_);
  double **xaB = block_matrix(aoccA_,noccB_);
  double **xBB = block_matrix(noccB_,noccB_);

  double *X = init_array(ndf_+3);

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_*(ndf_+3),1.0,A_p_aB[0],noccB_*(ndf_+3),
    B_p_aB[0],noccB_*(ndf_+3),0.0,xaa[0],aoccA_);

  energy += 2.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMV('t',noccA_*noccB_,ndf_+3,1.0,B_p_AB[0],ndf_+3,sAB[0],1,0.0,X,1);
  C_DGEMV('n',aoccA_*aoccA_,ndf_+3,1.0,A_p_aa[0],ndf_+3,X,1,0.0,xaa[0],1); 

  energy += 4.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,noccA_,1.0,sAB[0],noccB_,
      A_p_aA[a*noccA_],ndf_+3,0.0,C_p_aB[a*noccB_],ndf_+3);
  }

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_*(ndf_+3),1.0,C_p_aB[0],noccB_*(ndf_+3),
    B_p_aB[0],noccB_*(ndf_+3),0.0,xaa[0],aoccA_);

  energy -= 2.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMV('n',aoccA_*noccB_,ndf_+3,1.0,B_p_aB[0],ndf_+3,diagAA_,1,
    0.0,xaB[0],1);

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_,1.0,saB[0],noccB_,xaB[0],noccB_,
    0.0,xaa[0],aoccA_); 

  energy += 4.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  memset(&(xaB[0][0]),'\0',sizeof(double)*aoccA_*noccB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',aoccA_,noccB_,ndf_+3,1.0,A_p_Aa[a*aoccA_],ndf_+3,
      B_p_AB[a*noccB_],ndf_+3,1.0,xaB[0],noccB_);
  }

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_,1.0,saB[0],noccB_,xaB[0],noccB_,
    0.0,xaa[0],aoccA_);

  energy -= 2.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMV('n',aoccA_*noccB_,ndf_+3,1.0,A_p_aB[0],ndf_+3,diagBB_,1,
    0.0,xaB[0],1);

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_,1.0,xaB[0],noccB_,saB[0],noccB_,
    0.0,xaa[0],aoccA_);
  
  energy += 4.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMM('N','T',aoccA_,noccB_,noccB_*(ndf_+3),1.0,A_p_aB[0],noccB_*(ndf_+3),
    B_p_BB[0],noccB_*(ndf_+3),0.0,xaB[0],noccB_);

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_,1.0,xaB[0],noccB_,saB[0],noccB_,
    0.0,xaa[0],aoccA_);
  
  energy -= 2.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,sAB[0],noccB_,sAB[0],noccB_,
    0.0,xBB[0],noccB_);

  C_DGEMV('t',noccB_*noccB_,ndf_+3,1.0,B_p_BB[0],ndf_+3,xBB[0],1,0.0,X,1);

  C_DGEMV('n',aoccA_*aoccA_,ndf_+3,1.0,A_p_aa[0],ndf_+3,X,1,0.0,xaa[0],1);

  energy -= 4.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMV('n',noccB_*noccB_,ndf_+3,1.0,B_p_BB[0],ndf_+3,diagAA_,1,0.0,
    xBB[0],1);

  C_DGEMM('N','N',aoccA_,noccB_,noccB_,1.0,saB[0],noccB_,xBB[0],noccB_,
    0.0,xaB[0],noccB_);

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_,1.0,xaB[0],noccB_,saB[0],noccB_,
    0.0,xaa[0],aoccA_);

  energy -= 4.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMM('N','T',aoccA_,noccA_,noccB_,1.0,saB[0],noccB_,sAB[0],noccB_,
    0.0,xaA[0],noccA_);

  C_DGEMV('n',aoccA_*noccA_,ndf_+3,1.0,A_p_aA[0],ndf_+3,diagBB_,1,0.0,
    yaA[0],1);

  C_DGEMM('N','T',aoccA_,aoccA_,noccA_,1.0,xaA[0],noccA_,yaA[0],noccA_,
    0.0,xaa[0],aoccA_);

  energy -= 8.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  C_DGEMM('N','N',aoccA_,noccB_*(ndf_+3),noccB_,1.0,saB[0],noccB_,B_p_BB[0],
    noccB_*(ndf_+3),0.0,C_p_aB[0],noccB_*(ndf_+3));

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',noccA_,ndf_+3,noccB_,1.0,sAB[0],noccB_,C_p_aB[a*noccB_],
      ndf_+3,0.0,C_p_aA[a*noccA_],ndf_+3);
  }

  C_DGEMM('N','T',aoccA_,aoccA_,noccA_*(ndf_+3),1.0,C_p_aA[0],noccA_*(ndf_+3),
    A_p_aA[0],noccA_*(ndf_+3),0.0,xaa[0],aoccA_);

  energy += 4.0*C_DDOT(aoccA_*aoccA_,xaa[0],1,paa[0],1);

  free(X);

  free_block(C_p_aB);
  free_block(C_p_aA);
  free_block(paa);
  free_block(xaa);
  free_block(xaA);
  free_block(yaA);
  free_block(xaB);
  free_block(xBB);

  free_block(A_p_aB);
  free_block(B_p_AB);
  free_block(B_p_aB);
  free_block(A_p_aA);
  free_block(A_p_Aa);
  free_block(A_p_aa);
  free_block(B_p_BB);
  free_block(sAB);
  free_block(saB);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_2       = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::exch102_k11u_2()
{
  double energy=0.0;

  double **pbb = block_matrix(aoccB_,aoccB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"pBB Density Matrix",(char *) pbb[0],
    sizeof(double)*aoccB_*aoccB_);

  double **A_p_AB = get_AB_ints(1);
  double **A_p_Ab = get_AB_ints(1,0,foccB_);

  double **B_p_Ab = get_AB_ints(2,0,foccB_);

  double **A_p_AA = get_AA_ints(1);

  double **B_p_bB = get_BB_ints(1,foccB_,0);
  double **B_p_Bb = get_BB_ints(1,0,foccB_);
  double **B_p_bb = get_BB_ints(1,foccB_,foccB_);

  double **sAB = block_matrix(noccA_,noccB_);
  double **sAb = block_matrix(noccA_,aoccB_);

  for (int a=0; a<noccA_; a++) {
    C_DCOPY(noccB_,sAB_[a],1,sAB[a],1);
  }

  for (int a=0; a<noccA_; a++) {
    C_DCOPY(aoccB_,&(sAB_[a][foccB_]),1,sAb[a],1);
  }

  double **C_p_Ab = block_matrix(noccA_*aoccB_,ndf_+3);
  double **C_p_bB = block_matrix(aoccB_*noccB_,ndf_+3);
  double **xbb = block_matrix(aoccB_,aoccB_);
  double **xbB = block_matrix(aoccB_,noccB_);
  double **ybB = block_matrix(aoccB_,noccB_);
  double **xAb = block_matrix(noccA_,aoccB_);
  double **xAA = block_matrix(noccA_,noccA_);

  double *X = init_array(ndf_+3);

  memset(&(xbb[0][0]),'\0',sizeof(double)*aoccB_*aoccB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',aoccB_,aoccB_,ndf_+3,1.0,B_p_Ab[a*aoccB_],ndf_+3,
      A_p_Ab[a*aoccB_],ndf_+3,1.0,xbb[0],aoccB_);
  }

  energy += 2.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  C_DGEMV('t',noccA_*noccB_,ndf_+3,1.0,A_p_AB[0],ndf_+3,sAB[0],1,0.0,X,1);
  C_DGEMV('n',aoccB_*aoccB_,ndf_+3,1.0,B_p_bb[0],ndf_+3,X,1,0.0,xbb[0],1); 

  energy += 4.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',noccA_,ndf_+3,noccB_,1.0,sAB[0],noccB_,
      B_p_bB[b*noccB_],ndf_+3,0.0,C_p_Ab[b],aoccB_*(ndf_+3));
  }

  memset(&(xbb[0][0]),'\0',sizeof(double)*aoccB_*aoccB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',aoccB_,aoccB_,ndf_+3,1.0,C_p_Ab[a*aoccB_],ndf_+3,
      A_p_Ab[a*aoccB_],ndf_+3,1.0,xbb[0],aoccB_);
  }

  energy -= 2.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  C_DGEMV('n',noccA_*aoccB_,ndf_+3,1.0,A_p_Ab[0],ndf_+3,diagBB_,1,
    0.0,xAb[0],1);

  C_DGEMM('T','N',aoccB_,aoccB_,noccA_,1.0,sAb[0],aoccB_,xAb[0],aoccB_,
    0.0,xbb[0],aoccB_); 

  energy += 4.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  memset(&(xAb[0][0]),'\0',sizeof(double)*noccA_*aoccB_);

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('N','T',noccA_,aoccB_,ndf_+3,1.0,A_p_AB[b],noccB_*(ndf_+3),
      B_p_Bb[b*aoccB_],ndf_+3,1.0,xAb[0],aoccB_);
  }

  C_DGEMM('T','N',aoccB_,aoccB_,noccA_,1.0,sAb[0],aoccB_,xAb[0],aoccB_,
    0.0,xbb[0],aoccB_);

  energy -= 2.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  C_DGEMV('n',noccA_*aoccB_,ndf_+3,1.0,B_p_Ab[0],ndf_+3,diagAA_,1,
    0.0,xAb[0],1);

  C_DGEMM('T','N',aoccB_,aoccB_,noccA_,1.0,xAb[0],aoccB_,sAb[0],aoccB_,
    0.0,xbb[0],aoccB_);
  
  energy += 4.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  memset(&(xAb[0][0]),'\0',sizeof(double)*noccA_*aoccB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',noccA_,aoccB_,ndf_+3,1.0,A_p_AA[a*noccA_],ndf_+3,
      B_p_Ab[a*aoccB_],ndf_+3,1.0,xAb[0],aoccB_);
  }

  C_DGEMM('T','N',aoccB_,aoccB_,noccA_,1.0,xAb[0],aoccB_,sAb[0],aoccB_,
    0.0,xbb[0],aoccB_);
  
  energy -= 2.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,sAB[0],noccB_,sAB[0],noccB_,
    0.0,xAA[0],noccA_);

  C_DGEMV('t',noccA_*noccA_,ndf_+3,1.0,A_p_AA[0],ndf_+3,xAA[0],1,0.0,X,1);

  C_DGEMV('n',aoccB_*aoccB_,ndf_+3,1.0,B_p_bb[0],ndf_+3,X,1,0.0,xbb[0],1);

  energy -= 4.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  C_DGEMV('n',noccA_*noccA_,ndf_+3,1.0,A_p_AA[0],ndf_+3,diagBB_,1,0.0,
    xAA[0],1);

  C_DGEMM('N','N',noccA_,aoccB_,noccA_,1.0,xAA[0],noccA_,sAb[0],aoccB_,
    0.0,xAb[0],aoccB_);

  C_DGEMM('T','N',aoccB_,aoccB_,noccA_,1.0,xAb[0],aoccB_,sAb[0],aoccB_,
    0.0,xbb[0],aoccB_);

  energy -= 4.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  C_DGEMM('T','N',aoccB_,noccB_,noccA_,1.0,sAb[0],aoccB_,sAB[0],noccB_,
    0.0,xbB[0],noccB_);

  C_DGEMV('n',aoccB_*noccB_,ndf_+3,1.0,B_p_bB[0],ndf_+3,diagAA_,1,0.0,
    ybB[0],1);

  C_DGEMM('N','T',aoccB_,aoccB_,noccB_,1.0,xbB[0],noccB_,ybB[0],noccB_,
    0.0,xbb[0],aoccB_);

  energy -= 8.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  C_DGEMM('T','N',aoccB_,noccA_*(ndf_+3),noccA_,1.0,sAb[0],aoccB_,A_p_AA[0],
    noccA_*(ndf_+3),0.0,C_p_Ab[0],noccA_*(ndf_+3));

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('T','N',noccB_,ndf_+3,noccA_,1.0,sAB[0],noccB_,C_p_Ab[b*noccA_],
      ndf_+3,0.0,C_p_bB[b*noccB_],ndf_+3);
  }

  C_DGEMM('N','T',aoccB_,aoccB_,noccB_*(ndf_+3),1.0,C_p_bB[0],noccB_*(ndf_+3),
    B_p_bB[0],noccB_*(ndf_+3),0.0,xbb[0],aoccB_);

  energy += 4.0*C_DDOT(aoccB_*aoccB_,xbb[0],1,pbb[0],1);

  free(X);

  free_block(C_p_Ab);
  free_block(C_p_bB);
  free_block(pbb);
  free_block(xbb);
  free_block(xbB);
  free_block(ybB);
  free_block(xAb);
  free_block(xAA);

  free_block(B_p_Ab);
  free_block(A_p_AB);
  free_block(A_p_Ab);
  free_block(B_p_bB);
  free_block(B_p_Bb);
  free_block(B_p_bb);
  free_block(A_p_AA);
  free_block(sAB);
  free_block(sAb);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_2       = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::exch120_k11u_3()
{
  double energy=0.0;

  double **temp_thetaARAR = block_matrix(aoccA_*nvirA_,aoccA_*nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARAR Amplitudes",
    (char *) temp_thetaARAR[0],sizeof(double)*aoccA_*nvirA_*aoccA_*nvirA_);
  antisym(temp_thetaARAR,aoccA_,nvirA_);

  double **thetaRRAA = block_matrix(nvirA_*nvirA_,aoccA_*aoccA_);

  for(int a1=0,a1r1=0; a1<aoccA_; a1++) {
    for(int r1=0; r1<nvirA_; r1++,a1r1++) {
      for(int a2=0,a2r2=0; a2<aoccA_; a2++) {
        for(int r2=0; r2<nvirA_; r2++,a2r2++) {
          int a1a2 = a1*aoccA_+a2;
          int r1r2 = r1*nvirA_+r2;
          thetaRRAA[r1r2][a1a2] = temp_thetaARAR[a1r1][a2r2];
  }}}}

  free_block(temp_thetaARAR);

  double **thetaRBAA = block_matrix(nvirA_*noccB_,aoccA_*aoccA_);

  for(int r1=0; r1<nvirA_; r1++) {
    C_DGEMM('T','N',noccB_,aoccA_*aoccA_,nvirA_,1.0,&(sAB_[noccA_][0]),
      nmoB_,&(thetaRRAA[r1*nvirA_][0]),aoccA_*aoccA_,0.0,
      &(thetaRBAA[r1*noccB_][0]),aoccA_*aoccA_);
  }

  free_block(thetaRRAA);

  double **temp_tARAR = block_matrix(aoccA_*nvirA_,aoccA_*nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARAR Amplitudes",
    (char *) temp_tARAR[0],sizeof(double)*aoccA_*nvirA_*aoccA_*nvirA_);

  double **tRRAA = block_matrix(nvirA_*nvirA_,aoccA_*aoccA_);

  for(int a1=0,a1r1=0; a1<aoccA_; a1++) {
    for(int r1=0; r1<nvirA_; r1++,a1r1++) {
      for(int a2=0,a2r2=0; a2<aoccA_; a2++) {
        for(int r2=0; r2<nvirA_; r2++,a2r2++) {
          int a1a2 = a1*aoccA_+a2;
          int r1r2 = r1*nvirA_+r2;
          tRRAA[r1r2][a1a2] = temp_tARAR[a1r1][a2r2];
  }}}}

  free_block(temp_tARAR);

  double **B_p_RB = get_RB_ints(1);
  double **B_p_RR = get_RR_ints(1);

  double *yRB = init_array(nvirA_*noccB_);
  double **zRB = block_matrix(nvirA_,nvirA_*
    noccB_);

  for (int r1=0; r1<nvirA_; r1++) {
    C_DGEMM('N','T',r1+1,nvirA_*noccB_,(ndf_+3),1.0,B_p_RR[r1*nvirA_],
      (ndf_+3),B_p_RB[0],(ndf_+3),0.0,zRB[0],nvirA_*noccB_);
    for (int r2=0; r2<=r1; r2++) {
      int r1r2 = r1*nvirA_+r2;
      C_DGEMM('N','T',nvirA_,noccB_,aoccA_*aoccA_,1.0,tRRAA[r2*nvirA_],
        aoccA_*aoccA_,thetaRBAA[r1*noccB_],aoccA_*aoccA_,0.0,yRB,noccB_);
      if (r1 != r2)
      C_DGEMM('N','T',nvirA_,noccB_,aoccA_*aoccA_,1.0,tRRAA[r1*nvirA_],
        aoccA_*aoccA_,thetaRBAA[r2*noccB_],aoccA_*aoccA_,1.0,yRB,noccB_);
      energy += 2.0*C_DDOT(nvirA_*noccB_,yRB,1,zRB[r2],1);
  }}

  free(yRB);
  free_block(zRB);
  free_block(B_p_RB);

  double **tRBAA = block_matrix(nvirA_*noccB_,aoccA_*aoccA_);

  for(int r1=0; r1<nvirA_; r1++) {
    C_DGEMM('T','N',noccB_,aoccA_*aoccA_,nvirA_,1.0,&(sAB_[noccA_][0]),
      nmoB_,&(tRRAA[r1*nvirA_][0]),aoccA_*aoccA_,0.0,&(tRBAA[r1*noccB_][0]),
      aoccA_*aoccA_);
  }

  free_block(tRRAA);

  double **xRR = block_matrix(nvirA_,nvirA_);
  double **yRR = block_matrix(nvirA_,nvirA_);

  C_DGEMM('N','T',nvirA_,nvirA_,aoccA_*aoccA_*noccB_,1.0,&(tRBAA[0][0]),
    aoccA_*aoccA_*noccB_,&(thetaRBAA[0][0]),aoccA_*aoccA_*noccB_,0.0,
    &(xRR[0][0]),nvirA_);

  C_DGEMV('n',nvirA_*nvirA_,(ndf_+3),1.0,B_p_RR[0],(ndf_+3),diagBB_,1,0.0,
    &(yRR[0][0]),1);

  energy += 4.0*C_DDOT(nvirA_*nvirA_,xRR[0],1,yRR[0],1);

  free_block(xRR);
  free_block(yRR);

  double **B_p_BB = get_BB_ints(1);

  double *yBB = init_array(noccB_*noccB_);
  double **zBB = block_matrix(nvirA_,noccB_*noccB_);

  for(int r1=0; r1<nvirA_; r1++) {
    C_DGEMM('N','T',r1+1,noccB_*noccB_,(ndf_+3),1.0,B_p_RR[r1*nvirA_],
      (ndf_+3),B_p_BB[0],(ndf_+3),0.0,zBB[0],noccB_*noccB_);
    for(int r2=0; r2<=r1; r2++) {
      int r1r2 = r1*nvirA_ + r2;
      C_DGEMM('N','T',noccB_,noccB_,aoccA_*aoccA_,1.0,&(tRBAA[r2*noccB_][0]),
        aoccA_*aoccA_,&(thetaRBAA[r1*noccB_][0]),aoccA_*aoccA_,0.0,yBB,noccB_);
      if (r1 != r2)
      C_DGEMM('N','T',noccB_,noccB_,aoccA_*aoccA_,1.0,&(tRBAA[r1*noccB_][0]),
        aoccA_*aoccA_,&(thetaRBAA[r2*noccB_][0]),aoccA_*aoccA_,1.0,yBB,noccB_);
      energy -= 2.0*C_DDOT(noccB_*noccB_,yBB,1,zBB[r2],1);
  }}

  free_block(tRBAA);
  free_block(thetaRBAA);
  free_block(B_p_BB);
  free_block(B_p_RR);
  free(yBB);
  free_block(zBB);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_3       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

double SAPT2::exch102_k11u_3()
{
  double energy=0.0;

  double **temp_thetaBSBS = block_matrix(aoccB_*nvirB_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tBSBS Amplitudes",
    (char *) temp_thetaBSBS[0],sizeof(double)*aoccB_*nvirB_*aoccB_*nvirB_);
  antisym(temp_thetaBSBS,aoccB_,nvirB_);

  double **thetaSSBB = block_matrix(nvirB_*nvirB_,aoccB_*aoccB_);

  for(int b1=0,b1s1=0; b1<aoccB_; b1++) {
    for(int s1=0; s1<nvirB_; s1++,b1s1++) {
      for(int b2=0,b2s2=0; b2<aoccB_; b2++) {
        for(int s2=0; s2<nvirB_; s2++,b2s2++) {
          int b1b2 = b1*aoccB_+b2;
          int s1s2 = s1*nvirB_+s2;
          thetaSSBB[s1s2][b1b2] = temp_thetaBSBS[b1s1][b2s2];
  }}}}

  free_block(temp_thetaBSBS);

  double **thetaSABB = block_matrix(nvirB_*noccA_,aoccB_*aoccB_);

  for(int s1=0; s1<nvirB_; s1++) {
    C_DGEMM('N','N',noccA_,aoccB_*aoccB_,nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
      &(thetaSSBB[s1*nvirB_][0]),aoccB_*aoccB_,0.0,
      &(thetaSABB[s1*noccA_][0]),aoccB_*aoccB_);
  }

  free_block(thetaSSBB);

  double **temp_tBSBS = block_matrix(aoccB_*nvirB_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tBSBS Amplitudes",
    (char *) temp_tBSBS[0],sizeof(double)*aoccB_*nvirB_*aoccB_*nvirB_);

  double **tSSBB = block_matrix(nvirB_*nvirB_,aoccB_*aoccB_);

  for(int b1=0,b1s1=0; b1<aoccB_; b1++) {
    for(int s1=0; s1<nvirB_; s1++,b1s1++) {
      for(int b2=0,b2s2=0; b2<aoccB_; b2++) {
        for(int s2=0; s2<nvirB_; s2++,b2s2++) {
          int b1b2 = b1*aoccB_+b2;
          int s1s2 = s1*nvirB_+s2;
          tSSBB[s1s2][b1b2] = temp_tBSBS[b1s1][b2s2];
  }}}}

  free_block(temp_tBSBS);

  double **B_AS_p = get_AS_ints(1);
  double **B_p_SA = block_matrix(noccA_*nvirB_,(ndf_+3));

  for(int a=0,as=0; a<noccA_; a++) {
    for(int s=0; s<nvirB_; s++,as++) {
      int sa = s*noccA_ + a;
      for(int P=0; P<(ndf_+3); P++) {
        B_p_SA[sa][P] = B_AS_p[as][P];
  }}}

  free_block(B_AS_p);

  double **B_p_SS = get_SS_ints(1);

  double *ySA = init_array(nvirB_*noccA_);
  double **zSA = block_matrix(nvirB_,nvirB_*
    noccA_);

  for (int s1=0; s1<nvirB_; s1++) {
    C_DGEMM('N','T',s1+1,nvirB_*noccA_,(ndf_+3),1.0,B_p_SS[s1*nvirB_],
      (ndf_+3),B_p_SA[0],(ndf_+3),0.0,zSA[0],nvirB_*noccA_);
    for (int s2=0; s2<=s1; s2++) {
      int s1s2 = s1*nvirA_+s2;
      C_DGEMM('N','T',nvirB_,noccA_,aoccB_*aoccB_,1.0,tSSBB[s2*nvirB_],
        aoccB_*aoccB_,thetaSABB[s1*noccA_],aoccB_*aoccB_,0.0,ySA,noccA_);
      if (s1 != s2)
      C_DGEMM('N','T',nvirB_,noccA_,aoccB_*aoccB_,1.0,tSSBB[s1*nvirB_],
        aoccB_*aoccB_,thetaSABB[s2*noccA_],aoccB_*aoccB_,1.0,ySA,noccA_);
      energy += 2.0*C_DDOT(nvirB_*noccA_,ySA,1,zSA[s2],1);
  }}

  free(ySA);
  free_block(zSA);
  free_block(B_p_SA);

  double **tSABB = block_matrix(nvirB_*noccA_,aoccB_*aoccB_);

  for(int s1=0; s1<nvirB_; s1++) {
    C_DGEMM('N','N',noccA_,aoccB_*aoccB_,nvirB_,1.0,&(sAB_[0][noccB_]),
      nmoB_,&(tSSBB[s1*nvirB_][0]),aoccB_*aoccB_,0.0,
      &(tSABB[s1*noccA_][0]),aoccB_*aoccB_);
  }

  free_block(tSSBB);

  double **xSS = block_matrix(nvirB_,nvirB_);
  double **ySS = block_matrix(nvirB_,nvirB_);

  C_DGEMM('N','T',nvirB_,nvirB_,aoccB_*aoccB_*noccA_,1.0,&(tSABB[0][0]),
    aoccB_*aoccB_*noccA_,&(thetaSABB[0][0]),aoccB_*aoccB_*noccA_,0.0,
    &(xSS[0][0]),nvirB_);

  C_DGEMV('n',nvirB_*nvirB_,(ndf_+3),1.0,B_p_SS[0],(ndf_+3),diagAA_,1,
    0.0,&(ySS[0][0]),1);

  energy += 4.0*C_DDOT(nvirB_*nvirB_,xSS[0],1,ySS[0],1);

  free_block(xSS);
  free_block(ySS);

  double **B_p_AA = get_AA_ints(1);

  double *yAA = init_array(noccA_*noccA_);
  double **zAA = block_matrix(nvirB_,noccA_*noccA_);

  for(int s1=0; s1<nvirB_; s1++) {
    C_DGEMM('N','T',s1+1,noccA_*noccA_,(ndf_+3),1.0,B_p_SS[s1*nvirB_],
      (ndf_+3),B_p_AA[0],(ndf_+3),0.0,zAA[0],noccA_*noccA_);
    for(int s2=0; s2<=s1; s2++) {
      int s1s2 = s1*nvirB_ + s2;
      C_DGEMM('N','T',noccA_,noccA_,aoccB_*aoccB_,1.0,&(tSABB[s2*noccA_][0]),
        aoccB_*aoccB_,&(thetaSABB[s1*noccA_][0]),aoccB_*aoccB_,0.0,yAA,noccA_);
      if (s1 != s2)
      C_DGEMM('N','T',noccA_,noccA_,aoccB_*aoccB_,1.0,&(tSABB[s1*noccA_][0]),
        aoccB_*aoccB_,&(thetaSABB[s2*noccA_][0]),aoccB_*aoccB_,1.0,yAA,noccA_);
      energy -= 2.0*C_DDOT(noccA_*noccA_,yAA,1,zAA[s2],1);
  }}

  free_block(tSABB);
  free_block(thetaSABB);
  free(yAA);
  free_block(zAA);
  free_block(B_p_AA);
  free_block(B_p_SS);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_3       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

double SAPT2::exch120_k11u_4()
{
  double energy=0.0;

  double *tARAR = init_array((long int) aoccA_*nvirA_*aoccA_*nvirA_);
  double *thetaARAR = init_array((long int) aoccA_*nvirA_*aoccA_*nvirA_);

  psio_->read_entry(PSIF_SAPT_AMPS,"tARAR Amplitudes",(char *) &(tARAR[0]),
    sizeof(double)*aoccA_*nvirA_*aoccA_*nvirA_);
  C_DCOPY(aoccA_*nvirA_*aoccA_*nvirA_,tARAR,1,thetaARAR,1);
  antisym(thetaARAR,aoccA_,nvirA_);

  ijkl_to_ikjl(tARAR,aoccA_,nvirA_,aoccA_,nvirA_);
  ijkl_to_ikjl(thetaARAR,aoccA_,nvirA_,aoccA_,nvirA_);

  double *AAAA = init_array((long int) aoccA_*aoccA_*aoccA_*aoccA_);

  C_DGEMM('N','T',aoccA_*aoccA_,aoccA_*aoccA_,nvirA_*nvirA_,1.0,thetaARAR,
    nvirA_*nvirA_,tARAR,nvirA_*nvirA_,0.0,AAAA,aoccA_*aoccA_);

  free(tARAR);
  free(thetaARAR);

  ijkl_to_ikjl(AAAA,aoccA_,aoccA_,aoccA_,aoccA_);

  double **B_p_AA = get_AA_ints(1,foccA_,foccA_);
  double **C_p_AA = block_matrix(aoccA_*aoccA_,ndf_+3);

  C_DGEMM('N','N',aoccA_*aoccA_,ndf_+3,aoccA_*aoccA_,1.0,AAAA,aoccA_*aoccA_,
    &(B_p_AA[0][0]),ndf_+3,0.0,&(C_p_AA[0][0]),ndf_+3);

  free(AAAA);
  free_block(B_p_AA);

  double **B_p_AB = get_AB_ints(2,foccA_,0);
  double **D_p_AA = block_matrix(aoccA_*aoccA_,ndf_+3);

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
      &(B_p_AB[a*noccB_][0]),ndf_+3,0.0,&(D_p_AA[a*aoccA_][0]),(ndf_+3));
  }

  energy += 2.0*C_DDOT((long int) aoccA_*aoccA_*(ndf_+3),C_p_AA[0],1,
    D_p_AA[0],1);

  free_block(B_p_AB);
  free_block(D_p_AA);

  double *X = init_array(ndf_+3);
  double **sAA = block_matrix(aoccA_,aoccA_);

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[foccA_][0]),nmoB_,0.0,&(sAA[0][0]),aoccA_);

  C_DGEMV('t',aoccA_*aoccA_,ndf_+3,1.0,C_p_AA[0],ndf_+3,sAA[0],1,0.0,X,1);

  energy += 4.0*C_DDOT(ndf_+3,X,1,diagBB_,1);

  free(X);
  free_block(sAA);

  double **B_p_BB = get_BB_ints(1);
  double **C_p_AB = block_matrix(aoccA_*noccB_,ndf_+3);

  C_DGEMM('N','N',aoccA_,noccB_*(ndf_+3),noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(B_p_BB[0][0]),noccB_*(ndf_+3),0.0,&(C_p_AB[0][0]),noccB_*(ndf_+3));

  free_block(B_p_BB);

  double **E_p_AA = block_matrix(aoccA_*aoccA_,ndf_+3);

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
      &(C_p_AB[a*noccB_][0]),ndf_+3,0.0,&(E_p_AA[a*aoccA_][0]),ndf_+3);
  }

  energy -= 2.0*C_DDOT((long int) aoccA_*aoccA_*(ndf_+3),C_p_AA[0],1,
    E_p_AA[0],1);

  free_block(C_p_AB);
  free_block(C_p_AA);
  free_block(E_p_AA);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_4       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

double SAPT2::exch102_k11u_4()
{
  double energy=0.0;

  double *tBSBS = init_array((long int) aoccB_*nvirB_*aoccB_*nvirB_);
  double *thetaBSBS = init_array((long int) aoccB_*nvirB_*aoccB_*nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"tBSBS Amplitudes",(char *) &(tBSBS[0]),
    sizeof(double)*aoccB_*nvirB_*aoccB_*nvirB_);
  C_DCOPY(aoccB_*nvirB_*aoccB_*nvirB_,tBSBS,1,thetaBSBS,1);
  antisym(thetaBSBS,aoccB_,nvirB_);

  ijkl_to_ikjl(tBSBS,aoccB_,nvirB_,aoccB_,nvirB_);
  ijkl_to_ikjl(thetaBSBS,aoccB_,nvirB_,aoccB_,nvirB_);

  double *BBBB = init_array((long int) aoccB_*aoccB_*aoccB_*aoccB_);

  C_DGEMM('N','T',aoccB_*aoccB_,aoccB_*aoccB_,nvirB_*nvirB_,1.0,thetaBSBS,
    nvirB_*nvirB_,tBSBS,nvirB_*nvirB_,0.0,BBBB,aoccB_*aoccB_);

  free(tBSBS);
  free(thetaBSBS);

  ijkl_to_ikjl(BBBB,aoccB_,aoccB_,aoccB_,aoccB_);

  double **B_p_BB = get_BB_ints(1,foccB_,foccB_);
  double **C_p_BB = block_matrix(aoccB_*aoccB_,ndf_+3);

  C_DGEMM('N','N',aoccB_*aoccB_,ndf_+3,aoccB_*aoccB_,1.0,BBBB,aoccB_*aoccB_,
    &(B_p_BB[0][0]),ndf_+3,0.0,&(C_p_BB[0][0]),ndf_+3);

  free(BBBB);
  free_block(B_p_BB);

  double **B_p_AB = get_AB_ints(1,0,foccB_);
  double **D_p_BB = block_matrix(aoccB_*aoccB_,ndf_+3);

  C_DGEMM('T','N',aoccB_,aoccB_*(ndf_+3),noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(B_p_AB[0][0]),aoccB_*(ndf_+3),0.0,&(D_p_BB[0][0]),aoccB_*(ndf_+3));

  energy += 2.0*C_DDOT((long int) aoccB_*aoccB_*(ndf_+3),C_p_BB[0],1,
    D_p_BB[0],1);

  free_block(B_p_AB);
  free_block(D_p_BB);

  double *X = init_array(ndf_+3);
  double **sBB = block_matrix(aoccB_,aoccB_);

  C_DGEMM('T','N',aoccB_,aoccB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][foccB_]),nmoB_,0.0,&(sBB[0][0]),aoccB_);

  C_DGEMV('t',aoccB_*aoccB_,(ndf_+3),1.0,C_p_BB[0],(ndf_+3),sBB[0],1,0.0,X,1);

  energy += 4.0*C_DDOT((ndf_+3),X,1,diagAA_,1);

  free(X);
  free_block(sBB);

  double **B_p_AA = get_AA_ints(1);
  double **C_p_BA = block_matrix(noccA_*aoccB_,ndf_+3);

  C_DGEMM('T','N',aoccB_,noccA_*(ndf_+3),noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(B_p_AA[0][0]),noccA_*(ndf_+3),0.0,&(C_p_BA[0][0]),noccA_*(ndf_+3));

  free_block(B_p_AA);

  double **E_p_BB = block_matrix(aoccB_*aoccB_,ndf_+3);

  for(int b=0; b<aoccB_; b++) {
    C_DGEMM('T','N',aoccB_,ndf_+3,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
      &(C_p_BA[b*noccA_][0]),ndf_+3,0.0,&(E_p_BB[b*aoccB_][0]),ndf_+3);
  }

  energy -= 2.0*C_DDOT((long int) aoccB_*aoccB_*(ndf_+3),C_p_BB[0],1,
    E_p_BB[0],1);

  free_block(C_p_BA);
  free_block(C_p_BB);
  free_block(E_p_BB);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_4       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

double SAPT2::exch120_k11u_5()
{
  double energy=0.0;

  double **theta_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"Theta AR Intermediates",
    (char *) theta_p_AR[0],sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  double **thetaARAR = block_matrix(aoccA_*nvirA_,aoccA_*nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARAR Amplitudes",
    (char *) thetaARAR[0],sizeof(double)*aoccA_*nvirA_*aoccA_*nvirA_);
  antisym(thetaARAR,aoccA_,nvirA_);

  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);

  C_DGEMM('N','N',aoccA_*nvirA_,ndf_+3,aoccA_*nvirA_,1.0,&(thetaARAR[0][0]),
    aoccA_*nvirA_,&(theta_p_AR[0][0]),ndf_+3,0.0,&(T_p_AR[0][0]),ndf_+3);

  free_block(theta_p_AR);
  free_block(thetaARAR);

  double **C_p_BR = block_matrix(noccB_*nvirA_,ndf_+3);

  C_DGEMM('T','N',noccB_,nvirA_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(T_p_AR[0][0]),nvirA_*(ndf_+3),0.0,&(C_p_BR[0][0]),nvirA_*(ndf_+3));

  double **B_p_RB = get_RB_ints(1);

  for (int r=0,rb=0; r<nvirA_; r++){
    for (int b=0; b<noccB_; b++,rb++){
      int br = b*nvirA_+r;
      energy += C_DDOT((ndf_+3),C_p_BR[br],1,B_p_RB[rb],1);
  }}

  free_block(B_p_RB);
  free_block(C_p_BR);

  double **C_p_AB = block_matrix(aoccA_*noccB_,ndf_+3);

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][0]),nmoB_,
      &(T_p_AR[a*nvirA_][0]),ndf_+3,0.0,&(C_p_AB[a*noccB_][0]),ndf_+3);
  }

  double **B_p_AB = get_AB_ints(2,foccA_,0);

  energy += C_DDOT((long int) aoccA_*noccB_*(ndf_+3),B_p_AB[0],1,C_p_AB[0],1);

  free_block(B_p_AB);

  double **C_p_BB = block_matrix(noccB_*noccB_,ndf_+3);

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(C_p_AB[0][0]),noccB_*(ndf_+3),0.0,&(C_p_BB[0][0]),noccB_*(ndf_+3));

  free_block(C_p_AB);

  double **B_p_BB = get_BB_ints(1);

  energy -= 2.0*C_DDOT((long int) noccB_*noccB_*(ndf_+3),B_p_BB[0],1,
    C_p_BB[0],1);

  free_block(B_p_BB);
  free_block(C_p_BB);

  double **xAR = block_matrix(aoccA_,nvirA_);
  double **yAR = block_matrix(aoccA_,nvirA_);

  C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(xAR[0][0]),nvirA_);

  C_DGEMV('n',aoccA_*nvirA_,ndf_+3,1.0,T_p_AR[0],ndf_+3,diagBB_,1,0.0,
    &(yAR[0][0]),1);

  energy += 4.0*C_DDOT(aoccA_*nvirA_,xAR[0],1,yAR[0],1);

  free_block(xAR);
  free_block(yAR);
  free_block(T_p_AR);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_5       = %18.12lf H\n",-2.0*energy);
    fflush(outfile);
  }

  return(-2.0*energy);
}

double SAPT2::exch102_k11u_5()
{
  double energy=0.0;

  double **theta_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"Theta BS Intermediates",
    (char *) theta_p_BS[0],sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  double **thetaBSBS = block_matrix(aoccB_*nvirB_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tBSBS Amplitudes",
    (char *) thetaBSBS[0],sizeof(double)*aoccB_*nvirB_*aoccB_*nvirB_);
  antisym(thetaBSBS,aoccB_,nvirB_);

  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('N','N',aoccB_*nvirB_,ndf_+3,aoccB_*nvirB_,1.0,&(thetaBSBS[0][0]),
    aoccB_*nvirB_,&(theta_p_BS[0][0]),ndf_+3,0.0,&(T_p_BS[0][0]),ndf_+3);

  free_block(theta_p_BS);
  free_block(thetaBSBS);

  double **C_p_AS = block_matrix(noccA_*nvirB_,ndf_+3);

  C_DGEMM('N','N',noccA_,nvirB_*(ndf_+3),aoccB_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(T_p_BS[0][0]),nvirB_*(ndf_+3),0.0,&(C_p_AS[0][0]),nvirB_*(ndf_+3));

  double **B_p_AS = get_AS_ints(1);

  energy += C_DDOT((long int) noccA_*nvirB_*(ndf_+3),C_p_AS[0],1,B_p_AS[0],1);

  free_block(B_p_AS);
  free_block(C_p_AS);

  double **C_p_BA = block_matrix(aoccB_*noccA_,ndf_+3);

  for(int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',noccA_,(ndf_+3),nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
      &(T_p_BS[b*nvirB_][0]),(ndf_+3),0.0,&(C_p_BA[b*noccA_][0]),(ndf_+3));
  }

  double **B_p_AB = get_AB_ints(1,0,foccB_);

  for(int a=0,ab=0; a<noccA_; a++) {
    for(int b=0; b<aoccB_; b++,ab++) {
      int ba = b*noccA_+a;
      energy += C_DDOT(ndf_+3,B_p_AB[ab],1,C_p_BA[ba],1);
  }}

  free_block(B_p_AB);

  double **C_p_AA = block_matrix(noccA_*noccA_,ndf_+3);

  C_DGEMM('N','N',noccA_,noccA_*(ndf_+3),aoccB_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(C_p_BA[0][0]),noccA_*(ndf_+3),0.0,&(C_p_AA[0][0]),noccA_*(ndf_+3));

  free_block(C_p_BA);

  double **B_p_AA = get_AA_ints(1);

  energy -= 2.0*C_DDOT((long int) noccA_*noccA_*(ndf_+3),B_p_AA[0],1,
    C_p_AA[0],1);

  free_block(B_p_AA);
  free_block(C_p_AA);

  double **xBS = block_matrix(aoccB_,nvirB_);
  double **yBS = block_matrix(aoccB_,nvirB_);

  C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,&(xBS[0][0]),nvirB_);

  C_DGEMV('n',aoccB_*nvirB_,ndf_+3,1.0,T_p_BS[0],ndf_+3,diagAA_,1,0.0,
    &(yBS[0][0]),1);

  energy += 4.0*C_DDOT(aoccB_*nvirB_,xBS[0],1,yBS[0],1);

  free_block(xBS);
  free_block(yBS);
  free_block(T_p_BS);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_5       = %18.12lf H\n",-2.0*energy);
    fflush(outfile);
  }

  return(-2.0*energy);
}

double SAPT2::exch120_k11u_6()
{
  double energy=0.0;

  double *T_ARAR = init_array((long int) aoccA_*nvirA_*aoccA_*nvirA_);

  double *tARAR = init_array((long int) aoccA_*nvirA_*aoccA_*nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARAR Amplitudes",
    (char *) tARAR,sizeof(double)*aoccA_*nvirA_*aoccA_*nvirA_);

  C_DGEMM('N','T',aoccA_*nvirA_,aoccA_*nvirA_,aoccA_*nvirA_,3.0,tARAR,
    aoccA_*nvirA_,tARAR,aoccA_*nvirA_,0.0,T_ARAR,aoccA_*nvirA_);

  antisym(tARAR,aoccA_,nvirA_);
  OVOpVp_to_OVpOpV(tARAR,aoccA_,nvirA_);

  C_DGEMM('N','T',aoccA_*nvirA_,aoccA_*nvirA_,aoccA_*nvirA_,1.0,tARAR,
    aoccA_*nvirA_,tARAR,aoccA_*nvirA_,1.0,T_ARAR,aoccA_*nvirA_);

  ijkl_to_ikjl(T_ARAR,aoccA_,nvirA_,aoccA_,nvirA_);

  free(tARAR);

  double **B_p_RR = get_RR_ints(1);
  double **T_p_AA = block_matrix(aoccA_*aoccA_,ndf_+3);

  C_DGEMM('N','N',aoccA_*aoccA_,ndf_+3,nvirA_*nvirA_,1.0,T_ARAR,nvirA_*nvirA_,
    &(B_p_RR[0][0]),ndf_+3,0.0,&(T_p_AA[0][0]),ndf_+3);

  free_block(B_p_RR);

  double **B_p_AA = get_AA_ints(1,foccA_,foccA_);
  double **T_p_RR = block_matrix(nvirA_*nvirA_,ndf_+3);

  C_DGEMM('T','N',nvirA_*nvirA_,ndf_+3,aoccA_*aoccA_,1.0,T_ARAR,nvirA_*nvirA_,
    &(B_p_AA[0][0]),ndf_+3,0.0,&(T_p_RR[0][0]),ndf_+3);

  free(T_ARAR);
  free_block(B_p_AA);

  double **B_p_AB = get_AB_ints(2,foccA_,0);
  double **C_p_BA = block_matrix(noccB_*aoccA_,ndf_+3);

  C_DGEMM('T','N',noccB_,aoccA_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(T_p_AA[0][0]),aoccA_*(ndf_+3),0.0,&(C_p_BA[0][0]),aoccA_*(ndf_+3));

  for(int a=0, ab=0; a<aoccA_; a++) {
    for(int b=0; b<noccB_; b++, ab++) {
      int ba = b*aoccA_+a;
      energy -= C_DDOT(ndf_+3,&(C_p_BA[ba][0]),1,&(B_p_AB[ab][0]),1);
  }}

  free_block(B_p_AB);

  double **B_p_BB = get_BB_ints(1);
  double **C_p_BB = block_matrix(noccB_*noccB_,ndf_+3);

  for(int b=0; b<noccB_; b++) {
    C_DGEMM('T','N',noccB_,ndf_+3,aoccA_,1.0,&(sAB_[foccA_][0]),nmoB_,
      &(C_p_BA[b*aoccA_][0]),ndf_+3,0.0,&(C_p_BB[b*noccB_][0]),ndf_+3);
  }

  energy += C_DDOT((long int) noccB_*noccB_*(ndf_+3),&(C_p_BB[0][0]),1,
    &(B_p_BB[0][0]),1);

  free_block(C_p_BB);
  free_block(C_p_BA);

  double **B_p_RB = get_RB_ints(1);
  double **C_p_BR = block_matrix(noccB_*nvirA_,ndf_+3);

  C_DGEMM('T','N',noccB_,nvirA_*(ndf_+3),nvirA_,1.0,&(sAB_[noccA_][0]),
    nmoB_,&(T_p_RR[0][0]),nvirA_*(ndf_+3),0.0,&(C_p_BR[0][0]),nvirA_*(ndf_+3));

  for(int r=0, rb=0; r<nvirA_; r++) {
    for(int b=0; b<noccB_; b++, rb++) {
      int br = b*nvirA_+r;
      energy -= C_DDOT(ndf_+3,B_p_RB[rb],1,C_p_BR[br],1);
  }}

  free_block(B_p_RB);

  double **D_p_BB = block_matrix(noccB_*noccB_,ndf_+3);

  for(int b=0; b<noccB_; b++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][0]),nmoB_,
      &(C_p_BR[b*nvirA_][0]),ndf_+3,0.0,&(D_p_BB[b*noccB_][0]),ndf_+3);
  }

  energy += C_DDOT((long int) noccB_*noccB_*(ndf_+3),&(D_p_BB[0][0]),1,
    &(B_p_BB[0][0]),1);

  free_block(B_p_BB);
  free_block(C_p_BR);
  free_block(D_p_BB);

  double **sAA = block_matrix(aoccA_,aoccA_);

  C_DGEMM('N','T',aoccA_,aoccA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[foccA_][0]),nmoB_,0.0,&(sAA[0][0]),aoccA_);

  double **sRR = block_matrix(nvirA_,nvirA_);

  C_DGEMM('N','T',nvirA_,nvirA_,noccB_,1.0,&(sAB_[noccA_][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(sRR[0][0]),nvirA_);

  double *X = init_array(ndf_+3);

  C_DGEMV('t',aoccA_*aoccA_,ndf_+3,1.0,T_p_AA[0],ndf_+3,&(sAA[0][0]),1,
    0.0,X,1);

  energy -= 2.0*C_DDOT(ndf_+3,X,1,diagBB_,1);

  C_DGEMV('t',nvirA_*nvirA_,ndf_+3,1.0,T_p_RR[0],ndf_+3,&(sRR[0][0]),1,
    0.0,X,1);

  energy -= 2.0*C_DDOT(ndf_+3,X,1,diagBB_,1);

  free(X);
  free_block(sAA);
  free_block(sRR);
  free_block(T_p_AA);
  free_block(T_p_RR);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_6       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

double SAPT2::exch102_k11u_6()
{
  double energy=0.0;

  double *T_BSBS = init_array((long int) aoccB_*nvirB_*aoccB_*nvirB_);
  
  double *tBSBS = init_array((long int) aoccB_*nvirB_*aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tBSBS Amplitudes",
    (char *) tBSBS,sizeof(double)*aoccB_*nvirB_*aoccB_*nvirB_);

  C_DGEMM('N','T',aoccB_*nvirB_,aoccB_*nvirB_,aoccB_*nvirB_,3.0,tBSBS,
    aoccB_*nvirB_,tBSBS,aoccB_*nvirB_,0.0,T_BSBS,aoccB_*nvirB_);

  antisym(tBSBS,aoccB_,nvirB_);
  OVOpVp_to_OVpOpV(tBSBS,aoccB_,nvirB_);

  C_DGEMM('N','T',aoccB_*nvirB_,aoccB_*nvirB_,aoccB_*nvirB_,1.0,tBSBS,
    aoccB_*nvirB_,tBSBS,aoccB_*nvirB_,1.0,T_BSBS,aoccB_*nvirB_);

  ijkl_to_ikjl(T_BSBS,aoccB_,nvirB_,aoccB_,nvirB_);

  free(tBSBS);

  double **B_p_SS = get_SS_ints(1);
  double **T_p_BB = block_matrix(aoccB_*aoccB_,ndf_+3);

  C_DGEMM('N','N',aoccB_*aoccB_,ndf_+3,nvirB_*nvirB_,1.0,T_BSBS,nvirB_*nvirB_,
    &(B_p_SS[0][0]),ndf_+3,0.0,&(T_p_BB[0][0]),ndf_+3);

  free_block(B_p_SS);

  double **B_p_BB = get_BB_ints(1,foccB_,foccB_);
  double **T_p_SS = block_matrix(nvirB_*nvirB_,ndf_+3);

  C_DGEMM('T','N',nvirB_*nvirB_,ndf_+3,aoccB_*aoccB_,1.0,T_BSBS,nvirB_*nvirB_,
    &(B_p_BB[0][0]),ndf_+3,0.0,&(T_p_SS[0][0]),ndf_+3);

  free(T_BSBS);
  free_block(B_p_BB);

  double **B_p_AB = get_AB_ints(1,0,foccB_);
  double **C_p_AB = block_matrix(noccA_*aoccB_,ndf_+3);

  C_DGEMM('N','N',noccA_,aoccB_*(ndf_+3),aoccB_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(T_p_BB[0][0]),aoccB_*(ndf_+3),0.0,&(C_p_AB[0][0]),aoccB_*(ndf_+3));

  energy -= C_DDOT((long int) noccA_*aoccB_*(ndf_+3),&(C_p_AB[0][0]),1,
    &(B_p_AB[0][0]),1);

  free_block(B_p_AB);

  double **B_p_AA = get_AA_ints(1);
  double **C_p_AA = block_matrix(noccA_*noccA_,ndf_+3);

  for(int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',noccA_,ndf_+3,aoccB_,1.0,&(sAB_[0][foccB_]),nmoB_,
      &(C_p_AB[a*aoccB_][0]),ndf_+3,0.0,&(C_p_AA[a*noccA_][0]),ndf_+3);
  }

  energy += C_DDOT((long int) noccA_*noccA_*(ndf_+3),&(C_p_AA[0][0]),1,
    &(B_p_AA[0][0]),1);

  free_block(C_p_AA);
  free_block(C_p_AB);

  double **B_p_AS = get_AS_ints(1);
  double **C_p_AS = block_matrix(noccA_*nvirB_,ndf_+3);

  C_DGEMM('N','N',noccA_,nvirB_*(ndf_+3),nvirB_,1.0,&(sAB_[0][noccB_]),
    nmoB_,&(T_p_SS[0][0]),nvirB_*(ndf_+3),0.0,&(C_p_AS[0][0]),nvirB_*(ndf_+3));

  energy -= C_DDOT((long int) noccA_*nvirB_*(ndf_+3),B_p_AS[0],1,C_p_AS[0],1);

  free_block(B_p_AS);

  double **D_p_AA = block_matrix(noccA_*noccA_,ndf_+3);

  for(int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',noccA_,ndf_+3,nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
      &(C_p_AS[a*nvirB_][0]),ndf_+3,0.0,&(D_p_AA[a*noccA_][0]),ndf_+3);
  }

  energy += C_DDOT((long int) noccA_*noccA_*(ndf_+3),&(D_p_AA[0][0]),1,
    &(B_p_AA[0][0]),1);

  free_block(B_p_AA);
  free_block(C_p_AS);
  free_block(D_p_AA);

  double **sBB = block_matrix(aoccB_,aoccB_);

  C_DGEMM('T','N',aoccB_,aoccB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][foccB_]),nmoB_,0.0,&(sBB[0][0]),aoccB_);

  double **sSS = block_matrix(nvirB_,nvirB_);

  C_DGEMM('T','N',nvirB_,nvirB_,noccA_,1.0,&(sAB_[0][noccB_]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,&(sSS[0][0]),nvirB_);

  double *X = init_array(ndf_+3);

  C_DGEMV('t',aoccB_*aoccB_,ndf_+3,1.0,T_p_BB[0],ndf_+3,&(sBB[0][0]),1,
    0.0,X,1);

  energy -= 2.0*C_DDOT(ndf_+3,X,1,diagAA_,1);

  C_DGEMV('t',nvirB_*nvirB_,ndf_+3,1.0,T_p_SS[0],(ndf_+3),&(sSS[0][0]),1,
    0.0,X,1);

  energy -= 2.0*C_DDOT(ndf_+3,X,1,diagAA_,1);

  free(X);
  free_block(sBB);
  free_block(sSS);
  free_block(T_p_BB);
  free_block(T_p_SS);

  if (debug_) {
    fprintf(outfile,"    Exch12_k11u_6       = %18.12lf H\n",-energy);
    fflush(outfile);
  }

  return(-energy);
}

}}
