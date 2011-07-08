#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::exch12()
{
  double e_exch111 = exch111();

  if (debug_) {
    fprintf(outfile,"    Exch111             = %18.12lf H\n",e_exch111);
  }
/*
  double e_exch120_k2u = exch110(PSIF_SAPT_AMPS,"Theta 2 AR Intermediates");

  if (debug_) {
    fprintf(outfile,"    Exch120 K2u         = %18.12lf H\n",e_exch120_k2u);
  }

  double e_exch102_k2u = exch101(PSIF_SAPT_AMPS,"Theta 2 BS Intermediates");

  if (debug_) {
    fprintf(outfile,"    Exch102 K2u         = %18.12lf H\n",e_exch102_k2u);
  }
*/
  double e_exch120_k2f = exch120_k2f();

  if (debug_) {
    fprintf(outfile,"    Exch120 K2f         = %18.12lf H\n",e_exch120_k2f);
  }

  double e_exch102_k2f = exch102_k2f();

  if (debug_) {
    fprintf(outfile,"    Exch102 K2f         = %18.12lf H\n",e_exch102_k2f);
  }

  double e_exch120_k11u = exch120_k11u_1();

  if (debug_) {
    fprintf(outfile,"    Exch120 K11u        = %18.12lf H\n",e_exch120_k11u);
  }

  double e_exch102_k11u = exch102_k11u_1();

  if (debug_) {
    fprintf(outfile,"    Exch102 K11u        = %18.12lf H\n",e_exch102_k11u);
  }

  e_exch12_ = e_exch111 + e_exch120_k2f + e_exch102_k2f + /*e_exch120_k2u +
    e_exch102_k2u +*/ e_exch120_k11u + e_exch102_k11u; 

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

  C_DGEMM('T','N',noccB_,nvirA_*(ndf_+3),noccA_,1.0,&(sAB_[0][0]),nmo_,
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

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),nvirA_,1.0,&(sAB_[noccA_][0]),nmo_,
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
    nmo_,0.0,zRB[0],noccB_);

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

  C_DGEMM('T','N',noccB_,noccB_,nvirA_,1.0,sAB_[noccA_],nmo_,zRB[0],
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

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,sAB_[0],nmo_,sAB_[0],nmo_,0.0,
    yBB[0],noccB_);

  energy -= 4.0*C_DDOT(noccB_*noccB_,xBB[0],1,yBB[0],1);

  free(X);
  free_block(xBB);
  free_block(yBB);
  free_block(B_p_BB);

  if (debug_) {
    fprintf(outfile,"\n    Exch12_k11u_1       = %18.12lf H\n",-energy);
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

  C_DGEMM('N','N',noccA_,nvirB_*(ndf_+3),noccB_,1.0,&(sAB_[0][0]),nmo_,
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
    C_DGEMM('N','N',noccA_,(ndf_+3),nvirB_,1.0,&(sAB_[0][noccB_]),nmo_,
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

  C_DGEMM('N','N',noccA_,nvirB_,nvirB_,1.0,&(sAB_[0][noccB_]),nmo_,pSS[0],
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
    nmo_,0.0,zAA[0],noccA_);

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

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,sAB_[0],nmo_,sAB_[0],nmo_,0.0,
    yAA[0],noccA_);

  energy -= 4.0*C_DDOT(noccA_*noccA_,xAA[0],1,yAA[0],1);

  free(X);
  free_block(xAA);
  free_block(yAA);
  free_block(B_p_AA);

  if (debug_) {
    fprintf(outfile,"\n    Exch12_k11u_1       = %18.12lf H\n",-energy);
  }

  return(-energy);
}

}}
