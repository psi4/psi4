#include "sapt2p3.h"

using namespace boost;

namespace psi { namespace sapt {

void SAPT2p3::exch_ind30()
{
  double **tAR = block_matrix(noccA_,nvirA_);
  double **vAR = block_matrix(noccA_,nvirA_);

  psio_->read_entry(PSIF_SAPT_AMPS,"Ind30 uAR Amplitudes", (char *) tAR[0],
    sizeof(double)*noccA_*nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"AR Exch-Ind Integrals", (char *) vAR[0],
    sizeof(double)*noccA_*nvirA_);

  double ex_1 = -2.0*C_DDOT(noccA_*nvirA_,tAR[0],1,vAR[0],1);

  free_block(tAR);
  free_block(vAR);

  double **tBS = block_matrix(noccB_,nvirB_);
  double **vBS = block_matrix(noccB_,nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"Ind30 uBS Amplitudes", (char *) tBS[0],
    sizeof(double)*noccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"BS Exch-Ind Integrals", (char *) vBS[0],
    sizeof(double)*noccB_*nvirB_);

  double ex_2 = -2.0*C_DDOT(noccB_*nvirB_,tBS[0],1,vBS[0],1);

  free_block(tBS);
  free_block(vBS);

  double **sAR = block_matrix(noccA_,nvirA_);
  double **sBS = block_matrix(noccB_,nvirB_);

  for (int a=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++) {
      sAR[a][r] = wBAR_[a][r] / (evalsA_[a] - evalsA_[r+noccA_]);
  }}

  for (int b=0; b<noccB_; b++) {
    for (int s=0; s<nvirB_; s++) {
      sBS[b][s] = wABS_[b][s] / (evalsB_[b] - evalsB_[s+noccB_]);
  }}

  double ex_3 = exch_ind30_1(sAR,sBS);
  double ex_4 = exch_ind30_2(sAR);
  double ex_5 = exch_ind30_3(sBS);

  free_block(sAR);
  free_block(sBS);

  e_exch_ind30_ = ex_1 + ex_2 + ex_3 + ex_4 + ex_5;

  if (debug_) {
    fprintf(outfile,"\n    Exch-Ind_1          = %18.12lf H\n",ex_1);
    fprintf(outfile,"    Exch-Ind_2          = %18.12lf H\n",ex_2);
    fprintf(outfile,"    Exch-Ind_3          = %18.12lf H\n",ex_3);
    fprintf(outfile,"    Exch-Ind_4          = %18.12lf H\n",ex_4);
    fprintf(outfile,"    Exch-Ind_5          = %18.12lf H\n",ex_5);
  }
  if (print_) {
    fprintf(outfile,"    Exch-Ind30          = %18.12lf H\n",e_exch_ind30_);
    fflush(outfile);
  }
}

double SAPT2p3::exch_ind30_1(double **sAR, double **sBS)
{
  double energy = 0.0;

  double **vARBS = block_matrix(noccA_*nvirA_,noccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",(char *) vARBS[0],
    sizeof(double)*noccA_*nvirA_*noccB_*nvirB_);

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      energy -= 2.0*sAR[a][r]*C_DDOT(noccB_*nvirB_,vARBS[ar],1,sBS[0],1);
  }}

  free_block(vARBS);

  double **xAR = block_matrix(noccA_,nvirA_);
  double **xBS = block_matrix(noccB_,nvirB_);

  C_DGEMM('N','T',noccA_,nvirA_,noccB_,1.0,sAB_[0],nmoB_,sAB_[noccA_],nmoB_,
    0.0,xAR[0],nvirA_);

  C_DGEMM('T','N',noccB_,nvirB_,noccA_,1.0,sAB_[0],nmoB_,&(sAB_[0][noccB_]),
    nmoB_,0.0,xBS[0],nvirB_);

  double **yAR = block_matrix(noccA_,nvirA_);
  double **yBS = block_matrix(noccB_,nvirB_);

  double **B_p_AR = get_AR_ints(1);
  double **B_p_BS = get_BS_ints(1);

  C_DGEMV('n',noccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,diagBB_,1,
    0.0,&(yAR[0][0]),1);

  C_DGEMV('n',noccB_*nvirB_,ndf_+3,1.0,&(B_p_BS[0][0]),ndf_+3,diagAA_,1,
    0.0,&(yBS[0][0]),1);

  energy += 8.0*C_DDOT(noccA_*nvirA_,xAR[0],1,sAR[0],1)*
    C_DDOT(noccB_*nvirB_,yBS[0],1,sBS[0],1);
  energy += 8.0*C_DDOT(noccA_*nvirA_,yAR[0],1,sAR[0],1)*
    C_DDOT(noccB_*nvirB_,xBS[0],1,sBS[0],1);

  free_block(B_p_AR);
  free_block(B_p_BS);
  free_block(xAR);
  free_block(xBS);
  free_block(yAR);
  free_block(yBS);

  return(energy);
}

double SAPT2p3::exch_ind30_2(double **sAR)
{
  double energy = 0.0;

  double **ssRB = block_matrix(nvirA_,noccB_);

  C_DGEMM('T','N',nvirA_,noccB_,noccA_,1.0,sAR[0],nvirA_,sAB_[0],nmoB_,
    0.0,ssRB[0],noccB_);

  double **A_p_AR = get_AR_ints(1);
  double **B_p_BB = get_BB_ints(1);
  double **B_p_RB = get_RB_ints(1);

  double **C_p_AB = block_matrix(noccA_*noccB_,ndf_+3);
  double **D_p_AB = block_matrix(noccA_*noccB_,ndf_+3);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,ssRB[0],noccB_,A_p_AR[a*nvirA_],
      ndf_+3,0.0,C_p_AB[a*noccB_],ndf_+3);
  }

  C_DGEMM('N','N',noccA_,noccB_*(ndf_+3),nvirA_,1.0,sAR[0],nvirA_,B_p_RB[0],
    noccB_*(ndf_+3),0.0,D_p_AB[0],noccB_*(ndf_+3));

  energy += 2.0*C_DDOT(noccA_*noccB_*(ndf_+3),C_p_AB[0],1,D_p_AB[0],1);

  free_block(C_p_AB);
  free_block(D_p_AB);

  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);

  C_DGEMV('t',noccA_*nvirA_,ndf_+3,1.0,A_p_AR[0],ndf_+3,sAR[0],1,0.0,X,1);
  C_DGEMV('t',nvirA_*noccB_,ndf_+3,1.0,B_p_RB[0],ndf_+3,ssRB[0],1,0.0,Y,1);

  energy -= 4.0*C_DDOT(ndf_+3,X,1,Y,1);

  double **xAB = block_matrix(noccA_,noccB_);
  double **xAR = block_matrix(noccA_,nvirA_);
  double **yAR = block_matrix(noccA_,nvirA_);

  C_DGEMM('N','N',noccA_,noccB_,nvirA_,1.0,sAR[0],nvirA_,sAB_[noccA_],nmoB_,
    0.0,xAB[0],noccB_);

  C_DGEMM('N','T',noccA_,nvirA_,noccB_,1.0,xAB[0],noccB_,ssRB[0],noccB_,
    0.0,xAR[0],nvirA_);

  C_DGEMV('n',noccA_*nvirA_,ndf_+3,1.0,A_p_AR[0],ndf_+3,diagBB_,1,
    0.0,yAR[0],1); 

  energy += 4.0*C_DDOT(noccA_*nvirA_,xAR[0],1,yAR[0],1);

  free_block(xAR);
  free_block(yAR);

  double **E_p_AB = block_matrix(noccA_*noccB_,ndf_+3);
  double **C_p_BB = block_matrix(noccB_*noccB_,ndf_+3);
 
  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,ssRB[0],noccB_,A_p_AR[a*nvirA_],
      ndf_+3,0.0,E_p_AB[a*noccB_],ndf_+3);
  }

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),noccA_,1.0,xAB[0],noccB_,E_p_AB[0],
    noccB_*(ndf_+3),0.0,C_p_BB[0],noccB_*(ndf_+3));

  energy -= 2.0*C_DDOT(noccB_*noccB_*(ndf_+3),C_p_BB[0],1,B_p_BB[0],1);
 
  free_block(xAB);
  free_block(E_p_AB);
  free_block(C_p_BB);

  double **xBB = block_matrix(noccB_,noccB_);

  C_DGEMM('T','N',noccB_,noccB_,nvirA_,1.0,ssRB[0],noccB_,sAB_[noccA_],nmoB_,
    0.0,xBB[0],noccB_);

  C_DGEMV('t',noccB_*noccB_,ndf_+3,1.0,B_p_BB[0],ndf_+3,xBB[0],1,0.0,Y,1);

  energy += 4.0*C_DDOT(ndf_+3,X,1,Y,1);

  free_block(xBB);
  free_block(ssRB);
  free(X);
  free(Y);

  free_block(A_p_AR);
  free_block(B_p_RB);
  free_block(B_p_BB);

  return(energy);
}

double SAPT2p3::exch_ind30_3(double **sBS)
{
  double energy = 0.0;

  double **ssAS = block_matrix(noccA_,nvirB_);

  C_DGEMM('N','N',noccA_,nvirB_,noccB_,1.0,sAB_[0],nmoB_,sBS[0],nvirB_,
    0.0,ssAS[0],nvirB_);

  double **A_p_AA = get_AA_ints(1);
  double **A_p_AS = get_AS_ints(1);
  double **B_p_BS = get_BS_ints(1);

  double **C_p_AB = block_matrix(noccA_*noccB_,ndf_+3);
  double **D_p_AB = block_matrix(noccA_*noccB_,ndf_+3);

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('N','N',noccA_,ndf_+3,nvirB_,1.0,ssAS[0],nvirB_,B_p_BS[b*nvirB_],
      ndf_+3,0.0,C_p_AB[b],noccB_*(ndf_+3));
  }

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',noccB_,ndf_+3,nvirB_,1.0,sBS[0],nvirB_,A_p_AS[a*nvirB_],
      ndf_+3,0.0,D_p_AB[a*noccB_],ndf_+3);
  }

  energy += 2.0*C_DDOT(noccA_*noccB_*(ndf_+3),C_p_AB[0],1,D_p_AB[0],1);

  free_block(C_p_AB);
  free_block(D_p_AB);

  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);

  C_DGEMV('t',noccB_*nvirB_,ndf_+3,1.0,B_p_BS[0],ndf_+3,sBS[0],1,0.0,X,1);
  C_DGEMV('t',noccA_*nvirB_,ndf_+3,1.0,A_p_AS[0],ndf_+3,ssAS[0],1,0.0,Y,1);

  energy -= 4.0*C_DDOT(ndf_+3,X,1,Y,1);

  double **xAB = block_matrix(noccA_,noccB_);
  double **xBS = block_matrix(noccB_,nvirB_);
  double **yBS = block_matrix(noccB_,nvirB_);

  C_DGEMM('N','T',noccA_,noccB_,nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,sBS[0],
    nvirB_,0.0,xAB[0],noccB_);

  C_DGEMM('T','N',noccB_,nvirB_,noccA_,1.0,xAB[0],noccB_,ssAS[0],nvirB_,
    0.0,xBS[0],nvirB_);

  C_DGEMV('n',noccB_*nvirB_,ndf_+3,1.0,B_p_BS[0],ndf_+3,diagAA_,1,
    0.0,yBS[0],1); 

  energy += 4.0*C_DDOT(noccB_*nvirB_,xBS[0],1,yBS[0],1);

  free_block(xBS);
  free_block(yBS);

  double **E_p_AB = block_matrix(noccA_*noccB_,ndf_+3);
  double **C_p_AA = block_matrix(noccA_*noccA_,ndf_+3);
 
  for (int b=0; b<noccB_; b++) {
    C_DGEMM('N','N',noccA_,ndf_+3,nvirB_,1.0,ssAS[0],nvirB_,B_p_BS[b*nvirB_],
      ndf_+3,0.0,E_p_AB[b*noccA_],ndf_+3);
  }

  C_DGEMM('N','N',noccA_,noccA_*(ndf_+3),noccB_,1.0,xAB[0],noccB_,E_p_AB[0],
    noccA_*(ndf_+3),0.0,C_p_AA[0],noccA_*(ndf_+3));

  energy -= 2.0*C_DDOT(noccA_*noccA_*(ndf_+3),C_p_AA[0],1,A_p_AA[0],1);
 
  free_block(xAB);
  free_block(E_p_AB);
  free_block(C_p_AA);

  double **xAA = block_matrix(noccA_,noccA_);

  C_DGEMM('N','T',noccA_,noccA_,nvirB_,1.0,ssAS[0],nvirB_,&(sAB_[0][noccB_]),
    nmoB_,0.0,xAA[0],noccA_);

  C_DGEMV('t',noccA_*noccA_,ndf_+3,1.0,A_p_AA[0],ndf_+3,xAA[0],1,0.0,Y,1);

  energy += 4.0*C_DDOT(ndf_+3,X,1,Y,1);

  free_block(xAA);
  free(X);
  free(Y);

  free_block(ssAS);
  free_block(A_p_AS);
  free_block(A_p_AA);
  free_block(B_p_BS);

  return(energy);
}

}}

