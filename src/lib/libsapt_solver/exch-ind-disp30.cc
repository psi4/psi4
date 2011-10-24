#include "sapt2p3.h"

using namespace boost;

namespace psi { namespace sapt {

void SAPT2p3::exch_ind_disp30()
{
  double **tAR = block_matrix(aoccA_,nvirA_);
  double **vAR = block_matrix(noccA_,nvirA_);

  psio_->read_entry(PSIF_SAPT_AMPS,"IndDisp30 uAR Amplitudes", (char *) tAR[0],
    sizeof(double)*aoccA_*nvirA_);
  psio_->read_entry(PSIF_SAPT_AMPS,"AR Exch-Ind Integrals", (char *) vAR[0],
    sizeof(double)*noccA_*nvirA_);

  double ex_1 = -2.0*C_DDOT(aoccA_*nvirA_,tAR[0],1,vAR[foccA_],1);

  free_block(tAR);
  free_block(vAR);

  double **tBS = block_matrix(aoccB_,nvirB_);
  double **vBS = block_matrix(noccB_,nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"IndDisp30 uBS Amplitudes", (char *) tBS[0],
    sizeof(double)*aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"BS Exch-Ind Integrals", (char *) vBS[0],
    sizeof(double)*noccB_*nvirB_);

  double ex_2 = -2.0*C_DDOT(aoccB_*nvirB_,tBS[0],1,vBS[foccB_],1);

  free_block(tBS);
  free_block(vBS);

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"IndDisp30 uARBS Amplitudes",(char *)
    tARBS[0],sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  double **vARBS = block_matrix(noccA_*nvirA_,noccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",(char *) vARBS[0],
    sizeof(double)*noccA_*nvirA_*noccB_*nvirB_);

  double ex_3 = 0.0;

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      int aarr = (a+foccA_)*nvirA_+r;
      ex_3 -= 2.0*C_DDOT(aoccB_*nvirB_,&(vARBS[aarr][foccB_*nvirB_]),1,
        tARBS[ar],1);
  }}

  free_block(tARBS);
  free_block(vARBS);

  double **sAR = block_matrix(aoccA_,nvirA_);
  double **sBS = block_matrix(aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++) {
      sAR[a][r] = wBAR_[a+foccA_][r] / (evalsA_[a+foccA_] - evalsA_[r+noccA_]);
  }}

  for (int b=0; b<aoccB_; b++) {
    for (int s=0; s<nvirB_; s++) {
      sBS[b][s] = wABS_[b+foccB_][s] / (evalsB_[b+foccB_] - evalsB_[s+noccB_]);
  }}

  double ex_4 = exch_ind_disp30_21(sAR);
  double ex_5 = exch_ind_disp30_12(sBS);

  free_block(sAR);
  free_block(sBS);

  e_exch_ind_disp30_ = ex_1 + ex_2 + ex_3 + ex_4 + ex_5;

  if (debug_) {
    fprintf(outfile,"\n    Exch-Ind-Disp_1     = %18.12lf H\n",ex_1);
    fprintf(outfile,"    Exch-Ind-Disp_2     = %18.12lf H\n",ex_2);
    fprintf(outfile,"    Exch-Ind-Disp_3     = %18.12lf H\n",ex_3);
    fprintf(outfile,"    Exch-Ind-Disp_4     = %18.12lf H\n",ex_4);
    fprintf(outfile,"    Exch-Ind-Disp_5     = %18.12lf H\n",ex_5);
  }
  if (print_) {
    fprintf(outfile,"    Exch-Ind-Disp30     = %18.12lf H\n",
      e_exch_ind_disp30_);
    fflush(outfile);
  }
}

double SAPT2p3::exch_ind_disp30_21(double **sAR)
{
  double energy = 0.0;

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *)
    tARBS[0],sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  double **tAS_RB = block_matrix(nvirA_,aoccB_);
  double **tRB_AS = block_matrix(aoccA_,nvirB_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0,bs=0; b<aoccB_; b++) {
        for (int s=0; s<nvirB_; s++,bs++) {
          tAS_RB[r][b] += tARBS[ar][bs]*sAB_[a+foccA_][s+noccB_];
          tRB_AS[a][s] += tARBS[ar][bs]*sAB_[r+noccA_][b+foccB_];
  }}}}

  double **B_p_AR = get_AR_ints(1,foccA_);
  double **B_p_RB = get_RB_ints(1,foccB_);

  double **xRS = block_matrix(nvirA_,nvirB_);
  double **C_p_AS = block_matrix(aoccA_*nvirB_,ndf_+3);

  C_DGEMM('T','N',nvirA_,nvirB_,aoccA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[foccA_][noccB_]),nmoB_,0.0,&(xRS[0][0]),nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',nvirB_,(ndf_+3),nvirA_,1.0,
      &(xRS[0][0]),nvirB_,&(B_p_AR[a*nvirA_][0]),
      (ndf_+3),0.0,&(C_p_AS[a*nvirB_][0]),(ndf_+3));
  }

  double **xRBS = block_matrix(nvirA_*aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','T',nvirA_*aoccB_,nvirB_,ndf_+3,1.0,&(B_p_RB[0][0]),ndf_+3,
      &(C_p_AS[a*nvirB_][0]),ndf_+3,0.0,&(xRBS[0][0]),nvirB_);
    energy += C_DDOT(nvirA_*aoccB_*nvirB_,tARBS[a*nvirA_],1,xRBS[0],1);
  }

  free_block(xRBS);
  free_block(C_p_AS);

  double **B_p_AB = block_matrix(aoccA_*aoccB_,ndf_+3);
  double **C_p_AB = block_matrix(aoccA_*aoccB_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',aoccB_,ndf_+3,nvirA_,1.0,&(tAS_RB[0][0]),aoccB_,
      &(B_p_AR[a*nvirA_][0]),ndf_+3,0.0,&(B_p_AB[a*aoccB_][0]),ndf_+3);
  }

  C_DGEMM('N','N',aoccA_,aoccB_*(ndf_+3),nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(B_p_RB[0][0]),aoccB_*(ndf_+3),0.0,&(C_p_AB[0][0]),aoccB_*(ndf_+3));

  energy += C_DDOT(aoccA_*aoccB_*(ndf_+3),B_p_AB[0],1,C_p_AB[0],1);

  free_block(B_p_AB);
  free_block(C_p_AB);

  double *X = init_array((ndf_+3));
  double *Y = init_array((ndf_+3));

  C_DGEMV('t',aoccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,sAR[0],1,
    0.0,X,1);

  C_DGEMV('t',nvirA_*aoccB_,ndf_+3,1.0,&(B_p_RB[0][0]),ndf_+3,tAS_RB[0],1,
    0.0,Y,1);

  energy -= 2.0*C_DDOT(ndf_+3,X,1,Y,1);

  free(X);
  free(Y);

  double **C_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('T','N',nvirB_,ndf_+3,nvirA_,1.0,&(xRS[0][0]),nvirB_,
      &(B_p_RB[b][0]),aoccB_*(ndf_+3),0.0,&(C_p_BS[b*nvirB_][0]),ndf_+3);
  }

  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  energy -= 2.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),C_p_BS[0],1,T_p_BS[0],1);

  free_block(xRS);
  free_block(B_p_RB);
  free_block(C_p_BS);
  free_block(T_p_BS);

  double **B_p_BS = get_BS_ints(1);

  double **xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[noccA_][foccB_]),nmoB_,0.0,&(xAB[0][0]),aoccB_);

  C_p_AS = block_matrix(aoccA_*nvirB_,ndf_+3);
  double **C_p_RB = block_matrix(aoccB_*nvirA_,ndf_+3);

  for (int r=0; r<nvirA_; r++) {
    C_DGEMM('T','N',aoccB_,ndf_+3,aoccA_,1.0,&(xAB[0][0]),aoccB_,
      &(B_p_AR[r][0]),nvirA_*(ndf_+3),0.0,&(C_p_RB[r*aoccB_][0]),ndf_+3);
  }

  C_DGEMM('N','N',aoccA_,nvirB_*(ndf_+3),noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(B_p_BS[0][0]),nvirB_*(ndf_+3),0.0,&(C_p_AS[0][0]),nvirB_*(ndf_+3));

  xRBS = block_matrix(nvirA_*aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','T',nvirA_*aoccB_,nvirB_,ndf_+3,1.0,&(C_p_RB[0][0]),ndf_+3,
      &(C_p_AS[a*nvirB_][0]),ndf_+3,0.0,&(xRBS[0][0]),nvirB_);
    energy -= C_DDOT(nvirA_*aoccB_*nvirB_,tARBS[a*nvirA_],1,xRBS[0],1);
  }

  free_block(xAB);
  free_block(xRBS);
  free_block(C_p_AS);
  free_block(C_p_RB);

  double **xRB = block_matrix(nvirA_,noccB_);

  C_DGEMM('T','N',nvirA_,noccB_,aoccA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[foccA_][0]),nmoB_,0.0,&(xRB[0][0]),noccB_);

  B_p_AB = block_matrix(aoccA_*noccB_,ndf_+3);
  C_p_AB = block_matrix(aoccA_*noccB_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,&(xRB[0][0]),noccB_,
      &(B_p_AR[a*nvirA_][0]),ndf_+3,0.0,&(B_p_AB[a*noccB_][0]),ndf_+3);
  }

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,nvirB_,1.0,&(tRB_AS[0][0]),nvirB_,
      &(B_p_BS[b*nvirB_][0]),ndf_+3,0.0,&(C_p_AB[b][0]),noccB_*(ndf_+3));
  }

  energy -= C_DDOT(aoccA_*noccB_*(ndf_+3),B_p_AB[0],1,C_p_AB[0],1);

  free_block(xRB);
  free_block(B_p_AB);
  free_block(C_p_AB);

  double **xAR = block_matrix(aoccA_,nvirA_);

  C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(xAR[0][0]),nvirA_);

  double **xAA = block_matrix(aoccA_,aoccA_);

  C_DGEMM('N','T',aoccA_,aoccA_,nvirA_,1.0,&(xAR[0][0]),nvirA_,
    &(sAR[0][0]),nvirA_,0.0,&(xAA[0][0]),aoccA_);

  double **xRR = block_matrix(nvirA_,nvirA_);

  C_DGEMM('T','N',nvirA_,nvirA_,aoccA_,1.0,&(xAR[0][0]),nvirA_,
    &(sAR[0][0]),nvirA_,0.0,&(xRR[0][0]),nvirA_);

  double **C_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);

  C_DGEMM('N','N',aoccA_,nvirA_*(ndf_+3),aoccA_,1.0,&(xAA[0][0]),aoccA_,
    &(B_p_AR[0][0]),nvirA_*(ndf_+3),0.0,&(C_p_AR[0][0]),nvirA_*(ndf_+3));

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',nvirA_,ndf_+3,nvirA_,1.0,&(xRR[0][0]),nvirA_,
      &(B_p_AR[a*nvirA_][0]),ndf_+3,1.0,&(C_p_AR[a*nvirA_][0]),ndf_+3);
  }

  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  energy += 2.0*C_DDOT(aoccA_*nvirA_*(ndf_+3),C_p_AR[0],1,T_p_AR[0],1);

  free_block(xAR);
  free_block(xAA);
  free_block(xRR);
  free_block(C_p_AR);
  free_block(T_p_AR);

  xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[noccA_][foccB_]),nmoB_,0.0,&(xAB[0][0]),aoccB_);

  double **xBB = block_matrix(aoccB_,noccB_);

  C_DGEMM('T','N',aoccB_,noccB_,aoccA_,1.0,
    &(xAB[0][0]),aoccB_,&(sAB_[foccA_][0]),nmoB_,
    0.0,&(xBB[0][0]),noccB_);

  T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  C_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('N','N',aoccB_,nvirB_*(ndf_+3),noccB_,1.0,&(xBB[0][0]),noccB_,
    &(B_p_BS[0][0]),nvirB_*(ndf_+3),0.0,&(C_p_BS[0][0]),nvirB_*(ndf_+3));

  energy += 2.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),C_p_BS[0],1,T_p_BS[0],1);

  free_block(xAB);
  free_block(xBB);
  free_block(C_p_BS);
  free_block(T_p_BS);

  double **xBS = block_matrix(noccB_,nvirB_);

  C_DGEMM('T','N',noccB_,nvirB_,aoccA_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(tRB_AS[0][0]),nvirB_,0.0,&(xBS[0][0]),nvirB_);

  X = init_array((ndf_+3));
  Y = init_array((ndf_+3));

  C_DGEMV('t',aoccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,sAR[0],1,
    0.0,X,1);

  C_DGEMV('t',noccB_*nvirB_,ndf_+3,1.0,&(B_p_BS[0][0]),ndf_+3,xBS[0],1,
    0.0,Y,1);

  energy += 2.0*C_DDOT((ndf_+3),X,1,Y,1);

  free(X);
  free(Y);
  free_block(xBS);
  free_block(B_p_BS);

  double **B_p_BB = get_BB_ints(1,0,foccB_);

  xAB = block_matrix(aoccA_,noccB_);

  C_DGEMM('N','N',aoccA_,noccB_,nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(xAB[0][0]),noccB_);

  B_p_AB = block_matrix(aoccA_*aoccB_,(ndf_+3));

  C_DGEMM('N','N',aoccA_,aoccB_*(ndf_+3),noccB_,1.0,&(xAB[0][0]),noccB_,
    &(B_p_BB[0][0]),aoccB_*(ndf_+3),0.0,&(B_p_AB[0][0]),aoccB_*(ndf_+3));

  xRB = block_matrix(nvirA_,aoccB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','T',nvirA_,aoccB_,ndf_+3,1.0,&(B_p_AR[a*nvirA_][0]),ndf_+3,
      &(B_p_AB[a*aoccB_][0]),(ndf_+3),1.0,&(xRB[0][0]),aoccB_);
  }

  energy -= C_DDOT(nvirA_*aoccB_,xRB[0],1,tAS_RB[0],1);

  free_block(xAB);
  free_block(xRB);
  free_block(B_p_AB);

  xRS = block_matrix(nvirA_,nvirB_);

  C_DGEMM('T','N',nvirA_,nvirB_,aoccA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[foccA_][noccB_]),nmoB_,0.0,&(xRS[0][0]),nvirB_);

  C_p_AS = block_matrix(aoccA_*nvirB_,ndf_+3);
  C_p_RB = block_matrix(aoccB_*nvirA_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',nvirB_,(ndf_+3),nvirA_,1.0,&(xRS[0][0]),nvirB_,
      &(B_p_AR[a*nvirA_][0]),ndf_+3,0.0,&(C_p_AS[a*nvirB_][0]),ndf_+3);
  }

  C_DGEMM('N','N',nvirA_,aoccB_*(ndf_+3),noccB_,1.0,&(sAB_[noccA_][0]),
    nmoB_,&(B_p_BB[0][0]),aoccB_*(ndf_+3),0.0,&(C_p_RB[0][0]),aoccB_*(ndf_+3));

  xRBS = block_matrix(nvirA_*aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','T',nvirA_*aoccB_,nvirB_,ndf_+3,1.0,&(C_p_RB[0][0]),ndf_+3,
      &(C_p_AS[a*nvirB_][0]),ndf_+3,0.0,&(xRBS[0][0]),nvirB_);
    energy -= C_DDOT(nvirA_*aoccB_*nvirB_,tARBS[a*nvirA_],1,xRBS[0],1);
  }

  free_block(xRS);
  free_block(xRBS);
  free_block(B_p_BB);
  free_block(C_p_AS);
  free_block(C_p_RB);

  B_p_BB = get_BB_ints(1,foccB_,0);

  xAB = block_matrix(aoccA_,noccB_);

  C_DGEMM('N','N',aoccA_,noccB_,nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(xAB[0][0]),noccB_);

  xBS = block_matrix(noccB_,nvirB_);

  C_DGEMM('T','N',noccB_,nvirB_,aoccA_,1.0,&(xAB[0][0]),noccB_,
    &(sAB_[foccA_][noccB_]),nmoB_,0.0,&(xBS[0][0]),nvirB_);

  C_p_BS = block_matrix(aoccB_*nvirB_,(ndf_+3));

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('T','N',nvirB_,ndf_+3,noccB_,1.0,&(xBS[0][0]),nvirB_,
      &(B_p_BB[b*noccB_][0]),ndf_+3,0.0,&(C_p_BS[b*nvirB_][0]),ndf_+3);
  }

  T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  energy += 2.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),C_p_BS[0],1,T_p_BS[0],1);

  free_block(xAB);
  free_block(xBS);
  free_block(C_p_BS);
  free_block(T_p_BS);

  X = init_array(ndf_+3);

  C_DGEMV('t',aoccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,sAR[0],1,
    0.0,X,1);

  xBB = block_matrix(aoccB_,noccB_);

  C_DGEMV('n',aoccB_*noccB_,ndf_+3,1.0,&(B_p_BB[0][0]),ndf_+3,X,1,
    0.0,xBB[0],1);

  xRB = block_matrix(nvirA_,aoccB_);

  C_DGEMM('N','T',nvirA_,aoccB_,noccB_,1.0,&(sAB_[noccA_][0]),nmoB_,
    &(xBB[0][0]),noccB_,0.0,&(xRB[0][0]),aoccB_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<aoccB_; b++) {
        energy += 2.0*xRB[r][b]*C_DDOT(nvirB_,&(tARBS[ar][b*nvirB_]),1,
          &(sAB_[a+foccA_][noccB_]),1);
  }}}

  free(X);
  free_block(xBB);
  free_block(xRB);
  free_block(B_p_BB);

  xAR = block_matrix(aoccA_,nvirA_);

  C_DGEMV('n',aoccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,diagBB_,1,
    0.0,xAR[0],1);

  xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[noccA_][foccB_]),nmoB_,0.0,&(xAB[0][0]),aoccB_);

  xRB = block_matrix(nvirA_,aoccB_);

  C_DGEMM('T','N',nvirA_,aoccB_,aoccA_,1.0,&(xAR[0][0]),nvirA_,&(xAB[0][0]),
    aoccB_,0.0,&(xRB[0][0]),aoccB_);

  energy += 2.0*C_DDOT(nvirA_*aoccB_,xRB[0],1,tAS_RB[0],1);

  free_block(xAB);
  free_block(xRB);

  xAA = block_matrix(aoccA_,aoccA_);

  C_DGEMM('N','T',aoccA_,aoccA_,nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(xAR[0][0]),nvirA_,0.0,&(xAA[0][0]),aoccA_);

  double **xAS = block_matrix(aoccA_,nvirB_);

  C_DGEMM('T','N',aoccA_,nvirB_,aoccA_,1.0,&(xAA[0][0]),aoccA_,
    &(sAB_[foccA_][noccB_]),nmoB_,0.0,&(xAS[0][0]),nvirB_);

  energy += 2.0*C_DDOT(aoccA_*nvirB_,xAS[0],1,tRB_AS[0],1);

  free_block(xAR);
  free_block(xAA);
  free_block(xAS);

  xAR = block_matrix(aoccA_,nvirA_);

  C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(xAR[0][0]),nvirA_);

  T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  X = init_array((ndf_+3));
  Y = init_array((ndf_+3));

  C_DGEMV('t',aoccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,sAR[0],1,
    0.0,X,1);

  C_DGEMV('t',aoccA_*nvirA_,ndf_+3,1.0,&(T_p_AR[0][0]),ndf_+3,xAR[0],1,
    0.0,Y,1);

  energy -= 4.0*C_DDOT(ndf_+3,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xAR);
  free_block(T_p_AR);

  xAR = block_matrix(aoccA_,nvirA_);

  C_DGEMV('n',aoccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,diagBB_,1,
    0.0,xAR[0],1);

  xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,&(sAR[0][0]),nvirA_,
    &(sAB_[noccA_][foccB_]),nmoB_,0.0,&(xAB[0][0]),aoccB_);

  xBS = block_matrix(aoccB_,nvirB_);

  C_DGEMM('T','N',aoccB_,nvirB_,aoccA_,1.0,&(xAB[0][0]),aoccB_,
    &(sAB_[foccA_][noccB_]),nmoB_,0.0,&(xBS[0][0]),nvirB_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      energy -= 4.0*xAR[a][r]*C_DDOT(aoccB_*nvirB_,&(tARBS[ar][0]),1,
        &(xBS[0][0]),1);
  }}

  free_block(tRB_AS);
  free_block(tAS_RB);
  free_block(xAR);
  free_block(xAB);
  free_block(xBS);
  free_block(B_p_AR);
  free_block(tARBS);

  return(2.0*energy);
}

double SAPT2p3::exch_ind_disp30_12(double **sBS)
{
  double energy = 0.0;

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *)
    tARBS[0],sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  double **tAS_RB = block_matrix(nvirA_,aoccB_);
  double **tRB_AS = block_matrix(aoccA_,nvirB_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0,bs=0; b<aoccB_; b++) {
        for (int s=0; s<nvirB_; s++,bs++) {
          tAS_RB[r][b] += tARBS[ar][bs]*sAB_[a+foccA_][s+noccB_];
          tRB_AS[a][s] += tARBS[ar][bs]*sAB_[r+noccA_][b+foccB_];
  }}}}

  double **B_p_BS = get_BS_ints(1,foccB_);
  double **B_p_AS = get_AS_ints(1,foccA_);

  double **xRS = block_matrix(nvirA_,nvirB_);
  double **C_p_RB = block_matrix(nvirA_*aoccB_,ndf_+3);

  C_DGEMM('N','N',nvirA_,nvirB_,aoccB_,1.0,&(sAB_[noccA_][foccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xRS[0][0]),nvirB_);

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',nvirA_,(ndf_+3),nvirB_,1.0,
      &(xRS[0][0]),nvirB_,&(B_p_BS[b*nvirB_][0]),
      (ndf_+3),0.0,&(C_p_RB[b][0]),aoccB_*(ndf_+3));
  }

  double **xRBS = block_matrix(nvirA_*aoccB_,
    nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','T',nvirA_*aoccB_,nvirB_,
      (ndf_+3),1.0,&(C_p_RB[0][0]),(ndf_+3),
      &(B_p_AS[a*nvirB_][0]),(ndf_+3),0.0,&(xRBS[0][0]),
      nvirB_);
    energy += C_DDOT(nvirA_*aoccB_*nvirB_,
      tARBS[a*nvirA_],1,xRBS[0],1);
  }

  free_block(xRBS);
  free_block(C_p_RB);

  double **B_p_AB = block_matrix(aoccA_*aoccB_,(ndf_+3));
  double **C_p_AB = block_matrix(aoccA_*aoccB_,(ndf_+3));

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',aoccA_,(ndf_+3),nvirB_,1.0,
      &(tRB_AS[0][0]),nvirB_,&(B_p_BS[b*nvirB_][0]),
      (ndf_+3),0.0,&(B_p_AB[b][0]),aoccB_*(ndf_+3));
  }

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',aoccB_,(ndf_+3),nvirB_,1.0,
      &(sBS[0][0]),nvirB_,&(B_p_AS[a*nvirB_][0]),
      (ndf_+3),0.0,&(C_p_AB[a*aoccB_][0]),(ndf_+3));
  }

  energy += C_DDOT(aoccA_*aoccB_*(ndf_+3),B_p_AB[0],1,C_p_AB[0],1);

  free_block(B_p_AB);
  free_block(C_p_AB);

  double *X = init_array((ndf_+3));
  double *Y = init_array((ndf_+3));

  C_DGEMV('t',aoccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),sBS[0],1,0.0,X,1);

  C_DGEMV('t',aoccA_*nvirB_,(ndf_+3),1.0,
    &(B_p_AS[0][0]),(ndf_+3),tRB_AS[0],1,0.0,Y,1);

  energy -= 2.0*C_DDOT((ndf_+3),X,1,Y,1);

  free(X);
  free(Y);

  double **C_p_AR = block_matrix(aoccA_*nvirA_,(ndf_+3));

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',nvirA_,(ndf_+3),nvirB_,1.0,
      &(xRS[0][0]),nvirB_,&(B_p_AS[a*nvirB_][0]),
      (ndf_+3),0.0,&(C_p_AR[a*nvirA_][0]),(ndf_+3));
  }

  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  energy -= 2.0*C_DDOT(aoccA_*nvirA_*(ndf_+3),C_p_AR[0],1,T_p_AR[0],1);

  free_block(xRS);
  free_block(B_p_AS);
  free_block(C_p_AR);
  free_block(T_p_AR);

  double **B_p_AR = get_AR_ints(1);

  double **xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,&(sAB_[foccA_][noccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xAB[0][0]),aoccB_);

  double **C_p_AS = block_matrix(aoccA_*nvirB_,(ndf_+3));
  C_p_RB = block_matrix(aoccB_*nvirA_,(ndf_+3));

  C_DGEMM('N','N',aoccA_,nvirB_*(ndf_+3),aoccB_,1.0,&(xAB[0][0]),aoccB_,
    &(B_p_BS[0][0]),nvirB_*(ndf_+3),0.0,&(C_p_AS[0][0]),nvirB_*(ndf_+3));

  for (int r=0; r<nvirA_; r++) {
    C_DGEMM('T','N',aoccB_,(ndf_+3),noccA_,1.0,
      &(sAB_[0][foccB_]),nmoB_,&(B_p_AR[r][0]),
      nvirA_*(ndf_+3),0.0,&(C_p_RB[r*aoccB_][0]),(ndf_+3));
  }

  xRBS = block_matrix(nvirA_*aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','T',nvirA_*aoccB_,nvirB_,
      (ndf_+3),1.0,&(C_p_RB[0][0]),(ndf_+3),
      &(C_p_AS[a*nvirB_][0]),(ndf_+3),0.0,&(xRBS[0][0]),nvirB_);
    energy -= C_DDOT(nvirA_*aoccB_*nvirB_,tARBS[a*nvirA_],1,xRBS[0],1);
  }

  free_block(xAB);
  free_block(xRBS);
  free_block(C_p_AS);
  free_block(C_p_RB);

  double **xAS = block_matrix(noccA_,nvirB_);

  C_DGEMM('N','N',noccA_,nvirB_,aoccB_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xAS[0][0]),nvirB_);

  B_p_AB = block_matrix(noccA_*aoccB_,(ndf_+3));
  C_p_AB = block_matrix(noccA_*aoccB_,(ndf_+3));

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',noccA_,(ndf_+3),nvirB_,1.0,
      &(xAS[0][0]),nvirB_,&(B_p_BS[b*nvirB_][0]),
      (ndf_+3),0.0,&(B_p_AB[b][0]),aoccB_*(ndf_+3));
  }

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',aoccB_,(ndf_+3),nvirA_,1.0,
      &(tAS_RB[0][0]),aoccB_,&(B_p_AR[a*nvirA_][0]),
      (ndf_+3),0.0,&(C_p_AB[a*aoccB_][0]),(ndf_+3));
  }

  energy -= C_DDOT(noccA_*aoccB_*(ndf_+3),B_p_AB[0],1,C_p_AB[0],1);

  free_block(xAS);
  free_block(B_p_AB);
  free_block(C_p_AB);

  double **xBS = block_matrix(aoccB_,nvirB_);

  C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,&(xBS[0][0]),nvirB_);

  double **xBB = block_matrix(aoccB_,aoccB_);

  C_DGEMM('N','T',aoccB_,aoccB_,nvirB_,1.0,
    &(xBS[0][0]),nvirB_,&(sBS[0][0]),nvirB_,
    0.0,&(xBB[0][0]),aoccB_);

  double **xSS = block_matrix(nvirB_,nvirB_);

  C_DGEMM('T','N',nvirB_,nvirB_,aoccB_,1.0,
    &(xBS[0][0]),nvirB_,&(sBS[0][0]),nvirB_,
    0.0,&(xSS[0][0]),nvirB_);

  double **C_p_BS = block_matrix(aoccB_*nvirB_,
    (ndf_+3));

  C_DGEMM('N','N',aoccB_,nvirB_*(ndf_+3),
    aoccB_,1.0,&(xBB[0][0]),aoccB_,&(B_p_BS[0][0]),
    nvirB_*(ndf_+3),0.0,&(C_p_BS[0][0]),nvirB_*(ndf_+3));

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',nvirB_,(ndf_+3),nvirB_,1.0,
      &(xSS[0][0]),nvirB_,&(B_p_BS[b*nvirB_][0]),
      (ndf_+3),1.0,&(C_p_BS[b*nvirB_][0]),(ndf_+3));
  }

  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  energy += 2.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),
    C_p_BS[0],1,T_p_BS[0],1);

  free_block(xBS);
  free_block(xBB);
  free_block(xSS);
  free_block(C_p_BS);
  free_block(T_p_BS);

  xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,&(sAB_[foccA_][noccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xAB[0][0]),aoccB_);

  double **xAA = block_matrix(aoccA_,noccA_);

  C_DGEMM('N','T',aoccA_,noccA_,aoccB_,1.0,
    &(xAB[0][0]),aoccB_,&(sAB_[0][foccB_]),nmoB_,
    0.0,&(xAA[0][0]),noccA_);

  T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  C_p_AR = block_matrix(aoccA_*nvirA_,(ndf_+3));

  C_DGEMM('N','N',aoccA_,nvirA_*(ndf_+3),
    noccA_,1.0,&(xAA[0][0]),noccA_,&(B_p_AR[0][0]),
    nvirA_*(ndf_+3),0.0,&(C_p_AR[0][0]),nvirA_*(ndf_+3));

  energy += 2.0*C_DDOT(aoccA_*nvirA_*(ndf_+3),C_p_AR[0],1,T_p_AR[0],1);

  free_block(xAB);
  free_block(xAA);
  free_block(C_p_AR);
  free_block(T_p_AR);

  double **xAR = block_matrix(noccA_,nvirA_);

  C_DGEMM('N','T',noccA_,nvirA_,aoccB_,1.0,
    &(sAB_[0][foccB_]),nmoB_,&(tAS_RB[0][0]),aoccB_,
    0.0,&(xAR[0][0]),nvirA_);

  X = init_array((ndf_+3));
  Y = init_array((ndf_+3));

  C_DGEMV('t',aoccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),sBS[0],1,0.0,X,1);

  C_DGEMV('t',noccA_*nvirA_,(ndf_+3),1.0,
    &(B_p_AR[0][0]),(ndf_+3),xAR[0],1,0.0,Y,1);

  energy += 2.0*C_DDOT((ndf_+3),X,1,Y,1);

  free(X);
  free(Y);
  free_block(xAR);
  free_block(B_p_AR);

  double **B_p_AA = get_AA_ints(1,0,foccA_);

  xAB = block_matrix(noccA_,aoccB_);

  C_DGEMM('N','T',noccA_,aoccB_,nvirB_,1.0,
    &(sAB_[0][noccB_]),nmoB_,&(sBS[0][0]),nvirB_,0.0,&(xAB[0][0]),
    aoccB_);

  double **B_p_BA = block_matrix(aoccA_*aoccB_,
    (ndf_+3));

  C_DGEMM('T','N',aoccB_,aoccA_*(ndf_+3),
    noccA_,1.0,&(xAB[0][0]),aoccB_,&(B_p_AA[0][0]),
    aoccA_*(ndf_+3),0.0,&(B_p_BA[0][0]),aoccA_*(ndf_+3));

  xAS = block_matrix(aoccA_,nvirB_);

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','T',aoccA_,nvirB_,(ndf_+3),
      1.0,&(B_p_BA[b*aoccA_][0]),(ndf_+3),
      &(B_p_BS[b*nvirB_][0]),(ndf_+3),1.0,&(xAS[0][0]),nvirB_);
  }

  energy -= C_DDOT(aoccA_*nvirB_,xAS[0],1,tRB_AS[0],1);

  free_block(xAB);
  free_block(xAS);
  free_block(B_p_BA);
  free_block(B_p_AA);

  B_p_AA = get_AA_ints(1,foccA_,0);

  xRS = block_matrix(nvirA_,nvirB_);

  C_DGEMM('N','N',nvirA_,nvirB_,aoccB_,1.0,&(sAB_[noccA_][foccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xRS[0][0]),nvirB_);

  C_p_AS = block_matrix(aoccA_*nvirB_,(ndf_+3));
  C_p_RB = block_matrix(aoccB_*nvirA_,(ndf_+3));

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',nvirA_,(ndf_+3),nvirB_,1.0,
      &(xRS[0][0]),nvirB_,&(B_p_BS[b*nvirB_][0]),
      (ndf_+3),0.0,&(C_p_RB[b][0]),aoccB_*(ndf_+3));
  }

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',nvirB_,(ndf_+3),noccA_,1.0,&(sAB_[0][noccB_]),nmoB_,
      &(B_p_AA[a*noccA_][0]),(ndf_+3),0.0,&(C_p_AS[a*nvirB_][0]),(ndf_+3));
  }

  xRBS = block_matrix(nvirA_*aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','T',nvirA_*aoccB_,nvirB_,
      (ndf_+3),1.0,&(C_p_RB[0][0]),(ndf_+3),
      &(C_p_AS[a*nvirB_][0]),(ndf_+3),0.0,&(xRBS[0][0]),nvirB_);
    energy -= C_DDOT(nvirA_*aoccB_*nvirB_,tARBS[a*nvirA_],1,xRBS[0],1);
  }

  free_block(xRS);
  free_block(xRBS);
  free_block(C_p_AS);
  free_block(C_p_RB);

  xAB = block_matrix(noccA_,aoccB_);

  C_DGEMM('N','T',noccA_,aoccB_,nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xAB[0][0]),aoccB_);

  xAR = block_matrix(noccA_,nvirA_);

  C_DGEMM('N','T',noccA_,nvirA_,aoccB_,1.0,
    &(xAB[0][0]),aoccB_,&(sAB_[noccA_][foccB_]),
    nmoB_,0.0,&(xAR[0][0]),nvirA_);

  C_p_AR = block_matrix(aoccA_*nvirA_,(ndf_+3));

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',nvirA_,(ndf_+3),noccA_,1.0,
      &(xAR[0][0]),nvirA_,&(B_p_AA[a*noccA_][0]),
      (ndf_+3),0.0,&(C_p_AR[a*nvirA_][0]),(ndf_+3));
  }

  T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  energy += 2.0*C_DDOT(aoccA_*nvirA_*(ndf_+3),C_p_AR[0],1,T_p_AR[0],1);

  free_block(xAB);
  free_block(xAR);
  free_block(C_p_AR);
  free_block(T_p_AR);

  X = init_array((ndf_+3));

  C_DGEMV('t',aoccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),sBS[0],1,0.0,X,1);

  xAA = block_matrix(aoccA_,noccA_);

  C_DGEMV('n',aoccA_*noccA_,(ndf_+3),1.0,
    &(B_p_AA[0][0]),(ndf_+3),X,1,0.0,xAA[0],1);

  xAS = block_matrix(aoccA_,nvirB_);

  C_DGEMM('N','N',aoccA_,nvirB_,noccA_,1.0,
    &(xAA[0][0]),noccA_,&(sAB_[0][noccB_]),
    nmoB_,0.0,&(xAS[0][0]),nvirB_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<aoccB_; b++) {
        energy += 2.0*sAB_[r+noccA_][b+foccB_]*
          C_DDOT(nvirB_,&(tARBS[ar][b*nvirB_]),1,
          &(xAS[a][0]),1);
  }}}

  free(X);
  free_block(xAA);
  free_block(xAS);

  xBS = block_matrix(aoccB_,nvirB_);

  C_DGEMV('n',aoccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),diagAA_,1,0.0,xBS[0],1);

  xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,
    &(sAB_[foccA_][noccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xAB[0][0]),
    aoccB_);

  xAS = block_matrix(aoccA_,nvirB_);

  C_DGEMM('N','N',aoccA_,nvirB_,aoccB_,1.0,
    &(xAB[0][0]),aoccB_,&(xBS[0][0]),nvirB_,0.0,
    &(xAS[0][0]),nvirB_);

  energy + 2.0*C_DDOT(aoccA_*nvirB_,xAS[0],1,tRB_AS[0],1);

  free_block(xAB);
  free_block(xAS);

  xBB = block_matrix(aoccB_,aoccB_);

  C_DGEMM('N','T',aoccB_,aoccB_,nvirB_,1.0,
    &(sBS[0][0]),nvirB_,&(xBS[0][0]),nvirB_,
    0.0,&(xBB[0][0]),aoccB_);

  double **xRB = block_matrix(nvirA_,aoccB_);

  C_DGEMM('N','N',nvirA_,aoccB_,aoccB_,1.0,
    &(sAB_[noccA_][foccB_]),nmoB_,&(xBB[0][0]),
    aoccB_,0.0,&(xRB[0][0]),aoccB_);

  energy += 2.0*C_DDOT(nvirA_*aoccB_,xRB[0],1,tAS_RB[0],1);

  free_block(xBS);
  free_block(xBB);
  free_block(xRB);
  free_block(B_p_AA);

  xBS = block_matrix(aoccB_,nvirB_);

  C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,&(xBS[0][0]),nvirB_);

  T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  X = init_array((ndf_+3));
  Y = init_array((ndf_+3));

  C_DGEMV('t',aoccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),sBS[0],1,0.0,X,1);

  C_DGEMV('t',aoccB_*nvirB_,(ndf_+3),1.0,
    &(T_p_BS[0][0]),(ndf_+3),xBS[0],1,0.0,Y,1);

  energy -= 4.0*C_DDOT((ndf_+3),X,1,Y,1);

  free(X);
  free(Y);
  free_block(xBS);
  free_block(T_p_BS);

  xBS = block_matrix(aoccB_,nvirB_);

  C_DGEMV('n',aoccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),diagAA_,1,0.0,xBS[0],1);

  xAB = block_matrix(aoccA_,aoccB_);

  C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,
    &(sAB_[foccA_][noccB_]),nmoB_,
    &(sBS[0][0]),nvirB_,0.0,&(xAB[0][0]),aoccB_);

  xAR = block_matrix(aoccA_,nvirA_);

  C_DGEMM('N','T',aoccA_,nvirA_,aoccB_,1.0,
    &(xAB[0][0]),aoccB_,&(sAB_[noccA_][foccB_]),
    nmoB_,0.0,&(xAR[0][0]),nvirA_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      energy -= 4.0*xAR[a][r]*C_DDOT(aoccB_*nvirB_,
        &(tARBS[ar][0]),1,&(xBS[0][0]),1);
  }}

  free_block(tRB_AS);
  free_block(tAS_RB);
  free_block(xAR);
  free_block(xAB);
  free_block(xBS);
  free_block(B_p_BS);
  free_block(tARBS);

  return(2.0*energy);
}

}}

