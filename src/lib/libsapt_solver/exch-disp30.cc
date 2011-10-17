#include "sapt2p3.h"

using namespace boost;

namespace psi { namespace sapt {

void SAPT2p3::exch_disp30()
{
  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"Disp30 uARBS Amplitudes",(char *)
    tARBS[0],sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  double **vARBS = block_matrix(noccA_*nvirA_,noccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",(char *) vARBS[0],
    sizeof(double)*noccA_*nvirA_*noccB_*nvirB_);

  double ex_1 = 0.0;

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      int aarr = (a+foccA_)*nvirA_+r;
      ex_1 -= 2.0*C_DDOT(aoccB_*nvirB_,&(vARBS[aarr][foccB_*nvirB_]),1,
        tARBS[ar],1);
  }}

  free_block(tARBS);
  free_block(vARBS);

  double ex_2 = exch_disp30_20();
  double ex_3 = exch_disp30_02();
  double ex_4 = exch_disp30_22();

  e_exch_disp30_ = ex_1 + ex_2 + ex_3 + ex_4;

  if (debug_) {
    fprintf(outfile,"\n    Exch-Disp_1         = %18.12lf H\n",ex_1);
    fprintf(outfile,"    Exch-Disp_2         = %18.12lf H\n",ex_2);
    fprintf(outfile,"    Exch-Disp_3         = %18.12lf H\n",ex_3);
    fprintf(outfile,"    Exch-Disp_4         = %18.12lf H\n",ex_4);
  }
  if (print_) {
    fprintf(outfile,"    Exch-Disp30         = %18.12lf H\n",
      e_exch_disp30_);
    fflush(outfile);
  }
}

double SAPT2p3::exch_disp30_20()
{
  double energy = 0.0;

  double **uARAR = block_matrix(aoccA_*nvirA_,aoccA_*nvirA_);
  double **B_p_AR = get_AR_ints(1,foccA_);
  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  C_DGEMM('N','T',aoccA_*nvirA_,aoccA_*nvirA_,ndf_+3,1.0,&(B_p_AR[0][0]),
    ndf_+3,&(T_p_AR[0][0]),ndf_+3,0.0,&(uARAR[0][0]),aoccA_*nvirA_);

  free_block(T_p_AR);

  for(int ar=0; ar<aoccA_*nvirA_; ar++) {
    for(int a1r1=0; a1r1<ar; a1r1++) {
      double tval = uARAR[ar][a1r1] + uARAR[a1r1][ar];
      uARAR[a1r1][ar] = tval;
      uARAR[ar][a1r1] = tval;
  }}

  C_DSCAL(aoccA_*nvirA_,2.0,&(uARAR[0][0]),
    aoccA_*nvirA_+1);

  for (int a=0, ar=0; a < aoccA_; a++) {
  for (int r=0; r < nvirA_; r++, ar++) {
    for (int aa=0, aarr=0; aa < aoccA_; aa++) {
    for (int rr=0; rr < nvirA_; rr++, aarr++) {
      double denom = evalsA_[a+foccA_]+evalsA_[aa+foccA_]-
        evalsA_[r+noccA_]-evalsA_[rr+noccA_];
      uARAR[ar][aarr] /= denom;
    }}
  }}

  double **U_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);

  C_DGEMM('N','N',aoccA_*nvirA_,ndf_+3,aoccA_*nvirA_,1.0,&(uARAR[0][0]),
    aoccA_*nvirA_,&(B_p_AR[0][0]),ndf_+3,0.0,&(U_p_AR[0][0]),ndf_+3);

  double *X = init_array(nvirA_);

  for(int a=0; a<aoccA_; a++) {
  for(int a1=0; a1<=a; a1++) {
    for(int r=0; r<nvirA_; r++) {
      int ar = a*nvirA_+r;
      int a1r = a1*nvirA_+r;
      C_DCOPY(nvirA_,&(uARAR[ar][a1*nvirA_]),1,X,1);
      C_DCOPY(nvirA_,&(uARAR[a1r][a*nvirA_]),1,
        &(uARAR[ar][a1*nvirA_]),1);
      C_DCOPY(nvirA_,X,1,&(uARAR[a1r][a*nvirA_]),1);
  }}}

  free(X);

  double **U_p_ApR = block_matrix(aoccA_*nvirA_,ndf_+3);

  C_DGEMM('N','N',aoccA_*nvirA_,ndf_+3,aoccA_*nvirA_,1.0,&(uARAR[0][0]),
    aoccA_*nvirA_,&(B_p_AR[0][0]),ndf_+3,0.0,&(U_p_ApR[0][0]),ndf_+3);

  free_block(B_p_AR);
  free_block(uARAR);

  double **B_p_RB = get_RB_ints(1);
  double **X_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);

  for(int r=0; r<nvirA_; r++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
      &(B_p_RB[r*noccB_][0]),ndf_+3,0.0,&(X_p_AR[r][0]),nvirA_*(ndf_+3));
  }

  energy = C_DDOT(aoccA_*nvirA_*(ndf_+3),&(U_p_ApR[0][0]),1,&(X_p_AR[0][0]),1);

  energy -= 2.0*C_DDOT(aoccA_*nvirA_*(ndf_+3),&(U_p_AR[0][0]),1,
    &(X_p_AR[0][0]),1);

  free_block(B_p_RB);
  free_block(X_p_AR);

  double **xAR = block_matrix(aoccA_,nvirA_);
  double **yAR = block_matrix(aoccA_,nvirA_);

  C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(xAR[0][0]),nvirA_);

  C_DGEMV('n',aoccA_*nvirA_,ndf_+3,1.0,&(U_p_ApR[0][0]),ndf_+3,diagBB_,1,
    0.0,yAR[0],1);

  energy += 2.0*C_DDOT(aoccA_*nvirA_,&(xAR[0][0]),1,
    &(yAR[0][0]),1);

  C_DGEMV('n',aoccA_*nvirA_,ndf_+3,1.0,&(U_p_AR[0][0]),(ndf_+3),diagBB_,1,
    0.0,yAR[0],1);

  energy -= 4.0*C_DDOT(aoccA_*nvirA_,&(xAR[0][0]),1,&(yAR[0][0]),1);

  free_block(xAR);
  free_block(yAR);

  double **A_p_AB = block_matrix(aoccA_*noccB_,ndf_+3);
  double **A_p_BB = block_matrix(noccB_*noccB_,ndf_+3);

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][0]),nmoB_,
      &(U_p_ApR[a*nvirA_][0]),ndf_+3,0.0,&(A_p_AB[a*noccB_][0]),ndf_+3);
  }

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),aoccA_,-1.0,&(sAB_[foccA_][0]),nmoB_,
    &(A_p_AB[0][0]),noccB_*(ndf_+3),0.0,&(A_p_BB[0][0]),noccB_*(ndf_+3));

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][0]),nmoB_,
      &(U_p_AR[a*nvirA_][0]),ndf_+3,0.0,&(A_p_AB[a*noccB_][0]),ndf_+3);
  }

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),aoccA_,2.0,&(sAB_[foccA_][0]),nmoB_,
    &(A_p_AB[0][0]),noccB_*(ndf_+3),1.0,&(A_p_BB[0][0]),noccB_*(ndf_+3));

  double **B_p_BB = get_BB_ints(1);

  energy += C_DDOT(noccB_*noccB_*(ndf_+3),&(A_p_BB[0][0]),1,&(B_p_BB[0][0]),1);

  free_block(A_p_AB);
  free_block(A_p_BB);
  free_block(U_p_AR);
  free_block(U_p_ApR);
  free_block(B_p_BB);

  return(4.0*energy);
}

double SAPT2p3::exch_disp30_02()
{
  double energy = 0.0;

  double **uBSBS = block_matrix(aoccB_*nvirB_,aoccB_*nvirB_);
  double **B_p_BS = get_BS_ints(1,foccB_);
  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  C_DGEMM('N','T',aoccB_*nvirB_,aoccB_*nvirB_,ndf_+3,1.0,&(B_p_BS[0][0]),
    ndf_+3,&(T_p_BS[0][0]),ndf_+3,0.0,&(uBSBS[0][0]),aoccB_*nvirB_);

  free_block(T_p_BS);

  for(int bs=0; bs<aoccB_*nvirB_; bs++) {
    for(int b1s1=0; b1s1<bs; b1s1++) {
      double tval = uBSBS[bs][b1s1] + uBSBS[b1s1][bs];
      uBSBS[b1s1][bs] = tval;
      uBSBS[bs][b1s1] = tval;
  }}

  C_DSCAL(aoccB_*nvirB_,2.0,&(uBSBS[0][0]),aoccB_*nvirB_+1);

  for (int b=0, bs=0; b < aoccB_; b++) {
  for (int s=0; s < nvirB_; s++, bs++) {
    for (int bb=0, bbss=0; bb < aoccB_; bb++) {
    for (int ss=0; ss < nvirB_; ss++, bbss++) {
      double denom = evalsB_[b+foccB_]+evalsB_[bb+foccB_]-
        evalsB_[s+noccB_]-evalsB_[ss+noccB_];
      uBSBS[bs][bbss] /= denom;
    }}
  }}

  double **U_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('N','N',aoccB_*nvirB_,(ndf_+3),aoccB_*nvirB_,1.0,&(uBSBS[0][0]),
    aoccB_*nvirB_,&(B_p_BS[0][0]),ndf_+3,0.0,&(U_p_BS[0][0]),ndf_+3);

  double *X = init_array(nvirB_);

  for(int b=0; b<aoccB_; b++) {
  for(int b1=0; b1<=b; b1++) {
    for(int s=0; s<nvirB_; s++) {
      int bs = b*nvirB_+s;
      int b1s = b1*nvirB_+s;
      C_DCOPY(nvirB_,&(uBSBS[bs][b1*nvirB_]),1,X,1);
      C_DCOPY(nvirB_,&(uBSBS[b1s][b*nvirB_]),1,
        &(uBSBS[bs][b1*nvirB_]),1);
      C_DCOPY(nvirB_,X,1,&(uBSBS[b1s][b*nvirB_]),1);
  }}}

  free(X);

  double **U_p_BpS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('N','N',aoccB_*nvirB_,ndf_+3,aoccB_*nvirB_,1.0,&(uBSBS[0][0]),
    aoccB_*nvirB_,&(B_p_BS[0][0]),ndf_+3,0.0,&(U_p_BpS[0][0]),ndf_+3);

  free_block(B_p_BS);
  free_block(uBSBS);

  double **B_p_AS = get_AS_ints(1);
  double **X_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('T','N',aoccB_,nvirB_*(ndf_+3),noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(B_p_AS[0][0]),nvirB_*(ndf_+3),0.0,&(X_p_BS[0][0]),nvirB_*(ndf_+3));

  energy = C_DDOT(aoccB_*nvirB_*(ndf_+3),&(U_p_BpS[0][0]),1,&(X_p_BS[0][0]),1);

  energy -= 2.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),&(U_p_BS[0][0]),1,
    &(X_p_BS[0][0]),1);

  free_block(B_p_AS);
  free_block(X_p_BS);

  double **xBS = block_matrix(aoccB_,nvirB_);
  double **yBS = block_matrix(aoccB_,nvirB_);

  C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,&(xBS[0][0]),nvirB_);

  C_DGEMV('n',aoccB_*nvirB_,ndf_+3,1.0,&(U_p_BpS[0][0]),ndf_+3,diagAA_,1,
    0.0,yBS[0],1);

  energy += 2.0*C_DDOT(aoccB_*nvirB_,&(xBS[0][0]),1,&(yBS[0][0]),1);

  C_DGEMV('n',aoccB_*nvirB_,ndf_+3,1.0,
    &(U_p_BS[0][0]),ndf_+3,diagAA_,1,0.0,yBS[0],1);

  energy -= 4.0*C_DDOT(aoccB_*nvirB_,&(xBS[0][0]),1,&(yBS[0][0]),1);

  free_block(xBS);
  free_block(yBS);

  double **A_p_BA = block_matrix(aoccB_*noccA_,ndf_+3);
  double **A_p_AA = block_matrix(noccA_*noccA_,ndf_+3);

  for(int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',noccA_,(ndf_+3),nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
      &(U_p_BpS[b*nvirB_][0]),(ndf_+3),0.0,&(A_p_BA[b*noccA_][0]),(ndf_+3));
  }

  C_DGEMM('N','N',noccA_,noccA_*(ndf_+3),aoccB_,-1.0,&(sAB_[0][foccB_]),nmoB_,
    &(A_p_BA[0][0]),noccA_*(ndf_+3),0.0,&(A_p_AA[0][0]),noccA_*(ndf_+3));

  for(int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',noccA_,(ndf_+3),nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
      &(U_p_BS[b*nvirB_][0]),(ndf_+3),0.0,&(A_p_BA[b*noccA_][0]),(ndf_+3));
  }

  C_DGEMM('N','N',noccA_,noccA_*(ndf_+3),aoccB_,2.0,&(sAB_[0][foccB_]),nmoB_,
    &(A_p_BA[0][0]),noccA_*(ndf_+3),1.0,&(A_p_AA[0][0]),noccA_*(ndf_+3));

  double **B_p_AA = get_AA_ints(1);

  energy += C_DDOT(noccA_*noccA_*(ndf_+3),&(A_p_AA[0][0]),1,&(B_p_AA[0][0]),1);

  free_block(A_p_BA);
  free_block(A_p_AA);
  free_block(U_p_BS);
  free_block(U_p_BpS);
  free_block(B_p_AA);

  return(4.0*energy);
}

double SAPT2p3::exch_disp30_22()
{
  double energy = 0.0;

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  double **tAS_RB = block_matrix(nvirA_,aoccB_);
  double **tRB_AS = block_matrix(aoccA_,nvirB_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0,bs=0; b<aoccB_; b++) {
        for (int s=0; s<nvirB_; s++,bs++) {
          tAS_RB[r][b] += tARBS[ar][bs]*sAB_[a+foccA_][s+noccB_];
          tRB_AS[a][s] += tARBS[ar][bs]*sAB_[r+noccA_][b+foccB_];
  }}}}

  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T AR Intermediates",(char *) T_p_AR[0],
    sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(PSIF_SAPT_AMPS,"T BS Intermediates",(char *) T_p_BS[0],
    sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  double **B_p_AR = get_AR_ints(0,foccA_);
  double **B_p_BS = get_BS_ints(0,foccB_);

  double **X_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  double **xBB = block_matrix(aoccB_,aoccB_);
  double **xSS = block_matrix(nvirB_,nvirB_);

  C_DGEMM('T','N',aoccB_,aoccB_,nvirA_,1.0,&(sAB_[noccA_][foccB_]),nmoB_,
    &(tAS_RB[0][0]),aoccB_,0.0,&(xBB[0][0]),aoccB_);

  C_DGEMM('N','N',aoccB_,nvirB_*(ndf_+3),aoccB_,2.0,&(xBB[0][0]),aoccB_,
    &(B_p_BS[0][0]),nvirB_*(ndf_+3),0.0,&(X_p_BS[0][0]),nvirB_*(ndf_+3));

  C_DGEMM('T','N',nvirB_,nvirB_,aoccA_,1.0,&(sAB_[foccA_][noccB_]),nmoB_,
    &(tRB_AS[0][0]),nvirB_,0.0,&(xSS[0][0]),nvirB_);

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',nvirB_,ndf_+3,nvirB_,2.0,&(xSS[0][0]),nvirB_,
      &(B_p_BS[b*nvirB_][0]),ndf_+3,1.0,&(X_p_BS[b*nvirB_][0]),ndf_+3);
  }

  energy += C_DDOT(aoccB_*nvirB_*(ndf_+3),&(T_p_BS[0][0]),1,&(X_p_BS[0][0]),1);

  free_block(xBB);
  free_block(xSS);
  free_block(X_p_BS);

  double **X_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);

  double **xAA = block_matrix(aoccA_,aoccA_);
  double **xRR = block_matrix(nvirA_,nvirA_);

  C_DGEMM('N','T',aoccA_,aoccA_,nvirB_,1.0,&(sAB_[foccA_][noccB_]),nmoB_,
    &(tRB_AS[0][0]),nvirB_,0.0,&(xAA[0][0]),aoccA_);

  C_DGEMM('N','N',aoccA_,nvirA_*(ndf_+3),aoccA_,2.0,&(xAA[0][0]),aoccA_,
    &(B_p_AR[0][0]),nvirA_*(ndf_+3),0.0,&(X_p_AR[0][0]),nvirA_*(ndf_+3));

  C_DGEMM('N','T',nvirA_,nvirA_,aoccB_,1.0,&(sAB_[noccA_][foccB_]),nmoB_,
    &(tAS_RB[0][0]),aoccB_,0.0,&(xRR[0][0]),nvirA_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',nvirA_,ndf_+3,nvirA_,2.0,&(xRR[0][0]),nvirA_,
      &(B_p_AR[a*nvirA_][0]),ndf_+3,1.0,&(X_p_AR[a*nvirA_][0]),ndf_+3);
  }

  energy += C_DDOT(aoccA_*nvirA_*(ndf_+3),&(T_p_AR[0][0]),1,&(X_p_AR[0][0]),1);

  free_block(xAA);
  free_block(xRR);
  free_block(X_p_AR);

  double **A_p_AB = block_matrix(aoccA_*aoccB_,ndf_+3);
  double **B_p_AB = block_matrix(aoccA_*aoccB_,ndf_+3);

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',aoccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][foccB_]),nmoB_,
      &(T_p_AR[a*nvirA_][0]),ndf_+3,0.0,&(A_p_AB[a*aoccB_][0]),ndf_+3);
  }

  for(int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,nvirB_,1.0,&(sAB_[foccA_][noccB_]),nmoB_,
      &(T_p_BS[b*nvirB_][0]),ndf_+3,0.0,&(B_p_AB[b][0]),aoccB_*(ndf_+3));
  }

  energy -= 4.0*C_DDOT(aoccA_*aoccB_*(ndf_+3),&(A_p_AB[0][0]),1,
    &(B_p_AB[0][0]),1);

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',aoccB_,ndf_+3,nvirA_,1.0,&(tAS_RB[0][0]),aoccB_,
      &(B_p_AR[a*nvirA_][0]),ndf_+3,0.0,&(A_p_AB[a*aoccB_][0]),ndf_+3);
  }

  for(int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,nvirB_,1.0,&(tRB_AS[0][0]),nvirB_,
      &(B_p_BS[b*nvirB_][0]),ndf_+3,0.0,&(B_p_AB[b][0]),aoccB_*(ndf_+3));
  }

  energy -= C_DDOT(aoccA_*aoccB_*(ndf_+3),&(A_p_AB[0][0]),1,&(B_p_AB[0][0]),1);

  free_block(A_p_AB);
  free_block(B_p_AB);
  free_block(T_p_AR);
  free_block(T_p_BS);

  double **tABRS = block_matrix(aoccA_*aoccB_,nvirA_*nvirB_);

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<aoccB_; b++) {
        int ab = a*aoccB_+b;
        C_DCOPY(nvirB_,&(tARBS[ar][b*nvirB_]),1,&(tABRS[ab][r*nvirB_]),1);
  }}}

  free_block(tARBS);

  double **vABRS = block_matrix(aoccA_*aoccB_,nvirA_*nvirB_);

  for(int a=0,ab=0; a<aoccA_; a++) {
    for(int b=0; b<aoccB_; b++,ab++) {
      C_DGEMM('N','T',nvirA_,nvirB_,ndf_+3,1.0,&(B_p_AR[a*nvirA_][0]),ndf_+3,
        &(B_p_BS[b*nvirB_][0]),ndf_+3,0.0,&(vABRS[ab][0]),nvirB_);
  }}

  free_block(B_p_AR);
  free_block(B_p_BS);

  double **xAR = block_matrix(aoccA_,nvirA_);
  double **ABAB = block_matrix(aoccA_*aoccB_,aoccA_*aoccB_);

  for(int a=0,ab=0; a<aoccA_; a++) {
    for(int b=0; b<aoccB_; b++,ab++) {
      C_DGEMM('N','T',aoccA_,nvirA_,nvirB_,1.0,&(sAB_[foccA_][noccB_]),nmoB_,
        &(tABRS[ab][0]),nvirB_,0.0,&(xAR[0][0]),nvirA_);
      C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,&(xAR[0][0]),nvirA_,
        &(sAB_[noccA_][foccB_]),nmoB_,0.0,&(ABAB[ab][0]),aoccB_);
  }}

  free_block(xAR);

  double **xABRS = block_matrix(aoccA_*aoccB_,nvirA_*nvirB_);

  C_DGEMM('T','N',aoccA_*aoccB_,nvirA_*nvirB_,aoccA_*aoccB_,1.0,&(ABAB[0][0]),
    aoccA_*aoccB_,&(vABRS[0][0]),nvirA_*nvirB_,0.0,&(xABRS[0][0]),
    nvirA_*nvirB_);

  energy -= C_DDOT((long int) aoccA_*aoccB_*nvirA_*nvirB_,&(tABRS[0][0]),1,
    &(xABRS[0][0]),1);

  free_block(tABRS);
  free_block(ABAB);
  free_block(vABRS);
  free_block(xABRS);
  free_block(tAS_RB);
  free_block(tRB_AS);

  return(2.0*energy);
}

}}
