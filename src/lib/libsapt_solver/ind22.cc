#include "sapt2.h"

using namespace boost;

namespace psi { namespace sapt {

void SAPT2::ind22()
{
  double e_ind220 = ind220();

  if (debug_) {
    fprintf(outfile,"    Ind220              = %18.12lf H\n",e_ind220);
    fflush(outfile);
  }

  double e_ind202 = ind202();

  if (debug_) {
    fprintf(outfile,"    Ind202              = %18.12lf H\n\n",e_ind202);
    fflush(outfile);
  }

  e_ind22_ = e_ind220 + e_ind202;
  e_exch_ind22_ = e_ind22_*(e_exch_ind20_/e_ind20_);

  if (print_) {
    fprintf(outfile,"    Ind22               = %18.12lf H\n",e_ind22_);
    fflush(outfile);
  }
}

double SAPT2::ind220()
{
  double **iAR = block_matrix(aoccA_,nvirA_);

  for (int a=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++) {
      iAR[a][r] = wBAR_[a+foccA_][r]/(evalsA_[a+foccA_] - evalsA_[r+noccA_]);
  }}

  double **iBS = block_matrix(aoccB_,nvirB_);

  for (int b=0; b<aoccB_; b++) {
    for (int s=0; s<nvirB_; s++) {
      iBS[b][s] = wABS_[b+foccB_][s]/(evalsB_[b+foccB_] - evalsB_[s+noccB_]);
  }}

  double energy = 0.0;

  energy += ind220_1(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_AMPS,"tARAR Amplitudes",iAR,wBAA_,wBRR_,
    foccA_,noccA_,nvirA_,evalsA_);

  energy += ind220_2(PSIF_SAPT_AMPS,"T2 AR Amplitudes",iAR,wBAA_,wBRR_,
    foccA_,noccA_,nvirA_);

  energy += ind220_3(PSIF_SAPT_AMPS,"pAA Density Matrix","pRR Density Matrix",
    iAR,wBAR_,foccA_,noccA_,nvirA_);

  energy += ind220_4(PSIF_SAPT_AMPS,"Theta AR Intermediates",
    PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",iAR,foccA_,noccA_,nvirA_);

  energy += ind220_5(PSIF_SAPT_AMPS,"t2ARAR Amplitudes",iAR,foccA_,noccA_,
    nvirA_,evalsA_);

  energy += ind220_6(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_AMPS,"tARAR Amplitudes",iAR,foccA_,
    noccA_,nvirA_);

  energy += ind220_7(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_SAPT_AMPS,
    "T2 AR Amplitudes","pAA Density Matrix","pRR Density Matrix",iBS,
    foccA_,noccA_,nvirA_,foccB_,noccB_,nvirB_);

  free_block(iAR);
  free_block(iBS);

  return(energy);
}

double SAPT2::ind202()
{
  double **iAR = block_matrix(aoccA_,nvirA_);

  for (int a=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++) {
      iAR[a][r] = wBAR_[a+foccA_][r]/(evalsA_[a+foccA_] - evalsA_[r+noccA_]);
  }}

  double **iBS = block_matrix(aoccB_,nvirB_);

  for (int b=0; b<aoccB_; b++) {
    for (int s=0; s<nvirB_; s++) {
      iBS[b][s] = wABS_[b+foccB_][s]/(evalsB_[b+foccB_] - evalsB_[s+noccB_]);
  }}

  double energy = 0.0;

  energy += ind220_1(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AMPS,"tBSBS Amplitudes",iBS,wABB_,wASS_,
    foccB_,noccB_,nvirB_,evalsB_);

  energy += ind220_2(PSIF_SAPT_AMPS,"T2 BS Amplitudes",iBS,wABB_,wASS_,
    foccB_,noccB_,nvirB_);

  energy += ind220_3(PSIF_SAPT_AMPS,"pBB Density Matrix","pSS Density Matrix",
    iBS,wABS_,foccB_,noccB_,nvirB_);

  energy += ind220_4(PSIF_SAPT_AMPS,"Theta BS Intermediates",
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",iBS,foccB_,noccB_,nvirB_);

  energy += ind220_5(PSIF_SAPT_AMPS,"t2BSBS Amplitudes",iBS,foccB_,noccB_,
    nvirB_,evalsB_);

  energy += ind220_6(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AMPS,"tBSBS Amplitudes",iBS,foccB_,
    noccB_,nvirB_);

  energy += ind220_7(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_SAPT_AMPS,
    "T2 BS Amplitudes","pBB Density Matrix","pSS Density Matrix",iAR,
    foccB_,noccB_,nvirB_,foccA_,noccA_,nvirA_);

  free_block(iAR);
  free_block(iBS);

  return(energy);
}

double SAPT2::ind220_1(int intfile, const char *AAlabel, const char *ARlabel, 
  const char *RRlabel, int ampfile, const char *tlabel, double **iAR,
  double **wBAA, double **wBRR, int foccA, int noccA, int nvirA, 
  double *evalsA)
{
  int aoccA = noccA - foccA;

  double **C_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  C_DGEMM('N','N',aoccA,nvirA*(ndf_+3),nvirA,1.0,iAR[0],nvirA,B_p_RR[0],
    nvirA*(ndf_+3),0.0,C_p_AR[0],nvirA*(ndf_+3));
 
  free_block(B_p_RR);
 
  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);

  for (int a=0; a<aoccA; a++) {
  C_DGEMM('T','N',nvirA,ndf_+3,aoccA,-1.0,iAR[0],nvirA,B_p_AA[a*aoccA],ndf_+3,
    1.0,C_p_AR[a*nvirA],ndf_+3);
  }

  free_block(B_p_AA);

  double **xARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,C_p_AR[0],ndf_+3,
    B_p_AR[0],ndf_+3,0.0,xARAR[0],aoccA*nvirA);

  free_block(B_p_AR);
  free_block(C_p_AR);

  double **tARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);

  C_DGEMM('N','N',aoccA,nvirA*aoccA*nvirA,aoccA,-1.0,&(wBAA[foccA][foccA]),
    noccA,tARAR[0],nvirA*aoccA*nvirA,1.0,xARAR[0],nvirA*aoccA*nvirA);

  C_DGEMM('N','T',aoccA*nvirA*aoccA,nvirA,nvirA,1.0,tARAR[0],nvirA,wBRR[0],
    nvirA,1.0,xARAR[0],nvirA);

  free_block(tARAR);

  symmetrize(xARAR[0],aoccA,nvirA);

  double **yARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  C_DCOPY((long int) aoccA*nvirA*aoccA*nvirA,xARAR[0],1,yARAR[0],1);
  antisym(yARAR,aoccA,nvirA);

  for (int a=0, ar=0; a<aoccA; a++){
    for (int r=0; r<nvirA; r++, ar++){
      for (int aa=0, aarr=0; aa<aoccA; aa++){
        for (int rr=0; rr<nvirA; rr++, aarr++){
          xARAR[ar][aarr] /= evalsA[a+foccA]+evalsA[aa+foccA]-evalsA[r+noccA]-
            evalsA[rr+noccA];
  }}}}

  double energy = C_DDOT((long int) aoccA*nvirA*aoccA*nvirA,xARAR[0],1,
    yARAR[0],1);

  free_block(xARAR);
  free_block(yARAR);

  if (debug_) {
    fprintf(outfile,"\n    Ind22_1             = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::ind220_2(int ampfile, const char *tlabel, double **iAR,
  double **wBAA, double **wBRR, int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **tAR = block_matrix(aoccA,nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tAR[0],sizeof(double)*aoccA*nvirA);

  double **zAR = block_matrix(aoccA,nvirA);

  C_DGEMM('N','T',aoccA,nvirA,nvirA,1.0,iAR[0],nvirA,wBRR[0],nvirA,
    0.0,zAR[0],nvirA);

  C_DGEMM('N','N',aoccA,nvirA,aoccA,-1.0,&(wBAA[foccA][foccA]),noccA,
    iAR[0],nvirA,1.0,zAR[0],nvirA);

  double energy = 4.0*C_DDOT((long int) aoccA*nvirA,tAR[0],1,zAR[0],1);

  free_block(tAR);
  free_block(zAR);

  if (debug_) {
    fprintf(outfile,"    Ind22_2             = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::ind220_3(int ampfile, const char *AAlabel, const char *RRlabel,
  double **iAR, double **wBAR, int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **pAA = block_matrix(aoccA,aoccA);
  double **pRR = block_matrix(nvirA,nvirA);

  psio_->read_entry(ampfile,AAlabel,(char *) pAA[0],
    sizeof(double)*aoccA*aoccA);
  psio_->read_entry(ampfile,RRlabel,(char *) pRR[0],
    sizeof(double)*nvirA*nvirA);

  double **xAA = block_matrix(aoccA,aoccA);
  double **xRR = block_matrix(nvirA,nvirA);

  C_DGEMM('N','T',aoccA,aoccA,nvirA,1.0,iAR[0],nvirA,wBAR[foccA],nvirA,
    0.0,xAA[0],aoccA); 
  C_DGEMM('T','N',nvirA,nvirA,aoccA,1.0,iAR[0],nvirA,wBAR[foccA],nvirA,
    0.0,xRR[0],nvirA); 

  double energy = 0.0;

  energy -= 2.0*C_DDOT(aoccA*aoccA,pAA[0],1,xAA[0],1);
  energy -= 2.0*C_DDOT(nvirA*nvirA,pRR[0],1,xRR[0],1);

  free_block(pAA);
  free_block(pRR);
  free_block(xAA);
  free_block(xRR);

  if (debug_) {
    fprintf(outfile,"    Ind22_3             = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::ind220_4(int ampfile, const char *thetalabel, int intfile, 
  const char *ARlabel, double **iAR, int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **xAA = block_matrix(aoccA,aoccA);
  double **xRR = block_matrix(nvirA,nvirA);

  C_DGEMM('N','T',aoccA,aoccA,nvirA,1.0,iAR[0],nvirA,iAR[0],nvirA,
    0.0,xAA[0],aoccA);
  C_DGEMM('T','N',nvirA,nvirA,aoccA,1.0,iAR[0],nvirA,iAR[0],nvirA,
    0.0,xRR[0],nvirA);

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **C_p_AR = block_matrix(aoccA*nvirA,ndf_+3);

  C_DGEMM('N','N',aoccA,nvirA*(ndf_+3),aoccA,1.0,xAA[0],aoccA,B_p_AR[0],
    nvirA*(ndf_+3),0.0,C_p_AR[0],nvirA*(ndf_+3));

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','N',nvirA,ndf_+3,nvirA,1.0,xRR[0],nvirA,B_p_AR[a*nvirA],
      ndf_+3,1.0,C_p_AR[a*nvirA],ndf_+3);
  }

  free_block(xAA);
  free_block(xRR);
  free_block(B_p_AR);

  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double energy = -2.0*C_DDOT(aoccA*nvirA*(ndf_+3),C_p_AR[0],1,T_p_AR[0],1);

  free_block(C_p_AR);
  free_block(T_p_AR);

  if (debug_) {
    fprintf(outfile,"    Ind22_4             = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::ind220_5(int ampfile, const char *tlabel, double **iAR, 
  int foccA, int noccA, int nvirA, double *evalsA)
{
  int aoccA = noccA - foccA;

  double **tARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  antisym(tARAR,aoccA,nvirA);

  for (int a=0, ar=0; a<aoccA; a++){
    for (int r=0; r<nvirA; r++, ar++){
      for (int aa=0, aarr=0; aa<aoccA; aa++){
        for (int rr=0; rr<nvirA; rr++, aarr++){
          tARAR[ar][aarr] *= evalsA[a+foccA]+evalsA[aa+foccA]-evalsA[r+noccA]-
            evalsA[rr+noccA];
  }}}}

  double **xAR = block_matrix(aoccA,nvirA);

  C_DGEMV('n',aoccA*nvirA,aoccA*nvirA,1.0,tARAR[0],aoccA*nvirA,iAR[0],1,
    0.0,xAR[0],1);

  double energy = 2.0*C_DDOT(aoccA*nvirA,xAR[0],1,iAR[0],1);

  free_block(xAR);
  free_block(tARAR);

  if (debug_) {
    fprintf(outfile,"    Ind22_5             = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::ind220_6(int intfile, const char *AAlabel, const char *ARlabel,
  const char *RRlabel, int ampfile, const char *tlabel, double **iAR, 
  int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **gARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,2.0,B_p_AR[0],ndf_+3,
    B_p_AR[0],ndf_+3,0.0,gARAR[0],aoccA*nvirA);

  free_block(B_p_AR);

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  for (int a=0, ar=0; a<aoccA; a++){
    for (int r=0; r<nvirA; r++, ar++){
      C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-1.0,B_p_AA[a*aoccA],ndf_+3,
        B_p_RR[r*nvirA],ndf_+3,1.0,gARAR[ar],nvirA);
  }}

  free_block(B_p_AA);
  free_block(B_p_RR);

  double **xAR = block_matrix(aoccA,nvirA);
  double **yAR = block_matrix(aoccA,nvirA);

  C_DGEMV('n',aoccA*nvirA,aoccA*nvirA,1.0,gARAR[0],aoccA*nvirA,iAR[0],1,
    0.0,xAR[0],1);

  free_block(gARAR);

  double **tARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  antisym(tARAR,aoccA,nvirA);

  C_DGEMV('n',aoccA*nvirA,aoccA*nvirA,1.0,tARAR[0],aoccA*nvirA,iAR[0],1,
    0.0,yAR[0],1);

  free_block(tARAR);

  double energy = -4.0*C_DDOT(aoccA*nvirA,xAR[0],1,yAR[0],1);

  free_block(xAR);
  free_block(yAR);

  if (debug_) {
    fprintf(outfile,"    Ind22_6             = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

double SAPT2::ind220_7(int AAfile, const char *AAlabel, const char *ARlabel,
  const char *RRlabel, int BBfile, const char *BSlabel, int ampfile, 
  const char *tlabel, const char *pAAlabel, const char *pRRlabel, 
  double **iBS, int foccA, int noccA, int nvirA, int foccB, int noccB, 
  int nvirB)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **pAA = block_matrix(aoccA,aoccA);
  double **tAR = block_matrix(aoccA,nvirA);
  double **pRR = block_matrix(nvirA,nvirA);

  psio_->read_entry(ampfile,pAAlabel,(char *) pAA[0],
    sizeof(double)*aoccA*aoccA);
  psio_->read_entry(ampfile,tlabel,(char *) tAR[0],
    sizeof(double)*aoccA*nvirA);
  psio_->read_entry(ampfile,pRRlabel,(char *) pRR[0],
    sizeof(double)*nvirA*nvirA);

  double *W = init_array(ndf_+3);
  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);
  double *Z = init_array(ndf_+3);

  double **B_p_AA = get_DF_ints(AAfile,AAlabel,foccA,noccA,foccA,noccA);

  C_DGEMV('t',aoccA*aoccA,ndf_+3,1.0,B_p_AA[0],ndf_+3,pAA[0],1,0.0,W,1);

  free_block(B_p_AA); 

  double **B_p_RR = get_DF_ints(AAfile,RRlabel,0,nvirA,0,nvirA);

  C_DGEMV('t',nvirA*nvirA,ndf_+3,1.0,B_p_RR[0],ndf_+3,pRR[0],1,0.0,X,1);

  free_block(B_p_RR); 

  double **B_p_AR = get_DF_ints(AAfile,ARlabel,foccA,noccA,0,nvirA);

  C_DGEMV('t',aoccA*nvirA,ndf_+3,1.0,B_p_AR[0],ndf_+3,tAR[0],1,0.0,Y,1);

  free_block(B_p_AR); 

  double **B_p_BS = get_DF_ints(BBfile,BSlabel,foccB,noccB,0,nvirB);

  C_DGEMV('t',aoccB*nvirB,ndf_+3,1.0,B_p_BS[0],ndf_+3,iBS[0],1,0.0,Z,1);

  free_block(B_p_BS); 

  double energy = 0.0;

  energy -= 8.0*C_DDOT(ndf_+3,W,1,Z,1);
  energy += 8.0*C_DDOT(ndf_+3,X,1,Z,1);
  energy += 16.0*C_DDOT(ndf_+3,Y,1,Z,1);

  free(W);
  free(X);
  free(Y);
  free(Z);
  free_block(pAA);
  free_block(pRR);
  free_block(tAR);

  if (debug_) {
    fprintf(outfile,"    Ind22_7             = %18.12lf H\n",energy);
    fflush(outfile);
  }

  return(energy);
}

}}

