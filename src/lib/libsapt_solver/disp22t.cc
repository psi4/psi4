#include "sapt2p.h"

namespace psi { namespace sapt {

void SAPT2p::disp22t()
{
  if (print_) {
    fprintf(outfile,"\n");
  }

  double e_disp220t;

  if (nat_orbs_) {
    e_disp220t = disp220t(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR NO RI Integrals","RR NO RI Integrals",PSIF_SAPT_BB_DF_INTS,
      "BS NO RI Integrals",PSIF_SAPT_AMPS,"tARAR NO Amplitudes",foccA_,
      noccA_,no_nvirA_,foccB_,noccB_,no_nvirB_,no_evalsA_,no_evalsB_);
  }
  else {
    e_disp220t = disp220t(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR RI Integrals","RR RI Integrals",PSIF_SAPT_BB_DF_INTS,
      "BS RI Integrals",PSIF_SAPT_AMPS,"tARAR Amplitudes",foccA_,
      noccA_,nvirA_,foccB_,noccB_,nvirB_,evalsA_,evalsB_);
  }

  if (print_) {
    fprintf(outfile,"\n    Disp220 (T)         = %18.12lf H\n\n",e_disp220t);
    fflush(outfile);
  }

  double e_disp202t;

  if (nat_orbs_) {
    e_disp202t = disp220t(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS NO RI Integrals","SS NO RI Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR NO RI Integrals",PSIF_SAPT_AMPS,"tBSBS NO Amplitudes",foccB_,
      noccB_,no_nvirB_,foccA_,noccA_,no_nvirA_,no_evalsB_,no_evalsA_);
  }
  else {
    e_disp202t = disp220t(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS RI Integrals","SS RI Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR RI Integrals",PSIF_SAPT_AMPS,"tBSBS Amplitudes",foccB_,
      noccB_,nvirB_,foccA_,noccA_,nvirA_,evalsB_,evalsA_);
  }

  if (print_) {
    fprintf(outfile,"\n    Disp202 (T)         = %18.12lf H\n\n",e_disp202t);
    fflush(outfile);
  }

  e_disp22t_ = e_disp220t + e_disp202t;

  if (print_) {
    fprintf(outfile,"    Disp22 (T)          = %18.12lf H\n",e_disp22t_);
    fflush(outfile);
  }

  if (nat_orbs_) {
    double scale = e_disp20_/e_no_disp20_;
    e_disp220t *= scale;
    e_disp202t *= scale;
    e_est_disp22t_ = e_disp220t + e_disp202t;

    if (print_) {
      fprintf(outfile,"\n    Est. Disp220 (T)    = %18.12lf H\n",e_disp220t);
      fprintf(outfile,"    Est. Disp202 (T)    = %18.12lf H\n\n",e_disp202t);
      fprintf(outfile,"    Est. Disp22 (T)     = %18.12lf H\n",e_est_disp22t_);
      fflush(outfile);
    }
  }
}

double SAPT2p::disp220t(int AAfile, const char *AAlabel, const char *ARlabel,
  const char *RRlabel, int BBfile, const char *BSlabel, int ampfile, 
  const char *tlabel, int foccA, int noccA, int nvirA, int foccB, int noccB, 
  int nvirB, double *evalsA, double *evalsB)
{
  double energy = 0.0;

  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **wARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);

  double **vbsAA = block_matrix(aoccA,aoccA);
  double **vbsRR = block_matrix(nvirA,nvirA);
  double **vARAA = block_matrix(aoccA*nvirA,aoccA*aoccA);

  double **tARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  double **tbsAR = block_matrix(aoccA,nvirA);

  double **B_p_AA = get_DF_ints(AAfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_AR = get_DF_ints(AAfile,ARlabel,foccA,noccA,0,nvirA);
  double **B_p_RR = get_DF_ints(AAfile,RRlabel,0,nvirA,0,nvirA);
  double *B_p_bs = init_array(ndf_+3);

  double **C_p_AR = block_matrix(aoccA*nvirA,ndf_+3);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*aoccA,ndf_+3,1.0,&(B_p_AR[0][0]),
    ndf_+3,&(B_p_AA[0][0]),ndf_+3,0.0,&(vARAA[0][0]),aoccA*aoccA);

  psio_address next_DF_BS = PSIO_ZERO;

  time_t start = time(NULL);
  time_t stop;

  for(int b=0,bs=0; b<aoccB; b++) {
  for(int s=0; s<nvirB; s++,bs++) {

    next_DF_BS = psio_get_address(PSIO_ZERO,
      sizeof(double)*(b+foccB)*nvirB*(ndf_+3) + sizeof(double)*s*(ndf_+3));
    psio_->read(BBfile,BSlabel,(char *) &(B_p_bs[0]),sizeof(double)*(ndf_+3),
      next_DF_BS,&next_DF_BS);

    C_DGEMV('n',aoccA*nvirA,ndf_+3,1.0,B_p_AR[0],ndf_+3,B_p_bs,1,
      0.0,tbsAR[0],1);

    for(int a=0,ar=0; a<aoccA; a++) {
    for(int r=0; r<nvirA; r++,ar++) {
      double denom = evalsA[a+foccA]+evalsB[b+foccB]
        -evalsA[r+noccA]-evalsB[s+noccB];
      tbsAR[a][r] /= denom;
    }}

    C_DGEMV('n',aoccA*aoccA,ndf_+3,1.0,B_p_AA[0],ndf_+3,B_p_bs,1,
      0.0,vbsAA[0],1);
    C_DGEMV('n',nvirA*nvirA,ndf_+3,1.0,B_p_RR[0],ndf_+3,B_p_bs,1,
      0.0,vbsRR[0],1);

    C_DGEMM('N','N',aoccA*nvirA*aoccA,nvirA,nvirA,1.0,&(tARAR[0][0]),nvirA,
      &(vbsRR[0][0]),nvirA,0.0,&(wARAR[0][0]),nvirA);
    C_DGEMM('N','N',aoccA,nvirA*aoccA*nvirA,aoccA,-1.0,&(vbsAA[0][0]),aoccA,
      &(tARAR[0][0]),nvirA*aoccA*nvirA,1.0,&(wARAR[0][0]),nvirA*aoccA*nvirA);
    C_DGEMM('N','N',aoccA*nvirA*aoccA,nvirA,aoccA,-1.0,&(vARAA[0][0]),aoccA,
      &(tbsAR[0][0]),nvirA,1.0,&(wARAR[0][0]),nvirA);
    C_DGEMM('N','N',aoccA,nvirA*(ndf_+3),nvirA,1.0,&(tbsAR[0][0]),nvirA,
      &(B_p_RR[0][0]),nvirA*(ndf_+3),0.0,&(C_p_AR[0][0]),nvirA*(ndf_+3));
    C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,&(B_p_AR[0][0]),ndf_+3,
      &(C_p_AR[0][0]),ndf_+3,1.0,&(wARAR[0][0]),aoccA*nvirA);

    for(int a=0,ar=0; a<aoccA; a++) {
    for(int r=0; r<nvirA; r++,ar++) {
      for(int a1=0,a1r1=0; a1<aoccA; a1++) {
      for(int r1=0; r1<nvirA; r1++,a1r1++) {
        int a1r = a1*nvirA+r;
        int ar1 = a*nvirA+r1;
        double tval1 = wARAR[ar][a1r1]+wARAR[a1r1][ar];
        double tval2 = wARAR[a1r][ar1]+wARAR[ar1][a1r];
        double denom = evalsA[a+foccA]+evalsA[a1+foccA]+evalsB[b+foccB]-
                  evalsA[r+noccA]-evalsA[r1+noccA]-evalsB[s+noccB];
        energy += ((4.0*tval1-2.0*tval2)*tval1)/denom;
      }}
    }}
   }
  stop = time(NULL);
  if (print_) {
    fprintf(outfile,"    (i = %3d of %3d) %10ld seconds\n",b+1,aoccB,
      stop-start);
  } fflush(outfile); }

  free(B_p_bs);
  free_block(wARAR);
  free_block(vbsAA);
  free_block(vbsRR);
  free_block(vARAA);
  free_block(tARAR);
  free_block(tbsAR);
  free_block(B_p_AA);
  free_block(B_p_AR);
  free_block(B_p_RR);
  free_block(C_p_AR);

  return(energy);
}

}}
