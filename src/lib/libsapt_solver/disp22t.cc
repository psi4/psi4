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

void SAPT2p::disp22tccd()
{
  if (print_) {
    fprintf(outfile,"\n");
  }

  if (nat_orbs_) {
    natural_orbitalify_ccd();
  }

  double e_disp220t;


  if (nat_orbs_) {
    e_disp220t = disp220tccd(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals","RR NO RI Integrals",PSIF_SAPT_BB_DF_INTS,
      "BS NO RI Integrals",PSIF_SAPT_CCD,"T ARAR Natorb Amplitudes","T BSAR Natorb Amplitudes", 
      no_evalsA_,no_evalsB_,noccA_,no_nvirA_,foccA_,noccB_,no_nvirB_,foccB_);
  }
  else {
    e_disp220t = disp220tccd(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      PSIF_SAPT_AA_DF_INTS,"AR RI Integrals","RR RI Integrals",PSIF_SAPT_BB_DF_INTS,
      "BS RI Integrals",PSIF_SAPT_CCD,"T ARAR Amplitudes","T BSAR Amplitudes", 
      evalsA_,evalsB_,noccA_,nvirA_,foccA_,noccB_,nvirB_,foccB_);
  }

  if (print_) {
    fprintf(outfile,"\n    Disp220 (T)         = %18.12lf H\n\n",e_disp220t);
    fflush(outfile);
  }

  double e_disp202t;

  if (nat_orbs_) {
    e_disp220t = disp220tccd(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals","SS NO RI Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR NO RI Integrals",PSIF_SAPT_CCD,"T BSBS Natorb Amplitudes","T ARBS Natorb Amplitudes", 
      no_evalsB_,no_evalsA_,noccB_,no_nvirB_,foccB_,noccA_,no_nvirA_,foccA_);
  }
  else {
    e_disp220t = disp220tccd(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      PSIF_SAPT_BB_DF_INTS,"BS RI Integrals","SS RI Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR RI Integrals",PSIF_SAPT_CCD,"T BSBS Amplitudes","T ARBS Amplitudes", 
      evalsB_,evalsA_,noccB_,nvirB_,foccB_,noccA_,nvirA_,foccA_);
  }

  if (print_) {
    fprintf(outfile,"\n    Disp202 (T)         = %18.12lf H\n\n",e_disp202t);
    fflush(outfile);
  }

  e_disp22t_ccd_ = e_disp220t + e_disp202t;

  if (print_) {
    fprintf(outfile,"    Disp22 (T)          = %18.12lf H\n",e_disp22t_ccd_);
    fflush(outfile);
  }

  if (nat_orbs_) {
    double scale = e_disp20_/e_no_disp20_;
    e_disp220t *= scale;
    e_disp202t *= scale;
    e_est_disp22t_ccd_ = e_disp220t + e_disp202t;

    if (print_) {
      fprintf(outfile,"\n    Est. Disp220 (T)    = %18.12lf H\n",e_disp220t);
      fprintf(outfile,"    Est. Disp202 (T)    = %18.12lf H\n\n",e_disp202t);
      fprintf(outfile,"    Est. Disp22 (T)     = %18.12lf H\n",e_est_disp22t_ccd_);
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

  time_t start = time(NULL);
  time_t stop;

  for(int b=0,bs=0; b<aoccB; b++) {
  for(int s=0; s<nvirB; s++,bs++) {

    psio_address next_DF_BS = psio_get_address(PSIO_ZERO,
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

double SAPT2p::disp220tccd(int AAnum, char *AA_label, int Rnum, char *AR_label, 
  char *RR_label, int BBnum, char *BS_label, int ampnum, char *tarar, 
  char *tbsar, double *evalsA, double *evalsB, int noccA, int nvirA, 
  int foccA, int noccB, int nvirB, int foccB)
{
  double energy = 0.0;

  noccA -= foccA;
  noccB -= foccB;

  double **w_ARAR = block_matrix(noccA*nvirA,noccA*nvirA);

  double **v_bsAA = block_matrix(noccA,noccA);
  double **v_bsRR = block_matrix(nvirA,nvirA);
  double **v_ARAA = block_matrix(noccA*nvirA,noccA*noccA);

  double **B_p_AA = get_DF_ints_nongimp(AAnum,AA_label,foccA,noccA+foccA,
    foccA,noccA+foccA);
  double **B_p_AR = get_DF_ints_nongimp(Rnum,AR_label,foccA,noccA+foccA,0,nvirA);
  double **B_p_RR = get_DF_ints_nongimp(Rnum,RR_label,0,nvirA,0,nvirA);
  double *B_p_bs = init_array(ndf_);

  double **t_bsAR = block_matrix(noccA,nvirA);
  double **t_ARAR;

  psio_address next_ARAR;

  if (ampnum == PSIF_SAPT_CCD)  {
    t_ARAR = block_matrix(noccA*nvirA,noccA*nvirA);
    psio_->read_entry(ampnum,tarar,(char *) t_ARAR[0],noccA*nvirA*noccA*
      nvirA*(ULI) sizeof(double));
  }
  else if (ampnum)  {
    t_ARAR = block_matrix(noccA*nvirA,noccA*nvirA);
    for(int a=0,ar=0; a<noccA; a++) {
    for(int r=0; r<nvirA; r++,ar++) {
      next_ARAR = psio_get_address(PSIO_ZERO,
        ((foccA*nvirA+ar)*(noccA+foccA)*nvirA+foccA*nvirA)*sizeof(double));
      psio_->read(ampnum,tarar,(char *) t_ARAR[ar],noccA*nvirA*
        (ULI) sizeof(double),next_ARAR,&next_ARAR);
    }}
  }
  else {
    t_ARAR = block_matrix(noccA*nvirA,noccA*nvirA);
    C_DGEMM('N','T',noccA*nvirA,noccA*nvirA,ndf_,1.0,
      &(B_p_AR[0][0]),ndf_,&(B_p_AR[0][0]),ndf_,0.0,
      &(t_ARAR[0][0]),noccA*nvirA);

    for(int a=0,ar=0; a<noccA; a++) {
    for(int r=0; r<nvirA; r++,ar++) {
      for(int a1=0,a1r1=0; a1<noccA; a1++) {
      for(int r1=0; r1<nvirA; r1++,a1r1++) {
        double denom = evalsA[a+foccA]+evalsA[a1+foccA]-
                  evalsA[r+noccA+foccA]-evalsA[r1+noccA+foccA];
        t_ARAR[ar][a1r1] /= denom;
      }}
    }} 
  }

  double **C_p_AR = block_matrix(noccA*nvirA,ndf_);

  C_DGEMM('N','T',noccA*nvirA,noccA*noccA,ndf_,1.0,&(B_p_AR[0][0]),
    ndf_,&(B_p_AA[0][0]),ndf_,0.0,&(v_ARAA[0][0]),
    noccA*noccA);

  psio_address next_BSAR;

  time_t start = time(NULL);
  time_t stop;
  
  for(int b=0,bs=0; b<noccB; b++) {
  for(int s=0; s<nvirB; s++,bs++) {
  
    psio_address next_DF_BS = psio_get_address(PSIO_ZERO,((foccB + b)*nvirB + s)*
      (ndf_+3)*(ULI) sizeof(double));
    psio_->read(BBnum,BS_label,(char *) &(B_p_bs[0]),sizeof(double)*
      ndf_,next_DF_BS,&next_DF_BS);

    if (ampnum == PSIF_SAPT_CCD) {
      next_BSAR = psio_get_address(PSIO_ZERO,bs*noccA*nvirA*sizeof(double));
      psio_->read(ampnum,tbsar,(char *) t_bsAR[0],sizeof(double)*
        noccA*nvirA,next_BSAR,&next_BSAR);
    }
    else if (ampnum) {
      next_BSAR = psio_get_address(PSIO_ZERO,
        ((foccB*nvirB+bs)*(noccA+foccA)*nvirA+foccA*nvirA)*sizeof(double));
      psio_->read(ampnum,tbsar,(char *) t_bsAR[0],sizeof(double)*
        noccA*nvirA,next_BSAR,&next_BSAR);
    }
    else {
      C_DGEMV('n',noccA*nvirA,ndf_,1.0,B_p_AR[0],ndf_,
        B_p_bs,1,0.0,t_bsAR[0],1);

      for(int a=0; a<noccA; a++) {
      for(int r=0; r<nvirA; r++) {
        double denom = evalsA[a+foccA]+evalsB[b+foccB]-
                  evalsA[r+noccA+foccA]-evalsB[s+noccB+foccB];
        t_bsAR[a][r] /= denom;
      }}
    }

    C_DGEMV('n',noccA*noccA,ndf_,1.0,B_p_AA[0],ndf_,
      B_p_bs,1,0.0,v_bsAA[0],1);
    C_DGEMV('n',nvirA*nvirA,ndf_,1.0,B_p_RR[0],ndf_,
      B_p_bs,1,0.0,v_bsRR[0],1);

    C_DGEMM('N','N',noccA*nvirA*noccA,nvirA,nvirA,1.0,&(t_ARAR[0][0]),nvirA,
            &(v_bsRR[0][0]),nvirA,0.0,&(w_ARAR[0][0]),nvirA);
    C_DGEMM('N','N',noccA,nvirA*noccA*nvirA,noccA,-1.0,&(v_bsAA[0][0]),noccA,
            &(t_ARAR[0][0]),nvirA*noccA*nvirA,1.0,&(w_ARAR[0][0]),
            nvirA*noccA*nvirA);
    C_DGEMM('N','N',noccA*nvirA*noccA,nvirA,noccA,-1.0,&(v_ARAA[0][0]),noccA,
            &(t_bsAR[0][0]),nvirA,1.0,&(w_ARAR[0][0]),nvirA);
    C_DGEMM('N','N',noccA,nvirA*ndf_,nvirA,1.0,&(t_bsAR[0][0]),
            nvirA,&(B_p_RR[0][0]),nvirA*ndf_,0.0,&(C_p_AR[0][0]),
            nvirA*ndf_);
    C_DGEMM('N','T',noccA*nvirA,noccA*nvirA,ndf_,1.0,
            &(B_p_AR[0][0]),ndf_,&(C_p_AR[0][0]),ndf_,
            1.0,&(w_ARAR[0][0]),noccA*nvirA);

    for(int a=0,ar=0; a<noccA; a++) {
    for(int r=0; r<nvirA; r++,ar++) {
      for(int a1=0,a1r1=0; a1<noccA; a1++) {
      for(int r1=0; r1<nvirA; r1++,a1r1++) {
        int a1r = a1*nvirA+r;
        int ar1 = a*nvirA+r1;
        double tval1 = w_ARAR[ar][a1r1]+w_ARAR[a1r1][ar];
        double tval2 = w_ARAR[a1r][ar1]+w_ARAR[ar1][a1r];
        double denom = evalsA[a+foccA]+evalsA[a1+foccA]+evalsB[b+foccB]-
                  evalsA[r+noccA+foccA]-evalsA[r1+noccA+foccA]-
                  evalsB[s+noccB+foccB];
        energy += ((4.0*tval1-2.0*tval2)*tval1)/denom;
      }}
    }}
   }
  stop = time(NULL);
    fprintf(outfile,"    (i = %3d of %3d) %10ld seconds\n",b+1,noccB,stop-start);
  fflush(outfile);
  }

  free(B_p_bs);
  free_block(w_ARAR);
  free_block(v_bsAA);
  free_block(v_bsRR);
  free_block(v_ARAA);
  free_block(t_ARAR);
  free_block(t_bsAR);
  free_block(B_p_AA);
  free_block(B_p_AR);
  free_block(B_p_RR);
  free_block(C_p_AR);

  return(energy);
}

}}
