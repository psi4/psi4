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

void SAPT2p::disp22t()
{
  if (print_) {
    outfile->Printf("\n");
  }

  double e_disp220t;

  if (nat_orbs_t3_) {
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
    outfile->Printf("\n    Disp220 (T)         = %18.12lf [Eh]\n\n",e_disp220t);
    
  }

  double e_disp202t;

  if (nat_orbs_t3_) {
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
    outfile->Printf("\n    Disp202 (T)         = %18.12lf [Eh]\n\n",e_disp202t);
    
  }

  e_disp22t_ = e_disp220t + e_disp202t;

  if (print_) {
    outfile->Printf("    Disp22 (T)          = %18.12lf [Eh]\n",e_disp22t_);
    
  }

  if (nat_orbs_t3_) {
    double scale = e_disp20_/e_no_disp20_;
    e_disp220t *= scale;
    e_disp202t *= scale;
    e_est_disp22t_ = e_disp220t + e_disp202t;

    if (print_) {
      outfile->Printf("\n    Est. Disp220 (T)    = %18.12lf [Eh]\n",e_disp220t);
      outfile->Printf("    Est. Disp202 (T)    = %18.12lf [Eh]\n\n",e_disp202t);
      outfile->Printf("    Est. Disp22 (T)     = %18.12lf [Eh]\n",e_est_disp22t_);
      
    }
  }
}

void SAPT2p::disp22tccd()
{
  if (print_) {
    outfile->Printf("\n");
  }

  if (nat_orbs_t3_) {
    natural_orbitalify_ccd();
  }

  double e_disp220t;


  if (nat_orbs_t3_) {
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
    outfile->Printf("\n    Disp220 (T)         = %18.12lf [Eh]\n\n",e_disp220t);
    
  }

  double e_disp202t;

  if (nat_orbs_t3_) {
    e_disp202t = disp220tccd(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals","SS NO RI Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR NO RI Integrals",PSIF_SAPT_CCD,"T BSBS Natorb Amplitudes","T ARBS Natorb Amplitudes",
      no_evalsB_,no_evalsA_,noccB_,no_nvirB_,foccB_,noccA_,no_nvirA_,foccA_);
  }
  else {
    e_disp202t = disp220tccd(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      PSIF_SAPT_BB_DF_INTS,"BS RI Integrals","SS RI Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR RI Integrals",PSIF_SAPT_CCD,"T BSBS Amplitudes","T ARBS Amplitudes",
      evalsB_,evalsA_,noccB_,nvirB_,foccB_,noccA_,nvirA_,foccA_);
  }

  if (print_) {
    outfile->Printf("\n    Disp202 (T)         = %18.12lf [Eh]\n\n",e_disp202t);
    
  }

  e_disp22t_ccd_ = e_disp220t + e_disp202t;

  if (print_) {
    outfile->Printf("    Disp22 (T)          = %18.12lf [Eh]\n",e_disp22t_ccd_);
    
  }

  if (nat_orbs_t3_) {
    double scale = e_disp20_/e_no_disp20_;
    e_disp220t *= scale;
    e_disp202t *= scale;
    e_est_disp22t_ccd_ = e_disp220t + e_disp202t;

    if (print_) {
      outfile->Printf("\n    Est. Disp220 (T)    = %18.12lf [Eh]\n",e_disp220t);
      outfile->Printf("    Est. Disp202 (T)    = %18.12lf [Eh]\n\n",e_disp202t);
      outfile->Printf("    Est. Disp22 (T)     = %18.12lf [Eh]\n",e_est_disp22t_ccd_);
      
    }
  }
}

void SAPT2p::natural_orbitalify_ccd()
{
  int occA = noccA_ - foccA_;
  int occB = noccB_ - foccB_;

  double **tARAR = block_matrix(occA*nvirA_,occA*nvirA_);

  psio_->read_entry(PSIF_SAPT_CCD,"T ARAR Amplitudes",(char *) tARAR[0],
    occA*nvirA_*occA*nvirA_*(ULI) sizeof(double));

  double **tARAr = block_matrix(occA*nvirA_,occA*no_nvirA_);

  C_DGEMM('N','N',occA*nvirA_*occA,no_nvirA_,nvirA_,
    1.0,tARAR[0],nvirA_,no_CA_[0],
    no_nvirA_,0.0,tARAr[0],no_nvirA_);

  free_block(tARAR);
  double **tArAr = block_matrix(occA*no_nvirA_,occA*no_nvirA_);

  for (int a=0; a<occA; a++) {
    C_DGEMM('T','N',no_nvirA_,occA*no_nvirA_,nvirA_,
      1.0,no_CA_[0],no_nvirA_,
      tARAr[a*nvirA_],occA*no_nvirA_,0.0,
      tArAr[a*no_nvirA_],occA*no_nvirA_);
  }

  free_block(tARAr);

  psio_->write_entry(PSIF_SAPT_CCD,"T ARAR Natorb Amplitudes",(char *)
    tArAr[0],occA*no_nvirA_*occA*no_nvirA_*(ULI) sizeof(double));

  free_block(tArAr);

  double **tBSBS = block_matrix(occB*nvirB_,occB*nvirB_);

  psio_->read_entry(PSIF_SAPT_CCD,"T BSBS Amplitudes",(char *) tBSBS[0],
    occB*nvirB_*occB*nvirB_*(ULI) sizeof(double));

  double **tBSBs = block_matrix(occB*nvirB_,occB*no_nvirB_);

  C_DGEMM('N','N',occB*nvirB_*occB,no_nvirB_,nvirB_,
    1.0,tBSBS[0],nvirB_,no_CB_[0],
    no_nvirB_,0.0,tBSBs[0],no_nvirB_);

  free_block(tBSBS);
  double **tBsBs = block_matrix(occB*no_nvirB_,occB*no_nvirB_);

  for (int b=0; b<occB; b++) {
    C_DGEMM('T','N',no_nvirB_,occB*no_nvirB_,nvirB_,
      1.0,no_CB_[0],no_nvirB_,
      tBSBs[b*nvirB_],occB*no_nvirB_,0.0,
      tBsBs[b*no_nvirB_],occB*no_nvirB_);
  }

  free_block(tBSBs);

  psio_->write_entry(PSIF_SAPT_CCD,"T BSBS Natorb Amplitudes",(char *)
    tBsBs[0],occB*no_nvirB_*occB*no_nvirB_*(ULI) sizeof(double));

  free_block(tBsBs);

  double **tARBS = block_matrix(occA*nvirA_,occB*nvirB_);

  psio_->read_entry(PSIF_SAPT_CCD,"T ARBS Amplitudes",(char *) tARBS[0],
    occA*nvirA_*occB*nvirB_*(ULI) sizeof(double));

  double **tARBs = block_matrix(occA*nvirA_,occB*no_nvirB_);

  C_DGEMM('N','N',occA*nvirA_*occB,no_nvirB_,nvirB_,
    1.0,tARBS[0],nvirB_,no_CB_[0],
    no_nvirB_,0.0,tARBs[0],no_nvirB_);

  free_block(tARBS);
  double **tArBs = block_matrix(occA*no_nvirA_,occB*no_nvirB_);

  for (int a=0; a<occA; a++) {
    C_DGEMM('T','N',no_nvirA_,occB*no_nvirB_,nvirA_,
      1.0,no_CA_[0],no_nvirA_,
      tARBs[a*nvirA_],occB*no_nvirB_,0.0,
      tArBs[a*no_nvirA_],occB*no_nvirB_);
  }

  free_block(tARBs);

  double **tBsAr = block_matrix(occB*no_nvirB_,occA*no_nvirA_);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<no_nvirA_; r1++,a1r1++) {
    for(int b1=0,b1s1=0; b1<occB; b1++) {
    for(int s1=0; s1<no_nvirB_; s1++,b1s1++) {
      tBsAr[b1s1][a1r1] = tArBs[a1r1][b1s1];
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,"T ARBS Natorb Amplitudes",(char *)
    tArBs[0],occA*no_nvirA_*occB*no_nvirB_*(ULI) sizeof(double));
  psio_->write_entry(PSIF_SAPT_CCD,"T BSAR Natorb Amplitudes",(char *)
    tBsAr[0],occA*no_nvirA_*occB*no_nvirB_*(ULI) sizeof(double));

  free_block(tArBs);
  free_block(tBsAr);
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
    outfile->Printf("    (i = %3d of %3d) %10ld seconds\n",b+1,aoccB,
      stop-start);
  }  }

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

double SAPT2p::disp220tccd(int AAnum, const char *AA_label, int Rnum, const char *AR_label,
  const char *RR_label, int BBnum, const char *BS_label, int ampnum, const char *tarar,
  const char *tbsar, double *evalsA, double *evalsB, int noccA, int nvirA,
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
    outfile->Printf("    (i = %3d of %3d) %10ld seconds\n",b+1,noccB,stop-start);
  
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
