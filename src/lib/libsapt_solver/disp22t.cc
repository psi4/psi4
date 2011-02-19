/* This function calculates the Disp20 energy */

#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <cstring>
#include <iostream>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt2p.h"

namespace psi { namespace sapt {

void SAPT2p::disp22t()
{
  double d220t,d202t;

  if (params_.print) {
    fprintf(outfile,"Begining Disp220(T) Calculation\n\n");
    fflush(outfile);
  }

  if (params_.nat_orbs) {
    psio_->open(PSIF_SAPT_TEMP,0);
    natural_orbitalify_triples(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR NO Integrals","RR NO Integrals",PSIF_SAPT_BB_DF_INTS,
      "BS NO Integrals","T ARAR Amplitudes","T BSAR Amplitudes",
      no_info_.evalsA,no_info_.evalsB,no_info_.CA,no_info_.CB,
      calc_info_.noccA,calc_info_.nvirA,params_.foccA,no_info_.nvirA,
      calc_info_.noccB,calc_info_.nvirB,params_.foccB,no_info_.nvirB);
    d220t = disp220t(PSIF_SAPT_TEMP,"AA RI Integrals","AR NO Integrals",
      "RR NO Integrals",PSIF_SAPT_TEMP,"BS NO Integrals",PSIF_SAPT_TEMP,
      "T ARAR Amplitudes","T BSAR Amplitudes",no_info_.evalsA,
      no_info_.evalsB,calc_info_.noccA,no_info_.nvirA,params_.foccA,
      calc_info_.noccB,no_info_.nvirB,params_.foccB);
    psio_->close(PSIF_SAPT_TEMP,0);
  }
  else if (params_.foccA || params_.foccB) {
    psio_->open(PSIF_SAPT_TEMP,0);
    fzn_triples(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
      "RR RI Integrals",PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_SAPT_AMPS,
      "T ARAR Amplitudes","T BSAR Amplitudes",calc_info_.noccA,
      calc_info_.nvirA,params_.foccA,calc_info_.noccB,calc_info_.nvirB,
      params_.foccB);
    d220t = disp220t(PSIF_SAPT_TEMP,"AA RI Integrals","AR RI Integrals",
      "RR RI Integrals",PSIF_SAPT_TEMP,"BS RI Integrals",PSIF_SAPT_TEMP,
      "T ARAR Amplitudes","T BSAR Amplitudes",calc_info_.evalsA,
      calc_info_.evalsB,calc_info_.noccA,calc_info_.nvirA,params_.foccA,
      calc_info_.noccB,calc_info_.nvirB,params_.foccB);
    psio_->close(PSIF_SAPT_TEMP,0);
  }
  else {
    d220t = disp220t(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
      "RR RI Integrals",PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_SAPT_AMPS,
      "T ARAR Amplitudes","T BSAR Amplitudes",calc_info_.evalsA,
      calc_info_.evalsB,calc_info_.noccA,calc_info_.nvirA,params_.foccA,
      calc_info_.noccB,calc_info_.nvirB,params_.foccB);
  }

  if (params_.print) {
    fprintf(outfile,"\ndisp220(T)         = %18.12lf  H\n\n",d220t);
    fflush(outfile);
  }

  if (params_.print) {
    fprintf(outfile,"Begining Disp202(T) Calculation\n\n");
    fflush(outfile);
  }

  if (params_.nat_orbs) {
    psio_->open(PSIF_SAPT_TEMP,0);
    natural_orbitalify_triples(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS NO Integrals","SS NO Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR NO Integrals","T BSBS Amplitudes","T ARBS Amplitudes",
      no_info_.evalsB,no_info_.evalsA,no_info_.CB,no_info_.CA,
      calc_info_.noccB,calc_info_.nvirB,params_.foccB,no_info_.nvirB,
      calc_info_.noccA,calc_info_.nvirA,params_.foccA,no_info_.nvirA);
    d202t = disp220t(PSIF_SAPT_TEMP,"BB RI Integrals","BS NO Integrals",
      "SS NO Integrals",PSIF_SAPT_TEMP,"AR NO Integrals",PSIF_SAPT_TEMP,
      "T BSBS Amplitudes","T ARBS Amplitudes",no_info_.evalsB,
      no_info_.evalsA,calc_info_.noccB,no_info_.nvirB,params_.foccB,
      calc_info_.noccA,no_info_.nvirA,params_.foccA);
    psio_->close(PSIF_SAPT_TEMP,0);
  }
  else if (params_.foccA || params_.foccB) {
    psio_->open(PSIF_SAPT_TEMP,0);
    fzn_triples(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
      "SS RI Integrals",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_SAPT_AMPS,
      "T BSBS Amplitudes","T ARBS Amplitudes",calc_info_.noccB,
      calc_info_.nvirB,params_.foccB,calc_info_.noccA,calc_info_.nvirA,
      params_.foccA);
    d202t = disp220t(PSIF_SAPT_TEMP,"BB RI Integrals","BS RI Integrals",
      "SS RI Integrals",PSIF_SAPT_TEMP,"AR RI Integrals",PSIF_SAPT_TEMP,
      "T BSBS Amplitudes","T ARBS Amplitudes",calc_info_.evalsB,
      calc_info_.evalsA,calc_info_.noccB,calc_info_.nvirB,params_.foccB,
      calc_info_.noccA,calc_info_.nvirA,params_.foccA);
    psio_->close(PSIF_SAPT_TEMP,0);
  }
  else {
    d202t = disp220t(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
      "SS RI Integrals",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_SAPT_AMPS,
      "T BSBS Amplitudes","T ARBS Amplitudes",calc_info_.evalsB,
      calc_info_.evalsA,calc_info_.noccB,calc_info_.nvirB,params_.foccB,
      calc_info_.noccA,calc_info_.nvirA,params_.foccA);
  }

  if (params_.print) {
    fprintf(outfile,"\ndisp202(T)         = %18.12lf  H\n\n",d202t);
    fflush(outfile);
  }

  if (params_.nat_orbs) {
    no_info_.disp20 = nat_orb_disp20();
    double tval = results_.disp20/no_info_.disp20;
    d220t *= tval;
    d202t *= tval;
    if (params_.print) {
      fprintf(outfile,"est. disp220(T)    = %18.12lf  H\n",d220t);
      fprintf(outfile,"est. disp202(T)    = %18.12lf  H\n\n",d202t);
      fflush(outfile);
    }
  }

  results_.disp22t = d220t + d202t;
}

double SAPT2p::disp220t(int AAnum, char *AA_label, char *AR_label, 
  char *RR_label, int BBnum, char *BS_label, int ampnum, char *tarar, 
  char *tbsar, double *evalsA, double *evalsB, int noccA, int nvirA, 
  int foccA, int noccB, int nvirB, int foccB)
{
  double energy;

  noccA -= foccA;
  noccB -= foccB;

  double **w_ARAR = block_matrix(noccA*nvirA,noccA*nvirA);

  double **v_bsAA = block_matrix(noccA,noccA);
  double **v_bsRR = block_matrix(nvirA,nvirA);
  double **v_ARAA = block_matrix(noccA*nvirA,noccA*noccA);

  double **t_ARAR = read_IJKL(ampnum,tarar,noccA*nvirA,noccA*nvirA);
  double **t_bsAR = block_matrix(noccA,nvirA);

  double **B_p_AA = get_DF_ints(AAnum,AA_label,noccA*noccA);
  double **B_p_AR = get_DF_ints(AAnum,AR_label,noccA*nvirA);
  double **B_p_RR = get_DF_ints(AAnum,RR_label,nvirA*nvirA);
  double *B_p_bs = init_array(calc_info_.nrio);

  double **C_p_AR = block_matrix(noccA*nvirA,calc_info_.nrio);

  C_DGEMM('N','T',noccA*nvirA,noccA*noccA,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
    calc_info_.nrio,&(B_p_AA[0][0]),calc_info_.nrio,0.0,&(v_ARAA[0][0]),
    noccA*noccA);

  psio_address next_DF_BS = PSIO_ZERO;
  psio_address next_t_BSAR = PSIO_ZERO;

  time_t start = time(NULL);
  time_t stop;
  
  for(int b=0,bs=0; b<noccB; b++) {
  for(int s=0; s<nvirB; s++,bs++) {
  
    psio_->read(BBnum,BS_label,(char *) &(B_p_bs[0]),sizeof(double)*
      calc_info_.nrio,next_DF_BS,&next_DF_BS);
    psio_->read(ampnum,tbsar,(char *) t_bsAR[0],sizeof(double)*
      noccA*nvirA,next_t_BSAR,&next_t_BSAR);

    C_DGEMV('n',noccA*noccA,calc_info_.nrio,1.0,B_p_AA[0],calc_info_.nrio,
      B_p_bs,1,0.0,v_bsAA[0],1);
    C_DGEMV('n',nvirA*nvirA,calc_info_.nrio,1.0,B_p_RR[0],calc_info_.nrio,
      B_p_bs,1,0.0,v_bsRR[0],1);

    C_DGEMM('N','N',noccA*nvirA*noccA,nvirA,nvirA,1.0,&(t_ARAR[0][0]),nvirA,
            &(v_bsRR[0][0]),nvirA,0.0,&(w_ARAR[0][0]),nvirA);
    C_DGEMM('N','N',noccA,nvirA*noccA*nvirA,noccA,-1.0,&(v_bsAA[0][0]),noccA,
            &(t_ARAR[0][0]),nvirA*noccA*nvirA,1.0,&(w_ARAR[0][0]),
            nvirA*noccA*nvirA);
    C_DGEMM('N','N',noccA*nvirA*noccA,nvirA,noccA,-1.0,&(v_ARAA[0][0]),noccA,
            &(t_bsAR[0][0]),nvirA,1.0,&(w_ARAR[0][0]),nvirA);
    C_DGEMM('N','N',noccA,nvirA*calc_info_.nrio,nvirA,1.0,&(t_bsAR[0][0]),
            nvirA,&(B_p_RR[0][0]),nvirA*calc_info_.nrio,0.0,&(C_p_AR[0][0]),
            nvirA*calc_info_.nrio);
    C_DGEMM('N','T',noccA*nvirA,noccA*nvirA,calc_info_.nrio,1.0,
            &(B_p_AR[0][0]),calc_info_.nrio,&(C_p_AR[0][0]),calc_info_.nrio,
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
  fprintf(outfile,"(i = %3d) %10ld seconds\n",b,stop-start);
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

void SAPT2p::fzn_triples(int AAnum, char *AA_label, char *AR_label,
  char *RR_label, int BBnum, char *BS_label, int ampnum, char *tarar,
  char *tbsar, int noccA, int nvirA, int foccA, int noccB, int nvirB, 
  int foccB)
{
  psio_address next_psio;
  double **B_p_AA = block_matrix(noccA*noccA,calc_info_.nrio);

  psio_->read_entry(AAnum,AA_label,(char *) &(B_p_AA[0][0]),sizeof(double)*
      calc_info_.nrio*noccA*(ULI) noccA);

  next_psio = PSIO_ZERO;
  for(int a=foccA; a<noccA; a++) {
    int aa = a*noccA + foccA;
    psio_->write(PSIF_SAPT_TEMP,AA_label,(char *) &(B_p_AA[aa][0]),
      (noccA-foccA)*calc_info_.nrio*(ULI) sizeof(double),next_psio,&next_psio);
  }

  free_block(B_p_AA);

  double **B_p_AR = block_matrix(noccA*nvirA,calc_info_.nrio);

  psio_->read_entry(AAnum,AR_label,(char *) &(B_p_AR[0][0]),sizeof(double)*
      calc_info_.nrio*noccA*(ULI) nvirA);

  next_psio = PSIO_ZERO;
  for(int a=foccA; a<noccA; a++) {
    int ar = a*nvirA;
    psio_->write(PSIF_SAPT_TEMP,AR_label,(char *) &(B_p_AR[ar][0]),
      nvirA*calc_info_.nrio*(ULI) sizeof(double),next_psio,&next_psio);
  }

  free_block(B_p_AR);

  double **B_p_RR = block_matrix(nvirA*nvirA,calc_info_.nrio);

  psio_->read_entry(AAnum,RR_label,(char *) &(B_p_RR[0][0]),sizeof(double)*
      calc_info_.nrio*nvirA*(ULI) nvirA);

  psio_->write_entry(PSIF_SAPT_TEMP,RR_label,(char *) &(B_p_RR[0][0]),
    sizeof(double)*calc_info_.nrio*nvirA*(ULI) nvirA);

  free_block(B_p_RR);

  double **B_p_BS = block_matrix(noccB*nvirB,calc_info_.nrio);

  psio_->read_entry(BBnum,BS_label,(char *) &(B_p_BS[0][0]),sizeof(double)*
      calc_info_.nrio*noccB*(ULI) nvirB);

  next_psio = PSIO_ZERO;
  for(int b=foccB; b<noccB; b++) {
    int bs = b*nvirB;
    psio_->write(PSIF_SAPT_TEMP,BS_label,(char *) &(B_p_BS[bs][0]),
      nvirB*calc_info_.nrio*(ULI) sizeof(double),next_psio,&next_psio);
  }

  free_block(B_p_BS);

  double **tARAR = block_matrix(noccA*nvirA,noccA*nvirA);

  psio_->read_entry(ampnum,tarar,(char *) &(tARAR[0][0]),sizeof(double)*
      noccA*nvirA*noccA*(ULI) nvirA);

  next_psio = PSIO_ZERO;
  for(int a=foccA; a<noccA; a++) {
    for(int r=0; r<nvirA; r++) {
      int ar = a*nvirA+r;
      int aarr = foccA*nvirA;
      psio_->write(PSIF_SAPT_TEMP,tarar,(char *) &(tARAR[ar][aarr]),
        (noccA-foccA)*nvirA*(ULI) sizeof(double),next_psio,&next_psio);
  }}

  free_block(tARAR);

  double **tBSAR = block_matrix(noccB*nvirB,noccA*nvirA);

  psio_->read_entry(ampnum,tbsar,(char *) &(tBSAR[0][0]),sizeof(double)*
      noccB*nvirB*noccA*(ULI) nvirA);

  next_psio = PSIO_ZERO;
  for(int b=foccB; b<noccB; b++) {
    for(int s=0; s<nvirB; s++) {
      int bs = b*nvirB+s;
      int ar = foccA*nvirA;
      psio_->write(PSIF_SAPT_TEMP,tbsar,(char *) &(tBSAR[bs][ar]),
        (noccA-foccA)*nvirA*(ULI) sizeof(double),next_psio,&next_psio);
  }}

  free_block(tBSAR);

}

void SAPT2p::natural_orbitalify_triples(int AAnum, char *AA_label, 
  char *AR_label, char *RR_label, int BBnum, char *BS_label, char *tarar,
  char *tbsar, double *evalsA, double *evalsB, double **mo2noA, 
  double **mo2noB, int noccA, int nvirA, int foccA, int novirA, int noccB, 
  int nvirB, int foccB, int novirB)
{
  psio_address next_psio;
  double **B_p_AA = block_matrix(noccA*noccA,calc_info_.nrio);
  
  psio_->read_entry(AAnum,AA_label,(char *) &(B_p_AA[0][0]),sizeof(double)*
      calc_info_.nrio*noccA*(ULI) noccA);

  next_psio = PSIO_ZERO;
  for(int a=foccA; a<noccA; a++) {
    int aa = a*noccA + foccA;
    psio_->write(PSIF_SAPT_TEMP,AA_label,(char *) &(B_p_AA[aa][0]),
      (noccA-foccA)*calc_info_.nrio*(ULI) sizeof(double),next_psio,&next_psio);
  }

  free_block(B_p_AA);

  double **B_p_RR = block_matrix(novirA*novirA,calc_info_.nrio);

  psio_->read_entry(AAnum,RR_label,(char *) &(B_p_RR[0][0]),sizeof(double)*
      calc_info_.nrio*novirA*(ULI) novirA);

  psio_->write_entry(PSIF_SAPT_TEMP,RR_label,(char *) &(B_p_RR[0][0]),
    sizeof(double)*calc_info_.nrio*novirA*(ULI) novirA);

  free_block(B_p_RR);

  double **C_p_AR = block_matrix((noccA-foccA)*novirA,calc_info_.nrio);

  next_psio = psio_get_address(PSIO_ZERO,foccA*novirA*calc_info_.nrio*(ULI) 
    sizeof(double));
  psio_->read(AAnum,AR_label,(char *) &(C_p_AR[0][0]),(noccA-foccA)*
    novirA*calc_info_.nrio*(ULI) sizeof(double),next_psio,&next_psio);

  psio_->write_entry(PSIF_SAPT_TEMP,AR_label,(char *) &(C_p_AR[0][0]),
    (noccA-foccA)*novirA*calc_info_.nrio*(ULI) sizeof(double));

  double **C_p_BS = block_matrix((noccB-foccB)*novirB,calc_info_.nrio);

  next_psio = psio_get_address(PSIO_ZERO,foccB*novirB*calc_info_.nrio*(ULI)
    sizeof(double));
  psio_->read(BBnum,BS_label,(char *) &(C_p_BS[0][0]),(noccB-foccB)*
    novirB*calc_info_.nrio*(ULI) sizeof(double),next_psio,&next_psio);

  psio_->write_entry(PSIF_SAPT_TEMP,BS_label,(char *) &(C_p_BS[0][0]),
    (noccB-foccB)*novirB*calc_info_.nrio*(ULI) sizeof(double));

  double **tARAR = block_matrix((noccA-foccA)*novirA,(noccA-foccA)*novirA);

  C_DGEMM('N','T',(noccA-foccA)*novirA,(noccA-foccA)*novirA,calc_info_.nrio,
    1.0,&(C_p_AR[0][0]),calc_info_.nrio,&(C_p_AR[0][0]),calc_info_.nrio,0.0,
    &(tARAR[0][0]),(noccA-foccA)*novirA);

  for (int a=0, ar=0; a < (noccA-foccA); a++) {
  for (int r=0; r < novirA; r++, ar++) {
    for (int aa=0, aarr=0; aa < (noccA-foccA); aa++) {
    for (int rr=0; rr < novirA; rr++, aarr++) {
      double denom = evalsA[a+foccA]+evalsA[aa+foccA]-
        evalsA[r+noccA]-evalsA[rr+noccA];
      tARAR[ar][aarr] /= denom;
    }}
  }}

  psio_->write_entry(PSIF_SAPT_TEMP,tarar,(char *) &(tARAR[0][0]),
    (noccA-foccA)*novirA*(noccA-foccA)*novirA*(ULI) sizeof(double));

  free_block(tARAR);

  double **tBSAR = block_matrix((noccB-foccB)*novirB,(noccA-foccA)*novirA);

  C_DGEMM('N','T',(noccB-foccB)*novirB,(noccA-foccA)*novirA,calc_info_.nrio,
    1.0,&(C_p_BS[0][0]),calc_info_.nrio,&(C_p_AR[0][0]),calc_info_.nrio,0.0,
    &(tBSAR[0][0]),(noccA-foccA)*novirA);

  for (int b=0, bs=0; b < (noccB-foccB); b++) {
  for (int s=0; s < novirB; s++, bs++) {
    for (int a=0, ar=0; a < (noccA-foccA); a++) {
    for (int r=0; r < novirA; r++, ar++) {
      double denom = evalsA[a+foccA]+evalsB[b+foccB]-
        evalsA[r+noccA]-evalsB[s+noccB];
      tBSAR[bs][ar] /= denom;
    }}
  }}

  psio_->write_entry(PSIF_SAPT_TEMP,tbsar,(char *) &(tBSAR[0][0]),
    (noccB-foccB)*novirB*(noccA-foccA)*novirA*(ULI) sizeof(double));

  free_block(tBSAR);

  free_block(C_p_AR);
  free_block(C_p_BS);
}

double SAPT2p::nat_orb_disp20()
{
  double energy = 0.0;

  double **B_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    (char *) &(B_p_AR[0][0]),calc_info_.noccA*calc_info_.nvirA*
    calc_info_.nrio*(ULI) sizeof(double));

  double **C_p_AR = block_matrix(calc_info_.noccA*no_info_.nvirA,
    calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',no_info_.nvirA,calc_info_.nrio,calc_info_.nvirA,1.0,
      &(no_info_.CA[calc_info_.noccA][calc_info_.noccA]),
      calc_info_.noccA+no_info_.nvirA,B_p_AR[a*calc_info_.nvirA],
      calc_info_.nrio,0.0,C_p_AR[a*no_info_.nvirA],calc_info_.nrio);
  }

  free_block(B_p_AR);

  double **B_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    (char *) &(B_p_BS[0][0]),calc_info_.noccB*calc_info_.nvirB*
    calc_info_.nrio*(ULI) sizeof(double));

  double **C_p_BS = block_matrix(calc_info_.noccB*no_info_.nvirB,
    calc_info_.nrio);

  for(int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',no_info_.nvirB,calc_info_.nrio,calc_info_.nvirB,1.0,
      &(no_info_.CB[calc_info_.noccB][calc_info_.noccB]),
      calc_info_.noccB+no_info_.nvirB,B_p_BS[b*calc_info_.nvirB],
      calc_info_.nrio,0.0,C_p_BS[b*no_info_.nvirB],calc_info_.nrio);
  }

  free_block(B_p_BS);

  double **ARBS = block_matrix(calc_info_.noccA*no_info_.nvirA,
    calc_info_.noccB*no_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccA*no_info_.nvirA,calc_info_.noccB*
    no_info_.nvirB,calc_info_.nrio,1.0,&(C_p_AR[0][0]),calc_info_.nrio,
    &(C_p_BS[0][0]),calc_info_.nrio,0.0,&(ARBS[0][0]),calc_info_.noccB*
    no_info_.nvirB);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < no_info_.nvirA; r++, ar++) {
    for (int b=0, bs=0; b < calc_info_.noccB; b++) {
    for (int s=0; s < no_info_.nvirB; s++, bs++) {
      double denom = no_info_.evalsA[a]+no_info_.evalsB[b]-
        no_info_.evalsA[r+calc_info_.noccA]-
        no_info_.evalsB[s+calc_info_.noccB];
      double tval = ARBS[ar][bs];
      energy += tval*tval/denom;
    }}
  }}

  free_block(C_p_AR);
  free_block(C_p_BS);
  free_block(ARBS);

  return(4.0*energy);
}

}}
