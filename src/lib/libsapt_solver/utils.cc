#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt2b.h"

namespace psi { namespace sapt {

void SAPT::zero_disk(int file, char *array, char *zero, int nri, int ijmax)
{
  psio_address next_PSIF = PSIO_ZERO;

  for (int ij=0; ij<ijmax; ij++) {
    psio_->write(file,array,zero,sizeof(double)*(ULI) nri,next_PSIF,&next_PSIF);
  }
}

double** SAPT::get_DF_ints(int filenum, char *label, int length)
{
  double **A = block_matrix(length,ribasis_->nbf()+3);
  psio_->read_entry(filenum,label,(char *) A[0],
                  sizeof(double)*length*(ULI) (ribasis_->nbf()+3));
  return(A);
}

double** SAPT2B::get_AA_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA*calc_info_.noccA,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.noccA*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int a=0; a<calc_info_.noccA; a++){
      int aa = a*calc_info_.noccA+a;
      A[aa][calc_info_.nrio-3] = 1.0;
      A[aa][calc_info_.nrio-1] = enuc;
      for (int ap=0; ap<calc_info_.noccA; ap++){
        int aap = a*calc_info_.noccA+ap;
        A[aap][calc_info_.nrio-2] = NB*calc_info_.VBAA[a][ap];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_diag_AA_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA,calc_info_.nrio);

  psio_address next_PSIF = PSIO_ZERO;
  for (int a=0; a<calc_info_.noccA; a++){
    psio_->read(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
      &(A[a][0]),sizeof(double)*(ULI) calc_info_.nrio,next_PSIF,&next_PSIF);
    next_PSIF = psio_get_address(next_PSIF,calc_info_.noccA*calc_info_.nrio*
      (ULI) sizeof(double));
    if (dress) {
      A[a][calc_info_.nrio-3] = 1.0;
            A[a][calc_info_.nrio-2] = NB*calc_info_.VBAA[a][a];
      A[a][calc_info_.nrio-1] = enuc;
    }
  }

  return(A);
}

double** SAPT2B::get_BB_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccB*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccB*calc_info_.noccB*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int b=0; b<calc_info_.noccB; b++){
      int bb = b*calc_info_.noccB+b;
      A[bb][calc_info_.nrio-2] = 1.0;
      A[bb][calc_info_.nrio-1] = enuc;
      for (int bp=0; bp<calc_info_.noccB; bp++){
        int bbp = b*calc_info_.noccB+bp;
        A[bbp][calc_info_.nrio-3] = NA*calc_info_.VABB[b][bp];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_diag_BB_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccB,calc_info_.nrio);

  psio_address next_PSIF = PSIO_ZERO;
  for (int b=0; b<calc_info_.noccB; b++){
    psio_->read(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
      &(A[b][0]),sizeof(double)*(ULI) calc_info_.nrio,next_PSIF,&next_PSIF);
    next_PSIF = psio_get_address(next_PSIF,calc_info_.noccB*calc_info_.nrio*
      (ULI) sizeof(double));
    if (dress) {
      A[b][calc_info_.nrio-3] = NA*calc_info_.VABB[b][b];
      A[b][calc_info_.nrio-2] = 1.0;
      A[b][calc_info_.nrio-1] = enuc;
    }
  }

  return(A);
}

double** SAPT2B::get_AB_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.noccB*
    (ULI) calc_info_.nrio);

  if (dress==1) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int b=0; b<calc_info_.noccB; b++){
        int ab = a*calc_info_.noccB+b;
        A[ab][calc_info_.nrio-3] = calc_info_.S_AB[a][b];
        A[ab][calc_info_.nrio-2] = NB*calc_info_.VBAB[a][b];
        A[ab][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][b];
      }
    }
  }
  else if (dress==2) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int b=0; b<calc_info_.noccB; b++){
        int ab = a*calc_info_.noccB+b;
        A[ab][calc_info_.nrio-3] = NA*calc_info_.VAAB[a][b];
        A[ab][calc_info_.nrio-2] = calc_info_.S_AB[a][b];
        A[ab][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][b];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_AS_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.nvirB*
    (ULI) calc_info_.nrio);

  if (dress==1) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int s=0; s<calc_info_.nvirB; s++){
        int as = a*calc_info_.nvirB+s;
        A[as][calc_info_.nrio-3] = calc_info_.S_AB[a][s+calc_info_.noccB];
        A[as][calc_info_.nrio-2] = NB*calc_info_.VBAB[a][s+calc_info_.noccB];
        A[as][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][s+calc_info_.noccB];
      }
    }
  }
  else if (dress==2) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int s=0; s<calc_info_.nvirB; s++){
        int as = a*calc_info_.nvirB+s;
        A[as][calc_info_.nrio-3] = NA*calc_info_.VAAB[a][s+calc_info_.noccB];
        A[as][calc_info_.nrio-2] = calc_info_.S_AB[a][s+calc_info_.noccB];
        A[as][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][s+calc_info_.noccB];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_RB_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.nvirA*calc_info_.noccB*
    (ULI) calc_info_.nrio);

    if (dress == 1) {
    for (int r=0; r<calc_info_.nvirA; r++){
      for (int b=0; b<calc_info_.noccB; b++){
        int rb = r*calc_info_.noccB+b;
        A[rb][calc_info_.nrio-3] = NA*calc_info_.VAAB[r+calc_info_.noccA][b];
        A[rb][calc_info_.nrio-2] = calc_info_.S_AB[r+calc_info_.noccA][b];
        A[rb][calc_info_.nrio-1] = enuc*calc_info_.S_AB[r+calc_info_.noccA][b];
      }
    }
  }
  else if (dress == 2) {
    for (int r=0; r<calc_info_.nvirA; r++){
      for (int b=0; b<calc_info_.noccB; b++){
        int rb = r*calc_info_.noccB+b;
        A[rb][calc_info_.nrio-3] = calc_info_.S_AB[r+calc_info_.noccA][b];
        A[rb][calc_info_.nrio-2] = NB*calc_info_.VBAB[r+calc_info_.noccA][b];
        A[rb][calc_info_.nrio-1] = enuc*calc_info_.S_AB[r+calc_info_.noccA][b];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_AR_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirA*calc_info_.noccA,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.nvirA*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int r=0; r<calc_info_.nvirA; r++){
        int ar = a*calc_info_.nvirA+r;
        A[ar][calc_info_.nrio-2] = NB*calc_info_.VBAA[a][r+calc_info_.noccA];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_BS_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirB*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccB*calc_info_.nvirB*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int b=0; b<calc_info_.noccB; b++){
      for (int s=0; s<calc_info_.nvirB; s++){
        int bs = b*calc_info_.nvirB+s;
        A[bs][calc_info_.nrio-3] = NA*calc_info_.VABB[b][s+calc_info_.noccB];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_RR_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirA*calc_info_.nvirA,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.nvirA*calc_info_.nvirA*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int r=0; r<calc_info_.nvirA; r++){
      int rr = r*calc_info_.nvirA+r;
      A[rr][calc_info_.nrio-3] = 1.0;
      A[rr][calc_info_.nrio-1] = enuc;
      for (int rp=0; rp<calc_info_.nvirA; rp++){
        int rrp = r*calc_info_.nvirA+rp;
        A[rrp][calc_info_.nrio-2] =
          NB*calc_info_.VBAA[r+calc_info_.noccA][rp+calc_info_.noccA];
      }
    }
  }

  return(A);

}

double** SAPT2B::get_SS_ints(int dress)
{

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirB*calc_info_.nvirB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.nvirB*calc_info_.nvirB*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int s=0; s<calc_info_.nvirB; s++){
      int ss = s*calc_info_.nvirB+s;
      A[ss][calc_info_.nrio-2] = 1.0;
      A[ss][calc_info_.nrio-1] = enuc;
      for (int sp=0; sp<calc_info_.nvirB; sp++){
        int ssp = s*calc_info_.nvirB+sp;
        A[ssp][calc_info_.nrio-3] =
          NA*calc_info_.VABB[s+calc_info_.noccB][sp+calc_info_.noccB];
      }
    }
  }

  return(A);

}

double **SAPT::read_IJKL(int filenum, char *label, int length_IJ,
  int length_KL)
{
  double **A = block_matrix(length_IJ,length_KL);

  psio_->read_entry(filenum,label,(char *) &(A[0][0]),
                  sizeof(double)*length_IJ*length_KL);

  return(A);
}

void SAPT::write_IJKL(double **A, int filenum, char *label, int length_IJ,
  int length_KL)
{
  psio_->write_entry(filenum,label,(char *) &(A[0][0]),
                  sizeof(double)*length_IJ*length_KL);

  free_block(A);
}

double **SAPT::IJKL_ints(int IJfile, char *IJlabel, int IJlength, int KLfile,
  char *KLlabel, int KLlength)
{
  double **IJKL = block_matrix(IJlength, KLlength);
  double **DF_p_IJ = get_DF_ints(IJfile, IJlabel, IJlength);
  double **DF_p_KL = get_DF_ints(KLfile, KLlabel, KLlength);

  C_DGEMM('N','T',IJlength,KLlength,ribasis_->nbf()+3,1.0,&(DF_p_IJ[0][0]),
    ribasis_->nbf()+3,&(DF_p_KL[0][0]),ribasis_->nbf()+3,0.0,&(IJKL[0][0]),
    KLlength);

  free_block(DF_p_IJ);
  free_block(DF_p_KL);

  return(IJKL);
}

double **SAPT::IJIJ_ints(int IJfile, char *IJlabel, int IJlength)
{
  double **IJIJ = block_matrix(IJlength, IJlength);
  double **DF_p_IJ = get_DF_ints(IJfile, IJlabel, IJlength);

  C_DGEMM('N','T',IJlength,IJlength,ribasis_->nbf()+3,1.0,&(DF_p_IJ[0][0]),
    ribasis_->nbf()+3,&(DF_p_IJ[0][0]),ribasis_->nbf()+3,0.0,&(IJIJ[0][0]),
    IJlength);

  free_block(DF_p_IJ);

  return(IJIJ);
}

}}
