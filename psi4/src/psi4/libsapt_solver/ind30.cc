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

#include "sapt2p3.h"



namespace psi { namespace sapt {

void SAPT2p3::ind30()
{
  double **tAR = block_matrix(noccA_,nvirA_);

  psio_->read_entry(PSIF_SAPT_AMPS,"Ind30 uAR Amplitudes", (char *) tAR[0],
    sizeof(double)*noccA_*nvirA_);

  double indA_B = 2.0*C_DDOT(noccA_*nvirA_,tAR[0],1,wBAR_[0],1);

  free_block(tAR);

  double **tBS = block_matrix(noccB_,nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"Ind30 uBS Amplitudes", (char *) tBS[0],
    sizeof(double)*noccB_*nvirB_);

  double indB_A = 2.0*C_DDOT(noccB_*nvirB_,tBS[0],1,wABS_[0],1);

  free_block(tBS);

  e_ind30_ = indA_B + indB_A;

  if (debug_) {
    outfile->Printf("\n    Ind30_1             = %18.12lf [Eh]\n",indA_B);
    outfile->Printf("    Ind30_2             = %18.12lf [Eh]\n",indB_A);
  }
  if (print_) {
    outfile->Printf("    Ind30               = %18.12lf [Eh]\n",e_ind30_);

  }
}

void SAPT2p3::ind30r()
{
  double indA_B = ind30r_1(CHFA_,CHFB_,wBAA_,wBRR_,PSIF_SAPT_AA_DF_INTS,
    "AA RI Integrals","AR RI Integrals","RR RI Integrals",
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",noccA_,nvirA_,noccB_,nvirB_);

  double indB_A = ind30r_1(CHFB_,CHFA_,wABB_,wASS_,PSIF_SAPT_BB_DF_INTS,
    "BB RI Integrals","BS RI Integrals","SS RI Integrals",
    PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",noccB_,nvirB_,noccA_,nvirA_);

  e_ind30r_ = indA_B + indB_A;

  if (debug_) {
    outfile->Printf("\n    Ind30_1             = %18.12lf [Eh]\n",indA_B);
    outfile->Printf("    Ind30_2             = %18.12lf [Eh]\n",indB_A);
  }
  if (print_) {
    outfile->Printf("    Ind30,r             = %18.12lf [Eh]\n",e_ind30r_);

  }
}

double SAPT2p3::ind30r_1(double **cAR, double **cBS, double **wBAA,
  double **wBRR, int intfileA, const char *AAlabel, const char *ARlabel,
  const char *RRlabel, int intfileB, const char *BSlabel, int noccA,
  int nvirA, int noccB, int nvirB)
{
  double energy = 0.0;

  double **xAR = block_matrix(noccA,nvirA);

  C_DGEMM('N','N',noccA,nvirA,nvirA,1.0,cAR[0],nvirA,wBRR[0],nvirA,
    0.0,xAR[0],nvirA);

  C_DGEMM('N','N',noccA,nvirA,noccA,-1.0,wBAA[0],noccA,cAR[0],nvirA,
    1.0,xAR[0],nvirA);

  energy = 2.0*C_DDOT(noccA*nvirA,cAR[0],1,xAR[0],1);

  free_block(xAR);

  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);

  double **B_p_BS = get_DF_ints(intfileB,BSlabel,0,noccB,0,nvirB);

  C_DGEMV('t',noccB*nvirB,ndf_+3,1.0,B_p_BS[0],ndf_+3,cBS[0],1,0.0,Y,1);

  free_block(B_p_BS);

  double **B_p_AR = get_DF_ints(intfileA,ARlabel,0,noccA,0,nvirA);

  C_DGEMV('t',noccA*nvirA,ndf_+3,1.0,B_p_AR[0],ndf_+3,cAR[0],1,0.0,X,1);

  energy += 8.0*C_DDOT(ndf_+3,X,1,Y,1);

  double **xAA = block_matrix(noccA,noccA);
  double **xRR = block_matrix(nvirA,nvirA);

  C_DGEMM('N','T',noccA,noccA,nvirA,1.0,cAR[0],nvirA,cAR[0],nvirA,
    0.0,xAA[0],noccA);

  C_DGEMM('T','N',nvirA,nvirA,noccA,1.0,cAR[0],nvirA,cAR[0],nvirA,
    0.0,xRR[0],nvirA);

  double **B_p_RR = get_DF_ints(intfileA,RRlabel,0,nvirA,0,nvirA);

  C_DGEMV('t',nvirA*nvirA,ndf_+3,1.0,B_p_RR[0],ndf_+3,xRR[0],1,0.0,Y,1);

  energy += 8.0*C_DDOT(ndf_+3,X,1,Y,1);

  double **C_p_AR = block_matrix(noccA*nvirA,ndf_+3);

  C_DGEMM('N','N',noccA,nvirA*(ndf_+3),nvirA,1.0,cAR[0],nvirA,B_p_RR[0],
    nvirA*(ndf_+3),0.0,C_p_AR[0],nvirA*(ndf_+3));

  free_block(B_p_RR);

  double **D_p_AR = block_matrix(noccA*nvirA,ndf_+3);

  for (int a=0; a<noccA; a++) {
    C_DGEMM('N','N',nvirA,ndf_+3,nvirA,1.0,xRR[0],nvirA,C_p_AR[a*nvirA],
      ndf_+3,0.0,D_p_AR[a*nvirA],ndf_+3);
  }

  energy -= 4.0*C_DDOT(noccA*nvirA*(ndf_+3),B_p_AR[0],1,D_p_AR[0],1);

  free_block(C_p_AR);
  free_block(D_p_AR);

  double **B_p_AA = get_DF_ints(intfileA,AAlabel,0,noccA,0,noccA);

  C_DGEMV('t',noccA*noccA,ndf_+3,1.0,B_p_AA[0],ndf_+3,xAA[0],1,0.0,Y,1);

  energy -= 8.0*C_DDOT(ndf_+3,X,1,Y,1);

  double **C_p_AA = block_matrix(noccA*noccA,ndf_+3);
  double **D_p_AA = block_matrix(noccA*noccA,ndf_+3);

  for (int a=0; a<noccA; a++) {
    C_DGEMM('N','N',noccA,ndf_+3,nvirA,1.0,cAR[0],nvirA,B_p_AR[a*nvirA],
      ndf_+3,0.0,C_p_AA[a*noccA],ndf_+3);
  }

  C_DGEMM('N','N',noccA,noccA*(ndf_+3),noccA,1.0,xAA[0],noccA,C_p_AA[0],
    noccA*(ndf_+3),0.0,D_p_AA[0],noccA*(ndf_+3));

  energy += 4.0*C_DDOT(noccA*noccA*(ndf_+3),B_p_AA[0],1,D_p_AA[0],1);

  free(X);
  free(Y);
  free_block(xAA);
  free_block(xRR);
  free_block(B_p_AA);
  free_block(C_p_AA);
  free_block(D_p_AA);
  free_block(B_p_AR);

  return(energy);
}

}}
