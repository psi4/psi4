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

#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::exch11()
{
  double e_exch110 = exch110(PSIF_SAPT_AMPS,"Theta AR Intermediates");

  if (debug_) {
    outfile->Printf("    Exch110             = %18.12lf [Eh]\n",e_exch110);
    
  }

  double e_exch101 = exch101(PSIF_SAPT_AMPS,"Theta BS Intermediates");

  if (debug_) {
    outfile->Printf("    Exch101             = %18.12lf [Eh]\n\n",e_exch101);
    
  }

  e_exch11_ = e_exch110 + e_exch101; 


  if (print_) {
    outfile->Printf("    Exch11              = %18.12lf [Eh]\n",e_exch11_);
    
  }
}

double SAPT2::exch110(int ampfile, const char *thetalabel)
{
  double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0;

  double **T_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
     sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  double **B_p_AB = get_AB_ints(2,foccA_,0);
  double **C_p_AB = block_matrix(aoccA_*noccB_,ndf_+3);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',noccB_,ndf_+3,nvirA_,1.0,&(sAB_[noccA_][0]),nmoB_,
      T_p_AR[a*nvirA_],ndf_+3,0.0,C_p_AB[a*noccB_],ndf_+3);
  } 

  e1 -= 2.0*C_DDOT((long int) aoccA_*noccB_*(ndf_+3),C_p_AB[0],1,
    B_p_AB[0],1);

  free_block(B_p_AB);

  double **C_p_BB = block_matrix(noccB_*noccB_,ndf_+3);

  C_DGEMM('T','N',noccB_,noccB_*(ndf_+3),aoccA_,1.0,&(sAB_[foccA_][0]),nmoB_,
    C_p_AB[0],noccB_*(ndf_+3),0.0,C_p_BB[0],noccB_*(ndf_+3));

  free_block(C_p_AB);

  double **B_p_BB = get_BB_ints(1);

  e2 += 4.0*C_DDOT((long int) noccB_*noccB_*(ndf_+3),B_p_BB[0],1,
    C_p_BB[0],1);

  free_block(B_p_BB);
  free_block(C_p_BB);

  double **B_p_RB = get_RB_ints(1);
  double **C_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);

  for (int r=0; r<nvirA_; r++) {
    C_DGEMM('N','N',aoccA_,ndf_+3,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
      B_p_RB[r*noccB_],ndf_+3,0.0,C_p_AR[r],nvirA_*(ndf_+3));
  }

  e3 -= 2.0*C_DDOT(aoccA_*nvirA_*(ndf_+3),T_p_AR[0],1,C_p_AR[0],1);

  free_block(B_p_RB);
  free_block(C_p_AR);

  double **xAR = block_matrix(aoccA_,nvirA_);
  double **yAR = block_matrix(aoccA_,nvirA_);

  C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,xAR[0],nvirA_);

  C_DGEMV('n',aoccA_*nvirA_,ndf_+3,1.0,T_p_AR[0],ndf_+3,diagBB_,1,0.0,
    yAR[0],1);

  e4 -= 8.0*C_DDOT(aoccA_*nvirA_,xAR[0],1,yAR[0],1);

  free_block(xAR);
  free_block(yAR);
  free_block(T_p_AR);

  if (debug_) {
    outfile->Printf("\n    Exch11_1            = %18.12lf [Eh]\n",e1);
    outfile->Printf("    Exch11_2            = %18.12lf [Eh]\n",e2);
    outfile->Printf("    Exch11_3            = %18.12lf [Eh]\n",e3);
    outfile->Printf("    Exch11_4            = %18.12lf [Eh]\n",e4);
    
  }

  return(e1+e2+e3+e4);
}

double SAPT2::exch101(int ampfile, const char *thetalabel)
{
  double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0;

  double **T_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_BS[0],
     sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  double **B_p_AB = get_AB_ints(1,0,foccB_);
  double **C_p_AB = block_matrix(noccA_*aoccB_,ndf_+3);

  for (int b=0; b<aoccB_; b++) {
    C_DGEMM('N','N',noccA_,ndf_+3,nvirB_,1.0,&(sAB_[0][noccB_]),nmoB_,
      T_p_BS[b*nvirB_],ndf_+3,0.0,C_p_AB[b],aoccB_*(ndf_+3));
  } 

  e1 -= 2.0*C_DDOT((long int) noccA_*aoccB_*(ndf_+3),C_p_AB[0],1,
    B_p_AB[0],1);

  free_block(B_p_AB);

  double **C_p_AA = block_matrix(noccA_*noccA_,ndf_+3);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',noccA_,ndf_+3,aoccB_,1.0,&(sAB_[0][foccB_]),nmoB_,
      C_p_AB[a*aoccB_],ndf_+3,0.0,C_p_AA[a*noccA_],ndf_+3);
  }

  free_block(C_p_AB);

  double **B_p_AA = get_AA_ints(1);

  e2 += 4.0*C_DDOT((long int) noccA_*noccA_*(ndf_+3),B_p_AA[0],1,
    C_p_AA[0],1);

  free_block(B_p_AA);
  free_block(C_p_AA);

  double **B_p_AS = get_AS_ints(1);
  double **C_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  C_DGEMM('T','N',aoccB_,nvirB_*(ndf_+3),noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    B_p_AS[0],nvirB_*(ndf_+3),0.0,C_p_BS[0],nvirB_*(ndf_+3));

  e3 -= 2.0*C_DDOT(aoccB_*nvirB_*(ndf_+3),T_p_BS[0],1,C_p_BS[0],1);

  free_block(B_p_AS);
  free_block(C_p_BS);

  double **xBS = block_matrix(aoccB_,nvirB_);
  double **yBS = block_matrix(aoccB_,nvirB_);

  C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,xBS[0],nvirB_);

  C_DGEMV('n',aoccB_*nvirB_,ndf_+3,1.0,T_p_BS[0],ndf_+3,diagAA_,1,0.0,
    yBS[0],1);

  e4 -= 8.0*C_DDOT(aoccB_*nvirB_,xBS[0],1,yBS[0],1);

  free_block(xBS);
  free_block(yBS);
  free_block(T_p_BS);

  if (debug_) {
    outfile->Printf("\n    Exch11_1            = %18.12lf [Eh]\n",e1);
    outfile->Printf("    Exch11_2            = %18.12lf [Eh]\n",e2);
    outfile->Printf("    Exch11_3            = %18.12lf [Eh]\n",e3);
    outfile->Printf("    Exch11_4            = %18.12lf [Eh]\n",e4);
    
  }

  return(e1+e2+e3+e4);
}

}}
